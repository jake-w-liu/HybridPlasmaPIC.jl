# Rank-owned kinetic ions and algebraic isothermal electrons; no serial grid/state.
struct DistributedHybrid{T,G,R,P,M,F}
    g::G
    physics::M
    nfloor::T
    fields::F
    owners::Vector{P}
    _species_scalars::Vector{Tuple{T,T}}
    _population_counts::Vector{Int}
    initialized::Base.RefValue{Bool}
    primed::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    _candidate::F
    _k1::NTuple{3,R}
    _k2::NTuple{3,R}
    _k3::NTuple{3,R}
    _k4::NTuple{3,R}
    _Btmp::NTuple{3,R}
    _lapJ::NTuple{3,R}
    _floor_count::Base.RefValue{Int}
    _staged::Vector{P}
    _midpoints::Vector{NTuple{3,Vector{T}}}
    _vpred::Vector{NTuple{3,Vector{T}}}
    _layout::_PencilParticleLayout
end

_pencil_hybrid_field_arrays(f) = (f.B..., f.E..., f.n, f.ui..., f.pe, f.J..., f.gradp..., f.ninv)
_pencil_hybrid_arrays(model) = (_pencil_hybrid_field_arrays(model.fields)...,
    _pencil_hybrid_field_arrays(model._candidate)..., model._k1..., model._k2...,
    model._k3..., model._k4..., model._Btmp..., model._lapJ...)
_pencil_hybrid_metadata(model) = (model._species_scalars, model._population_counts,
    model._layout.by_rank, model._layout.lookup, model._layout.partitions..., model._layout.starts...)
_pencil_hybrid_particle_work(model) = Tuple(a for rows in (model._midpoints, model._vpred) for row in rows for a in row)

function _pencil_hybrid_mode(::Type{T}) where {T}
    Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() || throw(
        ArgumentError("distributed hybrid requires nearest rounding and gradual underflow"))
    return nothing
end

function _pencil_hybrid_physics(physics, ::Type{T}) where {T}
    physics isa HybridPlasmaPIC.HybridModel{HybridPlasmaPIC.IsothermalElectrons{T}} ||
        throw(ArgumentError("distributed hybrid requires isothermal electrons at grid precision"))
    isfinite(physics.closure.Te) && HybridPlasmaPIC._stored_nonnegative(physics.closure.Te) ||
        throw(ArgumentError("distributed hybrid Te must be finite and nonnegative"))
    all(HybridPlasmaPIC._particle_stored_zero, (physics.η, physics.ηH, physics.de2)) ||
        throw(ArgumentError("distributed hybrid requires original η = ηH = de2 = 0"))
    return HybridPlasmaPIC._require_finite_positive_real("distributed hybrid nfloor", physics.nfloor, T)
end

function _pencil_hybrid_population(ps)
    _pencil_pic_population(ps)
    HybridPlasmaPIC._stored_positive(ps.q) ||
        throw(ArgumentError("distributed hybrid requires positive ion charge"))
    return nothing
end

# Representation only: restore does not inspect old field/clock/count values.
function _pencil_hybrid_storage(model::DistributedHybrid)
    ns = length(model.owners)
    ns > 0 && length(model._species_scalars) == ns && axes(model._population_counts, 1) == Base.OneTo(ns) ||
        throw(ArgumentError("distributed hybrid cached species/count storage is invalid"))
    model.initialized !== model.primed && model.last_dt !== model.time && model.step !== model._floor_count ||
        throw(ArgumentError("distributed hybrid phase/clock references must be independent"))
    _pencil_particle_layout_storage(model._layout, model.g)
    arrays = _pencil_hybrid_arrays(model)
    for (i, a) in enumerate(arrays)
        _require_physical_field("hybrid storage", a, model.g, i)
        for j = 1:(i-1)
            Base.mightalias(a, arrays[j]) &&
                throw(ArgumentError("distributed hybrid fields and scratch must be independent"))
        end
    end
    for rows in (model._midpoints, model._vpred, model._staged)
        isempty(rows) || length(rows) == ns ||
            throw(ArgumentError("distributed hybrid particle work rows must match species"))
    end
    metadata = (_pencil_hybrid_metadata(model)..., _pencil_hybrid_particle_work(model)...)
    for (i, a) in enumerate(metadata)
        _aliases_workspace(a, model.g) &&
            throw(ArgumentError("distributed hybrid bindings must not alias Fourier storage"))
        for b in arrays
            Base.mightalias(a, b) &&
                throw(ArgumentError("distributed hybrid bindings must not alias fields"))
        end
        for j = 1:(i-1)
            Base.mightalias(a, metadata[j]) &&
                throw(ArgumentError("distributed hybrid bindings and particle work must be independent"))
        end
    end
    return nothing
end

function _pencil_hybrid_species(model::DistributedHybrid{T,G,R,P}, species;
    bound = false, values = true, factors = true, counts = bound) where {T,G,R,P}
    species isa Vector{<:HybridPlasmaPIC.ParticleSet{3,T}} ||
        throw(ArgumentError("distributed hybrid species must be an ordinary Vector of matching 3D ParticleSets"))
    Base.require_one_based_indexing(species)
    length(species) == length(model.owners) || throw(ArgumentError("distributed hybrid species count changed"))
    arrays = (_pencil_hybrid_arrays(model)..., _pencil_hybrid_metadata(model)...,
        _pencil_hybrid_particle_work(model)...)
    for (s, ps) in enumerate(species)
        ps isa P || throw(ArgumentError("distributed hybrid requires Vector-backed particles"))
        HybridPlasmaPIC._check_particle_storage_axes(ps)
        HybridPlasmaPIC._require_particle_ordinary_stores(ps, "distributed hybrid")
        HybridPlasmaPIC._particle_components_mightalias(ps) &&
            throw(ArgumentError("distributed hybrid particle components must be independent"))
        for j = 1:(s-1)
            HybridPlasmaPIC._particle_sets_mightalias(ps, species[j]) &&
                throw(ArgumentError("distributed hybrid populations must be independent"))
        end
        for a in arrays
            HybridPlasmaPIC._mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("distributed hybrid particles must not alias model storage"))
        end
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            _aliases_workspace(a, model.g) &&
                throw(ArgumentError("distributed hybrid particles must not alias Fourier storage"))
        end
        for staged in model._staged
            HybridPlasmaPIC._particle_sets_mightalias(ps, staged) &&
                throw(ArgumentError("distributed hybrid particles must not alias staged storage"))
        end
        bound && ps !== model.owners[s] && throw(ArgumentError("distributed hybrid owner/order changed; reinitialize"))
        bound && factors && !isequal((ps.q, ps.m), model._species_scalars[s]) &&
            throw(ArgumentError("distributed hybrid factors changed; reinitialize"))
        counts && model._population_counts[s] != length(ps.weight) &&
            throw(ArgumentError("distributed hybrid local particle count changed; reinitialize"))
        values && _pencil_hybrid_population(ps)
    end
    return nothing
end

function _pencil_hybrid_bound_config(model::DistributedHybrid{T}) where {T}
    isequal(model.nfloor, _pencil_hybrid_physics(model.physics, T)) ||
        throw(ArgumentError("distributed hybrid converted density floor changed"))
    length(model._species_scalars) == length(model.owners) ||
        throw(ArgumentError("distributed hybrid cached species count is invalid"))
    for (q, m) in model._species_scalars
        HybridPlasmaPIC._validated_charge_mass(q, m, T)
        HybridPlasmaPIC._stored_positive(q) || throw(ArgumentError("distributed hybrid cached ion charge must be positive"))
    end
    return (T = T, n = model.g.n, L = model.g.L, dx = model.g.dx,
        Te = model.physics.closure.Te, original_nfloor = model.physics.nfloor, nfloor = model.nfloor,
        η = model.physics.η, ηH = model.physics.ηH, de2 = model.physics.de2,
        species_scalars = Tuple(model._species_scalars), owned_ranges = Tuple(model._layout.by_rank))
end

function _pencil_hybrid_density(n, nfloor)
    all(value -> isfinite(value) && value >= nfloor, n) || throw(
        ArgumentError("distributed hybrid represented density must be finite and at least nfloor; clipping is unsupported"))
    return nothing
end

function _pencil_hybrid_finite_fields(fields, nfloor)
    foreach(_finite_magnitude, _pencil_hybrid_field_arrays(fields))
    _pencil_hybrid_density(parent(fields.n), nfloor)
    return nothing
end

_pencil_hybrid_phase_values(initialized::Bool, primed::Bool, last_dt, time, step::Int) =
    _pencil_pic_phase_values(initialized, primed, last_dt, time, step)
_pencil_hybrid_phase(model) = _pencil_hybrid_phase_values(model.initialized[], model.primed[],
    model.last_dt[], model.time[], model.step[])

function _pencil_hybrid_agreement(model::DistributedHybrid{T}, species, context; ready = true) where {T}
    return _agree_configuration(model.g.fft.comm, context) do
        _pencil_hybrid_mode(T)
        _pencil_hybrid_storage(model)
        _pencil_hybrid_species(model, species; bound = ready)
        if ready
            model.initialized[] || throw(ArgumentError("initialize distributed hybrid before use"))
            _pencil_hybrid_phase(model)
            _pencil_hybrid_finite_fields(model.fields, model.nfloor)
        end
        (_pencil_hybrid_bound_config(model), Tuple((ps.q, ps.m) for ps in species), ready ?
            (model.initialized[], model.primed[], model.last_dt[], model.time[], model.step[]) : ())
    end
end

"""
    distributed_hybrid(grid, species, physics::HybridModel{<:IsothermalElectrons})

Collectively allocate periodic 3D CPU Float32/64 positive-ion CIC hybrid fields
on actual physical pencils. Te has grid precision; original η/ηH/de2 are zero.
All populations are mobile and use original q/m. Seed `model.fields.B`, then
`init!` with physical full-time velocities. Density must remain at least nfloor;
no floor, Poisson neutralizer or replicated field/population is introduced.

Calls require nearest rounding, gradual underflow, faithful Base alias reporting
and exclusive workspace use. The communicator/plan remains caller-owned. Ready
operations enforce cached counts; same-count external replacement is unsupported
and cannot be detected. See `docs/DISTRIBUTED_HYBRID_IMPLEMENTATION.md` for phase,
range, transaction, diagnostics and explicitly deferred execution acceptance.
"""
function HybridPlasmaPIC.distributed_hybrid(g::DistributedFourierGrid{T}, species, physics) where {T}
    context = "Distributed hybrid construction"
    layout = HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "ownership construction") do
        _pencil_particle_layout(g)
    end
    parameters = _agree_configuration(g.fft.comm, context) do
        _pencil_hybrid_mode(T)
        nf = _pencil_hybrid_physics(physics, T)
        species isa Vector{<:HybridPlasmaPIC.ParticleSet{3,T}} && !isempty(species) ||
            throw(ArgumentError("species must be a nonempty ordinary Vector of matching 3D ParticleSets"))
        Base.require_one_based_indexing(species)
        all(ps -> ps isa _PencilParticle{T}, species) ||
            throw(ArgumentError("distributed hybrid requires Vector-backed particles"))
        foreach(_pencil_hybrid_population, species)
        (T, g.n, g.L, g.dx, physics, nf, Tuple((ps.q, ps.m) for ps in species), Tuple(layout.by_rank))
    end
    return HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "model allocation") do
        vector() = ntuple(_ -> HybridPlasmaPIC.distributed_field(g), 3)
        fields() = (B = vector(), E = vector(), n = HybridPlasmaPIC.distributed_field(g),
            ui = vector(), pe = HybridPlasmaPIC.distributed_field(g), J = vector(),
            gradp = vector(), ninv = HybridPlasmaPIC.distributed_field(g))
        P = _PencilParticle{T}
        model = DistributedHybrid(g, physics, parameters[6], fields(), P[ps for ps in species],
            Tuple{T,T}[(ps.q, ps.m) for ps in species], Int[length(ps.weight) for ps in species],
            Ref(false), Ref(false), Ref(zero(T)), Ref(zero(T)), Ref(0), fields(),
            vector(), vector(), vector(), vector(), vector(), vector(), Ref(0),
            P[], NTuple{3,Vector{T}}[], NTuple{3,Vector{T}}[], layout)
        _pencil_hybrid_storage(model)
        _pencil_hybrid_species(model, species)
        model
    end
end

function _pencil_hybrid_stages!(model::DistributedHybrid{T}, species, context) where {T}
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "particle and field staging") do
        HybridPlasmaPIC._pic_prepare_species_stages!(model._staged, species)
        for ps in model._staged
            HybridPlasmaPIC.apply_periodic!(ps, ntuple(_ -> zero(T), 3), model.g.L)
        end
        for rows in (model._midpoints, model._vpred)
            while length(rows) < length(species)
                push!(rows, ntuple(_ -> T[], 3))
            end
            resize!(rows, length(species))
            for (ps, row) in zip(model._staged, rows), a in row
                resize!(a, length(ps.weight))
            end
        end
        for (out, source) in zip(_pencil_hybrid_field_arrays(model._candidate),
            _pencil_hybrid_field_arrays(model.fields))
            copyto!(parent(out), parent(source))
        end
        foreach(a -> fill!(parent(a), zero(T)), model._lapJ)
        model._staged
    end
end

function _pencil_hybrid_moments!(model::DistributedHybrid{T}, species, context; midpoints = nothing) where {T}
    f = model._candidate
    _pencil_deposit_cic!(f.n, species, model.g, model._layout, zero(T), 0, context; midpoints)
    for c = 1:3
        _pencil_deposit_cic!(f.ui[c], species, model.g, model._layout, zero(T), c, context; midpoints)
    end
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "strict density and ion normalization") do
        _pencil_hybrid_density(parent(f.n), model.nfloor)
        _finite_vector_magnitude(f.ui)
        HybridPlasmaPIC._normalize_ion_moments!(parent(f.n), map(parent, f.ui), model.nfloor)
        _finite_vector_magnitude(f.ui)
        nothing
    end
end

function _pencil_hybrid_ohm_prep!(model::DistributedHybrid{T}, context) where {T}
    f = model._candidate
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "isothermal pressure") do
        _pencil_hybrid_density(parent(f.n), model.nfloor)
        HybridPlasmaPIC.electron_pressure!(parent(f.pe), parent(f.n), model.physics.closure)
        _finite_magnitude(f.pe)
    end
    SpectralOperators.gradient!(f.gradp, f.pe, model.g)
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "pressure gradient and reciprocal") do
        _finite_vector_magnitude(f.gradp)
        HybridPlasmaPIC._ohm_ninv!(parent(f.ninv), parent(f.n), model.nfloor, model._floor_count)
        iszero(model._floor_count[]) || throw(ArgumentError("distributed hybrid density clipping is unsupported"))
        _finite_magnitude(f.ninv)
        nothing
    end
end

function _pencil_hybrid_ohm!(model::DistributedHybrid{T}, B, context; differentiated = false) where {T}
    f = model._candidate
    SpectralOperators.curl!(f.J, B, model.g)
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "Ohm field admission") do
        _finite_vector_magnitude(f.J)
        risk = HybridPlasmaPIC._ohm_field_rows!(map(parent, f.E), map(parent, f.ui),
            map(parent, B), map(parent, f.J), map(parent, model._lapJ), map(parent, f.gradp),
            parent(f.ninv), T(model.physics.η), T(model.physics.ηH))
        _finite_vector_magnitude(f.E)
        differentiated && risk && throw(ArgumentError(
            "distributed hybrid original-product Faraday range recovery is unsupported; " *
            "a finite completed local Ohm field does not settle its differentiated source"))
        nothing
    end
end

function _pencil_hybrid_rhs!(out, B, model, context)
    _pencil_hybrid_ohm!(model, B, context; differentiated = true)
    SpectralOperators.curl!(out, model._candidate.E, model.g)
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "Faraday rate admission") do
        HybridPlasmaPIC._faraday_negate!(map(parent, out))
        all(a -> all(isfinite, parent(a)), out) || throw(ArgumentError(
            "distributed hybrid requires a finite materialized Faraday rate; " *
            "completion through a later duration is unsupported"))
        nothing
    end
end

# Each local affine publication is agreed before the next collective RHS.
function _pencil_hybrid_rk4!(model, h, context)
    B = model._candidate.B
    half, full, final = HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "RK coefficient preparation") do
        HybridPlasmaPIC._RK4Affine(h, 2), HybridPlasmaPIC._RK4Affine(h, 1), HybridPlasmaPIC._RK4Affine(h, 6)
    end
    _pencil_hybrid_rhs!(model._k1, B, model, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "RK second stage") do
        for c = 1:3
            tmp, base, rate = parent(model._Btmp[c]), parent(B[c]), parent(model._k1[c])
            @. tmp = half(base, rate)
        end
        _finite_vector_magnitude(model._Btmp)
    end
    _pencil_hybrid_rhs!(model._k2, model._Btmp, model, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "RK third stage") do
        for c = 1:3
            tmp, base, rate = parent(model._Btmp[c]), parent(B[c]), parent(model._k2[c])
            @. tmp = half(base, rate)
        end
        _finite_vector_magnitude(model._Btmp)
    end
    _pencil_hybrid_rhs!(model._k3, model._Btmp, model, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "RK fourth stage") do
        for c = 1:3
            tmp, base, rate = parent(model._Btmp[c]), parent(B[c]), parent(model._k3[c])
            @. tmp = full(base, rate)
        end
        _finite_vector_magnitude(model._Btmp)
    end
    _pencil_hybrid_rhs!(model._k4, model._Btmp, model, context)
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "RK final field") do
        for c = 1:3
            base = parent(B[c])
            k1, k2, k3, k4 = parent(model._k1[c]), parent(model._k2[c]), parent(model._k3[c]), parent(model._k4[c])
            @. base = final(base, k1, k2, k3, k4)
        end
        _finite_vector_magnitude(B)
        nothing
    end
end

function _pencil_hybrid_predictor!(model::DistributedHybrid{T}, species, h, context) where {T}
    f = model._candidate
    electric = _pencil_gather(f.E, species, model.g, model._layout, context)
    magnetic = _pencil_gather(f.B, species, model.g, model._layout, context)
    predicted = HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "integer velocity predictor") do
        for (ps, vp, Ep, Bp) in zip(species, model._vpred, electric, magnetic)
            HybridPlasmaPIC._predict_half_kick!(vp, ps.v, Ep, Bp,
                HybridPlasmaPIC._validated_boris_factor(ps), h / 2, length(ps.weight))
            all(a -> all(isfinite, a), vp) ||
                throw(ArgumentError("distributed hybrid integer velocity predictor must be finite"))
        end
        [HybridPlasmaPIC.ParticleSet{3,T}(ps.x, vp, ps.weight, ps.id, ps.tag, ps.q, ps.m)
            for (ps, vp) in zip(species, model._vpred)]
    end
    _pencil_hybrid_moments!(model, predicted, context)
    _pencil_hybrid_ohm_prep!(model, context)
    return _pencil_hybrid_ohm!(model, f.B, context)
end

function _pencil_hybrid_final(model::DistributedHybrid{T}, species, prepared, context, dt, clock; reset = false) where {T}
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "publication preparation") do
        _pencil_hybrid_species(model, prepared; counts = false)
        _pencil_owned_particles(prepared, model.g, model._layout)
        _pencil_hybrid_finite_fields(model._candidate, model.nfloor)
        _pencil_hybrid_phase_values(true, !reset, reset ? zero(T) : dt, clock[1], clock[2])
        length(prepared) == length(species) || throw(ArgumentError("prepared hybrid species count changed"))
        (scalars = Tuple{T,T}[(ps.q, ps.m) for ps in prepared], counts = Int[length(ps.weight) for ps in prepared])
    end
end

# Complete ordinary stores only, after the final all-rank preparation barrier.
function _pencil_hybrid_commit!(model, species, prepared, binding, dt, clock; reset = false)
    for (ps, replacement) in zip(species, prepared)
        HybridPlasmaPIC._publish_particle_replacement!(ps, replacement)
    end
    for (out, source) in zip(_pencil_hybrid_field_arrays(model.fields), _pencil_hybrid_field_arrays(model._candidate))
        copyto!(parent(out), parent(source))
    end
    copyto!(model.owners, species)
    copyto!(model._species_scalars, binding.scalars)
    copyto!(model._population_counts, binding.counts)
    model.initialized[] = true
    model.primed[] = !reset
    model.last_dt[] = reset ? zero(dt) : dt
    model.time[], model.step[] = clock
    return model
end

"""
    init!(model::DistributedHybrid, species)

Collectively initialize from finite seeded B and physical full-time velocities.
Private periodic migration, q-weighted moments and canonical isothermal Ohm fields
complete before all eighteen physical fields, particles and reset phase publish.
Old E/derived field values may be stale. No velocity unpriming is inferred.
"""
function HybridPlasmaPIC.init!(model::DistributedHybrid{T}, species) where {T}
    context = "Distributed hybrid initialization"
    _pencil_hybrid_agreement(model, species, context; ready = false)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "seed magnetic field") do
        _finite_vector_magnitude(model.fields.B)
    end
    staged = _pencil_hybrid_stages!(model, species, context)
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    _pencil_hybrid_moments!(model, prepared, context)
    _pencil_hybrid_ohm_prep!(model, context)
    _pencil_hybrid_ohm!(model, model._candidate.B, context)
    clock = (zero(T), 0)
    binding = _pencil_hybrid_final(model, species, prepared, context, zero(T), clock; reset = true)
    return _pencil_hybrid_commit!(model, species, prepared, binding, zero(T), clock; reset = true)
end

"""
    step!(model::DistributedHybrid, species, dt; NB=1)

Collectively advance original-factor ions and NB frozen-midpoint RK4 magnetic
subcycles. Carried E is recentered with the canonical forward half-kick predictor.
The timestep history preserves the variable-duration half-phase. Nonzero steps
require finite resolved stages; differentiated original-product Ohm range events
are explicitly unsupported. No nonlinear stability certificate is implied.
Recoverable failure changes private work only. True stored-zero dt validates,
leaves work/phase unchanged and does not consume initial priming.
"""
function HybridPlasmaPIC.step!(model::DistributedHybrid{T}, species, dt; NB = 1) where {T}
    context = "Distributed hybrid step"
    _pencil_hybrid_agreement(model, species, context)
    h, subcycles, clock = _agree_configuration(model.g.fft.comm, context) do
        dt isa Real && NB isa Integer || throw(ArgumentError("distributed hybrid requires real dt and integer NB"))
        nB = HybridPlasmaPIC._require_positive_intlike("NB", NB)
        h = HybridPlasmaPIC._validated_step_dt(T, dt, nB; min_NB = 1, name = context)
        if HybridPlasmaPIC._particle_stored_zero(h)
            h, nB, (model.time[], model.step[])
        else
            time, step = HybridPlasmaPIC._validated_step_clock(model.time[], model.step[], h)
            time > model.time[] || throw(ArgumentError("distributed hybrid time must advance"))
            h, nB, (time, step)
        end
    end
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "current particle ownership") do
        _pencil_owned_particles(species, model.g, model._layout)
    end
    HybridPlasmaPIC._particle_stored_zero(h) && return model
    staged = _pencil_hybrid_stages!(model, species, context)
    electric = _pencil_gather(model._candidate.E, staged, model.g, model._layout, context)
    magnetic = _pencil_gather(model._candidate.B, staged, model.g, model._layout, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "canonical ion kick and midpoint drift") do
        kick_dt = model.primed[] ? model.last_dt[] / 2 + h / 2 : h
        for (ps, mid, Ep, Bp) in zip(staged, model._midpoints, electric, magnetic)
            factor = HybridPlasmaPIC._validated_boris_factor(ps)
            qm = HybridPlasmaPIC._boris_native_factor(factor)
            ws = (Ep = Ep, Bp = Bp, xmid = mid)
            HybridPlasmaPIC._push_hybrid_gathered!(ps, ws, factor, qm, h / 2, h, kick_dt,
                !model.primed[], model.g.L)
            _pencil_hybrid_population(ps)
        end
    end
    _pencil_hybrid_moments!(model, staged, context; midpoints = model._midpoints)
    _pencil_hybrid_ohm_prep!(model, context)
    for _ = 1:subcycles
        _pencil_hybrid_rk4!(model, h / subcycles, context)
    end
    _pencil_hybrid_moments!(model, staged, context)
    _pencil_hybrid_ohm_prep!(model, context)
    _pencil_hybrid_ohm!(model, model._candidate.B, context)
    _pencil_hybrid_predictor!(model, staged, h, context)
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    binding = _pencil_hybrid_final(model, species, prepared, context, h, clock)
    return _pencil_hybrid_commit!(model, species, prepared, binding, h, clock)
end

# Local log maps stay local. Only exact rational enclosure endpoints cross ranks.
function _pencil_hybrid_free_energy(model::DistributedHybrid{T}, mechanical, context) where {T}
    comm = model.g.fft.comm
    terms, odd = HybridPlasmaPIC._mpi_collective_stage(comm, context, "electron log preparation") do
        terms = HybridPlasmaPIC._isothermal_log_terms(parent(model.fields.n), model.physics.closure.Te, model.g.dx)
        odd = false
        for (base, coefficient) in terms
            base > 1 && !iszero(coefficient) || throw(ArgumentError("invalid hybrid log term"))
            if base != 2
                isodd(base) && coefficient > 0 ||
                    throw(ArgumentError("hybrid endpoint odd log terms must have positive coefficients"))
                odd = true
            end
        end
        terms, odd
    end
    any_odd = MPI.Allreduce(Int(odd), MPI.MAX, comm) != 0
    logtwo = _pencil_exact_scalar_sum(comm, context, "log2 coefficient preparation", "log2 coefficient accumulation") do
        get(terms, UInt64(2), zero(_PencilCharge))
    end
    if !any_odd && iszero(logtwo)
        return HybridPlasmaPIC._mpi_collective_stage(comm, context, "zero electron energy publication") do
            zero(T), HybridPlasmaPIC._from_exact(T, mechanical)
        end
    end
    bits = precision(T) + 16
    while true
        local_bounds = HybridPlasmaPIC._mpi_collective_stage(comm, context, "local electron energy enclosure") do
            HybridPlasmaPIC._isothermal_log_bounds(terms, bits)
        end
        lower = _pencil_exact_scalar_sum(comm, context, "lower energy preparation", "lower energy accumulation") do
            local_bounds[1]
        end
        upper = _pencil_exact_scalar_sum(comm, context, "upper energy preparation", "upper energy accumulation") do
            local_bounds[2]
        end
        result = HybridPlasmaPIC._mpi_collective_stage(comm, context, "global energy rounding") do
            fe_lo, fe_hi = HybridPlasmaPIC._from_exact(T, lower), HybridPlasmaPIC._from_exact(T, upper)
            total_lo = HybridPlasmaPIC._from_exact(T, mechanical + lower)
            total_hi = HybridPlasmaPIC._from_exact(T, mechanical + upper)
            done = isequal(fe_lo, fe_hi) && isequal(total_lo, total_hi)
            (done = done, electron = fe_lo, total = total_lo, bits = done ? bits : Base.checked_mul(bits, 2))
        end
        result.done && return result.electron, result.total
        bits = result.bits
    end
end

"""
    hybrid_diagnostics(model::DistributedHybrid, species)

Collective original-factor particle totals, represented electron/ion charge and
volume/deposition ledgers, magnetic/free-electron reservoirs and represented divB.
Local canonical log maps are enclosed before global exact scalar endpoint sums;
no rank-rounded energy or global density/log dictionary is summed or replicated.
Variable precision and rational communication have unmeasured cost. Complete
reservoirs outside T can report Inf through canonical rounding. E² is absent.
Velocities retain the nominal negative last_dt/2 offset; fields/predictor moments
are integer-time. No exact discrete conservation or relative FFT bound is claimed.
Physical state is unchanged; private diagnostic scratch may change.
"""
function HybridPlasmaPIC.hybrid_diagnostics(model::DistributedHybrid{T}, species) where {T}
    context = "Distributed hybrid diagnostics"
    _pencil_hybrid_agreement(model, species, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "particle ownership") do
        _pencil_owned_particles(species, model.g, model._layout)
    end
    charge = _pencil_exact_scalar_sum(model.g.fft.comm, context, "ion charge preparation", "ion charge accumulation") do
        HybridPlasmaPIC._pic_collection_charges(model.g, zero(T), species).particle
    end
    grid_n = _pencil_exact_scalar_sum(model.g.fft.comm, context, "grid density preparation", "grid density accumulation") do
        HybridPlasmaPIC._exact_float_sum(parent(model.fields.n))
    end
    kinetic = _pencil_exact_scalar_sum(model.g.fft.comm, context, "kinetic preparation", "kinetic accumulation") do
        HybridPlasmaPIC._kinetic_energy_rational(species)
    end
    magnetic = _pencil_exact_scalar_sum(model.g.fft.comm, context, "magnetic preparation", "magnetic accumulation") do
        HybridPlasmaPIC._field_energy_rational((), map(parent, model.fields.B), one(_PencilCharge), model.g.dx)
    end
    local_momentum = HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "momentum preparation") do
        HybridPlasmaPIC._momentum_rational(species)
    end
    momentum = ntuple(3) do c
        _pencil_exact_scalar_sum(model.g.fft.comm, context, "momentum serialization", "momentum accumulation") do
            local_momentum[c]
        end
    end
    count = _pencil_exact_scalar_sum(model.g.fft.comm, context, "count preparation", "count accumulation") do
        sum((BigInt(length(ps.weight)) for ps in species); init = BigInt(0)) // BigInt(1)
    end
    mechanical = HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "mechanical energy preparation") do
        kinetic + magnetic
    end
    electron, total = _pencil_hybrid_free_energy(model, mechanical, context)
    SpectralOperators.divergence!(model._k1[1], model.fields.B, model.g)
    local_stats = HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "diagnostic local extrema") do
        _finite_magnitude(model._k1[1]), minimum(parent(model.fields.n)), maximum(parent(model.fields.n))
    end
    max_divB = MPI.Allreduce(local_stats[1], MPI.MAX, model.g.fft.comm)
    density_min = MPI.Allreduce(local_stats[2], MPI.MIN, model.g.fft.comm)
    density_max = MPI.Allreduce(local_stats[3], MPI.MAX, model.g.fft.comm)
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "diagnostic publication") do
        Q = HybridPlasmaPIC._exact_rational
        volume = prod(Q, model.g.dx)
        declared_volume, grid_volume = prod(Q, model.g.L), prod(BigInt, model.g.n) * volume
        electron_charge = -grid_n * volume
        offset = model.primed[] ? -model.last_dt[] / 2 : zero(T)
        (time = model.time[], step = model.step[], particle_count = numerator(count),
            particle_charge = charge, represented_electron_charge = electron_charge,
            particle_grid_charge_residual = charge + electron_charge,
            deposition_discrepancy = -electron_charge - charge,
            declared_volume = declared_volume, stored_grid_volume = grid_volume,
            volume_discrepancy = grid_volume - declared_volume,
            momentum = map(value -> HybridPlasmaPIC._from_exact(T, value), momentum),
            kinetic = HybridPlasmaPIC._from_exact(T, kinetic), magnetic = HybridPlasmaPIC._from_exact(T, magnetic),
            electron_free = electron, total = total, density_min = density_min, density_max = density_max,
            max_divB = max_divB, velocity_time_offset = offset, field_time_offset = zero(T))
    end
end

"""
    recommended_dt(model::DistributedHybrid; NB=1, safety=0.8, d_i=1, Ω_ci=1)

Collective existing RK4 magnetic advisory on actual grid spacings. This is not a
nonlinear hybrid/particle/ion-acoustic stability certificate or automatic limiter.
"""
function HybridPlasmaPIC.recommended_dt(model::DistributedHybrid{T}; NB = 1,
    safety = 0.8, d_i = 1.0, Ω_ci = 1.0) where {T}
    context = "Distributed hybrid timestep advisory"
    _pencil_hybrid_agreement(model, model.owners, context)
    return first(_agree_configuration(model.g.fft.comm, context) do
        NB isa Integer || throw(ArgumentError("distributed hybrid NB must be an integer"))
        nB = HybridPlasmaPIC._require_positive_intlike("NB", NB)
        value = HybridPlasmaPIC._recommended_dt_geometry(model.g, T; NB = nB,
            integrator = :rk4, safety, d_i, Ω_ci, η = 0.0, ηH = 0.0)
        (value,)
    end)
end
