# Actual physical-pencil electromagnetic PIC; no replicated fields or serial state.
struct DistributedEMPIC{T,G,R,P,S,F}
    g::G
    E::NTuple{3,R}
    B::NTuple{3,R}
    J::NTuple{3,R}
    rho_n::R
    rho_np1::R
    c::T
    background_charge::T
    neutrality::Symbol
    charge_rtol::T
    grid_dc::Symbol
    owners::Vector{P}
    _species_scalars::Vector{Tuple{T,T}}
    _population_counts::Vector{Int}
    initialized::Base.RefValue{Bool}
    primed::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    _candidate::F
    _Bc::NTuple{3,R}
    _Jraw::NTuple{3,R}
    _curl::NTuple{3,R}
    _source_spectral::S
    _staged::Vector{P}
    _midpoints::Vector{NTuple{3,Vector{T}}}
    _layout::_PencilParticleLayout
    _cfl::T
end

_pencil_empic_physical(model) =
    (E = model.E, B = model.B, J = model.J, rho_n = model.rho_n, rho_np1 = model.rho_np1)
_pencil_empic_field_arrays(fields) =
    (fields.E..., fields.B..., fields.J..., fields.rho_n, fields.rho_np1)
_pencil_empic_arrays(model) = (
    _pencil_empic_field_arrays(_pencil_empic_physical(model))...,
    _pencil_empic_field_arrays(model._candidate)...,
    model._Bc...,
    model._Jraw...,
    model._curl...,
)
_pencil_empic_metadata(model) = (
    model._species_scalars,
    model._population_counts,
    model._layout.by_rank,
    model._layout.lookup,
    model._layout.partitions...,
    model._layout.starts...,
)
_pencil_empic_midpoints(model) = Tuple(a for row in model._midpoints for a in row)

function _pencil_empic_mode(::Type{T}) where {T}
    Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() || throw(
        ArgumentError(
            "distributed electromagnetic PIC requires nearest rounding and gradual underflow",
        ),
    )
    return nothing
end

# Representation only: load may replace stale public values, phase and counts.
function _pencil_empic_storage(model::DistributedEMPIC{T}) where {T}
    ns = length(model.owners)
    length(model._species_scalars) == ns && axes(model._population_counts, 1) == Base.OneTo(ns) ||
        throw(ArgumentError("distributed electromagnetic cached species/count storage is invalid"))
    model.initialized !== model.primed && model.last_dt !== model.time || throw(
        ArgumentError("distributed electromagnetic phase and clock references must be independent"),
    )
    _pencil_particle_layout_storage(model._layout, model.g)
    arrays = _pencil_empic_arrays(model)
    for (i, a) in enumerate(arrays)
        _require_physical_field("electromagnetic storage", a, model.g, i)
        for j = 1:(i-1)
            Base.mightalias(a, arrays[j]) && throw(
                ArgumentError("distributed electromagnetic fields and scratch must be independent"),
            )
        end
    end
    source = model._source_spectral
    source isa PencilArrays.PencilArray{Complex{T},3} &&
        parent(source) isa Array{Complex{T},3} &&
        PencilArrays.pencil(source) === PencilArrays.pencil(model.g.spectral) &&
        axes(source) == axes(model.g.spectral) || throw(
        ArgumentError(
            "distributed electromagnetic source spectrum must use the grid's spectral pencil",
        ),
    )
    _aliases_workspace(source, model.g) && throw(
        ArgumentError(
            "distributed electromagnetic source spectrum must be independent of Fourier storage",
        ),
    )
    all(a -> !Base.mightalias(source, a), arrays) || throw(
        ArgumentError("distributed electromagnetic source spectrum must be independent of fields"),
    )
    isempty(model._midpoints) ||
        length(model._midpoints) == ns ||
        throw(ArgumentError("distributed electromagnetic midpoint rows must match species"))
    isempty(model._staged) ||
        length(model._staged) == ns ||
        throw(ArgumentError("distributed electromagnetic staged rows must match species"))
    metadata, midpoints = _pencil_empic_metadata(model), _pencil_empic_midpoints(model)
    for (i, a) in enumerate((metadata..., midpoints...))
        _aliases_workspace(a, model.g) && throw(
            ArgumentError("distributed electromagnetic bindings must not alias Fourier storage"),
        )
        for b in (arrays..., source)
            Base.mightalias(a, b) && throw(
                ArgumentError(
                    "distributed electromagnetic bindings must not alias fields or spectra",
                ),
            )
        end
        previous = (metadata..., midpoints...)
        for j = 1:(i-1)
            Base.mightalias(a, previous[j]) && throw(
                ArgumentError(
                    "distributed electromagnetic bindings and midpoints must be independent",
                ),
            )
        end
    end
    return nothing
end

function _pencil_empic_species(
    model::DistributedEMPIC{T,G,R,P},
    species;
    bound = false,
    values = true,
    factors = true,
    counts = bound,
) where {T,G,R,P}
    species isa Vector{<:HybridPlasmaPIC.ParticleSet{3,T}} || throw(
        ArgumentError(
            "distributed electromagnetic species must be an ordinary Vector of matching 3D ParticleSets",
        ),
    )
    Base.require_one_based_indexing(species)
    length(species) == length(model.owners) ||
        throw(ArgumentError("distributed electromagnetic species count changed"))
    arrays = (
        _pencil_empic_arrays(model)...,
        model._source_spectral,
        _pencil_empic_metadata(model)...,
        _pencil_empic_midpoints(model)...,
    )
    for (s, ps) in enumerate(species)
        ps isa P ||
            throw(ArgumentError("distributed electromagnetic PIC requires Vector-backed particles"))
        HybridPlasmaPIC._check_particle_storage_axes(ps)
        HybridPlasmaPIC._require_particle_ordinary_stores(ps, "distributed electromagnetic PIC")
        HybridPlasmaPIC._particle_components_mightalias(ps) && throw(
            ArgumentError("distributed electromagnetic particle components must be independent"),
        )
        for previous = 1:(s-1)
            HybridPlasmaPIC._particle_sets_mightalias(ps, species[previous]) &&
                throw(ArgumentError("distributed electromagnetic populations must be independent"))
        end
        for a in arrays
            HybridPlasmaPIC._mightalias_particle_storage(a, ps) && throw(
                ArgumentError("distributed electromagnetic particles must not alias model storage"),
            )
        end
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            _aliases_workspace(a, model.g) && throw(
                ArgumentError(
                    "distributed electromagnetic particles must not alias Fourier storage",
                ),
            )
        end
        for staged in model._staged
            HybridPlasmaPIC._particle_sets_mightalias(ps, staged) && throw(
                ArgumentError(
                    "distributed electromagnetic particles must not alias staged storage",
                ),
            )
        end
        bound &&
            ps !== model.owners[s] &&
            throw(ArgumentError("distributed electromagnetic owner/order changed; reinitialize"))
        bound &&
            factors &&
            !isequal((ps.q, ps.m), model._species_scalars[s]) &&
            throw(ArgumentError("distributed electromagnetic factors changed; reinitialize"))
        counts &&
            model._population_counts[s] != length(ps.weight) &&
            throw(
                ArgumentError(
                    "distributed electromagnetic local particle count changed; reinitialize",
                ),
            )
        values && _pencil_pic_population(ps)
    end
    return nothing
end

function _pencil_empic_bound_config(model::DistributedEMPIC{T}) where {T}
    HybridPlasmaPIC._require_finite_positive_real("c", model.c, T)
    isequal(model._cfl, _distributed_maxwell_cfl(model.g, model.c)) ||
        throw(ArgumentError("distributed electromagnetic cached CFL does not match the grid/speed"))
    model.neutrality in (:require, :neutralize) && model.grid_dc === :project ||
        throw(ArgumentError("distributed electromagnetic source policies are invalid"))
    isfinite(model.background_charge) &&
        isfinite(model.charge_rtol) &&
        HybridPlasmaPIC._stored_nonnegative(model.charge_rtol) ||
        throw(ArgumentError("distributed electromagnetic background or tolerance is invalid"))
    length(model._species_scalars) == length(model.owners) ||
        throw(ArgumentError("distributed electromagnetic cached species count is invalid"))
    for (q, m) in model._species_scalars
        HybridPlasmaPIC._validated_charge_mass(q, m, T)
    end
    return (
        T = T,
        n = model.g.n,
        L = model.g.L,
        dx = model.g.dx,
        c = model.c,
        background_charge = model.background_charge,
        neutrality = model.neutrality,
        charge_rtol = model.charge_rtol,
        grid_dc = model.grid_dc,
        species_scalars = Tuple(model._species_scalars),
        owned_ranges = Tuple(model._layout.by_rank),
    )
end

function _pencil_empic_finite_fields(fields)
    _finite_vector_magnitude(fields.E)
    _finite_vector_magnitude(fields.B)
    _finite_vector_magnitude(fields.J)
    _finite_magnitude(fields.rho_n)
    _finite_magnitude(fields.rho_np1)
    return nothing
end

function _pencil_empic_phase_values(
    initialized::Bool,
    primed::Bool,
    last_dt::T,
    time::T,
    step::Int,
    cfl::T,
) where {T}
    _pencil_pic_phase_values(initialized, primed, last_dt, time, step)
    primed &&
        !(last_dt < cfl) &&
        throw(
            ArgumentError("distributed electromagnetic timestep history must be strictly sub-CFL"),
        )
    return nothing
end

function _pencil_empic_phase_fields(J, rho_n, rho_np1, step::Int)
    if step == 0
        all(a -> all(HybridPlasmaPIC._particle_stored_zero, a), J) || throw(
            ArgumentError("unstepped distributed electromagnetic current must be stored zero"),
        )
        all(isequal(a, b) for (a, b) in zip(rho_n, rho_np1)) || throw(
            ArgumentError(
                "unstepped distributed electromagnetic density endpoints must be identical",
            ),
        )
    end
    return nothing
end

function _pencil_empic_phase(model)
    _pencil_empic_phase_values(
        model.initialized[],
        model.primed[],
        model.last_dt[],
        model.time[],
        model.step[],
        model._cfl,
    )
    _pencil_empic_phase_fields(
        map(parent, model.J),
        parent(model.rho_n),
        parent(model.rho_np1),
        model.step[],
    )
    return nothing
end

function _pencil_empic_agreement(
    model::DistributedEMPIC{T},
    species,
    context;
    ready = true,
) where {T}
    return _agree_configuration(model.g.fft.comm, context) do
        _pencil_empic_mode(T)
        _pencil_empic_storage(model)
        _pencil_empic_species(model, species; bound = ready)
        if ready
            model.initialized[] ||
                throw(ArgumentError("initialize distributed electromagnetic PIC before use"))
            _pencil_empic_phase(model)
            _pencil_empic_finite_fields(_pencil_empic_physical(model))
        end
        (
            _pencil_empic_bound_config(model),
            Tuple((ps.q, ps.m) for ps in species),
            ready ?
            (model.initialized[], model.primed[], model.last_dt[], model.time[], model.step[]) : (),
        )
    end
end

"""
    distributed_empic(grid, species; grid_dc, c=5, background_charge=0,
                      neutrality=:require, charge_rtol=64eps(T))

Collectively allocate periodic 3D CPU Float32/64 electromagnetic PIC on actual
physical pencils. All ordered Vector-backed populations are mobile, nonrelativistic
and CIC, with one common particle/field step. Every rank supplies every species,
including empty local populations; q/m and configuration agree, local counts need
not. Original finite signed charge and positive mass feed the canonical mover.
Explicit `grid_dc=:project` separates physical neutrality from represented-grid DC.

Seed owned E/B, then call `init_empic!` with physical full-time velocities. The
model owns no communicator or plan lifetime. Calls are collective, nonconcurrent,
nearest/gradual only. Migration replaces component leaves while owner objects stay
stable. Ready operations reject changed local counts; same-count external mutation
cannot be detected and is unsupported. Reinitialize after external population
changes. See `docs/DISTRIBUTED_EMPIC_IMPLEMENTATION.md` for phase, source, range,
publication and pending runtime/checkpoint acceptance.
"""
function HybridPlasmaPIC.distributed_empic(
    g::DistributedFourierGrid{T},
    species;
    grid_dc = nothing,
    c = 5,
    background_charge = zero(T),
    neutrality = :require,
    charge_rtol = 64eps(T),
) where {T}
    context = "Distributed electromagnetic PIC"
    layout = HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "ownership construction") do
        _pencil_particle_layout(g)
    end
    parameters = _agree_configuration(g.fft.comm, context) do
        _pencil_empic_mode(T)
        species isa Vector{<:HybridPlasmaPIC.ParticleSet{3,T}} || throw(
            ArgumentError("species must be an ordinary Vector of matching 3D ParticleSets"),
        )
        Base.require_one_based_indexing(species)
        all(ps -> ps isa _PencilParticle{T}, species) || throw(
            ArgumentError("distributed electromagnetic PIC requires Vector-backed particles"),
        )
        foreach(_pencil_pic_population, species)
        neutrality in (:require, :neutralize) ||
            throw(ArgumentError("neutrality must be :require or :neutralize"))
        grid_dc === :project || throw(ArgumentError("grid_dc must explicitly be :project"))
        c isa Real && background_charge isa Real && charge_rtol isa Real ||
            throw(ArgumentError("speed, background charge and tolerance must be real"))
        speed = HybridPlasmaPIC._require_finite_positive_real("c", c, T)
        background =
            HybridPlasmaPIC._require_finite_real("background_charge", background_charge, T)
        rtol = HybridPlasmaPIC._require_finite_nonnegative_real("charge_rtol", charge_rtol, T)
        T,
        g.n,
        g.L,
        g.dx,
        speed,
        background,
        neutrality,
        rtol,
        grid_dc,
        Tuple((ps.q, ps.m) for ps in species),
        Tuple(layout.by_rank)
    end
    return HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "model allocation") do
        vector() = ntuple(_ -> HybridPlasmaPIC.distributed_field(g), 3)
        fields() = (
            E = vector(),
            B = vector(),
            J = vector(),
            rho_n = HybridPlasmaPIC.distributed_field(g),
            rho_np1 = HybridPlasmaPIC.distributed_field(g),
        )
        physical = fields()
        P = _PencilParticle{T}
        model = DistributedEMPIC(
            g,
            physical.E,
            physical.B,
            physical.J,
            physical.rho_n,
            physical.rho_np1,
            parameters[5],
            parameters[6],
            parameters[7],
            parameters[8],
            parameters[9],
            P[ps for ps in species],
            Tuple{T,T}[(ps.q, ps.m) for ps in species],
            Int[length(ps.weight) for ps in species],
            Ref(false),
            Ref(false),
            Ref(zero(T)),
            Ref(zero(T)),
            Ref(0),
            fields(),
            vector(),
            vector(),
            vector(),
            HybridPlasmaPIC.distributed_fft_output(g.fft),
            P[],
            NTuple{3,Vector{T}}[],
            layout,
            _distributed_maxwell_cfl(g, parameters[5]),
        )
        _pencil_empic_storage(model)
        _pencil_empic_species(model, species)
        model
    end
end

HybridPlasmaPIC.maxwell_cfl(model::DistributedEMPIC) = model._cfl

# Finite, upward-rounded headroom for the shared represented modal constraint.
function _pencil_empic_projection_gain(
    g::DistributedFourierGrid{T},
    divisor::T,
    endpoints::Bool,
) where {T}
    minimum_wave = T(Inf)
    for axis in g.kvec, wave in axis
        iszero(wave) || (minimum_wave = min(minimum_wave, abs(wave)))
    end
    Q = HybridPlasmaPIC._exact_rational
    exact = Q(4)
    isinf(minimum_wave) || (exact += (endpoints ? 2 : 1) / (Q(divisor) * Q(minimum_wave)))
    gain = HybridPlasmaPIC._from_exact(T, exact)
    isfinite(gain) && Q(gain) < exact && (gain = nextfloat(gain))
    isfinite(gain) || throw(
        ArgumentError(
            "distributed electromagnetic projection requires a finite FFT multiplier bound; " *
            "this does not classify the complete physical field range",
        ),
    )
    return gain
end

function _pencil_empic_source_forward!(
    model::DistributedEMPIC{T},
    rho_n,
    rho_np1::Nothing,
    scale,
    context,
) where {T}
    return _forward_field!(model._source_spectral, rho_n, model.g, scale, context)
end

# Retain both finite endpoints through normalization, not a materialized T delta.
function _pencil_empic_source_forward!(
    model::DistributedEMPIC{T},
    rho_n,
    rho_np1,
    scale,
    context,
) where {T}
    g = model.g
    HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        context,
        "density difference normalization",
    ) do
        Q = HybridPlasmaPIC._exact_rational
        denominator = Q(scale)
        destination, preceding, following = parent(g.physical), parent(rho_n), parent(rho_np1)
        for I in eachindex(destination, preceding, following)
            delta = Q(following[I]) - Q(preceding[I])
            normalized = HybridPlasmaPIC._from_exact(T, delta / denominator)
            isfinite(normalized) ||
                throw(ArgumentError("normalized electromagnetic density difference is not finite"))
            !iszero(delta) &&
                HybridPlasmaPIC._particle_stored_zero(normalized) &&
                throw(
                    ArgumentError(
                        "electromagnetic density difference loses nonzero range during FFT normalization",
                    ),
                )
            destination[I] = Complex{T}(normalized)
        end
    end
    return HybridPlasmaPIC.distributed_fft_forward!(model._source_spectral, g.fft, g.physical)
end

# Gauss/current share transform orchestration, not a copied projection kernel.
function _pencil_empic_constraint!(
    model::DistributedEMPIC{T},
    out,
    raw,
    rho_n,
    rho_np1,
    divisor::T,
    ::Val{current},
    context,
) where {T,current}
    g = model.g
    magnitude, gain =
        HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "constraint preflight") do
            local_max = max(_finite_vector_magnitude(raw), _finite_magnitude(rho_n))
            current && (local_max = max(local_max, _finite_magnitude(rho_np1)))
            local_max, _pencil_empic_projection_gain(g, divisor, current)
        end
    scale = _agreed_scale(g, magnitude, gain)
    spectra = (g.spectral, g.spectral_work...)
    for c = 1:3
        _forward_field!(spectra[c], raw[c], g, scale, context)
    end
    _pencil_empic_source_forward!(model, rho_n, rho_np1, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "constraint modes") do
        all(a -> all(isfinite, parent(a)), (spectra..., model._source_spectral)) || throw(
            ArgumentError(
                "electromagnetic constraint requires finite represented Fourier coefficients",
            ),
        )
        components = map(PencilArrays.global_view, spectra)
        source = PencilArrays.global_view(model._source_spectral)
        dc = CartesianIndex(1, 1, 1)
        for I in CartesianIndices(source)
            wave = ntuple(d -> g.kvec[d][I[d]], 3)
            if current
                HybridPlasmaPIC._empic_current_mode!(components, wave, source, divisor, I, dc)
            elseif !HybridPlasmaPIC._spectral_mode_is_null(wave)
                original = ntuple(d -> components[d][I], 3)
                corrected =
                    HybridPlasmaPIC._spectral_project_mode(original, wave, source[I], divisor)
                for c = 1:3
                    components[c][I] = corrected[c]
                end
            end
        end
        all(a -> all(isfinite, parent(a)), spectra) || throw(
            ArgumentError("complete electromagnetic projected Fourier coefficient is not finite"),
        )
    end
    staging = _staging_buffers(g)
    for c = 1:3
        _inverse_field!(staging[c], spectra[c], g, scale, context)
    end
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "constraint output validation") do
        _finite_vector_magnitude(staging)
    end
    return _publish_components!(out, staging)
end

function _pencil_empic_stages!(model::DistributedEMPIC{T}, species, context) where {T}
    return HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "particle and field staging",
    ) do
        HybridPlasmaPIC._pic_prepare_species_stages!(model._staged, species)
        for ps in model._staged
            HybridPlasmaPIC.apply_periodic!(ps, ntuple(_ -> zero(T), 3), model.g.L)
        end
        while length(model._midpoints) < length(species)
            push!(model._midpoints, ntuple(_ -> T[], 3))
        end
        resize!(model._midpoints, length(species))
        for (ps, mid) in zip(model._staged, model._midpoints), a in mid
            resize!(a, length(ps.weight))
        end
        for (out, source) in zip(
            _pencil_empic_field_arrays(model._candidate),
            _pencil_empic_field_arrays(_pencil_empic_physical(model)),
        )
            copyto!(parent(out), parent(source))
        end
        model._staged
    end
end

function _pencil_empic_final(model::DistributedEMPIC{T}, species, prepared, context) where {T}
    return HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "publication preparation",
    ) do
        _pencil_empic_species(model, prepared; counts = false)
        _pencil_owned_particles(prepared, model.g, model._layout)
        _pencil_empic_finite_fields(model._candidate)
        length(prepared) == length(species) ||
            throw(ArgumentError("prepared electromagnetic species count changed"))
        (
            scalars = Tuple{T,T}[(ps.q, ps.m) for ps in prepared],
            counts = Int[length(ps.weight) for ps in prepared],
        )
    end
end

# No allocation/resize/validation occurs after this complete publication barrier.
function _pencil_empic_commit!(
    model,
    species,
    prepared,
    binding;
    reset = false,
    dt = model.last_dt[],
    clock = (model.time[], model.step[]),
)
    for (ps, replacement) in zip(species, prepared)
        HybridPlasmaPIC._publish_particle_replacement!(ps, replacement)
    end
    for (out, source) in zip(
        _pencil_empic_field_arrays(_pencil_empic_physical(model)),
        _pencil_empic_field_arrays(model._candidate),
    )
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
    init_empic!(model::DistributedEMPIC, species)

Collectively initialize from seeded finite E/B and physical full-time velocities.
Private copies wrap/migrate, form density and correct resolved Gauss modes while
preserving seed transverse/null E and B. Zero current, identical endpoint density
and reset phase/count bindings publish only after complete collective success.
"""
function HybridPlasmaPIC.init_empic!(model::DistributedEMPIC{T}, species) where {T}
    context = "Distributed electromagnetic initialization"
    _pencil_empic_agreement(model, species, context; ready = false)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "seed field validation") do
        _finite_vector_magnitude(model.E)
        _finite_vector_magnitude(model.B)
    end
    staged = _pencil_empic_stages!(model, species, context)
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    _pencil_pic_charges(model, prepared, context)
    fields = model._candidate
    _pencil_deposit_cic!(
        fields.rho_n,
        prepared,
        model.g,
        model._layout,
        model.background_charge,
        0,
        context,
    )
    _pencil_empic_constraint!(
        model,
        fields.E,
        fields.E,
        fields.rho_n,
        nothing,
        one(T),
        Val(false),
        context,
    )
    HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "initialized field completion",
    ) do
        copyto!(parent(fields.rho_np1), parent(fields.rho_n))
        foreach(a -> fill!(parent(a), zero(T)), fields.J)
        _pencil_empic_phase_fields(
            map(parent, fields.J),
            parent(fields.rho_n),
            parent(fields.rho_np1),
            0,
        )
    end
    binding = _pencil_empic_final(model, species, prepared, context)
    return _pencil_empic_commit!(
        model,
        species,
        prepared,
        binding;
        reset = true,
        dt = zero(T),
        clock = (zero(T), 0),
    )
end

function _pencil_empic_curl!(model, source, context)
    SpectralOperators.curl!(model._curl, source, model.g)
    HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "materialized curl admission",
    ) do
        all(a -> all(isfinite, parent(a)), model._curl) || throw(
            ArgumentError(
                "distributed electromagnetic PIC requires a finite materialized curl; " *
                "completion through a later timestep factor is not supported",
            ),
        )
    end
    return nothing
end

function _pencil_empic_magnetic!(out, old, curl, previous, following, backward = false)
    for c = 1:3
        HybridPlasmaPIC._em_magnetic_update!(
            parent(out[c]),
            parent(old[c]),
            parent(curl[c]),
            nothing,
            previous,
            following,
            backward,
        )
    end
    _finite_vector_magnitude(out)
    return nothing
end

# Includes original private first-step priming; public phase remains unchanged.
function _pencil_empic_advance_candidate!(
    model::DistributedEMPIC{T},
    species,
    h,
    next_time,
    next_step,
    context,
) where {T}
    staged = _pencil_empic_stages!(model, species, context)
    fields = model._candidate
    _pencil_deposit_cic!(
        fields.rho_n,
        staged,
        model.g,
        model._layout,
        model.background_charge,
        0,
        context,
    )
    electric = _pencil_gather(fields.E, staged, model.g, model._layout, context)
    _pencil_empic_curl!(model, fields.E, context)
    previous = model.primed[] ? model.last_dt[] : h
    if !model.primed[]
        magnetic = _pencil_gather(fields.B, staged, model.g, model._layout, context)
        HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "leapfrog priming") do
            for (ps, Ep, Bp) in zip(staged, electric, magnetic)
                HybridPlasmaPIC._prime_leapfrog!(
                    ps.v,
                    Ep,
                    Bp,
                    HybridPlasmaPIC._validated_boris_factor(ps),
                    h / 2,
                    length(ps.weight),
                    false,
                    model.c,
                )
            end
            _pencil_empic_magnetic!(fields.B, fields.B, model._curl, h, zero(T), true)
        end
    end
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "magnetic centering") do
        _pencil_empic_magnetic!(model._Bc, fields.B, model._curl, previous, zero(T))
    end
    magnetic = _pencil_gather(model._Bc, staged, model.g, model._layout, context)
    HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "canonical kick and midpoint drift",
    ) do
        for (ps, mid, Ep, Bp) in zip(staged, model._midpoints, electric, magnetic)
            HybridPlasmaPIC._push_empic_velocities!(
                ps,
                Ep,
                Bp,
                previous / 2 + h / 2,
                false,
                model.c,
            )
            HybridPlasmaPIC._drift_record_mid!(ps, mid, h, model.g.L)
            HybridPlasmaPIC.apply_periodic!(ps, ntuple(_ -> zero(T), 3), model.g.L)
            _pencil_pic_population(ps)
        end
    end
    for c = 1:3
        _pencil_deposit_cic!(
            model._Jraw[c],
            staged,
            model.g,
            model._layout,
            zero(T),
            c,
            context;
            midpoints = model._midpoints,
        )
    end
    _pencil_deposit_cic!(
        fields.rho_np1,
        staged,
        model.g,
        model._layout,
        model.background_charge,
        0,
        context,
    )
    _pencil_empic_constraint!(
        model,
        fields.J,
        model._Jraw,
        fields.rho_n,
        fields.rho_np1,
        h,
        Val(true),
        context,
    )
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "magnetic update") do
        _pencil_empic_magnetic!(fields.B, fields.B, model._curl, previous, h)
    end
    _pencil_empic_curl!(model, fields.B, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "electric update") do
        HybridPlasmaPIC._empic_electric_update!(
            map(parent, fields.E),
            map(parent, model._curl),
            map(parent, fields.J),
            model.c,
            h,
            model.c * model.c,
            nothing,
        )
        _pencil_empic_finite_fields(fields)
        _pencil_empic_phase_values(true, true, h, next_time, next_step, model._cfl)
    end
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    return prepared
end

"""
    step_empic!(model::DistributedEMPIC, species, dt)

Collectively advance all mobile populations and staggered Maxwell fields by one
nonnegative converted dt. Nonzero dt must advance the finite clock and be strictly
sub-CFL. Existing original-factor prime/Boris/midpoint drift and scalar Maxwell
owners are shared. Exact routed CIC moments precede one represented-mode current
constraint; periodic ownership migration precedes final atomic publication.
Recoverable failure preserves all public particles/fields/phase/count bindings;
private work may change. True stored-zero dt validates and leaves work unchanged.
"""
function HybridPlasmaPIC.step_empic!(model::DistributedEMPIC{T}, species, dt) where {T}
    context = "Distributed electromagnetic step"
    _pencil_empic_agreement(model, species, context)
    h, next_time, next_step = _agree_configuration(model.g.fft.comm, context) do
        dt isa Real || throw(ArgumentError("distributed electromagnetic dt must be real"))
        h = HybridPlasmaPIC._require_resolved_dt(
            context,
            HybridPlasmaPIC._validated_nonnegative_dt(T, dt; name = context),
        )
        if HybridPlasmaPIC._particle_stored_zero(h)
            h, model.time[], model.step[]
        else
            h < model._cfl || throw(
                ArgumentError(
                    "dt must be strictly below the distributed spectral Maxwell CFL limit",
                ),
            )
            time, step = HybridPlasmaPIC._validated_step_clock(model.time[], model.step[], h)
            time > model.time[] ||
                throw(ArgumentError("distributed electromagnetic time must advance"))
            h, time, step
        end
    end
    HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "current particle ownership",
    ) do
        _pencil_owned_particles(species, model.g, model._layout)
    end
    HybridPlasmaPIC._particle_stored_zero(h) && return model
    _pencil_pic_charges(model, species, context)
    prepared = _pencil_empic_advance_candidate!(model, species, h, next_time, next_step, context)
    binding = _pencil_empic_final(model, species, prepared, context)
    return _pencil_empic_commit!(
        model,
        species,
        prepared,
        binding;
        dt = h,
        clock = (next_time, next_step),
    )
end

# Private represented-source filtering for diagnostics; physical rho is retained.
function _pencil_empic_filter_source!(
    model::DistributedEMPIC{T},
    out,
    rho_n,
    rho_np1,
    keep_dc::Bool,
    context,
) where {T}
    g = model.g
    magnitude = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        context,
        "diagnostic source preparation",
    ) do
        value = _finite_magnitude(rho_n)
        rho_np1 === nothing ? value : max(value, _finite_magnitude(rho_np1))
    end
    scale = _agreed_scale(g, magnitude, rho_np1 === nothing ? one(T) : T(2))
    _pencil_empic_source_forward!(model, rho_n, rho_np1, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "diagnostic null-mode filtering") do
        all(isfinite, parent(model._source_spectral)) || throw(
            ArgumentError("electromagnetic diagnostic source requires finite Fourier coefficients"),
        )
        modes = PencilArrays.global_view(model._source_spectral)
        dc = CartesianIndex(1, 1, 1)
        for I in CartesianIndices(modes)
            wave = ntuple(d -> g.kvec[d][I[d]], 3)
            HybridPlasmaPIC._spectral_mode_is_null(wave) &&
                !(keep_dc && I == dc) &&
                (modes[I] = zero(Complex{T}))
        end
    end
    _inverse_field!(g.staging, model._source_spectral, g, scale, context)
    HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "diagnostic source validation") do
        _finite_magnitude(g.staging)
    end
    copyto!(parent(out), parent(g.staging))
    return nothing
end

function _pencil_empic_residuals(model::DistributedEMPIC{T}, context) where {T}
    g = model.g
    SpectralOperators.divergence!(model._curl[1], model.E, g)
    SpectralOperators.divergence!(model._curl[2], model.B, g)
    SpectralOperators.divergence!(model._curl[3], model.J, g)
    rho = model.step[] == 0 ? model.rho_n : model.rho_np1
    _pencil_empic_filter_source!(model, model._candidate.rho_np1, rho, nothing, false, context)
    if model.step[] > 0
        _pencil_empic_filter_source!(
            model,
            model._candidate.rho_n,
            model.rho_n,
            model.rho_np1,
            true,
            context,
        )
    end
    maxima = HybridPlasmaPIC._mpi_collective_stage(
        g.fft.comm,
        context,
        "diagnostic residual preparation",
    ) do
        _finite_vector_magnitude(model._curl)
        Q = HybridPlasmaPIC._exact_rational
        divE, divB, divJ = map(parent, model._curl)
        density, resolved = parent(rho), parent(model._candidate.rho_np1)
        gauss, unresolved, continuity, source =
            zero(_PencilCharge), zero(_PencilCharge), zero(_PencilCharge), zero(_PencilCharge)
        for I in eachindex(density, resolved, divE)
            gauss = max(gauss, abs(Q(divE[I]) - Q(resolved[I])))
            unresolved = max(unresolved, abs(Q(density[I]) - Q(resolved[I])))
        end
        if model.step[] > 0
            delta, h = parent(model._candidate.rho_n), Q(model.last_dt[])
            for I in eachindex(delta, divJ)
                value = Q(delta[I])
                continuity = max(continuity, abs(value + h * Q(divJ[I])))
                source = max(source, abs(value))
            end
        end
        HybridPlasmaPIC._from_exact(T, gauss),
        HybridPlasmaPIC._from_exact(T, unresolved),
        HybridPlasmaPIC._from_exact(T, continuity),
        HybridPlasmaPIC._from_exact(T, source),
        maximum(abs, divB)
    end
    # Rounding is monotone for these nonnegative maxima; no field is gathered.
    result = ntuple(i -> MPI.Allreduce(maxima[i], MPI.MAX, g.fft.comm), 5)
    return (
        gauss_residual_max = result[1],
        unresolved_density_max = result[2],
        continuity_residual_max = model.step[] == 0 ? nothing : result[3],
        continuity_source_max = model.step[] == 0 ? nothing : result[4],
        max_divB = result[5],
    )
end

"""
    electromagnetic_diagnostics(model::DistributedEMPIC, species)

Collective original-factor charge/count/momentum/energy and represented-grid DC,
volume and deposition ledgers. E/positions are full-time, velocity and B carry
nominal negative last_dt/2 offsets. Absolute resolved Gauss/continuity diagnostics
retain their declared finite FFT/materialized-derivative boundaries; unresolved
charge/DC and source size are reported separately. No exact conservation or
cancellation-relative FFT accuracy is asserted. Diagnostics change private work
only. All numerical/compiler/collective/resource acceptance remains deferred.
"""
function HybridPlasmaPIC.electromagnetic_diagnostics(model::DistributedEMPIC{T}, species) where {T}
    context = "Distributed electromagnetic diagnostics"
    _pencil_empic_agreement(model, species, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "particle ownership") do
        _pencil_owned_particles(species, model.g, model._layout)
    end
    charges = _pencil_pic_charges(model, species, context)
    grid_n = _pencil_exact_scalar_sum(
        model.g.fft.comm,
        context,
        "old grid charge preparation",
        "old grid charge accumulation",
    ) do
        HybridPlasmaPIC._exact_float_sum(parent(model.rho_n))
    end
    grid_np1 = _pencil_exact_scalar_sum(
        model.g.fft.comm,
        context,
        "new grid charge preparation",
        "new grid charge accumulation",
    ) do
        HybridPlasmaPIC._exact_float_sum(parent(model.rho_np1))
    end
    kinetic = _pencil_exact_scalar_sum(
        model.g.fft.comm,
        context,
        "kinetic preparation",
        "kinetic accumulation",
    ) do
        HybridPlasmaPIC._kinetic_energy_rational(species)
    end
    local_momentum =
        HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "momentum preparation") do
            HybridPlasmaPIC._momentum_rational(species)
        end
    momentum = ntuple(3) do c
        _pencil_exact_scalar_sum(
            model.g.fft.comm,
            context,
            "momentum serialization",
            "momentum accumulation",
        ) do
            local_momentum[c]
        end
    end
    field = _pencil_exact_scalar_sum(
        model.g.fft.comm,
        context,
        "field energy preparation",
        "field energy accumulation",
    ) do
        HybridPlasmaPIC._field_energy_rational(
            map(parent, model.E),
            map(parent, model.B),
            HybridPlasmaPIC._exact_rational(model.c)^2,
            model.g.dx,
        )
    end
    count = _pencil_exact_scalar_sum(
        model.g.fft.comm,
        context,
        "particle count preparation",
        "particle count accumulation",
    ) do
        sum((BigInt(length(ps.weight)) for ps in species); init = BigInt(0)) // BigInt(1)
    end
    residuals = _pencil_empic_residuals(model, context)
    return HybridPlasmaPIC._mpi_collective_stage(
        model.g.fft.comm,
        context,
        "diagnostic publication",
    ) do
        Q = HybridPlasmaPIC._exact_rational
        volume = prod(Q, model.g.dx)
        declared_volume, grid_volume = prod(Q, model.g.L), prod(BigInt, model.g.n) * volume
        grid_charge = (model.step[] == 0 ? grid_n : grid_np1) * volume
        physical_compensator = model.neutrality === :neutralize ? -charges.net : zero(_PencilCharge)
        offset = model.primed[] ? -model.last_dt[] / 2 : zero(T)
        merge(
            (
                time = model.time[],
                step = model.step[],
                particle_count = numerator(count),
                particle_charge = charges.particle,
                background_charge = charges.background,
                physical_neutrality_residual = charges.net,
                physical_policy_compensator = physical_compensator,
                deposited_grid_charge = grid_charge,
                grid_dc_compensator = -grid_charge,
                grid_dc_charge_change = (grid_np1 - grid_n) * volume,
                declared_volume = declared_volume,
                stored_grid_volume = grid_volume,
                background_volume_discrepancy = Q(model.background_charge) *
                                                (grid_volume - declared_volume),
                deposition_discrepancy = grid_charge - Q(model.background_charge) * grid_volume -
                                         charges.particle,
                momentum = map(value -> HybridPlasmaPIC._from_exact(T, value), momentum),
                kinetic = HybridPlasmaPIC._from_exact(T, kinetic),
                field = HybridPlasmaPIC._from_exact(T, field),
                total = HybridPlasmaPIC._from_exact(T, kinetic + field),
                velocity_time_offset = offset,
                magnetic_time_offset = offset,
            ),
            residuals,
        )
    end
end
