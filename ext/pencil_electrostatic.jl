struct DistributedElectrostaticPIC{T,G,R,P}
    g::G
    E::NTuple{3,R}
    rho::R
    background_charge::T
    neutrality::Symbol
    charge_rtol::T
    grid_dc::Symbol
    owners::Vector{P}
    _species_scalars::Vector{Tuple{T,T}}
    initialized::Base.RefValue{Bool}
    primed::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    _candidate_E::NTuple{3,R}
    _candidate_rho::R
    _staged::Vector{P}
    _layout::_PencilParticleLayout
end

function _pencil_es_mode(::Type{T}) where {T}
    Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() || throw(
        ArgumentError("distributed electrostatic PIC requires nearest rounding and gradual underflow"),
    )
    return nothing
end


function _pencil_es_storage(model)
    model.initialized !== model.primed && model.last_dt !== model.time ||
        throw(ArgumentError("distributed electrostatic phase and clock references must be independent"))
    length(model._species_scalars) == length(model.owners) ||
        throw(ArgumentError("distributed electrostatic cached species count is invalid"))
    _pencil_particle_layout_storage(model._layout, model.g)
    fields = (model.E..., model.rho, model._candidate_E..., model._candidate_rho)
    for (i, a) in enumerate(fields)
        _require_physical_field("electrostatic storage", a, model.g, i)
        for j = 1:(i-1)
            Base.mightalias(a, fields[j]) &&
                throw(ArgumentError("distributed electrostatic fields and candidates must be independent"))
        end
    end
    return nothing
end

function _pencil_es_species(model::DistributedElectrostaticPIC{T,G,R,P}, species; bound = false, values = true, factors = true) where {T,G,R,P}
    species isa Vector{<:HybridPlasmaPIC.ParticleSet{3,T}} ||
        throw(ArgumentError("distributed electrostatic species must be an ordinary Vector of matching 3D ParticleSets"))
    Base.require_one_based_indexing(species)
    length(species) == length(model.owners) ||
        throw(ArgumentError("distributed electrostatic species count changed"))
    fields = (model.E..., model.rho, model._candidate_E..., model._candidate_rho)
    for (s, ps) in enumerate(species)
        ps isa P || throw(ArgumentError("distributed electrostatic PIC requires Vector-backed particle storage"))
        HybridPlasmaPIC._check_particle_storage_axes(ps)
        HybridPlasmaPIC._require_particle_ordinary_stores(ps, "distributed electrostatic PIC")
        HybridPlasmaPIC._particle_components_mightalias(ps) &&
            throw(ArgumentError("distributed electrostatic particle components must be independent"))
        for previous = 1:(s-1)
            HybridPlasmaPIC._particle_sets_mightalias(ps, species[previous]) &&
                throw(ArgumentError("distributed electrostatic populations must be independent"))
        end
        for field in fields
            HybridPlasmaPIC._mightalias_particle_storage(field, ps) &&
                throw(ArgumentError("distributed electrostatic particles must not alias field storage"))
        end
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            _aliases_workspace(a, model.g) &&
                throw(ArgumentError("distributed electrostatic particles must not alias Fourier storage or metadata"))
        end
        for staged in model._staged
            HybridPlasmaPIC._particle_sets_mightalias(ps, staged) &&
                throw(ArgumentError("distributed electrostatic particles must not alias staged storage"))
        end
        if bound
            ps === model.owners[s] && (!factors || isequal((ps.q, ps.m), model._species_scalars[s])) ||
                throw(ArgumentError("distributed electrostatic owners or factors changed; reinitialize or refresh"))
        end
        values && _pencil_pic_population(ps)
    end
    return nothing
end


_pencil_es_phase(model) = _pencil_pic_phase_values(
    model.initialized[], model.primed[], model.last_dt[], model.time[], model.step[],
)

# Cached binding, not replaceable destination values: checkpoint restore may
# validate its immutable model contract while old particle scalars are stale.
function _pencil_es_bound_config(model::DistributedElectrostaticPIC{T}) where {T}
    model.neutrality in (:require, :neutralize) && model.grid_dc === :project ||
        throw(ArgumentError("distributed electrostatic source policies are invalid"))
    isfinite(model.background_charge) && isfinite(model.charge_rtol) &&
        HybridPlasmaPIC._stored_nonnegative(model.charge_rtol) ||
        throw(ArgumentError("distributed electrostatic background or charge tolerance is invalid"))
    length(model._species_scalars) == length(model.owners) ||
        throw(ArgumentError("distributed electrostatic cached species count is invalid"))
    for (q, m) in model._species_scalars
        HybridPlasmaPIC._validated_charge_mass(q, m, T)
    end
    return (
        T = T, n = model.g.n, L = model.g.L, dx = model.g.dx,
        background_charge = model.background_charge, neutrality = model.neutrality,
        charge_rtol = model.charge_rtol, grid_dc = model.grid_dc,
        species_scalars = Tuple(model._species_scalars), owned_ranges = Tuple(model._layout.by_rank),
    )
end

_pencil_es_config(model, species) = (_pencil_es_bound_config(model), Tuple((ps.q, ps.m) for ps in species))

function _pencil_es_agreement(model::DistributedElectrostaticPIC{T}, species, context; ready = true, factors = true) where {T}
    return _agree_configuration(model.g.fft.comm, context) do
        _pencil_es_mode(T)
        _pencil_es_storage(model)
        _pencil_es_species(model, species; bound = ready, factors = factors)
        if ready
            _pencil_es_phase(model)
            model.initialized[] || throw(ArgumentError("initialize distributed electrostatic PIC before use"))
            _finite_vector_magnitude(model.E)
            _finite_magnitude(model.rho)
        end
        (
            _pencil_es_config(model, species), ready ? (
                model.initialized[], model.primed[], model.last_dt[], model.time[], model.step[],
            ) : (),
        )
    end
end

"""
    distributed_electrostatic_pic(grid, species; grid_dc, background_charge=0,
                                 neutrality=:require, charge_rtol=64eps(T))

Collectively allocate a periodic 3D CPU Float32/64 rank-owned electrostatic PIC
collection. Each rank supplies the same ordered species (empty local populations
are valid), using independent Vector-backed ParticleSets with finite signed q,
positive finite m and nonnegative finite weights. Initialize with `init_espic!`.

`grid_dc=:project` is required explicitly: physical neutrality is checked against
original particle/background charge, whereas represented deposited-grid charge
is separately projected by the canonical Poisson zero-mode policy and reported.
CIC charge contributions retain original factors through exact rank routing and
one density rounding; this is a new map, not serial normal-density bit identity.
Gather and leapfrog reuse canonical native owners. Component arrays may be
replaced on migration while ParticleSet owners remain stable. No global field or
population is replicated. The caller owns its communicator/grid lifetime and must
exclude concurrent use. See `docs/DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md`.
"""
function HybridPlasmaPIC.distributed_electrostatic_pic(
    g::DistributedFourierGrid{T}, species; grid_dc = nothing, background_charge = zero(T),
    neutrality = :require, charge_rtol = 64eps(T),
) where {T}
    context = "Distributed electrostatic PIC"
    layout = HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "ownership construction") do
        _pencil_particle_layout(g)
    end
    parameters = _agree_configuration(g.fft.comm, context) do
        _pencil_es_mode(T)
        species isa Vector{<:HybridPlasmaPIC.ParticleSet{3,T}} ||
            throw(ArgumentError("species must be an ordinary Vector of matching 3D ParticleSets"))
        Base.require_one_based_indexing(species)
        all(ps -> ps isa _PencilParticle{T}, species) ||
            throw(ArgumentError("distributed electrostatic PIC requires Vector-backed particle storage"))
        foreach(_pencil_pic_population, species)
        neutrality in (:require, :neutralize) ||
            throw(ArgumentError("neutrality must be :require or :neutralize"))
        grid_dc === :project || throw(ArgumentError("grid_dc must explicitly be :project"))
        background_charge isa Real && charge_rtol isa Real ||
            throw(ArgumentError("background_charge and charge_rtol must be real"))
        background = HybridPlasmaPIC._require_finite_real("background_charge", background_charge, T)
        rtol = HybridPlasmaPIC._require_finite_nonnegative_real("charge_rtol", charge_rtol, T)
        T, g.n, g.L, g.dx, background, neutrality, rtol, grid_dc,
        Tuple((ps.q, ps.m) for ps in species), Tuple(layout.by_rank)
    end
    return HybridPlasmaPIC._mpi_collective_stage(g.fft.comm, context, "model allocation") do
        P = _PencilParticle{T}
        model = DistributedElectrostaticPIC(
            g, ntuple(_ -> HybridPlasmaPIC.distributed_field(g), 3), HybridPlasmaPIC.distributed_field(g),
            parameters[5], parameters[6], parameters[7], parameters[8], P[ps for ps in species],
            Tuple{T,T}[(ps.q, ps.m) for ps in species], Ref(false), Ref(false), Ref(zero(T)),
            Ref(zero(T)), Ref(0), ntuple(_ -> HybridPlasmaPIC.distributed_field(g), 3),
            HybridPlasmaPIC.distributed_field(g), P[], layout,
        )
        _pencil_es_storage(model)
        _pencil_es_species(model, species)
        model
    end
end


function _pencil_es_fields!(model, species, context)
    _pencil_pic_charges(model, species, context)
    _pencil_deposit_cic!(model._candidate_rho, species, model.g, model._layout, model.background_charge, 0, context)
    # Explicit stored policy: do not weaken the standalone strict :require API.
    HybridPlasmaPIC.distributed_poisson_field!(model._candidate_E, model._candidate_rho, model.g; neutrality = :neutralize)
    return nothing
end

function _pencil_es_stages!(model, species, context)
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "particle staging") do
        HybridPlasmaPIC._pic_prepare_species_stages!(model._staged, species)
        for ps in model._staged
            HybridPlasmaPIC.apply_periodic!(ps, ntuple(_ -> zero(eltype(ps)), 3), model.g.L)
        end
        model._staged
    end
end

function _pencil_es_commit!(model, species, prepared, scalars; reset = false, dt = model.last_dt[], clock = (model.time[], model.step[]))
    for (ps, replacement) in zip(species, prepared)
        HybridPlasmaPIC._publish_particle_replacement!(ps, replacement)
    end
    for c = 1:3
        copyto!(parent(model.E[c]), parent(model._candidate_E[c]))
    end
    copyto!(parent(model.rho), parent(model._candidate_rho))
    copyto!(model.owners, species)
    copyto!(model._species_scalars, scalars)
    model.initialized[] = true
    model.primed[] = reset ? false : model.primed[]
    model.last_dt[] = reset ? zero(dt) : dt
    model.time[], model.step[] = clock
    return model
end

function _pencil_es_final(model::DistributedElectrostaticPIC{T}, species, prepared, context) where {T}
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "publication preparation") do
        length(prepared) == length(species) || throw(ArgumentError("incoming species count changed"))
        foreach(_pencil_pic_population, prepared)
        _pencil_owned_particles(prepared, model.g, model._layout)
        _finite_vector_magnitude(model._candidate_E)
        _finite_magnitude(model._candidate_rho)
        # Typed replacements and bindings have their complete final lengths now;
        # no public array is resized or appended during publication.
        Tuple{T,T}[(ps.q, ps.m) for ps in prepared]
    end
end

"""
    init_espic!(model::DistributedElectrostaticPIC, species)

Collectively wrap/redistribute supplied full-time particles, form the common
field and publish together. Initialization resets time/step and leapfrog phase;
call `poisson_E!(model)` to refresh sources without resetting an existing phase.
"""
function HybridPlasmaPIC.init_espic!(model::DistributedElectrostaticPIC{T}, species) where {T}
    context = "Distributed electrostatic initialization"
    _pencil_es_agreement(model, species, context; ready = false)
    staged = _pencil_es_stages!(model, species, context)
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    _pencil_es_fields!(model, prepared, context)
    scalars = _pencil_es_final(model, species, prepared, context)
    return _pencil_es_commit!(model, species, prepared, scalars; reset = true, dt = zero(T), clock = (zero(T), 0))
end

"Collectively refresh bound particles/density/E, including redistribution, without changing the velocity half level or clock."
function HybridPlasmaPIC.poisson_E!(model::DistributedElectrostaticPIC)
    context = "Distributed electrostatic refresh"
    species = model.owners
    _pencil_es_agreement(model, species, context; factors = false)
    staged = _pencil_es_stages!(model, species, context)
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    _pencil_es_fields!(model, prepared, context)
    scalars = _pencil_es_final(model, species, prepared, context)
    _pencil_es_commit!(model, species, prepared, scalars)
    return model.E
end

# Includes canonical private first-step priming, but no public phase store.
function _pencil_es_advance_candidate!(model, species, h, context)
    staged = _pencil_es_stages!(model, species, context)
    gathers = _pencil_gather(model.E, staged, model.g, model._layout, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "canonical particle push") do
        for (ps, field) in zip(staged, gathers)
            HybridPlasmaPIC._espic_push!(ps, field, model.g, h, model.primed[], model.last_dt[])
        end
    end
    prepared = _pencil_migrate(staged, model.g, model._layout, context)
    _pencil_es_fields!(model, prepared, context)
    return prepared
end

"""
    step_espic!(model::DistributedElectrostaticPIC, species, dt)

Collectively gather old E, apply the shared staggered mover, route complete
particles to their new owners and solve the new charge/E. All public particles,
fields and phase/clocks commit after final collective success. Recoverable local
failure preserves them; scratch may change. MPI loss/fatal resource failure is
excluded. Zero converted dt validates then leaves the complete state unchanged.
"""
function HybridPlasmaPIC.step_espic!(model::DistributedElectrostaticPIC{T}, species, dt) where {T}
    context = "Distributed electrostatic step"
    _pencil_es_agreement(model, species, context)
    h, next_time, next_step = _agree_configuration(model.g.fft.comm, context) do
        dt isa Real || throw(ArgumentError("distributed electrostatic dt must be real"))
        h = HybridPlasmaPIC._require_resolved_dt(context, HybridPlasmaPIC._validated_nonnegative_dt(T, dt; name = context))
        if HybridPlasmaPIC._particle_stored_zero(h)
            h, model.time[], model.step[]
        else
            next_time, next_step = HybridPlasmaPIC._validated_step_clock(model.time[], model.step[], h)
            next_time > model.time[] || throw(ArgumentError("distributed electrostatic time must advance"))
            h, next_time, next_step
        end
    end
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "current particle ownership") do
        _pencil_owned_particles(species, model.g, model._layout)
    end
    HybridPlasmaPIC._particle_stored_zero(h) && return model
    prepared = _pencil_es_advance_candidate!(model, species, h, context)
    scalars = _pencil_es_final(model, species, prepared, context)
    _pencil_es_commit!(model, species, prepared, scalars; dt = h, clock = (next_time, next_step))
    model.primed[] = true
    return model
end

"""
    electrostatic_diagnostics(model::DistributedElectrostaticPIC, species)

Collective exact physical/grid charge ledgers, global particle count, momentum,
kinetic and field energy, time/step and stored half-velocity offset. Charges/volumes
remain rational; global energies and momentum components round once after exact
scalar reduction. Momentum uses the stored carried velocity half level. No field
or population is gathered. The ledgers separate physical-policy and represented-grid
compensators; these are descriptions of the same projection, not two additions.
"""
function HybridPlasmaPIC.electrostatic_diagnostics(model::DistributedElectrostaticPIC{T}, species) where {T}
    context = "Distributed electrostatic diagnostics"
    _pencil_es_agreement(model, species, context)
    HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "particle ownership") do
        _pencil_owned_particles(species, model.g, model._layout)
    end
    charges = _pencil_pic_charges(model, species, context)
    grid_sum = _pencil_exact_scalar_sum(model.g.fft.comm, context, "grid charge preparation", "grid charge accumulation") do
        HybridPlasmaPIC._exact_float_sum(parent(model.rho))
    end
    kinetic = _pencil_exact_scalar_sum(model.g.fft.comm, context, "kinetic preparation", "kinetic accumulation") do
        HybridPlasmaPIC._kinetic_energy_rational(species)
    end
    local_momentum = HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "momentum preparation") do
        HybridPlasmaPIC._momentum_rational(species)
    end
    momentum = ntuple(3) do c
        _pencil_exact_scalar_sum(model.g.fft.comm, context, "momentum component serialization", "momentum component accumulation") do
            local_momentum[c]
        end
    end
    field = _pencil_exact_scalar_sum(model.g.fft.comm, context, "field energy preparation", "field energy accumulation") do
        HybridPlasmaPIC._field_energy_rational(map(parent, model.E), (), zero(_PencilCharge), model.g.dx)
    end
    count = _pencil_exact_scalar_sum(model.g.fft.comm, context, "particle count preparation", "particle count accumulation") do
        sum((BigInt(length(ps.weight)) for ps in species); init = BigInt(0)) // BigInt(1)
    end
    return HybridPlasmaPIC._mpi_collective_stage(model.g.fft.comm, context, "diagnostic publication") do
        Q = HybridPlasmaPIC._exact_rational
        cell_volume = prod(Q, model.g.dx)
        declared_volume, grid_volume = prod(Q, model.g.L), prod(BigInt, model.g.n) * cell_volume
        grid_charge = grid_sum * cell_volume
        physical_compensator = model.neutrality === :neutralize ? -charges.net : zero(_PencilCharge)
        (
            time = model.time[], step = model.step[], particle_count = numerator(count),
            particle_charge = charges.particle, background_charge = charges.background,
            physical_neutrality_residual = charges.net,
            physical_policy_compensator = physical_compensator,
            deposited_grid_charge = grid_charge, grid_dc_compensator = -grid_charge,
            declared_volume = declared_volume, stored_grid_volume = grid_volume,
            background_volume_discrepancy = Q(model.background_charge) * (grid_volume - declared_volume),
            deposition_discrepancy = grid_charge - Q(model.background_charge) * grid_volume - charges.particle,
            kinetic = HybridPlasmaPIC._from_exact(T, kinetic), field = HybridPlasmaPIC._from_exact(T, field),
            total = HybridPlasmaPIC._from_exact(T, kinetic + field),
            momentum = map(value -> HybridPlasmaPIC._from_exact(T, value), momentum),
            velocity_time_offset = model.primed[] ? -model.last_dt[] / 2 : zero(T),
        )
    end
end
