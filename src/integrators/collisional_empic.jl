# Field-first fixed-step composition on the canonical EMPIC collection candidate.

"""
    CollisionalEMPICStepper(inner::EMPIC, species::AbstractVector;
        dt, rng, coulomb_log, epsilon0, weight_to_number, u_floor=1e-3)

First-order, fixed-step CPU composition of synchronous electromagnetic PIC and
cell-local nonrelativistic [`collide_coulomb!`](@ref). Supports Float32/64 in
D=1:3 with ordinary Vector particle leaves, all populations mobile, `n_sub=1`,
zero background charge, `neutrality=:require` and no species schedule owner,
even for opted-in all-fast `n_sub=1`. `weight_to_number` explicitly
converts macro-weights to real-particle counts. The collision `epsilon0` and
`weight_to_number` must resolve to the same positive value: with these stored
weights, the Maxwell charge/current coefficient is `weight_to_number/epsilon0`.
Nearest rounding and gradual
underflow are required; the operation never changes the caller's arithmetic mode.

Call `init!` with physical full-time velocities before stepping. Each nonzero
`step!(state, species, dt)` first forms the canonical EMPIC candidate, then
scatters its carried half-time velocities at the new integer-time positions.
Positions and charge snapshots do not change during scattering. The current
remains the actual preceding drift source; neither Maxwell nor current is
recomputed after scattering. This is a field-then-collision Lie split, not a
symmetric split or an exact discrete energy law. Unequal macro-weights inherit
the collision owner's momentum/energy conservation in expectation.

Particle objects/order/counts, q/m, grid, field/collision configuration and dt
are fixed. Weights and provenance are unchanged by evolution. External edits
retain the underlying EMPIC phase and consistent-source responsibilities. A
zero dt validates the completed state and consumes no priming or RNG draws;
all-zero Coulomb logarithms skip collision and all RNG copying. Relativistic,
fixed-population, source-injection and independent-species schedules are outside
this first scope.

The inner state is borrowed exclusively. Fields, particles, clocks and RNG stage
until the whole operation succeeds; expected input/arithmetic/sampling errors
publish none of them. Private candidates, FFT buffers and collision scratch may
change on rejection and are reset before reuse. Ordinary writes exclude resource
exhaustion, OS-protected memory, custom setters and concurrent external mutation.
The explicit Xoshiro or MersenneTwister RNG is deeply detached from its caller.
Initialization preserves its current stream. Checkpoints/archives retain full
advanced RNG state and require the same Julia runtime/platform.
"""
mutable struct CollisionalEMPICStepper{D,T,S,R,L,SH,G<:FourierGrid{D,T}}
    const inner::S
    const grid::G
    const species::Vector{_ElectrostaticHostParticle{D,T}}
    const identities::Vector{Tuple{Int,T,T}}
    const c::T
    const shape::SH
    const charge_rtol::T
    const dt::T
    const coulomb_log::L
    const log_binding::L
    const epsilon0::T
    const weight_to_number::T
    const u_floor::T
    const rng::R
    const candidate_rng::R
    const workspace::CoulombCollisionWorkspace{T}
    const busy::Threads.Atomic{Bool}
    phase::Symbol
    source_time::T
    source_count::Int
end

function _validate_collisional_empic_mode()
    _boris_accuracy_mode() || throw(
        ArgumentError("CollisionalEMPICStepper requires nearest rounding and gradual underflow"),
    )
    return nothing
end

function _validate_collisional_empic_model(inner::EMPIC{D,T}) where {D,T}
    T in (Float32, Float64) ||
        throw(ArgumentError("CollisionalEMPICStepper requires Float32 or Float64"))
    work = _empic_collection(inner)
    work.schedule === nothing || throw(
        ArgumentError("CollisionalEMPICStepper requires an unscheduled generic EMPIC collection"),
    )
    !inner.relativistic ||
        throw(ArgumentError("CollisionalEMPICStepper requires nonrelativistic EMPIC"))
    inner.n_sub == 1 || throw(ArgumentError("CollisionalEMPICStepper requires n_sub = 1"))
    _particle_stored_zero(inner.n0) ||
        throw(ArgumentError("CollisionalEMPICStepper requires zero background charge"))
    work.neutrality === :require ||
        throw(ArgumentError("CollisionalEMPICStepper requires neutrality = :require"))
    all(work.mobile) ||
        throw(ArgumentError("CollisionalEMPICStepper requires every population mobile"))
    inner.shape isa Union{NGP,CIC,TSC} ||
        throw(ArgumentError("CollisionalEMPICStepper requires NGP, CIC or TSC"))
    _require_finite_positive_real("c", inner.c, T)
    return work
end

function _validate_collisional_empic_population_storage(species, grid::FourierGrid{D,T}) where {D,T}
    _espic_validate_collection_type(species, grid)
    for ps in species
        ps isa _ElectrostaticHostParticle{D,T} || throw(
            ArgumentError("CollisionalEMPICStepper requires ordinary Vector particle storage"),
        )
    end
    return nothing
end

function CollisionalEMPICStepper(
    inner::EMPIC{D,T},
    species::AbstractVector;
    dt::Real,
    rng,
    coulomb_log,
    epsilon0::Real,
    weight_to_number::Real,
    u_floor::Real = 1e-3,
) where {D,T}
    _validate_collisional_empic_mode()
    work = _validate_collisional_empic_model(inner)
    _validate_collisional_empic_population_storage(species, inner.g)
    _empic_collection_validate_species(inner, species)
    _empic_collection_validate_phase(inner)
    _particle_stored_zero(inner.time[]) &&
        inner.step[] == 0 &&
        _particle_stored_zero(inner.last_dt[]) &&
        _particle_stored_zero(inner.last_sub_dt[]) ||
        throw(ArgumentError("construct collisional EMPIC from an unadvanced inner state"))
    h = _require_resolved_dt("collisional EMPIC dt", _require_finite_positive_real("dt", dt, T))
    epsilon = _require_finite_positive_real("epsilon0", epsilon0, T)
    number_scale = _require_finite_positive_real("weight_to_number", weight_to_number, T)
    isequal(epsilon, number_scale) || throw(
        ArgumentError("collisional EMPIC requires equal resolved epsilon0 and weight_to_number"),
    )
    uf = _require_finite_positive_real("u_floor", u_floor, T)
    _validate_coulomb_logs(coulomb_log, length(species), T)
    logs = coulomb_log isa Real ? T(coulomb_log) : Matrix{T}(coulomb_log)
    rng isa Union{Random.Xoshiro,Random.MersenneTwister} ||
        throw(ArgumentError("CollisionalEMPICStepper requires Xoshiro or MersenneTwister"))
    _validate_collisional_rng(rng)
    _validate_checkpoint_configuration((species = species,), inner, (rng,))
    populations = Vector{_ElectrostaticHostParticle{D,T}}(species)
    owned_rng = deepcopy(rng)
    st = CollisionalEMPICStepper(
        inner,
        inner.g,
        populations,
        Tuple{Int,T,T}[(nparticles(ps), ps.q, ps.m) for ps in populations],
        inner.c,
        inner.shape,
        work.charge_rtol,
        h,
        logs,
        copy(logs),
        epsilon,
        number_scale,
        uf,
        owned_rng,
        deepcopy(owned_rng),
        CoulombCollisionWorkspace(populations, inner.g),
        Threads.Atomic{Bool}(false),
        :uninitialized,
        zero(T),
        0,
    )
    _validate_collisional_empic_binding(st, populations)
    return st
end

function _validate_collisional_empic_rngs(st::CollisionalEMPICStepper)
    _validate_collisional_rng_pair(st.rng, st.candidate_rng, "collisional EMPIC")
    sources = (
        species = st.species,
        identities = st.identities,
        workspace = st.workspace,
        logs = (st.coulomb_log, st.log_binding),
    )
    _validate_checkpoint_configuration(sources, st.inner, (st.rng,))
    _validate_checkpoint_configuration(sources, st.inner, (st.candidate_rng,))
    return nothing
end

function _validate_collisional_empic_binding(
    st::CollisionalEMPICStepper{D,T},
    species;
    values::Bool = true,
) where {D,T}
    _validate_collisional_empic_mode()
    st.phase in (:uninitialized, :ready) ||
        throw(ArgumentError("collisional EMPIC state is in flight or has an invalid phase"))
    _validate_collisional_empic_population_storage(species, st.grid)
    length(species) == length(st.species) ||
        throw(ArgumentError("collisional EMPIC species count is fixed"))
    for a in eachindex(st.species)
        ps = species[a]
        ps === st.species[a] ||
            throw(ArgumentError("collisional EMPIC species identities and order are fixed"))
        nparticles(ps) == st.identities[a][1] ||
            throw(ArgumentError("collisional EMPIC particle counts are fixed"))
        !values ||
            isequal((nparticles(ps), ps.q, ps.m), st.identities[a]) ||
            throw(ArgumentError("collisional EMPIC charges and masses are fixed"))
    end
    st.inner.g === st.grid || throw(ArgumentError("collisional EMPIC grid is fixed"))
    work = _validate_collisional_empic_model(st.inner)
    isequal(st.inner.c, st.c) &&
        isequal(st.inner.shape, st.shape) &&
        isequal(work.charge_rtol, st.charge_rtol) ||
        throw(ArgumentError("collisional EMPIC field configuration is fixed"))
    isequal(st.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional EMPIC Coulomb logarithms are fixed"))
    _require_resolved_dt("collisional EMPIC dt", _require_finite_positive_real("dt", st.dt, T))
    _require_finite_positive_real("epsilon0", st.epsilon0, T)
    _require_finite_positive_real("weight_to_number", st.weight_to_number, T)
    _require_finite_positive_real("u_floor", st.u_floor, T)
    isequal(st.epsilon0, st.weight_to_number) || throw(
        ArgumentError("collisional EMPIC requires equal resolved epsilon0 and weight_to_number"),
    )
    _validate_coulomb_logs(st.coulomb_log, length(st.species), T)
    charges = if values
        result = _empic_collection_validate_species(st.inner, st.species)
        _validate_local_coulomb_species(st.species, st.grid)
        result
    else
        _empic_collection_validate_storage(st.inner, st.species)
        nothing
    end
    # Mobility is readonly configuration, including on Julia versions whose
    # BitArray dataids do not expose the packed backing array.
    for ps in st.species
        _mightalias_particle_storage(work.mobile.chunks, ps) &&
            throw(ArgumentError("collisional EMPIC particles must not alias mobility storage"))
    end
    if st.coulomb_log isa AbstractArray
        _particle_arrays_mightalias(st.coulomb_log, st.log_binding) &&
            throw(ArgumentError("collisional EMPIC fixed logarithm storage must be independent"))
    end
    _validate_checkpoint_configuration(
        (
            species = st.species,
            identities = st.identities,
            workspace = st.workspace,
            rngs = (st.rng, st.candidate_rng),
        ),
        st.inner,
        (st.coulomb_log, st.log_binding),
    )
    _validate_coulomb_workspace(st.workspace, st.species, st.grid, st.coulomb_log)
    _validate_checkpoint_configuration(
        (
            species = st.species,
            identities = st.identities,
            logs = (st.coulomb_log, st.log_binding),
            rngs = (st.rng, st.candidate_rng),
        ),
        st.inner,
        (st.workspace,),
    )
    _validate_collisional_empic_rngs(st)
    return charges
end

function _validate_collisional_empic_boundary(st::CollisionalEMPICStepper)
    st.phase === :ready || throw(ArgumentError("initialize collisional EMPIC before use"))
    inner = st.inner
    work = _empic_collection(inner)
    work.initialized[] || throw(ArgumentError("collisional EMPIC inner state must be initialized"))
    _empic_collection_bound(work, st.species)
    _empic_collection_validate_phase(inner)
    isequal(st.source_time, inner.time[]) && st.source_count == inner.step[] ||
        throw(ArgumentError("collisional EMPIC source ledger does not match its field clock"))
    initial = st.source_count == 0
    (initial ? _particle_stored_zero(st.source_time) : true) &&
        (initial ? _particle_stored_zero(inner.last_dt[]) : isequal(inner.last_dt[], st.dt)) &&
        (initial ? _particle_stored_zero(inner.last_sub_dt[]) : isequal(inner.last_sub_dt[], st.dt)) ||
        throw(ArgumentError("collisional EMPIC history does not match its fixed dt"))
    return nothing
end

function _acquire_collisional_empic!(st::CollisionalEMPICStepper)
    Threads.atomic_cas!(st.busy, false, true) &&
        throw(ArgumentError("concurrent or reentrant collisional EMPIC operation"))
    return nothing
end

function init!(st::CollisionalEMPICStepper, species::AbstractVector)
    _acquire_collisional_empic!(st)
    try
        _validate_collisional_empic_binding(st, species)
        # The inner initializer is already a complete staged publication. Only
        # typed assignments follow its first public store; RNG is untouched.
        init_empic!(st.inner, st.species)
        st.source_time = st.inner.time[]
        st.source_count = st.inner.step[]
        st.phase = :ready
        return st
    finally
        st.busy[] = false
    end
end

function step!(st::CollisionalEMPICStepper{D,T}, species::AbstractVector, dt::Real) where {D,T}
    h = _require_resolved_dt(
        "collisional EMPIC step!",
        _validated_nonnegative_dt(T, dt; name = "collisional EMPIC step!"),
    )
    _particle_stored_zero(h) &&
        !iszero(dt) &&
        throw(ArgumentError("nonzero collisional EMPIC dt must not underflow at grid precision"))
    _acquire_collisional_empic!(st)
    try
        charges = _validate_collisional_empic_binding(st, species)
        _validate_collisional_empic_boundary(st)
        _empic_collection_require_finite_fields(st.inner)
        _particle_stored_zero(h) && return st
        isequal(h, st.dt) || throw(ArgumentError("collisional EMPIC step! requires its fixed dt"))
        next_time, next_step = _validated_step_clock(st.inner.time[], st.inner.step[], h, 1)
        active = _collisional_active(st.coulomb_log)
        active && copy!(st.candidate_rng, st.rng)
        previous_phase = st.phase
        st.phase = :inflight
        completed = false
        try
            work = _empic_collection(st.inner)
            _empic_collection_advance_candidate!(
                st.inner,
                st.species,
                work,
                charges,
                h,
                h,
                next_time,
                next_step,
            )
            if active
                collide_coulomb!(
                    work.staged,
                    st.grid,
                    h;
                    coulomb_log = st.coulomb_log,
                    epsilon0 = st.epsilon0,
                    weight_to_number = st.weight_to_number,
                    u_floor = st.u_floor,
                    rng = st.candidate_rng,
                    workspace = st.workspace,
                )
                for ps in work.staged
                    _require_finite_particle_velocities(ps, "collisional EMPIC output")
                end
                _validate_collisional_rng(st.candidate_rng)
            end
            # Physical candidate validation is already complete: Coulomb only
            # changes the staged velocities checked above. The normalized stock
            # RNG copy and ordinary fixed-size publication cannot reject valid
            # preflighted storage; no validation/resize/arithmetic follows.
            active && copy!(st.rng, st.candidate_rng)
            _empic_collection_publish!(st.inner, work; particles = true)
            st.source_time = next_time
            st.source_count = next_step
            completed = true
            return st
        finally
            st.phase = completed ? :ready : previous_phase
        end
    finally
        st.busy[] = false
    end
end

function electromagnetic_diagnostics(st::CollisionalEMPICStepper, species::AbstractVector)
    _acquire_collisional_empic!(st)
    try
        _validate_collisional_empic_binding(st, species)
        _validate_collisional_empic_boundary(st)
        return merge(
            electromagnetic_diagnostics(st.inner, st.species),
            (
                split_policy = :field_then_coulomb,
                source_time = st.source_time,
                source_count = st.source_count,
            ),
        )
    finally
        st.busy[] = false
    end
end
