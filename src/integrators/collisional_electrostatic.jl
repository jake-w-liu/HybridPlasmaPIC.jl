# Field-first fixed-step composition on the canonical ES collection candidate.

"""
    CollisionalElectrostaticPIC(inner::ElectrostaticPIC, species::AbstractVector;
        dt, rng, coulomb_log, epsilon0, weight_to_number, u_floor=1e-3)

Fixed-step CPU composition of periodic electrostatic PIC and cell-local
nonrelativistic [`collide_coulomb!`](@ref), for Float32/64 in D=1:3. Uses the
generic all-mobile CIC collection, ordinary Vector particle leaves, zero
background charge and `neutrality=:require`. Both normalization keywords are
required and must resolve to the same positive value: the Poisson source
coefficient is `weight_to_number/epsilon0` in these stored particle units.
Nearest rounding and gradual underflow are required and never changed internally.

Call `init!` with physical full-time velocities. Each positive step forms the
canonical ES kick/drift/Poisson candidate, then scatters its carried half-time
velocities at the new integer-time positions. Scattering changes neither charge
nor positions, so E/ne are retained without another Poisson solve. This is a
field-then-collision first-order Lie split, not a symmetric split or an exact
discrete energy law. Unequal macro-weights inherit conservation in expectation.

Particle identities/order/counts, q/m, field/collision configuration and dt are
fixed. Evolution leaves weights/provenance unchanged. Zero dt validates without
priming/draws; all-zero logs skip collisions and both RNG copies. The wrapper's
`source_time` and `source_count` are its sole clock and completed-source ledger;
the inner ES state has no clock. It is borrowed exclusively: an extra same-dt
bare-inner advance after priming cannot be detected from its phase metadata.
External edits retain consistent-field and declared velocity-phase obligations.

Fields, particles, phase, clock and advanced RNG publish only after successful
private staging. Expected input/arithmetic/sampling failures preserve public
state; private FFT/particle/collision/RNG scratch may change and resets on reuse.
Ordinary publication excludes resource exhaustion, protected memory, custom
setters and concurrent external mutation. Use independent states/grids/populations
for concurrent calls. The explicit Xoshiro or MersenneTwister is deeply detached;
init preserves its stream, and restart retains its full same-runtime state.
"""
mutable struct CollisionalElectrostaticPIC{D,T,S,R,L,G<:FourierGrid{D,T}}
    const inner::S
    const grid::G
    const species::Vector{_ElectrostaticHostParticle{D,T}}
    const identities::Vector{Tuple{Int,T,T}}
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

function _validate_collisional_espic_mode()
    _boris_accuracy_mode() || throw(
        ArgumentError("CollisionalElectrostaticPIC requires nearest rounding and gradual underflow"),
    )
    return nothing
end

function _validate_collisional_espic_model(inner::ElectrostaticPIC{D,T}) where {D,T}
    T in (Float32, Float64) ||
        throw(ArgumentError("CollisionalElectrostaticPIC requires Float32 or Float64"))
    work = _espic_collection(inner)
    _particle_stored_zero(inner.n0) ||
        throw(ArgumentError("CollisionalElectrostaticPIC requires zero background charge"))
    work.neutrality === :require ||
        throw(ArgumentError("CollisionalElectrostaticPIC requires neutrality = :require"))
    return work
end

function _validate_collisional_espic_population_storage(species, grid::FourierGrid{D,T}) where {D,T}
    _espic_validate_collection_type(species, grid)
    for ps in species
        ps isa _ElectrostaticHostParticle{D,T} ||
            throw(ArgumentError("CollisionalElectrostaticPIC requires ordinary Vector particle storage"))
    end
    return nothing
end

function CollisionalElectrostaticPIC(
    inner::ElectrostaticPIC{D,T},
    species::AbstractVector;
    dt::Real,
    rng,
    coulomb_log,
    epsilon0::Real,
    weight_to_number::Real,
    u_floor::Real = 1e-3,
) where {D,T}
    _validate_collisional_espic_mode()
    work = _validate_collisional_espic_model(inner)
    _validate_collisional_espic_population_storage(species, inner.g)
    _espic_validate_species(inner, species)
    !inner.primed[] && _particle_stored_zero(inner.last_dt[]) ||
        throw(ArgumentError("construct collisional ES from an unprimed inner state"))
    h = _require_resolved_dt("collisional ES dt", _require_finite_positive_real("dt", dt, T))
    epsilon = _require_finite_positive_real("epsilon0", epsilon0, T)
    number_scale = _require_finite_positive_real("weight_to_number", weight_to_number, T)
    isequal(epsilon, number_scale) || throw(
        ArgumentError("collisional ES requires equal resolved epsilon0 and weight_to_number"),
    )
    uf = _require_finite_positive_real("u_floor", u_floor, T)
    _validate_coulomb_logs(coulomb_log, length(species), T)
    logs = coulomb_log isa Real ? T(coulomb_log) : Matrix{T}(coulomb_log)
    rng isa Union{Random.Xoshiro,Random.MersenneTwister} ||
        throw(ArgumentError("CollisionalElectrostaticPIC requires Xoshiro or MersenneTwister"))
    _validate_collisional_rng(rng)
    _validate_checkpoint_configuration((species = species,), inner, (rng,))
    populations = Vector{_ElectrostaticHostParticle{D,T}}(species)
    owned_rng = deepcopy(rng)
    st = CollisionalElectrostaticPIC(
        inner,
        inner.g,
        populations,
        Tuple{Int,T,T}[(nparticles(ps), ps.q, ps.m) for ps in populations],
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
    _validate_collisional_espic_binding(st, populations)
    return st
end

function _require_collisional_espic_metadata_independent(array, g)
    for metadata in (g.midx..., g.kfull..., g.ik..., g.kvec...)
        _particle_arrays_mightalias(array, metadata) &&
            throw(ArgumentError("collisional ES writes must not alias readonly grid metadata"))
    end
    return nothing
end

function _validate_collisional_espic_role_storage(st::CollisionalElectrostaticPIC)
    inner, g = st.inner, st.grid
    work = _espic_collection(inner)
    writers = (
        inner.E..., inner.ne, work.E..., work.ne, work.density, inner.Ep...,
        g.cbuf, g.tbuf, g.abuf, g.sbuf,
    )
    for a in eachindex(writers)
        _require_collisional_espic_metadata_independent(writers[a], g)
        for b = 1:(a-1)
            _particle_arrays_mightalias(writers[a], writers[b]) &&
                throw(ArgumentError("collisional ES field, gather and FFT storage must be independent"))
        end
    end
    for (a, ps) in enumerate(work.staged)
        _particle_components_mightalias(ps) &&
            throw(ArgumentError("collisional ES candidate components must be independent"))
        for b = 1:(a-1)
            _particle_sets_mightalias(ps, work.staged[b]) &&
                throw(ArgumentError("collisional ES candidate populations must be independent"))
        end
        for array in writers
            _mightalias_particle_storage(array, ps) &&
                throw(ArgumentError("collisional ES candidates must not alias solver storage"))
        end
    end
    # Particles, collision scratch and RNG caches also write during step/restore.
    for populations in (st.species, work.staged), ps in populations
        for array in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            _require_collisional_espic_metadata_independent(array, g)
        end
    end
    for population in st.workspace.populations
        for array in (population.velocity..., _coulomb_index_arrays(population)...)
            _require_collisional_espic_metadata_independent(array, g)
        end
    end
    for rng in (st.rng, st.candidate_rng), array in _collisional_rng_buffers(rng)
        _require_collisional_espic_metadata_independent(array, g)
    end
    return nothing
end

function _validate_collisional_espic_rngs(st::CollisionalElectrostaticPIC)
    _validate_collisional_rng_pair(st.rng, st.candidate_rng, "collisional ES")
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

function _validate_collisional_espic_binding(
    st::CollisionalElectrostaticPIC{D,T}, species; values::Bool = true,
) where {D,T}
    _validate_collisional_espic_mode()
    st.phase in (:uninitialized, :ready) ||
        throw(ArgumentError("collisional ES state is in flight or has an invalid phase"))
    _validate_collisional_espic_population_storage(species, st.grid)
    length(species) == length(st.species) ||
        throw(ArgumentError("collisional ES species count is fixed"))
    for a in eachindex(st.species)
        ps = species[a]
        ps === st.species[a] ||
            throw(ArgumentError("collisional ES species identities and order are fixed"))
        nparticles(ps) == st.identities[a][1] ||
            throw(ArgumentError("collisional ES particle counts are fixed"))
        !values || isequal((nparticles(ps), ps.q, ps.m), st.identities[a]) ||
            throw(ArgumentError("collisional ES charges and masses are fixed"))
    end
    st.inner.g === st.grid || throw(ArgumentError("collisional ES grid is fixed"))
    work = _validate_collisional_espic_model(st.inner)
    isequal(work.charge_rtol, st.charge_rtol) ||
        throw(ArgumentError("collisional ES field configuration is fixed"))
    isequal(st.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional ES Coulomb logarithms are fixed"))
    _require_resolved_dt("collisional ES dt", _require_finite_positive_real("dt", st.dt, T))
    _require_finite_positive_real("epsilon0", st.epsilon0, T)
    _require_finite_positive_real("weight_to_number", st.weight_to_number, T)
    _require_finite_positive_real("u_floor", st.u_floor, T)
    isequal(st.epsilon0, st.weight_to_number) || throw(
        ArgumentError("collisional ES requires equal resolved epsilon0 and weight_to_number"),
    )
    _validate_coulomb_logs(st.coulomb_log, length(st.species), T)
    if values
        _espic_validate_species(st.inner, st.species)
        _validate_local_coulomb_species(st.species, st.grid)
    else
        _espic_validate_storage(st.inner, st.species)
    end
    _espic_validate_owner_types(work, st.species)
    if st.coulomb_log isa AbstractArray
        _particle_arrays_mightalias(st.coulomb_log, st.log_binding) &&
            throw(ArgumentError("collisional ES fixed logarithm storage must be independent"))
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
    _validate_collisional_espic_rngs(st)
    _validate_collisional_espic_role_storage(st)
    return nothing
end

function _validate_collisional_espic_boundary(st::CollisionalElectrostaticPIC{D,T}) where {D,T}
    st.phase === :ready || throw(ArgumentError("initialize collisional ES before use"))
    inner = st.inner
    _espic_collection(inner).initialized[] ||
        throw(ArgumentError("collisional ES inner state must be initialized"))
    _espic_require_bound_species(_espic_collection(inner), st.species)
    _validate_checkpoint_clock_state((time = st.source_time, step = st.source_count), T)
    (
        st.source_count == 0 ?
        _particle_stored_zero(st.source_time) && !inner.primed[] && _particle_stored_zero(inner.last_dt[]) :
        st.source_time > zero(T) && inner.primed[] && isequal(inner.last_dt[], st.dt)
    ) || throw(ArgumentError("collisional ES clock and carried phase do not match its fixed dt"))
    return nothing
end

function _acquire_collisional_espic!(st::CollisionalElectrostaticPIC)
    Threads.atomic_cas!(st.busy, false, true) &&
        throw(ArgumentError("concurrent or reentrant collisional ES operation"))
    return nothing
end

function init!(st::CollisionalElectrostaticPIC{D,T}, species::AbstractVector) where {D,T}
    _acquire_collisional_espic!(st)
    try
        _validate_collisional_espic_binding(st, species)
        # Canonical init publishes fields before rebinding owners. Reserve all
        # reference capacity first; only compatible ordinary writes then follow.
        sizehint!(_espic_collection(st.inner).owners, length(st.species); shrink = false)
        init_espic!(st.inner, st.species)
        st.source_time = zero(T)
        st.source_count = 0
        st.phase = :ready
        return st
    finally
        st.busy[] = false
    end
end

function step!(st::CollisionalElectrostaticPIC{D,T}, species::AbstractVector, dt::Real) where {D,T}
    h = _require_resolved_dt(
        "collisional ES step!", _validated_nonnegative_dt(T, dt; name = "collisional ES step!"),
    )
    _particle_stored_zero(h) && !iszero(dt) &&
        throw(ArgumentError("nonzero collisional ES dt must not underflow at grid precision"))
    _acquire_collisional_espic!(st)
    try
        _validate_collisional_espic_binding(st, species)
        _validate_collisional_espic_boundary(st)
        all(a -> all(isfinite, a), st.inner.E) ||
            throw(ArgumentError("collisional ES field must be finite"))
        _particle_stored_zero(h) && return st
        isequal(h, st.dt) || throw(ArgumentError("collisional ES step! requires its fixed dt"))
        next_time, next_count = _validated_step_clock(st.source_time, st.source_count, h)
        next_time > st.source_time ||
            throw(ArgumentError("collisional ES clock must advance at the stored precision"))
        active = _collisional_active(st.coulomb_log)
        active && copy!(st.candidate_rng, st.rng)
        previous_phase = st.phase
        st.phase = :inflight
        completed = false
        try
            work = _espic_collection(st.inner)
            _espic_collection_advance_candidate!(st.inner, st.species, work, h)
            if active
                collide_coulomb!(
                    work.staged, st.grid, h;
                    coulomb_log = st.coulomb_log,
                    epsilon0 = st.epsilon0,
                    weight_to_number = st.weight_to_number,
                    u_floor = st.u_floor,
                    rng = st.candidate_rng,
                    workspace = st.workspace,
                )
                for ps in work.staged
                    _require_finite_particle_velocities(ps, "collisional ES output")
                end
                _validate_collisional_rng(st.candidate_rng)
            end
            # All candidate arithmetic, dimensions and RNG structure passed.
            # Fixed-size stock copying and typed assignments cannot discover a
            # new expected input error after the first public RNG/physical store.
            active && copy!(st.rng, st.candidate_rng)
            _espic_collection_publish!(st.inner, st.species, work, h)
            st.source_time = next_time
            st.source_count = next_count
            completed = true
            return st
        finally
            st.phase = completed ? :ready : previous_phase
        end
    finally
        st.busy[] = false
    end
end

function electrostatic_diagnostics(st::CollisionalElectrostaticPIC, species::AbstractVector)
    _acquire_collisional_espic!(st)
    try
        _validate_collisional_espic_binding(st, species)
        _validate_collisional_espic_boundary(st)
        return merge(
            electrostatic_diagnostics(st.inner, st.species),
            (split_policy = :field_then_coulomb, source_time = st.source_time, source_count = st.source_count),
        )
    finally
        st.busy[] = false
    end
end
