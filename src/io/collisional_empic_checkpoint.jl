# Completed stochastic EMPIC boundaries extend the existing checkpoint owners.
# The nested generic collection state is the only physical particle/field payload.
const _COLLISIONAL_EMPIC_CHECKPOINT_FIELDS = (
    :integrator,
    :schema,
    :split_policy,
    :collision_position_policy,
    :collision_velocity_policy,
    :inner,
    :dt,
    :coulomb_log,
    :epsilon0,
    :weight_to_number,
    :u_floor,
    :rng,
    :rng_algorithm,
    :runtime,
    :phase,
    :initialized,
    :source_time,
    :source_count,
)

function _checkpoint_state(st::CollisionalEMPICStepper, species::AbstractVector)
    _validate_collisional_empic_binding(st, species)
    _validate_collisional_empic_boundary(st)
    _validate_collisional_empic_rngs(st)
    return (
        integrator = :collisional_empic,
        schema = 1,
        split_policy = 1,
        collision_position_policy = :new_integer,
        collision_velocity_policy = :carried_fast_half,
        inner = _checkpoint_state(st.inner, st.species),
        dt = st.dt,
        coulomb_log = copy(st.coulomb_log),
        epsilon0 = st.epsilon0,
        weight_to_number = st.weight_to_number,
        u_floor = st.u_floor,
        rng = deepcopy(st.rng),
        rng_algorithm = typeof(st.rng),
        runtime = _collisional_runtime(),
        phase = :ready,
        initialized = true,
        source_time = st.source_time,
        source_count = st.source_count,
    )
end

function _validate_checkpoint_state(s, st::CollisionalEMPICStepper{D,T}) where {D,T}
    _validate_checkpoint_container(s, _COLLISIONAL_EMPIC_CHECKPOINT_FIELDS)
    s.integrator === :collisional_empic && s.schema === 1 && s.split_policy === 1 ||
        throw(ArgumentError("unsupported collisional EMPIC checkpoint schema or split policy"))
    s.collision_position_policy === :new_integer &&
        s.collision_velocity_policy === :carried_fast_half ||
        throw(ArgumentError("collisional EMPIC checkpoint collision phase does not match"))
    isequal(s.runtime, _collisional_runtime()) ||
        throw(ArgumentError("collisional EMPIC checkpoint requires the same Julia runtime and platform"))
    s.phase === :ready && s.initialized === true ||
        throw(ArgumentError("collisional EMPIC checkpoint must be an initialized completed boundary"))
    _validate_checkpoint_state(s.inner, st.inner)
    s.inner.integrator === :empic_collection &&
        s.inner.particle_layout === :collection &&
        s.inner.schema_version === 1 &&
        s.inner.initialized === true &&
        s.inner.relativistic === false &&
        s.inner.n_sub === 1 &&
        all(s.inner.mobile) &&
        _particle_stored_zero(s.inner.n0) &&
        s.inner.neutrality === :require ||
        throw(ArgumentError("collisional EMPIC checkpoint requires its fixed generic electromagnetic collection"))
    isequal(s.inner.c, st.c) &&
        isequal(s.inner.shape, st.shape) &&
        isequal(s.inner.charge_rtol, st.charge_rtol) ||
        throw(ArgumentError("collisional EMPIC checkpoint field configuration must match the destination"))
    for name in (:dt, :epsilon0, :weight_to_number, :u_floor)
        value = getproperty(s, name)
        value isa T && isequal(value, getfield(st, name)) ||
            throw(ArgumentError("collisional EMPIC checkpoint $name must match the destination"))
    end
    # Stored weights deposit with unit Maxwell coefficient: the collision
    # permittivity must equal the physical count conversion in these units.
    isfinite(s.epsilon0) && s.epsilon0 > zero(T) &&
        isequal(s.epsilon0, s.weight_to_number) ||
        throw(ArgumentError("collisional EMPIC checkpoint requires equal positive epsilon0 and weight_to_number"))
    typeof(s.coulomb_log) === typeof(st.coulomb_log) && isequal(s.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional EMPIC checkpoint Coulomb logarithms must match the destination"))
    s.source_time isa T &&
        isequal(s.source_time, s.inner.time) &&
        s.source_count isa Int &&
        s.source_count == s.inner.step ||
        throw(ArgumentError("collisional EMPIC checkpoint source ledger must match its field clock"))
    (
        s.source_count == 0 ?
        iszero(s.source_time) && iszero(s.inner.last_dt) && iszero(s.inner.last_sub_dt) :
        isequal(s.inner.last_dt, st.dt) && isequal(s.inner.last_sub_dt, st.dt)
    ) || throw(ArgumentError("collisional EMPIC checkpoint has inconsistent fixed timestep history"))
    length(s.inner.species) == length(st.identities) ||
        throw(ArgumentError("collisional EMPIC checkpoint species count must match the destination"))
    for a in eachindex(st.identities)
        p = s.inner.species[a]
        isequal((length(p.weight), p.q, p.m), st.identities[a]) || throw(
            ArgumentError("collisional EMPIC checkpoint particle counts, charges and masses must match"),
        )
        all(array -> array isa Vector, (p.x..., p.v..., p.weight, p.id, p.tag)) ||
            throw(ArgumentError("collisional EMPIC checkpoint requires standard Vector particle storage"))
    end
    all(array -> array isa Array, _empic_collection_checkpoint_arrays(st.inner, s.inner.fields)) ||
        throw(ArgumentError("collisional EMPIC checkpoint requires standard Array field storage"))
    s.rng_algorithm === typeof(st.rng) && typeof(s.rng) === typeof(st.rng) ||
        throw(ArgumentError("collisional EMPIC checkpoint RNG algorithm must match the destination"))
    _validate_collisional_empic_rngs(st)
    _validate_collisional_rng(s.rng)
    # An RNG can retain mutable seed metadata. Reject links before deepcopy can
    # hide them, including links to an uninitialized destination's particle leaves.
    _validate_checkpoint_configuration(
        (
            inner = s.inner,
            logs = s.coulomb_log,
            live_species = st.species,
            identities = st.identities,
            workspace = st.workspace,
            live_logs = st.coulomb_log,
            log_binding = st.log_binding,
            live_rng = st.rng,
            candidate_rng = st.candidate_rng,
        ),
        st.inner,
        (s.rng,),
    )
    return nothing
end

function _restore_checkpoint_state!(st::CollisionalEMPICStepper, species::AbstractVector, s)
    _validate_collisional_empic_binding(st, species; values = false)
    _validate_checkpoint_state(s, st)
    # Fixed ordinary Vector destinations and all three RNG structural states
    # have passed preflight. Detach the entire incoming graph before publication.
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    # The canonical nested restore completes its own validation, detachment and
    # capacity reservation before its first public store. Counts cannot change.
    # Stock RNG copy! then only copies the already validated fixed-size state and
    # caches, plus typed metadata; it cannot discover a size mismatch afterward.
    _restore_checkpoint_state!(st.inner, st.species, saved.inner)
    copy!(st.rng, saved.rng)
    st.source_time = saved.source_time
    st.source_count = saved.source_count
    st.phase = :ready
    return st
end

function save_checkpoint(
    path::AbstractString,
    st::CollisionalEMPICStepper,
    species::AbstractVector,
)
    _acquire_collisional_empic!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy[] = false
    end
end

function load_checkpoint!(
    st::CollisionalEMPICStepper,
    species::AbstractVector,
    path::AbstractString,
)
    _acquire_collisional_empic!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy[] = false
    end
end
