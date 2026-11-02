# Complete ES–Coulomb boundaries use the existing nested ES physical state.
const _COLLISIONAL_ESPIC_CHECKPOINT_FIELDS = (
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

function _checkpoint_state(st::CollisionalElectrostaticPIC, species::AbstractVector)
    _validate_collisional_espic_binding(st, species)
    _validate_collisional_espic_boundary(st)
    _validate_collisional_espic_rngs(st)
    return (
        integrator = :collisional_espic,
        schema = 1,
        split_policy = 1,
        collision_position_policy = :new_integer,
        collision_velocity_policy = :carried_half,
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

function _validate_checkpoint_state(s, st::CollisionalElectrostaticPIC{D,T}) where {D,T}
    _validate_checkpoint_container(s, _COLLISIONAL_ESPIC_CHECKPOINT_FIELDS)
    s.integrator === :collisional_espic && s.schema === 1 && s.split_policy === 1 ||
        throw(ArgumentError("unsupported collisional ES checkpoint schema or split policy"))
    s.collision_position_policy === :new_integer &&
        s.collision_velocity_policy === :carried_half ||
        throw(ArgumentError("collisional ES checkpoint collision phase does not match"))
    isequal(s.runtime, _collisional_runtime()) ||
        throw(ArgumentError("collisional ES checkpoint requires the same Julia runtime and platform"))
    s.phase === :ready && s.initialized === true ||
        throw(ArgumentError("collisional ES checkpoint must be an initialized completed boundary"))
    _validate_checkpoint_state(s.inner, st.inner)
    s.inner.integrator === :espic_collection &&
        s.inner.particle_layout === :collection &&
        s.inner.schema_version === 1 &&
        s.inner.initialized === true &&
        _particle_stored_zero(s.inner.n0) &&
        s.inner.neutrality === :require &&
        isequal(s.inner.charge_rtol, st.charge_rtol) ||
        throw(ArgumentError("collisional ES checkpoint requires its fixed electrostatic collection policy"))
    for name in (:dt, :epsilon0, :weight_to_number, :u_floor)
        value = getproperty(s, name)
        value isa T && isequal(value, getfield(st, name)) ||
            throw(ArgumentError("collisional ES checkpoint $name must match the destination"))
    end
    # Poisson and field energy use unit coefficients for stored weights.
    isfinite(s.epsilon0) && s.epsilon0 > zero(T) &&
        isequal(s.epsilon0, s.weight_to_number) ||
        throw(ArgumentError("collisional ES checkpoint requires equal positive epsilon0 and weight_to_number"))
    typeof(s.coulomb_log) === typeof(st.coulomb_log) && isequal(s.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional ES checkpoint Coulomb logarithms must match the destination"))
    # ES has no inner clock. The wrapper clock is authoritative; do not infer
    # time from count*h, whose rounding history differs from stored additions.
    _validate_checkpoint_clock_state((time = s.source_time, step = s.source_count), T)
    (
        s.source_count == 0 ?
        _particle_stored_zero(s.source_time) && !s.inner.primed &&
        _particle_stored_zero(s.inner.last_dt) :
        s.source_time > zero(T) && s.inner.primed && isequal(s.inner.last_dt, st.dt)
    ) || throw(ArgumentError("collisional ES checkpoint source clock and velocity phase disagree"))
    length(s.inner.species) == length(st.identities) &&
        axes(s.inner.species, 1) == Base.OneTo(length(st.identities)) ||
        throw(ArgumentError("collisional ES checkpoint requires matching ordered population records"))
    for a in eachindex(st.identities)
        p = s.inner.species[a]
        isequal((length(p.weight), p.q, p.m), st.identities[a]) || throw(
            ArgumentError("collisional ES checkpoint particle counts, charges and masses must match"),
        )
        all(array -> array isa Vector, (p.x..., p.v..., p.weight, p.id, p.tag)) ||
            throw(ArgumentError("collisional ES checkpoint requires standard Vector particle storage"))
    end
    fields = (s.inner.E..., s.inner.ne)
    for (i, array) in enumerate(fields)
        array isa Array ||
            throw(ArgumentError("collisional ES checkpoint requires standard Array field storage"))
        for j = 1:(i-1)
            _particle_arrays_mightalias(array, fields[j]) &&
                throw(ArgumentError("collisional ES checkpoint physical fields must be independent"))
        end
    end
    s.rng_algorithm === typeof(st.rng) && typeof(s.rng) === typeof(st.rng) ||
        throw(ArgumentError("collisional ES checkpoint RNG algorithm must match the destination"))
    _validate_collisional_espic_rngs(st)
    _validate_collisional_rng(s.rng)
    # Reject retained incoming/live RNG links before deepcopy can hide them.
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

function _restore_checkpoint_state!(st::CollisionalElectrostaticPIC, species::AbstractVector, s)
    _validate_collisional_espic_binding(st, species; values = false)
    _validate_checkpoint_state(s, st)
    # Detach every incoming leaf, including RNG metadata, before publication.
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    # Canonical ES restore detaches again, validates destinations and reserves
    # owner capacity before its first store. Fixed matching Vector counts avoid
    # particle resizing. All stock RNG state/cache sizes passed preflight.
    _restore_checkpoint_state!(st.inner, st.species, saved.inner)
    copy!(st.rng, saved.rng)
    st.source_time = saved.source_time
    st.source_count = saved.source_count
    st.phase = :ready
    return st
end

function save_checkpoint(
    path::AbstractString,
    st::CollisionalElectrostaticPIC,
    species::AbstractVector,
)
    _acquire_collisional_espic!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy[] = false
    end
end

function load_checkpoint!(
    st::CollisionalElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    _acquire_collisional_espic!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy[] = false
    end
end
