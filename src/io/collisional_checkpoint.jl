# Versioned stochastic composition state; publication remains owned by checkpoint.jl.
const _COLLISIONAL_CHECKPOINT_FIELDS = (
    :integrator,
    :schema,
    :inner,
    :fields,
    :floor_count,
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

_collisional_runtime() =
    (julia = VERSION, word_size = Sys.WORD_SIZE, arch = Sys.ARCH, kernel = Sys.KERNEL)

function _checkpoint_state(st::CollisionalHybridStepper, species::AbstractVector)
    _validate_collisional_binding(st, species)
    _validate_collisional_boundary(st)
    return (
        integrator = :collisional_hybrid,
        schema = 1,
        inner = _checkpoint_state(st.inner, st.species),
        fields = _hybrid_field_state(st.inner.fields),
        floor_count = st.inner.fields.floor_count[],
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

function _validate_checkpoint_state(s, st::CollisionalHybridStepper{D,T}) where {D,T}
    _validate_checkpoint_container(s, _COLLISIONAL_CHECKPOINT_FIELDS)
    s.integrator === :collisional_hybrid && s.schema === 1 ||
        throw(ArgumentError("unsupported collisional checkpoint schema"))
    isequal(s.runtime, _collisional_runtime()) ||
        throw(ArgumentError("collisional checkpoint requires the same Julia runtime and platform"))
    s.phase === :ready && s.initialized === true ||
        throw(ArgumentError("collisional checkpoint must be an initialized completed boundary"))
    _validate_checkpoint_state(s.inner, st.inner)
    haskey(s.inner, :particle_layout) && s.inner.particle_layout === :collection ||
        throw(ArgumentError("collisional checkpoint requires an ordered species collection"))
    isequal(s.inner.model, st.model) && isequal(s.inner.shape, st.shape) ||
        throw(ArgumentError("collisional checkpoint model and shape must match the destination"))
    for name in (:dt, :epsilon0, :weight_to_number, :u_floor)
        value = getproperty(s, name)
        value isa T && isequal(value, getfield(st, name)) ||
            throw(ArgumentError("collisional checkpoint $name must match the destination"))
    end
    typeof(s.coulomb_log) === typeof(st.coulomb_log) && isequal(s.coulomb_log, st.log_binding) ||
        throw(ArgumentError("collisional checkpoint Coulomb logarithms must match the destination"))
    s.source_time isa T &&
        isequal(s.source_time, s.inner.time) &&
        s.source_count isa Int &&
        s.source_count == s.inner.step ||
        throw(ArgumentError("collisional checkpoint source ledger must match its field clock"))
    (
        s.source_count == 0 ? iszero(s.source_time) && iszero(s.inner.last_dt) :
        isequal(s.inner.last_dt, st.dt)
    ) || throw(ArgumentError("collisional checkpoint has inconsistent timestep history"))
    length(s.inner.species) == length(st.identities) ||
        throw(ArgumentError("collisional checkpoint species count must match the destination"))
    for a in eachindex(st.identities)
        p = s.inner.species[a]
        isequal((length(p.weight), p.q, p.m), st.identities[a]) || throw(
            ArgumentError("collisional checkpoint particle counts, charges and masses must match"),
        )
        for array in (p.x..., p.v..., p.weight, p.id, p.tag)
            array isa Vector || throw(
                ArgumentError("collisional checkpoint requires standard Vector particle storage"),
            )
        end
    end
    _validate_checkpoint_fields(s.fields, _hybrid_field_state(st.inner.fields))
    for components in values(s.fields)
        arrays = components isa Tuple ? components : (components,)
        for array in arrays
            array isa Array ||
                throw(ArgumentError("collisional checkpoint requires standard Array field storage"))
            all(isfinite, array) ||
                throw(ArgumentError("collisional checkpoint fields must be finite"))
        end
    end
    isequal(s.fields.B, s.inner.B) && isequal(s.fields.E, s.inner.E) ||
        throw(ArgumentError("collisional checkpoint inner and public fields disagree"))
    s.floor_count isa Int && s.floor_count >= 0 ||
        throw(ArgumentError("collisional checkpoint density floor count must be nonnegative"))
    s.rng_algorithm === typeof(st.rng) && typeof(s.rng) === typeof(st.rng) ||
        throw(ArgumentError("collisional checkpoint RNG algorithm must match the destination"))
    _validate_collisional_rng(s.rng)
    # RNG metadata may be mutable (notably MersenneTwister.seed). Reuse the
    # canonical graph boundary to forbid links into either incoming or live state.
    _validate_checkpoint_configuration(
        (inner = s.inner, fields = s.fields, logs = s.coulomb_log, live_rng = st.rng),
        st.inner,
        (s.rng,),
    )
    return nothing
end

function _restore_checkpoint_state!(st::CollisionalHybridStepper, species::AbstractVector, s)
    _validate_collisional_binding(st, species)
    _validate_checkpoint_state(s, st)
    populations = st.species
    count = _validate_hybrid_checkpoint_restore(st.inner, populations, s.inner)
    # Allocate and detach every incoming object before publishing values. The
    # canonical particle/field restore then copies into the existing public arrays.
    saved = deepcopy(s)
    metadata = [(copy(p.weight), copy(p.id), copy(p.tag)) for p in populations]
    _resize_hybrid_particle_workspaces!(st.inner, count)
    _save_collisional_transaction!(st)
    previous_phase = st.phase
    previous_ledger = (st.source_time, st.source_count)
    st.phase = :inflight
    completed = false
    try
        _restore_validated_hybrid_checkpoint!(st.inner, populations, saved.inner, count)
        _restore_checkpoint_fields!(_hybrid_field_state(st.inner.fields), saved.fields)
        st.inner.fields.floor_count[] = saved.floor_count
        copy!(st.rng, saved.rng)
        st.source_time = saved.source_time
        st.source_count = saved.source_count
        completed = true
        return st
    finally
        if !completed
            _restore_collisional_transaction!(st)
            for (p, (weight, id, tag)) in zip(populations, metadata)
                copyto!(p.weight, weight)
                copyto!(p.id, id)
                copyto!(p.tag, tag)
            end
            st.source_time, st.source_count = previous_ledger
        end
        st.phase = completed ? :ready : previous_phase
    end
end

function save_checkpoint(
    path::AbstractString,
    st::CollisionalHybridStepper,
    species::AbstractVector,
)
    _acquire_collisional!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy[] = false
    end
end

function load_checkpoint!(
    st::CollisionalHybridStepper,
    species::AbstractVector,
    path::AbstractString,
)
    _acquire_collisional!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy[] = false
    end
end
