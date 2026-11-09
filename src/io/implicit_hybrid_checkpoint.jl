# This family persists its synchronous physical state, never an explicit-stepper
# surrogate or a nonlinear solver history/initial guess.
const _IMPLICIT_HYBRID_CHECKPOINT_FIELDS = (
    :integrator,
    :schema,
    :configuration,
    :phase,
    :initialized,
    :species,
    :fields,
    :floor_count,
    :time,
    :step,
    :last_dt,
    :last_report,
)

function _implicit_hybrid_checkpoint_keys(s, expected)
    _validate_checkpoint_container(s, expected)
    length(keys(s)) == length(expected) ||
        throw(ArgumentError("implicit hybrid checkpoint contains unsupported fields"))
    return nothing
end

function _implicit_hybrid_checkpoint_report(
    r,
    st::ImplicitHybridStepper{D,T},
    time,
    step,
    h,
) where {D,T}
    _implicit_hybrid_phase_values(true, time, step, h, r)
    _implicit_hybrid_report_policy(st, r)
    # These are history metadata checks, not a re-solve of an unavailable old
    # interval or evidence of the original calculation's numerical accuracy.
    return nothing
end

function _checkpoint_state(st::ImplicitHybridStepper, species::AbstractVector)
    _implicit_hybrid_binding(st, species)
    _implicit_hybrid_boundary(st)
    configuration = _implicit_hybrid_configuration(st)
    state = (
        integrator = :implicit_hybrid_midpoint,
        schema = st.model.closure isa EvolvingElectrons ? 2 : 1,
        configuration = configuration,
        phase = :synchronous_integer,
        initialized = true,
        species = [_checkpoint_particle_state(ps) for ps in species],
        fields = _hybrid_field_state(st.fields),
        floor_count = st.fields.floor_count[],
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        last_report = st.last_report,
    )
    # The evolving receipt is a known exact numeric record, validated by its
    # report owner. It is not opaque configuration or retained storage metadata.
    metadata =
        st.model.closure isa EvolvingElectrons ? configuration : (configuration, st.last_report)
    _validate_checkpoint_configuration(state, st, metadata)
    # Capture itself is detached so the shared writer's validation observes the
    # same independent ordinary arrays that a deserialized payload must own.
    return deepcopy(state)
end

function _validate_checkpoint_state(s, st::ImplicitHybridStepper{D,T}) where {D,T}
    _implicit_hybrid_binding(st, st.owners; values = false, counts = false)
    _implicit_hybrid_checkpoint_keys(s, _IMPLICIT_HYBRID_CHECKPOINT_FIELDS)
    s.integrator === :implicit_hybrid_midpoint &&
        s.schema === (st.model.closure isa EvolvingElectrons ? 2 : 1) &&
        s.phase === :synchronous_integer &&
        s.initialized === true ||
        throw(ArgumentError("unsupported implicit hybrid checkpoint family, schema or phase"))
    configuration = _implicit_hybrid_configuration(st)
    typeof(s.configuration) === typeof(configuration) && isequal(s.configuration, configuration) ||
        throw(
            ArgumentError(
                "implicit hybrid checkpoint fixed configuration must match the destination",
            ),
        )
    s.time isa T && s.step isa Int && s.last_dt isa T ||
        throw(ArgumentError("implicit hybrid checkpoint clock/history types differ"))
    _implicit_hybrid_checkpoint_report(s.last_report, st, s.time, s.step, s.last_dt)
    s.floor_count isa Int && s.floor_count == 0 ||
        throw(ArgumentError("implicit hybrid checkpoint requires an inactive density floor"))
    references = _hybrid_field_state(st.fields)
    _implicit_hybrid_checkpoint_keys(s.fields, keys(references))
    _validate_checkpoint_fields(s.fields, references)
    field_arrays =
        Tuple(a for value in values(s.fields) for a in (value isa Tuple ? value : (value,)))
    all(a -> a isa Array{T,D} && all(isfinite, a), field_arrays) || throw(
        ArgumentError("implicit hybrid checkpoint fields must be finite ordinary matching Arrays"),
    )
    all(n -> n >= T(st.model.nfloor), s.fields.n) || throw(
        ArgumentError("implicit hybrid checkpoint density must satisfy the strict floor bound"),
    )
    st.model.closure isa EvolvingElectrons && _validate_evolving_pressure(
        s.fields.pe,
        st.g;
        context = "implicit hybrid checkpoint electron pressure",
    )
    s.species isa Vector && length(s.species) == length(st.owners) ||
        throw(ArgumentError("implicit hybrid checkpoint requires its ordered population count"))
    populations = _ElectrostaticHostParticle{D,T}[]
    for (index, p) in enumerate(s.species)
        _implicit_hybrid_checkpoint_keys(p, _CHECKPOINT_PARTICLE_FIELDS)
        ps = _checkpoint_species(p, Val(D), T)
        ps isa _ElectrostaticHostParticle{D,T} &&
            isequal((ps.q, ps.m), st._species_scalars[index]) || throw(
            ArgumentError(
                "implicit hybrid checkpoint requires Vector leaves and fixed cached charges/masses",
            ),
        )
        push!(populations, ps)
    end
    _validate_hybrid_species(populations)
    arrays = (s.species, field_arrays..., _implicit_hybrid_particle_arrays(populations)...)
    live = (
        _implicit_hybrid_field_arrays(st.fields)...,
        st.owners,
        st._species_scalars,
        st._population_counts,
        _implicit_hybrid_particle_arrays(st.owners)...,
        _implicit_hybrid_work_arrays(st.work)...,
        st.g.cbuf,
        st.g.tbuf,
        st.g.abuf,
        st.g.sbuf,
        st.g.midx...,
        st.g.kfull...,
        st.g.ik...,
        st.g.kvec...,
    )
    for i in eachindex(arrays)
        for j = 1:(i-1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(
                ArgumentError(
                    "implicit hybrid checkpoint physical storage roles must be independent",
                ),
            )
        end
        for a in live
            _particle_arrays_mightalias(arrays[i], a) && throw(
                ArgumentError(
                    "implicit hybrid checkpoint must be detached from live public/private storage",
                ),
            )
        end
    end
    metadata =
        st.model.closure isa EvolvingElectrons ? s.configuration : (s.configuration, s.last_report)
    _validate_checkpoint_configuration(s, st, metadata)
    return nothing
end

function _restore_checkpoint_state!(
    st::ImplicitHybridStepper{D,T},
    species::AbstractVector,
    s,
) where {D,T}
    # Old physical/qm/count/clock/phase values are replaceable. Fixed owner
    # identity, ordinary storage, cached binding, grid and solver policy are not.
    _implicit_hybrid_binding(st, species; values = false, counts = false)
    _validate_checkpoint_state(s, st)
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    incoming =
        _ElectrostaticHostParticle{D,T}[_checkpoint_species(p, Val(D), T) for p in saved.species]
    prepared, counts = _implicit_hybrid_replacements(incoming)
    fields = HybridFields{D,T}(st.g.n)
    _restore_checkpoint_fields!(_hybrid_field_state(fields), saved.fields)
    fields.floor_count[] = saved.floor_count
    # Reserve all stage leaves, gather/nonlinear vectors and concrete GMRES
    # storage before the first public store. This does not reconstruct physics.
    _implicit_hybrid_reserve!(st, incoming)
    return _implicit_hybrid_commit!(
        st,
        prepared,
        counts,
        fields,
        true,
        saved.time,
        saved.step,
        saved.last_dt,
        saved.last_report,
    )
end

"""
    save_checkpoint(path, st::ImplicitHybridStepper, species)

Atomically save `:implicit_hybrid_midpoint` (isothermal schema1 or evolving
pressure schema2): complete public
fields, ordered particles, synchronous clock/history, exact fixed physics and
solver policy, and the last accepted scalar report. Nonlinear/FFT scratch is
omitted. Only initialized completed boundaries can be saved.
"""
function save_checkpoint(path::AbstractString, st::ImplicitHybridStepper, species::AbstractVector)
    _acquire_implicit_hybrid!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy = false
    end
end

"""
    load_checkpoint!(st::ImplicitHybridStepper, species, path)

Validate matching fixed configuration and complete incoming state, detach and
prepare every replacement and solver buffer, then publish fields/particles and
history verbatim. Stale destination physical values and phase are allowed. No
half-time migration, initialization, redeposition or implicit solve is performed.
"""
function load_checkpoint!(st::ImplicitHybridStepper, species::AbstractVector, path::AbstractString)
    _acquire_implicit_hybrid!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy = false
    end
end

function archive_run(
    path::AbstractString,
    st::ImplicitHybridStepper,
    species::AbstractVector{<:ParticleSet};
    kwargs...,
)
    _acquire_implicit_hybrid!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
    finally
        st.busy = false
    end
end

"Verify the run-container checksum and perform the same prepared synchronous restore."
function load_archive!(st::ImplicitHybridStepper, species::AbstractVector, path::AbstractString)
    _acquire_implicit_hybrid!(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, species, loaded.state)
    finally
        st.busy = false
    end
end
