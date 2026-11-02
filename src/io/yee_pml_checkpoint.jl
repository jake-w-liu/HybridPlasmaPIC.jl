const _YEE_PML_CHECKPOINT_FIELDS = (
    :integrator,
    :schema_version,
    :pml_policy,
    :profile_policy,
    :inner,
    :memory,
    :dt,
    :cells,
    :peak_rate,
    :profile,
    :initialized,
    :phase,
    :ledger_step,
    _YEE_PML_LEDGER_FIELDS...,
)

function _checkpoint_state(model::YeePML)
    _yp_validate_config(model)
    _yp_validate_storage(model)
    _yp_boundary(model)
    return (
        integrator = :yee_ade_pml,
        schema_version = 1,
        pml_policy = :causal_directional_kappa1_alpha0_midpoint_kdk_original_h_v1,
        profile_policy = :index_staggered_polynomial_v1,
        inner = _checkpoint_state(model.inner),
        memory = (E = model._memory.E, B = model._memory.B),
        dt = model.dt,
        cells = model.cells,
        peak_rate = model.peak_rate,
        profile = model.profile,
        initialized = true,
        phase = :ready,
        ledger_step = model.ledger_step,
        initial_energy = model.initial_energy,
        last_exchange = model.last_exchange,
        total_exchange = model.total_exchange,
        last_field_defect = model.last_field_defect,
        total_field_defect = model.total_field_defect,
    )
end

function _yp_checkpoint_memory(s, model::YeePML{D,T}) where {D,T}
    _validate_checkpoint_container(s, (:E, :B))
    length(keys(s)) == 2 ||
        throw(ArgumentError("PML checkpoint contains unsupported memory families"))
    s.E isa NTuple{3,NTuple{2,Union{Nothing,Array{T,D}}}} &&
        s.B isa NTuple{3,NTuple{2,Union{Nothing,Array{T,D}}}} || throw(
        ArgumentError(
            "PML checkpoint memories must use ordinary matching CPU arrays or absent slots",
        ),
    )
    memory = _YeePMLMemory{D,T}(s.E, s.B)
    _yp_memory_storage(memory, model.inner.g, model.cells)
    return memory
end

function _validate_checkpoint_state(s, model::YeePML{D,T}) where {D,T}
    _yp_validate_config(model; cache = false)
    _validate_checkpoint_container(s, _YEE_PML_CHECKPOINT_FIELDS)
    length(keys(s)) == length(_YEE_PML_CHECKPOINT_FIELDS) ||
        throw(ArgumentError("PML checkpoint contains unsupported fields"))
    s.integrator === :yee_ade_pml &&
        s.schema_version === 1 &&
        s.pml_policy === :causal_directional_kappa1_alpha0_midpoint_kdk_original_h_v1 &&
        s.profile_policy === :index_staggered_polynomial_v1 ||
        throw(ArgumentError("unsupported PML checkpoint schema or policy"))
    _validate_checkpoint_state(s.inner, model.inner)
    for name in (:dt, :cells, :peak_rate, :profile)
        saved, expected = getproperty(s, name), getfield(model, name)
        typeof(saved) === typeof(expected) && isequal(saved, expected) ||
            throw(ArgumentError("PML checkpoint $name does not match destination"))
    end
    memory = _yp_checkpoint_memory(s.memory, model)
    rates = _yee_profile_axis_cache(model.inner.g, s.cells, s.peak_rate, s.profile, one(_YeeLayerQ))
    _yp_memory_values(memory, model.inner.g, rates, s.inner.step)
    _yee_check_independent_arrays((s.inner.fields.E..., s.inner.fields.B..., _yp_arrays(memory)...))
    fields = YeeFields{D,T}(s.inner.fields.E, s.inner.fields.B)
    energy = _yee_layer_energy(fields, model.inner)
    _yp_validate_budget(s, energy, s.inner.time, s.inner.step, s.dt, T)
    return nothing
end

function _restore_checkpoint_state!(model::YeePML, s)
    # Old field/memory/clock/ledger/cache values are replaceable. Only fixed
    # representation/configuration and destination storage gate this restore.
    _yp_validate_config(model; cache = false)
    _yp_validate_storage(model; check_pec = false, cache = false)
    _validate_checkpoint_state(s, model)
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, model)
    memory = _yp_checkpoint_memory(saved.memory, model)
    rates = _yee_profile_axis_cache(
        model.inner.g,
        saved.cells,
        saved.peak_rate,
        saved.profile,
        one(_YeeLayerQ),
    )
    # The whole source graph is detached. Canonical nested restore completes its
    # own detachment/checks before its first store. Fixed-size copies and typed
    # already-prepared assignments only follow that physical publication.
    _restore_checkpoint_state!(model.inner, saved.inner)
    _yp_copy_memory!(model._memory, memory)
    model._axis_rate = rates
    model.initial_energy = saved.initial_energy
    model.last_exchange = saved.last_exchange
    model.total_exchange = saved.total_exchange
    model.last_field_defect = saved.last_field_defect
    model.total_field_defect = saved.total_field_defect
    model.ledger_step = saved.ledger_step
    model.initialized = true
    model.phase = :ready
    return model
end

"""
    save_checkpoint(path, model::YeePML)

Atomically save the distinct causal directional ADE PML schema, fixed staggered
profile/timestep, synchronized vacuum fields/clock, every live memory and exact
signed exchange/field-defect ledgers. Ready count-zero memories must all be
stored zero; advanced memories remain actual state. Derived caches/candidates
are omitted. No scalar-absorber, forcing, CPML tuning or particle state is inferred.
"""
function save_checkpoint(path::AbstractString, model::YeePML)
    _yp_acquire!(model)
    try
        return _save_checkpoint(path, model, _checkpoint_state(model))
    finally
        model.busy[] = false
    end
end

"""
    load_checkpoint!(model::YeePML, path)

Validate matching configuration and complete incoming PML history, detach all
fields/memories together, rebuild derived caches and restore verbatim. Old
replaceable destination values may be stale. Physical/memory arrays and clock
Refs retain identities; no zero-memory migration, field projection or history
reconstruction occurs. Trusted Serialization/ordinary-copy/exclusive-access
premises apply. Runtime continuation acceptance is deferred.
"""
function load_checkpoint!(model::YeePML, path::AbstractString)
    _yp_acquire!(model)
    try
        return _restore_checkpoint_state!(model, deserialize(path))
    finally
        model.busy[] = false
    end
end

"""
    archive_run(path, model::YeePML; rng_seed, kwargs...)

Archive the same complete ADE PML state through the checksummed run container.
The seed is provenance only. Default units/boundary metadata describe the actual
vacuum/PEC/periodic grid and the distinct directional kappa1/alpha0 PML; no RNG,
source waveform or measured reflection result is supplied.
"""
function archive_run(
    path::AbstractString,
    model::YeePML;
    normalization::AbstractString = "epsilon0=1;mu0=1/c^2;directional_ADE_PML_kappa1_alpha0",
    boundary_desc::AbstractString = join(
        map(p -> p ? "periodic" : "pec", model.inner.g.periodic),
        ",",
    ) * ";directional_ADE_PML",
    kwargs...,
)
    _yp_acquire!(model)
    try
        return _archive_run(
            path,
            model,
            _checkpoint_state(model);
            normalization,
            boundary_desc,
            kwargs...,
        )
    finally
        model.busy[] = false
    end
end
