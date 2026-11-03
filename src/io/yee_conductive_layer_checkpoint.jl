const _YEE_LAYER_CHECKPOINT_FIELDS = (
    :integrator,
    :schema_version,
    :layer_policy,
    :profile_policy,
    :inner,
    :dt,
    :cells,
    :peak_rate,
    :profile,
    :initialized,
    :phase,
    :ledger_step,
    _YEE_LAYER_LEDGER_FIELDS...,
)

function _checkpoint_state(model::YeeConductiveLayer)
    _yee_layer_validate_config(model)
    _yee_layer_validate_storage(model)
    _yee_layer_boundary(model)
    return (
        integrator = :yee_conductive_absorber,
        schema_version = 1,
        layer_policy = :scalar_matched_midpoint_strang_original_dt_v1,
        profile_policy = :index_staggered_polynomial_v1,
        inner = _checkpoint_state(model.inner),
        dt = model.dt,
        cells = model.cells,
        peak_rate = model.peak_rate,
        profile = model.profile,
        initialized = true,
        phase = :ready,
        ledger_step = model.ledger_step,
        initial_energy = model.initial_energy,
        last_loss = model.last_loss,
        total_loss = model.total_loss,
        last_work = model.last_work,
        total_work = model.total_work,
        last_defect = model.last_defect,
        total_defect = model.total_defect,
    )
end

function _validate_checkpoint_state(s, model::YeeConductiveLayer{D,T}) where {D,T}
    _yee_layer_mode()
    _validate_checkpoint_container(s, _YEE_LAYER_CHECKPOINT_FIELDS)
    length(keys(s)) == length(_YEE_LAYER_CHECKPOINT_FIELDS) ||
        throw(ArgumentError("conductive layer checkpoint contains unsupported fields"))
    s.integrator === :yee_conductive_absorber &&
        s.schema_version === 1 &&
        s.layer_policy === :scalar_matched_midpoint_strang_original_dt_v1 &&
        s.profile_policy === :index_staggered_polynomial_v1 ||
        throw(ArgumentError("unsupported conductive layer checkpoint schema or policy"))
    _validate_checkpoint_state(s.inner, model.inner)
    for name in (:dt, :cells, :peak_rate, :profile)
        saved, expected = getproperty(s, name), getfield(model, name)
        typeof(saved) === typeof(expected) && isequal(saved, expected) ||
            throw(ArgumentError("conductive layer checkpoint $name does not match destination"))
    end
    _yee_layer_config(model.inner.g, s.dt, s.cells, s.peak_rate, s.profile)
    s.dt < maxwell_cfl(model) ||
        throw(ArgumentError("conductive layer checkpoint timestep must be strictly sub-CFL"))
    fields = YeeFields{D,T}(s.inner.fields.E, s.inner.fields.B)
    energy = _yee_layer_energy(fields, model.inner)
    _yee_layer_validate_budget(s, energy, s.inner.time, s.inner.step, s.dt, T)
    return nothing
end

function _restore_checkpoint_state!(model::YeeConductiveLayer, s)
    # Current physical, clock and ledger values will be replaced. Derived cache
    # contents are omitted from the payload and rebuilt before any public store.
    _yee_layer_validate_config(model; cache = false)
    _yee_layer_validate_storage(model; check_pec = false, cache = false)
    _validate_checkpoint_state(s, model)
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, model)
    cache =
        _yee_layer_axis_cache(model.inner.g, saved.dt, saved.cells, saved.peak_rate, saved.profile)
    # Canonical nested restore detaches again and completes all finite/shape/
    # alias/PEC checks before its first ordinary field copy. Only typed stores
    # of already prepared metadata and caches follow that publication.
    _restore_checkpoint_state!(model.inner, saved.inner)
    model._axis_r = cache
    model.initial_energy = saved.initial_energy
    model.last_loss = saved.last_loss
    model.total_loss = saved.total_loss
    model.last_work = saved.last_work
    model.total_work = saved.total_work
    model.last_defect = saved.last_defect
    model.total_defect = saved.total_defect
    model.ledger_step = saved.ledger_step
    model.initialized = true
    model.phase = :ready
    return model
end

"""
    save_checkpoint(path, model::YeeConductiveLayer)

Atomically save the distinct scalar-conductivity policy, fixed profile/timestep,
canonical synchronized vacuum fields/clock and all exact loss/work/defect ledgers.
Only initialized completed boundaries save. Derived axis caches and staging are
not physical checkpoint state. This schema does not represent PML memory or a
particle/source waveform. Runtime/restart acceptance is deferred.
"""
function save_checkpoint(path::AbstractString, model::YeeConductiveLayer)
    _yee_layer_acquire!(model)
    try
        return _save_checkpoint(path, model, _checkpoint_state(model))
    finally
        model.busy[] = false
    end
end

"""
    load_checkpoint!(model::YeeConductiveLayer, path)

Restore matching grid/boundaries/speed, fixed conductivity configuration and
complete synchronized fields/budget. Incoming nested fields and exact ledger
identity validate and detach before publication. Old destination field/clock/
ledger values may be stale; storage/configuration must match. Physical arrays/
Refs retain identities, derived caches rebuild, and no field projection occurs.
Trusted Serialization and the ordinary-copy/exclusive-access premises apply.
"""
function load_checkpoint!(model::YeeConductiveLayer, path::AbstractString)
    _yee_layer_acquire!(model)
    try
        return _restore_checkpoint_state!(model, deserialize(path))
    finally
        model.busy[] = false
    end
end

"""
    archive_run(path, model::YeeConductiveLayer; rng_seed, kwargs...)

Archive the same complete field-layer payload through the existing checksummed
container. The seed is provenance only; no RNG or forcing sequence is invented.
The default boundary description retains actual periodic/PEC faces and the
separate scalar-conductivity layer. This is not a PML archive.
"""
function archive_run(
    path::AbstractString,
    model::YeeConductiveLayer;
    normalization::AbstractString = "epsilon0=1;mu0=1/c^2;scalar_matched_conductivity",
    boundary_desc::Union{Nothing,AbstractString} = nothing,
    kwargs...,
)
    boundary_desc === nothing && (
        boundary_desc =
            join(map(p -> p ? "periodic" : "pec", model.inner.g.periodic), ",")::String *
            ";scalar_conductivity"
    )
    _yee_layer_acquire!(model)
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
