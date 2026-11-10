# The dynamic payload owns changing topology/global accounting and delegates
# each complete physical epoch to the existing static schema/validation owner.
const _DYNAMIC_COMPOSITE_ES_FIELDS = (
    :integrator,
    :schema_version,
    :configuration,
    :initialized,
    :current,
    :time,
    :step,
    :last_dt,
    :epoch_start_time,
    :epoch_start_step,
    :regrid_count,
    :initial_energy,
    :physical_defect,
    :mesh_work,
    :last_regrid_report,
    :history,
)
const _DCE_HISTORY_FIELDS =
    (:epoch_start_time, :epoch_start_step, :global_time, :global_step, :state)

function _dce_static_arrays(s, ::Type{T}) where {T}
    _sce_checkpoint_keys(s, _STATIC_COMPOSITE_ES_CHECKPOINT_FIELDS)
    s.E isa Vector{T} && s.phi isa Vector{T} && s.last_impulse isa Vector{_CompositeQ} ||
        throw(ArgumentError("dynamic epoch fields/history require ordinary exact-type Vectors"))
    s.species isa Vector ||
        throw(ArgumentError("dynamic epoch requires ordinary ordered particle records"))
    populations = _ElectrostaticHostParticle{1,T}[]
    arrays = Any[s.E, s.phi, s.last_impulse]
    for record in s.species
        _sce_checkpoint_keys(record, _CHECKPOINT_PARTICLE_FIELDS)
        ps = _checkpoint_species(record, Val(1), T)
        _sce_population_storage([ps], T)
        push!(populations, ps)
        append!(arrays, _sce_particle_arrays(ps))
    end
    if s.last_report !== nothing
        s.last_report isa StaticCompositeESReport{T} ||
            throw(ArgumentError("dynamic epoch has an unsupported physical step report"))
        push!(arrays, s.last_report.interface_crossings)
    end
    return arrays, populations
end

function _dce_history_arrays(history, ::Type{T}) where {T}
    _sce_checkpoint_keys(history, _DCE_HISTORY_FIELDS)
    arrays, _ = _dce_static_arrays(history.state, T)
    return arrays
end

function _dce_epoch_destination(st::DynamicCompositeElectrostaticPIC{T}, saved) where {T}
    _, populations = _dce_static_arrays(saved, T)
    expected = _sce_configuration(st._state)
    _sce_checkpoint_keys(saved.configuration, keys(expected))
    canonical = merge(saved.configuration, (refinements = (),))
    wanted = st._configuration.canonical
    typeof(canonical) === typeof(wanted) && isequal(canonical, wanted) ||
        throw(ArgumentError("dynamic epoch domain/source/factor or numerical policy differs"))
    saved.configuration.refinements isa Tuple ||
        throw(ArgumentError("dynamic epoch refinements must be an immutable canonical tuple"))
    mesh = StaticAMRHierarchy1D(st.mesh.nbase, st.mesh.L, saved.configuration.refinements)
    private = _dce_copy_species(populations, T)
    destination = StaticCompositeElectrostaticPIC(
        mesh,
        private;
        shape = st.shape,
        background_charge = st.background_charge,
        left_plate_charge = st.left_plate_charge,
    )
    # This checks the ORIGINAL static record graph and represented source map,
    # using independent destination particles/geometry rather than its own arrays.
    _sce_checkpoint_data(saved, destination)
    return destination
end

function _dce_prepare_epoch(st::DynamicCompositeElectrostaticPIC, saved)
    destination = _dce_epoch_destination(st, saved)
    _restore_checkpoint_state!(destination, _dce_private_species(destination), saved)
    return destination
end

function _dce_clone_epoch(st::DynamicCompositeElectrostaticPIC)
    saved = _checkpoint_state(st._state, _dce_private_species(st._state))
    return _dce_prepare_epoch(st, saved)
end

function _dce_validate_history(
    history,
    st::DynamicCompositeElectrostaticPIC{T},
    current,
    time,
    step,
    dt,
) where {T}
    history === nothing &&
        throw(ArgumentError("retained dynamic physical history is missing"))
    _sce_checkpoint_keys(history, _DCE_HISTORY_FIELDS)
    history.epoch_start_time isa T &&
        history.epoch_start_step isa Int &&
        history.global_time isa T &&
        history.global_step isa Int &&
        history.global_step == step &&
        isequal(history.global_time, time) ||
        throw(ArgumentError("retained dynamic physical step has invalid global/epoch provenance"))
    destination = _dce_epoch_destination(st, history.state)
    saved = history.state
    saved.step > 0 && isequal(saved.last_dt, dt) ||
        throw(ArgumentError("retained dynamic physical history must be a real fixed-duration step"))
    _dce_epoch_clock(
        history.global_time,
        history.global_step,
        dt,
        history.epoch_start_time,
        history.epoch_start_step,
        saved.time,
        saved.step,
        saved.last_dt,
    )
    _, populations = _dce_static_arrays(saved, T)
    _dce_species_equal(populations, _dce_private_species(current)) || throw(
        ArgumentError("regridding changed the retained physical-step particle inventory/state"),
    )
    # The saved old E/phi/current/report were validated against THEIR OWN mesh and
    # detached saved particles by _sce_checkpoint_data, never current new leaves.
    destination.mesh.nbase == current.mesh.nbase && isequal(destination.mesh.L, current.mesh.L) ||
        throw(ArgumentError("retained dynamic physical step has a different physical domain"))
    return nothing
end

function _dce_prepared_roles(
    st::DynamicCompositeElectrostaticPIC{T},
    current,
    species,
    history,
) where {T}
    _sce_population_storage(species, T)
    length(species) == length(st._factors) &&
        isequal(_dce_configuration(current), st._configuration) ||
        throw(ArgumentError("prepared dynamic state changes immutable configuration"))
    private = _dce_private_species(current)
    _sce_ready(current, private)
    for (actual, prepared, factor) in zip(species, private, st._factors)
        isequal((actual.q, actual.m), factor) &&
            length(actual.weight) == length(prepared.weight) &&
            isequal(actual.weight, prepared.weight) &&
            isequal(actual.id, prepared.id) &&
            isequal(actual.tag, prepared.tag) ||
            throw(ArgumentError("prepared dynamic state changes particle sources/provenance"))
    end
    writable, readonly = _sce_arrays(current, private)
    arrays = vcat(writable, readonly)
    for ps in species
        append!(arrays, _sce_particle_arrays(ps))
    end
    history === nothing || append!(arrays, _dce_history_arrays(history, T))
    _dce_independent(arrays)
    return nothing
end

function _dce_checkpoint_data(saved, st::DynamicCompositeElectrostaticPIC{T}) where {T}
    # Old physical scalar/particle values may be stale. The bound ordinary owner
    # graph, view identities, geometry caches and fixed configuration remain real.
    _dce_storage(st, st._owners; counts = false, values = false)
    _sce_checkpoint_keys(saved, _DYNAMIC_COMPOSITE_ES_FIELDS)
    saved.integrator === :dynamic_composite_es &&
        saved.schema_version === 1 &&
        saved.initialized === true &&
        typeof(saved.configuration) === typeof(st._configuration) &&
        isequal(saved.configuration, st._configuration) ||
        throw(ArgumentError("unsupported dynamic composite checkpoint family/schema/configuration"))
    incoming, _ = _dce_static_arrays(saved.current, T)
    saved.history === nothing || append!(incoming, _dce_history_arrays(saved.history, T))
    _dce_independent(incoming)
    live = _dce_owned_arrays(st, st._owners)
    for array in incoming, destination in live
        _particle_arrays_mightalias(array, destination) && throw(
            ArgumentError("dynamic checkpoint must be detached from live/private/geometry storage"),
        )
    end
    current = _dce_prepare_epoch(st, saved.current)
    _dce_phase(
        st,
        current,
        saved.history,
        saved.time,
        saved.step,
        saved.last_dt,
        saved.epoch_start_time,
        saved.epoch_start_step,
        saved.regrid_count,
        saved.initial_energy,
        saved.physical_defect,
        saved.mesh_work,
        saved.last_regrid_report,
    )
    return current
end

function _validate_checkpoint_state(saved, st::DynamicCompositeElectrostaticPIC)
    _dce_checkpoint_data(saved, st)
    return nothing
end

function _dce_snapshot_locked(st::DynamicCompositeElectrostaticPIC, species)
    _dce_ready(st, species)
    record = (
        integrator = :dynamic_composite_es,
        schema_version = 1,
        configuration = st._configuration,
        initialized = true,
        current = _checkpoint_state(st._state, _dce_private_species(st._state)),
        time = st.time,
        step = st.step,
        last_dt = st.last_dt,
        epoch_start_time = st.epoch_start_time,
        epoch_start_step = st.epoch_start_step,
        regrid_count = st.regrid_count,
        initial_energy = st.initial_energy,
        physical_defect = st.physical_defect,
        mesh_work = st.mesh_work,
        last_regrid_report = st.last_regrid_report,
        history = st._history,
    )
    saved = deepcopy(record)
    _validate_checkpoint_state(saved, st)
    return saved
end

function _checkpoint_state(st::DynamicCompositeElectrostaticPIC, species::AbstractVector)
    _dce_acquire(st)
    try
        return _dce_snapshot_locked(st, species)
    finally
        _dce_release(st)
    end
end

function _dce_restore_locked!(st::DynamicCompositeElectrostaticPIC{T}, species, record) where {T}
    _dce_storage(st, species; counts = false, values = false)
    # Fully check original links before deepcopy can hide an alias to old state.
    _dce_checkpoint_data(record, st)
    saved = deepcopy(record)
    current = _dce_checkpoint_data(saved, st)
    private = _dce_private_species(current)
    replacements = _dce_copy_species(private, T)
    owners = _ElectrostaticHostParticle{1,T}[ps for ps in species]
    components = Tuple[_sce_particle_arrays(ps) for ps in replacements]
    counts = Int[length(ps.weight) for ps in replacements]
    for (destination, prepared) in zip(species, replacements)
        typeof(destination) === typeof(prepared) ||
            throw(ArgumentError("dynamic checkpoint replacement particle storage types differ"))
    end
    _dce_prepared_roles(st, current, replacements, saved.history)
    for (destination, prepared) in zip(species, replacements)
        _publish_particle_replacement!(destination, prepared)
        destination.q, destination.m = prepared.q, prepared.m
    end
    return _dce_publish!(
        st,
        current,
        owners,
        components,
        counts,
        saved.history,
        saved.time,
        saved.step,
        saved.last_dt,
        saved.epoch_start_time,
        saved.epoch_start_step,
        saved.regrid_count,
        saved.initial_energy,
        saved.physical_defect,
        saved.mesh_work,
        saved.last_regrid_report,
    )
end

function _restore_checkpoint_state!(
    st::DynamicCompositeElectrostaticPIC,
    species::AbstractVector,
    saved,
)
    _dce_acquire(st)
    try
        return _dce_restore_locked!(st, species, saved)
    finally
        _dce_release(st)
    end
end

"Atomically save every dynamic physical epoch, prior-step mesh/history and exact ledger."
function save_checkpoint(
    path::AbstractString,
    st::DynamicCompositeElectrostaticPIC,
    species::AbstractVector,
)
    _dce_acquire(st)
    try
        return _save_checkpoint(path, st, _dce_snapshot_locked(st, species))
    finally
        _dce_release(st)
    end
end

"Restore the saved topology, particles, fields, global/epoch phase and ledgers verbatim."
function load_checkpoint!(
    st::DynamicCompositeElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    _dce_acquire(st)
    try
        return _dce_restore_locked!(st, species, deserialize(path))
    finally
        _dce_release(st)
    end
end

function archive_run(
    path::AbstractString,
    st::DynamicCompositeElectrostaticPIC,
    species::AbstractVector;
    boundary_desc::AbstractString = "closed isolated dynamic composite 1D",
    kwargs...,
)
    _dce_acquire(st)
    try
        return _archive_run(path, st, _dce_snapshot_locked(st, species); boundary_desc, kwargs...)
    finally
        _dce_release(st)
    end
end

function load_archive!(
    st::DynamicCompositeElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    _dce_acquire(st)
    try
        loaded = load_run(path)
        return _dce_restore_locked!(st, species, loaded.state)
    finally
        _dce_release(st)
    end
end
