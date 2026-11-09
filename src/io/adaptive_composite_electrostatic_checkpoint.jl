const _ACE_PHASE_FIELDS = (
    :time,
    :step,
    :last_dt,
    :epoch_start_time,
    :epoch_start_step,
    :regrid_count,
    :initial_energy,
    :field_defect,
    :kick_defect,
    :trajectory_defect,
    :physical_defect,
    :mesh_work,
    :last_regrid_report,
)
const _ACE_CHECKPOINT_FIELDS = (
    :integrator,
    :schema_version,
    :configuration,
    :phase,
    :initialized,
    :current,
    _ACE_PHASE_FIELDS...,
    :history,
)
const _ACE_PHYSICAL_FIELDS = (:refinements, :species, :phi, :last_solve_report, :last_step_report)
const _ACE_PERIODIC_PHYSICAL_FIELDS = (_ACE_PHYSICAL_FIELDS..., :background_density, :transport)
_ace_checkpoint_kind(st) =
    st.mesh.periodic ? :periodic_adaptive_composite_es : :adaptive_composite_es
const _ACE_HISTORY_FIELDS =
    (:epoch_start_time, :epoch_start_step, :global_time, :global_step, :last_dt, :state)
function _ace_keys(record, names)
    _validate_checkpoint_container(record, names)
    keys(record) == names ||
        throw(ArgumentError("adaptive composite checkpoint has noncanonical or unsupported fields"))
    return nothing
end
function _ace_particle_record(ps)
    return (
        x = map(copy, ps.x),
        v = map(copy, ps.v),
        weight = copy(ps.weight),
        id = copy(ps.id),
        tag = copy(ps.tag),
        q = ps.q,
        m = ps.m,
    )
end
function _ace_physical_payload(state)
    core = (
        refinements = state.mesh.refinements,
        species = [_ace_particle_record(ps) for ps in state.species],
        phi = copy(state.phi),
        last_solve_report = deepcopy(state.solve_report),
        last_step_report = deepcopy(state.step_report),
    )
    return state.mesh.periodic ?
           merge(
        core,
        (
            background_density = deepcopy(state.background_density),
            transport = deepcopy(state.transport),
        ),
    ) : core
end
# Structural role inventory deliberately ignores numerical values: old destination
# particles/fields/phase may be stale when a valid incoming record replaces them.
function _ace_payload_arrays(record, st::AdaptiveCompositeElectrostaticPIC{D,T}) where {D,T}
    _ace_keys(record, st.mesh.periodic ? _ACE_PERIODIC_PHYSICAL_FIELDS : _ACE_PHYSICAL_FIELDS)
    record.refinements isa Tuple &&
        record.species isa Vector &&
        length(record.species) == length(st._factors) &&
        record.phi isa Vector{T} ||
        throw(ArgumentError("adaptive composite physical payload layout differs"))
    arrays = Any[record.phi]
    for ps in record.species
        _ace_keys(ps, _CHECKPOINT_PARTICLE_FIELDS)
        ps.x isa NTuple{D,Vector{T}} &&
            ps.v isa NTuple{3,Vector{T}} &&
            ps.weight isa Vector{T} &&
            ps.id isa Vector{UInt64} &&
            ps.tag isa Vector{UInt32} &&
            ps.q isa T &&
            ps.m isa T || throw(
            ArgumentError(
                "adaptive composite particle record needs ordinary exact-type vectors/scalars",
            ),
        )
        n = length(ps.weight)
        for array in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            length(array) == n && axes(array, 1) == Base.OneTo(n) ||
                throw(DimensionMismatch("adaptive composite saved particle component sizes differ"))
            push!(arrays, array)
        end
    end
    record.last_solve_report isa AdaptiveCompositeESSolveReport{T} ||
        throw(ArgumentError("adaptive composite saved solve receipt type differs"))
    if record.last_step_report !== nothing
        record.last_step_report isa AdaptiveCompositeESStepReport{D,T} ||
            throw(ArgumentError("adaptive composite saved step receipt type differs"))
        push!(arrays, record.last_step_report.crossings)
    end
    if st.mesh.periodic
        record.background_density isa _CompositeQ ||
            throw(ArgumentError("periodic physical background density requires original Q"))
        if record.transport !== nothing
            record.transport isa _PeriodicCompositeTransport{D,T} ||
                throw(ArgumentError("periodic saved transport type differs"))
            append!(arrays, _pacc_transport_arrays(record.transport))
        end
    end
    return arrays
end
function _ace_history_arrays(history, st)
    _ace_keys(history, _ACE_HISTORY_FIELDS)
    return _ace_payload_arrays(history.state, st)
end
function _ace_prepare_payload(st::AdaptiveCompositeElectrostaticPIC{D,T}, record) where {D,T}
    _ace_independent(_ace_payload_arrays(record, st))
    all(a -> a isa Tuple{Int,NTuple{D,Int}}, record.refinements) ||
        throw(ArgumentError("adaptive composite saved cell addresses require canonical Int tuples"))
    mesh = CartesianAMRHierarchy(
        st.mesh.nbase,
        st.mesh.L,
        record.refinements;
        periodic = st.mesh.periodic,
    )
    isequal(record.refinements, mesh.refinements) ||
        throw(ArgumentError("adaptive composite saved refinement order is not canonical"))
    populations = _ElectrostaticHostParticle{D,T}[]
    for particle in record.species
        push!(populations, _checkpoint_species(particle, Val(D), T))
    end
    _ace_values(mesh, populations, st._factors)
    private = _ace_copy(populations, Val(D), T)
    state = _AdaptiveCompositePhysical{D,T}(
        mesh,
        _acq_data(mesh),
        private,
        copy(record.phi),
        _acc_bins(private, mesh),
        deepcopy(record.last_solve_report),
        deepcopy(record.last_step_report),
        st.mesh.periodic ? deepcopy(record.transport) : nothing,
        st.mesh.periodic ? deepcopy(record.background_density) : nothing,
    )
    _ace_physical_check(st, state)
    return state
end
function _ace_validate_history(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    history,
    current,
    phase,
) where {D,T}
    history isa NamedTuple ||
        throw(ArgumentError("adaptive composite empty epoch lost its preceding physical payload"))
    _ace_history_arrays(history, st)
    et, en, t, n, h = history.epoch_start_time,
    history.epoch_start_step,
    history.global_time,
    history.global_step,
    history.last_dt
    et isa T &&
        t isa T &&
        h isa T &&
        en isa Int &&
        n isa Int &&
        isfinite(et) &&
        isfinite(t) &&
        isfinite(h) &&
        _stored_nonnegative(et) &&
        _stored_positive(t) &&
        _stored_positive(h) &&
        0 <= en < n &&
        et < t &&
        h <= t &&
        isequal(t, phase.time) &&
        n == phase.step &&
        isequal(h, phase.last_dt) &&
        (iszero(en) ? _particle_stored_zero(et) : _stored_positive(et) && h <= et) ||
        throw(ArgumentError("adaptive composite retained history epoch/global clock differs"))
    previous = _ace_prepare_payload(st, history.state)
    _ace_step_report(st, previous, t, n, h, et, en)
    _ace_species_equal(previous.species, current.species) || throw(
        ArgumentError(
            "adaptive composite retained ending particles differ from unchanged current particles",
        ),
    )
    return nothing
end
function _ace_checkpoint_data(saved, st::AdaptiveCompositeElectrostaticPIC)
    _ace_storage(st, st._owners; counts = false)
    _ace_keys(saved, _ACE_CHECKPOINT_FIELDS)
    saved.integrator === _ace_checkpoint_kind(st) &&
        saved.schema_version === 1 &&
        saved.phase === :synchronous_integer &&
        saved.initialized === true &&
        typeof(saved.configuration) === typeof(st._configuration) &&
        isequal(saved.configuration, st._configuration) || throw(
        ArgumentError(
            "unsupported adaptive composite checkpoint family/schema/phase/configuration",
        ),
    )
    incoming = _ace_payload_arrays(saved.current, st)
    saved.history === nothing || append!(incoming, _ace_history_arrays(saved.history, st))
    _ace_independent(incoming)
    live_arrays = _ace_live_arrays(st, st._owners)
    for array in incoming, live in live_arrays
        _particle_arrays_mightalias(array, live) && throw(
            ArgumentError(
                "adaptive composite incoming checkpoint aliases live/private/geometry/history storage",
            ),
        )
    end
    current = _ace_prepare_payload(st, saved.current)
    phase =
        NamedTuple{_ACE_PHASE_FIELDS}(Tuple(getproperty(saved, name) for name in _ACE_PHASE_FIELDS))
    _ace_phase_check(st, current, saved.history, phase)
    return current
end
function _validate_checkpoint_state(saved, st::AdaptiveCompositeElectrostaticPIC)
    _ace_checkpoint_data(saved, st)
    return nothing
end
function _ace_snapshot_locked(st, species)
    _ace_ready(st, species)
    saved = deepcopy(
        merge(
            (
                integrator = _ace_checkpoint_kind(st),
                schema_version = 1,
                configuration = st._configuration,
                phase = :synchronous_integer,
                initialized = true,
                current = _ace_physical_payload(st._state),
            ),
            _ace_phase(st),
            (history = st._history,),
        ),
    )
    _ace_checkpoint_data(saved, st)
    return saved
end
function _checkpoint_state(st::AdaptiveCompositeElectrostaticPIC, species::AbstractVector)
    _ace_acquire(st)
    try
        return _ace_snapshot_locked(st, species)
    finally
        _ace_release(st)
    end
end
function _ace_restore_locked!(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    species,
    saved,
) where {D,T}
    _ace_storage(st, species; counts = false)
    # The original graph is inspected before deepcopy can conceal shared roles.
    _ace_checkpoint_data(saved, st)
    detached = deepcopy(saved)
    current = _ace_checkpoint_data(detached, st)
    history = detached.history
    phase = NamedTuple{_ACE_PHASE_FIELDS}(
        Tuple(getproperty(detached, name) for name in _ACE_PHASE_FIELDS),
    )
    prepared = _ace_copy(current.species, Val(D), T)
    public_phi = copy(current.phi)
    owners = _ElectrostaticHostParticle{D,T}[ps for ps in species]
    components = Tuple[_sce_particle_arrays(ps) for ps in prepared]
    counts = Int[length(ps.weight) for ps in prepared]
    for (destination, replacement) in zip(species, prepared)
        typeof(destination) === typeof(replacement) ||
            throw(ArgumentError("adaptive composite restore particle replacement type differs"))
    end
    _ace_prepared(st, current, prepared, public_phi, history)
    for (destination, replacement) in zip(species, prepared)
        _publish_particle_replacement!(destination, replacement)
        destination.q, destination.m = replacement.q, replacement.m
    end
    return _ace_publish!(st, current, public_phi, owners, components, counts, history, phase)
end
function _restore_checkpoint_state!(
    st::AdaptiveCompositeElectrostaticPIC,
    species::AbstractVector,
    saved,
)
    _ace_acquire(st)
    try
        return _ace_restore_locked!(st, species, saved)
    finally
        _ace_release(st)
    end
end
"Atomically save the complete adaptive field/particles/global phase, old-mesh history and work ledger."
function save_checkpoint(
    path::AbstractString,
    st::AdaptiveCompositeElectrostaticPIC,
    species::AbstractVector,
)
    _ace_acquire(st)
    try
        return _save_checkpoint(path, st, _ace_snapshot_locked(st, species))
    finally
        _ace_release(st)
    end
end
function load_checkpoint!(
    st::AdaptiveCompositeElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    _ace_acquire(st)
    try
        return _ace_restore_locked!(st, species, deserialize(path))
    finally
        _ace_release(st)
    end
end
function archive_run(
    path::AbstractString,
    st::AdaptiveCompositeElectrostaticPIC,
    species::AbstractVector;
    boundary_desc::AbstractString = st.mesh.periodic ? "periodic Cartesian composite torus" :
                                    "grounded Cartesian composite box; particle-bound cloud tails",
    normalization::AbstractString = "epsilon0=1",
    kwargs...,
)
    _ace_acquire(st)
    try
        return _archive_run(
            path,
            st,
            _ace_snapshot_locked(st, species);
            boundary_desc,
            normalization,
            kwargs...,
        )
    finally
        _ace_release(st)
    end
end
function load_archive!(
    st::AdaptiveCompositeElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    _ace_acquire(st)
    try
        record = load_run(path)
        return _ace_restore_locked!(st, species, record.state)
    finally
        _ace_release(st)
    end
end
