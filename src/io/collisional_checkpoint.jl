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

function _collisional_checkpoint1_state(st::CollisionalHybridStepper, species::AbstractVector)
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

function _validate_collisional_checkpoint1(s, st::CollisionalHybridStepper{D,T}) where {D,T}
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
    _collisional_model_matches(s.inner.model, st.model) && isequal(s.inner.shape, st.shape) ||
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
        isequal(s.inner.last_dt, st.dt) && st.dt <= s.source_time
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
    if _evolves_pressure(st.model.closure)
        # The full-field restore follows the canonical inner restore. It must
        # not overwrite the pressure that the inner validator just admitted.
        isequal(s.fields.pe, s.inner.pe) ||
            throw(ArgumentError("collisional checkpoint inner and public pressure disagree"))
    end
    if _transports_electrons(st.model.closure)
        # The inner validator owns policy/type/range checks. The wrapper stores
        # additional public fields and particles, which must not share force storage.
        for force in s.inner.Eforce
            for components in values(s.fields),
                array in (components isa Tuple ? components : (components,))

                _particle_arrays_mightalias(force, array) &&
                    throw(ArgumentError("collisional checkpoint force must not alias saved fields"))
            end
            for ps in s.inner.species, array in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
                _particle_arrays_mightalias(force, array) && throw(
                    ArgumentError("collisional checkpoint force must not alias saved particles"),
                )
            end
        end
    end
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
    previous_report = st.thermal_report
    st.phase = :inflight
    completed = false
    try
        _restore_validated_hybrid_checkpoint!(st.inner, populations, saved.inner, count)
        _restore_checkpoint_fields!(_hybrid_field_state(st.inner.fields), saved.fields)
        st.inner.fields.floor_count[] = saved.floor_count
        copy!(st.rng, saved.rng)
        st.source_time = saved.source_time
        st.source_count = saved.source_count
        if st.electron_ion !== nothing
            st.thermal_report = saved.thermal_report
            st.thermal_work.pending = nothing
        end
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
            st.thermal_report = previous_report
        end
        st.phase = completed ? :ready : previous_phase
    end
end

const _COLLISIONAL_THERMAL_CHECKPOINT_FIELDS =
    (_COLLISIONAL_CHECKPOINT_FIELDS..., :electron_ion, :thermal_report)

const _COLLISIONAL_CGL_POLICY_FIELDS = (:cgl_policy_version, :cgl_force_floor, :cgl_force_operator)
const _COLLISIONAL_CGL_CHECKPOINT_FIELDS =
    (_COLLISIONAL_CHECKPOINT_FIELDS..., _COLLISIONAL_CGL_POLICY_FIELDS...)
_collisional_cgl_policy(::Type{T}) where {T} = (
    cgl_policy_version = 1,
    cgl_force_floor = T(1e-12),
    cgl_force_operator = :global_invariants_trial_B_tensor,
)

function _checkpoint_state(st::CollisionalHybridStepper, species::AbstractVector)
    common = _collisional_checkpoint1_state(st, species)
    if st.model.closure isa CGLElectrons
        return merge(common, (schema = 4,), _collisional_cgl_policy(eltype(st.inner.fields.n)))
    end
    st.electron_ion === nothing && return common
    _validate_electron_ion_report(
        st.thermal_report,
        st.electron_ion,
        st.inner,
        st.species,
        st.inner.fields,
        st.source_time,
        st.source_count,
        st.dt,
        st.epsilon0,
        st.weight_to_number,
    )
    # The payload is subsequently validated by the common transport. Detach the
    # receipt now so that capture never needs to waive incoming/live links.
    return merge(
        common,
        (
            schema = st.electron_ion.policy == 1 ? 2 : 3,
            electron_ion = _electron_ion_configuration(st.electron_ion),
            thermal_report = deepcopy(st.thermal_report),
        ),
    )
end

function _validate_checkpoint_state(s, st::CollisionalHybridStepper)
    if st.model.closure isa CGLElectrons
        return _validate_collisional_cgl_checkpoint(s, st)
    end
    if s isa NamedTuple && any(k -> haskey(s, k), _COLLISIONAL_CGL_POLICY_FIELDS)
        throw(ArgumentError("non-CGL collisional state cannot carry CGL schema fields"))
    end
    if st.electron_ion !== nothing
        return _validate_collisional_thermal_checkpoint(s, st)
    end
    if s isa NamedTuple
        any(k -> haskey(s, k), (:electron_ion, :thermal_report)) &&
            throw(ArgumentError("schema-1 collisional state cannot carry thermal schema fields"))
    end
    return _validate_collisional_checkpoint1(s, st)
end

function _validate_collisional_cgl_checkpoint(s, st::CollisionalHybridStepper{D,T}) where {D,T}
    _electron_ion_record(s, _COLLISIONAL_CGL_CHECKPOINT_FIELDS)
    s.schema isa Int && s.schema == 4 ||
        throw(ArgumentError("collisional CGL requires checkpoint schema 4"))
    st.electron_ion === nothing && st.thermal_work === nothing && st.thermal_report === nothing ||
        throw(ArgumentError("collisional CGL cannot carry scalar thermal exchange state"))
    _validate_collisional_model(st.inner)
    expected = _collisional_cgl_policy(T)
    s.cgl_policy_version isa Int &&
        s.cgl_policy_version === expected.cgl_policy_version &&
        s.cgl_force_floor isa T &&
        isequal(s.cgl_force_floor, expected.cgl_force_floor) &&
        s.cgl_force_operator === expected.cgl_force_operator ||
        throw(ArgumentError("collisional CGL force policy must match the destination"))
    legacy = merge(NamedTuple{_COLLISIONAL_CHECKPOINT_FIELDS}(s), (schema = 1,))
    _validate_collisional_checkpoint1(legacy, st)
    _validate_collisional_cgl_coefficients(s.inner.model.closure)
    return nothing
end

# Known closed rational value fields use the canonical integer value codec;
# their private BigInt pointers are not stored process-local resources. All
# unrelated container/metadata edges retain the canonical graph policy.
function _validate_thermal_checkpoint_storage(s, st)
    owned = IdDict{Any,Bool}()
    arrays = AbstractArray[]
    function remember(value)
        (ismutable(value) || value isa AbstractArray) || return nothing
        haskey(owned, value) && return nothing
        owned[value] = true
        value isa AbstractArray && push!(arrays, value)
        return nothing
    end
    legacy = NamedTuple{_COLLISIONAL_CHECKPOINT_FIELDS}(s)
    _register_checkpoint_configuration_owned!(remember, remember, owned, legacy, st.inner)
    _visit_checkpoint_objects(
        remember,
        (st.transaction, st.workspace, st.thermal_work, st.electron_ion, st.rng, st.thermal_report),
    )
    graph = IdDict{Any,NamedTuple{(:logical, :edges),Tuple{Bool,Vector{Tuple{Any,Bool}}}}}()
    function record(value, logical, edges)
        value isa Union{_ElectronIonQ,BigInt,Regex} && empty!(edges)
        graph[value] = (logical = logical, edges = edges)
        return nothing
    end
    function independent(value)
        haskey(owned, value) && throw(
            ArgumentError(
                "incoming thermal receipt must not link to saved or live simulation storage",
            ),
        )
        if value isa AbstractArray
            any(a -> _particle_arrays_mightalias(value, a), arrays) &&
                throw(ArgumentError("incoming thermal receipt must not alias saved or live arrays"))
        end
        return nothing
    end
    _visit_checkpoint_objects(independent, s.thermal_report; record)
    _validate_checkpoint_configuration_graph(s.thermal_report, graph, independent)
    # Enumerate logical mutable roles rather than letting the graph's visited
    # set erase a repeated cell/ion-array role.
    if st.electron_ion.policy == 2
        _validate_electron_ion_overlap_roles(s.thermal_report)
    elseif s.thermal_report isa NamedTuple && haskey(s.thermal_report, :cells)
        cells = s.thermal_report.cells
        cells isa Vector{NamedTuple} ||
            throw(ArgumentError("thermal receipt cells must be a standard Vector"))
        roles = AbstractArray[cells]
        for cell in cells
            _electron_ion_record(cell, _ELECTRON_ION_CELL_KEYS)
            cell.ions isa Vector{NamedTuple} ||
                throw(ArgumentError("thermal receipt ions must be a standard Vector"))
            for other in roles
                _particle_arrays_mightalias(cell.ions, other) &&
                    throw(ArgumentError("thermal receipt array roles must be independent"))
            end
            push!(roles, cell.ions)
        end
    end
    return nothing
end

function _validate_collisional_thermal_checkpoint(s, st::CollisionalHybridStepper{D,T}) where {D,T}
    st.electron_ion.policy == 2 && return _validate_collisional_overlap_checkpoint(s, st)
    _electron_ion_record(s, _COLLISIONAL_THERMAL_CHECKPOINT_FIELDS)
    s.schema isa Int && s.schema == 2 ||
        throw(ArgumentError("enabled electron-ion exchange requires collisional schema 2"))
    expected = _electron_ion_configuration(st.electron_ion)
    typeof(s.electron_ion) === typeof(expected) && isequal(s.electron_ion, expected) || throw(
        ArgumentError("electron-ion checkpoint configuration must match the bound destination"),
    )
    legacy = merge(NamedTuple{_COLLISIONAL_CHECKPOINT_FIELDS}(s), (schema = 1,))
    _validate_collisional_checkpoint1(legacy, st)
    _validate_thermal_checkpoint_storage(s, st)
    # Borrow validated incoming arrays only for private endpoint inspection;
    # public destination values and its old phase do not enter this report.
    populations =
        [ParticleSet{D,T}(p.x, p.v, p.weight, p.id, p.tag, p.q, p.m) for p in s.inner.species]
    _validate_electron_ion_report(
        s.thermal_report,
        st.electron_ion,
        st.inner,
        populations,
        s.fields,
        s.source_time,
        s.source_count,
        s.dt,
        s.epsilon0,
        s.weight_to_number,
    )
    return nothing
end

function _validate_collisional_overlap_checkpoint(s, st::CollisionalHybridStepper{D,T}) where {D,T}
    _electron_ion_record(s, _COLLISIONAL_THERMAL_CHECKPOINT_FIELDS)
    s.schema isa Int && s.schema == 3 ||
        throw(ArgumentError("CIC/TSC thermal overlap requires collisional schema 3"))
    expected = _electron_ion_configuration(st.electron_ion)
    typeof(s.electron_ion) === typeof(expected) && isequal(s.electron_ion, expected) || throw(
        ArgumentError("thermal overlap checkpoint configuration must match the bound destination"),
    )
    legacy = merge(NamedTuple{_COLLISIONAL_CHECKPOINT_FIELDS}(s), (schema = 1,))
    _validate_collisional_checkpoint1(legacy, st)
    _validate_thermal_checkpoint_storage(s, st)
    populations =
        [ParticleSet{D,T}(p.x, p.v, p.weight, p.id, p.tag, p.q, p.m) for p in s.inner.species]
    _validate_electron_ion_report(
        s.thermal_report,
        st.electron_ion,
        st.inner,
        populations,
        s.fields,
        s.source_time,
        s.source_count,
        s.dt,
        s.epsilon0,
        s.weight_to_number,
    )
    return nothing
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

"Verify the archive checksum and use the existing locked stochastic restore."
function load_archive!(st::CollisionalHybridStepper, species::AbstractVector, path::AbstractString)
    _acquire_collisional!(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, species, loaded.state)
    finally
        st.busy[] = false
    end
end
