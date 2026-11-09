# Collection payloads extend the existing checkpoint/archive owners. Workspaces,
# gathers and FFT scratch remain recomputable and are not serialized.
const _EMPIC_COLLECTION_CHECKPOINT_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :n0,
    :c,
    :shape,
    :n_sub,
    :particle_layout,
    :schema_version,
    :relativistic,
    :mobile,
    :neutrality,
    :charge_rtol,
    :initialized,
    :species,
    :fields,
    :time,
    :step,
    :last_dt,
    :last_sub_dt,
)

const _EMPIC_SPECIES_SCHEDULE_CHECKPOINT_FIELDS =
    (:schedule_policy, :species_stride, :last_species_dt)
const _EMPIC_EVENT_CHECKPOINT_FIELDS = (:phase_policy, :event_proposal, :event_phase)

_empic_collection_checkpoint_integrator(::EMPIC) = :empic_collection
_empic_collection_checkpoint_integrator(::EMPIC1D) = :empic1d_collection

_empic_collection_checkpoint_arrays(::EMPIC, fields) =
    (fields.E..., fields.B..., fields.J..., fields.rho_n, fields.rho_np1)
_empic_collection_checkpoint_arrays(::EMPIC1D, fields) =
    (fields.Ex, fields.Ey, fields.Bz, fields.Jx, fields.Jy, fields.rho_n, fields.rho_np1)

function _validate_empic_collection_checkpoint(s, st::_ElectromagneticCheckpointStepper)
    normalized = _normalize_empic_collection_checkpoint(s, st)
    return _validate_normalized_empic_collection_checkpoint(normalized, st, st.g)
end

# Schema3 is deliberately closed: none of its pending source phase may be
# discarded by a legacy destination or inferred from a dividing history.
function _reject_empic_event_checkpoint_fields(s)
    _validate_checkpoint_container(s, (:schema_version,))
    _reject_empic_coupled_checkpoint_fields(s)
    any(k -> haskey(s, k), _EMPIC_EVENT_CHECKPOINT_FIELDS) &&
        throw(ArgumentError("persistent event checkpoint requires an event-schedule destination"))
    return nothing
end

function _normalize_empic_event_checkpoint(s)
    _reject_empic_coupled_checkpoint_fields(s)
    _validate_checkpoint_container(
        s,
        (:schema_version, :schedule_policy, _EMPIC_EVENT_CHECKPOINT_FIELDS...),
    )
    s.schema_version === 3 && s.schedule_policy === 2 && s.phase_policy === :persistent_event ||
        throw(ArgumentError("persistent schedules require schema3/policy2 without phase migration"))
    return s
end

function _empic_event_checkpoint_phase(wire, ::Type{T}) where {T}
    wire isa NamedTuple &&
        keys(wire) == (:time, :primed, :interval, :last, :next, :report) &&
        _empic_event_q(wire.time) &&
        wire.primed isa Bool &&
        wire.interval isa Vector{T} &&
        wire.last isa Vector{_EMPICEventQ} &&
        wire.next isa Vector{_EMPICEventQ} ||
        throw(ArgumentError("checkpoint persistent phase layout differs"))
    return _EMPICEventPhase{T}(
        wire.time,
        wire.primed,
        wire.interval,
        wire.last,
        wire.next,
        wire.report,
    )
end

function _empic_event_checkpoint_proposal_layout(p, ::Type{T}) where {T}
    p isa NamedTuple &&
        keys(p) == (:cap, :adaptive, :minimum, :limits, :budget) &&
        p.cap isa Vector{T} &&
        p.adaptive isa BitVector &&
        p.minimum isa Vector{T} &&
        p.budget isa Int &&
        p.limits isa NamedTuple &&
        keys(p.limits) == (:rotation, :electric, :travel, :growth) &&
        all(v -> v isa T, values(p.limits)) ||
        throw(ArgumentError("checkpoint persistent proposal layout differs"))
    return nothing
end

function _empic_event_checkpoint_arrays(s, populations, phase)
    p = s.event_proposal
    # Keep repeated semantic roles: a seen-object walk would erase a malformed
    # repeated leaf before the pairwise independence check could observe it.
    arrays = AbstractArray[
        s.mobile,
        s.species,
        p.cap,
        p.adaptive,
        p.minimum,
        phase.interval,
        phase.last,
        phase.next,
    ]
    for value in values(s.fields)
        append!(arrays, value isa Tuple ? value : (value,))
    end
    for ps in populations
        append!(arrays, _empic_event_arrays(ps))
    end
    append!(arrays, _empic_event_report_arrays(phase.report))
    for (i, a) in enumerate(arrays), j = 1:(i-1)
        b = arrays[j]
        left, right = _checkpoint_packed_chunks(a), _checkpoint_packed_chunks(b)
        # Packed flags also own their chunks; this catches a saved UInt64
        # particle-ID role retaining the same buffer across Julia versions.
        aliased =
            _particle_arrays_mightalias(a, b) ||
            (left !== nothing && _particle_arrays_mightalias(left, b)) ||
            (right !== nothing && _particle_arrays_mightalias(a, right)) ||
            (left !== nothing && right !== nothing && _particle_arrays_mightalias(left, right))
        aliased &&
            throw(ArgumentError("checkpoint persistent array/container roles must be independent"))
    end
    return arrays
end

function _empic_event_checkpoint_live_links(s, st, species)
    live = IdDict{Any,Bool}()
    arrays = AbstractArray[]
    function remember(value)
        # Closed Q fields are numeric value codecs; shared BigInt scalar values
        # neither own a second array role nor retain an executable resource.
        value isa Union{_EMPICEventQ,BigInt} && return nothing
        # Immutable value records (control/proposal tuples and phase carriers)
        # cannot retain executable resources on their own. Their mutable fields
        # are visited separately, so only mutable objects and arrays are live.
        (ismutable(value) || value isa AbstractArray) || return nothing
        live[value] = true
        value isa AbstractArray && push!(arrays, value)
        return nothing
    end
    _visit_checkpoint_objects(remember, (st, species); backing = remember)
    function independent(value)
        value isa Union{_EMPICEventQ,BigInt} && return nothing
        if ismutable(value) || value isa AbstractArray
            haskey(live, value) && throw(
                ArgumentError(
                    "checkpoint persistent source must not retain live solver objects ($(typeof(value)))",
                ),
            )
            if value isa AbstractArray
                any(a -> _particle_arrays_mightalias(value, a), arrays) && throw(
                    ArgumentError("checkpoint persistent source must not alias live solver storage"),
                )
            end
        end
        return nothing
    end
    _visit_checkpoint_objects(independent, s; backing = independent)
    return nothing
end

function _empic_event_checkpoint_graph(s)
    graph = IdDict{Any,NamedTuple{(:logical, :edges),Tuple{Bool,Vector{Tuple{Any,Bool}}}}}()
    function record(value, logical, edges)
        # Q is a closed numeric receipt. Its BigInts are values, not retained
        # native limbs or identity-bearing configuration. Other wrapper and
        # metadata edges keep the canonical serialization/resource policy.
        value isa Union{_EMPICEventQ,BigInt,Regex} && empty!(edges)
        graph[value] = (logical = logical, edges = edges)
        return nothing
    end
    _visit_checkpoint_objects(_ -> nothing, s; record)
    _validate_checkpoint_configuration_graph(s, graph, _ -> nothing)
    return nothing
end

function _validate_empic_event_checkpoint(
    s,
    st,
    g::FourierGrid{D,T};
    destination = nothing,
) where {D,T}
    _normalize_empic_event_checkpoint(s)
    required = (
        _EMPIC_COLLECTION_CHECKPOINT_FIELDS...,
        :schedule_policy,
        _EMPIC_EVENT_CHECKPOINT_FIELDS...,
        (st isa EMPIC1D ? (:polarization,) : ())...,
    )
    _validate_checkpoint_container(s, required)
    length(keys(s)) == length(required) &&
        s.D isa Int &&
        s.T === T &&
        s.ncell isa typeof(g.n) &&
        s.L isa typeof(g.L) &&
        s.mobile isa BitVector &&
        s.species isa Vector &&
        all(i -> isassigned(s.species, i), eachindex(s.species)) ||
        throw(ArgumentError("checkpoint persistent envelope must have the closed schema3 layout"))
    for record in s.species
        record isa NamedTuple && keys(record) == _CHECKPOINT_PARTICLE_FIELDS ||
            throw(ArgumentError("checkpoint persistent particle record layout differs"))
    end
    populations = _validate_empic_collection_checkpoint_envelope(s, st, g, 3)
    keys(s.fields) == keys(_empic_event_fields(st)) ||
        throw(ArgumentError("checkpoint persistent field inventory differs"))
    _validate_checkpoint_fields(s.fields, _empic_event_fields(st))
    for value in values(s.fields), array in (value isa Tuple ? value : (value,))
        array isa Array{T,D} && all(isfinite, array) ||
            throw(ArgumentError("checkpoint persistent fields require finite ordinary grid arrays"))
    end
    _empic_event_checkpoint_proposal_layout(s.event_proposal, T)
    phase = _empic_event_checkpoint_phase(s.event_phase, T)
    phase.report === nothing || _empic_event_report_layout(st, phase.report)
    _empic_event_checkpoint_arrays(s, populations, phase)
    _empic_event_checkpoint_graph(s)
    # Capture has its own live references by design. Restore must inspect the
    # original source graph before any copying or proposal normalization.
    destination === nothing || _empic_event_checkpoint_live_links(s, st, destination)
    p = s.event_proposal
    proposal = _empic_event_proposal(s.mobile, p.cap, p.adaptive, p.minimum, p.limits, p.budget, T)
    s.n_sub <= proposal.budget ||
        throw(ArgumentError("checkpoint base field partition exceeds its saved event budget"))
    candidate = _empic_event_zero_phase(T, length(s.mobile))
    prepared = _EMPICEventSchedule{T}(proposal, phase, candidate, nothing, nothing)
    _empic_collection_validate_schedule(st, prepared)
    _empic_event_phase_values(
        st,
        phase,
        populations,
        s.initialized,
        s.time,
        s.step,
        s.last_dt,
        s.last_sub_dt;
        speed = s.c,
        fields = s.fields,
    )
    return nothing
end

function _empic_event_checkpoint_state(st, common)
    schedule = _empic_collection(st).schedule
    phase = schedule.phase
    state = merge(
        common,
        (
            schema_version = 3,
            schedule_policy = 2,
            phase_policy = :persistent_event,
            fields = _empic_event_fields(st),
            event_proposal = schedule.proposal,
            event_phase = (
                time = phase.time,
                primed = phase.primed,
                interval = phase.interval,
                last = phase.last,
                next = phase.next,
                report = phase.report,
            ),
        ),
    )
    _validate_empic_event_checkpoint(state, st, st.g)
    saved = deepcopy(state)
    _validate_empic_event_checkpoint(saved, st, st.g)
    return saved
end

function _restore_empic_event_checkpoint!(st, species, s, g::FourierGrid{D,T}) where {D,T}
    work = _empic_collection(st)
    _validate_empic_event_checkpoint(s, st, g; destination = species)
    _empic_collection_validate_storage(st, species)
    length(species) == length(s.species) ||
        throw(ArgumentError("checkpoint population count does not match destination"))
    all(ps -> ps isa eltype(work.owners), species) || throw(
        ArgumentError(
            "population storage type changed; construct a new electromagnetic collection",
        ),
    )
    for (ps, record) in zip(species, s.species)
        _validate_checkpoint_destination(ps, length(record.weight))
    end
    saved = deepcopy(s)
    _validate_empic_event_checkpoint(saved, st, g)
    populations =
        ParticleSet{D,T}[_checkpoint_species(record, Val(D), T) for record in saved.species]
    p = saved.event_proposal
    proposal =
        _empic_event_proposal(saved.mobile, p.cap, p.adaptive, p.minimum, p.limits, p.budget, T)
    phase = _empic_event_checkpoint_phase(saved.event_phase, T)
    candidate = _empic_event_zero_phase(T, length(saved.mobile))
    configuration = _empic_event_configuration(
        st;
        c = saved.c,
        shape = saved.shape,
        relativistic = saved.relativistic,
    )
    clock = (saved.time, saved.step, saved.last_dt, saved.last_sub_dt)
    shadow =
        saved.initialized ?
        _empic_event_shadow_values(st, populations, saved.fields, species, configuration, clock) :
        nothing
    # Capacity and every new authoritative/private object are prepared before
    # ordinary stores, with the same resource/concurrency boundary as old IO.
    for (ps, record) in zip(species, saved.species)
        count = length(record.weight)
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            if length(a) != count
                a isa Vector ||
                    throw(ArgumentError("only Vector checkpoint destinations may resize"))
                sizehint!(a, count)
            end
        end
    end
    sizehint!(work.owners, length(species))
    for (ps, record) in zip(species, saved.species)
        _restore_checkpoint_particles!(ps, record)
    end
    _restore_checkpoint_fields!(_empic_event_fields(st), saved.fields)
    st.c, st.shape, st.n_sub, st.relativistic =
        saved.c, saved.shape, saved.n_sub, saved.relativistic
    st.time[], st.step[] = saved.time, saved.step
    st.last_dt[], st.last_sub_dt[] = saved.last_dt, saved.last_sub_dt
    schedule = work.schedule
    schedule.proposal, schedule.phase, schedule.candidate = proposal, phase, candidate
    schedule.shadow, schedule.candidate_shadow = shadow, nothing
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = saved.initialized
    return st
end

"Verify the archive checksum and restore an existing electromagnetic collection verbatim."
function load_archive!(
    st::_ElectromagneticCheckpointStepper,
    species::AbstractVector{<:ParticleSet},
    path::AbstractString,
)
    _empic_collection(st)
    loaded = load_run(path)
    return _restore_checkpoint_state!(st, species, loaded.state)
end

_empic_collection_checkpoint_schema(st) =
    _empic_collection(st).schedule isa _EMPICCoupledSchedule ? 4 :
    _empic_collection(st).schedule isa _EMPICEventSchedule ? 3 :
    (_empic_collection(st).schedule === nothing ? 1 : 2)

function _normalize_empic_collection_checkpoint(s, st::EMPIC)
    _empic_collection(st).schedule isa _EMPICCoupledSchedule &&
        return _normalize_empic_coupled_checkpoint(s)
    _empic_collection(st).schedule isa _EMPICEventSchedule &&
        return _normalize_empic_event_checkpoint(s)
    _reject_empic_event_checkpoint_fields(s)
    _validate_checkpoint_container(s, (:schema_version,))
    if _empic_collection(st).schedule === nothing
        any(k -> haskey(s, k), _EMPIC_SPECIES_SCHEDULE_CHECKPOINT_FIELDS) && throw(
            ArgumentError(
                "an unscheduled generic checkpoint cannot contain species schedule fields",
            ),
        )
        return s
    end
    return _normalize_empic_species_schedule_checkpoint(s, st)
end

_normalize_empic_collection_checkpoint(s, st::EMPIC1D) =
    _normalize_empic_species_schedule_checkpoint(s, st)

function _normalize_empic_species_schedule_checkpoint(s, st)
    _empic_collection(st).schedule isa _EMPICCoupledSchedule &&
        return _normalize_empic_coupled_checkpoint(s)
    _empic_collection(st).schedule isa _EMPICEventSchedule &&
        return _normalize_empic_event_checkpoint(s)
    _reject_empic_event_checkpoint_fields(s)
    _validate_checkpoint_container(s, (:schema_version,))
    if s.schema_version === 1
        any(k -> haskey(s, k), _EMPIC_SPECIES_SCHEDULE_CHECKPOINT_FIELDS) && throw(
            ArgumentError("schema-1 electromagnetic checkpoint cannot contain schedule fields"),
        )
        _validate_empic_collection_checkpoint_envelope(s, st, st.g, 1)
        T = typeof(st.n0)
        stride = Int[mobile ? 1 : 0 for mobile in s.mobile]
        history = T[s.step > 0 && mobile ? s.last_sub_dt : zero(T) for mobile in s.mobile]
        return merge(
            s,
            (
                schema_version = 2,
                schedule_policy = 1,
                species_stride = stride,
                last_species_dt = history,
            ),
        )
    end
    return s
end

function _validate_normalized_empic_collection_checkpoint(
    s,
    st::_ElectromagneticCheckpointStepper,
    g::FourierGrid{D,T},
) where {D,T}
    _empic_collection(st).schedule isa _EMPICCoupledSchedule &&
        return _validate_empic_coupled_checkpoint(s, st, g)
    _empic_collection(st).schedule isa _EMPICEventSchedule &&
        return _validate_empic_event_checkpoint(s, st, g)
    version = _empic_collection_checkpoint_schema(st)
    populations = _validate_empic_collection_checkpoint_envelope(s, st, g, version)
    if version == 2
        _validate_checkpoint_container(s, _EMPIC_SPECIES_SCHEDULE_CHECKPOINT_FIELDS)
        s.schedule_policy === 1 ||
            throw(ArgumentError("unsupported electromagnetic checkpoint schedule policy"))
        _empic_validate_species_schedule(s.mobile, s.n_sub, s.species_stride)
        _empic_validate_species_history(
            s.mobile,
            s.last_species_dt,
            s.step,
            s.last_dt,
            s.last_sub_dt,
            T,
        )
        metadata = (s.species_stride, s.last_species_dt)
        fields = _empic_collection_checkpoint_arrays(st, s.fields)
        for (i, a) in enumerate(metadata)
            for b in (fields..., metadata[1:(i-1)]...)
                Base.mightalias(a, b) && throw(
                    ArgumentError(
                        "checkpoint schedule storage must be independent of fields and history",
                    ),
                )
            end
            for ps in populations
                _mightalias_particle_storage(a, ps) &&
                    throw(ArgumentError("checkpoint schedule storage must not alias particles"))
            end
        end
    end
    return nothing
end

function _validate_empic_collection_checkpoint_envelope(
    s,
    st::_ElectromagneticCheckpointStepper,
    ::FourierGrid{D,T},
    schema_version::Int,
) where {D,T}
    _validate_checkpoint_container(s, _EMPIC_COLLECTION_CHECKPOINT_FIELDS)
    if st isa EMPIC1D
        _validate_checkpoint_container(s, (:polarization,))
        s.polarization === :ex_ey_bz ||
            throw(ArgumentError("unsupported EMPIC1D collection checkpoint polarization"))
    end
    _validate_checkpoint_grid(s, st.g)
    s.integrator === _empic_collection_checkpoint_integrator(st) &&
        s.particle_layout === :collection &&
        s.schema_version === schema_version ||
        throw(ArgumentError("unsupported electromagnetic collection checkpoint schema"))
    work = _empic_collection(st)
    s.n0 isa T && isequal(s.n0, st.n0) ||
        throw(ArgumentError("checkpoint background_charge does not match destination"))
    s.neutrality === work.neutrality ||
        throw(ArgumentError("checkpoint neutrality policy does not match destination"))
    s.charge_rtol isa T && isequal(s.charge_rtol, work.charge_rtol) ||
        throw(ArgumentError("checkpoint charge_rtol does not match destination"))
    s.mobile isa AbstractVector{Bool} &&
        axes(s.mobile, 1) == Base.OneTo(length(work.mobile)) &&
        s.mobile == work.mobile ||
        throw(ArgumentError("checkpoint per-population mobility does not match destination"))
    s.relativistic isa Bool || throw(ArgumentError("checkpoint relativistic flag must be Bool"))
    typeof(s.shape) === typeof(st.shape) && s.shape isa Union{NGP,CIC,TSC} ||
        throw(ArgumentError("checkpoint shape does not match destination"))
    s.c isa T || throw(ArgumentError("checkpoint c precision does not match destination"))
    _require_finite_positive_real("checkpoint c", s.c, T)
    s.n_sub isa Int && s.n_sub >= 1 ||
        throw(ArgumentError("checkpoint n_sub must be a positive Int"))
    s.initialized isa Bool || throw(ArgumentError("checkpoint initialized flag must be Bool"))
    _validate_checkpoint_time_state(s, T)
    s.last_sub_dt isa T &&
        isfinite(s.last_sub_dt) &&
        zero(T) <= s.last_sub_dt <= s.last_dt &&
        (s.step == 0 ? iszero(s.last_sub_dt) : s.last_sub_dt > zero(T)) ||
        throw(ArgumentError("checkpoint fast timestep history is inconsistent"))
    _validate_empic_clock_phase(s.time, s.step, s.last_dt, s.last_sub_dt)
    !s.initialized &&
        s.step != 0 &&
        throw(ArgumentError("advanced electromagnetic checkpoint must be initialized"))
    s.species isa AbstractVector && axes(s.species, 1) == Base.OneTo(length(s.species)) ||
        throw(ArgumentError("checkpoint species must be an ordered one-based vector"))
    populations = [_checkpoint_species(record, Val(D), T) for record in s.species]
    # Concrete records give a concretely typed vector. An empty collection of
    # untyped records infers Vector{Any}; keep the typed empty vector the
    # collection validators require, as the original typed comprehension did.
    eltype(populations) <: ParticleSet{D,T} || (populations = ParticleSet{D,T}[])
    _empic_collection_validate_populations(st, populations, s.c, s.relativistic)
    _validate_checkpoint_fields(s.fields, _empic_checkpoint_fields(st))
    arrays = _empic_collection_checkpoint_arrays(st, s.fields)
    for (i, a) in enumerate(arrays)
        all(isfinite, a) || throw(ArgumentError("checkpoint electromagnetic fields must be finite"))
        for j = 1:(i-1)
            Base.mightalias(a, arrays[j]) &&
                throw(ArgumentError("checkpoint electromagnetic fields must be independent"))
        end
        for ps in populations
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("checkpoint particles must not alias fields"))
        end
    end
    return populations
end

_checkpoint_state(st::_ElectromagneticCheckpointStepper, species::AbstractVector{<:ParticleSet}) =
    _empic_collection_checkpoint_state(st, species, st.g)

function _empic_collection_checkpoint_common(st, species, ::FourierGrid{D,T}) where {D,T}
    work = _empic_collection(st)
    state = (
        integrator = _empic_collection_checkpoint_integrator(st),
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        n0 = st.n0,
        c = st.c,
        shape = st.shape,
        n_sub = st.n_sub,
        particle_layout = :collection,
        schema_version = 1,
        relativistic = st.relativistic,
        mobile = copy(work.mobile),
        neutrality = work.neutrality,
        charge_rtol = work.charge_rtol,
        initialized = work.initialized[],
        species = [_checkpoint_particle_state(ps) for ps in species],
        fields = _empic_checkpoint_fields(st),
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        last_sub_dt = st.last_sub_dt[],
    )
    st isa EMPIC1D && (state = merge(state, (polarization = :ex_ey_bz,)))
    return state
end

function _empic_collection_checkpoint_state(
    st::_ElectromagneticCheckpointStepper,
    species::AbstractVector{<:ParticleSet},
    ::FourierGrid{D,T},
) where {D,T}
    _empic_collection(st).schedule isa _EMPICCoupledSchedule &&
        return _empic_coupled_checkpoint_state(st, species)
    work = _empic_collection(st)
    _empic_collection_validate_species(st, species)
    _empic_collection_validate_phase(st)
    work.initialized[] && _empic_collection_bound(work, species)
    state = _empic_collection_checkpoint_common(st, species, st.g)
    if work.schedule isa _EMPICEventSchedule
        return _empic_event_checkpoint_state(st, state)
    end
    if work.schedule !== nothing
        return merge(
            state,
            (
                schema_version = 2,
                schedule_policy = 1,
                species_stride = copy(work.schedule.stride),
                last_species_dt = copy(work.schedule.last_species_dt),
            ),
        )
    end
    return state
end

_restore_checkpoint_state!(
    st::_ElectromagneticCheckpointStepper,
    species::AbstractVector{<:ParticleSet},
    s,
) = _restore_empic_collection_checkpoint!(st, species, s, st.g)

function _restore_empic_collection_checkpoint!(
    st::_ElectromagneticCheckpointStepper,
    species::AbstractVector{<:ParticleSet},
    s,
    ::FourierGrid{D,T},
) where {D,T}
    _empic_collection(st).schedule isa _EMPICCoupledSchedule &&
        return _restore_empic_coupled_checkpoint!(st, species, s, st.g)
    _empic_collection(st).schedule isa _EMPICEventSchedule &&
        return _restore_empic_event_checkpoint!(st, species, s, st.g)
    work = _empic_collection(st)
    s = _normalize_empic_collection_checkpoint(s, st)
    _validate_normalized_empic_collection_checkpoint(s, st, st.g)
    _empic_collection_validate_storage(st, species)
    length(species) == length(s.species) ||
        throw(ArgumentError("checkpoint population count does not match destination"))
    all(ps -> ps isa eltype(work.owners), species) || throw(
        ArgumentError(
            "population storage type changed; construct a new electromagnetic collection state",
        ),
    )
    saved_populations = [_checkpoint_species(record, Val(D), T) for record in s.species]
    saved_fields = _empic_collection_checkpoint_arrays(st, s.fields)
    for (ps, record) in zip(species, s.species)
        _validate_checkpoint_destination(ps, length(record.weight))
        for saved in saved_populations
            _particle_sets_mightalias(ps, saved) && throw(
                ArgumentError("checkpoint source and destination particles must be independent"),
            )
        end
        for a in saved_fields
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("checkpoint fields must not alias destination particles"))
        end
    end
    # As in the collisional checkpoint owner, detach every incoming leaf before
    # publication. In-memory payload fields can be a permutation of live fields;
    # copying one directly would overwrite a later source leaf.
    saved = deepcopy(s)
    _validate_normalized_empic_collection_checkpoint(saved, st, st.g)
    stride =
        work.schedule === nothing ? nothing :
        _empic_species_stride_values(saved.mobile, saved.n_sub, saved.species_stride)
    history = work.schedule === nothing ? nothing : T[value for value in saved.last_species_dt]
    # Reserve ordinary Vector capacity before any values or lengths publish.
    # Fixed strided views have already passed the exact-count preflight above.
    for (ps, record) in zip(species, saved.species)
        count = length(record.weight)
        for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            if length(a) != count
                a isa Vector ||
                    throw(ArgumentError("only Vector checkpoint destinations may resize"))
                sizehint!(a, count)
            end
        end
    end
    sizehint!(work.owners, length(species))
    for (ps, record) in zip(species, saved.species)
        _restore_checkpoint_particles!(ps, record)
    end
    _restore_checkpoint_fields!(_empic_checkpoint_fields(st), saved.fields)
    st.c, st.shape, st.n_sub = saved.c, saved.shape, saved.n_sub
    if work.schedule !== nothing
        work.schedule.stride = stride
        copyto!(work.schedule.last_species_dt, history)
    end
    st.relativistic = saved.relativistic
    st.time[], st.step[] = saved.time, saved.step
    st.last_dt[], st.last_sub_dt[] = saved.last_dt, saved.last_sub_dt
    empty!(work.owners)
    append!(work.owners, species)
    work.initialized[] = saved.initialized
    return st
end
