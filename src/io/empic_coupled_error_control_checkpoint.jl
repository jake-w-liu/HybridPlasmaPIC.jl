# Schema4 preserves the full-time controller and its own historical admission.
# The shared serializer/archive owners remain the only byte-format owners.
const _EMPIC_COUPLED_CHECKPOINT_FIELDS =
    (:coupled_control,:coupled_proposal,:coupled_phase,:spectral_signature)
function _reject_empic_coupled_checkpoint_fields(s)
    _validate_checkpoint_container(s,(:schema_version,))
    (s.schema_version===4 || any(k -> haskey(s,k),_EMPIC_COUPLED_CHECKPOINT_FIELDS) ||
        (haskey(s,:phase_policy) && s.phase_policy===:full_time_coupled)) &&
        throw(ArgumentError("coupled full-time checkpoint requires a coupled destination without phase migration"))
    return nothing
end
function _normalize_empic_coupled_checkpoint(s)
    _validate_checkpoint_container(s,(:schema_version,:schedule_policy,:phase_policy,_EMPIC_COUPLED_CHECKPOINT_FIELDS...))
    s.schema_version===4 && s.schedule_policy===3 && s.phase_policy===:full_time_coupled ||
        throw(ArgumentError("coupled full-time checkpoints require schema4/policy3 without phase migration"))
    return s
end
function _ec_checkpoint_schedule(st)
    schedule=_empic_collection(st).schedule
    schedule isa _EMPICCoupledSchedule || throw(ArgumentError("construct a coupled full-time checkpoint destination"))
    return schedule
end
_ec_checkpoint_phase_wire(phase) = NamedTuple{_EC_PHASE_KEYS}(Tuple(getfield(phase,k) for k in _EC_PHASE_KEYS))
function _ec_checkpoint_phase(wire,::Type{T}) where {T}
    wire isa NamedTuple && keys(wire)==_EC_PHASE_KEYS || throw(ArgumentError("coupled saved phase has noncanonical fields"))
    all(_ec_q,(wire.time,wire.elapsed,wire.initial_energy,wire.field_defect,wire.particle_defect,wire.coupling_defect)) &&
        wire.interval isa Vector{T} && wire.last isa Vector{_EMPICCoupledQ} && wire.next isa Vector{_EMPICCoupledQ} &&
        wire.next_h isa T && (wire.report===nothing || wire.report isa NamedTuple) ||
        throw(ArgumentError("coupled saved phase requires exact canonical scalar/vector types"))
    _ec_identical(wire.admission_control,_ec_control(wire.admission_control,T)) ||
        throw(ArgumentError("coupled saved historical admission control is not canonical"))
    _ec_costs_check(wire.costs)
    _ec_closed(wire)
    return _EMPICCoupledPhase{T}(Tuple(getproperty(wire,k) for k in _EC_PHASE_KEYS)...)
end
function _ec_checkpoint_context(st,g,s,populations)
    p=s.coupled_proposal
    Model=st isa EMPIC ? EMPIC : EMPIC1D
    context=Model(g,populations;mobile=s.mobile,background_charge=s.n0,
        neutrality=s.neutrality,charge_rtol=s.charge_rtol,c=s.c,shape=s.shape,
        n_sub=s.n_sub,relativistic=s.relativistic,species_dt=p.cap,
        species_adaptive=p.adaptive,species_dt_min=p.minimum,schedule_limits=p.limits,
        event_budget=p.budget,error_control=s.coupled_control)
    _ec_restore_fields!(context,s.fields)
    _ec_sync_magnetic!(context)
    context.time[],context.step[]=s.time,s.step
    context.last_dt[],context.last_sub_dt[]=s.last_dt,s.last_sub_dt
    return context
end
# This helper is called only while the destination owner is held. It validates
# the ORIGINAL graph before any constructor/copy can hide repeated roles.
function _ec_prepare_checkpoint(s,st,g::FourierGrid{D,T};destination=nothing) where {D,T}
    _ec_checkpoint_schedule(st);_ec_mode(T)
    g===st.g || throw(ArgumentError("coupled checkpoint must use its destination grid"))
    _normalize_empic_coupled_checkpoint(s)
    required=(_EMPIC_COLLECTION_CHECKPOINT_FIELDS...,(st isa EMPIC1D ? (:polarization,) : ())...,
        :schedule_policy,:phase_policy,_EMPIC_COUPLED_CHECKPOINT_FIELDS...)
    keys(s)==required && s.integrator===_empic_collection_checkpoint_integrator(st) &&
        s.particle_layout===:collection && (st isa EMPIC || s.polarization===:ex_ey_bz) &&
        s.D isa Int && s.T===T && s.ncell isa typeof(g.n) && s.L isa typeof(g.L) &&
        s.mobile isa BitVector && s.species isa Vector && all(i -> isassigned(s.species,i),eachindex(s.species)) ||
        throw(ArgumentError("coupled checkpoint envelope has noncanonical fields or types"))
    for record in s.species
        record isa NamedTuple && keys(record)==_CHECKPOINT_PARTICLE_FIELDS ||
            throw(ArgumentError("coupled saved particle record has noncanonical fields"))
    end
    s.fields isa NamedTuple && keys(s.fields)==keys(_empic_event_fields(st)) ||
        throw(ArgumentError("coupled saved field inventory differs"))
    _ec_independent(_ec_role_arrays((s.mobile,s.species,s.fields,s.coupled_control,s.coupled_proposal,s.coupled_phase)))
    _empic_event_checkpoint_graph(s)
    destination===nothing || _empic_event_checkpoint_live_links(s,st,destination)
    _ec_closed((s.species,s.fields,s.coupled_control,s.coupled_proposal,s.coupled_phase))
    _empic_event_checkpoint_proposal_layout(s.coupled_proposal,T)
    control=_ec_control(s.coupled_control,T)
    _ec_identical(control,s.coupled_control) || throw(ArgumentError("coupled saved future control is not canonical"))
    p=s.coupled_proposal
    proposal=_empic_event_proposal(s.mobile,p.cap,p.adaptive,p.minimum,p.limits,p.budget,T)
    _ec_identical(proposal,p) || throw(ArgumentError("coupled saved nominal proposal is not canonical"))
    s.n_sub isa Int && 1<=s.n_sub<=proposal.budget || throw(ArgumentError("coupled saved base/event budget differs"))
    phase=_ec_checkpoint_phase(s.coupled_phase,T)
    s.spectral_signature isa String && s.spectral_signature==_ec_spectral_signature(g) ||
        throw(ArgumentError("coupled saved spectral operator signature differs"))
    populations=_validate_empic_collection_checkpoint_envelope(s,st,g,4)
    # The constructor only allocates/adopts saved scalar policy; no init/solve/
    # projection is invoked. Incoming primaries are copied verbatim to context.
    context=_ec_checkpoint_context(st,g,s,populations)
    populations=_ec_snapshot_check(context,(species=s.species,fields=s.fields))
    context.collection.schedule.phase=phase
    context.collection.initialized[]=s.initialized
    clock=(s.time,s.step,s.last_dt,s.last_sub_dt)
    _ec_identical(phase.next_h,_ec_future_hint(phase,control)) ||
        throw(ArgumentError("coupled saved next interval differs from its historical recommendation and future bounds"))
    if !s.initialized
        iszero(s.step) && all(iszero,(phase.time,phase.elapsed,phase.initial_energy,phase.field_defect,
            phase.particle_defect,phase.coupling_defect,phase.interval...,phase.last...,phase.next...)) &&
            all(iszero,values(phase.costs)) && phase.report===nothing ||
            throw(ArgumentError("uninitialized coupled seeds must retain an empty zero-time history"))
    end
    meter=_ec_meter(phase.admission_control)
    _ec_validate_phase(context,phase,populations,s.fields,s.initialized,clock,meter)
    return (context=context,populations=populations,proposal=proposal,control=control,phase=phase)
end
function _validate_empic_coupled_checkpoint(s,st,g;destination=nothing)
    schedule=_ec_checkpoint_schedule(st);_ec_acquire(schedule)
    try
        _ec_prepare_checkpoint(s,st,g;destination)
        return nothing
    finally
        _ec_release(schedule)
    end
end
function _empic_coupled_checkpoint_state_locked(st,species)
    schedule=_ec_checkpoint_schedule(st);work=st.collection
    T=eltype(st.rho_n);D=length(st.g.n)
    _empic_collection_validate_storage(st,species)
    fields=_empic_event_fields(st);wire=_ec_checkpoint_phase_wire(schedule.phase)
    # Inspect actual flags (including packed chunks), public views and every
    # private phase/report role before the common builder copies any leaves.
    _ec_independent(_ec_role_arrays((work.mobile,species,fields,schedule.control,schedule.proposal,wire)))
    graph=(mobile=work.mobile,species=[_checkpoint_particle_state(ps) for ps in species],
        fields=fields,control=schedule.control,proposal=schedule.proposal,phase=wire)
    _empic_event_checkpoint_graph(graph)
    _ec_closed((fields,schedule.control,schedule.proposal,wire))
    if work.initialized[]
        _ec_ready(st,species,_ec_meter(schedule.phase.admission_control))
    else
        _empic_collection_validate_species(st,species)
        _ec_phase_scalar_check(st,schedule.phase,length(species),false,
            (st.time[],st.step[],st.last_dt[],st.last_sub_dt[]))
    end
    ordinary=_ec_copy_species(species,Val(D),T)
    common=_empic_collection_checkpoint_common(st,species,st.g)
    original=merge(common,(schema_version=4,schedule_policy=3,phase_policy=:full_time_coupled,
        species=NamedTuple[_checkpoint_particle_state(ps) for ps in ordinary],fields=fields,
        coupled_control=schedule.control,coupled_proposal=schedule.proposal,coupled_phase=wire,
        spectral_signature=schedule.signature))
    saved=deepcopy(original)
    _ec_prepare_checkpoint(saved,st,st.g)
    return saved
end
function _empic_coupled_checkpoint_state(st,species)
    schedule=_ec_checkpoint_schedule(st);_ec_acquire(schedule)
    try
        return _empic_coupled_checkpoint_state_locked(st,species)
    finally
        _ec_release(schedule)
    end
end
function _ec_checkpoint_copy_stage!(destination,prepared)
    # Types, shapes and ordinary independent stores were admitted before commit.
    foreach(copyto!,_empic_transaction_fields(destination),_empic_transaction_fields(prepared))
    destination.c,destination.shape,destination.relativistic,destination.n_sub=
        prepared.c,prepared.shape,prepared.relativistic,prepared.n_sub
    destination.time[],destination.step[]=prepared.time[],prepared.step[]
    destination.last_dt[],destination.last_sub_dt[]=prepared.last_dt[],prepared.last_sub_dt[]
    return nothing
end
function _restore_empic_coupled_checkpoint_locked!(st,species,s,g::FourierGrid{D,T}) where {D,T}
    work=st.collection;schedule=_ec_checkpoint_schedule(st)
    # This inspects original aliases/live links before deepcopy. Stale numerical
    # destination values are deliberately not passed through _ec_ready.
    _empic_collection_validate_storage(st,species)
    _ec_independent(_ec_role_arrays((work.mobile,species)))
    _ec_prepare_checkpoint(s,st,g;destination=species)
    length(species)==length(s.species) && all(ps -> ps isa eltype(work.owners),species) ||
        throw(ArgumentError("coupled checkpoint destination population layout differs"))
    for (ps,record) in zip(species,s.species)
        _validate_checkpoint_destination(ps,length(record.weight))
    end
    saved=deepcopy(s)
    admitted=_ec_prepare_checkpoint(saved,st,g)
    context,phase,proposal,control=admitted.context,admitted.phase,admitted.proposal,admitted.control
    candidate=_ec_zero_phase(T,length(species),control)
    clock=(saved.time,saved.step,saved.last_dt,saved.last_sub_dt)
    shadow=saved.initialized ? _ec_shadow(context,admitted.populations,saved.fields,species,clock) : nothing
    arrays=vcat(_ec_role_arrays((proposal,_ec_checkpoint_phase_wire(phase),_ec_checkpoint_phase_wire(candidate))),
        _empic_event_shadow_arrays(shadow))
    _ec_independent(arrays)
    for array in arrays,ps in species
        _mightalias_particle_storage(array,ps) && throw(ArgumentError("coupled restored phase/report/shadow aliases public particles"))
    end
    # Reserve all permitted growth before the first ordinary value/length store.
    for (ps,record) in zip(species,saved.species)
        count=length(record.weight)
        for array in (ps.x...,ps.v...,ps.weight,ps.id,ps.tag)
            if length(array)!=count
                array isa Vector || throw(ArgumentError("only Vector coupled checkpoint destinations may resize"))
                sizehint!(array,count)
            end
        end
    end
    sizehint!(work.owners,length(species))
    # Public/stage objects and all their arrays retain identity. The independent
    # prepared context already carries full B in Bc/Bzc and reset scratch rows.
    for (ps,record) in zip(species,saved.species)
        _restore_checkpoint_particles!(ps,record)
    end
    _ec_checkpoint_copy_stage!(st,context)
    _ec_checkpoint_copy_stage!(work.stage,context)
    schedule.proposal,schedule.control=proposal,control
    schedule.phase,schedule.candidate=phase,candidate
    schedule.shadow,schedule.candidate_shadow=shadow,nothing
    schedule.signature=saved.spectral_signature
    empty!(work.owners);append!(work.owners,species)
    work.initialized[]=saved.initialized
    return st
end
function _restore_empic_coupled_checkpoint!(st,species,s,g)
    schedule=_ec_checkpoint_schedule(st);_ec_acquire(schedule)
    try
        return _restore_empic_coupled_checkpoint_locked!(st,species,s,g)
    finally
        _ec_release(schedule)
    end
end
