# This family has no projection/phase migration from spectral or closed models.
const _COMPATIBLE_YEE_CHECKPOINT_FIELDS = (
    :integrator,
    :schema,
    :configuration,
    :phase,
    :initialized,
    :species,
    :fields,
    :current,
    :charge_n,
    :charge_np1,
    :time,
    :step,
    :last_dt,
    :last_report,
)
function _pyee_checkpoint_keys(s, required)
    _validate_checkpoint_container(s, required)
    length(keys(s)) == length(required) ||
        throw(ArgumentError("compatible Yee checkpoint contains unsupported fields"))
    return nothing
end
function _pyee_checkpoint_data(s, st::ImplicitYeeEMPICStepper{T}) where {T}
    _pyee_binding(st, st.owners; values = false, counts = false, sources = false)
    _pyee_checkpoint_keys(s, _COMPATIBLE_YEE_CHECKPOINT_FIELDS)
    s.integrator === :implicit_periodic_yee_cloud_midpoint &&
        s.schema === 1 &&
        s.phase === :synchronous_integer &&
        s.initialized === true ||
        throw(ArgumentError("unsupported compatible Yee family/schema/phase"))
    configuration = _pyee_configuration(
        st.g,
        st.shape,
        st.c,
        st._species_scalars,
        st.controls,
        st.constraint_tolerances,
    )
    typeof(s.configuration) === typeof(configuration) && isequal(s.configuration, configuration) ||
        throw(ArgumentError("compatible Yee fixed grid/width/physics/solver policies must match"))
    _pyee_phase(st, s.initialized, s.time, s.step, s.last_dt, s.last_report)
    s.species isa Vector && length(s.species) == length(st.owners) ||
        throw(ArgumentError("compatible Yee checkpoint requires its ordered species count"))
    populations = _ElectrostaticHostParticle{1,T}[]
    for (index, record) in enumerate(s.species)
        _pyee_checkpoint_keys(record, _CHECKPOINT_PARTICLE_FIELDS)
        ps = _checkpoint_species(record, Val(1), T)
        ps isa _ElectrostaticHostParticle{1,T} &&
            isequal((ps.q, ps.m), st._species_scalars[index]) || throw(
            ArgumentError(
                "compatible Yee checkpoint requires ordinary matching populations/original factors",
            ),
        )
        push!(populations, ps)
    end
    _pyee_particles(st, populations)
    _pyee_checkpoint_keys(s.fields, (:E, :B))
    _validate_checkpoint_fields(s.fields, (E = st.fields.E, B = st.fields.B))
    all(a -> a isa Vector{T} && all(isfinite, a), (s.fields.E..., s.fields.B...)) || throw(
        ArgumentError("compatible Yee checkpoint fields must be finite ordinary matching Vectors"),
    )
    s.current isa NTuple{3,Vector{T}} &&
        all(a -> axes(a, 1) == Base.OneTo(st.g.n[1]) && all(isfinite, a), s.current) || throw(
        ArgumentError("compatible Yee checkpoint current must be finite ordinary matching Vectors"),
    )
    for charge in (s.charge_n, s.charge_np1)
        charge isa Vector{Rational{BigInt}} && axes(charge, 1) == Base.OneTo(st.g.n[1]) || throw(
            ArgumentError("compatible Yee checkpoint charge must be matching exact-Q Vectors"),
        )
        _pyee_finite(charge)
    end
    incoming = (
        s.species,
        s.fields.E...,
        s.fields.B...,
        s.current...,
        s.charge_n,
        s.charge_np1,
        _implicit_hybrid_particle_arrays(populations)...,
    )
    owned = _pyee_owned_arrays(st)
    for i in eachindex(incoming)
        for j = 1:(i-1)
            _particle_arrays_mightalias(incoming[i], incoming[j]) &&
                throw(ArgumentError("compatible Yee incoming physical roles must be independent"))
        end
        for a in owned
            _particle_arrays_mightalias(incoming[i], a) && throw(
                ArgumentError(
                    "compatible Yee incoming state must be detached from the destination",
                ),
            )
        end
    end
    # Configuration is a fixed scalar/type/policy tuple. Q source/report records
    # have their explicit numeric owner; shared BigInt values are not metadata.
    _validate_checkpoint_configuration(s, st, s.configuration)
    fields = YeeFields{1,T}(s.fields.E, s.fields.B)
    _pyee_field_values(fields, st.g)
    _pyee_charge!(st.work.charge1, populations, st.g, st.shape)
    isequal(st.work.charge1, s.charge_np1) ||
        throw(ArgumentError("compatible Yee incoming endpoint charge differs from its particles"))
    if s.step == 0
        all(a -> all(_particle_stored_zero, a), s.current) && isequal(s.charge_n, s.charge_np1) ||
            throw(
                ArgumentError(
                    "unstepped compatible Yee checkpoint requires zero current/equal sources",
                ),
            )
    end
    scale = _pyee_charge_scale(st, fields, populations, s.charge_np1)
    _pyee_constraints(
        st,
        fields,
        s.current,
        s.charge_np1,
        s.charge_np1,
        zero(T);
        continuity = false,
        scale = scale,
    )
    if s.step > 0
        actual = _pyee_constraints(
            st,
            fields,
            s.current,
            s.charge_n,
            s.charge_np1,
            s.last_dt;
            scale = s.last_report.charge_scale,
        )
        actual == (s.last_report.continuity, s.last_report.gauss) || throw(
            ArgumentError(
                "compatible Yee report differs from saved current/source/field constraints",
            ),
        )
    end
    return populations, fields
end
function _validate_checkpoint_state(s, st::ImplicitYeeEMPICStepper)
    _pyee_checkpoint_data(s, st)
    return nothing
end
function _checkpoint_state(st::ImplicitYeeEMPICStepper, species::AbstractVector)
    _pyee_binding(st, species)
    _pyee_boundary(st)
    _pyee_stage!(st, species, zero(st.c))
    _pyee_ready!(st)
    state = (
        integrator = :implicit_periodic_yee_cloud_midpoint,
        schema = 1,
        configuration = _pyee_configuration(
            st.g,
            st.shape,
            st.c,
            st._species_scalars,
            st.controls,
            st.constraint_tolerances,
        ),
        phase = :synchronous_integer,
        initialized = true,
        species = [_checkpoint_particle_state(ps) for ps in species],
        fields = (E = st.fields.E, B = st.fields.B),
        current = st.J,
        charge_n = st.charge_n,
        charge_np1 = st.charge_np1,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        last_report = st.last_report,
    )
    _validate_checkpoint_configuration(state, st, state.configuration)
    saved = deepcopy(state)
    _validate_checkpoint_state(saved, st)
    return saved
end
function _restore_checkpoint_state!(st::ImplicitYeeEMPICStepper, species::AbstractVector, s)
    _pyee_binding(st, species; values = false, counts = false, sources = false)
    _pyee_checkpoint_data(s, st)
    saved = deepcopy(s)
    incoming, fields = _pyee_checkpoint_data(saved, st)
    prepared, counts = _implicit_em_replacements(incoming)
    bindings = _implicit_em_bindings(incoming)
    _pyee_reserve!(st, incoming)
    return _pyee_commit!(
        st,
        prepared,
        counts,
        bindings,
        fields,
        saved.current,
        saved.charge_n,
        saved.charge_np1,
        saved.time,
        saved.step,
        saved.last_dt,
        saved.last_report,
    )
end

"""
    save_checkpoint(path, st::ImplicitYeeEMPICStepper, species)

Atomically save the distinct compatible periodic Yee family, all physical fields,
interval current, exact source snapshots, original-width/fixed policy and complete
synchronous clock/report. Raw checkpoints have no separate checksum.
"""
function save_checkpoint(path::AbstractString, st::ImplicitYeeEMPICStepper, species::AbstractVector)
    _acquire_pyee!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy = false
    end
end
"Validate original links, detach/revalidate and prepare complete physical/source/phase replacement before publication."
function load_checkpoint!(
    st::ImplicitYeeEMPICStepper,
    species::AbstractVector,
    path::AbstractString,
)
    _acquire_pyee!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy = false
    end
end
function archive_run(
    path::AbstractString,
    st::ImplicitYeeEMPICStepper,
    species::AbstractVector;
    kwargs...,
)
    _acquire_pyee!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
    finally
        st.busy = false
    end
end
"Verify run-container checksum and perform the same prepared restore without projection, priming or solve."
function load_archive!(st::ImplicitYeeEMPICStepper, species::AbstractVector, path::AbstractString)
    _acquire_pyee!(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, species, loaded.state)
    finally
        st.busy = false
    end
end
