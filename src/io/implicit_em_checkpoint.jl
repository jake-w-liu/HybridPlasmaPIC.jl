# Complete synchronous implicit EM state; no explicit leapfrog surrogate.
const _IMPLICIT_EM_CHECKPOINT_FIELDS = (
    :integrator, :schema, :configuration, :phase, :initialized,
    :species, :fields, :time, :step, :last_dt, :last_report,
)
function _implicit_em_checkpoint_keys(s, keys_required)
    _validate_checkpoint_container(s, keys_required)
    length(keys(s)) == length(keys_required) || throw(ArgumentError("implicit EM checkpoint contains unsupported fields"))
    return nothing
end
function _implicit_em_checkpoint_physical!(st::ImplicitEMPICStepper{D,T}, fields, species, step, h, report) where {D,T}
    Q, w = _exact_rational, st.work
    all(a -> all(isfinite, a), _implicit_em_field_arrays(fields)) || throw(ArgumentError("implicit EM checkpoint fields must be finite"))
    if step == 0
        all(a -> all(_particle_stored_zero, a), (fields.J..., fields.Jraw...)) && isequal(fields.rho_n, fields.rho_np1) ||
            throw(ArgumentError("unstepped implicit EM checkpoint requires zero current and equal charges"))
    end
    _implicit_em_set_constraint_scales!(st, fields, species)
    _implicit_em_charge!(w.midpoint.rho_n, st, species)
    mismatch = zero(Rational{BigInt})
    for I in eachindex(fields.rho_np1)
        mismatch = max(mismatch, abs(Q(w.midpoint.rho_n[I]) - Q(fields.rho_np1[I])))
    end
    mismatch / w.charge_scale <= Q(st.constraint_tolerances[1]) + Q(st.constraint_tolerances[2]) ||
        throw(ArgumentError("implicit EM checkpoint endpoint charge differs from its particles"))
    _implicit_em_constraint_values(st, fields, fields.rho_n, zero(T); continuity = false)
    if step > 0
        # Last current/rho are interval data. Use the saved actual reference
        # scales, not a reconstructed time or the destination's old values.
        w.charge_scale, w.divergence_scale = report.constraint_scales
        _implicit_em_constraint_values(st, fields, fields.rho_n, h)
    end
    return nothing
end
function _implicit_em_checkpoint_data(s, st::ImplicitEMPICStepper{D,T}) where {D,T}
    _implicit_em_binding(st, st.owners; values = false, counts = false, sources = false)
    _implicit_em_checkpoint_keys(s, _IMPLICIT_EM_CHECKPOINT_FIELDS)
    s.integrator === :implicit_empic_midpoint && s.schema === 1 &&
        s.phase === :synchronous_integer && s.initialized === true ||
        throw(ArgumentError("unsupported implicit EM checkpoint family/schema/phase"))
    config = _implicit_em_configuration(st.g, st.c, st._species_scalars, st.controls, st.constraint_tolerances)
    typeof(s.configuration) === typeof(config) && isequal(s.configuration, config) ||
        throw(ArgumentError("implicit EM checkpoint fixed configuration must match the destination"))
    _implicit_em_phase(st, s.initialized, s.time, s.step, s.last_dt, s.last_report)
    s.species isa Vector && length(s.species) == length(st.owners) ||
        throw(ArgumentError("implicit EM checkpoint requires its ordered species count"))
    populations = _ElectrostaticHostParticle{D,T}[]
    for (index, record) in enumerate(s.species)
        _implicit_em_checkpoint_keys(record, _CHECKPOINT_PARTICLE_FIELDS)
        ps = _checkpoint_species(record, Val(D), T)
        ps isa _ElectrostaticHostParticle{D,T} && isequal((ps.q, ps.m), st._species_scalars[index]) ||
            throw(ArgumentError("implicit EM checkpoint requires ordinary Vectors and fixed original factors"))
        push!(populations, ps)
    end
    _implicit_em_particles(st, populations)
    references = _implicit_em_field_state(st.fields)
    _implicit_em_checkpoint_keys(s.fields, keys(references))
    _validate_checkpoint_fields(s.fields, references)
    fields_in = Tuple(a for value in values(s.fields) for a in (value isa Tuple ? value : (value,)))
    all(a -> a isa Array{T,D} && all(isfinite, a), fields_in) ||
        throw(ArgumentError("implicit EM checkpoint requires finite ordinary matching field Arrays"))
    incoming = (s.species, fields_in..., _implicit_hybrid_particle_arrays(populations)...)
    owned = (_implicit_em_owned_arrays(st)..., st.g.midx..., st.g.kfull..., st.g.ik..., st.g.kvec...)
    for i in eachindex(incoming)
        for j = 1:(i-1)
            _particle_arrays_mightalias(incoming[i], incoming[j]) &&
                throw(ArgumentError("implicit EM checkpoint physical roles must be independent"))
        end
        for a in owned
            _particle_arrays_mightalias(incoming[i], a) &&
                throw(ArgumentError("implicit EM checkpoint must be detached from destination/grid storage"))
        end
    end
    # Configuration contains only fixed scalar/type/policy tuples. Known Q
    # report records are finite-validated by their owner, not opaque metadata.
    _validate_checkpoint_configuration(s, st, s.configuration)
    fields = _ImplicitEMFields(st.g)
    _restore_checkpoint_fields!(_implicit_em_field_state(fields), s.fields)
    _implicit_em_checkpoint_physical!(st, fields, populations, s.step, s.last_dt, s.last_report)
    return populations, fields
end
function _validate_checkpoint_state(s, st::ImplicitEMPICStepper)
    _implicit_em_checkpoint_data(s, st)
    return nothing
end
function _checkpoint_state(st::ImplicitEMPICStepper, species::AbstractVector)
    _implicit_em_binding(st, species)
    _implicit_em_boundary(st)
    _implicit_em_stage!(st, species, zero(st.c))
    _implicit_em_ready_constraints!(st)
    config = _implicit_em_configuration(st.g, st.c, st._species_scalars, st.controls, st.constraint_tolerances)
    state = (integrator = :implicit_empic_midpoint, schema = 1, configuration = config,
        phase = :synchronous_integer, initialized = true,
        species = [_checkpoint_particle_state(ps) for ps in species],
        fields = _implicit_em_field_state(st.fields), time = st.time[], step = st.step[],
        last_dt = st.last_dt[], last_report = st.last_report)
    _validate_checkpoint_configuration(state, st, config)
    saved = deepcopy(state)
    _validate_checkpoint_state(saved, st)
    return saved
end
function _restore_checkpoint_state!(st::ImplicitEMPICStepper, species::AbstractVector, s)
    # Stale physical/source-binding/count/qm/clock/report values are replaceable;
    # actual owner/storage/configuration and copy roles remain mandatory.
    _implicit_em_binding(st, species; values = false, counts = false, sources = false)
    _implicit_em_checkpoint_data(s, st)
    saved = deepcopy(s)
    incoming, fields = _implicit_em_checkpoint_data(saved, st)
    prepared, counts = _implicit_em_replacements(incoming)
    bindings = _implicit_em_bindings(incoming)
    _implicit_em_reserve!(st, incoming)
    return _implicit_em_commit!(st, prepared, counts, bindings, fields,
        saved.time, saved.step, saved.last_dt, saved.last_report)
end

"""
    save_checkpoint(path, st::ImplicitEMPICStepper, species)

Atomically save distinct schema1 synchronous implicit EM: full fields/current/
charge snapshots, ordered particles, fixed policy, actual clock and detached
complete solve/work report. Nonlinear history is omitted. Initialized completed
boundaries only. The raw checkpoint has no independent checksum.
"""
function save_checkpoint(path::AbstractString, st::ImplicitEMPICStepper, species::AbstractVector)
    _acquire_implicit_em!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy = false
    end
end
"""
    load_checkpoint!(st::ImplicitEMPICStepper, species, path)

Validate original links, detach/revalidate and prepare all physical replacements,
source bindings and solver capacity before one publication. Saved fields and
synchronous phase are authoritative; no public redeposition, projection, field
solve, priming or phase conversion is performed.
"""
function load_checkpoint!(st::ImplicitEMPICStepper, species::AbstractVector, path::AbstractString)
    _acquire_implicit_em!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy = false
    end
end
function archive_run(path::AbstractString, st::ImplicitEMPICStepper, species::AbstractVector; kwargs...)
    _acquire_implicit_em!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); kwargs...)
    finally
        st.busy = false
    end
end
"Verify the existing archive container/checksum, then perform the same prepared implicit EM restore."
function load_archive!(st::ImplicitEMPICStepper, species::AbstractVector, path::AbstractString)
    _acquire_implicit_em!(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, species, loaded.state)
    finally
        st.busy = false
    end
end
