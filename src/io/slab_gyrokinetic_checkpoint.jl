# Distinct physical-phase persistence; no ParticleSet/half-time surrogate.
const _SG_CHECKPOINT_KEYS = (:integrator, :schema, :configuration, :phase,
    :initialized, :populations, :rho_gc, :phi, :time, :step, :last_dt, :last_report)
const _SG_POPULATION_KEYS = (:X, :vpar, :mu, :weight, :id, :tag, :q, :m)

function _sg_configuration(st)
    return (dimension = 3, precision = eltype(st.phi), n = st.g.n, L = st.g.L,
        dx = st.g.dx, B0 = st.B0, Te = st.Te, e = st.e,
        chi = st.chi, alpha = st.alpha, ne0_display = st.ne0_display,
        counts = Tuple(st._counts), charges_masses = Tuple(st._qm),
        nring = st.nring, phases = copy(st._phases), shape = :tsc,
        background_policy = :exact_inventory,
        coefficient_policy = :rounded_equilibrium_susceptibilities,
        electron_policy = :linear_modified_adiabatic_yz_average,
        polarization_policy = :equilibrium_long_wavelength_full_k_perpendicular,
        normalization = :coherent_si_form_Te_energy,
        integrator_policy = :collective_rk4_four_stages_endpoint_refresh,
        field_rtol = st.field_rtol, field_atol = st.field_atol,
        field_norm = :normalized_modal_component_max,
        fixed_invariants = :particle_counts_charge_mass_mu_fullf_weight)
end
function _sg_checkpoint_keys(s, expected)
    s isa NamedTuple && keys(s) == expected ||
        throw(ArgumentError("gyrokinetic checkpoint fields/schema do not match"))
    return nothing
end
function _sg_checkpoint_population(s, st::SlabGyrokineticPIC{T}, k) where {T}
    _sg_checkpoint_keys(s, _SG_POPULATION_KEYS)
    s.X isa NTuple{3,Vector{T}} && s.vpar isa Vector{T} && s.mu isa Vector{T} &&
        s.weight isa Vector{T} && s.id isa Vector{UInt64} && s.tag isa Vector{UInt32} &&
        s.q isa T && s.m isa T || throw(ArgumentError("gyrokinetic checkpoint population storage/precision differs"))
    p = GyrocentreSet{T}(s.X, s.vpar, s.mu, s.weight, s.id, s.tag, s.q, s.m)
    _sg_population(p)
    nparticles(p) == st._counts[k] && isequal((p.q, p.m), st._qm[k]) &&
        isequal(p.mu, st._mu[k]) && isequal(p.weight, st._weight[k]) ||
        throw(ArgumentError("gyrokinetic checkpoint fixed population inventory/mu differs"))
    for d in 1:3, x in p.X[d]
        0 <= x < st.g.L[d] || throw(ArgumentError("gyrokinetic saved coordinates must be published periodic representatives"))
    end
    return p
end
function _sg_checkpoint_report(report::SlabGyrokineticFieldReport{T}) where {T}
    values = (report.residual, report.tolerance, report.source_norm, report.gauge,
        report.gauge_tolerance, report.ion_dc_mean, report.ion_dc_imaginary,
        report.ion_dc_inventory_defect, report.charge_defect, report.electron_response_amplitude)
    any(isnan, values) && throw(ArgumentError("gyrokinetic checkpoint receipt contains NaN"))
    all(>=(zero(T)), (report.residual, report.tolerance, report.source_norm,
        report.gauge, report.gauge_tolerance, report.electron_response_amplitude)) ||
        throw(ArgumentError("gyrokinetic checkpoint receipt norm is negative"))
    report.range_resolved === all(isfinite, values) ||
        throw(ArgumentError("gyrokinetic checkpoint receipt range flag disagrees with its displayed scalars"))
    return nothing
end
# The generic graph validator registers payload arrays as owned. Configuration
# itself must remain the independent graph being checked, so register only the
# physical payload (plus the complete destination/scratch/FFT graph it visits).
function _sg_checkpoint_configuration(s, st)
    physical = (populations = s.populations, rho_gc = s.rho_gc, phi = s.phi,
        last_report = s.last_report)
    return _validate_checkpoint_configuration(physical, st, s.configuration)
end
function _validate_checkpoint_state(s, st::SlabGyrokineticPIC{T}) where {T}
    _sg_bind(st, st.owners; values = false)
    _sg_checkpoint_keys(s, _SG_CHECKPOINT_KEYS)
    s.integrator === :slab_gyrokinetic_fullf && s.schema === 1 &&
        s.phase === :synchronous_integer && s.initialized === true ||
        throw(ArgumentError("unsupported gyrokinetic checkpoint family/schema/phase"))
    configuration = _sg_configuration(st)
    typeof(s.configuration) === typeof(configuration) && isequal(s.configuration, configuration) ||
        throw(ArgumentError("gyrokinetic checkpoint fixed physics/quadrature/field policy differs"))
    s.time isa T && s.step isa Int && s.last_dt isa T &&
        s.last_report isa SlabGyrokineticFieldReport{T} ||
        throw(ArgumentError("gyrokinetic checkpoint clock/receipt types differ"))
    _sg_phase(true, s.time, s.step, s.last_dt, s.last_report)
    _sg_checkpoint_report(s.last_report)
    s.populations isa Vector && length(s.populations) == length(st.owners) ||
        throw(ArgumentError("gyrokinetic checkpoint ordered population count differs"))
    populations = [_sg_checkpoint_population(p, st, k) for (k,p) in enumerate(s.populations)]
    s.rho_gc isa Array{T,3} && s.phi isa Array{T,3} &&
        size(s.rho_gc) == st.g.n && size(s.phi) == st.g.n &&
        all(isfinite, s.rho_gc) && all(isfinite, s.phi) ||
        throw(ArgumentError("gyrokinetic checkpoint fields must be finite matching ordinary Arrays"))
    arrays = Any[s.populations, s.rho_gc, s.phi, s.configuration.phases]
    for p in populations
        append!(arrays, _sg_arrays(p))
    end
    live = _sg_storage_arrays(st)
    push!(live, st.owners, st._components)
    _visit_checkpoint_objects((st.g.plan, st.g.iplan)) do value
        value isa AbstractArray && push!(live, value)
    end
    for i in eachindex(arrays)
        for j in 1:(i-1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) &&
                throw(ArgumentError("gyrokinetic incoming physical/configuration roles must be independent"))
        end
        for a in live
            _particle_arrays_mightalias(arrays[i], a) &&
                throw(ArgumentError("gyrokinetic checkpoint must be detached from live public/private/metadata storage"))
        end
    end
    _sg_checkpoint_configuration(s, st)
    # Original incoming roles are checked before any detach or scratch mutation.
    # These are validation calculations; the stored fields are not replaced by a solve.
    _sg_copy!(st.work.endpoint, populations)
    _sg_deposit!(st, st.work.endpoint)
    isequal(st.work.rho, s.rho_gc) ||
        throw(ArgumentError("gyrokinetic saved density differs from its full-f populations/quadrature"))
    copyto!(st.work.phi, s.phi)
    source_scale = _sg_fft!(st.work.source, s.rho_gc, st.g, st._gain)
    phi_scale = _sg_fft!(st.work.phihat, s.phi, st.g, one(T))
    _sg_field_receipt(st, source_scale, phi_scale)
    return nothing
end
function _checkpoint_state(st::SlabGyrokineticPIC, populations::AbstractVector)
    _sg_bind(st, populations)
    _sg_boundary(st)
    saved = (integrator = :slab_gyrokinetic_fullf, schema = 1,
        configuration = _sg_configuration(st), phase = :synchronous_integer,
        initialized = true,
        populations = [(X = map(copy, p.X), vpar = copy(p.vpar), mu = copy(p.mu),
            weight = copy(p.weight), id = copy(p.id), tag = copy(p.tag), q = p.q, m = p.m)
            for p in populations], rho_gc = copy(st.rho_gc), phi = copy(st.phi),
        time = st.time, step = st.step, last_dt = st.last_dt, last_report = st.last_report)
    _sg_checkpoint_configuration(saved, st)
    return saved
end
function _restore_checkpoint_state!(st::SlabGyrokineticPIC{T}, populations::AbstractVector, s) where {T}
    _sg_bind(st, populations; values = false)
    _validate_checkpoint_state(s, st)
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    incoming = [_sg_checkpoint_population(p, st, k) for (k,p) in enumerate(saved.populations)]
    _sg_copy!(st.work.endpoint, incoming)
    # All storage, values, graphs, scalar conversions and private computations
    # are complete. Fixed-length ordinary stores are the sole remaining action.
    for k in eachindex(st.owners), i in eachindex(st.owners[k].weight)
        @inbounds st.owners[k].mu[i] = incoming[k].mu[i]
        @inbounds st.owners[k].weight[i] = incoming[k].weight[i]
        @inbounds st.owners[k].id[i] = incoming[k].id[i]
        @inbounds st.owners[k].tag[i] = incoming[k].tag[i]
    end
    return _sg_publish!(st, st.work.endpoint, saved.rho_gc, saved.phi,
        saved.time, saved.step, saved.last_dt, saved.last_report)
end

"""Atomically save complete schema-1 full-f gyrocentre fields, invariants and synchronous phase."""
function save_checkpoint(path::AbstractString, st::SlabGyrokineticPIC,
    populations::AbstractVector = st.owners)
    _sg_acquire(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, populations))
    finally
        st.busy = false
    end
end
"""Validate/detach the complete incoming phase, then restore stored physical fields verbatim."""
function load_checkpoint!(st::SlabGyrokineticPIC, populations::AbstractVector, path::AbstractString)
    _sg_acquire(st)
    try
        return _restore_checkpoint_state!(st, populations, deserialize(path))
    finally
        st.busy = false
    end
end
function load_checkpoint!(st::SlabGyrokineticPIC, path::AbstractString)
    return load_checkpoint!(st, st.owners, path)
end
function archive_run(path::AbstractString, st::SlabGyrokineticPIC,
    populations::AbstractVector = st.owners;
    normalization::AbstractString = "coherent SI-form slab gyrokinetic (Te energy)", kwargs...)
    _sg_acquire(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, populations); normalization, kwargs...)
    finally
        st.busy = false
    end
end

"""Verify the archive checksum, then validate/detach and restore its complete physical phase."""
function load_archive!(st::SlabGyrokineticPIC, populations::AbstractVector, path::AbstractString)
    _sg_acquire(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, populations, loaded.state)
    finally
        st.busy = false
    end
end
function load_archive!(st::SlabGyrokineticPIC, path::AbstractString)
    return load_archive!(st, st.owners, path)
end
