# Distinct complete synchronous mapped phase; no Cartesian/half-time migration.
const _SH_CHECKPOINT_KEYS = (:integrator, :schema, :configuration, :phase, :initialized,
    :species, :beta, :fields, :time, :step, :last_dt, :last_report)
const _SH_PARTICLE_KEYS = (:x, :v, :weight, :id, :tag, :q, :m)

function _sh_configuration(st)
    return (dimension = 3, precision = eltype(st.beta[1]), n = st.g.n, L = st.g.L,
        dx = st.g.dx, amplitude = st.mesh.amplitude, tau = st.mesh.tau,
        metric = copy(st.mesh._s), map = :periodic_sinusoidal_shear_J1,
        metric_policy = :bounded_fraction_sincospi_original_factors,
        coordinate_policy = :logical_positions_cartesian_velocities,
        component_policy = :primary_contravariant_flux_physical_cartesian_fields,
        normalization = :Omega_ci_e1_mu0_1, quadrature = :product_stored_logical_dx,
        closure = st.closure, shape = :tsc_complete_original_axis_coefficients,
        moment_policy = :complete_collection_charge_cartesian_current,
        rk_policy = :coupled_four_stages_endpoint_refresh,
        differentiated_source_policy = :reject_ohm_source_risk,
        count = Tuple(st._counts), charges_masses = Tuple(st._qm), fixed_weights = true,
        div_rtol = st.div_rtol, div_atol = st.div_atol, refresh_rtol = st.refresh_rtol,
        divergence_norm = :normalized_modal_component_max,
        refresh_norm = :relative_block_max_no_absolute_tolerance)
end
function _sh_checkpoint_keys(s, expected)
    s isa NamedTuple && keys(s) == expected || throw(ArgumentError("sheared checkpoint keys/schema do not match"))
    return nothing
end
function _sh_checkpoint_configuration(s, st)
    physical = (species = s.species, beta = s.beta, fields = s.fields, last_report = s.last_report)
    return _validate_checkpoint_configuration(physical, st, s.configuration)
end
function _sh_checkpoint_population(s, st::ShearedPeriodicHybridPIC{T}, k) where {T}
    _sh_checkpoint_keys(s, _SH_PARTICLE_KEYS)
    s.x isa NTuple{3,Vector{T}} && s.v isa NTuple{3,Vector{T}} &&
        s.weight isa Vector{T} && s.id isa Vector{UInt64} && s.tag isa Vector{UInt32} &&
        s.q isa T && s.m isa T || throw(ArgumentError("sheared checkpoint requires matching ordinary population leaves"))
    p = ParticleSet{3,T}(s.x, s.v, s.weight, s.id, s.tag, s.q, s.m)
    _sh_particle(p)
    nparticles(p) == st._counts[k] && isequal((p.q, p.m), st._qm[k]) &&
        isequal(p.weight, st._weights[k]) || throw(ArgumentError("sheared checkpoint fixed population inventory differs"))
    for d in 1:3, x in p.x[d]
        0 <= x < st.g.L[d] || throw(ArgumentError("saved shear logical coordinates must be half-open periodic representatives"))
    end
    return p
end
function _sh_relative_block(saved, prepared, rtol)
    Q = _exact_rational
    left = saved isa Tuple ? saved : (saved,)
    right = prepared isa Tuple ? prepared : (prepared,)
    defect = scale = zero(Rational{BigInt})
    for (a, b) in zip(left, right), I in eachindex(a)
        defect = max(defect, abs(Q(a[I]) - Q(b[I])))
        scale = max(scale, abs(Q(a[I])), abs(Q(b[I])))
    end
    defect <= Q(rtol) * scale || throw(ArgumentError("saved sheared derived field failed its relative refresh gate"))
    return nothing
end
function _validate_checkpoint_state(s, st::ShearedPeriodicHybridPIC{T}) where {T}
    _sh_bind(st, st.owners; values = false)
    _sh_checkpoint_keys(s, _SH_CHECKPOINT_KEYS)
    s.integrator === :sheared_periodic_hybrid && s.schema === 1 &&
        s.phase === :synchronous_integer && s.initialized === true ||
        throw(ArgumentError("unsupported sheared checkpoint family/schema/phase"))
    configuration = _sh_configuration(st)
    typeof(s.configuration) === typeof(configuration) && isequal(s.configuration, configuration) ||
        throw(ArgumentError("sheared checkpoint map/closure/fixed configuration differs"))
    s.time isa T && s.step isa Int && s.last_dt isa T || throw(ArgumentError("sheared checkpoint clock types differ"))
    _sh_phase(s.time, s.step, s.last_dt, s.last_report)
    s.species isa Vector && length(s.species) == length(st.owners) ||
        throw(ArgumentError("sheared checkpoint ordered population count differs"))
    populations = [_sh_checkpoint_population(p, st, k) for (k, p) in enumerate(s.species)]
    s.beta isa NTuple{3,Array{T,3}} && typeof(s.fields) === typeof(st.fields) ||
        throw(ArgumentError("sheared checkpoint field precision/storage differs"))
    _sh_checkpoint_keys(s.fields, keys(st.fields))
    field_arrays = (s.beta..., _sh_field_arrays(s.fields)...)
    all(a -> size(a) == st.g.n && all(isfinite, a), field_arrays) &&
        all(>(0), s.fields.n) && all(>(0), s.fields.ninv) ||
        throw(ArgumentError("sheared saved fields must be finite matching arrays with positive density/reciprocal"))
    arrays = Any[s.species, s.configuration.metric, field_arrays...]
    for p in populations
        append!(arrays, _sh_particle_arrays(p))
    end
    live = _sh_storage_arrays(st)
    _visit_checkpoint_objects((st.g.plan, st.g.iplan)) do value
        value isa AbstractArray && push!(live, value)
    end
    for i in eachindex(arrays)
        for j in 1:(i - 1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(ArgumentError("incoming shear physical/configuration storage roles overlap"))
        end
        any(a -> _particle_arrays_mightalias(arrays[i], a), live) &&
            throw(ArgumentError("incoming shear state overlaps live public/private/readonly/plan storage"))
    end
    _sh_checkpoint_configuration(s, st)
    # Original graphs precede detach and scratch mutation. Refresh calculations
    # validate the incoming phase; they never replace its stored physical fields.
    _sh_copy!(st.work.endpoint, populations, s.beta)
    _sh_prepare!(st, st.work.endpoint)
    isequal(s.fields.n, st.work.fields.n) && isequal(s.fields.Ji, st.work.fields.Ji) ||
        throw(ArgumentError("saved shear density/current differ from their complete populations"))
    for name in keys(st.fields)
        name in (:n, :Ji) && continue
        _sh_relative_block(getproperty(s.fields, name), getproperty(st.work.fields, name), st.refresh_rtol)
    end
    return nothing
end
function _checkpoint_state(st::ShearedPeriodicHybridPIC, species::AbstractVector)
    _sh_bind(st, species)
    _sh_boundary(st)
    fields = NamedTuple{keys(st.fields)}(map(values(st.fields)) do value
        value isa Tuple ? map(copy, value) : copy(value)
    end)
    saved = (integrator = :sheared_periodic_hybrid, schema = 1,
        configuration = _sh_configuration(st), phase = :synchronous_integer, initialized = true,
        species = [(x = map(copy, p.x), v = map(copy, p.v), weight = copy(p.weight),
            id = copy(p.id), tag = copy(p.tag), q = p.q, m = p.m) for p in species],
        beta = map(copy, st.beta), fields = fields, time = st.time, step = st.step,
        last_dt = st.last_dt, last_report = st.last_report)
    _sh_checkpoint_configuration(saved, st)
    return saved
end
function _restore_checkpoint_state!(st::ShearedPeriodicHybridPIC, species::AbstractVector, incoming)
    _sh_bind(st, species; values = false)
    _validate_checkpoint_state(incoming, st)
    saved = deepcopy(incoming)
    _validate_checkpoint_state(saved, st)
    populations = [_sh_checkpoint_population(p, st, k) for (k, p) in enumerate(saved.species)]
    _sh_copy!(st.work.endpoint, populations, saved.beta)
    return _sh_publish!(st, st.work.endpoint, saved.fields, saved.time, saved.step,
        saved.last_dt, saved.last_report, populations)
end

"Atomically save a complete synchronous mapped phase through the canonical writer."
function save_checkpoint(path::AbstractString, st::ShearedPeriodicHybridPIC, species::AbstractVector = st.owners)
    _sh_acquire(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy = false
    end
end
function load_checkpoint!(st::ShearedPeriodicHybridPIC, species::AbstractVector, path::AbstractString)
    _sh_acquire(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy = false
    end
end
load_checkpoint!(st::ShearedPeriodicHybridPIC, path::AbstractString) = load_checkpoint!(st, st.owners, path)
function archive_run(path::AbstractString, st::ShearedPeriodicHybridPIC, species::AbstractVector = st.owners;
    normalization::AbstractString = "Omega_ci; periodic shear logical x, Cartesian v; e=mu0=1", kwargs...)
    _sh_acquire(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); normalization, kwargs...)
    finally
        st.busy = false
    end
end
function load_archive!(st::ShearedPeriodicHybridPIC, species::AbstractVector, path::AbstractString)
    _sh_acquire(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, species, loaded.state)
    finally
        st.busy = false
    end
end
load_archive!(st::ShearedPeriodicHybridPIC, path::AbstractString) = load_archive!(st, st.owners, path)
