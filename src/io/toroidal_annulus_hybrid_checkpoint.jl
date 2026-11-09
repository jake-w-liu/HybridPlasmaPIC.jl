# Closed synchronous annular phase; raw and archive transports share this owner.
const _TAH_CHECKPOINT_KEYS = (
    :integrator,
    :schema,
    :configuration,
    :phase,
    :initialized,
    :species,
    :beta,
    :fields,
    :time,
    :step,
    :last_dt,
    :initial_energy,
    :outward_work,
    :controller_work,
    :stage_defect,
    :rk_defect,
    :last_report,
    :last_endpoint,
)
const _TAH_PARTICLE_KEYS = (:x, :v, :weight, :id, :tag, :q, :m)

function _tah_checkpoint_keys(saved, expected)
    saved isa NamedTuple && keys(saved) == expected ||
        throw(ArgumentError("annular checkpoint keys/schema do not match"))
    return nothing
end

function _tah_checkpoint_population(saved, st::ToroidalAnnulusHybridPIC{T}, index) where {T}
    _tah_checkpoint_keys(saved, _TAH_PARTICLE_KEYS)
    saved.x isa NTuple{3,Vector{T}} && saved.v isa NTuple{3,Vector{T}} &&
        saved.weight isa Vector{T} && saved.id isa Vector{UInt64} &&
        saved.tag isa Vector{UInt32} && saved.q isa T && saved.m isa T ||
        throw(ArgumentError("annular checkpoint population precision/ordinary storage differs"))
    population = ParticleSet{3,T}(saved.x, saved.v, saved.weight,
        saved.id, saved.tag, saved.q, saved.m)
    _tah_particle_storage(population, T)
    _tah_particle_values(population, st.mesh)
    nparticles(population) == st._counts[index] &&
        isequal((population.q, population.m), st._qm[index]) &&
        isequal(population.weight, st._weights[index]) ||
        throw(ArgumentError("annular checkpoint fixed population factors/count/weights differ"))
    all(d -> all(x -> 0 <= x < 1, population.x[d]), 2:3) ||
        throw(ArgumentError("annular saved turns must be half-open representatives"))
    return population
end

function _tah_checkpoint_arrays(saved, populations)
    arrays = Any[saved.species, saved.beta..., _tah_field_arrays(saved.fields)...]
    for population in populations
        append!(arrays, _tah_particle_arrays(population))
    end
    return arrays
end

function _tah_checkpoint_live_graph(st)
    objects = IdDict{Any,Nothing}()
    arrays = AbstractArray[]
    function remember(value)
        (ismutable(value) || value isa AbstractArray) || return nothing
        if !haskey(objects, value)
            objects[value] = nothing
            value isa AbstractArray && push!(arrays, value)
        end
        return nothing
    end
    # Include nonlogical backing roles and every mutable scalar/metadata role,
    # including committed reports and retained plan children omitted from payload.
    record(value, logical, edges) = remember(value)
    _visit_checkpoint_objects(remember, st; backing = remember, record)
    foreach(remember, _tah_owned_arrays(st))
    return objects, arrays
end

function _tah_checkpoint_graph(saved, st, arrays)
    _tah_independent(arrays)
    live_objects, live_arrays = _tah_checkpoint_live_graph(st)
    function require_independent(value)
        (ismutable(value) || value isa AbstractArray) || return nothing
        haskey(live_objects, value) &&
            throw(ArgumentError("incoming annular state retains live mutable storage or metadata"))
        if value isa AbstractArray
            any(a -> _particle_arrays_mightalias(value, a), live_arrays) &&
                throw(ArgumentError("incoming annular state overlaps live public/private/geometry/plan storage"))
        end
        return nothing
    end
    # Historical rational receipts may share scalar values with each other;
    # they may never retain mutable scalar objects from the live destination.
    record(value, logical, edges) = require_independent(value)
    _visit_checkpoint_objects(require_independent, saved; backing = require_independent, record)
    return nothing
end

function _tah_checkpoint_layout(saved, st::ToroidalAnnulusHybridPIC{T}) where {T}
    _tah_storage(st, st.owners; values = false)
    _tah_checkpoint_keys(saved, _TAH_CHECKPOINT_KEYS)
    saved.integrator === :toroidal_annulus_hybrid && saved.schema === 1 &&
        saved.phase === :synchronous_integer && saved.initialized === true ||
        throw(ArgumentError("unsupported annular checkpoint family/schema/phase"))
    typeof(saved.configuration) === typeof(st._configuration) &&
        isequal(saved.configuration, st._configuration) ||
        throw(ArgumentError("annular checkpoint fixed geometry/source/closure/boundary policy differs"))
    saved.time isa T && saved.step isa Int && saved.last_dt isa T &&
        saved.last_report isa Union{Nothing,ToroidalAnnulusStepReport{T}} &&
        saved.last_endpoint isa ToroidalAnnulusStageReport{T} &&
        all(_sce_q_finite, (saved.initial_energy, saved.outward_work,
            saved.controller_work, saved.stage_defect, saved.rk_defect)) ||
        throw(ArgumentError("annular checkpoint phase/receipt/exact-ledger types differ"))
    saved.species isa Vector && length(saved.species) == length(st.owners) &&
        all(index -> isassigned(saved.species, index), eachindex(saved.species)) ||
        throw(ArgumentError("annular checkpoint ordered population layout differs"))
    populations = [_tah_checkpoint_population(p, st, index)
        for (index, p) in enumerate(saved.species)]
    saved.beta isa NTuple{3,Array{T,3}} && typeof(saved.fields) === typeof(st.fields) ||
        throw(ArgumentError("annular checkpoint primary/derived field precision/storage differs"))
    _tah_checkpoint_keys(saved.fields, keys(st.fields))
    all(a -> size(a) == st.mesh.n && all(isfinite, a),
        (saved.beta..., _tah_field_arrays(saved.fields)...)) &&
        all(>(zero(T)), saved.fields.n) && all(>(zero(T)), saved.fields.ninv) &&
        all(>=(zero(T)), saved.fields.pe) ||
        throw(ArgumentError("annular saved fields require finite matching arrays, positive density/reciprocal and nonnegative pressure"))
    _tah_checkpoint_graph(saved, st, _tah_checkpoint_arrays(saved, populations))
    return populations
end

function _tah_checkpoint_stage_graph(stage, saved, populations, st)
    candidates = _tah_stage_arrays(stage)
    _tah_independent(candidates)
    sources = _tah_checkpoint_arrays(saved, populations)
    _, live = _tah_checkpoint_live_graph(st)
    for candidate in candidates
        any(a -> _particle_arrays_mightalias(candidate, a), sources) &&
            throw(ArgumentError("annular prepared committed shadow overlaps saved source roles"))
        any(a -> _particle_arrays_mightalias(candidate, a), live) &&
            throw(ArgumentError("annular prepared committed shadow overlaps live storage"))
    end
    return nothing
end

function _tah_validate_checkpoint_locked(saved, st::ToroidalAnnulusHybridPIC{T}) where {T}
    populations = _tah_checkpoint_layout(saved, st)
    _tah_normal_flux(st, saved.beta)
    _tah_contact_fields(st, saved.fields)
    reservoirs = _tah_reservoirs(st, populations, saved.beta, saved.fields.pe)
    _tah_phase(st, reservoirs, saved.initialized, saved.time, saved.step,
        saved.last_dt, saved.initial_energy, saved.outward_work, saved.controller_work,
        saved.stage_defect, saved.rk_defect, saved.last_report, saved.last_endpoint)
    stage = _tah_make_stage(populations, saved.beta, saved.fields.pe, T)
    _tah_checkpoint_stage_graph(stage, saved, populations, st)
    refreshed = _tah_rhs!(st.work.endpoint_rate, st, stage, saved.last_dt)
    _tah_refresh(saved.fields, st.work.fields, st.refresh_rtol)
    _tah_endpoint_receipt(saved.last_endpoint, refreshed, saved.fields, st)
    return stage
end

# Canonical raw/archive writers invoke this under their caller's operation guard.
# Acquiring here would recursively reject save/archive after successful capture.
function _validate_checkpoint_state(saved, st::ToroidalAnnulusHybridPIC)
    _tah_validate_checkpoint_locked(saved, st)
    return nothing
end

function _tah_checkpoint_state_locked(st::ToroidalAnnulusHybridPIC, species)
    _tah_ready(st, species)
    # One complete detach also separates the BigInt children of exact reports
    # and ledgers. No captured mutable leaf remains a live simulation reference.
    saved = deepcopy((
        integrator = :toroidal_annulus_hybrid,
        schema = 1,
        configuration = st._configuration,
        phase = :synchronous_integer,
        initialized = st.initialized,
        species = [(x = p.x, v = p.v, weight = p.weight, id = p.id, tag = p.tag,
            q = p.q, m = p.m) for p in species],
        beta = st.beta,
        fields = st.fields,
        time = st.time,
        step = st.step,
        last_dt = st.last_dt,
        initial_energy = st.initial_energy,
        outward_work = st.outward_work,
        controller_work = st.controller_work,
        stage_defect = st.stage_defect,
        rk_defect = st.rk_defect,
        last_report = st.last_report,
        last_endpoint = st.last_endpoint,
    ))
    _tah_validate_checkpoint_locked(saved, st)
    return saved
end

function _tah_restore_checkpoint_locked!(st::ToroidalAnnulusHybridPIC, species, incoming)
    _tah_storage(st, species; values = false)
    # The original graph must reject before deepcopy can conceal its aliases.
    _tah_checkpoint_layout(incoming, st)
    saved = deepcopy(incoming)
    stage = _tah_validate_checkpoint_locked(saved, st)
    return _tah_publish!(st, stage, saved.fields, saved.time, saved.step, saved.last_dt,
        saved.initial_energy, saved.outward_work, saved.controller_work,
        saved.stage_defect, saved.rk_defect, saved.last_report, saved.last_endpoint;
        provenance = true)
end

function _checkpoint_state(st::ToroidalAnnulusHybridPIC, species::AbstractVector = st.owners)
    _tah_acquire(st)
    try
        return _tah_checkpoint_state_locked(st, species)
    finally
        _tah_release(st)
    end
end

function _restore_checkpoint_state!(st::ToroidalAnnulusHybridPIC,
    species::AbstractVector, incoming)
    _tah_acquire(st)
    try
        return _tah_restore_checkpoint_locked!(st, species, incoming)
    finally
        _tah_release(st)
    end
end
_restore_checkpoint_state!(st::ToroidalAnnulusHybridPIC, incoming) =
    _restore_checkpoint_state!(st, st.owners, incoming)

"Atomically save a complete synchronous annular phase through the canonical writer."
function save_checkpoint(path::AbstractString, st::ToroidalAnnulusHybridPIC,
    species::AbstractVector = st.owners)
    _tah_acquire(st)
    try
        return _save_checkpoint(path, st, _tah_checkpoint_state_locked(st, species))
    finally
        _tah_release(st)
    end
end

function load_checkpoint!(st::ToroidalAnnulusHybridPIC,
    species::AbstractVector, path::AbstractString)
    _tah_acquire(st)
    try
        return _tah_restore_checkpoint_locked!(st, species, deserialize(path))
    finally
        _tah_release(st)
    end
end
load_checkpoint!(st::ToroidalAnnulusHybridPIC, path::AbstractString) =
    load_checkpoint!(st, st.owners, path)

function archive_run(path::AbstractString, st::ToroidalAnnulusHybridPIC,
    species::AbstractVector = st.owners;
    rng_seed::Integer = 0,
    normalization::AbstractString = "Omega_ci; toroidal annulus logical x, Cartesian v; e=mu0=1",
    boundary_desc::AbstractString = "equipotential contacts; incoming electron-pressure reservoirs; stage-domain ions",
    kwargs...)
    _tah_acquire(st)
    try
        return _archive_run(path, st, _tah_checkpoint_state_locked(st, species);
            rng_seed, normalization, boundary_desc, kwargs...)
    finally
        _tah_release(st)
    end
end

function load_archive!(st::ToroidalAnnulusHybridPIC,
    species::AbstractVector, path::AbstractString)
    _tah_acquire(st)
    try
        loaded = load_run(path)
        return _tah_restore_checkpoint_locked!(st, species, loaded.state)
    finally
        _tah_release(st)
    end
end
load_archive!(st::ToroidalAnnulusHybridPIC, path::AbstractString) =
    load_archive!(st, st.owners, path)
