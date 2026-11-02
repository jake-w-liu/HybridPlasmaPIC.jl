# Complete synchronous open-Hall state. Plans, projection work and trial history
# remain private; no saved physical field or source is reconstructed on restore.
const _OPEN_HALL_CHECKPOINT_FIELDS = (
    :integrator, :schema, :configuration, :phase, :initialized, :species,
    :fields, :floor_count, :time, :step, :last_dt, :source, :rng,
    :rng_algorithm, :runtime, :last_report,
)

function _open_hall_checkpoint_keys(s, required)
    _validate_checkpoint_container(s, required)
    length(keys(s)) == length(required) ||
        throw(ArgumentError("open Hall checkpoint contains unsupported fields"))
    return nothing
end

_open_hall_checkpoint_source_arrays(source) =
    (source.alpha, source.next_index, source.initial_next_index,
        source.initial_totals, source.births, source.exits, source.last_births, source.last_exits)

function _open_hall_checkpoint_rng_metadata(st, s)
    # Register incoming state and both destination streams before checking only
    # opaque stock metadata. Registering the whole RNG as metadata against its
    # own recurrence buffers would reject every valid MersenneTwister.
    if s.rng === st.rng || s.rng === st.rng_candidate
        throw(ArgumentError("open Hall checkpoint RNG must be detached from destination streams"))
    end
    _validate_checkpoint_configuration_registered(_collisional_rng_metadata(s.rng)) do remember, backing, owned
        _open_hall_register_owned!(remember, backing, owned, st)
        for value in (s.species, s.fields, s.source, s.last_report, _collisional_rng_buffers(s.rng))
            _visit_checkpoint_objects(remember, value; backing)
        end
    end
    return nothing
end

function _open_hall_checkpoint_data(s, st::OpenHallHybrid{D,T}) where {D,T}
    # Existing values/phase/counts/RNG metadata are replaceable; fixed ordinary
    # storage, object identities, policy and stock copy destinations are not.
    _open_hall_binding(st, st.owners; values = false, counts = false)
    _open_hall_checkpoint_keys(s, _OPEN_HALL_CHECKPOINT_FIELDS)
    s.integrator === :open_hall_hybrid && s.schema === 1 &&
        s.phase === :synchronous_integer && s.initialized === true ||
        throw(ArgumentError("unsupported open Hall checkpoint family, schema or phase"))
    configuration = _open_hall_bound_configuration(st)
    typeof(s.configuration) === typeof(configuration) && isequal(s.configuration, configuration) ||
        throw(ArgumentError("open Hall checkpoint fixed configuration must match the destination"))
    runtime = _collisional_runtime()
    typeof(s.runtime) === typeof(runtime) && isequal(s.runtime, runtime) ||
        throw(ArgumentError("open Hall checkpoint requires the same Julia runtime and platform"))
    s.rng_algorithm === typeof(st.rng) && typeof(s.rng) === typeof(st.rng) ||
        throw(ArgumentError("open Hall checkpoint stock RNG algorithm must match the destination"))
    _validate_collisional_rng(s.rng)
    _open_hall_phase_values(st, s.initialized, s.time, s.step, s.last_dt, s.last_report)
    s.floor_count isa Int && s.floor_count == 0 ||
        throw(ArgumentError("open Hall checkpoint requires an inactive density floor"))

    s.species isa Vector && length(s.species) == length(st.owners) ||
        throw(ArgumentError("open Hall checkpoint requires its ordered population count"))
    populations = _ElectrostaticHostParticle{D,T}[]
    for (index, record) in enumerate(s.species)
        _open_hall_checkpoint_keys(record, _CHECKPOINT_PARTICLE_FIELDS)
        ps = _checkpoint_species(record, Val(D), T)
        ps isa _ElectrostaticHostParticle{D,T} &&
            isequal((ps.q, ps.m), st._species_scalars[index]) ||
            throw(ArgumentError("open Hall checkpoint requires ordinary Vectors and fixed original charges/masses"))
        push!(populations, ps)
    end
    _open_hall_particles(st, populations)
    _open_hall_source_values(st, s.source, populations, s.step; report = s.last_report)

    references = _hybrid_field_state(st.fields)
    _open_hall_checkpoint_keys(s.fields, keys(references))
    _validate_checkpoint_fields(s.fields, references)
    field_arrays = Tuple(a for value in values(s.fields)
        for a in (value isa Tuple ? value : (value,)))
    all(a -> a isa Array{T,D} && all(isfinite, a), field_arrays) ||
        throw(ArgumentError("open Hall checkpoint fields must be finite ordinary matching Arrays"))
    incoming = (s.species, field_arrays..., _open_hall_particle_arrays(populations)...,
        _open_hall_checkpoint_source_arrays(s.source)..., _collisional_rng_buffers(s.rng)...)
    destination = _open_hall_owned_arrays(st)
    for i in eachindex(incoming)
        for j = 1:(i-1)
            _particle_arrays_mightalias(incoming[i], incoming[j]) &&
                throw(ArgumentError("open Hall checkpoint physical/source/RNG arrays must be independent"))
        end
        for a in destination
            _particle_arrays_mightalias(incoming[i], a) &&
                throw(ArgumentError("open Hall checkpoint must be detached from destination public/private storage"))
        end
    end
    _open_hall_checkpoint_rng_metadata(st, s)
    # These are concrete scalar configuration/report/source records. Their Q
    # values have canonical finite-scalar validation; legitimate repeated BigInt
    # values are not arbitrary mutable configuration requiring identity recovery.
    fields = HybridFields{D,T}(st.op.n)
    _restore_checkpoint_fields!(_hybrid_field_state(fields), s.fields)
    fields.floor_count[] = s.floor_count
    _open_hall_fields(st, fields)
    _open_hall_direct_state(st, fields.B)
    return populations, fields
end

function _validate_checkpoint_state(s, st::OpenHallHybrid)
    _open_hall_checkpoint_data(s, st)
    return nothing
end

function _checkpoint_state(st::OpenHallHybrid, species::AbstractVector)
    _open_hall_boundary(st, species)
    state = (integrator = :open_hall_hybrid, schema = 1,
        configuration = _open_hall_bound_configuration(st), phase = :synchronous_integer,
        initialized = true,
        species = [_checkpoint_particle_state(ps) for ps in species],
        fields = _hybrid_field_state(st.fields), floor_count = st.fields.floor_count[],
        time = st.time[], step = st.step[], last_dt = st.last_dt[],
        source = _open_hall_source_state(st), rng = st.rng,
        rng_algorithm = typeof(st.rng), runtime = _collisional_runtime(),
        last_report = st.last_report)
    # Live metadata links were checked by the boundary owner before detachment.
    # One whole-graph copy preserves actual physical/source/stream values.
    saved = deepcopy(state)
    _validate_checkpoint_state(saved, st)
    return saved
end

function _restore_checkpoint_state!(st::OpenHallHybrid, species::AbstractVector, s)
    _open_hall_binding(st, species; values = false, counts = false)
    # Validate original links before deepcopy can hide a retained live reference.
    _open_hall_checkpoint_data(s, st)
    saved = deepcopy(s)
    populations, fields = _open_hall_checkpoint_data(saved, st)
    prepared, counts = _open_hall_replacements(st, populations)
    _open_hall_reserve!(st, populations)
    # All validation, direct constraints, allocation and stock-copy admission
    # precede the canonical commit's first public particle/field/source store.
    return _open_hall_commit!(st, prepared, counts, fields, saved.source, saved.rng,
        saved.initialized, saved.time, saved.step, saved.last_dt, saved.last_report)
end

"""
    save_checkpoint(path, st::OpenHallHybrid, species)

Atomically save distinct schema1 `:open_hall_hybrid`: complete physical fields,
ordered particles, synchronous phase, exact source/ID/energy histories, fixed
policies, detached last report and full stock RNG. Initialized completed
boundaries only. The raw checkpoint has no separate checksum; `archive_run`
uses the existing checksummed run container.
"""
function save_checkpoint(path::AbstractString, st::OpenHallHybrid, species::AbstractVector)
    _acquire_open_hall!(st)
    try
        return _save_checkpoint(path, st, _checkpoint_state(st, species))
    finally
        st.busy = false
    end
end

"""
    load_checkpoint!(st::OpenHallHybrid, species, path)

Validate matching policy and incoming physical/source/RNG state, detach and
prepare every replacement, then restore verbatim. Old destination values,
clock/report and particle counts are replaceable. No injection, redeposition,
field solve, projection, recentering or priming reconstructs saved state.
"""
function load_checkpoint!(st::OpenHallHybrid, species::AbstractVector, path::AbstractString)
    _acquire_open_hall!(st)
    try
        return _restore_checkpoint_state!(st, species, deserialize(path))
    finally
        st.busy = false
    end
end

"""
    archive_run(path, st::OpenHallHybrid, species; rng_seed, kwargs...)

Save the same complete state through the checksummed run archive. `rng_seed` is
caller provenance; continuation uses the saved full stream, not a reseeded RNG.
"""
function archive_run(path::AbstractString, st::OpenHallHybrid, species::AbstractVector;
    boundary_desc::AbstractString = "open Hall x reservoirs; periodic transverse", kwargs...)
    _acquire_open_hall!(st)
    try
        return _archive_run(path, st, _checkpoint_state(st, species); boundary_desc, kwargs...)
    finally
        st.busy = false
    end
end

"""
    load_archive!(st::OpenHallHybrid, species, path)

Verify the existing run-container schema and byte checksum, then use the same
strict prepared restore as `load_checkpoint!`. To inspect metadata without
restoring, use `load_archive(path)`.
"""
function load_archive!(st::OpenHallHybrid, species::AbstractVector, path::AbstractString)
    _acquire_open_hall!(st)
    try
        loaded = load_run(path)
        return _restore_checkpoint_state!(st, species, loaded.state)
    finally
        st.busy = false
    end
end
