# Synchronized field-only state; coupled closed-Yee schemas remain distinct.
const _MAXWELL_VACUUM_CHECKPOINT_FIELDS = (
    :integrator,
    :schema_version,
    :integrator_policy,
    :field_phase,
    :D,
    :T,
    :ncell,
    :L,
    :dx,
    :periodic,
    :c,
    :time,
    :step,
    :fields,
)

function _checkpoint_state(model::MaxwellVacuum{D,T}) where {D,T}
    _require_finite_positive_real("c", model.c, T)
    _yee_validate_storage(model)
    _yee_validate_fields(model.fields, model.g)
    _validate_checkpoint_clock_state((time = model.time[], step = model.step[]), T)
    return (
        integrator = :yee_maxwell_vacuum,
        schema_version = 1,
        integrator_policy = :synchronized_yee_kdk_native_normwise_v1,
        field_phase = :synchronized,
        D = D,
        T = T,
        ncell = model.g.n,
        L = model.g.L,
        dx = model.g.dx,
        periodic = model.g.periodic,
        c = model.c,
        time = model.time[],
        step = model.step[],
        fields = (E = model.fields.E, B = model.fields.B),
    )
end

function _validate_checkpoint_state(s, model::MaxwellVacuum{D,T}) where {D,T}
    _validate_checkpoint_container(s, _MAXWELL_VACUUM_CHECKPOINT_FIELDS)
    length(keys(s)) == length(_MAXWELL_VACUUM_CHECKPOINT_FIELDS) ||
        throw(ArgumentError("Maxwell vacuum checkpoint contains unsupported fields"))
    s.integrator === :yee_maxwell_vacuum && s.schema_version === 1 &&
        s.integrator_policy === :synchronized_yee_kdk_native_normwise_v1 &&
        s.field_phase === :synchronized ||
        throw(ArgumentError("unsupported Maxwell vacuum checkpoint schema or field policy"))
    s.D === D && s.T === T ||
        throw(ArgumentError("Maxwell vacuum checkpoint dimension and precision must match"))
    for (name, expected) in (
        (:ncell, model.g.n),
        (:L, model.g.L),
        (:dx, model.g.dx),
        (:periodic, model.g.periodic),
        (:c, model.c),
    )
        value = getproperty(s, name)
        typeof(value) === typeof(expected) && isequal(value, expected) ||
            throw(ArgumentError("Maxwell vacuum checkpoint $name does not match the destination"))
    end
    _require_finite_positive_real("c", s.c, T)
    _validate_checkpoint_clock_state(s, T)
    _validate_checkpoint_container(s.fields, (:E, :B))
    length(keys(s.fields)) == 2 ||
        throw(ArgumentError("Maxwell vacuum checkpoint contains unsupported field components"))
    s.fields.E isa NTuple{3,Array{T,D}} && s.fields.B isa NTuple{3,Array{T,D}} ||
        throw(ArgumentError("Maxwell vacuum checkpoint fields must use ordinary CPU arrays"))
    fields = YeeFields{D,T}(s.fields.E, s.fields.B)
    _yee_validate_fields(fields, model.g)
    _yee_check_independent_arrays((fields.E..., fields.B...))
    _yee_validate_pec(fields, model.g)
    return nothing
end

function _restore_checkpoint_state!(model::MaxwellVacuum{D,T}, s) where {D,T}
    _require_finite_positive_real("c", model.c, T)
    # Old field/clock values are replaced. Validate storage without requiring
    # those old fields to satisfy PEC or to be finite.
    _yee_validate_storage(model; check_pec = false)
    _validate_checkpoint_state(s, model)
    # In-memory payloads can permute live arrays. Detach every source together
    # before any copy can overwrite a later source component.
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, model)
    for c = 1:3
        copyto!(model.fields.E[c], saved.fields.E[c])
        copyto!(model.fields.B[c], saved.fields.B[c])
    end
    model.time[] = saved.time
    model.step[] = saved.step
    return model
end

"""
    save_checkpoint(path, model::MaxwellVacuum)

Atomically save synchronized Yee vacuum fields, grid/boundary/speed configuration
and clock in a distinct schema-1 payload. No particles, staging, RNG or prescribed
current sequence is stored. Existing periodic/PEC component shapes and finite
physical fields validate before serialization. Concurrent model mutation is
unsupported. Runtime/restart acceptance is deferred.
"""
function save_checkpoint(path::AbstractString, model::MaxwellVacuum)
    return _save_checkpoint(path, model, _checkpoint_state(model))
end

"""
    load_checkpoint!(model::MaxwellVacuum, path)

Restore the matching field-only grid, precision, boundaries and speed. Incoming
policy, clock, finite staggered fields, independence and PEC data validate and
detach before six ordinary copies and clock stores. Destination arrays/Refs retain
identity; old field/clock values may be stale. No field projection or scratch
restore occurs. Trusted Julia Serialization and the existing ordinary-copy/resource
and exclusive-access premises apply; runtime/restart acceptance remains deferred.
"""
function load_checkpoint!(model::MaxwellVacuum, path::AbstractString)
    return _restore_checkpoint_state!(model, deserialize(path))
end

"""
    archive_run(path, model::MaxwellVacuum; rng_seed, kwargs...)

Archive the same field-only checkpoint payload through the existing checksummed
run container. The seed is provenance only. Default normalization describes the
vacuum field units, and the boundary description follows the actual grid axes;
both may be supplied explicitly. No source waveform or random state is invented.
"""
function archive_run(
    path::AbstractString,
    model::MaxwellVacuum;
    normalization::AbstractString = "epsilon0=1;mu0=1/c^2",
    boundary_desc::AbstractString = join(map(p -> p ? "periodic" : "pec", model.g.periodic), ","),
    kwargs...,
)
    return _archive_run(
        path, model, _checkpoint_state(model);
        normalization, boundary_desc, kwargs...,
    )
end
