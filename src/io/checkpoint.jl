# checkpoint.jl — checkpoint/restart of the deterministic simulation state.
#
# Stepping (step!) uses no RNG, so the state needed for a bitwise-identical
# restart is: model and shape, particle phase space + provenance, canonical B,
# carried E, time/step counters, and the previous particle timestep. The frozen
# n+1/2 moments are scratch recomputed every step, so they are not stored.
# ponytail: stdlib Serialization (no extra dep). Swap for HDF5 when cross-version
# or cross-language portability is needed.

using Serialization

const _CHECKPOINT_REQUIRED_FIELDS = (
    :D,
    :T,
    :ncell,
    :L,
    :model,
    :shape,
    :x,
    :v,
    :weight,
    :id,
    :tag,
    :q,
    :m,
    :B,
    :E,
    :time,
    :step,
    :last_dt,
)

function _validate_checkpoint_container(s)
    s isa NamedTuple ||
        throw(ArgumentError("checkpoint file does not contain a valid checkpoint container"))
    missing = Symbol[k for k in _CHECKPOINT_REQUIRED_FIELDS if !hasproperty(s, k)]
    if !isempty(missing)
        missing_list = join(string.(missing), ", ")
        throw(ArgumentError("checkpoint is missing required fields: $missing_list"))
    end
    return nothing
end

function _validate_checkpoint_vector(name::Symbol, a, N::Int, ::Type{T}) where {T}
    a isa AbstractVector ||
        throw(ArgumentError("checkpoint $(name) is $(typeof(a)), expected AbstractVector{$T}"))
    eltype(a) === T || throw(
        ArgumentError(
            "checkpoint $(name) eltype $(eltype(a)) ≠ $T (would silently convert on restore)",
        ),
    )
    length(a) == N ||
        throw(ArgumentError("checkpoint $(name) length $(length(a)) ≠ particle count $N"))
    _check_particle_vector_axes(name, a, N)
    return nothing
end

function _validate_checkpoint_array(name::Symbol, a, reference::AbstractArray{T}) where {T}
    a isa AbstractArray ||
        throw(ArgumentError("checkpoint $(name) is $(typeof(a)), expected AbstractArray{$T}"))
    eltype(a) === T || throw(
        ArgumentError(
            "checkpoint $(name) eltype $(eltype(a)) ≠ $T (would silently convert on restore)",
        ),
    )
    axes(a) == axes(reference) ||
        throw(ArgumentError("checkpoint $(name) axes $(axes(a)) ≠ $(axes(reference))"))
    return nothing
end

function _validate_checkpoint_scalar_state(s, ::Type{T}) where {T}
    s.q isa T || throw(
        ArgumentError("checkpoint q type $(typeof(s.q)) ≠ $T (would silently convert on restore)"),
    )
    s.m isa T || throw(
        ArgumentError("checkpoint m type $(typeof(s.m)) ≠ $T (would silently convert on restore)"),
    )
    s.time isa T || throw(
        ArgumentError(
            "checkpoint time type $(typeof(s.time)) ≠ $T (would silently convert on restore)",
        ),
    )
    s.step isa Int || throw(
        ArgumentError(
            "checkpoint step type $(typeof(s.step)) ≠ Int (would fail or convert on restore)",
        ),
    )
    s.step >= 0 || throw(ArgumentError("checkpoint step must be nonnegative"))
    s.last_dt isa T || throw(ArgumentError("checkpoint last_dt type $(typeof(s.last_dt)) ≠ $T"))
    isfinite(s.last_dt) && s.last_dt >= zero(T) ||
        throw(ArgumentError("checkpoint last_dt must be finite and nonnegative"))
    (s.step == 0 ? iszero(s.last_dt) : s.last_dt > zero(T)) || throw(
        ArgumentError("checkpoint last_dt must be zero before stepping and positive afterward"),
    )
    return nothing
end

function _validate_checkpoint_particle_state(s, ::Val{D}) where {D}
    s.x isa Tuple || throw(ArgumentError("checkpoint positions must be stored as a tuple"))
    s.v isa Tuple || throw(ArgumentError("checkpoint velocities must be stored as a tuple"))
    length(s.x) == D ||
        throw(ArgumentError("checkpoint has $(length(s.x)) position arrays, expected $D"))
    length(s.v) == 3 ||
        throw(ArgumentError("checkpoint has $(length(s.v)) velocity arrays, expected 3"))
    s.weight isa AbstractVector ||
        throw(ArgumentError("checkpoint weight is $(typeof(s.weight)), expected AbstractVector"))
    N = length(s.weight)
    T = s.T
    _validate_checkpoint_vector(:weight, s.weight, N, T)
    for d = 1:D
        _validate_checkpoint_vector(Symbol(:x, d), s.x[d], N, T)
    end
    for c = 1:3
        _validate_checkpoint_vector(Symbol(:v, c), s.v[c], N, T)
    end
    _validate_checkpoint_vector(:id, s.id, N, UInt64)
    _validate_checkpoint_vector(:tag, s.tag, N, UInt32)
    return nothing
end

function _validate_checkpoint_field_state(s, st)
    s.B isa Tuple || throw(ArgumentError("checkpoint B fields must be stored as a tuple"))
    s.E isa Tuple || throw(ArgumentError("checkpoint E fields must be stored as a tuple"))
    length(s.B) == 3 || throw(ArgumentError("checkpoint has $(length(s.B)) B arrays, expected 3"))
    length(s.E) == 3 || throw(ArgumentError("checkpoint has $(length(s.E)) E arrays, expected 3"))
    for c = 1:3
        _validate_checkpoint_array(Symbol(:B, c), s.B[c], st.fields.B[c])
        _validate_checkpoint_array(Symbol(:E, c), s.E[c], st.fields.E[c])
    end
    return nothing
end

function _validate_checkpoint_state(s, st::HybridStepper{D,T}) where {D,T}
    _validate_checkpoint_container(s)
    s.D == D || throw(ArgumentError("checkpoint dimension $(s.D) ≠ $D"))
    s.T == T || throw(
        ArgumentError(
            "checkpoint eltype $(s.T) ≠ $T (would silently convert and break bitwise restart)",
        ),
    )
    s.ncell == st.g.n || throw(ArgumentError("checkpoint grid $(s.ncell) ≠ $(st.g.n)"))
    s.L == st.g.L || throw(ArgumentError("checkpoint box lengths $(s.L) ≠ $(st.g.L)"))
    typeof(s.model) === typeof(st.model) ||
        throw(ArgumentError("checkpoint model type $(typeof(s.model)) ≠ $(typeof(st.model))"))
    typeof(s.shape) === typeof(st.shape) ||
        throw(ArgumentError("checkpoint shape type $(typeof(s.shape)) ≠ $(typeof(st.shape))"))
    _validate_checkpoint_scalar_state(s, T)
    _validate_checkpoint_particle_state(s, Val(D))
    _validate_checkpoint_field_state(s, st)
    return nothing
end

function _checkpoint_state(st::HybridStepper{D,T}, ps::ParticleSet{D,T}) where {D,T}
    return (
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        shape = st.shape,
        x = ps.x,
        v = ps.v,
        weight = ps.weight,
        id = ps.id,
        tag = ps.tag,
        q = ps.q,
        m = ps.m,
        B = st.fields.B,
        E = st.fields.E,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
    )
end

function _validate_checkpoint_destination(ps::ParticleSet, N::Int)
    for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
        if length(a) != N && !applicable(resize!, a, N)
            throw(
                ArgumentError(
                    "checkpoint requires $N particles, but destination $(typeof(a)) cannot resize from $(length(a))",
                ),
            )
        end
    end
    return nothing
end

function _restore_checkpoint_state!(st::HybridStepper{D,T}, ps::ParticleSet{D,T}, s) where {D,T}
    _validate_checkpoint_state(s, st)
    N = length(s.weight)
    _validate_checkpoint_destination(ps, N)
    _resize_hybrid_particle_workspaces!(st, N)
    st.model = s.model
    st.shape = s.shape
    for d = 1:D
        length(ps.x[d]) == N || resize!(ps.x[d], N)
        copyto!(ps.x[d], s.x[d])
    end
    for c = 1:3
        length(ps.v[c]) == N || resize!(ps.v[c], N)
        copyto!(ps.v[c], s.v[c])
    end
    length(ps.weight) == N || resize!(ps.weight, N)
    copyto!(ps.weight, s.weight)
    length(ps.id) == N || resize!(ps.id, N)
    copyto!(ps.id, s.id)
    length(ps.tag) == N || resize!(ps.tag, N)
    copyto!(ps.tag, s.tag)
    ps.q = s.q
    ps.m = s.m
    for c = 1:3
        copyto!(st.fields.B[c], s.B[c])
        copyto!(st.fields.E[c], s.E[c])
    end
    st.time[] = s.time
    st.step[] = s.step
    st.last_dt[] = s.last_dt
    return st
end

"""
    save_checkpoint(path, stepper, ps)

Serialize the restartable state (model, shape, particles, B, carried E, time,
step, and previous timestep for the staggered particle velocities) to `path`.
"""
function save_checkpoint(
    path::AbstractString,
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
) where {D,T}
    state = _checkpoint_state(st, ps)
    _validate_checkpoint_state(state, st)
    serialize(path, state)
    return path
end

"""
    load_checkpoint!(stepper, ps, path)

Restore a checkpoint into a `stepper`/`ps` of matching dimension, precision,
grid, and concrete model/shape types (particle arrays are resized to match).
Destination arrays that do not support resizing must already have the saved
particle count; writable views are restored in place.
The checkpoint replaces the destination's model and shape parameters, including
custom mutable parameters, so restored evolution uses the saved physics.
Checkpoints without configuration or timestep-history fields are rejected. After loading, `step!`
continues bitwise-identically to the original run with the same stepping inputs.
"""
function load_checkpoint!(
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    path::AbstractString,
) where {D,T}
    return _restore_checkpoint_state!(st, ps, deserialize(path))
end
