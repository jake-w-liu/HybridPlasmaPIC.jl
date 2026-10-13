# checkpoint.jl — checkpoint/restart of the deterministic simulation state.
#
# Stepping (step!) uses no RNG, so the state needed for a bitwise-identical
# restart is: model and shape, particle phase space + provenance, canonical B,
# carried E, time/step counters, and the previous particle timestep. The frozen
# n+1/2 moments are scratch recomputed every step, so they are not stored.
# CAM-CL also retains its two magnetic copies, prior output, and subcycle phase.
# Electrostatic PIC retains E, density, background density, and priming/timestep
# state; the per-particle field gather buffers are scratch.
# ponytail: stdlib Serialization (no extra dep). Swap for HDF5 when cross-version
# or cross-language portability is needed.

using Serialization

const _HybridCheckpointStepper{D,T} = Union{HybridStepper{D,T},CAMCLStepper{D,T}}
const _ElectrostaticCheckpointStepper = Union{Electrostatic1D,ElectrostaticPIC}
const _CheckpointStepper = Union{_HybridCheckpointStepper,_ElectrostaticCheckpointStepper}

const _CHECKPOINT_REQUIRED_FIELDS = (
    :integrator,
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

const _ESPIC_CHECKPOINT_REQUIRED_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :n0,
    :x,
    :v,
    :weight,
    :id,
    :tag,
    :q,
    :m,
    :E,
    :ne,
    :primed,
    :last_dt,
)

function _validate_checkpoint_container(s, required = _CHECKPOINT_REQUIRED_FIELDS)
    s isa NamedTuple ||
        throw(ArgumentError("checkpoint file does not contain a valid checkpoint container"))
    missing = Symbol[k for k in required if !hasproperty(s, k)]
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

function _validate_checkpoint_particle_scalars(s, ::Type{T}) where {T}
    s.q isa T || throw(
        ArgumentError("checkpoint q type $(typeof(s.q)) ≠ $T (would silently convert on restore)"),
    )
    s.m isa T || throw(
        ArgumentError("checkpoint m type $(typeof(s.m)) ≠ $T (would silently convert on restore)"),
    )
    return nothing
end

function _validate_checkpoint_scalar_state(s, ::Type{T}) where {T}
    _validate_checkpoint_particle_scalars(s, T)
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

function _validate_checkpoint_field_components(name::Symbol, components, reference)
    components isa Tuple ||
        throw(ArgumentError("checkpoint $name fields must be stored as a tuple"))
    length(components) == 3 ||
        throw(ArgumentError("checkpoint has $(length(components)) $name arrays, expected 3"))
    for c = 1:3
        _validate_checkpoint_array(Symbol(name, c), components[c], reference[c])
    end
    return nothing
end

function _validate_checkpoint_field_state(s, st)
    _validate_checkpoint_field_components(:B, s.B, st.fields.B)
    _validate_checkpoint_field_components(:E, s.E, st.fields.E)
    return nothing
end

function _validate_checkpoint_integrator_state(s, st::_HybridCheckpointStepper{D,T}) where {D,T}
    expected = st isa CAMCLStepper ? :camcl : :rk4
    s.integrator === expected ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ $expected"))
    if st isa CAMCLStepper
        for name in (:Blead, :Blag, :Bout, :cl_h, :cl_live)
            hasproperty(s, name) || throw(ArgumentError("checkpoint is missing $name"))
        end
        for name in (:Blead, :Blag, :Bout)
            _validate_checkpoint_field_components(name, getproperty(s, name), getproperty(st, name))
        end
        s.cl_h isa T || throw(ArgumentError("checkpoint cl_h type $(typeof(s.cl_h)) ≠ $T"))
        s.cl_live isa Bool || throw(ArgumentError("checkpoint cl_live must be a Bool"))
        isfinite(s.cl_h) && s.cl_h >= zero(T) ||
            throw(ArgumentError("checkpoint cl_h must be finite and nonnegative"))
        s.cl_live &&
            iszero(s.cl_h) &&
            throw(ArgumentError("checkpoint live cyclic-leapfrog state requires positive cl_h"))
    end
    return nothing
end

function _validate_checkpoint_grid(s, g::FourierGrid{D,T}) where {D,T}
    s.D == D || throw(ArgumentError("checkpoint dimension $(s.D) ≠ $D"))
    s.T == T || throw(
        ArgumentError(
            "checkpoint eltype $(s.T) ≠ $T (would silently convert and break bitwise restart)",
        ),
    )
    s.ncell == g.n || throw(ArgumentError("checkpoint grid $(s.ncell) ≠ $(g.n)"))
    s.L == g.L || throw(ArgumentError("checkpoint box lengths $(s.L) ≠ $(g.L)"))
    return nothing
end

function _validate_checkpoint_state(s, st::_HybridCheckpointStepper{D,T}) where {D,T}
    _validate_checkpoint_container(s)
    _validate_checkpoint_grid(s, st.g)
    typeof(s.model) === typeof(st.model) ||
        throw(ArgumentError("checkpoint model type $(typeof(s.model)) ≠ $(typeof(st.model))"))
    typeof(s.shape) === typeof(st.shape) ||
        throw(ArgumentError("checkpoint shape type $(typeof(s.shape)) ≠ $(typeof(st.shape))"))
    _validate_checkpoint_scalar_state(s, T)
    _validate_checkpoint_particle_state(s, Val(D))
    _validate_checkpoint_field_state(s, st)
    _validate_checkpoint_integrator_state(s, st)
    return nothing
end

_espic_checkpoint_integrator(::Electrostatic1D) = :espic1d
_espic_checkpoint_integrator(::ElectrostaticPIC) = :espic

function _validate_checkpoint_state(s, st::_ElectrostaticCheckpointStepper)
    _validate_checkpoint_container(s, _ESPIC_CHECKPOINT_REQUIRED_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    expected = _espic_checkpoint_integrator(st)
    s.integrator === expected ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ $expected"))
    T = eltype(st.ne)
    _validate_checkpoint_particle_scalars(s, T)
    s.q == -one(T) && s.m == one(T) ||
        throw(ArgumentError("electrostatic checkpoint requires electrons with q = -1 and m = 1"))
    s.n0 isa T || throw(ArgumentError("checkpoint n0 type $(typeof(s.n0)) ≠ $T"))
    isfinite(s.n0) && s.n0 >= zero(T) ||
        throw(ArgumentError("checkpoint n0 must be finite and nonnegative"))
    isequal(s.n0, st.n0) || throw(ArgumentError("checkpoint background density $(s.n0) ≠ $(st.n0)"))
    s.primed isa Bool || throw(ArgumentError("checkpoint primed must be a Bool"))
    s.last_dt isa T || throw(ArgumentError("checkpoint last_dt type $(typeof(s.last_dt)) ≠ $T"))
    isfinite(s.last_dt) && s.last_dt >= zero(T) ||
        throw(ArgumentError("checkpoint last_dt must be finite and nonnegative"))
    (s.primed ? s.last_dt > zero(T) : iszero(s.last_dt)) ||
        throw(ArgumentError("checkpoint last_dt must be positive exactly when primed"))
    _validate_checkpoint_particle_state(s, Val(length(st.g.n)))
    _validate_checkpoint_array(:ne, s.ne, st.ne)
    if st isa Electrostatic1D
        _validate_checkpoint_array(:E, s.E, st.E)
    else
        _validate_checkpoint_field_components(:E, s.E, st.E)
    end
    return nothing
end

function _checkpoint_particle_state(ps::ParticleSet)
    return (x = ps.x, v = ps.v, weight = ps.weight, id = ps.id, tag = ps.tag, q = ps.q, m = ps.m)
end

function _checkpoint_state(st::_HybridCheckpointStepper{D,T}, ps::ParticleSet{D,T}) where {D,T}
    common = (;
        integrator = st isa CAMCLStepper ? :camcl : :rk4,
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        shape = st.shape,
        _checkpoint_particle_state(ps)...,
        B = st.fields.B,
        E = st.fields.E,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
    )
    if st isa CAMCLStepper
        return merge(
            common,
            (
                Blead = st.Blead,
                Blag = st.Blag,
                Bout = st.Bout,
                cl_h = st.cl_h[],
                cl_live = st.cl_live[],
            ),
        )
    end
    return common
end

function _checkpoint_state(st::_ElectrostaticCheckpointStepper, ps::ParticleSet{D,T}) where {D,T}
    return (;
        integrator = _espic_checkpoint_integrator(st),
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        n0 = st.n0,
        _checkpoint_particle_state(ps)...,
        E = st.E,
        ne = st.ne,
        primed = st.primed[],
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

function _restore_checkpoint_particles!(ps::ParticleSet{D}, s) where {D}
    N = length(s.weight)
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
    return ps
end

function _restore_checkpoint_state!(
    st::_HybridCheckpointStepper{D,T},
    ps::ParticleSet{D,T},
    s,
) where {D,T}
    _validate_checkpoint_state(s, st)
    N = length(s.weight)
    _validate_checkpoint_destination(ps, N)
    _resize_hybrid_particle_workspaces!(st, N)
    st.model = s.model
    st.shape = s.shape
    _restore_checkpoint_particles!(ps, s)
    for c = 1:3
        copyto!(st.fields.B[c], s.B[c])
        copyto!(st.fields.E[c], s.E[c])
    end
    st.time[] = s.time
    st.step[] = s.step
    st.last_dt[] = s.last_dt
    if st isa CAMCLStepper
        for name in (:Blead, :Blag, :Bout), c = 1:3
            copyto!(getproperty(st, name)[c], getproperty(s, name)[c])
        end
        st.cl_h[] = s.cl_h
        st.cl_live[] = s.cl_live
    end
    return st
end

function _restore_checkpoint_state!(
    st::_ElectrostaticCheckpointStepper,
    ps::ParticleSet{D,T},
    s,
) where {D,T}
    _validate_checkpoint_state(s, st)
    s.D == D && s.T === T || throw(
        ArgumentError("checkpoint particles do not match destination dimension and precision"),
    )
    N = length(s.weight)
    _validate_checkpoint_destination(ps, N)
    _resize_espic_particle_workspaces!(st, N)
    _restore_checkpoint_particles!(ps, s)
    copyto!(st.ne, s.ne)
    if st isa Electrostatic1D
        copyto!(st.E, s.E)
    else
        for c = 1:3
            copyto!(st.E[c], s.E[c])
        end
    end
    st.primed[] = s.primed
    st.last_dt[] = s.last_dt
    return st
end

"""
    save_checkpoint(path, stepper, ps)

Serialize the restartable state (model, shape, particles, B, carried E, time,
step, and previous timestep for the staggered particle velocities) to `path`.
Supports [`HybridStepper`](@ref) and [`CAMCLStepper`](@ref); CAM-CL checkpoints
also preserve the persistent cyclic-leapfrog field copies and subcycle phase.
For [`Electrostatic1D`](@ref) and [`ElectrostaticPIC`](@ref), save the electron
particles, electric field, number density, background density, and priming/timestep
state. These electrostatic solvers use CIC and an immobile ion background.
"""
function save_checkpoint(path::AbstractString, st::_CheckpointStepper, ps::ParticleSet)
    state = _checkpoint_state(st, ps)
    _validate_checkpoint_state(state, st)
    serialize(path, state)
    return path
end

"""
    load_checkpoint!(stepper, ps, path)

Restore a checkpoint into a `stepper`/`ps` of matching dimension, precision,
grid, and integrator (particle arrays are resized to match).
Destination arrays that do not support resizing must already have the saved
particle count; writable views are restored in place.
Hybrid destinations require matching concrete model/shape types; the checkpoint
replaces their parameters, including custom mutable parameters, so restored
evolution uses the saved physics.
Electrostatic destinations must have the same background density `n0`, which
is immutable; their saved fields, density, and leapfrog phase are restored directly.
Checkpoints without integrator, configuration, or timestep-history fields are rejected.
After loading, `step!`, `step_camcl!`, or `step_espic!` continues bitwise-identically to the
original run with the same stepping inputs.
"""
function load_checkpoint!(st::_CheckpointStepper, ps::ParticleSet, path::AbstractString)
    return _restore_checkpoint_state!(st, ps, deserialize(path))
end
