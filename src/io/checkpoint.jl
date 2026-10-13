# checkpoint.jl — checkpoint/restart of the deterministic simulation state.
#
# Stepping (step!) uses no RNG, so the state needed for a bitwise-identical
# restart is: model and shape, particle phase space + provenance, canonical B,
# carried E, time/step counters, and the previous particle timestep. The frozen
# n+1/2 moments are scratch recomputed every step, so they are not stored.
# CAM-CL also retains its two magnetic copies, prior output, and subcycle phase.
# Electrostatic PIC retains E, density, background density, and priming/timestep
# state; the per-particle field gather buffers are scratch.
# Electromagnetic PIC also retains each species separately, its fast/full-step
# phase, and the charge snapshots/current used by continuity diagnostics.
# ponytail: stdlib Serialization (no extra dep). Swap for HDF5 when cross-version
# or cross-language portability is needed.

using Serialization

const _HybridCheckpointStepper{D,T} = Union{HybridStepper{D,T},CAMCLStepper{D,T}}
const _ElectrostaticCheckpointStepper = Union{Electrostatic1D,ElectrostaticPIC}
const _ElectromagneticCheckpointStepper = Union{EMPIC1D,EMPIC}
const _CheckpointStepper = Union{
    _HybridCheckpointStepper,
    _ElectrostaticCheckpointStepper,
    _ElectromagneticCheckpointStepper,
}

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

const _EMPIC_CHECKPOINT_REQUIRED_FIELDS = (
    :integrator,
    :D,
    :T,
    :ncell,
    :L,
    :n0,
    :c,
    :shape,
    :relativistic,
    :mobile,
    :mi,
    :n_sub,
    :electrons,
    :ions,
    :fields,
    :time,
    :step,
    :last_dt,
    :last_sub_dt,
)

const _CHECKPOINT_PARTICLE_FIELDS = (:x, :v, :weight, :id, :tag, :q, :m)

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
    _validate_checkpoint_time_state(s, T)
    return nothing
end

function _validate_checkpoint_time_state(s, ::Type{T}) where {T}
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

_validate_checkpoint_particle_state(s, dims::Val) =
    _validate_checkpoint_particle_state(s, dims, s.T)

function _validate_checkpoint_particle_state(s, ::Val{D}, ::Type{T}) where {D,T}
    s.x isa Tuple || throw(ArgumentError("checkpoint positions must be stored as a tuple"))
    s.v isa Tuple || throw(ArgumentError("checkpoint velocities must be stored as a tuple"))
    length(s.x) == D ||
        throw(ArgumentError("checkpoint has $(length(s.x)) position arrays, expected $D"))
    length(s.v) == 3 ||
        throw(ArgumentError("checkpoint has $(length(s.v)) velocity arrays, expected 3"))
    s.weight isa AbstractVector ||
        throw(ArgumentError("checkpoint weight is $(typeof(s.weight)), expected AbstractVector"))
    N = length(s.weight)
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

function _checkpoint_species(s, dims::Val{D}, ::Type{T}) where {D,T}
    _validate_checkpoint_container(s, _CHECKPOINT_PARTICLE_FIELDS)
    _validate_checkpoint_particle_state(s, dims, T)
    _validate_checkpoint_particle_scalars(s, T)
    # A validated wrapper reuses the mover's species and relativistic checks
    # without copying particle arrays or converting saved scalar types.
    return ParticleSet{D,T}(s.x, s.v, s.weight, s.id, s.tag, s.q, s.m)
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

_empic_checkpoint_integrator(::EMPIC1D) = :empic1d
_empic_checkpoint_integrator(::EMPIC) = :empic
_empic_checkpoint_fields(st::EMPIC1D) = (
    Ex = st.Ex,
    Ey = st.Ey,
    Bz = st.Bz,
    Jx = st.Jx,
    Jy = st.Jy,
    rho_n = st.rho_n,
    rho_np1 = st.rho_np1,
)
_empic_checkpoint_fields(st::EMPIC) =
    (E = st.E, B = st.B, J = st.J, rho_n = st.rho_n, rho_np1 = st.rho_np1)

function _checkpoint_state(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet} = nothing,
)
    return (
        integrator = _empic_checkpoint_integrator(st),
        D = length(st.g.n),
        T = typeof(st.n0),
        ncell = st.g.n,
        L = st.g.L,
        n0 = st.n0,
        c = st.c,
        shape = st.shape,
        relativistic = st.relativistic,
        mobile = st.mobile,
        mi = st.mi,
        n_sub = st.n_sub,
        electrons = _checkpoint_particle_state(electrons),
        ions = ions === nothing ? nothing : _checkpoint_particle_state(ions),
        fields = _empic_checkpoint_fields(st),
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        last_sub_dt = st.last_sub_dt[],
    )
end

function _validate_checkpoint_state(s, st::_ElectromagneticCheckpointStepper)
    _validate_checkpoint_container(s, _EMPIC_CHECKPOINT_REQUIRED_FIELDS)
    _validate_checkpoint_grid(s, st.g)
    expected = _empic_checkpoint_integrator(st)
    s.integrator === expected ||
        throw(ArgumentError("checkpoint integrator $(s.integrator) ≠ $expected"))
    typeof(s.shape) === typeof(st.shape) ||
        throw(ArgumentError("checkpoint shape type $(typeof(s.shape)) ≠ $(typeof(st.shape))"))
    T = s.T
    for name in (:n0, :c, :mi)
        value = getproperty(s, name)
        value isa T || throw(ArgumentError("checkpoint $name type $(typeof(value)) ≠ $T"))
    end
    _require_finite_nonnegative_real("checkpoint n0", s.n0, T)
    _require_finite_positive_real("checkpoint c", s.c, T)
    _require_finite_positive_real("checkpoint mi", s.mi, T)
    s.relativistic isa Bool || throw(ArgumentError("checkpoint relativistic must be a Bool"))
    s.mobile isa Bool || throw(ArgumentError("checkpoint mobile must be a Bool"))
    s.mobile == st.mobile ||
        throw(ArgumentError("checkpoint ion mobility does not match destination"))
    s.n_sub isa Int && s.n_sub >= 1 ||
        throw(ArgumentError("checkpoint n_sub must be a positive Int"))
    _validate_checkpoint_time_state(s, T)
    isfinite(s.time) && s.time >= zero(T) ||
        throw(ArgumentError("checkpoint time must be finite and nonnegative"))
    s.last_sub_dt isa T ||
        throw(ArgumentError("checkpoint last_sub_dt type $(typeof(s.last_sub_dt)) ≠ $T"))
    isfinite(s.last_sub_dt) && zero(T) <= s.last_sub_dt <= s.last_dt ||
        throw(ArgumentError("checkpoint last_sub_dt must lie between zero and last_dt"))
    (s.step == 0 ? iszero(s.last_sub_dt) : s.last_sub_dt > zero(T)) || throw(
        ArgumentError("checkpoint last_sub_dt must be zero before stepping and positive afterward"),
    )
    dims = Val(s.D)
    electrons = _checkpoint_species(s.electrons, dims, T)
    _require_empic_electrons(electrons)
    s.relativistic && _require_subluminal_velocities(electrons, s.c, "checkpoint electrons")
    if s.mobile
        s.ions === nothing && throw(ArgumentError("mobile checkpoint requires ion particles"))
        ions = _checkpoint_species(s.ions, dims, T)
        _require_empic_ions(ions, s.mi)
        s.relativistic && _require_subluminal_velocities(ions, s.c, "checkpoint ions")
        _particle_sets_mightalias(electrons, ions) &&
            throw(ArgumentError("checkpoint electron and ion storage must be independent"))
    else
        s.ions === nothing ||
            throw(ArgumentError("immobile checkpoint cannot contain ion particles"))
    end
    fields = _empic_checkpoint_fields(st)
    _validate_checkpoint_container(s.fields, keys(fields))
    for (name, reference) in pairs(fields)
        if reference isa Tuple
            _validate_checkpoint_field_components(name, getproperty(s.fields, name), reference)
        else
            _validate_checkpoint_array(name, getproperty(s.fields, name), reference)
        end
    end
    return nothing
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
        for (destination, saved) in ((st.Blead, s.Blead), (st.Blag, s.Blag), (st.Bout, s.Bout)),
            c = 1:3

            copyto!(destination[c], saved[c])
        end
        st.cl_h[] = s.cl_h
        st.cl_live[] = s.cl_live
    end
    return st
end

function _resize_empic_checkpoint_workspaces!(st::EMPIC1D, ne, ni)
    _ensure_electron_buffers!(st, ne)
    _ensure_ion_buffers!(st, ni)
    return nothing
end

function _resize_empic_checkpoint_workspaces!(st::EMPIC, ne, ni)
    _ensure_empic_electron_buffers!(st, ne)
    _ensure_empic_ion_buffers!(st, ni)
    return nothing
end

_restore_checkpoint_state!(st::_ElectromagneticCheckpointStepper, electrons::ParticleSet, s) =
    _restore_checkpoint_state!(st, electrons, nothing, s)

function _restore_checkpoint_state!(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet},
    s,
)
    _validate_checkpoint_state(s, st)
    D, T = s.D, s.T
    electrons isa ParticleSet{D,T} || throw(
        ArgumentError("checkpoint electrons do not match destination dimension and precision"),
    )
    if ions === nothing
        s.mobile && throw(
            ArgumentError(
                "mobile checkpoint requires destination ions of matching dimension and precision",
            ),
        )
    else
        s.mobile || throw(ArgumentError("immobile checkpoint does not restore ion particles"))
        ions isa ParticleSet{D,T} || throw(
            ArgumentError(
                "mobile checkpoint requires destination ions of matching dimension and precision",
            ),
        )
        _particle_sets_mightalias(electrons, ions) &&
            throw(ArgumentError("destination electron and ion storage must be independent"))
    end
    ne = length(s.electrons.weight)
    ni = ions === nothing ? 0 : length(s.ions.weight)
    _validate_checkpoint_destination(electrons, ne)
    ions === nothing || _validate_checkpoint_destination(ions, ni)
    _resize_empic_checkpoint_workspaces!(st, ne, ni)
    for name in (:n0, :c, :shape, :relativistic, :mi, :n_sub)
        setproperty!(st, name, getproperty(s, name))
    end
    _restore_checkpoint_particles!(electrons, s.electrons)
    ions === nothing || _restore_checkpoint_particles!(ions, s.ions)
    fields = _empic_checkpoint_fields(st)
    saved_fields = NamedTuple{keys(fields)}(s.fields)
    map(values(fields), values(saved_fields)) do destination, saved
        if destination isa Tuple
            map(copyto!, destination, saved)
        else
            copyto!(destination, saved)
        end
    end
    st.time[] = s.time
    st.step[] = s.step
    st.last_dt[] = s.last_dt
    st.last_sub_dt[] = s.last_sub_dt
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
    save_checkpoint(path, empic, electrons, ions)

Serialize the restartable state (model, shape, particles, B, carried E, time,
step, and previous timestep for the staggered particle velocities) to `path`.
Supports [`HybridStepper`](@ref) and [`CAMCLStepper`](@ref); CAM-CL checkpoints
also preserve the persistent cyclic-leapfrog field copies and subcycle phase.
For [`Electrostatic1D`](@ref) and [`ElectrostaticPIC`](@ref), save the electron
particles, electric field, number density, background density, and priming/timestep
state. These electrostatic solvers use CIC and an immobile ion background.
[`EMPIC1D`](@ref) and [`EMPIC`](@ref) preserve their physics configuration,
electric/magnetic fields, current and charge snapshots, and full/fast timestep
history. Supply the ion `ParticleSet` exactly when `mobile=true`; electron and
ion storage must be independent. Particle workspaces are resized on restore;
scratch values are recomputed during stepping.
"""
function save_checkpoint(path::AbstractString, st::_CheckpointStepper, ps::ParticleSet)
    return _save_checkpoint(path, st, _checkpoint_state(st, ps))
end

function save_checkpoint(
    path::AbstractString,
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet},
)
    return _save_checkpoint(path, st, _checkpoint_state(st, electrons, ions))
end

function _save_checkpoint(path::AbstractString, st::_CheckpointStepper, state)
    _validate_checkpoint_state(state, st)
    serialize(path, state)
    return path
end

"""
    load_checkpoint!(stepper, ps, path)
    load_checkpoint!(empic, electrons, ions, path)

Restore a checkpoint into a `stepper`/`ps` of matching dimension, precision,
grid, and integrator (particle arrays are resized to match).
Destination arrays that do not support resizing must already have the saved
particle count; writable views are restored in place.
Hybrid destinations require matching concrete model/shape types; the checkpoint
replaces their parameters, including custom mutable parameters, so restored
evolution uses the saved physics.
Electrostatic destinations must have the same background density `n0`, which
is immutable; their saved fields, density, and leapfrog phase are restored directly.
Electromagnetic destinations require the same ion mobility and concrete shape
type; saved physics parameters replace the destination values. Supply independent
electron and ion destinations for mobile runs. Both species are validated before
either is changed.
Checkpoints without integrator, configuration, or timestep-history fields are rejected.
After loading, `step!`, `step_camcl!`, `step_espic!`, or `step_empic!` continues bitwise-identically to the
original run with the same stepping inputs.
"""
function load_checkpoint!(st::_CheckpointStepper, ps::ParticleSet, path::AbstractString)
    return _restore_checkpoint_state!(st, ps, deserialize(path))
end

function load_checkpoint!(
    st::_ElectromagneticCheckpointStepper,
    electrons::ParticleSet,
    ions::Union{Nothing,ParticleSet},
    path::AbstractString,
)
    return _restore_checkpoint_state!(st, electrons, ions, deserialize(path))
end
