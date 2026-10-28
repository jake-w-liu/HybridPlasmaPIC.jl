# Boris.jl — Boris pusher (from particles.jl)

# The plain Boris formula in the caller's arithmetic and historical operation
# order. It is shared by the guarded CPU kick, its wide repair and the Metal
# kernel; the intermediates are read by the CPU range guard only.
@inline function _boris_kick_parts(
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    qm::T,
    dt::T,
) where {T}
    h = qm * dt / 2
    hEx = h * Ex
    hEy = h * Ey
    hEz = h * Ez
    vmx = vx + hEx
    vmy = vy + hEy
    vmz = vz + hEz
    tx = h * Bx
    ty = h * By
    tz = h * Bz
    t2 = tx * tx + ty * ty + tz * tz
    f = 2 / (1 + t2)
    sx = f * tx
    sy = f * ty
    sz = f * tz
    # v' = v⁻ + v⁻ × t
    vpx = vmx + (vmy * tz - vmz * ty)
    vpy = vmy + (vmz * tx - vmx * tz)
    vpz = vmz + (vmx * ty - vmy * tx)
    # v⁺ = v⁻ + v' × s
    vnx = vmx + (vpy * sz - vpz * sy)
    vny = vmy + (vpz * sx - vpx * sz)
    vnz = vmz + (vpx * sy - vpy * sx)
    return (vnx + hEx, vny + hEy, vnz + hEz), (h, hEx, hEy, hEz, tx, ty, tz, t2, f)
end

# Native-only Boris kick: no range guard or fallback (device kernels).
@inline _boris_kick_native(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt) =
    first(_boris_kick_parts(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt))

# Native device contract. Device pushes evaluate _boris_kick_native without a
# wide fallback, and Metal flushes subnormals to zero. Before a uniform-field
# launch the host checks the velocity-independent part of the range guard:
# h = qm*dt/2 finite and normal (or an exact zero from a zero factor), h*E
# finite, and t = h*B with finite |t|² and f = 2/(1 + |t|²) >= floatmin.
# Velocity-dependent range loss on the device (subnormal or overflowing
# velocity terms, drifts) is unsupported; use CPU storage for such states.
function _require_native_uniform_boris(qm::T, dt::T, E::NTuple{3,T}, B::NTuple{3,T}) where {T}
    z = zero(T)
    _, (h, hEx, hEy, hEz, _, _, _, t2, f) = _boris_kick_parts(z, z, z, E..., B..., qm, dt)
    _kinematics_product_in_range(h, qm, dt) &
    isfinite(hEx) &
    isfinite(hEy) &
    isfinite(hEz) &
    isfinite(t2) &
    (f >= floatmin(T)) && return nothing
    throw(
        ArgumentError(
            "native device Boris push requires finite normal q/m*dt/2, finite q/m*dt*E/2 and |q/m*dt*B/2|² well inside the $T range; the device has no range fallback",
        ),
    )
end

# Range guard over the native intermediates. With h normal and no overflow, each
# underflowing product adds absolute error at most floatmin*eps/2, amplified at
# most by (1 + tau) through the rotation; the scale threshold keeps the total
# below eps*S/16. f = 2/(1 + t2) >= floatmin excludes t2 near overflow, where f
# and s lose precision. An intermediate overflow reaches an output as Inf or
# NaN because only f divides. Every term is combined with `&` (branch-free).
@inline function _boris_kick_in_range(
    out::NTuple{3,T},
    parts,
    vx::T,
    vy::T,
    vz::T,
    qm::T,
    dt::T,
) where {T}
    h, hEx, hEy, hEz, tx, ty, tz, t2, f = parts
    S = max(abs(vx), abs(vy), abs(vz), abs(hEx), abs(hEy), abs(hEz))
    tau = max(abs(tx), abs(ty), abs(tz))
    return isfinite(h) &
           ((abs(h) >= floatmin(T)) | (iszero(h) & (iszero(qm) | iszero(dt)))) &
           isfinite(t2) &
           (f >= floatmin(T)) &
           isfinite(out[1]) &
           isfinite(out[2]) &
           isfinite(out[3]) &
           (iszero(S) | (S >= floatmin(T) * T(256) * (one(T) + tau)))
end

"""
    boris_kick(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt) -> (vx′, vy′, vz′)

Boris velocity update over `dt` for charge-to-mass ratio `qm`: half electric
kick, magnetic rotation with `t = (qm dt/2) B`, half electric kick.

Ordinary inputs use the plain formula in its historical operation order. For
`Float16`, `Float32` and `Float64`, a guard over its intermediates detects
overflow or underflow (for example `qm*dt/2` overflowing in a zero field, an
overflowing `|t|²`, or an underflowed `qm*dt` that would drop a representable
kick). Only those evaluations repeat the same expression wider — `Float32` for
`Float16`, `Float64` for `Float32`, exact `Rational{BigInt}` for `Float64` —
rounded once. The result is within a few ulps of `max(|v|, |qm dt E/2|)`;
componentwise accuracy of a nearly cancelling component is not certified. A
result that exceeds the type is returned as `±Inf`, and non-finite inputs keep
IEEE propagation. Other number types use the plain formula.
"""
@inline function boris_kick(
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    qm::T,
    dt::T,
) where {T}
    out, parts = _boris_kick_parts(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
    T <: _RangeGuardedFloat || return out
    _boris_kick_in_range(out, parts, vx, vy, vz, qm, dt) && return out
    return _boris_kick_repair(out, vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
end

# The same expression in the wider format, rounded once. For M = floatmax(T)
# every intermediate is below 8*M^6 (|t| < M^3, |v⁻| < M^3, |v′| < 2*M^6 and
# |s| <= 1): 2^771 for Float32 in Float64 and 2^99 for Float16 in Float32, so
# nothing overflows. Inputs are multiples of nextfloat(zero(T)), so h, h*E, t
# and s are zero or normal in the wide format. Underflow is possible only in
# the rotation products v⁻×t (Float16) and v′×s and in |t|² (Float16, a
# relative 2^-148 effect on f). Those products reach the result additively or
# through |s| <= 1, below 2^-1070 (Float64) resp. 2^-145 (Float32) in absolute
# value and far below the smallest subnormal of T, so one final rounding keeps
# the normwise bound. Float64 inputs are evaluated exactly.
@noinline function _boris_kick_repair(
    out::NTuple{3,T},
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    qm::T,
    dt::T,
) where {T}
    inputs = (vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
    all(isfinite, inputs) || return out
    W = _kinematics_wide(T)
    return map(x -> _kinematics_round(T, x), _boris_kick_native(map(W, inputs)...))
end

@inline function _validated_qm(ps::ParticleSet{D,T}) where {D,T<:AbstractFloat}
    q = _require_finite_real("particle charge q", ps.q, T)
    m = _require_finite_positive_real("particle mass m", ps.m, T)
    qm = q / m
    isfinite(qm) || throw(ArgumentError("particle charge-to-mass ratio q/m must be finite"))
    return qm
end

@inline function _validated_uniform_field(
    name::AbstractString,
    field::NTuple{3},
    ::Type{T},
) where {T<:AbstractFloat}
    return ntuple(3) do c
        value = T(field[c])
        isfinite(value) || throw(ArgumentError("$(name)[$c] must be finite"))
        value
    end
end

function _require_subluminal_velocities(
    ps::ParticleSet{D,T},
    c::T,
    context::AbstractString,
) where {D,T<:AbstractFloat}
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        ux = vx[p]
        uy = vy[p]
        uz = vz[p]
        (isfinite(ux) && isfinite(uy) && isfinite(uz)) ||
            throw(ArgumentError("$context: particle $p velocity must be finite"))
        β2, valid = _subluminal_velocity_parts((ux, uy, uz), c)
        valid || throw(
            ArgumentError(
                "$context: particle $p is outside the supported subluminal range (rounded |v|²/c² = $β2)",
            ),
        )
    end
    return nothing
end

@inline function _require_finite_field_sample(
    name::AbstractString,
    value::Real,
    ::Type{T},
) where {T<:AbstractFloat}
    v = T(value)
    isfinite(v) || throw(ArgumentError("$name must contain only finite values"))
    return v
end

"""
    push_uniform!(ps, E::NTuple{3}, B::NTuple{3}, dt)

Advance every particle one Boris leapfrog step in a spatially uniform field:
v^{n-1/2} → v^{n+1/2}, then x^{n} → x^{n+1} = x^n + dt·v^{n+1/2} (only the D
spatial coordinates move). Returns `ps`.

`dt` must be finite and either zero or at least `2*floatmin(T)` in magnitude;
other values throw an `ArgumentError` before any particle changes. Kicks and
drifts follow the normwise range contract of [`boris_kick`](@ref): a drift
whose rounded `x + dt*v` overflows although the exact coordinate is finite is
evaluated wider and rounded once. A fused reduction over the batch certifies the
common case, which then runs the native vectorized loops (bitwise unchanged);
otherwise every particle takes the guarded evaluation.
"""
function push_uniform!(ps::ParticleSet{D,T}, E::NTuple{3}, B::NTuple{3}, dt::Real) where {D,T}
    qm = _validated_qm(ps)
    dtT = _require_resolved_dt("push_uniform!", _require_finite_real("dt", dt, T))
    electric = _validated_uniform_field("E", E, T)
    magnetic = _validated_uniform_field("B", B, T)
    indices = eachindex(ps.weight)
    h = qm * dtT / 2
    K = _abs3(h * electric[1], h * electric[2], h * electric[3])
    emax, bmax = _abs3(electric...), _abs3(magnetic...)
    # S_p = max(|v_p|, K): a kick K above the threshold covers every particle.
    covered = K >= floatmin(T) * T(256) * (one(T) + 2 * (abs(h) * bmax))
    x = all(d -> axes(ps.x[d]) == axes(ps.weight), 1:D) ? ps.x : ()
    vmax, smin, xmax =
        covered ? _uniform_batch_scales(ps.v, x, K, indices, Val(false)) :
        _uniform_batch_scales(ps.v, x, K, indices, Val(true))
    x === () && (xmax = _coordinate_batch_max(ps.x, T))
    kick, vbound = _boris_batch_certificate(h, qm, dtT, vmax, emax, bmax, smin)
    drift = kick && _coordinate_batch_certified(xmax, dtT, vbound)
    if kick
        _uniform_kicks!(ps, electric, magnetic, qm, dtT, Val(false))
    else
        _uniform_kicks!(ps, electric, magnetic, qm, dtT, Val(true))
    end
    drift ? _uniform_drifts!(ps, dtT, Val(false)) : _uniform_drifts!(ps, dtT, Val(true))
    return ps
end

function _uniform_kicks!(ps, electric, magnetic, qm, dtT, guard::Val)
    Ex, Ey, Ez = electric
    Bx, By, Bz = magnetic
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = _kick_eval(guard, vx[p], vy[p], vz[p], Ex, Ey, Ez, Bx, By, Bz, qm, dtT)
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
    end
    return ps
end

function _uniform_drifts!(ps::ParticleSet{D}, dtT, guard::Val) where {D}
    @inbounds for d = 1:D
        xd = ps.x[d]
        vd = ps.v[d]
        for p in eachindex(xd)
            xd[p] = _coordinate_eval(guard, xd[p], vd[p], dtT)
        end
    end
    return ps
end

function _check_gathered_motion_component(component, ps, name, c)
    for d = 1:length(ps.x)
        _particle_arrays_mightalias(component, ps.x[d]) &&
            throw(ArgumentError("$name[$c] must not alias particle position or velocity storage"))
    end
    for k = 1:3
        _particle_arrays_mightalias(component, ps.v[k]) &&
            throw(ArgumentError("$name[$c] must not alias particle position or velocity storage"))
    end
    return nothing
end
function _check_gathered_field_storage(ps::ParticleSet{D}, E, B, xmid) where {D}
    N = nparticles(ps)
    all(length(E[c]) == N for c = 1:3) ||
        throw(DimensionMismatch("E component length must match particle count"))
    all(length(B[c]) == N for c = 1:3) ||
        throw(DimensionMismatch("B component length must match particle count"))
    # Gathered inputs must stay fixed while positions and velocities publish.
    # Read-only weights/metadata and shared E/B inputs need no extra copies.
    # Literal component calls preserve specialization for mixed vector/view tuples.
    _check_gathered_motion_component(E[1], ps, "E", 1)
    _check_gathered_motion_component(E[2], ps, "E", 2)
    _check_gathered_motion_component(E[3], ps, "E", 3)
    _check_gathered_motion_component(B[1], ps, "B", 1)
    _check_gathered_motion_component(B[2], ps, "B", 2)
    _check_gathered_motion_component(B[3], ps, "B", 3)
    if xmid !== nothing
        length(xmid) == D ||
            throw(DimensionMismatch("xmid must have one component per spatial dimension"))
        all(xmid[d] isa AbstractVector for d = 1:D) ||
            throw(ArgumentError("xmid components must be vectors"))
        all(length(xmid[d]) == N for d = 1:D) ||
            throw(DimensionMismatch("xmid component length must match particle count"))
        for d = 1:D
            midpoint = xmid[d]
            _mightalias_particle_storage(midpoint, ps) &&
                throw(ArgumentError("xmid[$d] must not alias particle storage"))
            for c = 1:3
                (
                    _particle_arrays_mightalias(midpoint, E[c]) ||
                    _particle_arrays_mightalias(midpoint, B[c])
                ) && throw(ArgumentError("xmid[$d] must not alias gathered E or B storage"))
            end
            for k = 1:(d-1)
                _particle_arrays_mightalias(midpoint, xmid[k]) &&
                    throw(ArgumentError("xmid[$d] must not alias xmid[$k]"))
            end
        end
    end
    return nothing
end

"""
    push_gathered!(ps, E, B, dt; xmid=nothing)

Advance particles one Boris leapfrog step using per-particle gathered electric
and magnetic fields. `E` and `B` are 3-tuples of particle-length vectors. If
`xmid` is provided, it receives `x^n + (dt/2) v^{n+1/2}` before positions are
advanced to `x^{n+1}`. The `dt` validation and range contract are those of
[`push_uniform!`](@ref). Midpoint components must be writable vectors whose
storage does not alias any particle array, gathered field or other midpoint
component. Conflicts detected by the particle-storage alias predicate raise
`ArgumentError` before mutation on CPU and Metal. Read-only `E` and `B`
components may share storage with one another.
Gathered inputs must not overlap particle positions or velocities, which this
operation writes. Inputs may share unchanged particle weight or metadata storage.
"""
function push_gathered!(
    ps::ParticleSet{D,T},
    E::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    B::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    dt::Real;
    xmid = nothing,
) where {D,T}
    _check_gathered_field_storage(ps, E, B, xmid)
    qm = _validated_qm(ps)
    dtT = _require_resolved_dt("push_gathered!", _require_finite_real("dt", dt, T))
    indices = eachindex(ps.weight)
    hk = qm * dtT / 2
    scales = _boris_batch_scales(ps.v, E, B, hk, indices)
    # A non-finite field sample makes max|E| or max|B| non-finite; the original
    # validation then reports the first offending sample before any mutation.
    isfinite(scales[2]) & isfinite(scales[3]) || _require_finite_gathered_fields(E, B, indices, T)
    kick, vbound = _boris_batch_certificate(hk, qm, dtT, scales...)
    xmax = _coordinate_batch_max(ps.x, T)
    if kick && _coordinate_batch_certified(xmax, dtT, vbound)
        _gathered_push!(ps, E, B, qm, dtT, xmid, Val(false))
    else
        _gathered_push!(ps, E, B, qm, dtT, xmid, Val(true))
    end
    return ps
end

function _require_finite_gathered_fields(E, B, indices, ::Type{T}) where {T}
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    @inbounds for p in indices
        _require_finite_field_sample("E[1]", Ex[p], T)
        _require_finite_field_sample("E[2]", Ey[p], T)
        _require_finite_field_sample("E[3]", Ez[p], T)
        _require_finite_field_sample("B[1]", Bx[p], T)
        _require_finite_field_sample("B[2]", By[p], T)
        _require_finite_field_sample("B[3]", Bz[p], T)
    end
    return nothing
end

function _gathered_push!(ps::ParticleSet{D,T}, E, B, qm, dtT, xmid, guard::Val) where {D,T}
    h = dtT / T(2)
    vx, vy, vz = ps.v
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = _kick_eval(
            guard,
            vx[p],
            vy[p],
            vz[p],
            _require_finite_field_sample("E[1]", Ex[p], T),
            _require_finite_field_sample("E[2]", Ey[p], T),
            _require_finite_field_sample("E[3]", Ez[p], T),
            _require_finite_field_sample("B[1]", Bx[p], T),
            _require_finite_field_sample("B[2]", By[p], T),
            _require_finite_field_sample("B[3]", Bz[p], T),
            qm,
            dtT,
        )
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        for d = 1:D
            xmid !== nothing && (xmid[d][p] = _coordinate_eval(guard, ps.x[d][p], ps.v[d][p], h))
            ps.x[d][p] = _coordinate_eval(guard, ps.x[d][p], ps.v[d][p], dtT)
        end
    end
    return ps
end

# ---------------------------------------------------------------- boundaries
