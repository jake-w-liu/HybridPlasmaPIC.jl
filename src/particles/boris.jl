# Boris.jl — Boris pusher (from particles.jl)

# The plain Boris formula in the caller's arithmetic and historical operation
# order. It is shared by the guarded CPU kick, its wide repair and the Metal
# kernel; the optional compile-time trace supplies the CPU accuracy guard.
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
    ::Val{trace} = Val(false),
    ops::_ScalarOps = _DefaultScalarOps(),
) where {T,trace}
    n1 = _scalar_mul(ops, qm, dt)
    n2 = _scalar_div(ops, n1, 2)
    n3 = _scalar_mul(ops, n2, Ex)
    n4 = _scalar_mul(ops, n2, Ey)
    n5 = _scalar_mul(ops, n2, Ez)
    n6 = _scalar_add(ops, vx, n3)
    n7 = _scalar_add(ops, vy, n4)
    n8 = _scalar_add(ops, vz, n5)
    n9 = _scalar_mul(ops, n2, Bx)
    n10 = _scalar_mul(ops, n2, By)
    n11 = _scalar_mul(ops, n2, Bz)
    n12 = _scalar_mul(ops, n9, n9)
    n13 = _scalar_mul(ops, n10, n10)
    n14 = _scalar_mul(ops, n11, n11)
    if T <: _RangeGuardedFloat
        n15 = _scalar_add(ops, n12, n13)
        n16 = _scalar_add(ops, n15, n14)
    else
        # Generic arithmetic keeps its historical variadic operator dispatch.
        # Only an explicitly requested generic diagnostic trace adds a prefix.
        n15 = trace ? _scalar_add(ops, n12, n13) : nothing
        n16 = n12 + n13 + n14
    end
    n17 = _scalar_add(ops, 1, n16)
    n18 = _scalar_div(ops, 2, n17)
    n19 = _scalar_mul(ops, n18, n9)
    n20 = _scalar_mul(ops, n18, n10)
    n21 = _scalar_mul(ops, n18, n11)
    n22 = _scalar_mul(ops, n7, n11)
    n23 = _scalar_mul(ops, n8, n10)
    n24 = _scalar_sub(ops, n22, n23)
    n25 = _scalar_add(ops, n6, n24)
    n26 = _scalar_mul(ops, n8, n9)
    n27 = _scalar_mul(ops, n6, n11)
    n28 = _scalar_sub(ops, n26, n27)
    n29 = _scalar_add(ops, n7, n28)
    n30 = _scalar_mul(ops, n6, n10)
    n31 = _scalar_mul(ops, n7, n9)
    n32 = _scalar_sub(ops, n30, n31)
    n33 = _scalar_add(ops, n8, n32)
    n34 = _scalar_mul(ops, n29, n21)
    n35 = _scalar_mul(ops, n33, n20)
    n36 = _scalar_sub(ops, n34, n35)
    n37 = _scalar_add(ops, n6, n36)
    n38 = _scalar_mul(ops, n33, n19)
    n39 = _scalar_mul(ops, n25, n21)
    n40 = _scalar_sub(ops, n38, n39)
    n41 = _scalar_add(ops, n7, n40)
    n42 = _scalar_mul(ops, n25, n20)
    n43 = _scalar_mul(ops, n29, n19)
    n44 = _scalar_sub(ops, n42, n43)
    n45 = _scalar_add(ops, n8, n44)
    n46 = _scalar_add(ops, n37, n3)
    n47 = _scalar_add(ops, n41, n4)
    n48 = _scalar_add(ops, n45, n5)
    out = (n46, n47, n48)
    parts = (n2, n3, n4, n5, n9, n10, n11, n16, n18)
    trace && return out,
    parts,
    (
        n1,
        n2,
        n3,
        n4,
        n5,
        n6,
        n7,
        n8,
        n9,
        n10,
        n11,
        n12,
        n13,
        n14,
        n15,
        n16,
        n17,
        n18,
        n19,
        n20,
        n21,
        n22,
        n23,
        n24,
        n25,
        n26,
        n27,
        n28,
        n29,
        n30,
        n31,
        n32,
        n33,
        n34,
        n35,
        n36,
        n37,
        n38,
        n39,
        n40,
        n41,
        n42,
        n43,
        n44,
        n45,
        n46,
        n47,
        n48,
    )
    return out, parts
end

# Native-only Boris kick: no range guard or fallback (device kernels).
@inline _boris_kick_native(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt) =
    first(_boris_kick_parts(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt))

# Native device contract. Device pushes evaluate _boris_kick_native without a
# wide fallback, and Metal flushes subnormals to zero. Before a uniform-field
# launch the host admits original q,m only when qm is normal or q is true zero,
# then checks the velocity-independent part of the range guard:
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

For CPU `Float16`, `Float32` and `Float64`, the historical expression is used
when sufficient range and accuracy bounds prove the result. Unresolved finite
inputs evaluate that same expression exactly as `Rational{BigInt}`. A native
result already within four ulps at the exact Euclidean scale
`max(norm(v),norm(qm*dt*E/2))` keeps its bits when the old range guard admits it;
otherwise the exact result is rounded once. This is an absolute tolerance at
the operation scale, including nearly cancelling output components. Exact
nearest-even overflow returns the corresponding `±Inf`. Nonfinite inputs keep
IEEE propagation; other number types keep the plain formula. The explicit `qm`
argument means that stored number. CPU species pushes retain their original
stored charge and mass separately. Native device arithmetic retains its
separately documented explicit-qm domain.
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
    T <: _RangeGuardedFloat || return _boris_kick_native(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
    return _boris_kick_checked(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
end

@inline function _boris_kick_checked(
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    factor,
    dt::T,
) where {T}
    qm = _boris_native_factor(factor)
    out, parts, nodes = _boris_kick_parts(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt, Val(true))
    # Stored nonfinite operands keep the historical IEEE expression. A species
    # quotient can be nonfinite despite finite original q,m and must go cold.
    all(isfinite, (vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, dt)) || return out
    (!(factor isa _BorisSpeciesFactors) && !isfinite(qm)) && return out
    range = _boris_kick_in_range(out, parts, vx, vy, vz, qm, dt)
    range &&
        _boris_accuracy_scalar(
            out,
            parts,
            nodes,
            (vx, vy, vz),
            (Ex, Ey, Ez),
            (Bx, By, Bz),
            factor,
            dt,
        ) &&
        return out
    return _boris_kick_repair(out, vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, factor, dt, range)
end

# Exact replay uses the one canonical expression. No separately transcribed
# rotation is used for the fallback or its native-retention decision.
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
    factor,
    dt::T,
    retain_native::Bool = false,
) where {T}
    inputs = (vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, dt)
    all(isfinite, inputs) || return out
    qm = _boris_native_factor(factor)
    (!(factor isa _BorisSpeciesFactors) && !isfinite(qm)) && return out
    Q = _exact_rational
    exact, parts = _boris_kick_parts(
        Q(vx),
        Q(vy),
        Q(vz),
        Q(Ex),
        Q(Ey),
        Q(Ez),
        Q(Bx),
        Q(By),
        Q(Bz),
        _boris_exact_factor(factor),
        Q(dt),
    )
    scale2 = max(Q(vx)^2 + Q(vy)^2 + Q(vz)^2, parts[2]^2 + parts[3]^2 + parts[4]^2)
    if retain_native && _boris_exact_native_ok(out, exact, scale2)
        return out
    end
    return map(x -> _from_exact(T, x), exact)
end

function _boris_exact_native_ok(out::NTuple{3,T}, exact, scale2::Rational{BigInt}) where {T}
    Q = Rational{BigInt}
    minimum = exponent(floatmin(T))
    e = minimum
    if !iszero(scale2)
        n, d = numerator(scale2), denominator(scale2)
        b = ndigits(n; base = 2) - ndigits(d; base = 2)
        (b >= 0 ? n < d << b : n << -b < d) && (b -= 1)
        e = max(fld(b, 2), minimum)
    end
    tolerance = 4Q(2)^(e - precision(T) + 1)
    threshold = _exact_rational(floatmax(T)) + Q(2)^(exponent(floatmax(T)) - precision(T))
    return all(1:3) do c
        x = exact[c]
        abs(x) >= threshold ? isinf(out[c]) && signbit(out[c]) == (x < 0) :
        isfinite(out[c]) && abs(_exact_rational(out[c]) - x) <= tolerance
    end
end

@inline function _kick_eval(
    ::Val{false},
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    factor::_BorisSpeciesFactors{T},
    dt::T,
) where {T}
    return _boris_kick_native(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, factor.qm, dt)
end
@inline function _kick_eval(
    ::Val{true},
    vx::T,
    vy::T,
    vz::T,
    Ex::T,
    Ey::T,
    Ez::T,
    Bx::T,
    By::T,
    Bz::T,
    factor::_BorisSpeciesFactors{T},
    dt::T,
) where {T}
    return _boris_kick_checked(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, factor, dt)
end

@inline function _validated_charge_mass(q::Real, m::Real, ::Type{T}) where {T<:AbstractFloat}
    q = _require_finite_real("particle charge q", q, T)
    m = _require_finite_real("particle mass m", m, T)
    _stored_positive(m) || throw(ArgumentError("particle mass m must be finite and positive"))
    return q, m
end
@inline _validated_charge_mass(ps::ParticleSet{D,T}) where {D,T<:AbstractFloat} =
    _validated_charge_mass(ps.q, ps.m, T)
@inline function _validated_boris_factor(ps::ParticleSet{D,T}) where {D,T<:AbstractFloat}
    q, m = _validated_charge_mass(ps)
    qm = q / m
    T <: _RangeGuardedFloat && return _BorisSpeciesFactors(q, m, qm)
    isfinite(qm) || throw(ArgumentError("particle charge-to-mass ratio q/m must be finite"))
    return qm
end

# Native device kernels have no original-factor replay, and their scalar
# quotient argument must survive subnormal flushing. Keep CPU admission separate.
@inline function _validated_native_boris_qm(ps::ParticleSet{D,T}) where {D,T<:_RangeGuardedFloat}
    q, m = _validated_charge_mass(ps)
    _particle_stored_zero(q) && return q
    qm = q / m
    (isfinite(qm) && abs(qm) >= floatmin(T)) ||
        throw(
            ArgumentError("native device Boris push requires finite normal q/m for nonzero charge"),
        )
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
drifts follow the normwise range contract of [`boris_kick`](@ref). The drift
classifies nearest-even overflow from the exact stored operands, retaining
valid finite native bits and repairing nonfinite native results with one exact
rounding. A fused reduction over the batch certifies the
common case, which then runs the native vectorized loops (bitwise unchanged);
otherwise every particle takes the guarded evaluation.
"""
function push_uniform!(ps::ParticleSet{D,T}, E::NTuple{3}, B::NTuple{3}, dt::Real) where {D,T}
    factor = _validated_boris_factor(ps)
    qm = _boris_native_factor(factor)
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
    kick, vbound = _boris_batch_admission(ps, h, factor, dtT, vmax, emax, bmax, smin)
    drift = kick && _coordinate_batch_certified(ps, xmax, dtT, vbound)
    if kick
        _uniform_kicks!(ps, electric, magnetic, factor, dtT, Val(false))
    else
        _uniform_kicks!(ps, electric, magnetic, factor, dtT, Val(true))
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
    factor = _validated_boris_factor(ps)
    qm = _boris_native_factor(factor)
    dtT = _require_resolved_dt("push_gathered!", _require_finite_real("dt", dt, T))
    indices = eachindex(ps.weight)
    hk = qm * dtT / 2
    scales = _boris_batch_scales(ps.v, E, B, hk, indices)
    # A non-finite field sample makes max|E| or max|B| non-finite; the original
    # validation then reports the first offending sample before any mutation.
    isfinite(scales[2]) & isfinite(scales[3]) || _require_finite_gathered_fields(E, B, indices, T)
    kick, vbound = _boris_batch_admission(ps, hk, factor, dtT, scales...)
    xmax = _coordinate_batch_max(ps.x, T)
    if kick && _coordinate_batch_certified(ps, xmax, dtT, vbound)
        _gathered_push!(ps, E, B, factor, dtT, xmid, Val(false))
    else
        _gathered_push!(ps, E, B, factor, dtT, xmid, Val(true))
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
