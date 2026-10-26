# kinematics_range.jl — normwise range repair shared by the particle kinematics.
#
# Contract (docs/PARTICLE_KINEMATICS_RANGE_PLAN.md):
# - Ordinary inputs keep the native operation order, so results are bitwise
#   those of the plain formulas.
# - A cheap guard over the native intermediates flags intermediate overflow or
#   underflow. Only flagged Float16/Float32/Float64 evaluations repeat the same
#   expression in a wider format and round once to T. Other number types keep
#   native semantics, and non-finite inputs keep IEEE propagation.
# - Accuracy is a few ulps of the operation scale (|v|, |h E| for kicks and
#   |x| + |dt v| for drifts). Componentwise accuracy of a nearly cancelling
#   component is not certified.
# - A wide result that does not fit T is returned as ±Inf for the caller's
#   finiteness policy; leapfrog priming instead raises an ArgumentError.

const _RangeGuardedFloat = Union{Float16,Float32,Float64}

# Wider formats for flagged rational expressions. Each expression comment
# states why its wide evaluation cannot overflow; Float64 is evaluated exactly.
_kinematics_wide(::Type{Float16}) = Float32
_kinematics_wide(::Type{Float32}) = Float64
_kinematics_wide(::Type{Float64}) = Rational{BigInt}

# The single rounding back to T: an IEEE conversion or a correctly rounded
# exact conversion, including subnormal results and overflow to ±Inf.
@inline _kinematics_round(::Type{T}, x::AbstractFloat) where {T} = T(x)
_kinematics_round(::Type{T}, x::Rational{BigInt}) where {T} = _from_exact(T, x)
_kinematics_round(::Type{T}, x::BigFloat) where {T} = _from_exact(T, _exact_rational(x))

# A finite product that is normal, or an exact zero produced by a zero factor.
# An underflowed product of nonzero factors would silently drop a kick.
@inline _kinematics_product_in_range(p::T, a::T, b::T) where {T} =
    isfinite(p) & ((abs(p) >= floatmin(T)) | (iszero(p) & (iszero(a) | iszero(b))))

# Relativistic expressions contain square roots, so their flagged evaluations
# use 256-bit BigFloat. The setprecision/setrounding do-block scopes are
# task-local on the supported Julia versions (ScopedValue based, checked on
# 1.12.7 and 1.13.3); the global MPFR state is never changed. Only arithmetic
# and sqrt are used, which do not touch MPFR's shared constant caches.
const _KINEMATICS_BIGFLOAT_BITS = 256
_kinematics_bigfloat(f) =
    setprecision(BigFloat, _KINEMATICS_BIGFLOAT_BITS) do
        setrounding(f, BigFloat, RoundNearest)
    end

# A rounded gap 1 - β² is well conditioned when β² <= 1/2: the gap
# is at least 1/2, so errors in the positive norm sum are not magnified by
# cancellation. Faster particles resolve a scaled gap or retain exact-gap replay.
@inline _lorentz_gap_in_range(β2) = β2 <= oftype(β2, 0.5)

# The historical Lorentz factor is retained when its subtraction is well
# conditioned. Faster IEEE particles first resolve the scaled gap with the
# shared error-free/double-word arithmetic; an unresolved gap retains the
# existing exact-gap replay. The flag belongs to the factor, not the input
# speed, so a corrected factor may pass the ordinary momentum range guard.
@inline _native_lorentz_factor(β2) = one(β2) / sqrt(one(β2) - β2)

@inline function _conditioned_lorentz_factor(v::NTuple{3,T}, c::T, β2::T) where {T}
    native = _native_lorentz_factor(β2)
    (!(T <: _RangeGuardedFloat) || _lorentz_gap_in_range(β2)) && return native, true
    return _scaled_lorentz_factor(native, v, c)
end

@inline function _lorentz_gap_words(v::NTuple{3,T}, c::T) where {T}
    cf = abs(Float64(c))
    kc = exponent(cf)
    ĉ = ldexp(cf, -kc)
    v̂ = map(x -> ldexp(Float64(x), -kc), v)
    all(isfinite, v̂) || return ĉ, (0.0, 0.0), ntuple(_ -> 0.0, 8)
    c2 = _eft_two_product(ĉ, ĉ)
    p1 = _eft_two_product(v̂[1], v̂[1])
    p2 = _eft_two_product(v̂[2], v̂[2])
    p3 = _eft_two_product(v̂[3], v̂[3])
    squares = (c2[2], c2[1], -p1[2], -p1[1], -p2[2], -p2[1], -p3[2], -p3[1])
    return ĉ, _eft_sum_dd(squares), squares
end

@inline function _lorentz_factor_words(v::NTuple{3,T}, c::T) where {T}
    ĉ, gap, _ = _lorentz_gap_words(v, c)
    # At most a few minimum-subnormal quanta can be lost by scaled products
    # and their low parts. Above 2^-900 this is far below double-word relative
    # precision; sqrt corrections down to roughly 2^-1006 remain normal and
    # gamma is below roughly 2^450. Reject unresolved gaps instead of clipping.
    gap[1] >= 0x1p-900 || return (0.0, 0.0), false
    gamma = _dd_div((ĉ, 0.0), _dd_sqrt(gap))
    return gamma, all(isfinite, gamma)
end

@noinline function _scaled_lorentz_factor(native::T, v::NTuple{3,T}, c::T) where {T}
    words, resolved = _lorentz_factor_words(v, c)
    resolved || return native, false
    gamma = T(words[1] + words[2])
    return gamma, isfinite(gamma)
end

# 1 - |v|²/c² from the stored operands. The rounded test β² < 1 accepts speeds
# within rounding of c, where any rounded gap can vanish or change sign.
function _exact_lorentz_gap(vx, vy, vz, c)
    Q = _exact_rational
    c2 = Q(c)^2
    return (c2 - (Q(vx)^2 + Q(vy)^2 + Q(vz)^2)) / c2
end

# ---------------------------------------------------------------- drift

@inline _particle_coordinate_native(x, v, dt) = x + dt * v

"""
    _particle_coordinate(x, v, dt)

Particle drift `x + dt*v` (pass `dt/2` for a midpoint) in the native operation
order. Only a non-finite result from finite operands is repaired: additions of
subnormals are exact, so an underflowed product costs at most one subnormal
quantum, within the normwise `|x| + |dt v|` scale. In the wide format `dt*v` is
exact (twice T's precision fits the wider significand) and lies far inside its
range, so the repaired value is the sum rounded once in the wide format and
once to `T`. A coordinate that truly exceeds `T` stays `±Inf` for the caller.
"""
@inline function _particle_coordinate(x::T, v::T, dt::T) where {T}
    position = _particle_coordinate_native(x, v, dt)
    T <: _RangeGuardedFloat || return position
    isfinite(position) && return position
    return _particle_coordinate_repair(position, x, v, dt)
end

@noinline function _particle_coordinate_repair(position::T, x::T, v::T, dt::T) where {T}
    (isfinite(x) & isfinite(v) & isfinite(dt)) || return position
    W = _kinematics_wide(T)
    return _kinematics_round(T, _particle_coordinate_native(W(x), W(v), W(dt)))
end

# ---------------------------------------------------------------- Euler half-kicks

# Forward-Euler Lorentz half-kick u ∓ (h*qm)*(E + u×B) in the historical
# operation order: `:prime` backs a loaded velocity up a half step (leapfrog
# priming), `:predict` advances it (carried-E predictor). `hq` is the caller's
# loop-invariant product h*qm. Returns the result and the guarded intermediates.
@inline function _euler_kick_parts(
    ux::T,
    uy::T,
    uz::T,
    ex::T,
    ey::T,
    ez::T,
    bx::T,
    by::T,
    bz::T,
    hq::T,
    ::Val{mode},
) where {T,mode}
    ax = ex + (uy * bz - uz * by)
    ay = ey + (uz * bx - ux * bz)
    az = ez + (ux * by - uy * bx)
    kx = hq * ax
    ky = hq * ay
    kz = hq * az
    out = mode === :prime ? (ux - kx, uy - ky, uz - kz) : (ux + kx, uy + ky, uz + kz)
    return out, (ax, ay, az, kx, ky, kz)
end

# With h*qm normal (or an exact zero from a zero factor) and no overflow, each
# underflowing product adds absolute error at most floatmin*eps/2, and the
# products inside E + u×B are amplified by |h*qm|. Their total stays below
# floatmin*eps*(1 + |h*qm|), which the scale threshold keeps under eps*S/256.
@inline function _euler_kick_in_range(out::NTuple{3,T}, parts, u, h::T, qm::T, hq::T) where {T}
    T <: _RangeGuardedFloat || return true
    ax, ay, az, kx, ky, kz = parts
    S = max(abs(u[1]), abs(u[2]), abs(u[3]), abs(kx), abs(ky), abs(kz))
    return _kinematics_product_in_range(hq, h, qm) &
           isfinite(ax) &
           isfinite(ay) &
           isfinite(az) &
           isfinite(out[1]) &
           isfinite(out[2]) &
           isfinite(out[3]) &
           (iszero(S) | (S >= floatmin(T) * T(256) * (one(T) + abs(hq))))
end

@inline function _euler_kick(u::NTuple{3,T}, e, b, h::T, qm::T, hq::T, ::Val{mode}) where {T,mode}
    out, parts = _euler_kick_parts(u..., e..., b..., hq, Val(mode))
    _euler_kick_in_range(out, parts, u, h, qm, hq) && return out
    return _euler_kick_repair(out, u, e, b, h, qm, Val(mode))
end

# The same expression in the wider format. Every term is below
# 4*floatmax(T)^4 (2^514 for Float32 in Float64, 2^66 for Float16 in Float32).
# Inputs are multiples of nextfloat(zero(T)), so h*qm and the u×B products are
# exact and either zero or at least 2^-298 (2^-48), E + u×B is either zero or
# at least 2^-350 (2^-71), and the kick product is either zero or at least
# 2^-648 (2^-119), all normal; a sum landing in the subnormal range is exact.
# The wide evaluation therefore has no overflow and no lossy underflow.
# Missing planar components are exact zeros and leave the remaining terms exact.
@noinline function _euler_kick_repair(
    out::NTuple{3,T},
    u,
    e,
    b,
    h::T,
    qm::T,
    ::Val{mode},
) where {T,mode}
    all(isfinite, (u..., e..., b..., h, qm)) || return out
    W = _kinematics_wide(T)
    wide, _ = _euler_kick_parts(map(W, u)..., map(W, e)..., map(W, b)..., W(h) * W(qm), Val(mode))
    kicked = map(x -> _kinematics_round(T, x), wide)
    mode === :prime &&
        !all(isfinite, kicked) &&
        throw(
            ArgumentError(
                "leapfrog priming: the half-step backed-up velocity exceeds the range of $T; rescale the units or reduce dt",
            ),
        )
    return kicked
end

# Electrostatic priming of one velocity component, v - (h*qm)*E, in the
# historical order. Only the kick product can underflow; its error of at most
# floatmin*eps/2 stays below eps*S/512 above the scale threshold. In the wider
# format h*qm is exact and both products are zero or normal (at least 2^-447
# for Float32 in Float64, 2^-72 for Float16 in Float32) and below 2*floatmax^3.
@inline function _es_prime_component(v::T, E::T, h::T, qm::T, hq::T) where {T}
    k = hq * E
    out = v - k
    T <: _RangeGuardedFloat || return out
    S = max(abs(v), abs(k))
    _kinematics_product_in_range(hq, h, qm) &
    isfinite(out) &
    (iszero(S) | (S >= floatmin(T) * T(256))) && return out
    return _es_prime_component_repair(out, v, E, h, qm)
end

@noinline function _es_prime_component_repair(out::T, v::T, E::T, h::T, qm::T) where {T}
    all(isfinite, (v, E, h, qm)) || return out
    W = _kinematics_wide(T)
    primed = _kinematics_round(T, W(v) - (W(h) * W(qm)) * W(E))
    isfinite(primed) || throw(
        ArgumentError(
            "electrostatic priming: the half-step backed-up velocity exceeds the range of $T; rescale the units or reduce dt",
        ),
    )
    return primed
end

# ---------------------------------------------------------------- relativistic priming

# Guard for the momentum back-up u = γ0 v0 - h*qm*(E + v0×B) and v = u/γ(u).
# Underflowing products inside E + v0×B are amplified by |qm| and then |h|; the
# total stays below floatmin*eps*(1 + |h|)*(1 + |qm|), which the scale
# threshold keeps under eps*S/256 for S = max(|γ0 v0|, |h a|).
@inline function _rel_prime_in_range(
    out::NTuple{3,T},
    g,
    a,
    k,
    gm,
    h::T,
    qm::T,
    c2::T,
    gap_resolved::Bool,
) where {T}
    T <: _RangeGuardedFloat || return true
    tiny = floatmin(T)
    S = max(abs(g[1]), abs(g[2]), abs(g[3]), abs(k[1]), abs(k[2]), abs(k[3]))
    βout = map(x -> x / sqrt(c2), out)
    publication_resolved = _lorentz_gap_in_range(sum(abs2, βout))
    return gap_resolved &
           publication_resolved &
           isfinite(c2) &
           (c2 >= tiny) &
           isfinite(h) &
           ((abs(h) >= tiny) | iszero(h)) &
           isfinite(a[1]) &
           isfinite(a[2]) &
           isfinite(a[3]) &
           isfinite(gm) &
           isfinite(out[1]) &
           isfinite(out[2]) &
           isfinite(out[3]) &
           (iszero(S) | (S >= tiny * T(256) * (one(T) + abs(h)) * (one(T) + abs(qm))))
end

# Wide evaluation of the same back-up. (|u|/c)² can exceed the range of the
# next IEEE format for extreme c, so every flagged type uses 256-bit BigFloat
# with the exact Lorentz gap. Planar callers pass exact zeros for the absent
# components. The result |v| < c always fits T.
@noinline function _rel_prime_repair(out::NTuple{3,T}, v0, e, b, qm::T, h::T, c::T) where {T}
    all(isfinite, (v0..., e..., b..., qm, h, c)) || return out
    if T <: _RangeGuardedFloat
        candidate, resolved = _rel_prime_dd_attempt(out, v0, e, b, qm, h, c)
        resolved && return candidate
    end
    gap = _exact_lorentz_gap(v0..., c)
    gap > 0 || throw(ArgumentError("leapfrog priming: input velocity must satisfy |v| < c"))
    wide = _kinematics_bigfloat() do
        W = BigFloat
        g0 = inv(sqrt(W(gap)))
        vx, vy, vz = map(W, v0)
        ex, ey, ez = map(W, e)
        bx, by, bz = map(W, b)
        Qm, H = W(qm), W(h)
        ax = Qm * (ex + (vy * bz - vz * by))
        ay = Qm * (ey + (vz * bx - vx * bz))
        az = Qm * (ez + (vx * by - vy * bx))
        umx = g0 * vx - H * ax
        umy = g0 * vy - H * ay
        umz = g0 * vz - H * az
        gm = sqrt(one(W) + (umx * umx + umy * umy + umz * umz) / W(c)^2)
        (umx / gm, umy / gm, umz / gm)
    end
    return map(x -> _kinematics_round(T, x), wide)
end

# ---------------------------------------------------------------- batch certificates

# Batch owners run one fused reduction over their inputs (no calls, no writes)
# and admit the whole batch to the native loop when a global bound proves that
# every particle passes its per-particle guard; that loop is the historical code
# and vectorizes. Otherwise the guarded loop runs on the untouched inputs and
# returns the same native value for every particle its guard admits.
@inline _batch_sample(a::AbstractVector, p, ::Type{T}) where {T} = T(@inbounds a[p])
@inline _batch_sample(a::Number, p, ::Type{T}) where {T} = T(a)
@inline _batch_nonzero_min(m, s) = min(m, ifelse(iszero(s), oftype(s, Inf), s))
@inline _abs3(a, b, c) = max(abs(a), abs(b), abs(c))

# Kernels selected by a certificate: Val(false) is the native expression.
@inline _kick_eval(::Val{false}, args...) = _boris_kick_native(args...)
@inline _kick_eval(::Val{true}, args...) = boris_kick(args...)
@inline _rel_kick_eval(::Val{false}, args...) = _boris_rel_native(args...)
@inline _rel_kick_eval(::Val{true}, args...) = _boris_rel(args...)
@inline _coordinate_eval(::Val{false}, x, v, dt) = _particle_coordinate_native(x, v, dt)
@inline _coordinate_eval(::Val{true}, x, v, dt) = _particle_coordinate(x, v, dt)

# max|v|, max|E|, max|B| and the smallest nonzero Boris guard scale
# S_p = max(|v_p|, |h E_p|) (same h*E products). Fields are per-particle vectors
# or uniform numbers; samples are converted to T as the kick loops do.
@inline function _boris_batch_scales(v, E, B, h::T, indices) where {T}
    vmax = emax = bmax = zero(T)
    smin = T(Inf)
    @inbounds for p in indices
        ex, ey, ez = _batch_sample(E[1], p, T), _batch_sample(E[2], p, T), _batch_sample(E[3], p, T)
        bx, by, bz = _batch_sample(B[1], p, T), _batch_sample(B[2], p, T), _batch_sample(B[3], p, T)
        av = _abs3(_batch_sample(v[1], p, T), _batch_sample(v[2], p, T), _batch_sample(v[3], p, T))
        vmax = max(vmax, av)
        emax = max(emax, _abs3(ex, ey, ez))
        bmax = max(bmax, _abs3(bx, by, bz))
        smin = _batch_nonzero_min(smin, max(av, _abs3(h * ex, h * ey, h * ez)))
    end
    return vmax, emax, bmax, smin
end

# Uniform-field kick plus drift (push_uniform!): one fused pass over velocities
# and positions. With uniform fields S_p = max(|v_p|, K), K = max|h E_i|.
# With Val(false) the caller has shown K itself passes the scale threshold, so
# every S_p >= K does and only the maxima are reduced.
@inline function _uniform_batch_scales(
    v,
    x::NTuple{D},
    K::T,
    indices,
    ::Val{scale},
) where {D,T,scale}
    vmax = xmax = zero(T)
    smin = scale ? T(Inf) : K
    vx, vy, vz = v
    @inbounds for p in indices
        av = _abs3(vx[p], vy[p], vz[p])
        vmax = max(vmax, av)
        scale && (smin = _batch_nonzero_min(smin, max(av, K)))
        xmax = max(xmax, _abs_at(x, p, T))
    end
    return vmax, smin, xmax
end

# max_d |x[d][p]|, unrolled over the coordinate tuple (keeps the loop vectorizable).
@inline _abs_at(::Tuple{}, p, ::Type{T}) where {T} = zero(T)
@inline _abs_at(x::Tuple, p, ::Type{T}) where {T} =
    max(abs(@inbounds first(x)[p]), _abs_at(Base.tail(x), p, T))

# Every particle passes the Boris guard (boris.jl) when h is admissible and, with
# k = 2|h|max|E| >= |h E_i| and tau = 2|h|max|B| >= |t_i| (the factor 2 covers
# roundings): tau <= sqrt(floatmax)/4 gives |t|² <= floatmax/4, so |t|² is finite
# and f = 2/(1 + |t|²) >= 8/floatmax >= floatmin; |v⁻| <= max|v| + k,
# |v⁻×t| <= 2|v⁻|tau, |s| <= 1 and |v′×s| <= 2|v′| bound every intermediate and
# output by 10(max|v| + k)(1 + 2tau) < bound = 4(max|v| + k)(1 + 2tau)² <=
# floatmax/4; and tau_p <= tau, so S_p >= smin >= 256 floatmin (1 + tau) passes
# the scale test. Returns the certificate and a bound on the kicked |v|.
@inline function _boris_batch_certificate(
    h::T,
    qm::T,
    dt::T,
    vmax::T,
    emax::T,
    bmax::T,
    smin::T,
) where {T}
    T <: _RangeGuardedFloat || return true, zero(T)
    k = 2 * (abs(h) * emax)
    tau = 2 * (abs(h) * bmax)
    bound = 4 * (vmax + k) * (1 + 2tau)^2
    ok =
        _kinematics_product_in_range(h, qm, dt) &
        (tau <= sqrt(floatmax(T)) / 4) &
        (bound <= floatmax(T) / 4) &
        (smin >= floatmin(T) * T(256) * (one(T) + tau))
    return ok, 4bound
end

# Relativistic Boris: the rounded test β² < 1 leaves 1 - β² >= eps/2, so
# γ0 <= G/2 with G = 2 sqrt(2/eps) and |u| <= G max|v|; γ⁻ >= 1 only shrinks t,
# so the nonrelativistic bounds hold with w = G max|v| + k. |u⁻|²/c² and
# |u⁺|²/c² are below 3(2.5 bound)²/c² <= floatmax when bound²/c² <= floatmax/32,
# so γ⁻ and γ⁺ are finite and |v⁺| < c. γ0 >= 1 makes |u_p| >= |v_p|, so the
# velocity-based smin bounds the momentum scale from below. The additional
# vmax <= c/4 condition gives β² <= 3/16, so no native batch bypasses the
# Lorentz-gap conditioning guard.
@inline function _boris_rel_batch_certificate(
    h::T,
    qm::T,
    dt::T,
    c::T,
    vmax::T,
    emax::T,
    bmax::T,
    smin::T,
) where {T}
    T <: _RangeGuardedFloat || return true
    c2 = c * c
    k = 2 * (abs(h) * emax)
    tau = 2 * (abs(h) * bmax)
    bound = 4 * (2 * sqrt(2 / eps(T)) * vmax + k) * (1 + 2tau)^2
    return (vmax <= c / 4) &
           isfinite(c2) &
           (c2 >= floatmin(T)) &
           _kinematics_product_in_range(h, qm, dt) &
           (tau <= sqrt(floatmax(T)) / 4) &
           (bound <= floatmax(T) / 4) &
           ((bound * bound) / c2 <= floatmax(T) / 32) &
           (smin >= floatmin(T) * T(256) * (one(T) + tau))
end

@inline _batch_kick_magnitude(e, factor::Number) = abs(factor * e)
@inline _batch_kick_magnitude(e, factor::Tuple) = abs(factor[1] * (factor[2] * e))

# Euler half-kick scales: max|u|, max|E|, max|B| and the smallest nonzero lower
# bound σ_p of the guard scale S_p = max(|u_p|, |k_p|): |u_p| for a moving
# particle and, for u_p = 0, exactly |k_p| = |factor E_p| because then
# E + u×B = E ± 0. `factor` is h*qm, or (h, qm) for the relativistic h*(qm*E).
@inline function _euler_batch_scales(u, E, B, factor, indices, ::Type{T}) where {T}
    umax = emax = bmax = zero(T)
    smin = T(Inf)
    @inbounds for p in indices
        ex, ey, ez = _batch_sample(E[1], p, T), _batch_sample(E[2], p, T), _batch_sample(E[3], p, T)
        bx, by, bz = _batch_sample(B[1], p, T), _batch_sample(B[2], p, T), _batch_sample(B[3], p, T)
        au = _abs3(_batch_sample(u[1], p, T), _batch_sample(u[2], p, T), _batch_sample(u[3], p, T))
        ak = max(
            _batch_kick_magnitude(ex, factor),
            _batch_kick_magnitude(ey, factor),
            _batch_kick_magnitude(ez, factor),
        )
        umax = max(umax, au)
        emax = max(emax, _abs3(ex, ey, ez))
        bmax = max(bmax, _abs3(bx, by, bz))
        smin = _batch_nonzero_min(smin, ifelse(iszero(au), ak, au))
    end
    return umax, emax, bmax, smin
end

# |E + u×B| <= max|E| + 2 max|u| max|B| and every intermediate and result of
# u ∓ h qm (E + u×B) is below bound/4 <= floatmax/16 (the factor 4 covers the
# roundings); S_p >= σ_p >= smin passes the threshold of _euler_kick_in_range.
@inline function _euler_batch_certificate(h::T, qm::T, hq::T, umax, emax, bmax, smin) where {T}
    T <: _RangeGuardedFloat || return true
    bound = 4 * (umax + abs(hq) * (emax + 2 * (umax * bmax)))
    return _kinematics_product_in_range(hq, h, qm) &
           (bound <= floatmax(T) / 4) &
           (smin >= floatmin(T) * T(256) * (one(T) + abs(hq)))
end

# Relativistic back-up (guard _rel_prime_in_range): |a| <= A = 2|qm|(max|E| +
# 2 max|u| max|B|), |γ0 v0| <= G max|u| with G as above, |um| <= w/4 and
# |um|²/c² <= w²/c² <= floatmax/32, so a, γ and the result (|v| < c) are finite;
# γ0 >= 1 makes σ_p a lower bound of max(|γ0 v0|, |h a|). Requiring
# umax <= c/4 also proves the Lorentz gap is well conditioned for every row.
@inline function _rel_prime_batch_certificate(h::T, qm::T, c::T, umax, emax, bmax, smin) where {T}
    T <: _RangeGuardedFloat || return true
    c2 = c * c
    A = 2 * (abs(qm) * (emax + 2 * (umax * bmax)))
    w = 4 * (2 * sqrt(2 / eps(T)) * umax + abs(h) * A)
    return (umax <= c / 4) &
           isfinite(c2) &
           (c2 >= floatmin(T)) &
           isfinite(h) &
           ((abs(h) >= floatmin(T)) | iszero(h)) &
           (w <= min(c, floatmax(T) / 4)) &
           ((w * w) / c2 <= floatmax(T) / 32) &
           (smin >= floatmin(T) * T(256) * (one(T) + abs(h)) * (one(T) + abs(qm)))
end

# Electrostatic priming v - (h*qm)*E: the guard's exact per-component scale.
@inline function _es_prime_batch_certificate(v, E, h::T, qm::T, hq::T, C, indices) where {T}
    T <: _RangeGuardedFloat || return true
    vmax = emax = zero(T)
    smin = T(Inf)
    @inbounds for c = 1:C, p in indices
        a, e = v[c][p], E[c][p]
        vmax = max(vmax, abs(a))
        emax = max(emax, abs(e))
        smin = _batch_nonzero_min(smin, max(abs(a), abs(hq * e)))
    end
    return _kinematics_product_in_range(hq, h, qm) &
           (4 * (vmax + abs(hq) * emax) <= floatmax(T) / 4) &
           (smin >= floatmin(T) * T(256))
end

# max|a| over a tuple of arrays, each over its own indices (drift certificates).
@inline function _coordinate_batch_max(arrays::Tuple, ::Type{T}) where {T}
    amax = zero(T)
    for a in arrays
        @inbounds for p in eachindex(a)
            amax = max(amax, abs(a[p]))
        end
    end
    return amax
end

# x + dt*v cannot overflow when max|x| + |dt| max|v| <= floatmax/2 (the product
# and sum roundings add at most 2eps relative), so the native coordinate is the
# guarded one. Midpoints use dt/2 and are covered.
@inline _coordinate_batch_certified(xmax::T, dt::T, vbound::T) where {T} =
    !(T <: _RangeGuardedFloat) || (xmax + abs(dt) * vbound <= floatmax(T) / 2)

# Shock steppers gather fields inside their push loops with interpolation
# weights that are nonnegative and sum to one, so every sample is at most
# 2 max|grid| in magnitude. A particle at rest has the unknown scale |h E_p|, so
# it is certified only when the grid electric field vanishes identically.
function _shock_batch_certified(ps::ParticleSet{D,T}, Egrid::Tuple, Bgrid::Tuple, dt::T) where {D,T}
    T <: _RangeGuardedFloat || return true
    emax = 2 * _coordinate_batch_max(Egrid, T)
    bmax = 2 * _coordinate_batch_max(Bgrid, T)
    vmax = zero(T)
    vmin = vmin_nonzero = T(Inf)
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        av = _abs3(vx[p], vy[p], vz[p])
        vmax = max(vmax, av)
        vmin = min(vmin, av)
        vmin_nonzero = _batch_nonzero_min(vmin_nonzero, av)
    end
    smin = iszero(emax) ? vmin_nonzero : vmin
    kick, vbound = _boris_batch_certificate(one(T) * dt / 2, one(T), dt, vmax, emax, bmax, smin)
    return kick && _coordinate_batch_certified(_coordinate_batch_max(ps.x, T), dt, vbound)
end
