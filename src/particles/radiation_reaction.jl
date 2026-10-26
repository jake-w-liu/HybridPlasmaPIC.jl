# radiation_reaction.jl — Landau-Lifshitz radiation-reaction drag (Phase-4 forcing).
#
# A relativistic charged particle radiating in strong fields feels a back-reaction. This
# applies the dominant (γ²) term of the reduced Landau-Lifshitz force as a post-push
# velocity kick (uniform fields), analogous to how the collision operators act on velocity.
#
# Arithmetic for Float16/32/64 particles. The frozen-rate substep p ← p·exp(-α), with
# α = K·dt·γ·χ², is applied as a velocity factor
#
#     v ← F·v,   F = 1/sqrt(1 + r·(exp(2α) - 1)),   r = 1 - |v|²/c² = 1/γ²,
#
# which is algebraically identical to damping p = γv and recovering v = p/sqrt(1 + |p|²/c²),
# but every computed F lies in [0, 1], so the update never increases a stored speed (RR-003).
# Every particle is first evaluated natively (Float16 particles in Float32) with an
# a-posteriori bound on its rounding error; a particle whose bound exceeds the budget, or
# whose intermediates leave the floating-point range, is re-evaluated in Float64 (Float16 and
# Float32 particles), then in exponent-scaled double-double arithmetic (RR-002). An
# amplified invariant uncertainty or unresolved Lorentz gap is completed from the original
# stored operands with task-local exact integers. Hardware particles use no BigFloat;
# native/double-double paths allocate nothing.

# Absolute budget for the propagated rounding error of a native factor F, in eps(T). A native
# result is within (budget + 4)·eps(T)·|v| a priori (the rest covers evaluating F, v·F and a
# final narrowing); the bounds over-estimate the observed error several-fold (tests assert
# 8 ulps of |v| and observe below 4). The budget trades accuracy against the share of
# particles that take a wider evaluation.
const _RADIATION_BUDGET_EPS = 8

# ----------------------------------------------------------------------------------------
# Error-free transformations (no overflow; exact by Knuth/Dekker)
# ----------------------------------------------------------------------------------------

@inline function _radiation_two_sum(a::W, b::W) where {W<:AbstractFloat}
    s = a + b
    t = s - a
    return s, (a - (s - t)) + (b - t)
end

# Requires |a| ≥ |b| or a == 0.
@inline function _radiation_fast_two_sum(a::W, b::W) where {W<:AbstractFloat}
    s = a + b
    return s, b - (s - a)
end

@inline function _radiation_two_product(a::W, b::W) where {W<:AbstractFloat}
    p = a * b
    return p, fma(a, b, -p)
end

# ----------------------------------------------------------------------------------------
# Native evaluation with an a-posteriori rounding bound
# ----------------------------------------------------------------------------------------

# Loop invariants of the native evaluation in format W for particles stored in format T.
struct _RadiationNative{W<:AbstractFloat}
    E::NTuple{3,W}
    B::NTuple{3,W}
    Kdt::W
    c2::W          # fl(c²); c² = c2 + c2lo exactly
    c2lo::W
    half_c2::W
    invc2::W
    invc::W
    bound0::W      # u·(12.1|E|² + 8 floatmin): see _radiation_native_factor
    bound1::W      # u·18.2|B|²
    budget::W      # budget·eps(T)/1.09 (first-order slack, see below)
    usable::Bool   # loop-invariant range conditions of the error analysis hold
end

function _radiation_native(::Type{W}, ::Type{T}, E, B, K, dt, c) where {W,T}
    Ew = map(W, E)
    Bw = map(W, B)
    cw = W(c)
    u = eps(W) / 2
    Kdt = W(K) * W(dt)
    c2 = cw * cw
    c2lo = fma(cw, cw, -c2)
    E2 = Ew[1] * Ew[1] + Ew[2] * Ew[2] + Ew[3] * Ew[3]
    B2 = Bw[1] * Bw[1] + Bw[2] * Bw[2] + Bw[3] * Bw[3]
    bound0 = u * (W(12.1) * E2 + 8 * floatmin(W))
    bound1 = u * (W(18.2) * B2)
    budget = W(_RADIATION_BUDGET_EPS) * W(eps(T)) / W(1.09)
    # Kdt normal keeps its relative rounding u. c² ≥ 4·floatmin/eps keeps c2lo exact and the
    # squares' underflow (absolute u·floatmin each) below u²·c² in r. Finite bound terms
    # exclude |E|² or |B|² overflow; per-particle overflow is caught by the guard itself.
    usable =
        floatmin(W) <= Kdt <= floatmax(W) &&
        4 * floatmin(W) / eps(W) <= c2 <= floatmax(W) / 4 &&
        isfinite(bound0) &&
        isfinite(bound1)
    return _RadiationNative{W}(
        Ew,
        Bw,
        Kdt,
        c2,
        c2lo,
        c2 / 2,
        inv(c2),
        inv(cw),
        bound0,
        bound1,
        budget,
        usable,
    )
end

# Velocity factor F and whether its error bound passes. Error analysis (u = eps(W)/2, first
# order, standard model with an absolute u·floatmin per underflowing product):
# - r: for β² ≤ 1/2 the rounded c² - S has relative error ≤ 4u; otherwise c² - S is exact
#   (Sterbenz) and the exact residuals of S, the squares and c² are added, leaving
#   ≤ 4u + 16u²/r. With r ≥ 64u and the reciprocal of c², |δr| ≤ 7u in both branches.
# - α = (K·dt/sqrt(r))·χ²: relative error ≤ u + (|δr|/2 + u) + 2u = 7.5u besides χ².
# - χ² = |w|² - P with w = E + v×B and P = (v·E·(1/c))², squared at the force scale so an
#   underflow is not amplified by 1/c²: |Δw_i| ≤ 3.01u·M_i, M_i = |E_i| + |v_j B_k| +
#   |v_k B_j|, Σ M_i² ≤ 2|E|² + 6|v|²|B|² (Cauchy-Schwarz), and P has relative error 5u
#   plus 6u|v·E|Σ|v_i E_i|/c² ≤ 6u|E|². Underflowing products add u·floatmin each, which
#   after AM-GM costs one more u·|w|² and 2u·P. Hence
#   |Δχ²| ≤ u(|χ²| + 7.01|w|² + 7P + 12.04|E|² + 18.06|v|²|B|² + 8 floatmin).
#   With χ²⁺, |χ²| ≤ |w|² + P this gives errα ≤ G·u·(16|w|² + 8P + 12.1|E|² + 18.2|v|²|B|²
#   + 8 floatmin), G = K·dt/sqrt(r). The tight variant replaces the Cauchy-Schwarz scale by
#   the computed products (used only before the double-double evaluation).
# - |dF/dα| = (X + r)/(1 + X)^(3/2) with X = r(exp(2α) - 1); for errα ≤ 1/64 it varies by
#   at most exp(5/64) < 1.09 across the error interval, so |ΔF| ≤ budget·eps(T) follows from
#   1.09·(X + r)·errα ≤ budget·eps(T)·(1 + X)^(3/2). A clamped negative χ² (exact χ² ≥ 0) is
#   covered by the same bound.
# - Evaluating F: exp(2α) - 1 carries absolute error ≤ 2.2u·exp(2α), so the absolute error of
#   F from X stays ≤ 2.7u for every α (normwise; expm1's relative accuracy is not needed),
#   plus 2.5u for 1 + X, sqrt and the reciprocal, and u per component of v·F.
# Non-finite intermediates make errα, X or the comparison NaN/Inf, which fails the guard,
# except X = Inf with a finite errα: then the exact F < 1/sqrt(floatmax(W)) and F = 0 is
# within the budget.
@inline function _radiation_native_factor(
    v::NTuple{3,W},
    p::_RadiationNative{W},
    ::Val{tight},
) where {W<:AbstractFloat,tight}
    u = eps(W) / 2
    vx, vy, vz = v
    Ex, Ey, Ez = p.E
    Bx, By, Bz = p.B
    px, py, pz = vx * vx, vy * vy, vz * vz
    s1 = px + py
    S = s1 + pz
    N = p.c2 - S
    if S > p.half_c2
        t1 = (px - (s1 - (s1 - px))) + (py - (s1 - px))
        t2 = (s1 - (S - (S - s1))) + (pz - (S - s1))
        low = (t1 + t2) + ((fma(vx, vx, -px) + fma(vy, vy, -py)) + fma(vz, vz, -pz))
        N += p.c2lo - low
    end
    r = N * p.invc2
    # A true r below the residual's rounding bound (16u²) could round negative; the guard
    # rejects r < 64u, and the clamp only keeps sqrt and 1 + X in their domains meanwhile.
    rc = max(r, zero(W))
    a, b = vy * Bz, vz * By
    d, e = vz * Bx, vx * Bz
    f, g = vx * By, vy * Bx
    dx, dy, dz = a - b, d - e, f - g
    wx, wy, wz = Ex + dx, Ey + dy, Ez + dz
    qx, qy, qz = vx * Ex, vy * Ey, vz * Ez
    vE = (qx + qy) + qz
    w2 = wx * wx + wy * wy + wz * wz
    βE = vE * p.invc   # scale before squaring: (v·E)² itself may underflow for small c
    P = βE * βE
    χ2 = w2 - P
    G = p.Kdt / sqrt(rc)
    α = G * max(χ2, zero(W))
    X = rc * (exp(2α) - one(W))
    onepX = one(W) + X
    s = sqrt(onepX)
    F = inv(s)
    errα = if tight
        mx = abs(wx) + abs(dx) + abs(a) + abs(b)
        my = abs(wy) + abs(dy) + abs(d) + abs(e)
        mz = abs(wz) + abs(dz) + abs(f) + abs(g)
        Vc = (abs(qx) + abs(qy) + abs(qz)) * p.invc
        first =
            abs(χ2) +
            4 * w2 +
            W(2.01) * (abs(wx) * mx + abs(wy) * my + abs(wz) * mz) +
            7 * P +
            W(6.01) * (abs(βE) * Vc) +
            8 * floatmin(W)
        second = W(1.01) * (mx * mx + my * my + mz * mz) + W(9.1) * (Vc * Vc)
        W(7.5) * u * α + G * (u * first + u * u * second)
    else
        G * fma(16u, w2, fma(8u, P, fma(p.bound1, S, p.bound0)))
    end
    lhs = (X + r) * errα
    rhs = p.budget * (onepX * s)
    ok = (r >= 64u) & (errα <= W(1) / 64) & (lhs <= rhs)
    return F, ok
end

# ----------------------------------------------------------------------------------------
# Exponent-scaled double-double evaluation (Float64)
# ----------------------------------------------------------------------------------------

const _RadiationDD = NTuple{2,Float64}

@inline function _radiation_dd_add(a::_RadiationDD, b::Float64)
    s, e = _radiation_two_sum(a[1], b)
    return _radiation_fast_two_sum(s, e + a[2])
end

# Accurate double-word addition (relative error ≤ 3u² of the exact sum).
@inline function _radiation_dd_add(a::_RadiationDD, b::_RadiationDD)
    s, e = _radiation_two_sum(a[1], b[1])
    t, f = _radiation_two_sum(a[2], b[2])
    s, e = _radiation_fast_two_sum(s, e + t)
    return _radiation_fast_two_sum(s, e + f)
end

@inline _radiation_dd_sub(a::_RadiationDD, b::_RadiationDD) = _radiation_dd_add(a, (-b[1], -b[2]))

@inline function _radiation_dd_mul(a::_RadiationDD, b::_RadiationDD)
    p, e = _radiation_two_product(a[1], b[1])
    return _radiation_fast_two_sum(p, fma(a[1], b[2], fma(a[2], b[1], e)))
end

@inline function _radiation_dd_mul(a::_RadiationDD, b::Float64)
    p, e = _radiation_two_product(a[1], b)
    return _radiation_fast_two_sum(p, fma(a[2], b, e))
end

@inline function _radiation_dd_div(a::_RadiationDD, b::_RadiationDD)
    q1 = a[1] / b[1]
    rem = _radiation_dd_sub(a, _radiation_dd_mul(b, q1))
    q2 = rem[1] / b[1]
    rem = _radiation_dd_sub(rem, _radiation_dd_mul(b, q2))
    q3 = rem[1] / b[1]
    return _radiation_dd_add(_radiation_fast_two_sum(q1, q2), q3)
end

@inline function _radiation_dd_sqrt(a::_RadiationDD)   # a[1] > 0
    x = sqrt(a[1])
    p, e = _radiation_two_product(x, x)
    return _radiation_fast_two_sum(x, (((a[1] - p) - e) + a[2]) / (2x))
end

# Exact sum of Float64 values as a nonoverlapping expansion (Shewchuk's grow-expansion),
# then rounded to a double-double by accumulating from the smallest component.
@inline _radiation_grow(::Tuple{}, b::Float64) = (b,)
@inline function _radiation_grow(e::Tuple, b::Float64)
    q, h = _radiation_two_sum(b, first(e))
    return (h, _radiation_grow(Base.tail(e), q)...)
end
@inline _radiation_expansion(t::Tuple{Float64}) = t
@inline _radiation_expansion(t::Tuple) =
    _radiation_grow(_radiation_expansion(Base.front(t)), last(t))
@inline function _radiation_exact_sum(t::Tuple)
    acc = (0.0, 0.0)
    for x in _radiation_expansion(t)
        acc = _radiation_dd_add(acc, x)
    end
    return acc
end

struct _RadiationScaled
    E::NTuple{3,Float64}
    B::NTuple{3,Float64}
    hasE::Bool
    hasB::Bool
    eE::Int              # binary exponents of max|E_i| and max|B_i|
    eB::Int
    kc::Int              # c = ĉ·2^kc, ĉ ∈ [1, 2)
    c2::_RadiationDD     # ĉ² exactly
    invc2::_RadiationDD
    Kdt::_RadiationDD    # K·dt = Kdt·2^eKdt exactly (frexp mantissas)
    eKdt::Int
    # Concrete fields only: no implicit converting constructor.
    function _RadiationScaled(
        E::NTuple{3,Float64},
        B::NTuple{3,Float64},
        hasE::Bool,
        hasB::Bool,
        eE::Int,
        eB::Int,
        kc::Int,
        c2::_RadiationDD,
        invc2::_RadiationDD,
        Kdt::_RadiationDD,
        eKdt::Int,
    )
        return new(E, B, hasE, hasB, eE, eB, kc, c2, invc2, Kdt, eKdt)
    end
end

function _radiation_scaled(E, B, K, dt, c)
    Ef = map(Float64, E)
    Bf = map(Float64, B)
    Emax = max(abs(Ef[1]), abs(Ef[2]), abs(Ef[3]))
    Bmax = max(abs(Bf[1]), abs(Bf[2]), abs(Bf[3]))
    cf = Float64(c)
    kc = exponent(cf)
    ĉ = ldexp(cf, -kc)
    c2 = _radiation_two_product(ĉ, ĉ)
    mK, eK = frexp(Float64(K))
    md, ed = frexp(Float64(dt))
    return _RadiationScaled(
        Ef,
        Bf,
        !iszero(Emax),
        !iszero(Bmax),
        iszero(Emax) ? 0 : exponent(Emax),
        iszero(Bmax) ? 0 : exponent(Bmax),
        kc,
        c2,
        _radiation_dd_div((1.0, 0.0), c2),
        _radiation_two_product(mK, md),
        eK + ed,
    )
end

@inline function _radiation_scale_exact(value::Rational{BigInt}, power::Int)
    return power >= 0 ? value * (BigInt(1) << power) : value / (BigInt(1) << -power)
end

# The high word is the nearest Float64 to log(2); the low word is its rounded
# residual. Keeping both prevents a large binary exponent from magnifying the
# high word's rounding error near the cooling transition.
const _RADIATION_LOG2 = (0x1.62e42fefa39efp-1, 0x1.abc9e3b39803fp-56)

# Complete the invariant from the original stored operands when its field-scale
# double-double error could be amplified into a material velocity error. Integer
# arithmetic and the shared exact converters use no ambient MPFR state.
@noinline function _radiation_exact_factor(v::NTuple{3,Float64}, s::_RadiationScaled)
    velocity = map(_exact_rational, v)
    electric = map(_exact_rational, s.E)
    magnetic = map(_exact_rational, s.B)
    c2 = _radiation_scale_exact(_exact_rational(s.c2[1]) + _exact_rational(s.c2[2]), 2s.kc)
    r = (c2 - sum(abs2, velocity)) / c2
    force = ntuple(3) do i
        j, k = i % 3 + 1, (i + 1) % 3 + 1
        electric[i] + velocity[j] * magnetic[k] - velocity[k] * magnetic[j]
    end
    chi2 = sum(abs2, force) - sum(velocity .* electric)^2 / c2
    iszero(chi2) && return 1.0
    kdt = _radiation_scale_exact(_exact_rational(s.Kdt[1]) + _exact_rational(s.Kdt[2]), s.eKdt)
    alpha2 = (kdt * chi2)^2 / r
    # For finite binary64 operands, r > 2^-4196. At alpha >= 4096,
    # floatmax(Float64)*gamma*exp(-alpha) is below the minimum subnormal.
    alpha2 >= 4096^2 && return 0.0
    # |dF/dalpha| <= 1: this endpoint is below a Float64 normwise ulp.
    alpha2 < _exact_rational(0x1p-120) && return 1.0
    alpha = _sqrt_from_exact(Float64, alpha2)
    residual = alpha2 - _exact_rational(alpha)^2
    low = _from_exact(Float64, residual / (2 * _exact_rational(alpha)))
    rh = _from_exact(Float64, r)
    rl = _from_exact(Float64, r - _exact_rational(rh))
    if 2alpha <= 708
        em1 = expm1(2alpha)
        exponential = em1 + (1 + em1) * (2low)
        x = fma(rh, exponential, rl * exponential)
        return inv(sqrt(1 + x))
    end
    # Normalize the exact ratio before taking its logarithm, including ratios
    # too small for a Float64. Only the normalized mantissa enters libm.
    n, d = numerator(r), denominator(r)
    power = ndigits(n; base = 2) - ndigits(d; base = 2)
    (power >= 0 ? n < d << power : n << -power < d) && (power -= 1)
    mantissa = _radiation_scale_exact(r, -power)
    mh = _from_exact(Float64, mantissa)
    ml = _from_exact(Float64, mantissa - _exact_rational(mh))
    lnr = _radiation_dd_add(_radiation_dd_mul(_RADIATION_LOG2, Float64(power)), log(mh))
    lnr = _radiation_dd_add(lnr, ml / mh)
    argument = _radiation_dd_add(lnr, (2alpha, 2low))
    argument[1] < -60 * _RADIATION_LOG2[1] && return 1.0
    if argument[1] <= 50 * _RADIATION_LOG2[1]
        exponential = exp(argument[1]) * (1 + argument[2])
        return inv(sqrt(1 + (exponential - rh)))
    end
    # Here F < 2^-25. Omitting log1p((1-r)*exp(-argument)) changes
    # F by less than 2^-75, below the normwise velocity budget.
    return exp(-(argument[1] + argument[2]) / 2)
end

# Velocity factor for finite subluminal v and positive K, dt. Powers of two keep
# the ordinary invariant operations bounded. The DD certificate includes the
# coefficient/timestep and Lorentz amplification before discarding a small force;
# unresolved cases reconstruct the invariant exactly. The gap denominator has a
# separate sensitivity bound. The normalized-DD/IEEE assumptions, operation
# inventory and libm limits are recorded in docs/RADIATION_FORCE_SPAN_PLAN.md.
function _radiation_scaled_factor(v::NTuple{3,Float64}, s::_RadiationScaled)
    vmax = max(abs(v[1]), abs(v[2]), abs(v[3]))
    (iszero(vmax) || !(s.hasE || s.hasB)) && return 1.0
    v̂ = map(x -> ldexp(x, -s.kc), v)
    p1 = _radiation_two_product(v̂[1], v̂[1])
    p2 = _radiation_two_product(v̂[2], v̂[2])
    p3 = _radiation_two_product(v̂[3], v̂[3])
    N = _radiation_dd_sub(_radiation_dd_sub(_radiation_dd_sub(s.c2, p1), p2), p3)
    if !(N[1] > ldexp(s.c2[1], -40))
        # Partial sums can exceed N by 2^40; recover N exactly (the preflight guarantees
        # |v| < c exactly) up to the underflow noted below.
        N = _radiation_exact_sum((s.c2[2], s.c2[1], -p1[2], -p1[1], -p2[2], -p2[1], -p3[2], -p3[1]))
    end
    r = _radiation_dd_mul(N, s.invc2)
    # Below this range the discarded low square parts need not be small
    # relative to the Lorentz gap. Complete it from the stored operands.
    r[1] >= 0x1p-800 || return _radiation_exact_factor(v, s)
    ev = exponent(vmax)
    kF = if !s.hasB
        s.eE + 1
    elseif !s.hasE
        ev + s.eB + 2
    else
        max(s.eE + 1, ev + s.eB + 2)
    end
    # With |vtilde_i| < 2, |Etilde_i| < 1, |Btilde_i| < 1/2 and
    # 1 <= V2 < 12, the bounded double-double operations below have a
    # conservative 2^-80 invariant error in units of 2^(2kF). This allows
    # over 2^24 elementary error coefficients at 106-bit precision. Include
    # gamma and K*dt before accepting that error: an exponent <= 24 keeps
    # its optical-depth contribution below 2^-56. Large amplification must
    # retain even a small force component whose scaled square underflows.
    gain_exponent = s.eKdt + 2kF + cld(-exponent(r[1]), 2)
    gain_exponent <= 24 || return _radiation_exact_factor(v, s)
    # |ṽ_i| < 2, |Ẽ_i| < 1 and |B̃_i| < 1/2, so every scaled product is below 1.
    ṽ = map(x -> ldexp(x, -ev), v)
    Ẽ = map(x -> ldexp(x, -kF), s.E)
    B̃ = map(x -> ldexp(x, ev - kF), s.B)
    w = ntuple(3) do i
        j, k = i % 3 + 1, (i + 1) % 3 + 1
        a1, a2 = _radiation_two_product(ṽ[j], B̃[k])
        b1, b2 = _radiation_two_product(ṽ[k], B̃[j])
        acc = _radiation_dd_add(_radiation_dd_add((Ẽ[i], 0.0), a1), -b1)
        _radiation_dd_add(_radiation_dd_add(acc, a2), -b2)
    end
    q = ntuple(i -> _radiation_two_product(ṽ[i], Ẽ[i]), 3)
    D = _radiation_dd_add(_radiation_dd_add(q[1], q[2]), q[3])          # ṽ·Ẽ
    m = ntuple(i -> _radiation_two_product(ṽ[i], ṽ[i]), 3)
    V2 = _radiation_dd_add(_radiation_dd_add(m[1], m[2]), m[3])         # |ṽ|² ∈ [1, 12)
    # χ̃²·V2² = |V2·w̃ - D·ṽ|² + r·D²·V2: V2·w̃ - D·ṽ is V2 times w⊥ (w normal to v).
    n = ntuple(i -> _radiation_dd_sub(_radiation_dd_mul(V2, w[i]), _radiation_dd_mul(D, ṽ[i])), 3)
    χ2V4 = _radiation_dd_add(
        _radiation_dd_add(
            _radiation_dd_add(_radiation_dd_mul(n[1], n[1]), _radiation_dd_mul(n[2], n[2])),
            _radiation_dd_mul(n[3], n[3]),
        ),
        _radiation_dd_mul(r, _radiation_dd_mul(_radiation_dd_mul(D, D), V2)),
    )
    χ2V4[1] > 0 || return 1.0
    A = _radiation_dd_div(
        _radiation_dd_mul(s.Kdt, χ2V4),
        _radiation_dd_mul(_radiation_dd_mul(V2, V2), _radiation_dd_sqrt(r)),
    )
    A[1] > 0 || return 1.0  # underflow is below the certified absolute budget
    e = s.eKdt + 2kF
    ea = exponent(A[1]) + e
    ea >= 12 && return 0.0
    ea < -60 && return 1.0
    αh, αl = ldexp(A[1], e), ldexp(A[2], e)
    if 2αh <= 708
        E1 = expm1(2αh)
        Em = E1 + (1 + E1) * (2αl)
        X = fma(r[1], Em, r[2] * Em)
        X <= 0x1p50 && return inv(sqrt(1 + X))
    end
    # 1 + X = r·exp(2α)·(1 + q), q = exp(-2α)(1 - r)/r ≤ 2^-50 here.
    lnr = log(r[1]) + r[2] / r[1]
    qq = exp(log1p(-r[1]) - lnr - 2αh)
    return exp(-(lnr + 2αh + 2αl + log1p(qq)) / 2)
end

# ----------------------------------------------------------------------------------------
# Per-particle evaluation chain
# ----------------------------------------------------------------------------------------

# Float16 particles evaluate natively in Float32, then in Float64; Float32 particles natively,
# then in Float64; Float64 particles natively, then with the tight bound. All finish in the
# scaled double-double or exact-invariant evaluation. Inputs convert exactly into each
# wider format.
_radiation_native_format(::Type{Float16}) = Float32
_radiation_native_format(::Type{T}) where {T<:Union{Float32,Float64}} = T

struct _RadiationPlan{T,N,M}
    native::N                # _RadiationNative{_radiation_native_format(T)}
    wide::M                  # _RadiationNative{Float64} or nothing (Float64 particles)
    scaled::_RadiationScaled
end

function _radiation_plan(::Type{T}, E, B, K, dt, c) where {T<:Union{Float16,Float32,Float64}}
    native = _radiation_native(_radiation_native_format(T), T, E, B, K, dt, c)
    wide = T === Float64 ? nothing : _radiation_native(Float64, T, E, B, K, dt, c)
    return _RadiationPlan{T,typeof(native),typeof(wide)}(
        native,
        wide,
        _radiation_scaled(E, B, K, dt, c),
    )
end

@inline function _radiation_velocity(v::NTuple{3,T}, plan::_RadiationPlan{T}) where {T}
    p = plan.native
    if p.usable
        W = _radiation_native_format(T)
        vw = map(W, v)
        F, ok = _radiation_native_factor(vw, p, Val(false))
        ok && return map(x -> T(x * F), vw)
    end
    return _radiation_velocity_fallback(v, plan)
end

@noinline function _radiation_velocity_fallback(v::NTuple{3,T}, plan::_RadiationPlan{T}) where {T}
    vf = map(Float64, v)
    F, ok = if T === Float64
        p = plan.native
        p.usable ? _radiation_native_factor(vf, p, Val(true)) : (1.0, false)
    else
        p = plan.wide
        p.usable ? _radiation_native_factor(vf, p, Val(false)) : (1.0, false)
    end
    factor = ok ? F : _radiation_scaled_factor(vf, plan.scaled)
    return map(x -> T(x * factor), vf)
end

# Other element types (BigFloat and generic AbstractFloat) keep the original native momentum
# update at the ambient precision and rounding mode.
function _radiation_reaction_generic!(
    ps::ParticleSet{D,T},
    Ex,
    Ey,
    Ez,
    Bx,
    By,
    Bz,
    dtT,
    KT,
    cT,
) where {D,T}
    c2 = cT * cT
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(vx)
        ux = vx[p]
        uy = vy[p]
        uz = vz[p]
        v2 = ux * ux + uy * uy + uz * uz
        γ = one(T) / sqrt(one(T) - v2 / c2)
        wx = Ex + (uy * Bz - uz * By)                      # w = E + v×B
        wy = Ey + (uz * Bx - ux * Bz)
        wz = Ez + (ux * By - uy * Bx)
        vE = ux * Ex + uy * Ey + uz * Ez
        χ2 = (wx * wx + wy * wy + wz * wz) - (vE * vE) / c2 # radiation invariant
        χ2 > 0 || continue
        # dp/dt = F_rad = −K γ² χ² v = −(K γ χ²) p, so p decays by exactly exp(−K γ χ² dt)
        s = γ * exp(-KT * γ * χ2 * dtT)
        px = s * ux
        py = s * uy
        pz = s * uz
        γn = sqrt(one(T) + (px * px + py * py + pz * pz) / c2)
        vx[p] = px / γn
        vy[p] = py / γn
        vz[p] = pz / γn
    end
    return ps
end

function _radiation_reaction_ieee!(ps::ParticleSet{D,T}, E, B, dtT, KT, cT) where {D,T}
    plan = _radiation_plan(T, E, B, KT, dtT, cT)
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(vx)
        vx[p], vy[p], vz[p] = _radiation_velocity((vx[p], vy[p], vz[p]), plan)
    end
    return ps
end

"""
    apply_radiation_reaction!(ps::ParticleSet{D,T}, E, B, dt; K, c=1.0) -> ps

One substep of the **Landau-Lifshitz radiation-reaction** drag (dominant `γ²` term) in
spatially uniform fields `E, B` (3-tuples, per unit `q/m`). For each particle of velocity
`v`, `γ = 1/√(1 − |v|²/c²)`:

    w = E + v×B,   χ² = |w|² − (v·E/c)²,   F_rad = −K γ² χ² v

The relativistic momentum is `p = γ v`; `K` is the radiation-reaction coefficient
(∝ the classical electron radius). A gyrating particle in a magnetic field cools
synchrotron-like, `d(1/γ)/dt = K v⁴ B²/c² ≈ K c² B²` ultrarelativistically (so
`1/γ ≈ 1/γ₀ + K c² B² t`; this reduces to the familiar `K B²` rate only when `c = 1`).

The exponential substep `p ← p·exp(−α)`, `α = K γ χ² dt`, holds the initial drag rate
fixed. It approximates the nonlinear drag equation; it is exact only for that frozen
rate. For `Float16`/`Float32`/`Float64` particles it is applied as `v ← F v` with
`F = 1/√(1 + (1 − |v|²/c²)(e^{2α} − 1)) ∈ [0, 1]`, so no stored speed component
increases. Results are normwise accurate: within a few ulps of `|v|` of the exact
substep for the stored inputs. Particles whose native rounding bound or range fails
are re-evaluated in a wider format or exponent-scaled double-double arithmetic; amplified
invariant uncertainty is completed from the original operands with exact integers. Other
element types such as `BigFloat` use the native momentum update at the ambient
precision. `K ≥ 0`, `dt ≥ 0`, `c > 0`; these parameters and `E, B` must remain finite in
the particle precision. Inputs are validated before a `K=0` or `dt=0` no-op. Returns `ps`.
"""
function apply_radiation_reaction!(
    ps::ParticleSet{D,T},
    E::NTuple{3,<:Real},
    B::NTuple{3,<:Real},
    dt::Real;
    K::Real,
    c::Real = 1.0,
) where {D,T}
    K >= 0 || throw(ArgumentError("radiation-reaction coefficient K must be ≥ 0"))
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    KT = _require_finite_nonnegative_real("radiation-reaction coefficient K", K, T)
    dtT = _validated_nonnegative_dt(T, dt; name = "apply_radiation_reaction!")
    cT = _require_finite_positive_real("c", c, T)
    Ex, Ey, Ez = _validated_uniform_field("E", E, T)
    Bx, By, Bz = _validated_uniform_field("B", B, T)
    (iszero(KT) || iszero(dtT)) && return ps
    _require_subluminal_velocities(ps, cT, "apply_radiation_reaction!")
    if T <: Union{Float16,Float32,Float64}
        return _radiation_reaction_ieee!(ps, (Ex, Ey, Ez), (Bx, By, Bz), dtT, KT, cT)
    end
    return _radiation_reaction_generic!(ps, Ex, Ey, Ez, Bx, By, Bz, dtT, KT, cT)
end
