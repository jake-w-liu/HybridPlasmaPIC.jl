# inject.jl — flux-weighted particle injection at an open boundary (§11.4).
#
# Particles crossing an inflow boundary are NOT sampled from the volume
# Maxwellian: the inward flux distribution carries a factor of normal speed,
#   p(v_n | v_n>0) ∝ v_n · f_M(v_n),
# so a boundary sampler must reproduce the target number/momentum/energy fluxes.
# Here `flux_speed` samples the inward normal speed and `inject_face_1d!`
# appends a flux-correct, accumulator-metered batch each step.

# erf via Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7) — avoids a SpecialFunctions dep
@inline function _erf(x::T) where {T}
    s = sign(x)
    ax = abs(x)
    t = one(T) / (one(T) + T(0.3275911) * ax)
    y =
        one(T) -
        (
            ((((T(1.061405429) * t - T(1.453152027)) * t) + T(1.421413741)) * t - T(0.284496736)) *
            t + T(0.254829592)
        ) *
        t *
        exp(-ax^2)
    return s * y
end

# unnormalized ∫₀ˢ s' exp(−(s'−a)²/2σ²) ds'
@inline function _flux_integral(s::T, a::T, σ::T) where {T}
    u0 = -a / σ
    u1 = (s - a) / σ
    r2 = sqrt(T(2))
    return σ * (
        a * sqrt(T(π) / 2) * (_erf(u1 / r2) - _erf(u0 / r2)) +
        σ * (exp(-u0^2 / 2) - exp(-u1^2 / 2))
    )
end

# Dimensionless form of the same integral, with x=s/σ and d=a/σ. The
# common σ² factor cancels during inverse-CDF sampling; omitting it prevents
# recoverable overflow when both drift and thermal speed are large.
@inline function _flux_shape_integral(x::T, d::T) where {T}
    u0 = -d
    u1 = x - d
    r2 = sqrt(T(2))
    return d * sqrt(T(π) / 2) * (_erf(u1 / r2) - _erf(u0 / r2)) + exp(-(u0 * u0) / 2) -
           exp(-(u1 * u1) / 2)
end

@inline function _validated_flux_sampler_params(a::T, σ::T) where {T}
    aT = _require_finite_real("a", a, T)
    σT = _require_finite_nonnegative_real("σ", σ, T)
    return aT, σT
end

# Scaled complementary error function for x ≥ 0, using the Numerical Recipes
# minimax form. Unlike `1-erf(x)`, this retains relative accuracy in the tail.
@inline function _erfcx_positive(x::T) where {T}
    t = one(T) / (one(T) + T(0.5) * x)
    p =
        T(-1.26551223) +
        t * (
            T(1.00002368) +
            t * (
                T(0.37409196) +
                t * (
                    T(0.09678418) +
                    t * (
                        T(-0.18628806) +
                        t * (
                            T(0.27886807) +
                            t * (
                                T(-1.13520398) +
                                t * (T(1.48851587) + t * (T(-0.82215223) + t * T(0.17087277)))
                            )
                        )
                    )
                )
            )
        )
    return t * exp(p)
end

# H(b) = ∫₀∞ x exp(-b*x-x²/2) dx. The erfcx form is stable for moderate b;
# its subtraction loses digits for large b, where the optimally-truncated
# asymptotic series is rapidly accurate.
function _negative_flux_shape(b::T) where {T}
    if b < T(6)
        R = sqrt(T(π) / T(2)) * _erfcx_positive(b / sqrt(T(2)))
        return max(one(T) - b * R, zero(T))
    end
    invb = inv(b)
    invb2 = invb * invb
    term = invb2
    total = term
    for n = 1:64
        nextterm = term * T(2n + 1) * invb2
        nextterm >= term && break
        candidate = isodd(n) ? total - nextterm : total + nextterm
        total = candidate
        nextterm <= eps(T) * abs(total) && break
        term = nextterm
    end
    return max(total, zero(T))
end

@inline function _negative_flux_per_density(σ::T, b::T) where {T}
    normal_decay_limit = sqrt(-T(2) * log(floatmin(T)))
    if b <= normal_decay_limit
        decay = exp(-(b * b) / 2)
        if decay >= floatmin(T)
            # Retain the ordinary-range operation order and its established
            # accuracy whenever the Gaussian itself is normally representable.
            return σ / sqrt(T(2π)) * decay * _negative_flux_shape(b)
        end
    end

    # In the recoverable tail, exp(-b²/2) may underflow even though its
    # product with a large σ is finite. Split it into two exp(-b²/4)
    # factors and place the finite prefactor between them. Once a half-decay
    # is below the smallest subnormal, even floatmax(T) cannot rescue the full
    # decay, so returning zero is the correctly rounded result.
    half_decay_limit = T(2) * sqrt(-log(nextfloat(zero(T))))
    b <= half_decay_limit || return zero(T)
    half_decay = exp(-(b * b) / 4)
    iszero(half_decay) && return zero(T)
    prefactor = σ / sqrt(T(2π)) * _negative_flux_shape(b)
    return (prefactor * half_decay) * half_decay
end

"""
    flux_per_density(a, σ)

Inward number flux per unit upstream density for a drifting Maxwellian with
normal drift `a` (into the domain) and normal thermal speed `σ`:
`Γ/n0 = a/2·(1+erf(a/(σ√2))) + σ/√(2π)·exp(−a²/2σ²)`. (a=0 ⇒ σ/√(2π).)
"""
function flux_per_density(a::T, σ::T) where {T}
    aT, σT = _validated_flux_sampler_params(a, σ)
    σT == zero(T) && return max(aT, zero(T))
    d = aT / σT
    if aT < zero(T)
        b = -d
        return _negative_flux_per_density(σT, b)
    end
    return aT / 2 * (one(T) + _erf(d / sqrt(T(2)))) + σT / sqrt(T(2π)) * exp(-(d * d) / 2)
end

function _negative_flux_speed(rng, b::T, σ::T) where {T}
    tiny = nextfloat(zero(T))
    while true
        # Gamma(shape=2, rate=b) envelope: q(x) ∝ x*exp(-b*x).
        y = -log1p(-rand(rng, T)) - log1p(-rand(rng, T))
        x = max(y / b, tiny)
        rand(rng, T) < exp(-(x * x) / 2) || continue
        return max(σ * x, tiny)
    end
end

"""
    flux_speed(rng, a, σ)

Sample the inward normal speed `s>0` from `p(s) ∝ s·exp(−(s−a)²/2σ²)`. `a=0`
uses the exact Rayleigh inverse-CDF; otherwise inverse-CDF by bisection on the
closed-form (erf) cumulative.
"""
function flux_speed(rng, a::T, σ::T) where {T}
    aT, σT = _validated_flux_sampler_params(a, σ)
    σT == zero(T) && return max(aT, zero(T))
    if aT == 0
        # Rayleigh inverse CDF on 1−U ∈ (0,1]: rand ∈ [0,1) can return exactly 0
        # (probability 2⁻⁵³ per Float64 draw), and log(rand) = −Inf would emit an
        # infinite-speed particle; log1p(−U) keeps every sample finite (same law,
        # since 1−U is uniform on (0,1]).
        return σT * sqrt(-2 * log1p(-rand(rng, T)))
    end
    d = aT / σT
    if aT < zero(T)
        b = -d
        b >= one(T) && return _negative_flux_speed(rng, b, σT)
    end
    # If d overflowed, or d+14 rounds back to d, the entire thermal width is
    # below the representable spacing at a: every finite sample rounds to a.
    (!isfinite(d) || (aT > zero(T) && d + T(14) == d)) && return aT
    hi = d + T(14)
    Z = _flux_shape_integral(hi, d)
    (isfinite(Z) && Z > zero(T)) ||
        throw(ArgumentError("flux_speed: flux CDF normalization must be finite and positive"))
    U = rand(rng, T)
    lo = zero(T)
    _, drift_exponent = frexp(abs(d))
    for _ = 1:(64+max(0, drift_exponent))
        m = lo + (hi - lo) / 2
        (_flux_shape_integral(m, d) / Z < U) ? (lo = m) : (hi = m)
    end
    x = lo + (hi - lo) / 2
    speed = σT * x
    isfinite(speed) || throw(ArgumentError("flux_speed: sampled speed exceeds numeric range"))
    return max(speed, nextfloat(zero(T)))
end

@inline function _injection_increment(n0::Float64, flux::Float64, dt::Float64, weight::Float64)
    (iszero(n0) || iszero(flux) || iszero(dt)) && return 0.0
    nf, ne = frexp(n0)
    ff, fe = frexp(flux)
    df, de = frexp(dt)
    wf, we = frexp(weight)
    return ldexp((nf * ff * df) / wf, ne + fe + de - we)
end

"""
    inject_face_1d!(ps, rng, face_x, inward, n0, a, σ, ut, σt, dt, w, acc, nextid) -> Ninj

Append a flux-weighted batch of particles at the 1-D boundary `face_x`, moving
into the domain (`inward = +1` or `-1`). Normal velocity is `inward·flux_speed`;
transverse velocities are Maxwellian (drift `ut::NTuple{2}`, thermal `σt`).
The batch size is metered by the carried accumulator `acc` (a `Ref`) so the
mean injected number flux is exactly `n0·Γ/n0`; `nextid` (a `Ref`) supplies
unique particle ids. Returns the number injected this call.
"""
function inject_face_1d!(
    ps::ParticleSet{1,T},
    rng,
    face_x,
    inward::Integer,
    n0,
    a,
    σ,
    ut::NTuple{2},
    σt,
    dt,
    w,
    acc::Base.RefValue{Float64},
    nextid::Base.RefValue{UInt64},
) where {T}
    inward == 1 || inward == -1 || throw(ArgumentError("inward must be +1 or -1, got $inward"))
    face_xT = _require_finite_real("face_x", face_x, T)
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    aT = T(a)
    σT = T(σ)
    ut1 = _require_finite_real("ut[1]", ut[1], T)
    ut2 = _require_finite_real("ut[2]", ut[2], T)
    σtT = _require_finite_nonnegative_real("σt", σt, T)
    dtT = _validated_nonnegative_dt(T, dt; name = "inject_face_1d!")
    wT = _require_finite_positive_real("w", w, T)
    fpn0 = flux_per_density(aT, σT)
    acc0 = acc[]
    (isfinite(acc0) && acc0 >= 0) ||
        throw(ArgumentError("inject_face_1d!: acc[] must be finite and non-negative"))
    increment = _injection_increment(Float64(n0T), Float64(fpn0), Float64(dtT), Float64(wT))
    (isfinite(increment) && increment >= 0) ||
        throw(ArgumentError("inject_face_1d!: injected particle count must be finite"))
    acc1 = acc0 + increment
    (isfinite(acc1) && acc1 < Float64(typemax(Int))) ||
        throw(ArgumentError("inject_face_1d!: injected particle count exceeds Int capacity"))
    Ninj = floor(Int, acc1)
    if Ninj == 0
        acc[] = acc1
        return 0
    end

    firstid = nextid[]
    firstid > 0 || throw(ArgumentError("inject_face_1d!: nextid[] must be positive"))
    ninjU = UInt64(Ninj)
    ninjU <= typemax(UInt64) - firstid ||
        throw(ArgumentError("inject_face_1d!: particle id counter exhausted"))

    batch = ParticleSet{1,T}(Ninj; q = ps.q, m = ps.m)
    s_in = T(inward)
    @inbounds for k = 1:Ninj
        s = flux_speed(rng, aT, σT)
        (isfinite(s) && s >= 0) ||
            throw(ArgumentError("inject_face_1d!: flux sampler returned an invalid speed"))
        x = face_xT + s_in * s * dtT * rand(rng, T)   # fly-in within the swept slab
        vx = s_in * s
        vy = ut1 + σtT * randn(rng, T)
        vz = ut2 + σtT * randn(rng, T)
        (isfinite(x) && isfinite(vx) && isfinite(vy) && isfinite(vz)) ||
            throw(ArgumentError("inject_face_1d!: sampled particle state must be finite"))
        batch.x[1][k] = x
        batch.v[1][k] = vx
        batch.v[2][k] = vy
        batch.v[3][k] = vz
        batch.weight[k] = wT
        batch.id[k] = firstid + UInt64(k - 1)
        batch.tag[k] = UInt32(1)                                      # tag=1 ⇒ injected
    end
    append_particles!(ps, batch)
    acc[] = acc1 - Ninj
    nextid[] = firstid + ninjU
    return Ninj
end
