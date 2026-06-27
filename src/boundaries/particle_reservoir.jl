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

"""
    flux_per_density(a, σ)

Inward number flux per unit upstream density for a drifting Maxwellian with
normal drift `a` (into the domain) and normal thermal speed `σ`:
`Γ/n0 = a/2·(1+erf(a/(σ√2))) + σ/√(2π)·exp(−a²/2σ²)`. (a=0 ⇒ σ/√(2π).)
"""
function flux_per_density(a::T, σ::T) where {T}
    return a / 2 * (one(T) + _erf(a / (σ * sqrt(T(2))))) + σ / sqrt(T(2π)) * exp(-a^2 / (2σ^2))
end

"""
    flux_speed(rng, a, σ)

Sample the inward normal speed `s>0` from `p(s) ∝ s·exp(−(s−a)²/2σ²)`. `a=0`
uses the exact Rayleigh inverse-CDF; otherwise inverse-CDF by bisection on the
closed-form (erf) cumulative.
"""
function flux_speed(rng, a::T, σ::T) where {T}
    if a == 0
        return σ * sqrt(-2 * log(rand(rng, T)))
    end
    Z = _flux_integral(a + 14σ, a, σ)
    U = rand(rng, T)
    lo = zero(T)
    hi = a + 14σ
    for _ = 1:60
        m = (lo + hi) / 2
        (_flux_integral(m, a, σ) / Z < U) ? (lo = m) : (hi = m)
    end
    return (lo + hi) / 2
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
    aT = T(a)
    σT = T(σ)
    fpn0 = flux_per_density(aT, σT)
    acc[] += float(n0) * float(fpn0) * float(dt) / float(w)
    Ninj = floor(Int, acc[])
    acc[] -= Ninj
    Ninj == 0 && return 0
    s_in = T(inward)
    σtT = T(σt)
    @inbounds for _ = 1:Ninj
        s = flux_speed(rng, aT, σT)
        push!(ps.x[1], T(face_x) + s_in * s * T(dt) * rand(rng, T))   # fly-in within the swept slab
        push!(ps.v[1], s_in * s)
        push!(ps.v[2], T(ut[1]) + σtT * randn(rng, T))
        push!(ps.v[3], T(ut[2]) + σtT * randn(rng, T))
        push!(ps.weight, T(w))
        push!(ps.id, nextid[])
        nextid[] += one(UInt64)
        push!(ps.tag, UInt32(1))                                      # tag=1 ⇒ injected
    end
    return Ninj
end
