# abstract_integrator.jl — base type for time integrators (§5.3 integrators/).
abstract type AbstractIntegrator end

# Imaginary-axis stability limit |ω·Δt| of each magnetic mover (the largest ω·Δt
# that keeps the oscillatory whistler update non-growing): RK4 ≈ 2.83, leapfrog/CL
# ≈ 2, Crank–Nicolson unconditionally stable (∞). The :camcl bound is that of the
# PERSISTENT staggered cyclic leapfrog in camcl.jl (measured on the stiffest grid
# mode: bounded for ω·h ≤ 1.99, growing at ω·h = 2.0, non-finite at 2.05).
@inline function _integrator_stability_limit(integrator::Symbol)
    integrator === :rk4 && return 2.8
    integrator === :leapfrog && return 2.0
    integrator === :camcl && return 2.0      # persistent cyclic leapfrog, leapfrog bound
    integrator === :cn && return Inf         # Crank–Nicolson |M|≡1, unconditional
    throw(ArgumentError("integrator must be :rk4, :leapfrog, :camcl, or :cn, got :$integrator"))
end

# Real-axis stability limit ν·Δt_B of each magnetic mover for the explicitly
# integrated resistive/hyper-resistive decay rate ν(k) = η k² + ηH k⁴ (§6.2 Ohm
# terms η·J − ηH·∇²J): RK4's negative-real-axis interval is |ν·Δt| ≤ 2.785;
# Crank–Nicolson (A-stable) is unconditional. The leapfrog/CL family has an EMPTY
# real-axis interval in exact arithmetic (its substep map has det = 1, so any
# damping e^{−νh} of the physical mode grows the computational mode by e^{+νh});
# the persistent-CL implementation (camcl.jl) arrests that growth by re-syncing
# its two copies when they drift, measured bounded up to ν·h = 1.0 — we budget a
# conservative ν·h ≤ 0.25 per substep so re-syncs stay ≳ 10 steps apart and the
# resolved decay per substep stays accurate.
@inline function _integrator_real_axis_limit(integrator::Symbol)
    integrator === :rk4 && return 2.785
    integrator === :leapfrog && return 0.25
    integrator === :camcl && return 0.25
    integrator === :cn && return Inf
    throw(ArgumentError("integrator must be :rk4, :leapfrog, :camcl, or :cn, got :$integrator"))
end

# Whistler limit for the leapfrog/CL family when η/ηH damping is active: each
# drift-triggered re-sync of the CL copies acts like a restart of the leapfrog,
# whose per-restart amplification of a whistler mode grows steeply with ω·h
# (measured per-step |m|−1 with a forced restart EVERY step: 7e-6 at ω·h = 0.2,
# 2e-5 at 0.25, 6e-5 at 0.3, 4e-4 at 0.5, 8e-3 at 1.0). With damping present the
# re-syncs fire every ~5–30 steps, so capping ω·h at 0.8·0.3 = 0.24 keeps the
# injected growth ≲ 1e-5 per step; without damping the copies never drift and the
# full leapfrog bound 2.0 applies.
const _CL_DAMPED_WHISTLER_LIMIT = 0.3

"""
    recommended_dt(g; NB=1, integrator=:rk4, safety=0.8, d_i=1.0, Ω_ci=1.0,
                   η=0.0, ηH=0.0) -> Float64

Conservative recommended particle timestep from the whistler CFL (§10.3). The Hall
whistler branch `ω_W(k) = ½[√((k d_i)⁴+4(k d_i)²)+(k d_i)²]·Ω_ci` is stiffest at the
grid's largest representable wavenumber — the SPECTRAL-CORNER mode
`k_max = √(Σ_d (π/Δx_d)²)`, not the single-axis Nyquist `π/min(Δx)`: with the
background field along the corner direction `ω_W(k_max)` is up to `D`× the
single-axis value, so in 2-D/3-D the axis formula under-resolves the CFL and the
subcycle blows up at its "recommended" step (in 1-D the two coincide). The
magnetic substep must satisfy `ω_W(k_max)·Δt_B < C(integrator)`, so with `NB`
magnetic substeps per particle step (`Δt_B = Δt_p/NB`, §10.4) the recommended
particle step is

    Δt_p = safety · NB · C(integrator) / ω_W(k_max),

capped by the Boris gyro-accuracy limit `Ω_ci·Δt ≲ 0.3`. If the model carries
resistivity `η` or hyper-resistivity `ηH` (the same normalized coefficients passed
to [`HybridModel`](@ref)), the subcycle also integrates the REAL decay rate
`ν(k) = η k² + ηH k⁴` explicitly, so `Δt_p` is additionally capped by

    Δt_p ≤ safety · NB · C_real(integrator) / ν(k_max)

with the integrator's real-axis limit `C_real` (RK4 ≈ 2.785; CN unconditional;
the leapfrog/CL family tolerates essentially no real-axis damping — it gets a
conservative `ν·Δt_B ≤ 0.25` budget AND a reduced whistler constant `C = 0.3`,
see `_integrator_real_axis_limit`; prefer `:rk4` or `:cn` when `η`/`ηH` matter).
This reduces to the commonly-cited Hybrid-VPIC estimate
`Ω_ci Δt < (1/π)(Δx/d_i)²` at large `k_max·d_i` in 1-D. `integrator ∈ (:rk4,
:leapfrog, :camcl, :cn)`; `:cn` is unconditionally stable on the whistler and the
resistive decay, so only the gyro limit applies.
"""
function recommended_dt(
    g::FourierGrid{D,T};
    NB::Integer = 1,
    integrator::Symbol = :rk4,
    safety::Real = 0.8,
    d_i::Real = 1.0,
    Ω_ci::Real = 1.0,
    η::Real = 0.0,
    ηH::Real = 0.0,
) where {D,T}
    NB >= 1 || throw(ArgumentError("NB must be ≥ 1, got $NB"))
    (isfinite(safety) && 0 < safety <= 1) ||
        throw(ArgumentError("safety must be in (0, 1], got $safety"))
    (isfinite(d_i) && d_i > 0) || throw(ArgumentError("d_i must be finite and > 0, got $d_i"))
    (isfinite(Ω_ci) && Ω_ci > 0) || throw(ArgumentError("Ω_ci must be finite and > 0, got $Ω_ci"))
    (isfinite(η) && η >= 0) || throw(ArgumentError("η must be finite and ≥ 0, got $η"))
    (isfinite(ηH) && ηH >= 0) || throw(ArgumentError("ηH must be finite and ≥ 0, got $ηH"))
    kmax = sqrt(sum(dx -> abs2(π / dx), g.dx))          # spectral-corner |k| (stiffest mode)
    K = kmax * d_i
    ω_w = 0.5 * (sqrt(K^4 + 4 * K^2) + K^2) * Ω_ci      # whistler ω_W(k_max)
    ν = η * kmax^2 + ηH * kmax^4                        # explicit real-axis decay rate at k_max
    C = _integrator_stability_limit(integrator)
    if ν > 0 && (integrator === :leapfrog || integrator === :camcl)
        C = _CL_DAMPED_WHISTLER_LIMIT                   # damped CL: re-syncs inject restart error
    end
    dt_whistler = isfinite(C) ? safety * NB * C / ω_w : T(Inf)
    dt_gyro = safety * 0.3 / Ω_ci                       # Boris gyro-accuracy limit
    Creal = _integrator_real_axis_limit(integrator)
    dt_real = ν > 0 ? safety * NB * Creal / ν : T(Inf)  # real-axis (η/ηH) limit
    return float(min(dt_whistler, dt_gyro, dt_real))
end
