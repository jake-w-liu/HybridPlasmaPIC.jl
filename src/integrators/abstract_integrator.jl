# abstract_integrator.jl — base type for time integrators (§5.3 integrators/).
abstract type AbstractIntegrator end

# Imaginary-axis stability limit |ω·Δt| of each magnetic mover (the largest ω·Δt
# that keeps the oscillatory whistler update non-growing): RK4 ≈ 2.83, leapfrog/CL
# ≈ 2, Crank–Nicolson unconditionally stable (∞).
@inline function _integrator_stability_limit(integrator::Symbol)
    integrator === :rk4 && return 2.8
    integrator === :leapfrog && return 2.0
    integrator === :camcl && return 2.0      # cyclic leapfrog, leapfrog-class bound
    integrator === :cn && return Inf         # Crank–Nicolson |M|≡1, unconditional
    throw(ArgumentError("integrator must be :rk4, :leapfrog, :camcl, or :cn, got :$integrator"))
end

"""
    recommended_dt(g; NB=1, integrator=:rk4, safety=0.8, d_i=1.0, Ω_ci=1.0) -> Float64

Conservative recommended particle timestep from the whistler CFL (§10.3). The Hall
whistler branch `ω_W(k) = ½[√((k d_i)⁴+4(k d_i)²)+(k d_i)²]·Ω_ci` is stiffest at the
grid's maximum resolved wavenumber `k_max = π/min(Δx)`. The magnetic substep must
satisfy `ω_W(k_max)·Δt_B < C(integrator)`, so with `NB` magnetic substeps per
particle step (`Δt_B = Δt_p/NB`, §10.4) the recommended particle step is

    Δt_p = safety · NB · C(integrator) / ω_W(k_max),

capped by the Boris gyro-accuracy limit `Ω_ci·Δt ≲ 0.3`. This reduces to the
commonly-cited Hybrid-VPIC estimate `Ω_ci Δt < (1/π)(Δx/d_i)²` at large `k_max·d_i`.
`integrator ∈ (:rk4, :leapfrog, :camcl, :cn)`; `:cn` is unconditionally stable on the
whistler, so only the gyro limit applies.
"""
function recommended_dt(
    g::FourierGrid{D,T};
    NB::Integer = 1,
    integrator::Symbol = :rk4,
    safety::Real = 0.8,
    d_i::Real = 1.0,
    Ω_ci::Real = 1.0,
) where {D,T}
    NB >= 1 || throw(ArgumentError("NB must be ≥ 1, got $NB"))
    (isfinite(safety) && 0 < safety <= 1) ||
        throw(ArgumentError("safety must be in (0, 1], got $safety"))
    (isfinite(d_i) && d_i > 0) || throw(ArgumentError("d_i must be finite and > 0, got $d_i"))
    (isfinite(Ω_ci) && Ω_ci > 0) || throw(ArgumentError("Ω_ci must be finite and > 0, got $Ω_ci"))
    dxmin = minimum(g.dx)
    kmax = π / dxmin                                    # grid Nyquist (max resolved wavenumber)
    K = kmax * d_i
    ω_w = 0.5 * (sqrt(K^4 + 4 * K^2) + K^2) * Ω_ci      # whistler ω_W(k_max)
    C = _integrator_stability_limit(integrator)
    dt_whistler = isfinite(C) ? safety * NB * C / ω_w : T(Inf)
    dt_gyro = safety * 0.3 / Ω_ci                       # Boris gyro-accuracy limit
    return float(min(dt_whistler, dt_gyro))
end
