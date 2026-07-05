# RAYCON port: ICRF fast-wave ray in Alcator C-Mod with linear mode
# conversion — the reference scenario of the original MATLAB package (Tracy,
# Kaufman & Jaun), driven entirely in the package's Ω_ci-NORMALIZED units:
# lengths in d_i, wavenumbers in 1/d_i, frequencies in Ω_ci.
#
#   julia --project=. examples/raycon_cmod.jl

using HybridPlasmaPIC
using HybridPlasmaPIC.Raycon
using Printf

units = cmod_units()                 # B0 = 7.9 T, n0 = 1e20 m⁻³, proton mass
prob = cmod_parameters(units)        # the C-Mod case in normalized form
di = inertial_length(units)
Ωci = gyrofrequency(units)

@printf(
    "reference scales: d_i = %.3f cm, Ω_ci = %.3e rad/s (f_ci = %.1f MHz)\n",
    100di,
    Ωci,
    Ωci / 2π / 1e6,
)
@printf(
    "tokamak in d_i: R0 = %.1f, a = %.1f;  antenna ω = %.3f Ω_ci, kφ = %.3f d_i⁻¹\n",
    prob.eq.r0 / di,
    prob.eq.r0 * prob.eq.iaspr / di,
    prob.omega / Ωci,
    prob.kphi * di,
)

res = trace_rays(units, prob; s = 0.4, theta = 0.001, kr = -31.5 * di, kz = 0.0, sigma_span = 5e-2)

for (i, r) in enumerate(res.rays)
    @printf(
        "ray %d (parent %d): %5d points, status %-14s  end R = %.2f d_i  Z = %+.3f d_i\n",
        i,
        r.parent,
        length(r.sigma),
        r.status,
        r.y[1, end],
        r.y[2, end],
    )
end
for c in res.conversions
    cv = c.conversion
    @printf(
        "conversion on ray %d at σ=%.3e: η²=%.3f  τ=%.3f  |β|=%.3f  arg β=%+.3f rad\n",
        c.ray,
        c.sigma,
        cv.eta2,
        cv.tau,
        abs(cv.beta),
        angle(cv.beta),
    )
    @printf(
        "   saddle (R, Z, kR, kZ) = (%.2f d_i, %+.3f d_i, %.3f d_i⁻¹, %.3f d_i⁻¹)\n",
        cv.saddle...,
    )
end
isempty(res.conversions) && println("no mode conversion detected on this ray")
