# raycon_normalized.jl — Ω_ci-normalized interface to the RAYCON port.
#
# This is the UNIFIED, package-convention interface: every method here takes a
# `PlasmaUnits` (the package's SI↔normalized reference scales, verification/
# normalization.jl) as its first argument and speaks exclusively the hybrid
# solver's units — lengths in d_i, wavenumbers in 1/d_i, frequencies in Ω_ci,
# magnetic field in B0, densities in n0, temperatures in m_i·v_A²
# (= B0²/(μ0·n0), the pressure unit over the density unit). The ray parameter σ
# is invariant under the length rescaling (U is dimensionless), so σ values and
# spans are identical in both systems.
#
# Internally each call converts at the boundary and runs the SI core — the
# layer that is regression-pinned against the original MATLAB RAYCON — so the
# normalized interface is exactly as verified as the SI one: unit conversion is
# a pure multiplicative rescaling.
#
# Physics constants: by default the normalized constructor uses the PACKAGE's
# CODATA constants from `PlasmaUnits` (e, c, ε0 = 1/(μ0 c²), and `units.mi` as
# the `amass` reference mass) so results are self-consistent with the rest of
# HybridPlasmaPIC. Pass `cnst = RayconConstants()` to reproduce the upstream
# MATLAB constants instead (~1e-5 relative shifts; used by the parity tests).

"""
    RayconProblem(units::PlasmaUnits; r0, iaspr, elong, q0, b0=1.0,
                  amass, acharge, n0, na, nb, t0, ta, tb,
                  omega, kphi, model=:cld2x2, cnst=nothing)

Build a RAYCON problem from Ω_ci-NORMALIZED inputs (the package convention):
major radius `r0` [d_i], on-axis field `b0` [B0], per-species on-axis densities
`n0` [units.n0], on-axis temperatures `t0` [m_i·v_A²], antenna frequency
`omega` [Ω_ci], constant toroidal wavenumber `kphi` [1/d_i]. `amass` is in
units of the reference ion mass `units.mi`; `iaspr`, `elong`, `q0`, `acharge`,
`na`, `nb`, `ta`, `tb` are dimensionless as always. `cnst = nothing` (default)
derives the physics constants from `units` (package CODATA values); pass an
explicit [`RayconConstants`](@ref) for upstream-MATLAB constant parity.
"""
function RayconProblem(
    units::PlasmaUnits;
    r0::Real,
    iaspr::Real,
    elong::Real,
    q0::Real,
    b0::Real = 1.0,
    amass::AbstractVector{<:Real},
    acharge::AbstractVector{<:Real},
    n0::AbstractVector{<:Real},
    na::AbstractVector{<:Real},
    nb::AbstractVector{<:Real},
    t0::AbstractVector{<:Real},
    ta::AbstractVector{<:Real},
    tb::AbstractVector{<:Real},
    omega::Real,
    kphi::Real,
    model::Symbol = :cld2x2,
    cnst::Union{Nothing,RayconConstants} = nothing,
)
    di = inertial_length(units)
    Ωci = gyrofrequency(units)
    vA = alfven_speed(units)
    cnstv =
        cnst === nothing ?
        RayconConstants(units.c, units.e, units.mi, 1 / (units.mu0 * units.c^2)) : cnst
    (isfinite(omega) && omega > 0) ||
        throw(ArgumentError("omega must be finite and positive [Ω_ci]"))
    temp_scale = units.mi * vA^2 / (units.e * 1000)   # m_i·v_A² [J] → [keV]
    eq = SolovevEquilibrium(;
        b0 = b0 * units.B0,
        r0 = r0 * di,
        q0 = q0,
        iaspr = iaspr,
        elong = elong,
    )
    return RayconProblem(;
        eq,
        amass = collect(amass),
        acharge = collect(acharge),
        n0 = collect(n0) .* units.n0,
        na = collect(na),
        nb = collect(nb),
        t0 = collect(t0) .* temp_scale,
        ta = collect(ta),
        tb = collect(tb),
        freq = omega * Ωci / (2π),
        kphi = kphi / di,
        model,
        cnst = cnstv,
    )
end

"""
    cmod_parameters(units::PlasmaUnits; model=:cld2x2, kphi_di=nothing)

The C-Mod reference case expressed through the normalized interface. With
`cmod_units()` as the reference scales this reproduces `cmod_parameters()`
(up to the CODATA-vs-upstream constants; see the file header). `kphi_di`
defaults to the upstream −10 m⁻¹ converted to 1/d_i.
"""
function cmod_parameters(
    units::PlasmaUnits;
    model::Symbol = :cld2x2,
    kphi_di::Union{Nothing,Real} = nothing,
)
    di = inertial_length(units)
    Ωci = gyrofrequency(units)
    vA = alfven_speed(units)
    temp_scale = units.mi * vA^2 / (units.e * 1000)
    kphi = kphi_di === nothing ? -10.0 * di : Float64(kphi_di)
    return RayconProblem(
        units;
        r0 = 0.67 / di,
        iaspr = 0.22 / 0.67,
        elong = 1.6,
        q0 = 2.0,
        b0 = 7.9 / units.B0,
        amass = [1 / 1836, 2.0, 3.0] .* (1.6726e-27 / units.mi),
        acharge = [-1.0, 1.0, 2.0],
        n0 = [10.0, 5.2, 2.4] .* 1e19 ./ units.n0,
        na = [1.0, 0.7, 0.7],
        nb = [3.0, 3.0, 3.0],
        t0 = [3.0, 3.0, 3.0] ./ temp_scale,
        ta = [1.0, 1.0, 1.0],
        tb = [1.0, 1.0, 1.0],
        omega = 2π * 80.0e6 / Ωci,
        kphi,
        model,
    )
end

"""
    cmod_units() -> PlasmaUnits

Natural reference scales for the C-Mod case: `B0 = 7.9` T, `n0 = 1e20` m⁻³
(the on-axis electron density), proton reference mass.
"""
cmod_units() = PlasmaUnits(; n0 = 1.0e20, B0 = 7.9, mi = 1.6726e-27)

# normalized (r, z, kr, kz) → SI and back; σ, τ, β, η² are scale-invariant
@inline _state_to_SI(y, di::Float64) = [y[1] * di, y[2] * di, y[3] / di, y[4] / di]
@inline _state_to_norm(y, di::Float64) = [y[1] / di, y[2] / di, y[3] * di, y[4] * di]
@inline _tuple_to_norm(t::NTuple{4,Float64}, di::Float64) =
    (t[1] / di, t[2] / di, t[3] * di, t[4] * di)

function _normalize_conversion(cv::RayconConversion, di::Float64)
    return RayconConversion(
        _tuple_to_norm(cv.saddle, di),
        _tuple_to_norm(cv.incoming, di),
        _tuple_to_norm(cv.converted, di),
        _tuple_to_norm(cv.transmitted, di),
        cv.eta2,
        cv.tau,
        cv.beta,
        cv.eta2_estimate,
        cv.converged,
        cv.hyperbola_ok,
        cv.transmitted_ok,
    )
end

"""
    launch_ray(units::PlasmaUnits, prob; s, theta, kr, kz, m=0.0)

Normalized antenna launch: `kr`, `kz` in 1/d_i; returns `[r, z, kr, kz]` with
positions in d_i and wavenumbers in 1/d_i.
"""
function launch_ray(
    units::PlasmaUnits,
    prob::RayconProblem;
    s::Real,
    theta::Real,
    kr::Real,
    kz::Real,
    m::Real = 0.0,
)
    di = Float64(inertial_length(units))
    y = launch_ray(prob; s, theta, kr = kr / di, kz = kz / di, m)
    return _state_to_norm(y, di)
end

"""
    integrate_ray(units::PlasmaUnits, prob, y0, sigma0, sigma_end; kwargs...)

Normalized single-ray integration: `y0 = (r, z, kr, kz)` in (d_i, d_i, 1/d_i,
1/d_i); the returned `RayconTrace` carries positions in d_i and wavenumbers in
1/d_i. σ is scale-invariant, so spans and step controls keep their meaning
(`rtol`/`atol` act on the internal SI state — upstream tolerance parity).
"""
function integrate_ray(
    units::PlasmaUnits,
    prob::RayconProblem,
    y0::AbstractVector{<:Real},
    sigma0::Real,
    sigma_end::Real;
    kwargs...,
)
    length(y0) == 4 || throw(ArgumentError("state must be (r, z, kr, kz)"))
    di = Float64(inertial_length(units))
    tr = integrate_ray(prob, _state_to_SI(Float64.(y0), di), sigma0, sigma_end; kwargs...)
    yn = copy(tr.y)
    yn[1:2, :] ./= di
    yn[3:4, :] .*= di
    return RayconTrace(tr.sigma, yn, tr.mon2, tr.status, tr.zdot, tr.zddot)
end

"""
    trace_rays(units::PlasmaUnits, prob; s, theta, kr, kz, kwargs...)

Normalized mode-conversion ray tracing: identical to the SI
[`trace_rays`](@ref) but with `kr`, `kz` given in 1/d_i and every returned
trajectory and conversion point expressed in (d_i, 1/d_i). The dimensionless
conversion coefficients (η², τ, β) are unit-independent.
"""
function trace_rays(
    units::PlasmaUnits,
    prob::RayconProblem;
    s::Real,
    theta::Real,
    kr::Real,
    kz::Real,
    kwargs...,
)
    di = Float64(inertial_length(units))
    res = trace_rays(prob; s, theta, kr = kr / di, kz = kz / di, kwargs...)
    rays = map(res.rays) do r
        yn = copy(r.y)
        yn[1:2, :] ./= di
        yn[3:4, :] .*= di
        (; parent = r.parent, sigma = r.sigma, y = yn, status = r.status)
    end
    conversions = map(res.conversions) do c
        (; ray = c.ray, sigma = c.sigma, conversion = _normalize_conversion(c.conversion, di))
    end
    return (; rays, conversions)
end
