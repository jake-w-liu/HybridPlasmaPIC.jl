# raycon_types.jl — configuration structs (ports of initCnst/initPlasma/data.m).

"""
    RayconConstants()

Physical constants in SI units, with the EXACT values of the upstream
`initCnst.m` (not CODATA-refreshed) so that ported results are bitwise
comparable with the original MATLAB code.
"""
struct RayconConstants
    c::Float64      # speed of light [m/s]
    e::Float64      # elementary charge [C]
    mp::Float64     # proton mass [kg]
    eps0::Float64   # vacuum permittivity [F/m]
end

RayconConstants() = RayconConstants(2.9979e8, 1.6022e-19, 1.6726e-27, 8.8542e-12)

"""
    SolovevEquilibrium(; b0, r0, q0, iaspr, elong)

Analytic Solovev tokamak equilibrium (`solovev.m`): on-axis field `b0` [T],
major radius `r0` [m], on-axis safety factor `q0`, inverse aspect ratio
`iaspr = a/r0`, elongation `elong`. The edge flux normalization is
`psin = b0/(2 q0) · elong · (r0·iaspr)²` (from `ray.m 'auxval'`).
"""
struct SolovevEquilibrium
    b0::Float64
    r0::Float64
    q0::Float64
    iaspr::Float64
    elong::Float64
    psin::Float64
end

function SolovevEquilibrium(; b0::Real, r0::Real, q0::Real, iaspr::Real, elong::Real)
    for (name, v) in (("b0", b0), ("r0", r0), ("q0", q0), ("iaspr", iaspr), ("elong", elong))
        (isfinite(v) && v > 0) || throw(ArgumentError("$name must be finite and positive"))
    end
    psin = 0.5 * b0 / q0 * elong * (r0 * iaspr)^2
    return SolovevEquilibrium(b0, r0, q0, iaspr, elong, psin)
end

"minor-radius scale `ρ_a = elong·iaspr·r0` used to bracket flux-surface roots (mapFlux.m)"
rho_edge(eq::SolovevEquilibrium) = eq.elong * eq.iaspr * eq.r0

"""
    RayconProblem(; eq, amass, acharge, n0, na, nb, t0, ta, tb, freq, kphi,
                    model=:cld2x2, cnst=RayconConstants())

Full RAYCON problem definition (`initPlasma.m` + `data.m`): the equilibrium,
the plasma species — atomic masses `amass` [proton masses], charges `acharge`
[e], on-axis densities `n0` [m⁻³] with profile `n = n0·(1−na·s²)^nb`, on-axis
temperatures `t0` [keV] with profile `T = t0·(1−ta·s²)^nb` (upstream quirk:
the temperature exponent uses `nb`, not `tb`; `tb` is retained as dead
configuration for parity) — the antenna frequency `freq` [Hz], the constant
toroidal wavenumber `kphi` [1/m] (`kant(2)`), and the dispersion `model`
(`:cld2x2`, `:cld3x3`, or `:msw1x1`).
"""
struct RayconProblem
    cnst::RayconConstants
    eq::SolovevEquilibrium
    amass::Vector{Float64}
    acharge::Vector{Float64}
    n0::Vector{Float64}
    na::Vector{Float64}
    nb::Vector{Float64}
    t0::Vector{Float64}
    ta::Vector{Float64}
    tb::Vector{Float64}
    freq::Float64
    omega::Float64
    kphi::Float64
    model::Symbol
end

function RayconProblem(;
    eq::SolovevEquilibrium,
    amass::AbstractVector{<:Real},
    acharge::AbstractVector{<:Real},
    n0::AbstractVector{<:Real},
    na::AbstractVector{<:Real},
    nb::AbstractVector{<:Real},
    t0::AbstractVector{<:Real},
    ta::AbstractVector{<:Real},
    tb::AbstractVector{<:Real},
    freq::Real,
    kphi::Real,
    model::Symbol = :cld2x2,
    cnst::RayconConstants = RayconConstants(),
)
    ns = length(amass)
    ns >= 1 || throw(ArgumentError("at least one species is required"))
    for (name, v) in (
        ("acharge", acharge),
        ("n0", n0),
        ("na", na),
        ("nb", nb),
        ("t0", t0),
        ("ta", ta),
        ("tb", tb),
    )
        length(v) == ns || throw(
            DimensionMismatch(
                "$name must have one entry per species (got $(length(v)), expected $ns)",
            ),
        )
        all(isfinite, v) || throw(ArgumentError("$name must be finite"))
    end
    all(m -> isfinite(m) && m > 0, amass) || throw(ArgumentError("amass must be positive"))
    all(>=(0), n0) || throw(ArgumentError("densities n0 must be non-negative"))
    all(>(0), nb) || throw(
        ArgumentError("profile exponents nb must be positive (they divide the log-derivatives)"),
    )
    (isfinite(freq) && freq > 0) || throw(ArgumentError("freq must be finite and positive"))
    isfinite(kphi) || throw(ArgumentError("kphi must be finite"))
    model in (:cld2x2, :cld3x3, :msw1x1) ||
        throw(ArgumentError("model must be :cld2x2, :cld3x3 or :msw1x1 (got $model)"))
    return RayconProblem(
        cnst,
        eq,
        Float64.(collect(amass)),
        Float64.(collect(acharge)),
        Float64.(collect(n0)),
        Float64.(collect(na)),
        Float64.(collect(nb)),
        Float64.(collect(t0)),
        Float64.(collect(ta)),
        Float64.(collect(tb)),
        Float64(freq),
        2π * Float64(freq),
        Float64(kphi),
        model,
    )
end

"""
    cmod_parameters(; model=:cld2x2, kphi=-10.0)

The Alcator C-Mod ICRF minority mode-conversion case (`data.m` `'cmod'` with
the `main.m` overrides): B0 = 7.9 T, R0 = 0.67 m, e/D/He3 plasma, f = 80 MHz,
`kant = [-31.5, kphi, 0]`. This is the upstream reference scenario for
conversion between the fast magnetosonic and ion-hybrid waves.
"""
function cmod_parameters(; model::Symbol = :cld2x2, kphi::Real = -10.0)
    eq = SolovevEquilibrium(; b0 = 7.9, r0 = 0.67, q0 = 2.0, iaspr = 0.22 / 0.67, elong = 1.6)
    return RayconProblem(;
        eq,
        amass = [1 / 1836, 2.0, 3.0],
        acharge = [-1.0, 1.0, 2.0],
        n0 = [10.0, 5.2, 2.4] .* 1e19,
        na = [1.0, 0.7, 0.7],
        nb = [3.0, 3.0, 3.0],
        t0 = [3.0, 3.0, 3.0],
        ta = [1.0, 1.0, 1.0],
        tb = [1.0, 1.0, 1.0],
        freq = 80.0e6,
        kphi,
        model,
    )
end
