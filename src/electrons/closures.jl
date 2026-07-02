# hybrid.jl — the massless-electron hybrid field model in Ω_ci-normalized units
# (§6–7 of the checklist). Proton-only quasineutrality n_e = n_i = n; the
# electric field is algebraic (generalized Ohm's law), only B is dynamic.
#
# Normalized equations:
#   J = ∇×B
#   u_e = u_i − J/n
#   E = −u_i×B + (J×B)/n − ∇p_e/n + η J          (hyperresistivity optional, off here)
#   ∂B/∂t = −∇×E
# Boris with q/m = 1 advances protons in these E, B.

# ---------------------------------------------------------------- electron closures

abstract type ElectronClosure end

"Isothermal electrons: p_e = T_e · n  (T_e the normalized electron temperature)."
struct IsothermalElectrons{T} <: ElectronClosure
    Te::T
end

"Polytropic electrons: p_e = p_{e0} (n/n_0)^{γ_e}."
struct PolytropicElectrons{T} <: ElectronClosure
    pe0::T
    n0::T
    γ::T
end

function IsothermalElectrons(Te::Real)
    T = float(typeof(Te))
    return IsothermalElectrons{T}(_require_finite_nonnegative_real("Te", Te, T))
end

function PolytropicElectrons(pe0::Real, n0::Real, γ::Real)
    T = float(promote_type(typeof(pe0), typeof(n0), typeof(γ)))
    pe0T = _require_finite_nonnegative_real("pe0", pe0, T)
    n0T = _require_finite_positive_real("n0", n0, T)
    γT = _require_finite_positive_real("γ", γ, T)
    return PolytropicElectrons{T}(pe0T, n0T, γT)
end

electron_pressure!(pe, n, c::IsothermalElectrons) = (@. pe = c.Te * n; pe)
electron_pressure!(pe, n, c::PolytropicElectrons) = (@. pe = c.pe0 * (n / c.n0)^c.γ; pe)

"Adiabatic index of a closure (for the electron internal-energy budget)."
closure_gamma(::IsothermalElectrons) = 1.0
closure_gamma(c::PolytropicElectrons) = float(c.γ)

"""
    CGLElectrons(p_perp0, p_par0, n0, B0)

Anisotropic **double-adiabatic (Chew–Goldberger–Low)** electron closure. The gyrotropic
electron pressures follow the two CGL adiabatic invariants `p_⊥/(nB)=const` (frozen `μ`)
and `p_∥ B²/n³=const` (frozen `J_∥`), so from reference values `p_⊥0, p_∥0` at `n0, B0`:

    p_⊥(n,B) = C_⊥ · n · B,      C_⊥ = p_⊥0 / (n0 B0)
    p_∥(n,B) = C_∥ · n³ / B²,    C_∥ = p_∥0 B0² / n0³

The electron pressure is then the gyrotropic tensor `P_e = p_⊥ I + (p_∥ − p_⊥) bb`,
`b = B/|B|`, so the generalized Ohm's law carries `−(∇·P_e)/n` in place of the isotropic
`−∇p_e/n`. In the isotropic limit `p_⊥0 = p_∥0 = p_e0` this reduces exactly to a scalar
pressure. Unlike the scalar closures it depends on `B` (both magnitude and
direction), so the force `∇·P_e` is recomputed each B-subcycle stage rather than
frozen once per step — freezing the field direction would misalign the stress and
drive a numerical instability.

**Stability.** The pure double-adiabatic closure is well-posed and stable in the
firehose-/mirror-**stable** regime (mild anisotropy, `p_∥ > p_⊥`, or isotropic). In a
strongly **mirror-unstable** state (`β_⊥(p_⊥/p_∥ − 1) > 1`, i.e. `p_⊥ ≫ p_∥`) it is
ill-posed — the fluid mirror mode grows without a short-wavelength cutoff (no kinetic /
FLR stabilization) and the run blows up; that is a known property of the CGL model, not
a solver bug. Use it in the stable regime, or add a firehose/mirror pressure limiter or
a Landau-fluid closure for the unstable regime.
"""
struct CGLElectrons{T} <: ElectronClosure
    Cperp::T   # p_⊥ / (n B)
    Cpar::T    # p_∥ B² / n³
end

function CGLElectrons(p_perp0::Real, p_par0::Real, n0::Real, B0::Real)
    T = float(promote_type(typeof(p_perp0), typeof(p_par0), typeof(n0), typeof(B0)))
    pp0 = _require_finite_nonnegative_real("p_perp0", p_perp0, T)
    ppa0 = _require_finite_nonnegative_real("p_par0", p_par0, T)
    n0T = _require_finite_positive_real("n0", n0, T)
    B0T = _require_finite_positive_real("B0", B0, T)
    return CGLElectrons{T}(pp0 / (n0T * B0T), ppa0 * B0T^2 / n0T^3)
end

"Perpendicular CGL electron pressure `p_⊥ = C_⊥ n |B|` at density `n`, field magnitude `Bmag`."
@inline cgl_pperp(c::CGLElectrons{T}, n, Bmag) where {T} = c.Cperp * T(n) * T(Bmag)
"Parallel CGL electron pressure `p_∥ = C_∥ n³ / |B|²`."
@inline cgl_ppar(c::CGLElectrons{T}, n, Bmag) where {T} = c.Cpar * T(n)^3 / T(Bmag)^2

closure_gamma(::CGLElectrons) = 5 / 3      # effective γ; the anisotropic budget is separate

# CGL has no scalar pressure — the Ohm's-law path dispatches to anisotropic_pressure_force!
# instead of _ohm_prep!/electron_pressure!. This method makes accidental scalar use explicit.
electron_pressure!(::Any, ::Any, ::CGLElectrons) = throw(
    ArgumentError(
        "CGLElectrons is gyrotropic and has no scalar electron pressure; use anisotropic_pressure_force!",
    ),
)

"Whether a closure produces a gyrotropic (anisotropic) pressure needing `∇·P_e` in Ohm's law."
is_anisotropic(::ElectronClosure) = false
is_anisotropic(::CGLElectrons) = true

# ---------------------------------------------------------------- model + state
