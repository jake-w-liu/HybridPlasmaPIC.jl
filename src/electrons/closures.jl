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

electron_pressure!(pe, n, c::IsothermalElectrons) = (@. pe = c.Te * n; pe)
electron_pressure!(pe, n, c::PolytropicElectrons) = (@. pe = c.pe0 * (n / c.n0)^c.γ; pe)

"Adiabatic index of a closure (for the electron internal-energy budget)."
closure_gamma(::IsothermalElectrons) = 1.0
closure_gamma(c::PolytropicElectrons) = float(c.γ)

# ---------------------------------------------------------------- model + state
