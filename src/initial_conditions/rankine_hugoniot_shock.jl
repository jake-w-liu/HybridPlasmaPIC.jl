# shock.jl — independent coplanar MHD Rankine–Hugoniot solver (§19.1).
#
# Shock frame, normal = x, coplanar (B,u in the x–y plane). Given the upstream
# state and adiabatic index γ, solve for the downstream state across the shock.
# Construction: parametrize by compression X = ρ₂/ρ₁; mass, induction,
# tangential-momentum, and normal-momentum jumps are satisfied EXACTLY by the
# closed-form downstream(X); the energy jump is a single residual root-found for
# X. The solver returns the downstream state and the (normalized) residuals of
# all six conservation laws — independent of the time-domain solver, so it can
# serve as a verification oracle.

"Coplanar MHD state in the shock frame: density, normal/tangential velocity, pressure, normal/tangential B."
struct MHDState{T}
    ρ::T
    ux::T
    uy::T
    p::T
    Bn::T
    Bt::T
end

# closed-form downstream state for a trial compression X (exact for mass,
# induction, tangential & normal momentum)
function _downstream(up::MHDState{T}, X::T, μ0::T) where {T}
    G = up.ρ * up.ux                      # mass flux ρuₓ (conserved)
    a = up.Bn^2 / (μ0 * G)                # = v_Ax²/uₓ  (0 for a perpendicular/hydro shock)
    ρ2 = X * up.ρ
    u2x = up.ux / X
    B2t = up.Bt * (a - up.ux) / (a - up.ux / X)
    u2y = up.uy + up.Bn * (B2t - up.Bt) / (μ0 * G)
    p2 = up.p + up.ρ * up.ux^2 * (one(T) - one(T) / X) + (up.Bt^2 - B2t^2) / (2μ0)
    return MHDState{T}(ρ2, u2x, u2y, p2, up.Bn, B2t)
end

# normal energy flux F_E,n
function _energy_flux(s::MHDState{T}, γ::T, μ0::T) where {T}
    u2 = s.ux^2 + s.uy^2
    B2 = s.Bn^2 + s.Bt^2
    uB = s.ux * s.Bn + s.uy * s.Bt
    return s.ux * (T(0.5) * s.ρ * u2 + γ / (γ - one(T)) * s.p + B2 / μ0) - uB * s.Bn / μ0
end

# the six conservation-law normal fluxes (mass, Bn, induction, mom_n, mom_t, energy)
function _fluxes(s::MHDState{T}, γ::T, μ0::T) where {T}
    return (
        mass = s.ρ * s.ux,
        Bn = s.Bn,
        induction = s.ux * s.Bt - s.uy * s.Bn,
        mom_n = s.ρ * s.ux^2 + s.p + (s.Bt^2 - s.Bn^2) / (2μ0),
        mom_t = s.ρ * s.ux * s.uy - s.Bn * s.Bt / μ0,
        energy = _energy_flux(s, γ, μ0),
    )
end

"""
    rankine_hugoniot(up::MHDState, γ; μ0=1.0) -> (down, X, residuals)

Solve the coplanar MHD jump conditions for the compressive (X>1) downstream
state. Returns the downstream `MHDState`, the compression `X = ρ₂/ρ₁`, and a
NamedTuple of normalized residuals for all six conservation laws (≈0 at the
solution; energy is root-found, the rest are exact by construction). `X = 1`
(with the upstream returned) means no compressive shock was bracketed.
"""
function rankine_hugoniot(up::MHDState{T}, γ::Real; μ0::Real = 1.0) where {T}
    γT = T(γ)
    μ = T(μ0)
    up.ux > 0 || throw(ArgumentError("upstream normal velocity must be > 0 (shock frame)"))
    Fup = _fluxes(up, γT, μ)
    Renergy(X) = _energy_flux(_downstream(up, X, μ), γT, μ) - Fup.energy

    Xmax = (γT + one(T)) / (γT - one(T))          # strong-shock compression limit
    # scan (1, Xmax) for a sign change of the energy residual (X=1 is the trivial root)
    nscan = 512
    lo = one(T) + T(1e-7)
    hi = Xmax - T(1e-9)
    Xsol = one(T)
    found = false
    Rprev = Renergy(lo)
    Xprev = lo
    for i = 1:nscan
        Xc = lo + (hi - lo) * i / nscan
        Rc = Renergy(Xc)
        if Rprev == 0
            Xsol = Xprev
            found = true
            break
        elseif sign(Rc) != sign(Rprev)
            # bisect [Xprev, Xc]
            a, b, fa = Xprev, Xc, Rprev
            for _ = 1:200
                m = (a + b) / 2
                fm = Renergy(m)
                if fm == 0 || (b - a) < T(1e-14) * m
                    a = b = m
                    break
                end
                if sign(fm) != sign(fa)
                    b = m
                else
                    a = m
                    fa = fm
                end
            end
            Xsol = (a + b) / 2
            found = true
            break
        end
        Rprev = Rc
        Xprev = Xc
    end

    down = _downstream(up, Xsol, μ)
    Fdn = _fluxes(down, γT, μ)
    scale(name) = max(abs(getfield(Fup, name)), abs(getfield(Fdn, name)), eps(T))
    residuals = (
        mass = abs(Fdn.mass - Fup.mass) / scale(:mass),
        Bn = abs(Fdn.Bn - Fup.Bn) / scale(:Bn),
        induction = abs(Fdn.induction - Fup.induction) / scale(:induction),
        mom_n = abs(Fdn.mom_n - Fup.mom_n) / scale(:mom_n),
        mom_t = abs(Fdn.mom_t - Fup.mom_t) / scale(:mom_t),
        energy = abs(Fdn.energy - Fup.energy) / scale(:energy),
    )
    return (; down, X = found ? Xsol : one(T), residuals)
end

"""
    rh_branches(up, γ; μ0=1.0, nscan=2000) -> Vector{NamedTuple}

Branch tracking: find ALL compressive (X>1) Rankine–Hugoniot solutions by
scanning the energy-jump residual over (1, X_max] and bisecting every sign
change. Each returned `(; X, down, residuals)` satisfies all six jump conditions
(residuals ≈ 0). A perpendicular shock yields a single fast branch; oblique
upstream states can admit slow/intermediate/fast branches.
"""
function rh_branches(up::MHDState{T}, γ::Real; μ0::Real = 1.0, nscan::Int = 2000) where {T}
    γT = T(γ)
    μ = T(μ0)
    up.ux > 0 || throw(ArgumentError("upstream normal velocity must be > 0 (shock frame)"))
    Fup = _fluxes(up, γT, μ)
    R(X) = _energy_flux(_downstream(up, X, μ), γT, μ) - Fup.energy
    Xmax = (γT + one(T)) / (γT - one(T))
    lo = one(T) + T(1e-7)
    hi = Xmax - T(1e-9)
    branches = NamedTuple[]
    Rprev = R(lo)
    Xprev = lo
    for i = 1:nscan
        Xc = lo + (hi - lo) * i / nscan
        Rc = R(Xc)
        if Rprev != 0 && sign(Rc) != sign(Rprev)
            a, b, fa = Xprev, Xc, Rprev
            for _ = 1:200
                m = (a + b) / 2
                fm = R(m)
                if fm == 0 || (b - a) < T(1e-13) * m
                    a = b = m
                    break
                end
                sign(fm) != sign(fa) ? (b = m) : (a = m; fa = fm)
            end
            Xs = (a + b) / 2
            down = _downstream(up, Xs, μ)
            Fdn = _fluxes(down, γT, μ)
            scale(name) = max(abs(getfield(Fup, name)), abs(getfield(Fdn, name)), eps(T))
            res = (
                mass = abs(Fdn.mass - Fup.mass) / scale(:mass),
                Bn = abs(Fdn.Bn - Fup.Bn) / scale(:Bn),
                induction = abs(Fdn.induction - Fup.induction) / scale(:induction),
                mom_n = abs(Fdn.mom_n - Fup.mom_n) / scale(:mom_n),
                mom_t = abs(Fdn.mom_t - Fup.mom_t) / scale(:mom_t),
                energy = abs(Fdn.energy - Fup.energy) / scale(:energy),
            )
            push!(branches, (; X = Xs, down, residuals = res))
        end
        Rprev = Rc
        Xprev = Xc
    end
    return branches
end
