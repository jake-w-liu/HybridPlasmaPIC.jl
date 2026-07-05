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

function _require_valid_gamma(γ::Real, ::Type{T}) where {T}
    γT = T(γ)
    isfinite(γT) && γT > one(T) || throw(ArgumentError("γ must be finite and > 1"))
    return γT
end

function _require_valid_rh_inputs(up::MHDState{T}, μ0::Real) where {T}
    μ = T(μ0)
    isfinite(μ) && μ > zero(T) || throw(ArgumentError("μ0 must be finite and positive"))
    isfinite(up.ρ) && up.ρ > zero(T) ||
        throw(ArgumentError("upstream density must be finite and positive"))
    isfinite(up.ux) && up.ux > zero(T) ||
        throw(ArgumentError("upstream normal velocity must be finite and > 0 (shock frame)"))
    isfinite(up.uy) || throw(ArgumentError("upstream tangential velocity must be finite"))
    isfinite(up.p) && up.p >= zero(T) ||
        throw(ArgumentError("upstream pressure must be finite and non-negative"))
    isfinite(up.Bn) || throw(ArgumentError("upstream normal magnetic field must be finite"))
    isfinite(up.Bt) || throw(ArgumentError("upstream tangential magnetic field must be finite"))
    return μ
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

# absolute-magnitude scales (Σ|term|) of the six normal fluxes — the residual
# denominators. A flux can vanish by cancellation while its terms are O(1)
# (e.g. the induction and tangential-momentum fluxes across a switch-on shock),
# so a meaningful relative residual must be measured against the constituent
# term magnitude, not against the (near-zero) flux value itself.
function _flux_scales(s::MHDState{T}, γ::T, μ0::T) where {T}
    u2 = s.ux^2 + s.uy^2
    B2 = s.Bn^2 + s.Bt^2
    return (
        mass = abs(s.ρ * s.ux),
        Bn = abs(s.Bn),
        induction = abs(s.ux * s.Bt) + abs(s.uy * s.Bn),
        mom_n = abs(s.ρ * s.ux^2) + abs(s.p) + B2 / (2μ0),
        mom_t = abs(s.ρ * s.ux * s.uy) + abs(s.Bn * s.Bt) / μ0,
        energy = abs(s.ux) * (T(0.5) * s.ρ * u2 + γ / (γ - one(T)) * abs(s.p) + B2 / μ0) +
                 abs((s.ux * s.Bn + s.uy * s.Bt) * s.Bn) / μ0,
    )
end

# normalized residuals of the six conservation laws across the up → down jump:
# |F₂ − F₁| over the largest constituent-term magnitude of that flux on either
# side (≥ |F| by the triangle inequality, so never stricter than the plain
# |F|-normalized residual, and still meaningful for cancelling fluxes).
function _rh_residuals(up::MHDState{T}, down::MHDState{T}, γ::T, μ0::T) where {T}
    Fup = _fluxes(up, γ, μ0)
    Fdn = _fluxes(down, γ, μ0)
    Sup = _flux_scales(up, γ, μ0)
    Sdn = _flux_scales(down, γ, μ0)
    scale(name) = max(getfield(Sup, name), getfield(Sdn, name), eps(T))
    return (
        mass = abs(Fdn.mass - Fup.mass) / scale(:mass),
        Bn = abs(Fdn.Bn - Fup.Bn) / scale(:Bn),
        induction = abs(Fdn.induction - Fup.induction) / scale(:induction),
        mom_n = abs(Fdn.mom_n - Fup.mom_n) / scale(:mom_n),
        mom_t = abs(Fdn.mom_t - Fup.mom_t) / scale(:mom_t),
        energy = abs(Fdn.energy - Fup.energy) / scale(:energy),
    )
end

# Bt₁ = 0 bifurcation (switch-on shock). For a field-aligned upstream the
# induction + tangential-momentum jumps factor as Bt₂·(u₂ₓ − Bn²/(μ0 G)) = 0:
# besides the gasdynamic family (Bt₂ = 0, any X) there is an isolated branch
# pinned at u₂ₓ = Bn²/(μ0 G), i.e. X = u₁ₓ/u₂ₓ = M_An² (downstream exactly
# Alfvénic, M_An₂ = 1), where the energy jump closes in closed form:
#   Bt₂² = 2μ0(γ−1)·(u₁ₓ−u₂ₓ)/u₂ₓ · [γ/(γ−1)·(G·u₂ₓ − p₁) − G·(u₁ₓ+u₂ₓ)/2],
#   u₂y = u₁y + Bn·Bt₂/(μ0 G),   p₂ = p₁ + G·(u₁ₓ−u₂ₓ) − Bt₂²/(2μ0).
# Bt₂² > 0 delimits the switch-on window (1 < M_An² < (γ+1)/(γ−1) as β → 0).
# Sign convention: Bt₂ ≥ 0 — the jump conditions fix only Bt₂²; the switch-on
# tangential field direction is arbitrary. Returns `nothing` when the upstream
# is not (near-)field-aligned or the window is empty.
function _switch_on_downstream(up::MHDState{T}, γ::T, μ0::T) where {T}
    Bscale = max(abs(up.Bn), sqrt(μ0 * up.ρ) * up.ux)
    abs(up.Bt) <= 4 * eps(T) * Bscale || return nothing
    up.Bn == 0 && return nothing
    G = up.ρ * up.ux
    u2x = up.Bn^2 / (μ0 * G)
    X = up.ux / u2x                                    # = M_An²
    isfinite(X) && X > one(T) || return nothing
    Bt2sq =
        2μ0 * (γ - one(T)) * (up.ux - u2x) / u2x *
        (γ / (γ - one(T)) * (G * u2x - up.p) - G * (up.ux + u2x) / 2)
    Bt2sq > zero(T) || return nothing
    Bt2 = sqrt(Bt2sq)
    p2 = up.p + G * (up.ux - u2x) - Bt2sq / (2μ0)
    p2 > zero(T) || return nothing
    u2y = up.uy + up.Bn * Bt2 / (μ0 * G)
    return MHDState{T}(X * up.ρ, u2x, u2y, p2, up.Bn, Bt2), X
end

# Evolutionary (Alfvén-point) ordering: a physical fast-family shock keeps
# M_An ≥ 1 on both sides, a slow-family shock M_An ≤ 1 on both; a branch that
# CROSSES the Alfvén point (1→3, 1→4, 2→4) is non-evolutionary. Switch-on sits
# on the fast-family boundary (M_An₂ = 1 exactly), hence the tolerance.
function _alfven_ordered(up::MHDState{T}, down::MHDState{T}, μ0::T) where {T}
    up.Bn == 0 && return true                          # no normal Alfvén point
    M1sq = μ0 * up.ρ * up.ux^2 / up.Bn^2
    M2sq = μ0 * down.ρ * down.ux^2 / down.Bn^2
    tol = sqrt(eps(T))
    return M1sq >= one(T) ? M2sq >= one(T) - tol : M2sq <= one(T) + tol
end

"""
    rankine_hugoniot(up::MHDState, γ; μ0=1.0) -> (down, X, residuals)

Solve the coplanar MHD jump conditions for the compressive (X>1) downstream
state. Returns the downstream `MHDState`, the compression `X = ρ₂/ρ₁`, and a
NamedTuple of normalized residuals for all six conservation laws (≈0 at the
solution; energy is root-found, the rest are exact by construction). `X = 1`
(with the upstream returned) means no compressive shock was bracketed.

A field-aligned upstream (`Bt = 0`) is bifurcation-handled: inside the
switch-on window (`1 < M_An² < (γ+1)/(γ−1)` at low β) the gasdynamic root
crosses the Alfvén point and is non-evolutionary, so the evolutionary
switch-on branch (`X = M_An²`, `Bt₂² > 0` from the energy jump) is returned
instead. The switch-on tangential field direction is arbitrary; the sign
convention is `Bt₂ ≥ 0`.
"""
function rankine_hugoniot(up::MHDState{T}, γ::Real; μ0::Real = 1.0) where {T}
    γT = _require_valid_gamma(γ, T)
    μ = _require_valid_rh_inputs(up, μ0)
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
        elseif isfinite(Rprev) && isfinite(Rc) && sign(Rc) != sign(Rprev)
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
    # Bt₁ = 0 bifurcation: the X-parametrization has B2t ∝ Bt, so the switch-on
    # family is invisible to the scan above. Construct it explicitly and prefer
    # the evolutionary branch — inside the switch-on window the gasdynamic root
    # crosses the Alfvén point (a non-evolutionary 1→4 jump) and the switch-on
    # branch is the physical solution.
    so = _switch_on_downstream(up, γT, μ)
    if so !== nothing && (!found || !_alfven_ordered(up, down, μ))
        down, Xsol = so
        found = true
    end
    residuals = _rh_residuals(up, down, γT, μ)
    return (; down, X = found ? Xsol : one(T), residuals)
end

"""
    rh_branches(up, γ; μ0=1.0, nscan=2000) -> Vector{NamedTuple}

Branch tracking: find ALL compressive (X>1) Rankine–Hugoniot solutions by
scanning the energy-jump residual over (1, X_max] and bisecting every sign
change. Each returned `(; X, down, residuals)` satisfies all six jump conditions
(residuals ≈ 0). A perpendicular shock yields a single fast branch; oblique
upstream states can admit slow/intermediate/fast branches. For a field-aligned
upstream (`Bt = 0`) the scan's parametrization has `B2t ∝ Bt`, so the switch-on
branch (`X = M_An²`, `Bt₂ ≥ 0` fixed by the energy jump) is constructed
explicitly and included whenever the switch-on window is non-empty. Branches
are returned sorted by increasing `X`.
"""
function rh_branches(up::MHDState{T}, γ::Real; μ0::Real = 1.0, nscan::Int = 2000) where {T}
    γT = _require_valid_gamma(γ, T)
    μ = _require_valid_rh_inputs(up, μ0)
    nscan >= 1 || throw(ArgumentError("nscan must be positive"))
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
        if Rprev != 0 && isfinite(Rprev) && isfinite(Rc) && sign(Rc) != sign(Rprev)
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
            push!(branches, (; X = Xs, down, residuals = _rh_residuals(up, down, γT, μ)))
        end
        Rprev = Rc
        Xprev = Xc
    end
    # Bt₁ = 0 bifurcation: the parametrized scan has B2t ∝ Bt and cannot reach
    # the switch-on family — construct that branch explicitly (see
    # _switch_on_downstream) so the returned set really is ALL compressive roots.
    so = _switch_on_downstream(up, γT, μ)
    if so !== nothing
        sdown, sX = so
        push!(branches, (; X = sX, down = sdown, residuals = _rh_residuals(up, sdown, γT, μ)))
        sort!(branches; by = b -> b.X)
    end
    return branches
end
