# shock_sim.jl — 1D perpendicular collisionless hybrid shock.
#
# Non-periodic shock normal x via the SBP-(2,1) derivative; reflecting wall at
# x=0 (the piston, in the downstream frame); held inflow at x=Lx. Perpendicular
# geometry (B0 = B0 ẑ ⟂ normal) collapses the field to a single dynamic Bz:
#
#   J = (0, −∂x Bz, 0)
#   E_x = −u_y Bz − Bz(∂x Bz)/n − ∂x p_e/n         (cross-shock field)
#   E_y =  u_x Bz − η ∂x Bz                          (motional + resistive)
#   ∂t Bz = −∂x(u_x Bz) + η ∂²Bz   + inflow SAT toward B0
#
# Flux-conservative ⇒ Bz/n is advected ⇒ across a steady shock Bz2/Bz1 = n2/n1
# (cross-checked vs `rankine_hugoniot`). Formulation + BCs cross-checked by an
# independent design review; stability (SBP central + inflow SAT + small η) was
# validated by a field-only static-pulse pre-test before adding particles.

mutable struct PerpShock{T}
    s::SBP1D{T}
    x::Vector{T}                 # node coordinates
    Bz::Vector{T}                # the one dynamic field
    n::Vector{T}
    ux::Vector{T}
    uy::Vector{T}
    pe::Vector{T}
    Ex::Vector{T}
    Ey::Vector{T}
    Te::T
    γe::T
    η::T
    τ::T
    B0::T
    nfloor::T
    # workspaces
    wsum::Vector{T}
    mx::Vector{T}
    my::Vector{T}
    DBz::Vector{T}
    Dpe::Vector{T}
    Fb::Vector{T}
    DF::Vector{T}
    d2::Vector{T}
    k1::Vector{T}
    k2::Vector{T}
    k3::Vector{T}
    k4::Vector{T}
    tmp::Vector{T}
    closure::Symbol              # electron closure: :polytropic (pe=Te·n^γe) or :energy (Leroy eq 6)
end

"""
    PerpShock(N, Lx; Te, γe, η, τ, B0, nfloor)

Allocate a 1D perpendicular hybrid-shock field state on `N` SBP nodes over
`[0, Lx]`. `Te`/`γe` set the electron closure `p_e = Te·n^{γe}` (n0=1); `η` is
the resistive/numerical dissipation; `τ` the inflow SAT strength (≈ inflow
speed); `B0` the upstream field.
"""
function PerpShock(
    N::Integer,
    Lx::T;
    Te = 0.5,
    γe = 5 / 3,
    η = 0.01,
    τ = 3.0,
    B0 = 1.0,
    nfloor = 1e-6,
    closure::Symbol = :polytropic,
) where {T<:AbstractFloat}
    LxT = _require_finite_positive_real("Lx", Lx, T)
    TeT = _require_finite_nonnegative_real("Te", Te, T)
    γeT = _require_valid_gamma(γe, T)
    ηT = _require_finite_nonnegative_real("η", η, T)
    τT = _require_finite_real("τ", τ, T)
    B0T = _require_finite_real("B0", B0, T)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, T)
    (closure === :polytropic || closure === :energy) ||
        throw(ArgumentError("closure must be :polytropic or :energy, got :$closure"))
    s = SBP1D(N, LxT)
    x = collect(range(zero(T), LxT; length = N))
    z() = zeros(T, N)
    PerpShock{T}(
        s,
        x,
        fill(B0T, N),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        TeT,
        γeT,
        ηT,
        τT,
        B0T,
        nfloorT,
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        z(),
        closure,
    )
end

# Upstream particle injection at the inflow face (x=Lx). Without it the finite
# initial reservoir drifts to the wall and the upstream EMPTIES within a few
# Ω_ci⁻¹, so the shock cannot reach a sustained/quasi-stationary state (the field
# eventually collapses). Leroy et al. 1982 and Hellinger et al. 2002 both maintain
# a constant upstream flux n₁V₁; this reproduces that. Off by default so the
# validated transient SHK-002 behaviour is unchanged.
"""
    ShockInjector(rng; n0=1.0, drift, vthi, ut=(0,0), σt=vthi, weight)

Sustained flux-weighted upstream injection at the inflow face for [`step_shock!`].
`drift` is the inward (toward the wall) bulk speed (= the upstream drift `U0`),
`vthi` the normal thermal speed, `weight` the per-particle weight (use
[`shock_density_weight`](@ref) so the injected flux carries density `n0`). Carries
the flux accumulator and next-id counter across steps.
"""
mutable struct ShockInjector{T,R}
    rng::R
    n0::T
    drift::T
    vthi::T
    ut::NTuple{2,T}
    σt::T
    weight::T
    acc::Base.RefValue{Float64}
    nextid::Base.RefValue{UInt64}
end

function ShockInjector(
    rng::R;
    n0::Real = 1.0,
    drift::Real,
    vthi::Real,
    ut::NTuple{2,<:Real} = (0.0, 0.0),
    σt::Real = vthi,
    weight::Real,
    first_id::Integer = 1,
) where {R}
    T = float(promote_type(typeof(n0), typeof(drift), typeof(vthi), typeof(weight)))
    n0T = _require_finite_nonnegative_real("n0", n0, T)
    driftT = _require_finite_real("drift", drift, T)
    vthiT = _require_finite_nonnegative_real("vthi", vthi, T)
    σtT = _require_finite_nonnegative_real("σt", σt, T)
    wT = _require_finite_positive_real("weight", weight, T)
    first_id >= 1 || throw(ArgumentError("first_id must be ≥ 1, got $first_id"))
    return ShockInjector{T,R}(
        rng,
        n0T,
        driftT,
        vthiT,
        (T(ut[1]), T(ut[2])),
        σtT,
        wT,
        Ref(0.0),
        Ref(UInt64(first_id)),
    )
end

# Inject one step's worth of upstream plasma at the inflow face x=Lx, drifting
# inward (−x, toward the wall). Returns the number injected.
function _inject_upstream!(
    sh::PerpShock{T},
    ps::ParticleSet{1,T},
    inj::ShockInjector,
    dt::T,
) where {T}
    Lx = sh.x[end]
    return inject_face_1d!(
        ps,
        inj.rng,
        Lx,
        -1,                       # inward = −x (toward the wall)
        inj.n0,
        inj.drift,
        inj.vthi,
        inj.ut,
        inj.σt,
        dt,
        inj.weight,
        inj.acc,
        inj.nextid,
    )
end

# CIC stencil on the SBP node grid (node i at (i−1)dx), 1-based, FOLD/clamp at
# the non-periodic boundaries (no wrap). Conserves total deposited weight.
@inline function _cic_sbp(xp::T, dx::T, N::Int) where {T}
    sidx = xp / dx
    i0 = floor(Int, sidx)
    f = sidx - i0
    a = i0 + 1
    b = i0 + 2
    a = a < 1 ? 1 : (a > N ? N : a)
    b = b < 1 ? 1 : (b > N ? N : b)
    return a, b, one(T) - f, f
end

"Deposit number density n=Σw S/H and bulk velocity u=(Σw v S)/(Σw S) onto SBP nodes."
function deposit_moments!(sh::PerpShock{T}, ps::ParticleSet{1,T}) where {T}
    N = sh.s.n
    dx = sh.s.dx
    fill!(sh.wsum, zero(T))
    fill!(sh.mx, zero(T))
    fill!(sh.my, zero(T))
    xp = ps.x[1]
    vx = ps.v[1]
    vy = ps.v[2]
    w = ps.weight
    @inbounds for p in eachindex(w)
        a, b, wa, wb = _cic_sbp(xp[p], dx, N)
        ww = w[p]
        sh.wsum[a] += wa * ww
        sh.wsum[b] += wb * ww
        sh.mx[a] += wa * ww * vx[p]
        sh.mx[b] += wb * ww * vx[p]
        sh.my[a] += wa * ww * vy[p]
        sh.my[b] += wb * ww * vy[p]
    end
    nf = sh.nfloor
    @inbounds for i = 1:N
        vol = sh.s.H[i]
        ws = sh.wsum[i]
        sh.n[i] = ws / vol
        wsf = max(ws, nf * vol)
        sh.ux[i] = sh.mx[i] / wsf
        sh.uy[i] = sh.my[i] / wsf
    end
    return sh
end

"Electric field from the 1D perpendicular Ohm's law."
function compute_E!(sh::PerpShock{T}) where {T}
    N = sh.s.n
    η = sh.η
    nf = sh.nfloor
    # :polytropic — pe is an algebraic function of n; :energy — pe is a dynamical
    # field evolved by the energy equation (do not overwrite it here).
    if sh.closure === :polytropic
        @. sh.pe = sh.Te * sh.n^sh.γe
    end
    sbp_deriv!(sh.DBz, sh.Bz, sh.s)
    sbp_deriv!(sh.Dpe, sh.pe, sh.s)
    @inbounds for i = 1:N
        ninv = one(T) / max(sh.n[i], nf)
        sh.Ex[i] = -sh.uy[i] * sh.Bz[i] - sh.Bz[i] * sh.DBz[i] * ninv - sh.Dpe[i] * ninv
        sh.Ey[i] = sh.ux[i] * sh.Bz[i] - η * sh.DBz[i]
    end
    return sh
end

# narrow 2nd difference (damps the 2Δx mode that the wide D∘D would miss)
function _d2narrow!(d2::Vector{T}, f::Vector{T}, dx::T) where {T}
    n = length(f)
    fill!(d2, zero(T))
    @inbounds for i = 2:n-1
        d2[i] = (f[i+1] - 2f[i] + f[i-1]) / dx^2
    end
    return d2
end

# dBz = −∂x(ux Bz) + η ∂²Bz + SAT(inflow node N → B0)
function _bz_rhs!(dB::Vector{T}, Bz::Vector{T}, sh::PerpShock{T}) where {T}
    N = sh.s.n
    @. sh.Fb = sh.ux * Bz
    sbp_deriv!(sh.DF, sh.Fb, sh.s)
    _d2narrow!(sh.d2, Bz, sh.s.dx)
    @. dB = -sh.DF + sh.η * sh.d2
    dB[N] += -(sh.τ / sh.s.H[N]) * (Bz[N] - sh.B0)
    return dB
end

function _rk4_bz!(sh::PerpShock{T}, h::T) where {T}
    Bz = sh.Bz
    _bz_rhs!(sh.k1, Bz, sh)
    @. sh.tmp = Bz + h / 2 * sh.k1
    _bz_rhs!(sh.k2, sh.tmp, sh)
    @. sh.tmp = Bz + h / 2 * sh.k2
    _bz_rhs!(sh.k3, sh.tmp, sh)
    @. sh.tmp = Bz + h * sh.k3
    _bz_rhs!(sh.k4, sh.tmp, sh)
    @. Bz += h / 6 * (sh.k1 + 2sh.k2 + 2sh.k3 + sh.k4)
    return sh
end

@inline function _gather_sbp(field::Vector{T}, xp::T, dx::T, N::Int) where {T}
    a, b, wa, wb = _cic_sbp(xp, dx, N)
    return wa * field[a] + wb * field[b]
end

"Initialize carried E from the loaded particles (call once after loading)."
function init_shock!(sh::PerpShock{T}, ps::ParticleSet{1,T}) where {T}
    deposit_moments!(sh, ps)
    @. sh.pe = sh.Te * sh.n^sh.γe        # consistent IC (the :energy closure evolves from here)
    compute_E!(sh)
    return sh
end

"""
    step_shock!(sh, ps, dt; NB=2)

Advance the perpendicular shock one timestep: Boris push (E^n,B^n) → reflect at
the wall (x=0, lo-only) / absorb past the inflow (x=Lx) → deposit n,u_x →
subcycle Bz (SBP + SAT + η) → recompute E.
"""
function step_shock!(
    sh::PerpShock{T},
    ps::ParticleSet{1,T},
    dt::Real;
    NB::Integer = 2,
    injector::Union{Nothing,ShockInjector} = nothing,
) where {T}
    N = sh.s.n
    dx = sh.s.dx
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step_shock!")
    Lx = sh.x[end]
    xp = ps.x[1]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    # 1. push with carried E^n, B^n
    @inbounds for p in eachindex(ps.weight)
        Exp = _gather_sbp(sh.Ex, xp[p], dx, N)
        Eyp = _gather_sbp(sh.Ey, xp[p], dx, N)
        Bzp = _gather_sbp(sh.Bz, xp[p], dx, N)
        nx, ny, nz =
            boris_kick(vx[p], vy[p], vz[p], Exp, Eyp, zero(T), zero(T), zero(T), Bzp, one(T), dtT)
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        xp[p] += dtT * nx
    end
    # 2. reflecting wall at x=0 (lo only — DO NOT reflect the inflow at x=Lx)
    @inbounds for p in eachindex(ps.weight)
        if xp[p] < 0
            xp[p] = -xp[p]
            vx[p] = -vx[p]
        end
    end
    # absorb particles past the inflow boundary
    write = 0
    @inbounds for p in eachindex(ps.weight)
        if xp[p] <= Lx
            write += 1
            if write != p
                xp[write] = xp[p]
                for c = 1:3
                    ps.v[c][write] = ps.v[c][p]
                end
                ps.weight[write] = ps.weight[p]
                ps.id[write] = ps.id[p]
                ps.tag[write] = ps.tag[p]
            end
        end
    end
    if write < length(ps.weight)
        resize!(ps.x[1], write)
        for c = 1:3
            resize!(ps.v[c], write)
        end
        resize!(ps.weight, write)
        resize!(ps.id, write)
        resize!(ps.tag, write)
    end
    # 2b. sustained upstream injection at the inflow face (opt-in). Appends fresh
    #     upstream plasma so the reservoir does not deplete over the run; deposited
    #     in this step's moment pass below.
    injector === nothing || _inject_upstream!(sh, ps, injector, dtT)
    # 3. deposit moments
    deposit_moments!(sh, ps)
    # 4. subcycle Bz
    hb = dtT / NB
    for _ = 1:NB
        _rk4_bz!(sh, hb)
    end
    # 5. recompute E for the next push
    compute_E!(sh)
    return sh
end

# ============================================================================
# §11.3 Rankine–Hugoniot two-state shock — the wall-less, shock-REST-FRAME model
# Leroy et al. 1982 uses (vs the reflecting-wall §11.2 model above): no wall;
# upstream plasma flows IN at x=Lx and downstream plasma flows OUT at x=0, with
# the shock held roughly stationary by a two-ended flux balance, and BOTH ends
# carry a mean-field B BC (B0 upstream at x=Lx, B2 downstream at x=0). The
# downstream boundary is a THERMAL reservoir (exiting ions are reinserted with a
# downstream-thermal velocity), NOT a specular wall — which is what lets the
# self-consistent ion reflection / foot develop (Leroy's α).
# ============================================================================

# dBz with SAT mean-field BCs at BOTH ends: node N → B0 (upstream), node 1 → B_down.
function _bz_rhs_leroy!(
    dB::Vector{T},
    Bz::Vector{T},
    sh::PerpShock{T},
    B_down::T,
    τ_down::T,
) where {T}
    N = sh.s.n
    @. sh.Fb = sh.ux * Bz
    sbp_deriv!(sh.DF, sh.Fb, sh.s)
    _d2narrow!(sh.d2, Bz, sh.s.dx)
    @. dB = -sh.DF + sh.η * sh.d2
    dB[N] += -(sh.τ / sh.s.H[N]) * (Bz[N] - sh.B0)         # upstream → B0
    dB[1] += -(τ_down / sh.s.H[1]) * (Bz[1] - B_down)      # downstream → B_down
    return dB
end

function _rk4_bz_leroy!(sh::PerpShock{T}, h::T, B_down::T, τ_down::T) where {T}
    Bz = sh.Bz
    _bz_rhs_leroy!(sh.k1, Bz, sh, B_down, τ_down)
    @. sh.tmp = Bz + h / 2 * sh.k1
    _bz_rhs_leroy!(sh.k2, sh.tmp, sh, B_down, τ_down)
    @. sh.tmp = Bz + h / 2 * sh.k2
    _bz_rhs_leroy!(sh.k3, sh.tmp, sh, B_down, τ_down)
    @. sh.tmp = Bz + h * sh.k3
    _bz_rhs_leroy!(sh.k4, sh.tmp, sh, B_down, τ_down)
    @. Bz += h / 6 * (sh.k1 + 2sh.k2 + 2sh.k3 + sh.k4)
    return sh
end

# Leroy 1982 eq 6 — electron energy equation (advective form, our normalisation):
#   ∂t pe = −ux ∂x pe − γe pe ∂x ux + (γe−1) η (∂x Bz)²        [Jy = −∂x Bz]
# Advanced once per step (after the Bz subcycle), matching Leroy's scheme order
# (move ions → moments → field → pe → E). ux, ∂x ux and the Ohmic source are
# frozen across the RK4 stages (only ∂x pe varies); ∂x ux is held in sh.d2 and the
# Ohmic source in sh.Fb. The two boundary nodes carry fresh injected/reservoir
# plasma (not yet Ohmically heated) so are clamped to the polytropic value, and pe
# is floored positive. SBP central advection (no extra dissipation) — the resistive
# ramp it rides on is itself η-smoothed, which keeps the pe profile well-behaved.
function _pe_energy_rhs!(dpe::Vector{T}, pe::Vector{T}, sh::PerpShock{T}, γe::T) where {T}
    sbp_deriv!(sh.DF, pe, sh.s)                                   # ∂x pe (per stage)
    _d2narrow!(sh.Dpe, pe, sh.s.dx)                              # ∂² pe (electron heat conduction)
    # advection − compression + Ohmic heating + resistive-scale heat conduction.
    # The η∂²pe term is electron heat diffusion at the collisional (resistive) scale;
    # it also damps the 2Δx advection mode that bare SBP central differencing leaves
    # undamped (without it pe develops grid-scale spikes/undershoots). Leroy's eq 6
    # omits conduction, so this is the one numerically-required addition — kept at the
    # same η as the field so it is the minimal resistive-scale regularisation.
    @. dpe = -sh.ux * sh.DF - γe * pe * sh.d2 + sh.Fb + sh.η * sh.Dpe
    return dpe
end

function _rk4_pe_energy!(sh::PerpShock{T}, h::T) where {T}
    N = sh.s.n
    γe = sh.γe
    sbp_deriv!(sh.DBz, sh.Bz, sh.s)                               # ∂x Bz  (= −Jy)
    sbp_deriv!(sh.d2, sh.ux, sh.s)                                # ∂x ux  (frozen this step)
    @. sh.Fb = (γe - one(T)) * sh.η * sh.DBz^2                    # Ohmic source (frozen, ≥ 0)
    pe = sh.pe
    _pe_energy_rhs!(sh.k1, pe, sh, γe)
    @. sh.tmp = pe + h / 2 * sh.k1
    _pe_energy_rhs!(sh.k2, sh.tmp, sh, γe)
    @. sh.tmp = pe + h / 2 * sh.k2
    _pe_energy_rhs!(sh.k3, sh.tmp, sh, γe)
    @. sh.tmp = pe + h * sh.k3
    _pe_energy_rhs!(sh.k4, sh.tmp, sh, γe)
    @. pe += h / 6 * (sh.k1 + 2sh.k2 + 2sh.k3 + sh.k4)
    pf = sh.nfloor * sh.Te                                        # tiny positive floor
    @inbounds for i = 1:N
        pe[i] < pf && (pe[i] = pf)
    end
    pe[1] = sh.Te * sh.n[1]^γe                                    # downstream-reservoir BC
    pe[N] = sh.Te * sh.n[N]^γe                                    # upstream-inflow BC
    return sh
end

"""
    LeroyBoundary(rng; V1, vthi, V2, vth2, B_down, p_up, τ_down=V2)

Two-ended flux boundary for [`step_leroy_shock!`]. Ions are conserved: each ion
leaving the domain is reinserted — with probability `p_up` as fresh upstream
inflow at x=Lx (drift −V1, thermal `vthi`), otherwise as a downstream-reservoir
ion at x=0 (into the box, thermal `vth2`). `p_up = V2/flux_per_density(V2, vth2)`
makes the net upstream inflow equal the downstream outflow n₁V₁ = n₂V₂ in steady
state. `B_down` is the downstream field BC; `τ_down` the downstream field-SAT
strength.
"""
mutable struct LeroyBoundary{T,R}
    V1::T
    vthi::T
    V2::T
    vth2::T
    B_down::T
    τ_down::T
    p_up::T
    rng::R
end
function LeroyBoundary(
    rng::R;
    V1::Real,
    vthi::Real,
    V2::Real,
    vth2::Real,
    B_down::Real,
    p_up::Real,
    τ_down::Real = V2,
) where {R}
    T = float(promote_type(typeof(V1), typeof(V2), typeof(B_down)))
    (0 <= p_up <= 1) || throw(ArgumentError("p_up must be in [0,1], got $p_up"))
    return LeroyBoundary{T,R}(T(V1), T(vthi), T(V2), T(vth2), T(B_down), T(τ_down), T(p_up), rng)
end

"""
    step_leroy_shock!(sh, ps, dt; NB=2, bc::LeroyBoundary)

One step of the §11.3 wall-less shock-rest-frame perpendicular shock: Boris push →
two-ended particle recycle (no wall) → deposit → two-ended-BC Bz subcycle →
recompute E. Particle count is exactly conserved (recycle is in place).
"""
function step_leroy_shock!(
    sh::PerpShock{T},
    ps::ParticleSet{1,T},
    dt::Real;
    NB::Integer = 2,
    bc::LeroyBoundary,
) where {T}
    N = sh.s.n
    dx = sh.s.dx
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step_leroy_shock!")
    Lx = sh.x[end]
    xp = ps.x[1]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    rng = bc.rng
    # 1. push with carried E^n, B^n
    @inbounds for p in eachindex(ps.weight)
        Exp = _gather_sbp(sh.Ex, xp[p], dx, N)
        Eyp = _gather_sbp(sh.Ey, xp[p], dx, N)
        Bzp = _gather_sbp(sh.Bz, xp[p], dx, N)
        nx, ny, nz =
            boris_kick(vx[p], vy[p], vz[p], Exp, Eyp, zero(T), zero(T), zero(T), Bzp, one(T), dtT)
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        xp[p] += dtT * nx
    end
    # 2. two-ended recycle (NO wall): each exiting ion is reinserted in place, as
    #    fresh upstream inflow (prob p_up) or a downstream-reservoir ion (else).
    ε = Lx * T(1e-6)
    @inbounds for p in eachindex(ps.weight)
        if xp[p] < zero(T) || xp[p] > Lx
            if rand(rng, T) < bc.p_up                    # → upstream inflow at x=Lx
                xp[p] = Lx - ε
                vx[p] = -bc.V1 + bc.vthi * randn(rng, T)
                vy[p] = bc.vthi * randn(rng, T)
                vz[p] = bc.vthi * randn(rng, T)
            else                                         # → downstream reservoir at x=0 (into box)
                xp[p] = ε
                vx[p] = abs(bc.vth2 * randn(rng, T))
                vy[p] = bc.vth2 * randn(rng, T)
                vz[p] = bc.vth2 * randn(rng, T)
            end
        end
    end
    # 3. deposit
    deposit_moments!(sh, ps)
    # 4. subcycle Bz (whistler CFL) — and pe alongside it when the energy closure is
    #    on, so the Ohmic source rides the evolving Bz and both share the safe substep
    hb = dtT / NB
    energy = sh.closure === :energy
    for _ = 1:NB
        _rk4_bz_leroy!(sh, hb, bc.B_down, bc.τ_down)
        energy && _rk4_pe_energy!(sh, hb)
    end
    # 5. recompute E
    compute_E!(sh)
    return sh
end

"""
    shock_density_weight(n0, Lx, N)

Particle weight so a uniform load over the SBP domain `[0, Lx]` deposits density
`n0` (∫n dx = Σw exactly with the SBP norm as node volume): `w = n0·Lx/N`.
"""
function shock_density_weight(n0, Lx, N)
    isfinite(n0) && n0 >= 0 || throw(ArgumentError("n0 must be finite and non-negative"))
    isfinite(Lx) && Lx > 0 || throw(ArgumentError("Lx must be finite and positive"))
    Np = _require_positive_intlike("N", N)
    return n0 * Lx / Np
end
