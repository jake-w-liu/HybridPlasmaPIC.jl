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
) where {T<:AbstractFloat}
    LxT = _require_finite_positive_real("Lx", Lx, T)
    TeT = _require_finite_nonnegative_real("Te", Te, T)
    γeT = _require_valid_gamma(γe, T)
    ηT = _require_finite_nonnegative_real("η", η, T)
    τT = _require_finite_real("τ", τ, T)
    B0T = _require_finite_real("B0", B0, T)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, T)
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
        ws = sh.wsum[i]
        sh.n[i] = ws / sh.s.H[i]
        wsf = max(ws, nf)
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
    @. sh.pe = sh.Te * sh.n^sh.γe
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
    compute_E!(sh)
    return sh
end

"""
    step_shock!(sh, ps, dt; NB=2)

Advance the perpendicular shock one timestep: Boris push (E^n,B^n) → reflect at
the wall (x=0, lo-only) / absorb past the inflow (x=Lx) → deposit n,u_x →
subcycle Bz (SBP + SAT + η) → recompute E.
"""
function step_shock!(sh::PerpShock{T}, ps::ParticleSet{1,T}, dt::Real; NB::Integer = 2) where {T}
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

"""
    shock_density_weight(n0, Lx, N)

Particle weight so a uniform load over the SBP domain `[0, Lx]` deposits density
`n0` (∫n dx = Σw exactly with the SBP norm as node volume): `w = n0·Lx/N`.
"""
shock_density_weight(n0, Lx, N) = n0 * Lx / N
