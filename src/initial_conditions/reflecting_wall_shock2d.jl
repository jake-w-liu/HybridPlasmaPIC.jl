# shock2d.jl — 2D perpendicular collisionless hybrid shock (Phase 11).
#
# Shock normal x (non-periodic, SBP-(2,1) derivative + reflecting wall + SAT
# inflow); transverse y (periodic, Fourier). B0 = B0 ẑ is OUT of the (x,y)
# simulation plane, so ∇·B = ∂x Bx + ∂y By = 0 holds trivially (Bx=By=0) and the
# ONLY dynamic field is Bz(x,y):
#   J  = ∇×B = (∂y Bz, −∂x Bz, 0)
#   E_x = −u_y Bz + (J×B)_x/n − ∂x p_e/n + η J_x = −u_y Bz + J_y Bz/n − ∂x p_e/n + η J_x
#   E_y =  u_x Bz + (J×B)_y/n − ∂y p_e/n + η J_y =  u_x Bz − J_x Bz/n − ∂y p_e/n + η J_y
#   ∂t Bz = −(∇×E)_z = −(∂x E_y − ∂y E_x) + ν ∂xx Bz   (+ inflow SAT toward B0)
#
# Reduces to the 1D PerpShock when y-uniform; supports shock-front rippling x_s(y).
# The state carries one reusable transverse Fourier workspace for all y
# derivatives so field updates do not allocate FFT buffers every call.

mutable struct PerpShock2D{T,YW}
    sbp::SBP1D{T}
    nx::Int
    ny::Int
    Lx::T
    Ly::T
    ywork::YW
    dy::T
    x::Vector{T}
    y::Vector{T}
    Bz::Matrix{T}
    n::Matrix{T}
    ux::Matrix{T}
    uy::Matrix{T}
    uz::Matrix{T}
    pe::Matrix{T}
    Ex::Matrix{T}
    Ey::Matrix{T}
    Te::T
    γe::T
    η::T
    τ::T
    B0::T
    nfloor::T
    Jx::Matrix{T}
    Jy::Matrix{T}
    dpe_x::Matrix{T}
    dpe_y::Matrix{T}
    wsum::Matrix{T}
    mx::Matrix{T}
    my::Matrix{T}
    mz::Matrix{T}
    dEy_x::Matrix{T}
    dEx_y::Matrix{T}
    d2::Matrix{T}
    k1::Matrix{T}
    k2::Matrix{T}
    k3::Matrix{T}
    k4::Matrix{T}
    Btmp::Matrix{T}
end

"""
    PerpShock2D(nx, ny, Lx, Ly; Te, γe, η, τ, B0, nfloor)

2D perpendicular hybrid-shock field state: `nx` SBP nodes over `[0,Lx]` (shock
normal) × `ny` periodic Fourier nodes over `[0,Ly]` (transverse).
"""
function PerpShock2D(
    nx::Integer,
    ny::Integer,
    Lx::T,
    Ly::T;
    Te = 0.125,
    γe = 5 / 3,
    η = 0.02,
    τ = 3.0,
    B0 = 1.0,
    nfloor = 1e-6,
) where {T<:AbstractFloat}
    nxi = Int(nx)
    nyi = Int(ny)
    LxT = _require_finite_positive_real("Lx", Lx, T)
    LyT = _require_finite_positive_real("Ly", Ly, T)
    TeT = _require_finite_nonnegative_real("Te", Te, T)
    γeT = _require_valid_gamma(γe, T)
    ηT = _require_finite_nonnegative_real("η", η, T)
    τT = _require_finite_real("τ", τ, T)
    B0T = _require_finite_real("B0", B0, T)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, T)
    sbp = SBP1D(nxi, LxT)
    ywork = FourierDerivYWorkspace(nxi, nyi, LyT)
    x = collect(range(zero(T), LxT; length = nxi))
    y = [(j - 1) * LyT / nyi for j = 1:nyi]
    M() = zeros(T, nxi, nyi)
    PerpShock2D{T,typeof(ywork)}(
        sbp,
        nxi,
        nyi,
        LxT,
        LyT,
        ywork,
        LyT / nyi,
        x,
        y,
        fill(B0T, nxi, nyi),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        TeT,
        γeT,
        ηT,
        τT,
        B0T,
        nfloorT,
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
        M(),
    )
end

# CIC stencil along the SBP x-grid (clamped / FOLD) — same as the 1D shock.
@inline function _cicx(xp::T, dx::T, nx::Int) where {T}
    s = xp / dx
    i0 = floor(Int, s)
    f = s - i0
    a = i0 + 1
    b = i0 + 2
    a = a < 1 ? 1 : (a > nx ? nx : a)
    b = b < 1 ? 1 : (b > nx ? nx : b)
    return a, b, one(T) - f, f
end
# periodic CIC stencil along y
@inline function _cicy(yp::T, dy::T, ny::Int) where {T}
    s = yp / dy
    j0 = floor(Int, s)
    f = s - j0
    a = mod(j0, ny) + 1
    b = mod(j0 + 1, ny) + 1
    return a, b, one(T) - f, f
end

"Deposit n and bulk velocity (u_x,u_y,u_z) onto the (SBP-x × periodic-y) mesh."
function deposit_moments2d!(sh::PerpShock2D{T}, ps::ParticleSet{2,T}) where {T}
    fill!(sh.wsum, zero(T))
    fill!(sh.mx, zero(T))
    fill!(sh.my, zero(T))
    fill!(sh.mz, zero(T))
    dx = sh.sbp.dx
    xp = ps.x[1]
    yp = ps.x[2]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    w = ps.weight
    @inbounds for p in eachindex(w)
        ax, bx, wxa, wxb = _cicx(xp[p], dx, sh.nx)
        ay, by, wya, wyb = _cicy(yp[p], sh.dy, sh.ny)
        ww = w[p]
        for (ix, wx) in ((ax, wxa), (bx, wxb)), (iy, wy) in ((ay, wya), (by, wyb))
            g = wx * wy * ww
            sh.wsum[ix, iy] += g
            sh.mx[ix, iy] += g * vx[p]
            sh.my[ix, iy] += g * vy[p]
            sh.mz[ix, iy] += g * vz[p]
        end
    end
    nf = sh.nfloor
    H = sh.sbp.H
    dy = sh.dy
    @inbounds for j = 1:sh.ny, i = 1:sh.nx
        vol = H[i] * dy
        ws = sh.wsum[i, j]
        sh.n[i, j] = ws / vol
        wsf = max(ws, nf * vol)
        sh.ux[i, j] = sh.mx[i, j] / wsf
        sh.uy[i, j] = sh.my[i, j] / wsf
        sh.uz[i, j] = sh.mz[i, j] / wsf
    end
    return sh
end

@inline function _gather2d(F::Matrix{T}, xp::T, yp::T, dx::T, dy::T, nx::Int, ny::Int) where {T}
    ax, bx, wxa, wxb = _cicx(xp, dx, nx)
    ay, by, wya, wyb = _cicy(yp, dy, ny)
    return wxa * (wya * F[ax, ay] + wyb * F[ax, by]) + wxb * (wya * F[bx, ay] + wyb * F[bx, by])
end

"Electric field from the 2D perpendicular Ohm's law (uses current Bz, n, u)."
function compute_E2d!(sh::PerpShock2D{T}) where {T}
    fourier_deriv_y!(sh.Jx, sh.Bz, sh.ywork)         # J_x = ∂y Bz
    sbp_deriv_x!(sh.Jy, sh.Bz, sh.sbp)
    sh.Jy .*= -one(T)   # J_y = −∂x Bz
    @. sh.pe = sh.Te * sh.n^sh.γe
    sbp_deriv_x!(sh.dpe_x, sh.pe, sh.sbp)
    fourier_deriv_y!(sh.dpe_y, sh.pe, sh.ywork)
    η = sh.η
    nf = sh.nfloor
    @inbounds for k in eachindex(sh.Bz)
        bz = sh.Bz[k]
        inv = one(T) / max(sh.n[k], nf)
        sh.Ex[k] = -sh.uy[k] * bz + sh.Jy[k] * bz * inv - sh.dpe_x[k] * inv + η * sh.Jx[k]
        sh.Ey[k] = sh.ux[k] * bz - sh.Jx[k] * bz * inv - sh.dpe_y[k] * inv + η * sh.Jy[k]
    end
    return sh
end

# narrow ∂xx for x-dissipation (zero at x-boundaries)
function _d2x!(d2::Matrix{T}, f::Matrix{T}, dx::T) where {T}
    nx = size(f, 1)
    fill!(d2, zero(T))
    @inbounds for j in axes(f, 2), i = 2:nx-1
        d2[i, j] = (f[i+1, j] - 2f[i, j] + f[i-1, j]) / dx^2
    end
    return d2
end

# Magnetic flux is frozen to the ion bulk flow (motional induction): with the
# frozen n+1/2 ion velocity (u_x,u_y), Faraday of the motional field −u×B plus
# resistive Ohm term ηJ gives
#   ∂t Bz = −∂x(u_x Bz) − ∂y(u_y Bz) + η(∂xx + ∂yy)Bz   (+ inflow SAT toward B0).
# This is the 2D analog of the verified 1D flux form (stable advection — no stiff
# Hall whistler in the field advance; the Hall/pressure terms enter only the
# particle electric field via compute_E2d!). Reduces to the 1D flux-plus-resistivity form when y-uniform.
function _bz_rhs2d!(dB::Matrix{T}, Bz_trial::Matrix{T}, sh::PerpShock2D{T}) where {T}
    @. sh.Jx = sh.ux * Bz_trial                       # F_x = u_x Bz (scratch)
    @. sh.Jy = sh.uy * Bz_trial                       # F_y = u_y Bz (scratch)
    sbp_deriv_x!(sh.dEy_x, sh.Jx, sh.sbp)             # ∂x F_x
    fourier_deriv_y!(sh.dEx_y, sh.Jy, sh.ywork)       # ∂y F_y
    _d2x!(sh.d2, Bz_trial, sh.sbp.dx)
    fourier_deriv_y!(sh.Jx, Bz_trial, sh.ywork)
    fourier_deriv_y!(sh.Jy, sh.Jx, sh.ywork)
    @. dB = -(sh.dEy_x + sh.dEx_y) + sh.η * (sh.d2 + sh.Jy)
    τH = sh.τ / sh.sbp.H[sh.nx]
    @inbounds for j = 1:sh.ny
        dB[sh.nx, j] += -τH * (Bz_trial[sh.nx, j] - sh.B0)
    end
    return dB
end

function _rk4_bz2d!(sh::PerpShock2D{T}, h::T) where {T}
    B = sh.Bz
    _bz_rhs2d!(sh.k1, B, sh)
    @. sh.Btmp = B + h / 2 * sh.k1
    _bz_rhs2d!(sh.k2, sh.Btmp, sh)
    @. sh.Btmp = B + h / 2 * sh.k2
    _bz_rhs2d!(sh.k3, sh.Btmp, sh)
    @. sh.Btmp = B + h * sh.k3
    _bz_rhs2d!(sh.k4, sh.Btmp, sh)
    @. B += h / 6 * (sh.k1 + 2sh.k2 + 2sh.k3 + sh.k4)
    return sh
end

"Initialize the carried E from the loaded particles."
function init_shock2d!(sh::PerpShock2D{T}, ps::ParticleSet{2,T}) where {T}
    deposit_moments2d!(sh, ps)
    compute_E2d!(sh)
    return sh
end

"""
    step_shock2d!(sh, ps, dt; NB=2)

Advance the 2D perpendicular shock one step: Boris push (gathered E,B) →
reflect at x=0 / wrap y / absorb past x=Lx → deposit moments → subcycle Bz.
"""
function step_shock2d!(
    sh::PerpShock2D{T},
    ps::ParticleSet{2,T},
    dt::Real;
    NB::Integer = 2,
) where {T}
    dx = sh.sbp.dx
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step_shock2d!")
    Lx = sh.Lx
    Ly = sh.Ly
    xp = ps.x[1]
    yp = ps.x[2]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    @inbounds for p in eachindex(ps.weight)
        Exp = _gather2d(sh.Ex, xp[p], yp[p], dx, sh.dy, sh.nx, sh.ny)
        Eyp = _gather2d(sh.Ey, xp[p], yp[p], dx, sh.dy, sh.nx, sh.ny)
        Bzp = _gather2d(sh.Bz, xp[p], yp[p], dx, sh.dy, sh.nx, sh.ny)
        nvx, nvy, nvz =
            boris_kick(vx[p], vy[p], vz[p], Exp, Eyp, zero(T), zero(T), zero(T), Bzp, one(T), dtT)
        vx[p] = nvx
        vy[p] = nvy
        vz[p] = nvz
        xp[p] += dtT * nvx
        yp[p] += dtT * nvy
    end
    # boundaries: reflect wall x=0, wrap y, absorb x>Lx
    @inbounds for p in eachindex(ps.weight)
        if xp[p] < 0
            xp[p] = -xp[p]
            vx[p] = -vx[p]
        end
        yp[p] = mod(yp[p], Ly)
    end
    write = 0
    @inbounds for p in eachindex(ps.weight)
        if xp[p] <= Lx
            write += 1
            if write != p
                for d = 1:2
                    ps.x[d][write] = ps.x[d][p]
                end
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
        for d = 1:2
            resize!(ps.x[d], write)
        end
        for c = 1:3
            resize!(ps.v[c], write)
        end
        resize!(ps.weight, write)
        resize!(ps.id, write)
        resize!(ps.tag, write)
    end
    deposit_moments2d!(sh, ps)
    hb = dtT / NB
    for _ = 1:NB
        _rk4_bz2d!(sh, hb)
    end
    compute_E2d!(sh)
    return sh
end

"Particle weight for a uniform load over [0,Lx]×[0,Ly] depositing density n0."
function shock2d_density_weight(n0, Lx, Ly, Np)
    isfinite(n0) && n0 >= 0 || throw(ArgumentError("n0 must be finite and non-negative"))
    isfinite(Lx) && Lx > 0 || throw(ArgumentError("Lx must be finite and positive"))
    isfinite(Ly) && Ly > 0 || throw(ArgumentError("Ly must be finite and positive"))
    N = _require_positive_intlike("Np", Np)
    return n0 * Lx * Ly / N
end

"""
    shock_surface(sh) -> (xs, mean_xs, sigma_xs)

Shock-front position x_s(y) (steepest |∂Bz/∂x| per transverse column), its
transverse mean, and its standard deviation σ_s (the rippling amplitude).
"""
function shock_surface(sh::PerpShock2D{T}) where {T}
    xs = Vector{T}(undef, sh.ny)
    @inbounds for j = 1:sh.ny
        # locate the compressed peak (downstream), then scan OUTWARD (toward
        # upstream) for the first drop below the ramp midpoint = the shock front.
        # Starting from the peak avoids both the low wall node and the far
        # rarefaction front.
        ipk = 1
        bzmax = sh.Bz[1, j]
        for i = 2:sh.nx
            sh.Bz[i, j] > bzmax && (bzmax = sh.Bz[i, j]; ipk = i)
        end
        thresh = (sh.B0 + bzmax) / 2
        xf = sh.x[sh.nx]
        for i = ipk:sh.nx
            if sh.Bz[i, j] < thresh
                xf = sh.x[i]
                break
            end
        end
        xs[j] = xf
    end
    m = sum(xs) / sh.ny
    σ = sqrt(sum(abs2, xs .- m) / sh.ny)
    return xs, m, σ
end
