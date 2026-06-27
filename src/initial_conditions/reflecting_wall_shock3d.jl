# shock3d.jl — 3D perpendicular collisionless hybrid shock (Phase 11, 3-D).
#
# Generalizes the verified 1-D `PerpShock` / 2-D `PerpShock2D` to a full 3-D3V
# perpendicular shock. Shock normal x is non-periodic (SBP-(2,1) + reflecting
# wall at x=0 + SAT inflow at x=Lx); transverse y AND z are periodic (Fourier).
# B0 = B0 ẑ is PERPENDICULAR to the normal (along a simulated transverse axis),
# so unlike the 2-D case all three magnetic components are dynamic:
#   J  = ∇×B
#   E_push = −u×B + (J×B)/n − ∇p_e/n + ηJ        (full Ohm, for the Boris push)
#   E_ind  = −u×B + ηJ                            (motional + resistive, stable)
#   ∂t B   = −∇×E_ind  (+ inflow SAT at x=Lx toward B0 ẑ)
#
# DIV-B PRESERVATION: with ∂x via SBP (dim 1) and ∂y,∂z via Fourier (dims 2,3),
# the three derivative operators act on independent index ranges, so they
# commute EXACTLY (∂x∂y = ∂y∂x, etc.). Hence ∇·(∇×E) = 0 to machine precision and
# the induction step conserves ∇·B at its initial value (=0 for uniform B0). This
# is verified in test_shock3d.jl (max|∇·B| stays ~1e-12). The Hall term is kept
# out of the field advance (as in 2-D) — it would inject a stiff whistler; it
# enters only the particle electric field. Reduces to the 2-D shock when
# z-uniform, and to the 1-D shock when y- and z-uniform.
#
const _Arr3{T} = Array{T,3}
const _Vec3{T} = NTuple{3,Array{T,3}}

mutable struct PerpShock3D{T}
    sbp::SBP1D{T}
    nx::Int
    ny::Int
    nz::Int
    Lx::T
    Ly::T
    Lz::T
    dy::T
    dz::T
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    B::_Vec3{T}        # Bx, By, Bz
    n::_Arr3{T}
    u::_Vec3{T}        # ux, uy, uz (bulk ion velocity)
    pe::_Arr3{T}
    E::_Vec3{T}        # particle-push electric field
    J::_Vec3{T}
    Eind::_Vec3{T}     # induction electric field (scratch)
    k1::_Vec3{T}
    k2::_Vec3{T}
    k3::_Vec3{T}
    k4::_Vec3{T}
    Btmp::_Vec3{T}
    s::_Arr3{T}        # derivative scratch
    wsum::_Arr3{T}
    mx::_Arr3{T}
    my::_Arr3{T}
    mz::_Arr3{T}
    Te::T
    γe::T
    η::T
    τ::T
    B0::T
    nfloor::T
end

function _require_valid_positive_shock_ma(MA::Real, ::Type{T}) where {T}
    MAT = T(MA)
    isfinite(MAT) && MAT > zero(T) || throw(ArgumentError("MA must be finite and positive"))
    return MAT
end

"""
    PerpShock3D(nx, ny, nz, Lx, Ly, Lz; Te, γe, η, τ, B0, nfloor)

3-D perpendicular hybrid-shock field state: `nx` SBP nodes over `[0,Lx]` (shock
normal) × `ny`×`nz` periodic Fourier nodes over `[0,Ly]×[0,Lz]` (transverse).
B0 is along ẑ (a transverse axis), perpendicular to the shock normal x.
"""
function PerpShock3D(
    nx::Integer,
    ny::Integer,
    nz::Integer,
    Lx::T,
    Ly::T,
    Lz::T;
    Te = 0.125,
    γe = 5 / 3,
    η = 0.02,
    τ = 3.0,
    B0 = 1.0,
    nfloor = 1e-6,
) where {T<:AbstractFloat}
    LxT = _require_finite_positive_real("Lx", Lx, T)
    LyT = _require_finite_positive_real("Ly", Ly, T)
    LzT = _require_finite_positive_real("Lz", Lz, T)
    TeT = _require_finite_nonnegative_real("Te", Te, T)
    γeT = _require_valid_gamma(γe, T)
    ηT = _require_finite_nonnegative_real("η", η, T)
    τT = _require_finite_real("τ", τ, T)
    B0T = _require_finite_real("B0", B0, T)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, T)
    sbp = SBP1D(nx, LxT)
    x = collect(range(zero(T), LxT; length = nx))
    y = [(j - 1) * LyT / ny for j = 1:ny]
    z = [(k - 1) * LzT / nz for k = 1:nz]
    A() = zeros(T, nx, ny, nz)
    V() = (A(), A(), A())
    Bx = A()
    By = A()
    Bz = fill(B0T, nx, ny, nz)            # B0 along ẑ
    PerpShock3D{T}(
        sbp,
        nx,
        ny,
        nz,
        LxT,
        LyT,
        LzT,
        LyT / ny,
        LzT / nz,
        x,
        y,
        z,
        (Bx, By, Bz),
        A(),
        V(),
        A(),
        V(),
        V(),
        V(),
        V(),
        V(),
        V(),
        V(),
        V(),
        A(),
        A(),
        A(),
        A(),
        A(),
        TeT,
        γeT,
        ηT,
        τT,
        B0T,
        nfloorT,
    )
end

# ---- derivative operators on 3-D arrays -----------------------------------

"SBP first derivative along x (dim 1) of a 3-D array, fibre by fibre."
function _sbp_dx3d!(out::_Arr3{T}, f::_Arr3{T}, s::SBP1D{T}) where {T}
    @inbounds for k in axes(f, 3), j in axes(f, 2)
        sbp_deriv!(view(out, :, j, k), view(f, :, j, k), s)
    end
    return out
end

# Spectral first derivative along periodic dim `d` (2=y or 3=z), Nyquist zeroed.
function _fourier_d!(out::_Arr3{T}, f::_Arr3{T}, L::T, d::Int) where {T}
    (d == 2 || d == 3) || throw(ArgumentError("d must be 2 or 3"))
    nd = size(f, d)
    fh = Array{Complex{T}}(undef, size(f))
    @inbounds for i in eachindex(f)
        fh[i] = Complex{T}(f[i], zero(T))
    end
    plan_fft!(fh, d) * fh
    nyquist = iseven(nd) ? nd ÷ 2 + 1 : 0
    if d == 2
        @inbounds for k3 in axes(fh, 3), j in axes(fh, 2), i in axes(fh, 1)
            m = j - 1
            mp = m <= nd ÷ 2 ? m : m - nd
            kval = j == nyquist ? zero(T) : T(2π) * mp / L
            fh[i, j, k3] *= Complex{T}(zero(T), kval)
        end
    else
        @inbounds for k3 in axes(fh, 3), j in axes(fh, 2), i in axes(fh, 1)
            m = k3 - 1
            mp = m <= nd ÷ 2 ? m : m - nd
            kval = k3 == nyquist ? zero(T) : T(2π) * mp / L
            fh[i, j, k3] *= Complex{T}(zero(T), kval)
        end
    end
    plan_ifft!(fh, d) * fh
    @inbounds for i in eachindex(out, fh)
        out[i] = real(fh[i])
    end
    return out
end

# J = ∇×B, written into `J`, using `sh.s` as scratch.
function _curl3d!(J::_Vec3{T}, B::_Vec3{T}, sh::PerpShock3D{T}) where {T}
    Bx, By, Bz = B
    s = sh.s
    # Jx = ∂y Bz − ∂z By
    _fourier_d!(J[1], Bz, sh.Ly, 2)
    _fourier_d!(s, By, sh.Lz, 3)
    J[1] .-= s
    # Jy = ∂z Bx − ∂x Bz
    _fourier_d!(J[2], Bx, sh.Lz, 3)
    _sbp_dx3d!(s, Bz, sh.sbp)
    J[2] .-= s
    # Jz = ∂x By − ∂y Bx
    _sbp_dx3d!(J[3], By, sh.sbp)
    _fourier_d!(s, Bx, sh.Ly, 2)
    J[3] .-= s
    return J
end

"Magnetic-field divergence ∂x Bx + ∂y By + ∂z Bz (uses `sh.s`; returns a fresh array)."
function magnetic_divergence3d(sh::PerpShock3D{T}) where {T}
    div = similar(sh.B[1])
    _sbp_dx3d!(div, sh.B[1], sh.sbp)
    _fourier_d!(sh.s, sh.B[2], sh.Ly, 2)
    div .+= sh.s
    _fourier_d!(sh.s, sh.B[3], sh.Lz, 3)
    div .+= sh.s
    return div
end

# ---- deposition / gather (trilinear CIC) ----------------------------------

@inline function _cicz(zp::T, dz::T, nz::Int) where {T}
    s = zp / dz
    k0 = floor(Int, s)
    f = s - k0
    a = mod(k0, nz) + 1
    b = mod(k0 + 1, nz) + 1
    return a, b, one(T) - f, f
end

"Deposit n and bulk velocity (ux,uy,uz) onto the SBP-x × periodic-y,z mesh."
function deposit_moments3d!(sh::PerpShock3D{T}, ps::ParticleSet{3,T}) where {T}
    fill!(sh.wsum, zero(T))
    fill!(sh.mx, zero(T))
    fill!(sh.my, zero(T))
    fill!(sh.mz, zero(T))
    dx = sh.sbp.dx
    xp = ps.x[1]
    yp = ps.x[2]
    zp = ps.x[3]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    w = ps.weight
    @inbounds for p in eachindex(w)
        ax, bx, wxa, wxb = _cicx(xp[p], dx, sh.nx)
        ay, by, wya, wyb = _cicy(yp[p], sh.dy, sh.ny)
        az, bz, wza, wzb = _cicz(zp[p], sh.dz, sh.nz)
        ww = w[p]
        for (ix, wx) in ((ax, wxa), (bx, wxb)),
            (iy, wy) in ((ay, wya), (by, wyb)),
            (iz, wz) in ((az, wza), (bz, wzb))

            g = wx * wy * wz * ww
            sh.wsum[ix, iy, iz] += g
            sh.mx[ix, iy, iz] += g * vx[p]
            sh.my[ix, iy, iz] += g * vy[p]
            sh.mz[ix, iy, iz] += g * vz[p]
        end
    end
    nf = sh.nfloor
    H = sh.sbp.H
    vol_yz = sh.dy * sh.dz
    @inbounds for iz = 1:sh.nz, iy = 1:sh.ny, ix = 1:sh.nx
        ws = sh.wsum[ix, iy, iz]
        sh.n[ix, iy, iz] = ws / (H[ix] * vol_yz)
        wsf = max(ws, nf)
        sh.u[1][ix, iy, iz] = sh.mx[ix, iy, iz] / wsf
        sh.u[2][ix, iy, iz] = sh.my[ix, iy, iz] / wsf
        sh.u[3][ix, iy, iz] = sh.mz[ix, iy, iz] / wsf
    end
    return sh
end

@inline function _gather3d(
    F::_Arr3{T},
    xp::T,
    yp::T,
    zp::T,
    dx::T,
    dy::T,
    dz::T,
    nx::Int,
    ny::Int,
    nz::Int,
) where {T}
    ax, bx, wxa, wxb = _cicx(xp, dx, nx)
    ay, by, wya, wyb = _cicy(yp, dy, ny)
    az, bz, wza, wzb = _cicz(zp, dz, nz)
    s = zero(T)
    @inbounds for (ix, wx) in ((ax, wxa), (bx, wxb)),
        (iy, wy) in ((ay, wya), (by, wyb)),
        (iz, wz) in ((az, wza), (bz, wzb))

        s += wx * wy * wz * F[ix, iy, iz]
    end
    return s
end

# ---- electric field (full Ohm, for the particle push) ---------------------

"Particle-push electric field from the full 3-D hybrid Ohm's law."
function compute_E3d!(sh::PerpShock3D{T}) where {T}
    _curl3d!(sh.J, sh.B, sh)                      # J = ∇×B
    @. sh.pe = sh.Te * sh.n^sh.γe
    ux, uy, uz = sh.u
    Bx, By, Bz = sh.B
    Jx, Jy, Jz = sh.J
    # ∇p_e written into Eind as scratch, then folded in
    _sbp_dx3d!(sh.Eind[1], sh.pe, sh.sbp)
    _fourier_d!(sh.Eind[2], sh.pe, sh.Ly, 2)
    _fourier_d!(sh.Eind[3], sh.pe, sh.Lz, 3)
    dpx, dpy, dpz = sh.Eind
    η = sh.η
    nf = sh.nfloor
    @inbounds for q in eachindex(sh.n)
        inv = one(T) / max(sh.n[q], nf)
        bx = Bx[q]
        by = By[q]
        bz = Bz[q]
        jx = Jx[q]
        jy = Jy[q]
        jz = Jz[q]
        # E = −u×B + (J×B)/n − ∇p_e/n + ηJ
        sh.E[1][q] = -(uy[q] * bz - uz[q] * by) + (jy * bz - jz * by) * inv - dpx[q] * inv + η * jx
        sh.E[2][q] = -(uz[q] * bx - ux[q] * bz) + (jz * bx - jx * bz) * inv - dpy[q] * inv + η * jy
        sh.E[3][q] = -(ux[q] * by - uy[q] * bx) + (jx * by - jy * bx) * inv - dpz[q] * inv + η * jz
    end
    return sh
end

# ---- induction RHS + RK4 ---------------------------------------------------

# ∂t B = −∇×E_ind, E_ind = −u×B + ηJ, evaluated at trial field `Bt`.
function _b_rhs3d!(K::_Vec3{T}, Bt::_Vec3{T}, sh::PerpShock3D{T}) where {T}
    _curl3d!(sh.J, Bt, sh)                        # J = ∇×Bt
    ux, uy, uz = sh.u
    Bx, By, Bz = Bt
    Jx, Jy, Jz = sh.J
    η = sh.η
    @inbounds for q in eachindex(Bx)
        sh.Eind[1][q] = -(uy[q] * Bz[q] - uz[q] * By[q]) + η * Jx[q]
        sh.Eind[2][q] = -(uz[q] * Bx[q] - ux[q] * Bz[q]) + η * Jy[q]
        sh.Eind[3][q] = -(ux[q] * By[q] - uy[q] * Bx[q]) + η * Jz[q]
    end
    _curl3d!(K, sh.Eind, sh)                      # K = ∇×E_ind
    @inbounds for c = 1:3
        K[c] .*= -one(T)                          # ∂t B = −∇×E_ind
    end
    # SAT inflow at x = Lx face: relax B toward B0 ẑ
    τH = sh.τ / sh.sbp.H[sh.nx]
    B0vec = (zero(T), zero(T), sh.B0)
    @inbounds for c = 1:3, iz = 1:sh.nz, iy = 1:sh.ny
        K[c][sh.nx, iy, iz] += -τH * (Bt[c][sh.nx, iy, iz] - B0vec[c])
    end
    return K
end

@inline function _axpy_tuple!(out::_Vec3{T}, B::_Vec3{T}, a::T, k::_Vec3{T}) where {T}
    @inbounds for c = 1:3
        @. out[c] = B[c] + a * k[c]
    end
    return out
end

function _rk4_b3d!(sh::PerpShock3D{T}, h::T) where {T}
    B = sh.B
    _b_rhs3d!(sh.k1, B, sh)
    _axpy_tuple!(sh.Btmp, B, h / 2, sh.k1)
    _b_rhs3d!(sh.k2, sh.Btmp, sh)
    _axpy_tuple!(sh.Btmp, B, h / 2, sh.k2)
    _b_rhs3d!(sh.k3, sh.Btmp, sh)
    _axpy_tuple!(sh.Btmp, B, h, sh.k3)
    _b_rhs3d!(sh.k4, sh.Btmp, sh)
    @inbounds for c = 1:3
        @. B[c] += h / 6 * (sh.k1[c] + 2 * sh.k2[c] + 2 * sh.k3[c] + sh.k4[c])
    end
    return sh
end

# A-stable trapezoidal (Crank–Nicolson) field substep, solved by fixed-point
# (Picard) iteration. The field RHS is LINEAR in B (moments are frozen during the
# substep), so the trapezoid B^{n+1} = B^n + h/2(f(Bⁿ)+f(B^{n+1})) converges in a
# few iterations at the substep size. This is the conservative/semi-implicit
# integrator option used to compare integrators on an identical shock.
function _cn_b3d!(sh::PerpShock3D{T}, h::T; iters::Int = 4) where {T}
    B = sh.B
    _b_rhs3d!(sh.k1, B, sh)                       # f(Bⁿ)
    @inbounds for c = 1:3
        @. sh.Btmp[c] = B[c]                      # initial iterate = Bⁿ
    end
    for _ = 1:iters
        _b_rhs3d!(sh.k2, sh.Btmp, sh)            # f(iterate)
        @inbounds for c = 1:3
            @. sh.Btmp[c] = B[c] + h / 2 * (sh.k1[c] + sh.k2[c])
        end
    end
    @inbounds for c = 1:3
        @. B[c] = sh.Btmp[c]
    end
    return sh
end

# ---- public step / init ----------------------------------------------------

"Initialize the carried E from the loaded particles."
function init_shock3d!(sh::PerpShock3D{T}, ps::ParticleSet{3,T}) where {T}
    deposit_moments3d!(sh, ps)
    compute_E3d!(sh)
    return sh
end

"""
    step_shock3d!(sh, ps, dt; NB=2)

Advance the 3-D perpendicular shock one step: Boris push (gathered E,B) →
reflect at x=0 / wrap y,z / absorb past x=Lx → deposit moments → subcycle B
(RK4) → recompute the push E.
"""
function step_shock3d!(
    sh::PerpShock3D{T},
    ps::ParticleSet{3,T},
    dt::Real;
    NB::Integer = 2,
    field_method::Symbol = :rk4,
) where {T}
    field_method in (:rk4, :cn) || throw(ArgumentError("field_method must be :rk4 or :cn"))
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "step_shock3d!")
    dx = sh.sbp.dx
    dy = sh.dy
    dz = sh.dz
    Lx = sh.Lx
    Ly = sh.Ly
    Lz = sh.Lz
    nx = sh.nx
    ny = sh.ny
    nz = sh.nz
    xp = ps.x[1]
    yp = ps.x[2]
    zp = ps.x[3]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    @inbounds for p in eachindex(ps.weight)
        Exp = _gather3d(sh.E[1], xp[p], yp[p], zp[p], dx, dy, dz, nx, ny, nz)
        Eyp = _gather3d(sh.E[2], xp[p], yp[p], zp[p], dx, dy, dz, nx, ny, nz)
        Ezp = _gather3d(sh.E[3], xp[p], yp[p], zp[p], dx, dy, dz, nx, ny, nz)
        Bxp = _gather3d(sh.B[1], xp[p], yp[p], zp[p], dx, dy, dz, nx, ny, nz)
        Byp = _gather3d(sh.B[2], xp[p], yp[p], zp[p], dx, dy, dz, nx, ny, nz)
        Bzp = _gather3d(sh.B[3], xp[p], yp[p], zp[p], dx, dy, dz, nx, ny, nz)
        nvx, nvy, nvz = boris_kick(vx[p], vy[p], vz[p], Exp, Eyp, Ezp, Bxp, Byp, Bzp, one(T), dtT)
        vx[p] = nvx
        vy[p] = nvy
        vz[p] = nvz
        xp[p] += dtT * nvx
        yp[p] += dtT * nvy
        zp[p] += dtT * nvz
    end
    # boundaries: reflect wall x=0, wrap y,z, absorb x>Lx
    @inbounds for p in eachindex(ps.weight)
        if xp[p] < 0
            xp[p] = -xp[p]
            vx[p] = -vx[p]
        end
        yp[p] = mod(yp[p], Ly)
        zp[p] = mod(zp[p], Lz)
    end
    write = 0
    @inbounds for p in eachindex(ps.weight)
        if xp[p] <= Lx
            write += 1
            if write != p
                for d = 1:3
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
        for d = 1:3
            resize!(ps.x[d], write)
        end
        for c = 1:3
            resize!(ps.v[c], write)
        end
        resize!(ps.weight, write)
        resize!(ps.id, write)
        resize!(ps.tag, write)
    end
    deposit_moments3d!(sh, ps)
    hb = dtT / NB
    for _ = 1:NB
        field_method === :rk4 ? _rk4_b3d!(sh, hb) : _cn_b3d!(sh, hb)
    end
    compute_E3d!(sh)
    return sh
end

"Particle weight for a uniform load over [0,Lx]×[0,Ly]×[0,Lz] depositing density n0."
shock3d_density_weight(n0, Lx, Ly, Lz, Np) = n0 * Lx * Ly * Lz / Np

"""
    shock_surface3d(sh) -> (xs, mean_xs, sigma_xs)

Shock-front position `x_s(y,z)` per transverse column (the first drop below the
ramp midpoint, scanning outward from the compressed peak), its transverse mean,
and its standard deviation σ_s (the 3-D rippling amplitude). `xs` is `ny×nz`.
"""
function shock_surface3d(sh::PerpShock3D{T}) where {T}
    Bz = sh.B[3]
    xs = Matrix{T}(undef, sh.ny, sh.nz)
    @inbounds for kz = 1:sh.nz, jy = 1:sh.ny
        ipk = 1
        bzmax = Bz[1, jy, kz]
        for i = 2:sh.nx
            if Bz[i, jy, kz] > bzmax
                bzmax = Bz[i, jy, kz]
                ipk = i
            end
        end
        thresh = (sh.B0 + bzmax) / 2
        xf = sh.x[sh.nx]
        for i = ipk:sh.nx
            if Bz[i, jy, kz] < thresh
                xf = sh.x[i]
                break
            end
        end
        xs[jy, kz] = xf
    end
    m = sum(xs) / length(xs)
    σ = sqrt(sum(abs2, xs .- m) / length(xs))
    return xs, m, σ
end

# ---- driver ----------------------------------------------------------------

# Build a 3-D shock field state + drifting-ion load (shared by the driver and
# the campaign/restart helpers). Returns (sh, ps).
function _load_shock3d(;
    MA::Real,
    nx::Integer,
    ny::Integer,
    nz::Integer,
    Lx::Real,
    Ly::Real,
    Lz::Real,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    nppc::Integer = 8,
    seed::Integer = 1,
    db_turb::Real = 0.0,
)
    T = Float64
    LxT = _require_finite_positive_real("Lx", Lx, T)
    LyT = _require_finite_positive_real("Ly", Ly, T)
    LzT = _require_finite_positive_real("Lz", Lz, T)
    B0 = one(T)
    U0 = _require_valid_positive_shock_ma(MA, T)
    vth = _require_finite_nonnegative_real("vthi", vthi, T)
    db_turbT = _require_finite_nonnegative_real("db_turb", db_turb, T)
    sh = PerpShock3D(nx, ny, nz, LxT, LyT, LzT; Te, γe, η, τ = U0, B0 = B0)

    # Optional upstream turbulence: div-free transverse Alfvénic fluctuations
    # (δBz varies only in y, δBy only in z ⇒ ∇·δB = 0), summed over a few
    # harmonics with random phases and total rms amplitude ≈ db_turb·B0.
    if db_turbT > 0
        rngt = MersenneTwister(Int(seed) + 7919)
        nh = 3
        amp = db_turbT / sqrt(nh)
        for h = 1:nh
            φy = T(2π) * rand(rngt)
            φz = T(2π) * rand(rngt)
            ky = T(2π) * h / LyT
            kz = T(2π) * h / LzT
            @inbounds for kk = 1:nz, jj = 1:ny, ii = 1:nx
                sh.B[3][ii, jj, kk] += amp * sin(ky * sh.y[jj] + φy)
                sh.B[2][ii, jj, kk] += amp * sin(kz * sh.z[kk] + φz)
            end
        end
    end

    Np = Int(nppc) * Int(nx) * Int(ny) * Int(nz)
    ps = ParticleSet{3,T}(Np)
    rng = MersenneTwister(Int(seed))
    @inbounds for p = 1:Np
        ps.x[1][p] = LxT * rand(rng)
        ps.x[2][p] = LyT * rand(rng)
        ps.x[3][p] = LzT * rand(rng)
        ps.v[1][p] = -U0 + vth * randn(rng)
        ps.v[2][p] = vth * randn(rng)
        ps.v[3][p] = vth * randn(rng)
    end
    ps.weight .= shock3d_density_weight(one(T), LxT, LyT, LzT, Np)
    init_shock3d!(sh, ps)
    return sh, ps
end

"""
    run_perp_shock3d(; MA, nx=64, ny=16, nz=16, Lx=80, Ly=12, Lz=12, Te=0.125,
                       γe=5/3, vthi=0.35, η=0.02, nppc=8, nsteps=400, dt=0.02,
                       seed=1) -> NamedTuple

Run a 3-D perpendicular collisionless hybrid shock at upstream Alfvén Mach
number `MA` (upstream drift `U0 = MA·v_A`, `v_A=1` in Ω_ci units) and return
downstream / shock-front diagnostics:

  • `n2`            — downstream compression ρ₂/ρ₁ (upstream n₁=1),
  • `Bz2`           — downstream perpendicular field,
  • `frozen_ratio`  — `(Bz2/B0)/n2` (=1 when the field is frozen to the flow),
  • `X_rh`          — fluid RH compression at the realized Mach number,
  • `sigma_xs`      — transverse rippling amplitude of the front `x_s(y,z)`,
  • `mean_xs`       — mean front position,
  • `maxdivB`       — max |∇·B| (≈0 confirms the constraint is preserved),
  • `M_real`        — realized shock-frame Mach number.

Grid defaults are modest for fast verification; raise `nx,ny,nz,nppc,nsteps` for
a production case (the code is unchanged — only compute time scales).
"""
function run_perp_shock3d(;
    MA::Real,
    nx::Integer = 64,
    ny::Integer = 16,
    nz::Integer = 16,
    Lx::Real = 80.0,
    Ly::Real = 12.0,
    Lz::Real = 12.0,
    Te::Real = 0.125,
    γe::Real = 5 / 3,
    vthi::Real = 0.35,
    η::Real = 0.02,
    nppc::Integer = 8,
    nsteps::Integer = 400,
    dt::Real = 0.02,
    seed::Integer = 1,
    field_method::Symbol = :rk4,
    db_turb::Real = 0.0,
)
    nx >= 3 || throw(ArgumentError("nx must be at least 3"))
    ny >= 1 || throw(ArgumentError("ny must be positive"))
    nz >= 1 || throw(ArgumentError("nz must be positive"))
    nppc >= 1 || throw(ArgumentError("nppc must be positive"))
    nsteps >= 0 || throw(ArgumentError("nsteps must be non-negative"))
    T = Float64
    MAT = _require_valid_positive_shock_ma(MA, T)
    B0 = one(T)
    vA = one(T)
    U0 = MAT * vA
    vth = T(vthi)
    sh, ps = _load_shock3d(; MA, nx, ny, nz, Lx, Ly, Lz, Te, γe, vthi, η, nppc, seed, db_turb)

    for _ = 1:nsteps
        step_shock3d!(sh, ps, T(dt); NB = 2, field_method = field_method)
    end
    _require_all_finite("Bz", sh.B[3], "unstable 3-D run")

    # Downstream slab average. Use the transverse-averaged Bz/n profiles along x.
    # The shock front is the first outward drop of Bz below the ramp midpoint,
    # scanning from the compressed peak (robust: the far upstream depletes as ions
    # drift into the wall, so a global gradient locator is unreliable here). The
    # downstream slab runs from a few cells off the wall to just inside the front.
    nbar = dropdims(sum(sh.n; dims = (2, 3)); dims = (2, 3)) ./ (ny * nz)
    Bzbar = dropdims(sum(sh.B[3]; dims = (2, 3)); dims = (2, 3)) ./ (ny * nz)
    ipk = argmax(Bzbar)
    thr = (B0 + Bzbar[ipk]) / 2
    ifr = nx
    for i = ipk:nx
        if Bzbar[i] < thr
            ifr = i
            break
        end
    end
    xf = sh.x[ifr]
    ilo = findfirst(>(T(3)), sh.x)
    ilo === nothing && (ilo = 1)
    ihi = max(ilo, ifr - 1)
    slab = ilo:ihi
    n2 = sum(@view nbar[slab]) / length(slab)
    Bz2 = sum(@view Bzbar[slab]) / length(slab)

    # shock speed from mass conservation Vs = U0/(n2−1) (front tracking in 3-D is
    # noisy; the mass-flux estimate is the robust EOS-independent measure).
    Vs = isfinite(n2) && n2 > 1 ? U0 / (n2 - 1) : T(NaN)
    M_real = (U0 + Vs) / vA
    p1 = vth^2 + T(Te)
    # fluid RH oracle (only defined for a real compressive shock; NaN otherwise)
    X_rh =
        isfinite(M_real) && (U0 + Vs) > 0 ?
        rankine_hugoniot(MHDState(one(T), U0 + Vs, zero(T), p1, zero(T), B0), T(γe)).X : T(NaN)

    _, mean_xs, sigma_xs = shock_surface3d(sh)
    maxdivB = maximum(abs, magnetic_divergence3d(sh))

    return (;
        n2,
        Bz2,
        frozen_ratio = (Bz2 / B0) / n2,
        X_rh,
        sigma_xs,
        mean_xs,
        maxdivB,
        M_real,
        Vs,
        xf,
        sh,
    )
end
