# shock_diag.jl — diagnostics for the perpendicular hybrid shocks (PerpShock /
# PerpShock2D).
#
# Conventions shared with shock_sim.jl / shock2d.jl:
#   * shock normal is x; reflecting WALL at x=0 (piston/downstream side); held
#     INFLOW at x=Lx (upstream side). Upstream plasma streams in the −x direction
#     (bulk u_x < 0) toward the wall.
#   * perpendicular field B = B_z ẑ ⟂ normal, normalized units μ0 = 1, ion m = 1.
#   * fluid moments n, u_x, u_y carried on the SBP-x nodes.
#
# Energy densities (normalized, μ0=1): magnetic w_B = B_z²/2; the field/EM energy
# flux along +x is the Poynting flux S_x = (E×B)_x = E_y B_z (since B=B_zẑ,
# E_z=0). The fluid kinetic-energy flux along +x is F_K = (½ n u²) u_x with
# u² = u_x² + u_y². Both are reported as the NORMAL (+x) flux evaluated at the
# inflow node (x=Lx) and the wall node (x=0); a steady shock balances them.

# ---------------------------------------------------------------- 1D boundary energy flux

"""
    boundary_energy_flux(sh::PerpShock; ps=nothing)
        -> (; magnetic, kinetic, enthalpy, total)

Normal (+x) energy flux through the two x-boundaries of a 1-D perpendicular
shock — the FULL conserved flux, so the upstream/downstream entries balance for
a steady shock. Each entry is a `(inflow, wall)` pair (flux at the inflow node
x=Lx and the reflecting-wall node x=0):

* `magnetic`  — electromagnetic (Poynting) flux `S_x = (E×B)_x = E_y B_z`
  (μ0=1, B=B_zẑ, E_z=0).
* `kinetic`   — the ION kinetic-energy flux. With `ps` (the kinetic ions) it is
  the FULL flux `Σ_p w·½|v|²·v_x` (bulk + thermal) deposited at the boundary
  with the same CIC/H-norm as the moments; without `ps` it falls back to the
  bulk-only fluid estimate `½ n (u_x²+u_y²) u_x` (no ion thermal flux).
* `enthalpy`  — the electron ENTHALPY flux `γ_e/(γ_e−1)·p_e·u_x` = internal
  energy `p_e/(γ_e−1)·u_x` plus pressure work `p_e·u_x` (the term previously
  omitted). Non-finite for the isothermal limit γ_e→1 (no internal-energy
  invariant), as in `energy_budget`.
* `total`     — `magnetic .+ kinetic .+ enthalpy`.

Sign convention: positive = energy carried in +x. Inflowing upstream plasma
(u_x<0) carries energy in −x (toward the wall), so its inflow entries are
negative.
"""
function boundary_energy_flux(
    sh::PerpShock{T};
    ps::Union{Nothing,ParticleSet{1,T}} = nothing,
) where {T}
    N = sh.s.n
    γe = sh.γe
    # electron enthalpy flux  h_e = γe/(γe−1) · p_e · u_x  (internal energy + p·u work)
    enth(i) = (γe / (γe - one(T))) * sh.pe[i] * sh.ux[i]
    # magnetic / Poynting flux S_x = E_y B_z at each boundary
    mag_inflow = sh.Ey[N] * sh.Bz[N]
    mag_wall = sh.Ey[1] * sh.Bz[1]
    enth_inflow = enth(N)
    enth_wall = enth(1)
    if ps === nothing
        # bulk-only ion kinetic flux ½ n |u|² u_x (ion thermal flux needs the ions)
        bulk(i) = T(0.5) * sh.n[i] * (sh.ux[i]^2 + sh.uy[i]^2) * sh.ux[i]
        kin_inflow = bulk(N)
        kin_wall = bulk(1)
    else
        # full kinetic-ion energy flux Σ_p w ½|v|² v_x (bulk + thermal) at the boundary
        Fe = _ion_energy_flux(sh, ps)
        kin_inflow = Fe[N]
        kin_wall = Fe[1]
    end
    magnetic = (mag_inflow, mag_wall)
    kinetic = (kin_inflow, kin_wall)
    enthalpy = (enth_inflow, enth_wall)
    total = (mag_inflow + kin_inflow + enth_inflow, mag_wall + kin_wall + enth_wall)
    return (; magnetic, kinetic, enthalpy, total)
end

# Full kinetic-ion energy-flux density Σ_p w·½|v|²·v_x deposited on the SBP grid
# (the same CIC stencil + H-norm as deposit_moments!), returned as a length-N
# vector. CIC's partition of unity makes Σ_i H[i]·F[i] = Σ_p w·½|v|²·v_x exactly.
function _ion_energy_flux(sh::PerpShock{T}, ps::ParticleSet{1,T}) where {T}
    N = sh.s.n
    dx = sh.s.dx
    acc = zeros(T, N)
    xp = ps.x[1]
    vx = ps.v[1]
    vy = ps.v[2]
    vz = ps.v[3]
    w = ps.weight
    _require_finite_real_sequence("particle positions", xp)
    _require_finite_real_sequence("particle velocities vx", vx)
    _require_finite_real_sequence("particle velocities vy", vy)
    _require_finite_real_sequence("particle velocities vz", vz)
    _require_finite_real_sequence("particle weights", w)
    @inbounds for p in eachindex(w)
        a, b, wa, wb = _cic_sbp(xp[p], dx, N)
        keflux = T(0.5) * (vx[p]^2 + vy[p]^2 + vz[p]^2) * vx[p] * w[p]
        acc[a] += wa * keflux
        acc[b] += wb * keflux
    end
    @inbounds for i = 1:N
        acc[i] /= sh.s.H[i]
    end
    return acc
end

# ---------------------------------------------------------------- 2D surface spectrum

"""
    shock_surface_spectrum(sh::PerpShock2D) -> (; ky, Ps, xs, mean_xs)

Transverse power spectrum `P_s(k_y) = |FFT_y(x_s(y) − ⟨x_s⟩)|²` of the shock
front `x_s(y)` returned by [`shock_surface`](@ref). `ky` are the non-negative
angular wavenumbers `2π m / Ly` (m = 0 … ny÷2), `Ps` the folded one-sided power
at each `ky`, `xs` the raw front, and `mean_xs` its transverse mean. The mean is
removed so `Ps[1]` (the k_y=0 / DC bin) is ≈0 for a rippled-but-unbiased front.
"""
function shock_surface_spectrum(sh::PerpShock2D{T}) where {T}
    xs, m, _ = shock_surface(sh)
    ny = sh.ny
    f = xs .- m
    fh = fft(f)
    nk = ny ÷ 2 + 1
    Ps = zeros(T, nk)
    @inbounds for idx = 1:ny
        mm = idx - 1                          # 0-based mode index
        km = mm <= ny ÷ 2 ? mm : ny - mm       # fold to non-negative |m|
        km + 1 <= nk && (Ps[km+1] += abs2(fh[idx]))
    end
    ky = T[T(2π) * mm / sh.Ly for mm = 0:nk-1]
    return (; ky, Ps, xs, mean_xs = m)
end

# ---------------------------------------------------------------- 2D transverse coherence

"""
    transverse_coherence(sh::PerpShock2D) -> (; dy, C)

Normalized transverse autocorrelation `C_s(Δy)` of the shock front `x_s(y)`
(mean removed), for lags `Δy = 0, dy, …, (ny−1)·dy`:

    C_s(Δy) = ⟨ x̃_s(y) x̃_s(y+Δy) ⟩_y / ⟨ x̃_s(y)² ⟩_y ,

with `x̃_s = x_s − ⟨x_s⟩` and the average over the periodic y-direction (the
lag wraps modulo ny). `C[1]` (Δy=0) = 1 by construction; a flat front returns
all-NaN (zero variance). `dy[k] = (k−1)·sh.dy`.
"""
function transverse_coherence(sh::PerpShock2D{T}) where {T}
    xs, m, _ = shock_surface(sh)
    ny = sh.ny
    f = xs .- m
    var0 = sum(abs2, f) / ny
    C = Vector{T}(undef, ny)
    if var0 == 0
        fill!(C, T(NaN))
    else
        @inbounds for lag = 0:ny-1
            acc = zero(T)
            for j = 1:ny
                jj = mod(j - 1 + lag, ny) + 1
                acc += f[j] * f[jj]
            end
            C[lag+1] = (acc / ny) / var0
        end
    end
    dy = T[T(k) * sh.dy for k = 0:ny-1]
    return (; dy, C)
end

# ---------------------------------------------------------------- crossing logger

"""
    CrossingLogger{T}

Tracks particles crossing a moving control surface `x_surface` between successive
[`log_crossings!`](@ref) calls. For each known particle (keyed by `ps.id`) it
stores the previous signed offset `x − x_surface` and the previous kinetic energy
`½ m |v|²`. A crossing is a SIGN CHANGE of that offset between two calls; at each
crossing the kinetic-energy change since the previous call is accumulated.

Fields: `count` (total crossings logged), `gain` (Σ ΔKE over all crossings),
`prev_off`/`prev_ke` (per-id last-seen offset and kinetic energy).
"""
mutable struct CrossingLogger{T}
    count::Int
    gain::T
    prev_off::Dict{UInt64,T}
    prev_ke::Dict{UInt64,T}
end

"""
    CrossingLogger{T}() / CrossingLogger(T=Float64)

Create an empty crossing logger of element type `T`.
"""
CrossingLogger{T}() where {T} = CrossingLogger{T}(0, zero(T), Dict{UInt64,T}(), Dict{UInt64,T}())
CrossingLogger(::Type{T} = Float64) where {T} = CrossingLogger{T}()

@inline function _ke(ps::ParticleSet{D,T}, p::Int) where {D,T}
    vx = ps.v[1][p]
    vy = ps.v[2][p]
    vz = ps.v[3][p]
    return T(0.5) * ps.m * (vx * vx + vy * vy + vz * vz)
end

"""
    log_crossings!(logger, ps, x_surface) -> nnew

Record particles of `ps` whose signed offset `x − x_surface` changed sign since
the previous call (a crossing of the control surface). `x_surface` may be a
scalar (flat surface) or a per-particle vector `x_surface[p]` (e.g. a rippled
front sampled at each particle's y). Uses `ps.x[1]` as the normal coordinate.
Returns the number of NEW crossings recorded on this call; updates
`logger.count` and accumulates ΔKE (since the previous call) into `logger.gain`
for every crossing particle. Particles seen for the first time are registered
without counting a crossing (no previous offset to compare).
"""
function log_crossings!(logger::CrossingLogger{T}, ps::ParticleSet{D,T}, x_surface) where {D,T}
    x = ps.x[1]
    ids = ps.id
    _require_finite_real_sequence("particle positions", x)
    _require_finite_real_sequence("particle velocities vx", ps.v[1])
    _require_finite_real_sequence("particle velocities vy", ps.v[2])
    _require_finite_real_sequence("particle velocities vz", ps.v[3])
    nnew = 0
    surf = if x_surface isa AbstractVector
        length(x_surface) == length(ids) ||
            throw(DimensionMismatch("x_surface length must match the particle count"))
        _require_finite_real_sequence("x_surface", x_surface)
        p -> T(x_surface[p])
    else
        xsurf = _require_finite_real("x_surface", x_surface, T)
        _ -> xsurf
    end
    @inbounds for p in eachindex(ids)
        id = ids[p]
        off = x[p] - surf(p)
        ke = _ke(ps, p)
        prev = get(logger.prev_off, id, nothing)
        if prev !== nothing
            # a crossing is a strict sign change (an exact-zero touch is not a
            # crossing until the sign actually flips)
            if (prev > 0 && off < 0) || (prev < 0 && off > 0)
                logger.count += 1
                nnew += 1
                logger.gain += ke - logger.prev_ke[id]
            end
        end
        logger.prev_off[id] = off
        logger.prev_ke[id] = ke
    end
    return nnew
end

"Total number of surface crossings recorded by the logger."
crossing_count(logger::CrossingLogger) = logger.count

"Accumulated kinetic-energy change ΣΔKE over all logged crossings."
energy_gain(logger::CrossingLogger) = logger.gain

# ---------------------------------------------------------------- spurious reflection metric

"""
    boundary_reflection_fraction(sh, ps; ncells=3) -> frac

Fraction of particles within `ncells` grid cells of the INFLOW boundary (x=Lx)
that are moving back UPSTREAM (toward larger x, `v_x > 0`) — a spurious-reflection
metric for the inflow boundary, which should be small for a clean inflow.
`sh` may be a `PerpShock` (with `ps::ParticleSet{1}`) or a `PerpShock2D` (with
`ps::ParticleSet{2}`). The boundary band is `[Lx − ncells·dx, Lx]`. Returns 0
when no particle lies in the band.
"""
function boundary_reflection_fraction(
    sh::PerpShock{T},
    ps::ParticleSet{1,T};
    ncells::Integer = 3,
) where {T}
    return _refl_frac(ps.x[1], ps.v[1], sh.x[end], sh.s.dx, ncells)
end

function boundary_reflection_fraction(
    sh::PerpShock2D{T},
    ps::ParticleSet{2,T};
    ncells::Integer = 3,
) where {T}
    return _refl_frac(ps.x[1], ps.v[1], sh.Lx, sh.sbp.dx, ncells)
end

function _refl_frac(x::Vector{T}, vx::Vector{T}, Lx::T, dx::T, ncells::Integer) where {T}
    ncells >= 1 || throw(ArgumentError("ncells must be positive"))
    _require_finite_real_sequence("particle positions", x)
    _require_finite_real_sequence("particle velocities", vx)
    Lx = _require_finite_real("Lx", Lx, T)
    dx = _require_finite_positive_real("dx", dx, T)
    edge = Lx - T(ncells) * dx
    inband = 0
    back = 0
    @inbounds for p in eachindex(x)
        if x[p] >= edge && x[p] <= Lx
            inband += 1
            vx[p] > 0 && (back += 1)
        end
    end
    return inband == 0 ? zero(T) : T(back) / T(inband)
end

# ---------------------------------------------------------------- normal-incidence frame

"""
    normal_incidence_frame(u::NTuple{3}, B::NTuple{3}, n_hat::NTuple{3}) -> NTuple{3}

Normal-incidence-frame (NIF) boost velocity. The NIF is the frame, obtained by a
boost ALONG the shock surface (i.e. perpendicular to the shock normal `n_hat`),
in which the upstream bulk flow `u` has NO tangential component — the flow is
purely along the normal (normal incidence).

Convention: the boost velocity is the tangential part of `u`,

    V_NIF = u − (u·n̂) n̂              (n̂ normalized internally),

so the flow seen in the NIF is `u − V_NIF = (u·n̂) n̂`, which is parallel to `n̂`
(tangential component ≈ 0). `n_hat` need not be unit length; a zero `n_hat`
returns the zero vector (frame undefined).
"""
function normal_incidence_frame(u::NTuple{3,<:Real}, B::NTuple{3,<:Real}, n_hat::NTuple{3,<:Real})
    ux, uy, uz = _require_finite_point3("u", u, Float64)
    _require_finite_point3("B", B, Float64)
    nx, ny, nz = _require_finite_point3("n_hat", n_hat, Float64)
    n2 = nx * nx + ny * ny + nz * nz
    if n2 == 0
        z = zero(float(n2))
        return (z, z, z)
    end
    udotn = (ux * nx + uy * ny + uz * nz) / n2   # (u·n̂)/|n̂| in units of /|n̂|
    # tangential part of u = u − (u·n̂)n̂  (n̂ = n_hat/|n_hat|)
    return (ux - udotn * nx, uy - udotn * ny, uz - udotn * nz)
end

# --- shock-front locator (from diagnostics.jl) ---
"""
    shock_front(Bz, x) -> (x_s, width)

Shock-front position (steepest |∂Bz/∂x|) and ramp width
`(Bz_down − Bz_up) / max|∂Bz/∂x|`.
"""
function shock_front(Bz::AbstractVector{T}, x::AbstractVector{T}) where {T}
    n = length(Bz)
    n > 0 || throw(ArgumentError("Bz and x must be nonempty"))
    length(x) == n || throw(DimensionMismatch("Bz and x must have the same length"))
    n >= 2 || throw(ArgumentError("Bz and x must contain at least two samples"))
    _require_finite_real_sequence("Bz", Bz)
    _require_finite_real_sequence("x", x)
    gmax = zero(T)
    im = 1
    @inbounds for i = 2:n
        dx = x[i] - x[i-1]
        dx > zero(T) || throw(ArgumentError("x must be strictly increasing"))
        gx = abs(Bz[i] - Bz[i-1]) / dx
        gx > gmax && (gmax = gx; im = i)
    end
    bz_down = Bz[1]
    bz_up = Bz[end]                # wall side vs inflow side
    width = gmax > 0 ? abs(bz_down - bz_up) / gmax : T(NaN)
    return x[im], width
end
