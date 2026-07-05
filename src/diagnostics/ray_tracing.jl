# ray_tracing.jl — WKB / geometric-optics ray tracing for the hybrid wave branches.
#
# A supporting method (diagnostic), not a field solver: wave-packet trajectories
# x(t), k(t) are traced through a FROZEN, smoothly varying medium — either analytic
# n(x), B(x) functions or a snapshot of the simulation's own grid fields — by
# integrating the Hamiltonian ray equations of the warm Hall-MHD scalar dispersion
# relation. Everything is in the code's Ω_ci-normalized units (length d_i, time
# Ω_ci⁻¹, velocity v_A, B in B0, n in n0).
#
# Dispersion relation. With u = ω², the massless-electron quasineutral hybrid
# (warm Hall-MHD) model obeys the cubic
#
#     D(ω, k; n, B) = (u − A)(u² − P·u + Q) − R·u·(u − S) = 0
#
#     A = k∥² v_A²        P = k²(v_A² + c_s²)     Q = k² k∥² v_A² c_s²
#     R = k² k∥² B²/n²    S = k² c_s²
#     v_A² = B²/n         k∥² = (k·B)²/B²
#     c_s²(n) = γ_e T_e n^{γ_e−1} + γ_i T_i n^{γ_i−1}   (combined sound speed, HYB-005)
#
# whose three non-negative roots u are the slow, intermediate (shear-Alfvén /
# ion-cyclotron), and fast (magnetosonic / whistler) branches at arbitrary angle,
# Hall term included. RAY-001 verifies these roots against the independent HYB-006
# eigenvalue oracle (test/oracles/hybrid_dispersion_oracle.jl) and the closed-form
# parallel R/L modes ω = ±k²/2 + k√(1+k²/4) and perpendicular fast mode.
#
# Sound-speed parameters map from the electron closures: IsothermalElectrons(Te) →
# (Te, γe=1); PolytropicElectrons(pe0, n0, γ) → (Te=pe0/n0^γ, γe=γ); warm ions enter
# through the optional (Ti, γi). CGLElectrons is gyrotropic and has no scalar c_s,
# so it is not supported here.
#
# Ray equations. For a time-stationary medium ω is constant along a ray and
#
#     dx/dt = −(∂D/∂k)/(∂D/∂ω) = v_g        dk/dt = +(∂D/∂x)/(∂D/∂ω)
#
# ∂D/∂ω, ∂D/∂k, ∂D/∂n, ∂D/∂B are computed by COMPLEX-STEP differentiation of the
# single generic implementation of D (machine-precision; D is rational plus
# n^(γ−1), complex-analytic for Re n > 0), and ∂D/∂x follows by the chain rule
# through the medium gradients. Integration is classical fixed-step RK4; the
# relative residual |D|/scale is recorded every step as the standard ray-tracing
# accuracy diagnostic (it is conserved exactly by the continuous flow, so growth
# measures integration error).
#
# Rays always carry 3-component positions and wavevectors; a GridRayMedium{D}
# varies only along its D grid axes (the package's "positions carry D coordinates,
# velocities carry 3" convention, applied to rays). Media are periodic in the grid
# axes (CIC gather with wraparound, node i at (i−1)·dx exactly as deposit.jl);
# reported ray positions are left unwrapped so trajectories stay continuous.

# ---------------------------------------------------------------- dispersion core

const _RAY_BRANCHES = (:slow, :intermediate, :fast)

@inline _ray_branch_index(b::Symbol) =
    b === :slow ? 1 :
    b === :intermediate ? 2 :
    b === :fast ? 3 : throw(ArgumentError("branch must be :slow, :intermediate, or :fast (got $b)"))

# Generic (complex-capable) coefficient core. All five coefficients are
# non-negative for real physical inputs (n > 0, Te,Ti ≥ 0), which makes the cubic
# roots real and non-negative and gives a positive residual scale.
@inline function _hybrid_apqrs(kx, ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)
    b2 = Bx * Bx + By * By + Bz * Bz
    k2 = kx * kx + ky * ky + kz * kz
    kdB = kx * Bx + ky * By + kz * Bz
    kpar2 = (kdB * kdB) / b2
    va2 = b2 / n
    cs2 = γe * Te * n^(γe - 1) + γi * Ti * n^(γi - 1)
    A = kpar2 * va2
    P = k2 * (va2 + cs2)
    Q = k2 * kpar2 * va2 * cs2
    R = k2 * kpar2 * b2 / (n * n)
    S = k2 * cs2
    return A, P, Q, R, S
end

@inline function _hybrid_D(ω, kx, ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)
    A, P, Q, R, S = _hybrid_apqrs(kx, ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)
    u = ω * ω
    return (u - A) * (u * u - P * u + Q) - R * u * (u - S)
end

function _validated_cs_params(Te::Real, γe::Real, Ti::Real, γi::Real)
    return (
        _require_finite_nonnegative_real("Te", Te, Float64),
        _require_finite_positive_real("γe", γe, Float64),
        _require_finite_nonnegative_real("Ti", Ti, Float64),
        _require_finite_positive_real("γi", γi, Float64),
    )
end

function _validated_wave_point(k, n::Real, B)
    kk = _require_finite_point3("k", k, Float64)
    nn = _require_finite_positive_real("n", n, Float64)
    BB = _require_finite_point3("B", B, Float64)
    (BB[1]^2 + BB[2]^2 + BB[3]^2) > 0.0 ||
        throw(ArgumentError("B must be nonzero (k∥ is undefined in an unmagnetized plasma)"))
    return kk, nn, BB
end

"""
    hybrid_wave_dispersion(ω, k, n, B; Te=0.0, γe=1.0, Ti=0.0, γi=5/3) -> Float64

Scalar warm Hall-MHD dispersion function `D(ω, k; n, B)` (header of this file) in
Ω_ci-normalized units; `D = 0` on the slow/intermediate/fast branches. `k` and `B`
are 3-tuples, `n > 0`, and the sound speed is `c_s² = γe·Te·n^(γe−1) + γi·Ti·n^(γi−1)`.
"""
function hybrid_wave_dispersion(
    ω::Real,
    k::NTuple{3,<:Real},
    n::Real,
    B::NTuple{3,<:Real};
    Te::Real = 0.0,
    γe::Real = 1.0,
    Ti::Real = 0.0,
    γi::Real = 5 / 3,
)
    ωv = _require_finite_real("ω", ω, Float64)
    kk, nn, BB = _validated_wave_point(k, n, B)
    TeV, γeV, TiV, γiV = _validated_cs_params(Te, γe, Ti, γi)
    return _hybrid_D(ωv, kk..., nn, BB..., TeV, γeV, TiV, γiV)
end

# Real roots of u³ + c2·u² + c1·u + c0 via the companion matrix. The physical
# dispersion cubic has three real non-negative roots; a relative imaginary part
# beyond 1e-6 (far above double-root rounding ~√eps) means non-physical input.
function _cubic_roots_real(c2::Float64, c1::Float64, c0::Float64)
    C = zeros(3, 3)
    C[2, 1] = 1.0
    C[3, 2] = 1.0
    C[1, 3] = -c0
    C[2, 3] = -c1
    C[3, 3] = -c2
    z = eigvals(C)
    u1, u2, u3 = ntuple(3) do i
        zi = z[i]
        abs(imag(zi)) <= 1.0e-6 * max(1.0, abs(zi)) ||
            error("dispersion cubic produced a non-real root ($(zi)); non-physical parameters")
        max(real(zi), 0.0)
    end
    lo, mid, hi = min(u1, u2, u3), u1 + u2 + u3 - min(u1, u2, u3) - max(u1, u2, u3), max(u1, u2, u3)
    return lo, mid, hi
end

"""
    hybrid_wave_frequencies(k, n, B; Te=0.0, γe=1.0, Ti=0.0, γi=5/3) -> NTuple{3,Float64}

The three non-negative wave frequencies `ω ≥ 0` of the warm Hall-MHD dispersion
relation at wavevector `k` in a plasma with density `n` and magnetic field `B`
(Ω_ci-normalized units), sorted ascending: `(ω_slow, ω_intermediate, ω_fast)`.
Branch labels are by local frequency ordering (the standard ray-tracing caveat:
across parameter space the ordered roots may swap physical identity at branch
degeneracies). Verified against the HYB-006 eigenvalue oracle (RAY-001).
"""
function hybrid_wave_frequencies(
    k::NTuple{3,<:Real},
    n::Real,
    B::NTuple{3,<:Real};
    Te::Real = 0.0,
    γe::Real = 1.0,
    Ti::Real = 0.0,
    γi::Real = 5 / 3,
)
    kk, nn, BB = _validated_wave_point(k, n, B)
    TeV, γeV, TiV, γiV = _validated_cs_params(Te, γe, Ti, γi)
    A, P, Q, R, S = _hybrid_apqrs(kk..., nn, BB..., TeV, γeV, TiV, γiV)
    lo, mid, hi = _cubic_roots_real(-(A + P + R), Q + A * P + R * S, -A * Q)
    return (sqrt(lo), sqrt(mid), sqrt(hi))
end

"""
    hybrid_wavenumbers(ω, khat, n, B; Te=0.0, γe=1.0, Ti=0.0, γi=5/3)
        -> Vector{@NamedTuple{kmag, k, branch}}

All positive wavenumber magnitudes `kmag` solving `D(ω, kmag·k̂) = 0` along the
direction `khat` at prescribed frequency `ω > 0` — the standard way to launch rays
at a given frequency (e.g. from an antenna). Substituting `k = κ·k̂` makes `D` a
cubic in `w = κ²`; each positive real root is returned with the propagating
wavevector `k = kmag·k̂/|k̂|` and its branch label (`:slow`/`:intermediate`/`:fast`,
by verifying `hybrid_wave_frequencies` at that `k` reproduces `ω`). Sorted by
`kmag` ascending. Directions/frequencies with no propagating solution (evanescent)
return an empty vector.
"""
function hybrid_wavenumbers(
    ω::Real,
    khat::NTuple{3,<:Real},
    n::Real,
    B::NTuple{3,<:Real};
    Te::Real = 0.0,
    γe::Real = 1.0,
    Ti::Real = 0.0,
    γi::Real = 5 / 3,
)
    ωv = _require_finite_positive_real("ω", ω, Float64)
    kh = _require_finite_point3("khat", khat, Float64)
    khn = sqrt(kh[1]^2 + kh[2]^2 + kh[3]^2)
    khn > 0.0 || throw(ArgumentError("khat must be a nonzero direction"))
    u = (kh[1] / khn, kh[2] / khn, kh[3] / khn)
    _, nn, BB = _validated_wave_point(u, n, B)
    TeV, γeV, TiV, γiV = _validated_cs_params(Te, γe, Ti, γi)
    a, p, q, r, s = _hybrid_apqrs(u..., nn, BB..., TeV, γeV, TiV, γiV)
    uu = ωv * ωv
    # D(w) = (rs·u − aq)·w³ + u(q + ap − ru)·w² − (p+a)u²·w + u³,  w = κ²
    coeffs = (r * s * uu - a * q, uu * (q + a * p - r * uu), -(p + a) * uu * uu, uu^3)
    scale = max(abs(coeffs[1]), abs(coeffs[2]), abs(coeffs[3]), abs(coeffs[4]))
    lead = 1
    # trim (near-)vanishing leading coefficients: exactly zero at symmetric
    # geometries (e.g. perpendicular propagation), else the corresponding root is
    # far outside WKB validity
    while lead <= 3 && abs(coeffs[lead]) <= 1.0e-14 * scale
        lead += 1
    end
    deg = 4 - lead
    out = @NamedTuple{kmag::Float64, k::NTuple{3,Float64}, branch::Symbol}[]
    deg == 0 && return out                    # constant u³ > 0: no propagating root
    C = zeros(deg, deg)
    for i = 1:deg-1
        C[i+1, i] = 1.0
    end
    for i = 1:deg
        C[i, deg] = -coeffs[4-i+1] / coeffs[lead]   # -c_{i-1}/c_deg (ascending index)
    end
    for z in eigvals(C)
        abs(imag(z)) <= 1.0e-6 * max(1.0, abs(z)) || continue
        w = real(z)
        w > 0.0 || continue
        κ = sqrt(w)
        kvec = (κ * u[1], κ * u[2], κ * u[3])
        freqs = hybrid_wave_frequencies(kvec, nn, BB; Te = TeV, γe = γeV, Ti = TiV, γi = γiV)
        db = (abs(freqs[1] - ωv), abs(freqs[2] - ωv), abs(freqs[3] - ωv))
        bi = db[1] <= db[2] && db[1] <= db[3] ? 1 : (db[2] <= db[3] ? 2 : 3)
        # roots that fail to reproduce ω are numerical artifacts of the companion
        # solve near degeneracies; keep only verified propagating solutions
        db[bi] <= 1.0e-6 * ωv || continue
        push!(out, (; kmag = κ, k = kvec, branch = _RAY_BRANCHES[bi]))
    end
    sort!(out; by = t -> t.kmag)
    return out
end

# ---------------------------------------------------------------- derivatives

# Complex-step differentiation of _hybrid_D: exact to machine precision for the
# rational-plus-power dispersion function, no step-size cancellation (the step
# only needs to avoid underflow, so a relative 1e-20 is safe).
@inline _cstep_h(a::Float64) = 1.0e-20 * max(1.0, abs(a))

function _dispersion_derivs(
    ω::Float64,
    k::NTuple{3,Float64},
    n::Float64,
    B::NTuple{3,Float64},
    Te::Float64,
    γe::Float64,
    Ti::Float64,
    γi::Float64,
)
    kx, ky, kz = k
    Bx, By, Bz = B
    h = _cstep_h(ω)
    Dω = imag(_hybrid_D(complex(ω, h), kx, ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(kx)
    Dk1 = imag(_hybrid_D(ω, complex(kx, h), ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(ky)
    Dk2 = imag(_hybrid_D(ω, kx, complex(ky, h), kz, n, Bx, By, Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(kz)
    Dk3 = imag(_hybrid_D(ω, kx, ky, complex(kz, h), n, Bx, By, Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(n)
    Dn = imag(_hybrid_D(ω, kx, ky, kz, complex(n, h), Bx, By, Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(Bx)
    DB1 = imag(_hybrid_D(ω, kx, ky, kz, n, complex(Bx, h), By, Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(By)
    DB2 = imag(_hybrid_D(ω, kx, ky, kz, n, Bx, complex(By, h), Bz, Te, γe, Ti, γi)) / h
    h = _cstep_h(Bz)
    DB3 = imag(_hybrid_D(ω, kx, ky, kz, n, Bx, By, complex(Bz, h), Te, γe, Ti, γi)) / h
    Dval = _hybrid_D(ω, kx, ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)
    # positive residual scale: same monomials with all signs +, u = ω² > 0
    A, P, Q, R, S = _hybrid_apqrs(kx, ky, kz, n, Bx, By, Bz, Te, γe, Ti, γi)
    u = ω * ω
    scale = u^3 + (A + P + R) * u^2 + (Q + A * P + R * S) * u + A * Q
    residual = abs(Dval) / scale
    return (; Dω, Dk = (Dk1, Dk2, Dk3), Dn, DB = (DB1, DB2, DB3), residual)
end

"""
    wave_group_velocity(k, n, B; branch=:fast, Te=0.0, γe=1.0, Ti=0.0, γi=5/3)
        -> (; vg, ω)

Group velocity `v_g = ∂ω/∂k = −(∂D/∂k)/(∂D/∂ω)` (3-tuple, units of v_A) and
frequency of the selected `branch` (`:slow`, `:intermediate`, `:fast`) at a single
plasma point. Throws for degenerate points where the branch frequency vanishes or
`∂D/∂ω = 0`. Verified against the oracle's numerical `dω/dk` (RAY-003).
"""
function wave_group_velocity(
    k::NTuple{3,<:Real},
    n::Real,
    B::NTuple{3,<:Real};
    branch::Symbol = :fast,
    Te::Real = 0.0,
    γe::Real = 1.0,
    Ti::Real = 0.0,
    γi::Real = 5 / 3,
)
    bi = _ray_branch_index(branch)
    kk, nn, BB = _validated_wave_point(k, n, B)
    TeV, γeV, TiV, γiV = _validated_cs_params(Te, γe, Ti, γi)
    ω = hybrid_wave_frequencies(kk, nn, BB; Te = TeV, γe = γeV, Ti = TiV, γi = γiV)[bi]
    ω > 0.0 || throw(ArgumentError("branch $branch is degenerate (ω = 0) at this point"))
    d = _dispersion_derivs(ω, kk, nn, BB, TeV, γeV, TiV, γiV)
    d.Dω != 0.0 || throw(ArgumentError("∂D/∂ω = 0 at this point (branch degeneracy)"))
    vg = (-d.Dk[1] / d.Dω, -d.Dk[2] / d.Dω, -d.Dk[3] / d.Dω)
    all(isfinite, vg) || throw(ArgumentError("group velocity is not finite at this point"))
    return (; vg, ω)
end

# ---------------------------------------------------------------- media

"""
    RayMedium

A frozen, smoothly varying plasma background `(n(x), B(x))` plus sound-speed
parameters `(Te, γe, Ti, γi)` that rays are traced through. Concrete media:
[`AnalyticRayMedium`](@ref), [`GridRayMedium`](@ref).
"""
abstract type RayMedium end

"""
    AnalyticRayMedium(nfun, Bfun; Te=0.0, γe=1.0, Ti=0.0, γi=5/3, h=cbrt(eps()))

Analytic ray-tracing medium: `nfun(x, y, z) -> n` (density, > 0) and
`Bfun(x, y, z) -> (Bx, By, Bz)`, in normalized units. Medium gradients (needed for
`dk/dt`) are computed by central finite differences with per-axis step
`h·max(1, |xᵢ|)`; the default `h = cbrt(eps())` balances truncation against
round-off for smooth functions. The functions must be evaluable in an `h`-
neighborhood of every ray point.
"""
struct AnalyticRayMedium{FN,FB} <: RayMedium
    nfun::FN
    Bfun::FB
    Te::Float64
    γe::Float64
    Ti::Float64
    γi::Float64
    h::Float64
end

function AnalyticRayMedium(
    nfun,
    Bfun;
    Te::Real = 0.0,
    γe::Real = 1.0,
    Ti::Real = 0.0,
    γi::Real = 5 / 3,
    h::Real = cbrt(eps(Float64)),
)
    TeV, γeV, TiV, γiV = _validated_cs_params(Te, γe, Ti, γi)
    hV = _require_finite_positive_real("h", h, Float64)
    return AnalyticRayMedium{typeof(nfun),typeof(Bfun)}(nfun, Bfun, TeV, γeV, TiV, γiV, hV)
end

@inline _b3(Bv) = (Float64(Bv[1]), Float64(Bv[2]), Float64(Bv[3]))

# central-difference gradients of (n, B) along axis i at position pos
function _fd_axis(med::AnalyticRayMedium, pos::NTuple{3,Float64}, i::Int)
    hi = med.h * max(1.0, abs(pos[i]))
    pp = ntuple(j -> j == i ? pos[j] + hi : pos[j], 3)
    pm = ntuple(j -> j == i ? pos[j] - hi : pos[j], 3)
    np = Float64(med.nfun(pp...))
    nm = Float64(med.nfun(pm...))
    Bp = _b3(med.Bfun(pp...))
    Bm = _b3(med.Bfun(pm...))
    inv2h = 1.0 / (pp[i] - pm[i])                # exact spacing after rounding
    return (
        (np - nm) * inv2h,
        (Bp[1] - Bm[1]) * inv2h,
        (Bp[2] - Bm[2]) * inv2h,
        (Bp[3] - Bm[3]) * inv2h,
    )
end

function _medium_at(med::AnalyticRayMedium, x::Float64, y::Float64, z::Float64)
    n = Float64(med.nfun(x, y, z))
    B = _b3(med.Bfun(x, y, z))
    pos = (x, y, z)
    d1 = _fd_axis(med, pos, 1)
    d2 = _fd_axis(med, pos, 2)
    d3 = _fd_axis(med, pos, 3)
    gn = (d1[1], d2[1], d3[1])
    gB1 = (d1[2], d2[2], d3[2])
    gB2 = (d1[3], d2[3], d3[3])
    gB3 = (d1[4], d2[4], d3[4])
    return n, B, gn, (gB1, gB2, gB3)
end

"""
    GridRayMedium(g::FourierGrid{D}, n, B; Te=0.0, γe=1.0, Ti=0.0, γi=5/3)

Ray-tracing medium from a SNAPSHOT of grid fields on the periodic grid `g`
(`D = 1, 2, 3`): density `n::Array{T,D}` (positive everywhere) and magnetic field
`B = (Bx, By, Bz)` of `Array{T,D}`. The arrays are copied (later mutation of the
simulation fields does not affect the medium) and the medium gradients are
precomputed SPECTRALLY at construction; values and gradients are then
CIC-interpolated at ray positions with periodic wraparound (node `i` at
`(i−1)·dx`, the `deposit.jl`/`gather_at` convention). The medium varies only along
the grid's `D` axes; along missing axes rays propagate freely.

Pass the moment density from [`compute_moments!`](@ref) and the `HybridFields` `B`
components to trace rays through a hybrid snapshot.
"""
struct GridRayMedium{D,T<:AbstractFloat,G<:FourierGrid{D,T}} <: RayMedium
    g::G
    n::Array{T,D}
    B::NTuple{3,Array{T,D}}
    gradn::NTuple{D,Array{T,D}}
    gradB::NTuple{3,NTuple{D,Array{T,D}}}   # gradB[j][i] = ∂B_j/∂x_i
    Te::Float64
    γe::Float64
    Ti::Float64
    γi::Float64
end

function GridRayMedium(
    g::FourierGrid{D,T},
    n::Array{T,D},
    B::NTuple{3,Array{T,D}};
    Te::Real = 0.0,
    γe::Real = 1.0,
    Ti::Real = 0.0,
    γi::Real = 5 / 3,
) where {D,T<:AbstractFloat}
    size(n) == g.n || throw(DimensionMismatch("density size $(size(n)) ≠ grid size $(g.n)"))
    for j = 1:3
        size(B[j]) == g.n || throw(DimensionMismatch("B[$j] size $(size(B[j])) ≠ grid size $(g.n)"))
    end
    _require_all_finite("density", n, "GridRayMedium")
    minimum(n) > zero(T) ||
        throw(ArgumentError("density must be positive everywhere for ray tracing"))
    for j = 1:3
        _require_all_finite("B[$j]", B[j], "GridRayMedium")
    end
    TeV, γeV, TiV, γiV = _validated_cs_params(Te, γe, Ti, γi)
    nc = copy(n)
    Bc = (copy(B[1]), copy(B[2]), copy(B[3]))
    gradn = ntuple(i -> deriv(nc, g, i), D)
    gradB = ntuple(j -> ntuple(i -> deriv(Bc[j], g, i), D), 3)
    return GridRayMedium{D,T,typeof(g)}(g, nc, Bc, gradn, gradB, TeV, γeV, TiV, γiV)
end

# Per-axis CIC weights with periodic wrap; guards |x/dx| against Int overflow so a
# runaway (but still finite) ray reports NaN instead of throwing InexactError.
@inline function _cic_axis(x::Float64, dx::Float64, npts::Int)
    s = x / dx
    abs(s) <= 9.0e15 || return (1, 1, NaN)
    i0 = floor(Int, s)
    return (mod(i0, npts) + 1, mod(i0 + 1, npts) + 1, s - i0)
end

@inline function _cic(f::Array{T,1}, g::FourierGrid{1,T}, xs::NTuple{1,Float64}) where {T}
    ia, ib, w = _cic_axis(xs[1], Float64(g.dx[1]), g.n[1])
    isnan(w) && return NaN
    @inbounds return (1.0 - w) * Float64(f[ia]) + w * Float64(f[ib])
end

@inline function _cic(f::Array{T,2}, g::FourierGrid{2,T}, xs::NTuple{2,Float64}) where {T}
    ia, ib, wx = _cic_axis(xs[1], Float64(g.dx[1]), g.n[1])
    ja, jb, wy = _cic_axis(xs[2], Float64(g.dx[2]), g.n[2])
    (isnan(wx) || isnan(wy)) && return NaN
    @inbounds return (1.0 - wx) * ((1.0 - wy) * Float64(f[ia, ja]) + wy * Float64(f[ia, jb])) +
                     wx * ((1.0 - wy) * Float64(f[ib, ja]) + wy * Float64(f[ib, jb]))
end

@inline function _cic(f::Array{T,3}, g::FourierGrid{3,T}, xs::NTuple{3,Float64}) where {T}
    ia, ib, wx = _cic_axis(xs[1], Float64(g.dx[1]), g.n[1])
    ja, jb, wy = _cic_axis(xs[2], Float64(g.dx[2]), g.n[2])
    ka, kb, wz = _cic_axis(xs[3], Float64(g.dx[3]), g.n[3])
    (isnan(wx) || isnan(wy) || isnan(wz)) && return NaN
    @inbounds begin
        c00 = (1.0 - wx) * Float64(f[ia, ja, ka]) + wx * Float64(f[ib, ja, ka])
        c10 = (1.0 - wx) * Float64(f[ia, jb, ka]) + wx * Float64(f[ib, jb, ka])
        c01 = (1.0 - wx) * Float64(f[ia, ja, kb]) + wx * Float64(f[ib, ja, kb])
        c11 = (1.0 - wx) * Float64(f[ia, jb, kb]) + wx * Float64(f[ib, jb, kb])
    end
    return (1.0 - wz) * ((1.0 - wy) * c00 + wy * c10) + wz * ((1.0 - wy) * c01 + wy * c11)
end

function _medium_at(med::GridRayMedium{D}, x::Float64, y::Float64, z::Float64) where {D}
    xs = ntuple(i -> (x, y, z)[i], D)
    n = _cic(med.n, med.g, xs)
    B = (_cic(med.B[1], med.g, xs), _cic(med.B[2], med.g, xs), _cic(med.B[3], med.g, xs))
    gn = ntuple(i -> i <= D ? _cic(med.gradn[i], med.g, xs) : 0.0, 3)
    gB = ntuple(3) do j
        ntuple(i -> i <= D ? _cic(med.gradB[j][i], med.g, xs) : 0.0, 3)
    end
    return n, B, gn, gB
end

# ---------------------------------------------------------------- ray integration

@inline _cs_kwargs(med::RayMedium) = (med.Te, med.γe, med.Ti, med.γi)

# Ray-equation right-hand side at (x, k) for fixed ω. Returns status ≠ :ok instead
# of throwing so RK4 stages can terminate a trace cleanly mid-flight.
function _ray_rhs(
    med::RayMedium,
    x::NTuple{3,Float64},
    k::NTuple{3,Float64},
    ω::Float64,
    Bmin::Float64,
)
    z3 = (0.0, 0.0, 0.0)
    (all(isfinite, x) && all(isfinite, k)) ||
        return (; status = :nonfinite, vg = z3, kdot = z3, residual = NaN)
    n, B, gn, gB = _medium_at(med, x[1], x[2], x[3])
    (isfinite(n) && n > 0.0 && all(isfinite, B)) ||
        return (; status = :invalid_medium, vg = z3, kdot = z3, residual = NaN)
    (B[1]^2 + B[2]^2 + B[3]^2) > Bmin * Bmin ||
        return (; status = :invalid_medium, vg = z3, kdot = z3, residual = NaN)
    (all(isfinite, gn) && all(isfinite, gB[1]) && all(isfinite, gB[2]) && all(isfinite, gB[3])) ||
        return (; status = :invalid_medium, vg = z3, kdot = z3, residual = NaN)
    Te, γe, Ti, γi = _cs_kwargs(med)
    d = _dispersion_derivs(ω, k, n, B, Te, γe, Ti, γi)
    d.Dω == 0.0 && return (; status = :caustic, vg = z3, kdot = z3, residual = d.residual)
    vg = (-d.Dk[1] / d.Dω, -d.Dk[2] / d.Dω, -d.Dk[3] / d.Dω)
    kdot = ntuple(3) do i
        (d.Dn * gn[i] + d.DB[1] * gB[1][i] + d.DB[2] * gB[2][i] + d.DB[3] * gB[3][i]) / d.Dω
    end
    (all(isfinite, vg) && all(isfinite, kdot)) ||
        return (; status = :caustic, vg = z3, kdot = z3, residual = d.residual)
    return (; status = :ok, vg, kdot, residual = d.residual)
end

"""
    trace_ray(med::RayMedium, x0, k0; branch=:fast, dt, nsteps,
              Bmin=1e-8, residual_max=1e-2)
        -> (; t, x, k, vg, residual, ω, branch, status)

Trace one wave-packet ray from position `x0` with initial wavevector `k0` (both
3-tuples, normalized units) through the frozen medium `med`, on the dispersion
branch selected at launch (`:slow`, `:intermediate`, `:fast`). The launch
frequency `ω` is solved from `k0` and then held exactly constant (time-stationary
medium); `(x, k)` are advanced by classical RK4 with fixed step `dt` for `nsteps`
steps. To trace backward along a ray, launch with `-k0` (the relation is even in
`k`).

Returns trajectory matrices `x`, `k`, `vg` (each `3 × m`, columns = saved points,
`m ≤ nsteps + 1`), times `t`, and the relative dispersion residual
`|D(x, k, ω)| / scale` per point — the standard accuracy diagnostic (exactly
conserved by the continuous ray flow; growth measures RK4 error, so shrink `dt`
when it rises). `status` is `:ok`, or the reason the trace stopped early:
`:nonfinite` (state left the reals), `:invalid_medium` (`n ≤ 0`, `|B| ≤ Bmin`, or
non-finite medium), `:caustic` (`∂D/∂ω → 0`; branch degeneracy / mode-conversion
point, WKB breaks down), or `:residual` (residual exceeded `residual_max`). When a
final diagnostic evaluation fails, the last saved `vg`/`residual` column is NaN.
"""
function trace_ray(
    med::RayMedium,
    x0::NTuple{3,<:Real},
    k0::NTuple{3,<:Real};
    branch::Symbol = :fast,
    dt::Real,
    nsteps::Integer,
    Bmin::Real = 1.0e-8,
    residual_max::Real = 1.0e-2,
)
    bi = _ray_branch_index(branch)
    x = _require_finite_point3("x0", x0, Float64)
    k = _require_finite_point3("k0", k0, Float64)
    (k[1]^2 + k[2]^2 + k[3]^2) > 0.0 || throw(ArgumentError("k0 must be nonzero"))
    dtv = _require_finite_positive_real("dt", dt, Float64)
    ns = _require_positive_intlike("nsteps", nsteps)
    BminV = _require_finite_nonnegative_real("Bmin", Bmin, Float64)
    resmax = _require_finite_positive_real("residual_max", residual_max, Float64)

    n0, B0, _, _ = _medium_at(med, x[1], x[2], x[3])
    (isfinite(n0) && n0 > 0.0 && all(isfinite, B0)) ||
        throw(ArgumentError("medium is invalid at the launch point (n = $n0, B = $B0)"))
    (B0[1]^2 + B0[2]^2 + B0[3]^2) > BminV * BminV ||
        throw(ArgumentError("|B| ≤ Bmin at the launch point; rays require a magnetized medium"))
    Te, γe, Ti, γi = _cs_kwargs(med)
    ω = hybrid_wave_frequencies(k, n0, B0; Te, γe, Ti, γi)[bi]
    ω > 0.0 || throw(
        ArgumentError("branch $branch is degenerate (ω = 0) at the launch point; not traceable"),
    )

    xs = Matrix{Float64}(undef, 3, ns + 1)
    ks = Matrix{Float64}(undef, 3, ns + 1)
    vgs = Matrix{Float64}(undef, 3, ns + 1)
    ts = Vector{Float64}(undef, ns + 1)
    res = Vector{Float64}(undef, ns + 1)

    r1 = _ray_rhs(med, x, k, ω, BminV)
    r1.status === :ok ||
        throw(ArgumentError("ray equations are singular at the launch point ($(r1.status))"))
    m = 1
    ts[1] = 0.0
    xs[:, 1] .= x
    ks[:, 1] .= k
    vgs[:, 1] .= r1.vg
    res[1] = r1.residual
    status = :ok
    for step = 1:ns
        h2 = dtv / 2
        x2 = x .+ h2 .* r1.vg
        k2 = k .+ h2 .* r1.kdot
        r2 = _ray_rhs(med, x2, k2, ω, BminV)
        r2.status === :ok || (status = r2.status; break)
        x3 = x .+ h2 .* r2.vg
        k3 = k .+ h2 .* r2.kdot
        r3 = _ray_rhs(med, x3, k3, ω, BminV)
        r3.status === :ok || (status = r3.status; break)
        x4 = x .+ dtv .* r3.vg
        k4 = k .+ dtv .* r3.kdot
        r4 = _ray_rhs(med, x4, k4, ω, BminV)
        r4.status === :ok || (status = r4.status; break)
        x = x .+ (dtv / 6) .* (r1.vg .+ 2 .* r2.vg .+ 2 .* r3.vg .+ r4.vg)
        k = k .+ (dtv / 6) .* (r1.kdot .+ 2 .* r2.kdot .+ 2 .* r3.kdot .+ r4.kdot)
        rn = _ray_rhs(med, x, k, ω, BminV)
        m += 1
        ts[m] = step * dtv
        xs[:, m] .= x
        ks[:, m] .= k
        if rn.status === :ok
            vgs[:, m] .= rn.vg
            res[m] = rn.residual
        else
            vgs[:, m] .= NaN
            res[m] = NaN
            status = rn.status
            break
        end
        if rn.residual > resmax
            status = :residual
            break
        end
        r1 = rn
    end
    return (;
        t = ts[1:m],
        x = xs[:, 1:m],
        k = ks[:, 1:m],
        vg = vgs[:, 1:m],
        residual = res[1:m],
        ω,
        branch,
        status,
    )
end
