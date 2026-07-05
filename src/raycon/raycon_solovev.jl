# raycon_solovev.jl — Solovev equilibrium flux label and derivatives
# (port of solovev.m) plus flux-coordinate mapping (mapFlux.m, calcFlux.m).
#
# Coordinates (upstream convention, used everywhere):
#   r = ρ·cosθ + r0,   z = −ρ·sinθ,
#   ρ = √((r−r0)² + z²),   θ = atan2(r0−r, −z) + π/2.
# Flux:  ψ = fac·(r²z²/E² + ¼(r²−r0²)²),  fac = psin/(iaspr·r0²)²,
#        s = √(ψ/psin) − sflxa   (sflxa ≠ 0 only for root finding).

@inline _sol_coords(eq::SolovevEquilibrium, rho::Float64, theta::Float64) =
    (rho * cos(theta) + eq.r0, -rho * sin(theta))

"""
    solovev_flux(eq, rho, theta; sflxa=0.0, order=0)

Solovev flux label `s = √(ψ/ψn) − sflxa` at minor radius `rho`, poloidal angle
`theta`. `order=0` returns just `s`; `order=1` adds `dsdr, dsdz`; `order=2`
adds the second derivatives; `order=3` adds the third derivatives (the full
`solovev.m` 10-output form). Returns a NamedTuple.
"""
function solovev_flux(
    eq::SolovevEquilibrium,
    rho::Real,
    theta::Real;
    sflxa::Real = 0.0,
    order::Integer = 0,
)
    (isfinite(rho) && rho >= 0) || throw(ArgumentError("rho must be finite and ≥ 0"))
    isfinite(theta) || throw(ArgumentError("theta must be finite"))
    0 <= order <= 3 || throw(ArgumentError("order must be 0, 1, 2 or 3"))
    r, z = _sol_coords(eq, Float64(rho), Float64(theta))
    rsq = r^2
    zsq = z^2
    r0sq = eq.r0^2
    Esq = eq.elong^2
    psin = eq.psin
    fac = psin / (eq.iaspr * r0sq)^2
    psi = fac * (rsq * zsq / Esq + 0.25 * (rsq - r0sq)^2)
    sflx = sqrt(psi / psin) - Float64(sflxa)
    order == 0 && return (; sflx)

    # NB: derivative formulas use the UN-SHIFTED s (upstream always calls with
    # sflxa=0 when derivatives are requested); guard the axis singularity.
    sflxa == 0 || throw(ArgumentError("derivatives require sflxa = 0 (upstream convention)"))
    sflx > 0 ||
        throw(ArgumentError("solovev derivatives are singular on the magnetic axis (s = 0)"))
    dpdr = fac * r * (2 * zsq / Esq + (rsq - r0sq))
    dpdz = fac * z * (2 * rsq / Esq)
    dpdrsq = dpdr^2
    dpdzsq = dpdz^2
    dpds = 2 * sflx * psin
    dsdp = 1 / dpds
    dsdr = dsdp * dpdr
    dsdz = dsdp * dpdz
    order == 1 && return (; sflx, dsdr, dsdz)

    dpdr2 = fac * (2 * zsq / Esq + 3 * rsq - r0sq)
    dpdrz = fac * (4 * r * z / Esq)
    dpdz2 = fac * (2 * rsq / Esq)
    dsdr2 = dsdp * (dpdr2 - dsdp / sflx * dpdrsq)
    dsdrz = dsdp * (dpdrz - dsdp / sflx * dpdr * dpdz)
    dsdz2 = dsdp * (dpdz2 - dsdp / sflx * dpdzsq)
    order == 2 && return (; sflx, dsdr, dsdz, dsdr2, dsdrz, dsdz2)

    dpdr3 = fac * 6 * r
    dpdr2z = fac * 4 * z / Esq
    dpdrz2 = fac * 4 * r / Esq
    dpdz3 = 0.0
    dLNsdp = dsdp / sflx
    dsdrp = -dLNsdp * dsdp * dpdr
    dLNsdrp = dLNsdp * (dsdrp / dsdp - dsdr / sflx)
    dsdzp = -dLNsdp * dsdp * dpdz
    dLNsdzp = dLNsdp * (dsdzp / dsdp - dsdz / sflx)
    dsdr3 = dsdr2 * (dsdrp / dsdp) + dsdp * (dpdr3 - dLNsdrp * dpdrsq - 2 * dLNsdp * dpdr * dpdr2)
    dsdr2z = dsdr2 * (dsdzp / dsdp) + dsdp * (dpdr2z - dLNsdzp * dpdrsq - 2 * dLNsdp * dpdr * dpdrz)
    dsdrz2 = dsdz2 * (dsdrp / dsdp) + dsdp * (dpdrz2 - dLNsdrp * dpdzsq - 2 * dLNsdp * dpdz * dpdrz)
    dsdz3 = dsdz2 * (dsdzp / dsdp) + dsdp * (dpdz3 - dLNsdzp * dpdzsq - 2 * dLNsdp * dpdz * dpdz2)
    return (; sflx, dsdr, dsdz, dsdr2, dsdrz, dsdz2, dsdr3, dsdr2z, dsdrz2, dsdz3)
end

# Brent's method on a bracketing interval (stdlib-only; replaces MATLAB fzero
# with interval input). Tight tolerance is a deliberate upgrade over the
# upstream TolX=1e-4.
function _brent(f, a::Float64, b::Float64; xtol::Float64 = 1e-12, maxit::Int = 200)
    fa = f(a)
    fb = f(b)
    fa == 0 && return a
    fb == 0 && return b
    fa * fb < 0 || throw(ArgumentError("root not bracketed in [$a, $b]"))
    if abs(fa) < abs(fb)
        a, b, fa, fb = b, a, fb, fa
    end
    c, fc = a, fa
    d = b - a
    e = d
    for _ = 1:maxit
        if fb * fc > 0
            c, fc = a, fa
            d = b - a
            e = d
        end
        if abs(fc) < abs(fb)
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        end
        tol = 2 * eps(abs(b)) + xtol
        m = 0.5 * (c - b)
        (abs(m) <= tol || fb == 0) && return b
        if abs(e) < tol || abs(fa) <= abs(fb)
            d = m
            e = m
        else
            s = fb / fa
            if a == c
                p = 2 * m * s
                q = 1 - s
            else
                q = fa / fc
                rr = fb / fc
                p = s * (2 * m * q * (q - rr) - (b - a) * (rr - 1))
                q = (q - 1) * (rr - 1) * (s - 1)
            end
            p > 0 ? (q = -q) : (p = -p)
            if 2p < min(3 * m * q - abs(tol * q), abs(e * q))
                e = d
                d = p / q
            else
                d = m
                e = m
            end
        end
        a, fa = b, fb
        b += abs(d) > tol ? d : (m > 0 ? tol : -tol)
        fb = f(b)
    end
    error("Brent root finding did not converge")
end

# Scalar root find from a starting guess with automatic bracket expansion —
# the port of MATLAB fzero(f, x0) used by inittok/adjust_disp_m/dispertok'Trs'.
function _fzero_near(f, x0::Float64; maxexpand::Int = 60)
    f0 = f(x0)
    f0 == 0 && return x0
    isfinite(f0) || throw(ArgumentError("function not finite at the starting guess"))
    dx = x0 != 0 ? abs(x0) / 50 : 1 / 50   # MATLAB fzero's initial search step
    a = x0
    b = x0
    fa = f0
    fb = f0
    sq2 = sqrt(2.0)
    for _ = 1:maxexpand
        dx *= sq2
        a = x0 - dx
        fa = f(a)
        # hand Brent a bracket whose BOTH endpoints have known finite values:
        # [a, x0] resp. [x0, b] (the far endpoint of the other side may be
        # NaN/Inf when the search left the function domain on that side only)
        (isfinite(fa) && fa * f0 <= 0) && return _brent(f, a, x0)
        b = x0 + dx
        fb = f(b)
        (isfinite(fb) && fb * f0 <= 0) && return _brent(f, x0, b)
        (isfinite(fa) && isfinite(fb)) ||
            throw(ArgumentError("root search left the function domain"))
    end
    throw(ArgumentError("could not bracket a root near $x0"))
end

"""
    map_flux(eq, s, theta) -> (; rho, r, z)

Map flux coordinates `(s, θ)` to `(ρ, R, Z)` by solving `solovev_flux(ρ,θ) = s`
on `ρ ∈ (0, 1.5·ρa]` (port of `mapFlux.m`, with the root polished to 1e-12
instead of the upstream 1e-4).
"""
function map_flux(eq::SolovevEquilibrium, s::Real, theta::Real)
    (isfinite(s) && 0 < s) || throw(ArgumentError("s must be finite and positive"))
    isfinite(theta) || throw(ArgumentError("theta must be finite"))
    th = Float64(theta)
    f = rho -> solovev_flux(eq, rho, th; sflxa = s).sflx
    # ρa = elong·iaspr·r0 is the VERTICAL half-height of the s=1 surface; for
    # oblate shapes (elong < 1) the midplane radius of interior surfaces can
    # exceed 1.5·ρa (upstream's fixed bracket fails there), so expand the upper
    # bound geometrically until the target surface is bracketed
    hi = 1.5 * rho_edge(eq)
    nexpand = 0
    while f(hi) < 0 && nexpand < 20
        hi *= 2
        nexpand += 1
    end
    rho = _brent(f, 1e-12, hi)
    r, z = _sol_coords(eq, rho, th)
    return (; rho, r, z)
end

"""
    flux_surface_mesh(eq; ns=40, nt=45) -> (; s, theta, r, z)

Flux-surface mesh `(R, Z)[ns, nt+1]` over `s ∈ [1e-4, 1]`, `θ ∈ [0, 2π]`
(port of `calcFlux.m`), for deposition binning and plotting.
"""
function flux_surface_mesh(eq::SolovevEquilibrium; ns::Integer = 40, nt::Integer = 45)
    (ns >= 2 && nt >= 4) || throw(ArgumentError("mesh needs ns ≥ 2 and nt ≥ 4"))
    s = [0.0001 + (k - 1) * (0.999 / (ns - 1)) for k = 1:ns]
    theta = [2π * j / nt for j = 0:nt]
    r = Matrix{Float64}(undef, ns, nt + 1)
    z = Matrix{Float64}(undef, ns, nt + 1)
    for k = 1:ns, j = 1:nt+1
        m = map_flux(eq, s[k], theta[j])
        r[k, j] = m.r
        z[k, j] = m.z
    end
    return (; s, theta, r, z)
end
