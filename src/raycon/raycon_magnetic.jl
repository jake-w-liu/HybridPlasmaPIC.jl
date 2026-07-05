# raycon_magnetic.jl — magnetic field amplitude, flux-coordinate metric and
# Stix basis vectors with 1st/2nd derivatives (port of magnetic.m+dmagnetic.m).
#
# Basis (components in cylindrical (r, φ, z)): e_n = ∇s/|∇s| (normal), e_b = B̂
# (parallel), e_p = complement (e_b × e_n orientation as upstream); the wave
# vector is projected as k_{n,b,p} = kr·e_{n,b,p}r + kφ·e_{n,b,p}φ + kz·e_{n,b,p}z.
#
# The (dbdr, dbdz, dbdr2, dbdrz, dbdz2) block reproduces dmagnetic.m: central
# finite differences of b(R, Z) with the upstream-hardcoded step 1e-8 — kept
# (not replaced by analytic derivatives) for 1:1 parity with the original.

# b(R, Z) exactly as dmagnetic.m builds it (from 1st-order solovev data only).
function _bfield_at(eq::SolovevEquilibrium, R::Float64, Z::Float64)
    rho = sqrt((R - eq.r0)^2 + Z^2)
    th = atan(eq.r0 - R, -Z) + π / 2
    s = solovev_flux(eq, rho, th; order = 1)
    dpds = 2 * s.sflx * eq.psin
    h11 = s.dsdr^2 + s.dsdz^2
    gp2 = dpds^2 * h11
    t0 = eq.b0 * eq.r0
    return sqrt((t0 / R)^2 + gp2 / R^2)
end

# dmagnetic.m: FD derivatives of b in (R, Z). Upstream hardcodes step 1e-8,
# which puts ~4·eps·|b|/h² ≈ 10-20% roundoff noise on the SECOND derivatives;
# we use 1e-5 (documented deviation), reducing that noise to ~1e-5 relative
# while keeping truncation error negligible (b varies on the ~0.1 m scale).
function _dmagnetic(eq::SolovevEquilibrium, rho::Float64, theta::Float64)
    r = rho * cos(theta) + eq.r0
    z = -rho * sin(theta)
    dr = 1e-5
    dz = 1e-5
    b00 = _bfield_at(eq, r, z)
    bpr = _bfield_at(eq, r + dr, z)
    bmr = _bfield_at(eq, r - dr, z)
    bpz = _bfield_at(eq, r, z + dz)
    bmz = _bfield_at(eq, r, z - dz)
    bpp = _bfield_at(eq, r + dr, z + dz)
    bpm = _bfield_at(eq, r + dr, z - dz)
    bmp = _bfield_at(eq, r - dr, z + dz)
    bmm = _bfield_at(eq, r - dr, z - dz)
    dbdr = (bpr - bmr) / (2 * dr)
    dbdz = (bpz - bmz) / (2 * dz)
    dbdr2 = (bpr - 2 * b00 + bmr) / dr^2
    dbdz2 = (bpz - 2 * b00 + bmz) / dz^2
    dbdrz = (bpp - bpm - bmp + bmm) / (4 * dr * dz)
    return dbdr, dbdz, dbdr2, dbdrz, dbdz2
end

"""
    magnetic_geometry(eq, rho, theta) -> NamedTuple

Full local magnetic geometry at `(ρ, θ)` (the 76-output form of `magnetic.m`):
field amplitude `b` with flux-coordinate derivatives (`dbds`, `dbdt`, 2nd
order), flux label `sflx` and its (r,z)-derivatives to 3rd order, poloidal-
angle derivatives, the Stix basis vectors `e_n, e_b, e_p` (cylindrical
components `ener…epez`), and their 1st and 2nd derivatives with respect to
`(r, z)` (the field-curvature terms used by the eigenvalue dispersion).
"""
function magnetic_geometry(eq::SolovevEquilibrium, rho::Real, theta::Real)
    (isfinite(rho) && rho > 0) || throw(ArgumentError("rho must be finite and positive"))
    isfinite(theta) || throw(ArgumentError("theta must be finite"))
    rhov = Float64(rho)
    th = Float64(theta)
    cost = cos(th)
    sint = sin(th)
    rho2 = rhov^2
    r = rhov * cost + eq.r0
    s = solovev_flux(eq, rhov, th; order = 3)
    (; sflx, dsdr, dsdz, dsdr2, dsdrz, dsdz2, dsdr3, dsdr2z, dsdrz2, dsdz3) = s
    psin = eq.psin
    dpds = 2 * sflx * psin
    dpds2 = 2 * psin
    t0 = eq.b0 * eq.r0

    # poloidal angle derivatives
    dtdr = -sint / rhov
    dtdz = -cost / rhov
    dtdr2 = 2 * sint * cost / rho2
    dtdrz = (cost^2 - sint^2) / rho2
    dtdz2 = -2 * sint * cost / rho2

    # metric (1st order)
    jac = r / (dsdz * dtdr - dsdr * dtdz)
    drds = -jac / r * dtdz
    dzds = jac / r * dtdr
    drdt = jac / r * dsdz
    dzdt = -jac / r * dsdr
    h11 = dsdr^2 + dsdz^2
    dh11dr = 2 * (dsdr * dsdr2 + dsdz * dsdrz)
    dh11dz = 2 * (dsdr * dsdrz + dsdz * dsdz2)
    dh11ds = dh11dr * drds + dh11dz * dzds
    dh11dt = dh11dr * drdt + dh11dz * dzdt
    gp2 = dpds^2 * h11
    dgp2ds = dpds^2 * (dh11ds + 2 * h11 / sflx)
    dgp2dt = dpds^2 * dh11dt

    # field amplitude (1st derivatives)
    da2ds = dgp2ds
    da2dt = dgp2dt
    r2 = r^2
    bpol2 = gp2 / r2
    btor2 = (t0 / r)^2
    b = sqrt(btor2 + bpol2)
    dbds = da2ds / (2 * b * r2) - b * drds / r
    dbdt = da2dt / (2 * b * r2) - b * drdt / r
    bp = sqrt(bpol2 / btor2)

    # basis vectors
    gs = sqrt(h11)
    ener = dsdr / gs
    enef = 0.0
    enez = dsdz / gs
    eber = t0 * dsdz / (r * b * gs)
    ebef = -dpds * gs / (r * b)
    ebez = -t0 * dsdr / (r * b * gs)
    eper = dpds * dsdz / (r * b)
    epef = t0 / (r * b)
    epez = -dpds * dsdr / (r * b)

    # mixed second derivatives of (s, t) along flux coordinates
    dsdrt = dsdr2 * drdt + dsdrz * dzdt
    dsdzt = dsdrz * drdt + dsdz2 * dzdt
    dtdrs = dtdr2 * drds + dtdrz * dzds
    dtdrt = dtdr2 * drdt + dtdrz * dzdt
    dtdzs = dtdrz * drds + dtdz2 * dzds
    dtdzt = dtdrz * drdt + dtdz2 * dzdt

    # metric (2nd order)
    djacdr = jac / r * (1 - jac * (dsdrz * dtdr + dsdz * dtdr2 - dsdr2 * dtdz - dsdr * dtdrz))
    djacdz = jac / r * (-jac * (dsdz2 * dtdr + dsdz * dtdrz - dsdrz * dtdz - dsdr * dtdz2))
    djacds = djacdr * drds + djacdz * dzds
    djacdt = djacdr * drdt + djacdz * dzdt
    drds2 = -(djacds * dtdz - jac * drds * dtdz / r + jac * dtdzs) / r
    dzds2 = (djacds * dtdr - jac * drds * dtdr / r + jac * dtdrs) / r
    drdst = -(djacdt * dtdz - jac * drdt * dtdz / r + jac * dtdzt) / r
    dzdst = (djacdt * dtdr - jac * drdt * dtdr / r + jac * dtdrt) / r
    drdt2 = (djacdt * dsdz - jac * drdt * dsdz / r + jac * dsdzt) / r
    dzdt2 = -(djacdt * dsdr - jac * drdt * dsdr / r + jac * dsdrt) / r

    dh11dr2 = 2 * (dsdr2 * dsdr2 + dsdr * dsdr3 + dsdrz * dsdrz + dsdz * dsdr2z)
    dh11drz = 2 * (dsdrz * dsdr2 + dsdr * dsdr2z + dsdz2 * dsdrz + dsdz * dsdrz2)
    dh11dz2 = 2 * (dsdrz * dsdrz + dsdr * dsdrz2 + dsdz2 * dsdz2 + dsdz * dsdz3)

    # field-amplitude (r,z) derivatives — FD block (dmagnetic.m parity)
    dbdr, dbdz, dbdr2, dbdrz, dbdz2 = _dmagnetic(eq, rhov, th)
    dbds2 = drds2 * dbdr + dzds2 * dbdz + dbdr2 * drds^2 + dbdz2 * dzds^2 + 2 * drds * dzds * dbdrz
    dbdt2 = drdt2 * dbdr + dzdt2 * dbdz + dbdr2 * drdt^2 + dbdz2 * dzdt^2 + 2 * drdt * dzdt * dbdrz
    dbdst =
        drdst * dbdr +
        dzdst * dbdz +
        drds * drdt * dbdr2 +
        drds * dzdt * dbdrz +
        dzds * drdt * dbdrz +
        dzds * dzdt * dbdz2

    # ----- basis-vector first derivatives (field-curvature terms) -----
    zeta = 1 / (r * b)
    dzetadr = -zeta * (1 / r + dbdr / b)
    dzetadz = -zeta * dbdz / b
    mag = 1 / gs
    dmagdr = -0.5 * mag^3 * dh11dr
    dmagdz = -0.5 * mag^3 * dh11dz

    denerdr = mag * dsdr2 + dsdr * dmagdr
    denefdr = 0.0
    denezdr = mag * dsdrz + dsdz * dmagdr
    denerdz = mag * dsdrz + dsdr * dmagdz
    denefdz = 0.0
    denezdz = mag * dsdz2 + dsdz * dmagdz

    deberdr = t0 * (dsdrz * zeta * mag + dsdz * dzetadr * mag + dsdz * zeta * dmagdr)
    debefdr = -(
        dsdr * dpds2 * zeta * gs +
        dpds * dzetadr * gs +
        dpds * zeta * (dsdr * dsdr2 + dsdz * dsdrz) / gs
    )
    debezdr = -t0 * (dsdr2 * zeta * mag + dsdr * dzetadr * mag + dsdr * zeta * dmagdr)
    deberdz = t0 * (dsdz2 * zeta * mag + dsdz * dzetadz * mag + dsdz * zeta * dmagdz)
    debefdz = -(
        dsdz * dpds2 * zeta * gs +
        dpds * dzetadz * gs +
        dpds * zeta * (dsdz * dsdz2 + dsdr * dsdrz) / gs
    )
    debezdz = -t0 * (dsdrz * zeta * mag + dsdr * dzetadz * mag + dsdr * zeta * dmagdz)

    deperdr = dsdr * dpds2 * dsdz * zeta + dpds * dsdrz * zeta + dpds * dsdz * dzetadr
    depefdr = t0 * dzetadr
    depezdr = -dsdr * dpds2 * dsdr * zeta - dpds * dsdr2 * zeta - dpds * dsdr * dzetadr
    deperdz = dsdz * dpds2 * dsdz * zeta + dpds * dsdz2 * zeta + dpds * dsdz * dzetadz
    depefdz = t0 * dzetadz
    depezdz = -dsdz * dpds2 * dsdr * zeta - dpds * dsdrz * zeta - dpds * dsdr * dzetadz

    # ----- basis-vector second derivatives -----
    dzetadr2 = -dzetadr * (1 / r + dbdr / b) - zeta * (-1 / r^2 - (dbdr / b)^2 + dbdr2 / b)
    dzetadrz = -dzetadr * dbdz / b - zeta * (dbdrz / b - dbdr * dbdz / b^2)
    dzetadz2 = -dzetadz * dbdz / b - zeta * (dbdz2 / b - (dbdz / b)^2)
    dmagdr2 = 0.75 * mag^5 * dh11dr^2 - 0.5 * dh11dr2 * mag^3
    dmagdrz = 0.75 * mag^5 * dh11dr * dh11dz - 0.5 * dh11drz * mag^3
    dmagdz2 = 0.75 * mag^5 * dh11dz^2 - 0.5 * dh11dz2 * mag^3
    dgsdr = 0.5 * mag * dh11dr
    dgsdz = 0.5 * mag * dh11dz
    dgsdr2 = 0.5 * (dmagdr * dh11dr + mag * dh11dr2)
    dgsdrz = 0.5 * (dmagdr * dh11dz + mag * dh11drz)
    dgsdz2 = 0.5 * (dmagdz * dh11dz + mag * dh11dz2)
    dpdsr = dsdr * dpds2
    dpdsz = dsdz * dpds2
    dpdsr2 = dsdr2 * dpds2
    dpdsrz = dsdrz * dpds2
    dpdsz2 = dsdz2 * dpds2

    denerdr2 = mag * dsdr3 + 2 * dmagdr * dsdr2 + dsdr * dmagdr2
    denerdrz = mag * dsdr2z + dmagdz * dsdr2 + dmagdr * dsdrz + dsdr * dmagdrz
    denerdz2 = mag * dsdrz2 + 2 * dmagdz * dsdrz + dsdr * dmagdz2
    denefdr2 = 0.0
    denefdrz = 0.0
    denefdz2 = 0.0
    denezdr2 = dmagdr2 * dsdz + 2 * dmagdr * dsdrz + mag * dsdr2z
    denezdrz = dmagdrz * dsdz + dmagdr * dsdz2 + dmagdz * dsdrz + mag * dsdrz2
    denezdz2 = dmagdz2 * dsdz + 2 * dmagdz * dsdz2 + mag * dsdz3

    deberdr2 =
        t0 * (
            dsdr2z * zeta * mag +
            dsdrz * dzetadr * mag +
            dsdrz * zeta * dmagdr +
            dsdrz * dzetadr * mag +
            dsdz * dzetadr2 * mag +
            dsdz * dzetadr * dmagdr +
            dsdrz * zeta * dmagdr +
            dsdz * dzetadr * dmagdr +
            dsdz * zeta * dmagdr2
        )
    deberdrz =
        t0 * (
            dsdrz2 * zeta * mag +
            dsdz2 * dzetadr * mag +
            dsdz2 * zeta * dmagdr +
            dsdrz * dzetadz * mag +
            dsdz * dzetadrz * mag +
            dsdz * dzetadz * dmagdr +
            dsdrz * zeta * dmagdz +
            dsdz * dzetadr * dmagdz +
            dsdz * zeta * dmagdrz
        )
    deberdz2 =
        t0 * (
            dsdz3 * zeta * mag +
            dsdz2 * dzetadz * mag +
            dsdz2 * zeta * dmagdz +
            dsdz2 * dzetadz * mag +
            dsdz * dzetadz2 * mag +
            dsdz * dzetadz * dmagdz +
            dsdz2 * zeta * dmagdz +
            dsdz * dzetadz * dmagdz +
            dsdz * zeta * dmagdz2
        )

    debefdr2 = -(
        dpdsr2 * zeta * gs +
        dpdsr * dzetadr * gs +
        dpdsr * zeta * dgsdr +
        dpdsr * dzetadr * gs +
        dpds * dzetadr2 * gs +
        dpds * dzetadr * dgsdr +
        dpdsr * zeta * dgsdr +
        dpds * dzetadr * dgsdr +
        dpds * zeta * dgsdr2
    )
    debefdrz = -(
        dpdsrz * zeta * gs +
        dpdsz * dzetadr * gs +
        dpdsz * zeta * dgsdr +
        dpdsr * dzetadz * gs +
        dpds * dzetadrz * gs +
        dpds * dzetadz * dgsdr +
        dpdsr * zeta * dgsdz +
        dpds * dzetadr * dgsdz +
        dpds * zeta * dgsdrz
    )
    debefdz2 = -(
        dpdsz2 * zeta * gs +
        dpdsz * dzetadz * gs +
        dpdsz * zeta * dgsdz +
        dpdsz * dzetadz * gs +
        dpds * dzetadz2 * gs +
        dpds * dzetadz * dgsdz +
        dpdsz * zeta * dgsdz +
        dpds * dzetadz * dgsdz +
        dpds * zeta * dgsdz2
    )

    debezdr2 =
        -t0 * (
            dsdr3 * zeta * mag +
            dsdr2 * dzetadr * mag +
            dsdr2 * zeta * dmagdr +
            dsdr2 * dzetadr * mag +
            dsdr * dzetadr2 * mag +
            dsdr * dzetadr * dmagdr +
            dsdr2 * zeta * dmagdr +
            dsdr * dzetadr * dmagdr +
            dsdr * zeta * dmagdr2
        )
    debezdrz =
        -t0 * (
            dsdr2z * zeta * mag +
            dsdrz * dzetadr * mag +
            dsdrz * zeta * dmagdr +
            dsdr2 * dzetadz * mag +
            dsdr * dzetadrz * mag +
            dsdr * dzetadz * dmagdr +
            dsdr2 * zeta * dmagdz +
            dsdr * dzetadr * dmagdz +
            dsdr * zeta * dmagdrz
        )
    debezdz2 =
        -t0 * (
            dsdrz2 * zeta * mag +
            dsdrz * dzetadz * mag +
            dsdrz * zeta * dmagdz +
            dsdrz * dzetadz * mag +
            dsdr * dzetadz2 * mag +
            dsdr * dzetadz * dmagdz +
            dsdrz * zeta * dmagdz +
            dsdr * dzetadz * dmagdz +
            dsdr * zeta * dmagdz2
        )

    deperdr2 =
        dsdr2z * zeta * dpds +
        dsdrz * dzetadr * dpds +
        dsdrz * zeta * dpdsr +
        dsdrz * dzetadr * dpds +
        dsdz * dzetadr2 * dpds +
        dsdz * dzetadr * dpdsr +
        dsdrz * zeta * dpdsr +
        dsdz * dzetadr * dpdsr +
        dsdz * zeta * dpdsr2
    deperdrz =
        dsdrz2 * zeta * dpds +
        dsdz2 * dzetadr * dpds +
        dsdz2 * zeta * dpdsr +
        dsdrz * dzetadz * dpds +
        dsdz * dzetadrz * dpds +
        dsdz * dzetadz * dpdsr +
        dsdrz * zeta * dpdsz +
        dsdz * dzetadr * dpdsz +
        dsdz * zeta * dpdsrz
    deperdz2 =
        dsdz3 * zeta * dpds +
        dsdz2 * dzetadz * dpds +
        dsdz2 * zeta * dpdsz +
        dsdz2 * dzetadz * dpds +
        dsdz * dzetadz2 * dpds +
        dsdz * dzetadz * dpdsz +
        dsdz2 * zeta * dpdsz +
        dsdz * dzetadz * dpdsz +
        dsdz * zeta * dpdsz2

    depefdr2 = t0 * dzetadr2
    depefdrz = t0 * dzetadrz
    depefdz2 = t0 * dzetadz2

    depezdr2 = -(
        dsdr3 * zeta * dpds +
        dsdr2 * dzetadr * dpds +
        dsdr2 * zeta * dpdsr +
        dsdr2 * dzetadr * dpds +
        dsdr * dzetadr2 * dpds +
        dsdr * dzetadr * dpdsr +
        dsdr2 * zeta * dpdsr +
        dsdr * dzetadr * dpdsr +
        dsdr * zeta * dpdsr2
    )
    depezdrz = -(
        dsdr2z * zeta * dpds +
        dsdrz * dzetadr * dpds +
        dsdrz * zeta * dpdsr +
        dsdr2 * dzetadz * dpds +
        dsdr * dzetadrz * dpds +
        dsdr * dzetadz * dpdsr +
        dsdr2 * zeta * dpdsz +
        dsdr * dzetadr * dpdsz +
        dsdr * zeta * dpdsrz
    )
    depezdz2 = -(
        dsdrz2 * zeta * dpds +
        dsdrz * dzetadz * dpds +
        dsdrz * zeta * dpdsz +
        dsdrz * dzetadz * dpds +
        dsdr * dzetadz2 * dpds +
        dsdr * dzetadz * dpdsz +
        dsdrz * zeta * dpdsz +
        dsdr * dzetadz * dpdsz +
        dsdr * zeta * dpdsz2
    )

    return (;
        b,
        dbds,
        dbdt,
        bp,
        sflx,
        dsdr,
        dsdz,
        dtdr,
        dtdz,
        ener,
        enef,
        enez,
        eber,
        ebef,
        ebez,
        eper,
        epef,
        epez,
        dbds2,
        dbdst,
        dbdt2,
        dsdr2,
        dsdrz,
        dsdz2,
        dsdr3,
        dsdr2z,
        dsdrz2,
        dsdz3,
        dtdr2,
        dtdrz,
        dtdz2,
        denerdr,
        denefdr,
        denezdr,
        denerdz,
        denefdz,
        denezdz,
        deberdr,
        debefdr,
        debezdr,
        deberdz,
        debefdz,
        debezdz,
        deperdr,
        depefdr,
        depezdr,
        deperdz,
        depefdz,
        depezdz,
        denerdr2,
        denerdrz,
        denerdz2,
        denefdr2,
        denefdrz,
        denefdz2,
        denezdr2,
        denezdrz,
        denezdz2,
        deberdr2,
        deberdrz,
        deberdz2,
        debefdr2,
        debefdrz,
        debefdz2,
        debezdr2,
        debezdrz,
        debezdz2,
        deperdr2,
        deperdrz,
        deperdz2,
        depefdr2,
        depefdrz,
        depefdz2,
        depezdr2,
        depezdrz,
        depezdz2,
        # extras used internally by the dispersion layer
        dbdr,
        dbdz,
    )
end

"""
    plasma_profiles(prob, sflx) -> (; n, T, dLNnds, dLNnds2)

Per-species density/temperature profiles and logarithmic density derivatives
at flux label `s` (identical in dispertok.m/disp_eig.m):
`n = n0·(1−na·s²)^nb`, `T = t0·(1−ta·s²)^nb` (upstream quirk: `nb`, not `tb`),
`dLNn/ds = −2·s·na·nb/(1−na·s²)`, `dLNn/ds² = dLNn/ds / s − (dLNn/ds)²/nb`.
"""
function plasma_profiles(prob::RayconProblem, sflx::Real)
    ns = length(prob.amass)
    n = Vector{Float64}(undef, ns)
    T = Vector{Float64}(undef, ns)
    dLNnds = Vector{Float64}(undef, ns)
    dLNnds2 = Vector{Float64}(undef, ns)
    s = Float64(sflx)
    for k = 1:ns
        p = 1 - s^2 * prob.na[k]
        p > 0 || throw(DomainError(s, "flux label outside the plasma (1 − na·s² ≤ 0)"))
        pt = 1 - s^2 * prob.ta[k]
        pt > 0 || throw(DomainError(s, "flux label outside the temperature profile"))
        dLNnds[k] = -2 * s * prob.na[k] * prob.nb[k] / p
        # upstream writes dLNnds/s − dLNnds²/nb, which is 0/0 on the magnetic
        # axis; substitute the algebraically identical dLNnds/s = −2·na·nb/p
        dLNnds2[k] = -2 * prob.na[k] * prob.nb[k] / p - dLNnds[k]^2 / prob.nb[k]
        n[k] = prob.n0[k] * p^prob.nb[k]
        T[k] = prob.t0[k] * pt^prob.nb[k]
    end
    return (; n, T, dLNnds, dLNnds2)
end
