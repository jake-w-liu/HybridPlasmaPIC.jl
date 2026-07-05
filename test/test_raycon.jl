# RCN-001..011 — RAYCON port (src/raycon/): Solovev equilibrium, magnetic
# geometry, eigenvalue dispersion, Hamiltonian tracing, mode conversion.
#
# Verification layers:
#   1. finite-difference oracles for every analytic derivative;
#   2. exact identities (basis orthonormality, eigenpair residuals, symplectic
#      symmetry, complex-Γ identities, U conservation along rays);
#   3. machine-precision comparison against reference data dumped from the
#      ORIGINAL MATLAB RAYCON (test/reference/raycon_reference.txt, generated
#      by tools/raycon_reference.m) — runs when the file is present.

using HybridPlasmaPIC
using HybridPlasmaPIC.Raycon
using Test, LinearAlgebra

const RCN_PROB = cmod_parameters()
const RCN_EQ = RCN_PROB.eq

_rz(rho, th) = (rho * cos(th) + RCN_EQ.r0, -rho * sin(th))
_rho_th(r, z) = (sqrt((r - RCN_EQ.r0)^2 + z^2), atan(RCN_EQ.r0 - r, -z) + π / 2)

@testset "RCN-001 Solovev flux derivatives vs finite differences" begin
    for (rho, th) in ((0.05, 0.3), (0.12, 0.8), (0.2, 2.0), (0.16, -1.2))
        sd = solovev_flux(RCN_EQ, rho, th; order = 3)
        r, z = _rz(rho, th)
        sval(r_, z_) = solovev_flux(RCN_EQ, _rho_th(r_, z_)...).sflx
        h = 1e-6
        @test isapprox(sd.dsdr, (sval(r + h, z) - sval(r - h, z)) / 2h; rtol = 1e-7)
        @test isapprox(sd.dsdz, (sval(r, z + h) - sval(r, z - h)) / 2h; rtol = 1e-7)
        h2 = 1e-4
        fd_dsdr2 = (sval(r + h2, z) - 2 * sval(r, z) + sval(r - h2, z)) / h2^2
        fd_dsdz2 = (sval(r, z + h2) - 2 * sval(r, z) + sval(r, z - h2)) / h2^2
        fd_dsdrz =
            (
                sval(r + h2, z + h2) - sval(r + h2, z - h2) - sval(r - h2, z + h2) +
                sval(r - h2, z - h2)
            ) / (4 * h2^2)
        @test isapprox(sd.dsdr2, fd_dsdr2; rtol = 1e-5, atol = 1e-8)
        @test isapprox(sd.dsdz2, fd_dsdz2; rtol = 1e-5, atol = 1e-8)
        @test isapprox(sd.dsdrz, fd_dsdrz; rtol = 1e-5, atol = 1e-8)
        # 3rd derivatives: FD of the analytic 2nd derivatives
        s2r(r_, z_) = solovev_flux(RCN_EQ, _rho_th(r_, z_)...; order = 2).dsdr2
        s2z(r_, z_) = solovev_flux(RCN_EQ, _rho_th(r_, z_)...; order = 2).dsdz2
        @test isapprox(sd.dsdr3, (s2r(r + h, z) - s2r(r - h, z)) / 2h; rtol = 1e-5, atol = 1e-6)
        @test isapprox(sd.dsdr2z, (s2r(r, z + h) - s2r(r, z - h)) / 2h; rtol = 1e-5, atol = 1e-6)
        @test isapprox(sd.dsdz3, (s2z(r, z + h) - s2z(r, z - h)) / 2h; rtol = 1e-5, atol = 1e-6)
        @test isapprox(sd.dsdrz2, (s2z(r + h, z) - s2z(r - h, z)) / 2h; rtol = 1e-5, atol = 1e-6)
    end
end

@testset "RCN-002 flux mapping round trip" begin
    for s in (0.15, 0.4, 0.7, 0.95), th in (0.001, 0.8, 2.0, -1.2, 3.1)
        mp = map_flux(RCN_EQ, s, th)
        @test isapprox(solovev_flux(RCN_EQ, mp.rho, th).sflx, s; atol = 1e-10)
        @test isapprox(mp.r, mp.rho * cos(th) + RCN_EQ.r0; rtol = 1e-14)
        @test isapprox(mp.z, -mp.rho * sin(th); atol = 1e-14)
    end
    mesh = flux_surface_mesh(RCN_EQ; ns = 6, nt = 8)
    @test size(mesh.r) == (6, 9)
    @test all(isfinite, mesh.r) && all(isfinite, mesh.z)
    @test_throws ArgumentError map_flux(RCN_EQ, -0.1, 0.0)
end

@testset "RCN-003 magnetic geometry: identities and FD oracles" begin
    for (rho, th) in ((0.08, 0.5), (0.15, 2.2), (0.11, -0.9))
        g = magnetic_geometry(RCN_EQ, rho, th)
        en = [g.ener, g.enef, g.enez]
        eb = [g.eber, g.ebef, g.ebez]
        ep = [g.eper, g.epef, g.epez]
        for e in (en, eb, ep)
            @test isapprox(norm(e), 1.0; rtol = 1e-12)
        end
        @test abs(dot(en, eb)) < 1e-12
        @test abs(dot(en, ep)) < 1e-12
        @test abs(dot(eb, ep)) < 1e-12
        # dbds via FD along s at fixed theta
        h = 1e-7
        s0 = g.sflx
        bofs(s_) = magnetic_geometry(RCN_EQ, map_flux(RCN_EQ, s_, th).rho, th).b
        @test isapprox(g.dbds, (bofs(s0 + h) - bofs(s0 - h)) / 2h; rtol = 1e-5)
        # basis-vector first derivatives vs FD in (r, z)
        r, z = _rz(rho, th)
        hb = 1e-6
        comp(r_, z_, f) = f(magnetic_geometry(RCN_EQ, _rho_th(r_, z_)...))
        for (an, f) in (
            (g.denerdr, gg -> gg.ener),
            (g.denezdr, gg -> gg.enez),
            (g.deberdr, gg -> gg.eber),
            (g.debefdr, gg -> gg.ebef),
            (g.debezdr, gg -> gg.ebez),
            (g.deperdr, gg -> gg.eper),
            (g.depefdr, gg -> gg.epef),
            (g.depezdr, gg -> gg.epez),
        )
            fd = (comp(r + hb, z, f) - comp(r - hb, z, f)) / 2hb
            @test isapprox(an, fd; rtol = 2e-4, atol = 1e-6)
        end
        for (an, f) in (
            (g.denerdz, gg -> gg.ener),
            (g.denezdz, gg -> gg.enez),
            (g.deberdz, gg -> gg.eber),
            (g.debefdz, gg -> gg.ebef),
            (g.debezdz, gg -> gg.ebez),
            (g.deperdz, gg -> gg.eper),
            (g.depefdz, gg -> gg.epef),
            (g.depezdz, gg -> gg.epez),
        )
            fd = (comp(r, z + hb, f) - comp(r, z - hb, f)) / 2hb
            @test isapprox(an, fd; rtol = 2e-4, atol = 1e-6)
        end
        # basis-vector second derivatives vs FD of the analytic first derivatives
        for (an, f) in (
            (g.denerdr2, gg -> gg.denerdr),
            (g.denerdrz, gg -> gg.denerdz),
            (g.deberdr2, gg -> gg.deberdr),
            (g.deberdrz, gg -> gg.deberdz),
            (g.debefdr2, gg -> gg.debefdr),
            (g.depezdr2, gg -> gg.depezdr),
        )
            fd = (comp(r + hb, z, f) - comp(r - hb, z, f)) / 2hb
            @test isapprox(an, fd; rtol = 1e-3, atol = 1e-4)
        end
    end
end

@testset "RCN-004 dispersion tensor eigenpair and Stix layer" begin
    pts = (
        (0.80, 0.05, -31.5, 5.0),
        (0.75, -0.10, -20.0, 0.0),
        (0.62, 0.02, 10.0, 3.0),
        (0.55, 0.12, -31.5, 10.0),
    )
    st = stix_elements(RCN_PROB, 0.75, 0.0)
    @test st.S isa Float64 && isfinite(st.S) && isfinite(st.D) && isfinite(st.P)
    @test st.omc[1] < 0            # electron cyclotron frequency is signed
    fr = cyclotron_frequencies(RCN_PROB, 0.75, 0.0)
    @test length(fr) == 3 && all(isfinite, fr)
    # ion-ion hybrid frequency lies between the two ion cyclotron frequencies
    @test min(fr[2], fr[3]) < fr[1] < max(fr[2], fr[3])
    for y in pts
        U = dispersion_U(RCN_PROB, y)
        pol = polarization(RCN_PROB, y)
        m = conversion_monitors(RCN_PROB, y)
        @test isfinite(U) && isfinite(m.mon2) && m.mon1 == 0.0
        @test isapprox(norm(pol), 1.0; rtol = 1e-12)
        # eigenpair residual: DD·pol = U·pol for the assembled tensor
        core = HybridPlasmaPIC.Raycon._disp_core(RCN_PROB, Float64.(y)...)
        resid = norm(core.DD * pol .- U .* pol)
        @test resid <= 1e-10 * max(1.0, maximum(abs.(core.DD)))
    end
end

@testset "RCN-005 exact eigenvalue gradient vs finite differences" begin
    for model in (:cld2x2, :cld3x3)
        prob = cmod_parameters(; model)
        for y in ((0.80, 0.05, -31.5, 5.0), (0.62, 0.02, 10.0, 3.0), (0.75, -0.1, -20.0, 0.0))
            rhs = trajectory_rhs(prob, collect(y))
            @test all(isfinite, rhs)
            hx = 1e-7          # position step [m]
            hk = 1e-4          # wavenumber step [1/m]
            Uf(yy) = dispersion_U(prob, yy)
            dUdr = (Uf((y[1] + hx, y[2], y[3], y[4])) - Uf((y[1] - hx, y[2], y[3], y[4]))) / 2hx
            dUdz = (Uf((y[1], y[2] + hx, y[3], y[4])) - Uf((y[1], y[2] - hx, y[3], y[4]))) / 2hx
            dUdkr = (Uf((y[1], y[2], y[3] + hk, y[4])) - Uf((y[1], y[2], y[3] - hk, y[4]))) / 2hk
            dUdkz = (Uf((y[1], y[2], y[3], y[4] + hk)) - Uf((y[1], y[2], y[3], y[4] - hk))) / 2hk
            # dz/dσ = (dU/dkr, dU/dkz, −dU/dr, −dU/dz)
            @test isapprox(rhs[1], dUdkr; rtol = 2e-5, atol = 1e-12)
            @test isapprox(rhs[2], dUdkz; rtol = 2e-5, atol = 1e-12)
            @test isapprox(rhs[3], -dUdr; rtol = 2e-5, atol = 1e-8)
            @test isapprox(rhs[4], -dUdz; rtol = 2e-5, atol = 1e-8)
        end
    end
end

@testset "RCN-006 dU/dω vs frequency-perturbed problems" begin
    y = (0.75, -0.10, -20.0, 0.0)
    dom = dUdomega(RCN_PROB, y)
    h = RCN_PROB.freq * 1e-7
    mk(freq) = RayconProblem(;
        eq = RCN_EQ,
        amass = RCN_PROB.amass,
        acharge = RCN_PROB.acharge,
        n0 = RCN_PROB.n0,
        na = RCN_PROB.na,
        nb = RCN_PROB.nb,
        t0 = RCN_PROB.t0,
        ta = RCN_PROB.ta,
        tb = RCN_PROB.tb,
        freq,
        kphi = RCN_PROB.kphi,
        model = RCN_PROB.model,
    )
    fd =
        (dispersion_U(mk(RCN_PROB.freq + h), y) - dispersion_U(mk(RCN_PROB.freq - h), y)) /
        (2π * 2h)
    @test isapprox(dom, fd; rtol = 1e-5)
end

@testset "RCN-007 symplectic tangent map and Hessian symmetry" begin
    for model in (:cld2x2, :cld3x3)
        prob = cmod_parameters(; model)
        y20 = vcat([0.75, -0.10, -20.0, 0.0], vec(Matrix{Float64}(I, 4, 4)))
        rhs = trajectory_rhs(prob, y20)
        @test length(rhs) == 20
        @test all(isfinite, rhs)
        @test rhs[1:4] == trajectory_rhs(prob, y20[1:4])
        # dS = J·H·S with H = ∇∇U symmetric ⟹ H = −J·dS is symmetric for S = I
        J = [0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; -1.0 0.0 0.0 0.0; 0.0 -1.0 0.0 0.0]
        H = -J * reshape(rhs[5:20], 4, 4)
        @test maximum(abs.(H .- transpose(H))) <= 1e-8 * max(1.0, maximum(abs.(H)))
        # Hessian vs FD of the trajectory RHS: kr column AND kz column (the
        # kz column is the one the upstream dkz²/dkr·dkz swap bug corrupted)
        JH = J * H
        hk = 1e-3
        for (col, dy) in ((3, [0.0, 0.0, hk, 0.0]), (4, [0.0, 0.0, 0.0, hk]))
            r1 = trajectory_rhs(prob, [0.75, -0.10, -20.0, 0.0] .+ dy)
            r2 = trajectory_rhs(prob, [0.75, -0.10, -20.0, 0.0] .- dy)
            fdcol = (r1[1:4] .- r2[1:4]) ./ (2hk)
            for i = 1:4
                @test isapprox(JH[i, col], fdcol[i]; rtol = 5e-3, atol = 1e-6)
            end
        end
    end
end

@testset "RCN-008 complex gamma function" begin
    @test isapprox(real(cgamma(1.0)), 1.0; rtol = 1e-12)
    @test abs(imag(cgamma(1.0))) < 1e-14
    @test isapprox(real(cgamma(0.5)), sqrt(π); rtol = 1e-12)
    @test isapprox(real(cgamma(5.0)), 24.0; rtol = 1e-12)
    for yv in (0.25, 1.0, 2.7)
        @test isapprox(abs2(cgamma(im * yv)), π / (yv * sinh(π * yv)); rtol = 1e-11)
        @test isapprox(abs2(cgamma(1 + im * yv)), π * yv / sinh(π * yv); rtol = 1e-11)
    end
    for z in (0.3 + 0.7im, -1.4 + 2.2im, 2.0 - 3.0im)
        @test isapprox(cgamma(z + 1), z * cgamma(z); rtol = 1e-10)
        @test isapprox(cgamma(conj(z)), conj(cgamma(z)); rtol = 1e-12)
    end
end

@testset "RCN-009 ray tracing conserves U; launch sits on the surface" begin
    y0 = launch_ray(RCN_PROB; s = 0.4, theta = 0.001, kr = -31.5, kz = 0.0)
    U0 = dispersion_U(RCN_PROB, y0)
    # U scale: a unit change of N² changes U by O(1); require machine-level launch
    @test abs(U0) < 1e-8
    tr = integrate_ray(RCN_PROB, y0, 0.0, 1e-3)
    @test tr.status in (:end_of_span, :conversion_event)
    @test size(tr.y, 2) == length(tr.sigma) >= 10
    Umax = maximum(abs(dispersion_U(RCN_PROB, tr.y[:, i])) for i = 1:size(tr.y, 2))
    @test Umax < 1e-4                       # conserved at integrator tolerance
    # the ray stays inside the plasma
    smax = maximum(
        solovev_flux(RCN_EQ, _rho_th(tr.y[1, i], tr.y[2, i])...).sflx for i = 1:size(tr.y, 2)
    )
    @test smax < 1.0
end

@testset "RCN-010 cmod conversion run: τ, β, ray splitting" begin
    res = trace_rays(RCN_PROB; s = 0.4, theta = 0.001, kr = -31.5, kz = 0.0, sigma_span = 5e-2)
    @test length(res.rays) >= 1
    @test all(r -> all(isfinite, r.y), res.rays)
    # this launch deterministically hits the ion-ion hybrid layer: conversions
    # MUST be found (a silent zero would previously have passed vacuously)
    @test !isempty(res.conversions)
    for c in res.conversions
        cv = c.conversion
        @test cv.converged && cv.hyperbola_ok && cv.transmitted_ok
        @test 0.0 < cv.tau < 1.0
        @test isfinite(cv.eta2) && cv.eta2 > 0
        @test isfinite(abs(cv.beta))
        # transmitted ray launches on the dispersion surface
        @test abs(dispersion_U(RCN_PROB, collect(cv.transmitted))) < 1e-3
        # converted ray continues from the detection point
        @test cv.converted == cv.incoming
    end
    # every valid conversion queued exactly one transmitted ray
    @test length(res.rays) == 1 + length(res.conversions)

    # cap semantics: max_conversions = 0 disables splitting entirely
    res0 = trace_rays(
        RCN_PROB;
        s = 0.4,
        theta = 0.001,
        kr = -31.5,
        kz = 0.0,
        sigma_span = 5e-2,
        max_conversions = 0,
    )
    @test isempty(res0.conversions)
    @test length(res0.rays) == 1
    @test res0.rays[1].status === :end_of_span

    # cld3x3 traces without splitting instead of aborting mid-run
    res3 = trace_rays(
        cmod_parameters(; model = :cld3x3);
        s = 0.4,
        theta = 0.001,
        kr = -31.5,
        kz = 0.0,
        sigma_span = 5e-3,
    )
    @test length(res3.rays) == 1
    @test isempty(res3.conversions)
    @test res3.rays[1].status === :end_of_span
end

@testset "RCN-012 integrator/driver edge cases (review findings)" begin
    y0 = launch_ray(RCN_PROB; s = 0.4, theta = 0.001, kr = -31.5, kz = 0.0)
    # exhausting the step budget reports :max_steps even when the final
    # iteration is a rejected step
    trm = integrate_ray(
        RCN_PROB,
        y0,
        0.0,
        5e-2;
        initial_step = 2e-4,
        rtol = 1e-12,
        atol = 1e-14,
        max_steps = 1,
        detect_conversion = false,
    )
    @test trm.status === :max_steps
    # tiny σ spans integrate instead of instantly underflowing
    trt = integrate_ray(RCN_PROB, y0, 0.0, 1e-8; detect_conversion = false)
    @test trt.status === :end_of_span
    @test trt.sigma[end] ≈ 1e-8 rtol = 1e-6
    # leaving the plasma terminates cleanly
    yout = [RCN_EQ.r0 + 1.3 * rho_edge(RCN_EQ), 0.0, -31.5, 0.0]
    trl = integrate_ray(RCN_PROB, yout, 0.0, 1e-3; detect_conversion = false)
    @test trl.status === :left_domain
    # _fzero_near hands Brent a valid bracket when one search side leaves the
    # function domain (NaN) — regression for the dropped-transmitted-ray bug
    @test HybridPlasmaPIC.Raycon._fzero_near(x -> x < 0 ? NaN : x - 1.0, 0.5) ≈ 1.0 atol = 1e-10
    # plasma_profiles is finite on the magnetic axis (removable singularity)
    prof0 = plasma_profiles(RCN_PROB, 0.0)
    @test all(isfinite, prof0.dLNnds2)
    @test all(iszero, prof0.dLNnds)
    # oblate equilibria (elong < 1) map interior flux surfaces
    eqo = SolovevEquilibrium(; b0 = 1.0, r0 = 1.0, q0 = 1.0, iaspr = 0.3, elong = 0.5)
    mo = map_flux(eqo, 0.9, 0.0)
    @test isapprox(solovev_flux(eqo, mo.rho, 0.0).sflx, 0.9; atol = 1e-10)
    # nb = 0 profile exponents are rejected at construction
    @test_throws ArgumentError RayconProblem(;
        eq = RCN_EQ,
        amass = [1.0],
        acharge = [1.0],
        n0 = [1e19],
        na = [1.0],
        nb = [0.0],
        t0 = [3.0],
        ta = [1.0],
        tb = [1.0],
        freq = 8e7,
        kphi = -10.0,
    )
end

@testset "RCN-013 Ω_ci-normalized (PlasmaUnits) interface" begin
    units = cmod_units()
    di = inertial_length(units)
    Ωci = gyrofrequency(units)
    vA = alfven_speed(units)
    @test isapprox(di, vA / Ωci; rtol = 1e-14)

    # exact-parity construction: with the upstream constants the normalized
    # constructor reproduces the SI cmod problem to rounding
    temp_scale = units.mi * vA^2 / (units.e * 1000)
    pn = RayconProblem(
        units;
        r0 = 0.67 / di,
        iaspr = 0.22 / 0.67,
        elong = 1.6,
        q0 = 2.0,
        b0 = 7.9 / units.B0,
        amass = [1 / 1836, 2.0, 3.0] .* (1.6726e-27 / units.mi),
        acharge = [-1.0, 1.0, 2.0],
        n0 = [10.0, 5.2, 2.4] .* 1e19 ./ units.n0,
        na = [1.0, 0.7, 0.7],
        nb = [3.0, 3.0, 3.0],
        t0 = [3.0, 3.0, 3.0] ./ temp_scale,
        ta = [1.0, 1.0, 1.0],
        tb = [1.0, 1.0, 1.0],
        omega = 2π * 80.0e6 / Ωci,
        kphi = -10.0 * di,
        cnst = RayconConstants(),
    )
    ps = cmod_parameters()
    @test isapprox(pn.eq.r0, ps.eq.r0; rtol = 1e-13)
    @test isapprox(pn.eq.b0, ps.eq.b0; rtol = 1e-13)
    @test isapprox(pn.eq.psin, ps.eq.psin; rtol = 1e-13)
    @test all(isapprox.(pn.n0, ps.n0; rtol = 1e-13))
    @test all(isapprox.(pn.t0, ps.t0; rtol = 1e-12))
    @test all(isapprox.(pn.amass, ps.amass; rtol = 1e-13))
    @test isapprox(pn.freq, ps.freq; rtol = 1e-13)
    @test isapprox(pn.omega, ps.omega; rtol = 1e-13)
    @test isapprox(pn.kphi, ps.kphi; rtol = 1e-13)

    # the convenience preset takes the same path (package CODATA constants)
    pu = cmod_parameters(units)
    @test isapprox(pu.eq.r0, ps.eq.r0; rtol = 1e-13)
    @test isapprox(pu.freq, ps.freq; rtol = 1e-13)

    # normalized launch/trace equal the exactly-rescaled SI results
    y_n = launch_ray(units, pn; s = 0.4, theta = 0.001, kr = -31.5 * di, kz = 0.0)
    y_s = launch_ray(pn; s = 0.4, theta = 0.001, kr = -31.5, kz = 0.0)
    @test isapprox(y_n[1], y_s[1] / di; rtol = 1e-13)
    @test isapprox(y_n[2], y_s[2] / di; rtol = 1e-12, atol = 1e-15)
    @test isapprox(y_n[3], y_s[3] * di; rtol = 1e-12)
    @test isapprox(y_n[4], y_s[4] * di; rtol = 1e-12)

    trn = integrate_ray(units, pn, y_n, 0.0, 2e-3; detect_conversion = false)
    trs = integrate_ray(pn, y_s, 0.0, 2e-3; detect_conversion = false)
    @test trn.status === trs.status === :end_of_span
    @test length(trn.sigma) == length(trs.sigma)
    @test maximum(abs.(trn.y[1:2, :] .- trs.y[1:2, :] ./ di)) <=
          1e-9 * maximum(abs.(trs.y[1:2, :] ./ di))
    @test maximum(abs.(trn.y[3:4, :] .- trs.y[3:4, :] .* di)) <=
          1e-9 * maximum(abs.(trs.y[3:4, :] .* di))

    resn =
        trace_rays(units, pn; s = 0.4, theta = 0.001, kr = -31.5 * di, kz = 0.0, sigma_span = 2e-2)
    ress = trace_rays(pn; s = 0.4, theta = 0.001, kr = -31.5, kz = 0.0, sigma_span = 2e-2)
    @test length(resn.rays) == length(ress.rays)
    @test length(resn.conversions) == length(ress.conversions)
    for (cn, cs) in zip(resn.conversions, ress.conversions)
        @test isapprox(cn.conversion.tau, cs.conversion.tau; rtol = 1e-12)
        @test isapprox(cn.conversion.eta2, cs.conversion.eta2; rtol = 1e-12)
        @test isapprox(cn.conversion.saddle[1], cs.conversion.saddle[1] / di; rtol = 1e-12)
        @test isapprox(cn.conversion.saddle[3], cs.conversion.saddle[3] * di; rtol = 1e-12)
    end

    # dimensionless sanity of the normalized quantities
    @test isapprox(pn.eq.r0 / di, 0.67 / di; rtol = 1e-13)   # R0 ≈ 29-34 d_i for cmod
    @test 0.5 < pn.omega / Ωci < 1.0                          # below proton Ω_ci
end

@testset "RCN-011 validation errors" begin
    @test_throws ArgumentError SolovevEquilibrium(;
        b0 = -1.0,
        r0 = 0.67,
        q0 = 2.0,
        iaspr = 0.3,
        elong = 1.6,
    )
    @test_throws DimensionMismatch RayconProblem(;
        eq = RCN_EQ,
        amass = [1.0, 2.0],
        acharge = [-1.0],
        n0 = [1e19, 1e19],
        na = [1.0, 1.0],
        nb = [3.0, 3.0],
        t0 = [3.0, 3.0],
        ta = [1.0, 1.0],
        tb = [1.0, 1.0],
        freq = 8e7,
        kphi = -10.0,
    )
    @test_throws ArgumentError cmod_parameters(; model = :bogus)
    @test_throws ArgumentError solovev_flux(RCN_EQ, -0.1, 0.0)
    @test_throws ArgumentError dispersion_U(RCN_PROB, (NaN, 0.0, 1.0, 1.0))
    @test_throws ArgumentError trajectory_rhs(RCN_PROB, [1.0, 2.0, 3.0])
    @test_throws ArgumentError integrate_ray(RCN_PROB, [0.7, 0.0, -30.0, 0.0], 0.0, 0.0)
    p3 = cmod_parameters(; model = :cld3x3)
    @test_throws ArgumentError analyze_conversion(
        p3,
        [0.62, 0.02, 10.0, 3.0],
        [0.5, -0.3, 2000.0, 1000.0],
        [0.1, 0.2, -500.0, 300.0],
    )
    # msw1x1 provides the launch estimate but not the eigenvalue tracer
    pm = cmod_parameters(; model = :msw1x1)
    @test isfinite(msw_dispersion(pm, (0.75, 0.0, -20.0, 0.0)))
    @test_throws ArgumentError trajectory_rhs(pm, [0.75, 0.0, -20.0, 0.0])
end

# ---------------------------------------------------------------- MATLAB reference
const RCN_REF = joinpath(@__DIR__, "reference", "raycon_reference.txt")
if isfile(RCN_REF)
    # flat "name<TAB>v1 v2 ..." format converted from the MATLAB JSON dump
    refdata = Dict{String,Vector{Float64}}()
    for line in eachline(RCN_REF)
        isempty(strip(line)) && continue
        name, rest = split(line, '\t'; limit = 2)
        refdata[name] = parse.(Float64, split(rest))
    end
    @testset "RCN-REF MATLAB reference comparison" begin
        npts = Int(refdata["solovev/npts"][1])
        for k = 1:npts
            rho, th = refdata["solovev/$k/rho"][1], refdata["solovev/$k/theta"][1]
            sd = solovev_flux(RCN_EQ, rho, th; order = 3)
            for (fld, val) in (
                ("sflx", sd.sflx),
                ("dsdr", sd.dsdr),
                ("dsdz", sd.dsdz),
                ("dsdr2", sd.dsdr2),
                ("dsdrz", sd.dsdrz),
                ("dsdz2", sd.dsdz2),
                ("dsdr3", sd.dsdr3),
                ("dsdr2z", sd.dsdr2z),
                ("dsdrz2", sd.dsdrz2),
                ("dsdz3", sd.dsdz3),
            )
                @test isapprox(val, refdata["solovev/$k/$fld"][1]; rtol = 1e-9, atol = 1e-9)
            end
            g = magnetic_geometry(RCN_EQ, rho, th)
            for fld in (
                "b",
                "dbds",
                "dbdt",
                "bp",
                "ener",
                "enez",
                "eber",
                "ebef",
                "ebez",
                "eper",
                "epef",
                "epez",
            )
                @test isapprox(
                    getproperty(g, Symbol(fld)),
                    refdata["magnetic/$k/$fld"][1];
                    rtol = 1e-7,
                    atol = 1e-9,
                )
            end
            # first derivatives of basis vectors: dominated by analytic terms,
            # dbdr enters via FD (upstream step 1e-8 vs ours 1e-5) → 1e-5 slack
            for fld in (
                "denerdr",
                "denezdr",
                "deberdr",
                "debefdr",
                "debezdr",
                "deperdr",
                "depefdr",
                "depezdr",
                "denerdz",
                "denezdz",
                "deberdz",
                "debefdz",
                "debezdz",
                "deperdz",
                "depefdz",
                "depezdz",
            )
                @test isapprox(
                    getproperty(g, Symbol(fld)),
                    refdata["magnetic/$k/$fld"][1];
                    rtol = 1e-4,
                    atol = 1e-4,
                )
            end
            # 2nd b derivatives carry upstream FD noise (~10-20%): loose gate
            for fld in ("dbds2", "dbdst", "dbdt2")
                @test isapprox(
                    getproperty(g, Symbol(fld)),
                    refdata["magnetic/$k/$fld"][1];
                    rtol = 0.3,
                    atol = 1.0,
                )
            end
        end
        nde = Int(refdata["disp_eig/npts"][1])
        for k = 1:nde
            y = Tuple(refdata["disp_eig/$k/y"])
            @test isapprox(
                dispersion_U(RCN_PROB, y),
                refdata["disp_eig/$k/U"][1];
                rtol = 1e-6,
                atol = 1e-8,
            )
            @test isapprox(
                conversion_monitors(RCN_PROB, y).mon2,
                refdata["disp_eig/$k/mon"][2];
                rtol = 1e-6,
            )
            rhs = trajectory_rhs(RCN_PROB, collect(y))
            for i = 1:4
                @test isapprox(
                    rhs[i],
                    refdata["disp_eig/$k/trj"][i];
                    rtol = 1e-4,
                    atol = 1e-6 * maximum(abs.(refdata["disp_eig/$k/trj"])),
                )
            end
            p3 = cmod_parameters(; model = :cld3x3)
            @test isapprox(
                dispersion_U(p3, y),
                refdata["disp_eig_3x3/$k/U"][1];
                rtol = 1e-6,
                atol = 1e-8,
            )
        end
        ncg = Int(refdata["cgamma/npts"][1])
        zre = refdata["cgamma/z_re"]
        zim = refdata["cgamma/z_im"]
        gre = refdata["cgamma/g_re"]
        gim = refdata["cgamma/g_im"]
        for k = 1:ncg
            # skip the Γ poles (z = 0, −1, …): the reference encodes them as NaN
            (isnan(gre[k]) || isnan(gim[k])) && continue
            gv = cgamma(complex(zre[k], zim[k]))
            @test isapprox(real(gv), gre[k]; rtol = 1e-9, atol = 1e-12)
            @test isapprox(imag(gv), gim[k]; rtol = 1e-9, atol = 1e-12)
        end

        # ---- dispertok (det-U conversion layer) ----
        zdot = refdata["dispertok_zdot"]
        zddot = refdata["dispertok_zddot"]
        A = zdot[1] * zddot[3] + zdot[2] * zddot[4] - zdot[3] * zddot[1] - zdot[4] * zddot[2]
        eqv = zdot ./ sqrt(abs(A))
        epv = zddot ./ sqrt(abs(A))
        ndt = Int(refdata["dispertok/npts"][1])
        for k = 1:ndt
            y = Tuple(refdata["dispertok/$k/y"])
            core = HybridPlasmaPIC.Raycon._detU_core(RCN_PROB, Float64.(y); need2nd = true)
            @test isapprox(core.U, refdata["dispertok/$k/Dsp"][1]; rtol = 1e-8)
            @test isapprox(msw_dispersion(RCN_PROB, y), refdata["dispertok/$k/Msw"][1]; rtol = 1e-8)
            frq = cyclotron_frequencies(RCN_PROB, y[1], y[2])
            for i = 1:3
                @test isapprox(frq[i], refdata["dispertok/$k/Frq"][i]; rtol = 1e-9)
            end
            # dispertok 'Trj' = (−dU/dkr, −dU/dkz, +dU/dr, +dU/dz) in its sign
            # convention (opposite overall sign to the disp_eig tracer)
            trj = refdata["dispertok/$k/Trj"]
            @test isapprox(-core.dUdkr, trj[1]; rtol = 1e-6)
            @test isapprox(-core.dUdkz, trj[2]; rtol = 1e-6)
            @test isapprox(core.dUdr, trj[3]; rtol = 1e-6)
            @test isapprox(core.dUdz, trj[4]; rtol = 1e-6)
            # 'Mon' = [mon1, log|V|, η² estimate, yg1, yg3]; 'Sdl' = z + zinzst
            osc = HybridPlasmaPIC.Raycon._osculating(core, collect(eqv), collect(epv))
            mon = refdata["dispertok/$k/Mon"]
            @test isapprox(log(abs(core.V)), mon[2]; rtol = 1e-9)
            @test isapprox(osc.eta2, mon[3]; rtol = 1e-8)
            @test isapprox(y[1] + 2 * osc.zinzst[1], mon[4]; rtol = 1e-7)
            @test isapprox(y[3] + 2 * osc.zinzst[3], mon[5]; rtol = 1e-7)
            sdl = refdata["dispertok/$k/Sdl"]
            for i = 1:4
                @test isapprox(
                    y[i] + osc.zinzst[i],
                    sdl[i];
                    rtol = 1e-7,
                    atol = 1e-7 * max(abs(sdl[3]), abs(sdl[4])),
                )
            end
            # NOTE: 'Sgn' (∂U/∂ω) is deliberately NOT compared — upstream
            # carries a sign bug in dD12dom that the port fixes (see notes).
        end

        # ---- full conversion pipeline at the upstream saddle ----
        zst = Tuple(refdata["conversion/zst"])
        y0c = Tuple(refdata["conversion/y0"])
        coeff =
            HybridPlasmaPIC.Raycon._conversion_coefficients(RCN_PROB, Float64.(zst), Float64.(y0c))
        @test coeff.hyperbola_ok
        # power-iteration stopping (upstream tol 1e-4) limits agreement to ~1e-5
        @test isapprox(coeff.tau, refdata["conversion/tau"][1]; rtol = 1e-5)
        betaref = complex(refdata["conversion/beta_re"][1], refdata["conversion/beta_im"][1])
        @test isapprox(abs(coeff.beta), abs(betaref); rtol = 1e-5)
        # arg β is defined modulo π (eigenvector sign freedom)
        dphi = mod(angle(coeff.beta) - angle(betaref), π)
        @test min(dphi, π - dphi) < 1e-4
        @test coeff.transmitted_ok
        ytrs = refdata["conversion/ytrs"]
        for i = 1:4
            @test isapprox(
                coeff.transmitted[i],
                ytrs[i];
                rtol = 1e-6,
                atol = 1e-8 * max(abs(ytrs[3]), abs(ytrs[4])),
            )
        end

        # ---- antenna launch and traced ray ----
        y0ref = refdata["launch/y0"]
        y0j = launch_ray(RCN_PROB; s = 0.4, theta = 0.001, kr = -31.5, kz = 0.0)
        for i = 1:4
            @test isapprox(y0j[i], y0ref[i]; rtol = 1e-9, atol = 1e-10)
        end
        # MATLAB ode45 stopped at its conversion event; integrate to the same σ
        # (jsonencode emits matrices as nested ROW arrays → row-major decode)
        tend = refdata["ray/t"][end]
        npts = Int(refdata["ray/npts"][1])
        yref = permutedims(reshape(refdata["ray/y"], 4, npts))
        for frac in (0.5, 1.0)
            idx = clamp(round(Int, frac * npts), 2, npts)
            tr = integrate_ray(RCN_PROB, y0j, 0.0, refdata["ray/t"][idx]; detect_conversion = false)
            @test tr.status === :end_of_span
            # measured agreement vs the original ode45 run: ≲1e-5 absolute;
            # gate with an order of margin for platform variation
            @test abs(tr.y[1, end] - yref[idx, 1]) <= 1e-5      # R [m]
            @test abs(tr.y[2, end] - yref[idx, 2]) <= 1e-5      # Z [m]
            @test abs(tr.y[3, end] - yref[idx, 3]) <= 1e-3      # kR [1/m]
            @test abs(tr.y[4, end] - yref[idx, 4]) <= 1e-3      # kZ [1/m]
        end
    end
else
    @info "raycon MATLAB reference data not found — skipping RCN-REF (generate with tools/raycon_reference.m)"
end
