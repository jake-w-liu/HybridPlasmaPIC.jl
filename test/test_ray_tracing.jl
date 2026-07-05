# RAY-001..009 — WKB ray tracing for the hybrid wave branches
# (src/diagnostics/ray_tracing.jl).
#
# The dispersion core is verified against the independent HYB-006 eigenvalue
# oracle (shares no code path), group velocities against the oracle's numerical
# dω/dk, and the tracer against exact local invariants: straight rays in uniform
# media, ω and transverse-k conservation in stratified media, and the pointwise
# relation k(x) = ω/c(x) for the perpendicular fast wave in a 1-D density wave.

using HybridPlasmaPIC, Test, Random, LinearAlgebra

if !@isdefined(HybridDispersionOracle)
    include("oracles/hybrid_dispersion_oracle.jl")
end
using .HybridDispersionOracle

@testset "RAY-001 dispersion roots match the HYB-006 oracle" begin
    # closed forms first: parallel cold R/L (whistler / ion-cyclotron)
    for K in (0.3, 1.0, 2.5)
        f = hybrid_wave_frequencies((K, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0))
        ωR = K^2 / 2 + K * sqrt(1 + K^2 / 4)
        ωL = -K^2 / 2 + K * sqrt(1 + K^2 / 4)
        @test isapprox(f[3], ωR; rtol = 1e-10)
        @test isapprox(f[2], ωL; rtol = 1e-10)
        @test f[1] <= 1e-10
    end
    # perpendicular fast magnetosonic ω = k√(vA²+cs²) (isothermal electrons)
    for (K, Te) in ((1.0, 0.0), (1.7, 0.36), (2.0, 0.49))
        f = hybrid_wave_frequencies((K, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); Te)
        @test isapprox(f[3], K * sqrt(1 + Te); rtol = 1e-10)
    end
    # random oblique warm sweep vs the eigenvalue oracle (γe = γi = 1 makes the
    # combined sound speed constant, matching the oracle's cs parameter)
    rng = MersenneTwister(7)
    compared = 0
    for _ = 1:150
        k = (1.5 * randn(rng), 1.5 * randn(rng), 1.5 * randn(rng))
        B = (randn(rng), randn(rng), 1.0 + randn(rng))
        n = 0.4 + 2.0 * rand(rng)
        Te = 0.8 * rand(rng)
        Ti = 0.4 * rand(rng)
        sum(abs2, B) < 0.05 && continue
        sum(abs2, k) < 1e-4 && continue
        fo = HybridDispersionOracle.dispersion_frequencies(
            k,
            B;
            cs = sqrt(Te + Ti),
            n0 = n,
            rho0 = n,
        )
        length(fo) == 3 || continue
        f = hybrid_wave_frequencies(k, n, B; Te, γe = 1.0, Ti, γi = 1.0)
        compared += 1
        for b = 1:3
            @test isapprox(f[b], fo[b]; rtol = 1e-7, atol = 1e-8)
        end
    end
    @test compared > 100
end

@testset "RAY-002 dispersion function consistency and derivatives" begin
    rng = MersenneTwister(11)
    for _ = 1:25
        k = (1.2 * randn(rng), 1.2 * randn(rng), 1.2 * randn(rng))
        B = (0.3 * randn(rng), 0.3 * randn(rng), 1.0 + 0.3 * randn(rng))
        n = 0.5 + rand(rng)
        Te = 0.5 * rand(rng)
        γe = 1.0 + rand(rng)
        sum(abs2, k) < 1e-3 && continue
        sum(abs2, B) < 0.3 && continue
        f = hybrid_wave_frequencies(k, n, B; Te, γe)
        # the solved branch frequencies are roots of the (factored) D
        offroot = abs(hybrid_wave_dispersion(f[3] + 0.5, k, n, B; Te, γe))
        for b = 1:3
            f[b] > 1e-8 || continue
            @test abs(hybrid_wave_dispersion(f[b], k, n, B; Te, γe)) <= 1e-9 * (offroot + 1.0)
        end
        # complex-step derivatives agree with central differences of the public D
        ω = f[3]
        d = HybridPlasmaPIC._dispersion_derivs(
            ω,
            Float64.(k),
            Float64(n),
            Float64.(B),
            Te,
            γe,
            0.0,
            5 / 3,
        )
        h = 1e-6
        fd(g, a) = (g(a + h) - g(a - h)) / (2h)
        @test isapprox(d.Dω, fd(a -> hybrid_wave_dispersion(a, k, n, B; Te, γe), ω); rtol = 1e-5)
        @test isapprox(d.Dn, fd(a -> hybrid_wave_dispersion(ω, k, a, B; Te, γe), n); rtol = 1e-5)
        for j = 1:3
            kj = a -> hybrid_wave_dispersion(ω, ntuple(i -> i == j ? a : k[i], 3), n, B; Te, γe)
            Bj = a -> hybrid_wave_dispersion(ω, k, n, ntuple(i -> i == j ? a : B[i], 3); Te, γe)
            @test isapprox(d.Dk[j], fd(kj, k[j]); rtol = 1e-5, atol = 1e-6)
            @test isapprox(d.DB[j], fd(Bj, B[j]); rtol = 1e-5, atol = 1e-6)
        end
    end
end

@testset "RAY-003 group velocity matches the oracle dω/dk" begin
    for (khat, B, Te, branches) in (
        ((1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.0, (2, 3)),          # parallel cold: L, R
        ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.25, (3,)),           # perpendicular warm: fast
        ((0.8, 0.0, 0.6), (0.0, 0.0, 1.0), 0.16, (1, 2, 3)),      # oblique warm: all
    )
        for K in (0.4, 1.3), b in branches
            n = 1.2
            kvec = (K * khat[1], K * khat[2], K * khat[3])
            f = hybrid_wave_frequencies(kvec, n, B; Te, γe = 1.0)
            f[b] > 1e-6 || continue
            gv = wave_group_velocity(kvec, n, B; branch = (:slow, :intermediate, :fast)[b], Te)
            @test isapprox(gv.ω, f[b]; rtol = 1e-12)
            # the oracle indexes only its POSITIVE roots: shift past our ω≈0 branches
            bo = b - count(<=(1e-8), f[1:b-1])
            vo = HybridDispersionOracle.group_velocity(
                bo,
                khat,
                K,
                B;
                cs = sqrt(Te),
                n0 = n,
                rho0 = n,
            )
            vg_along = gv.vg[1] * khat[1] + gv.vg[2] * khat[2] + gv.vg[3] * khat[3]
            @test isapprox(vg_along, vo; rtol = 1e-5, atol = 1e-7)
        end
    end
end

@testset "RAY-004 uniform media give straight rays with constant k" begin
    k0 = (0.6, 0.2, 0.3)
    gv = wave_group_velocity(k0, 1.0, (0.9, 0.1, 0.4); branch = :fast, Te = 0.09)
    # analytic
    med = AnalyticRayMedium((x, y, z) -> 1.0, (x, y, z) -> (0.9, 0.1, 0.4); Te = 0.09)
    r = trace_ray(med, (0.5, -1.0, 2.0), k0; branch = :fast, dt = 0.05, nsteps = 200)
    @test r.status === :ok
    @test length(r.t) == 201
    @test isapprox(r.ω, gv.ω; rtol = 1e-12)
    for i = 1:3
        @test maximum(abs.(r.k[i, :] .- k0[i])) <= 1e-10
        @test isapprox(r.x[i, end], (0.5, -1.0, 2.0)[i] + 10.0 * gv.vg[i]; atol = 1e-8)
        @test maximum(abs.(r.vg[i, :] .- gv.vg[i])) <= 1e-10
    end
    @test maximum(r.residual) <= 1e-10
    # grid (1-D, uniform fields: spectral gradients vanish to round-off)
    g = HybridPlasmaPIC.SpectralOperators.FourierGrid((32,), (12.0,))
    med2 = GridRayMedium(g, fill(1.0, 32), (fill(0.9, 32), fill(0.1, 32), fill(0.4, 32)); Te = 0.09)
    r2 = trace_ray(med2, (0.5, -1.0, 2.0), k0; branch = :fast, dt = 0.05, nsteps = 200)
    @test r2.status === :ok
    for i = 1:3
        @test maximum(abs.(r2.k[i, :] .- k0[i])) <= 1e-9
        @test isapprox(r2.x[i, end], r.x[i, end]; atol = 1e-8)
    end
end

@testset "RAY-005 stratified perpendicular fast wave: k(x) = ω/c(x) exactly" begin
    L = 20.0
    nfun = x -> 1.0 + 0.3 * sin(2π * x / L)
    # cold: c(x) = vA(x) = 1/√n(x)
    med = AnalyticRayMedium((x, y, z) -> nfun(x), (x, y, z) -> (0.0, 0.0, 1.0))
    r = trace_ray(med, (0.0, 0.0, 0.0), (0.7, 0.0, 0.0); branch = :fast, dt = 0.01, nsteps = 1500)
    @test r.status === :ok
    @test isapprox(r.ω, 0.7; rtol = 1e-12)          # n(0) = 1, vA = 1
    for m = 1:25:length(r.t)
        @test isapprox(r.k[1, m], r.ω * sqrt(nfun(r.x[1, m])); rtol = 2e-6)
    end
    @test maximum(abs.(r.k[2, :])) <= 1e-12          # transverse k exactly conserved
    @test maximum(abs.(r.k[3, :])) <= 1e-12
    @test maximum(r.residual) <= 1e-8
    # ω is conserved: the local branch frequency at the endpoint equals launch ω
    nend = nfun(r.x[1, end])
    fend = hybrid_wave_frequencies((r.k[1, end], 0.0, 0.0), nend, (0.0, 0.0, 1.0))
    @test isapprox(fend[3], r.ω; rtol = 1e-6)

    # warm polytropic electrons: c(x)² = vA²(x) + γe Te n^{γe−1}
    Te, γe = 0.08, 2.0
    medw = AnalyticRayMedium((x, y, z) -> nfun(x), (x, y, z) -> (0.0, 0.0, 1.0); Te, γe)
    rw = trace_ray(medw, (0.0, 0.0, 0.0), (0.7, 0.0, 0.0); branch = :fast, dt = 0.01, nsteps = 1200)
    @test rw.status === :ok
    for m = 1:25:length(rw.t)
        nx = nfun(rw.x[1, m])
        c = sqrt(1.0 / nx + γe * Te * nx^(γe - 1))
        @test isapprox(rw.k[1, m], rw.ω / c; rtol = 2e-6)
    end

    # grid medium reproduces the same invariant to interpolation accuracy
    nx = 256
    g = HybridPlasmaPIC.SpectralOperators.FourierGrid((nx,), (L,))
    xg = [(i - 1) * L / nx for i = 1:nx]
    medg = GridRayMedium(g, nfun.(xg), (zeros(nx), zeros(nx), fill(1.0, nx)))
    rg = trace_ray(medg, (0.0, 0.0, 0.0), (0.7, 0.0, 0.0); branch = :fast, dt = 0.01, nsteps = 1500)
    @test rg.status === :ok
    for m = 1:50:length(rg.t)
        @test isapprox(rg.k[1, m], rg.ω * sqrt(nfun(rg.x[1, m])); rtol = 1e-3)
    end
end

@testset "RAY-006 x-stratified medium conserves transverse wavevector" begin
    L = 16.0
    med = AnalyticRayMedium(
        (x, y, z) -> 1.0 + 0.2 * sin(2π * x / L),
        (x, y, z) -> (0.8, 0.1 + 0.15 * cos(2π * x / L), 1.0);
        Te = 0.05,
    )
    k0 = (0.4, 0.3, 0.2)
    r = trace_ray(med, (1.0, 0.0, 0.0), k0; branch = :fast, dt = 0.02, nsteps = 800)
    @test r.status === :ok
    @test maximum(abs.(r.k[2, :] .- k0[2])) <= 1e-9   # ∂/∂y = ∂/∂z = 0 ⇒ ky, kz const
    @test maximum(abs.(r.k[3, :] .- k0[3])) <= 1e-9
    @test maximum(r.residual) <= 1e-6
end

@testset "RAY-007 hybrid_wavenumbers round trip" begin
    rng = MersenneTwister(23)
    found = 0
    for _ = 1:40
        khat = (randn(rng), randn(rng), randn(rng))
        sum(abs2, khat) < 1e-3 && continue
        B = (0.2 * randn(rng), 0.2 * randn(rng), 1.0)
        n = 0.6 + rand(rng)
        Te = 0.3 * rand(rng)
        K = 0.3 + 1.5 * rand(rng)
        nrm = sqrt(sum(abs2, khat))
        kvec = (K * khat[1] / nrm, K * khat[2] / nrm, K * khat[3] / nrm)
        f = hybrid_wave_frequencies(kvec, n, B; Te)
        for b = 1:3
            f[b] > 1e-3 || continue
            sols = hybrid_wavenumbers(f[b], khat, n, B; Te)
            hit = findall(
                s ->
                    isapprox(s.kmag, K; rtol = 1e-8) &&
                        s.branch === (:slow, :intermediate, :fast)[b],
                sols,
            )
            @test length(hit) == 1
            found += 1
            # every returned solution satisfies the dispersion relation
            for s in sols
                Dref = abs(hybrid_wave_dispersion(f[b] + 0.3, s.k, n, B; Te))
                @test abs(hybrid_wave_dispersion(f[b], s.k, n, B; Te)) <= 1e-7 * (Dref + 1.0)
            end
        end
    end
    @test found > 40
end

@testset "RAY-008 validation and termination" begin
    @test_throws ArgumentError hybrid_wave_frequencies((1.0, 0.0, 0.0), 0.0, (1.0, 0.0, 0.0))
    @test_throws ArgumentError hybrid_wave_frequencies((1.0, 0.0, 0.0), -1.0, (1.0, 0.0, 0.0))
    @test_throws ArgumentError hybrid_wave_frequencies((1.0, 0.0, 0.0), 1.0, (0.0, 0.0, 0.0))
    @test_throws ArgumentError hybrid_wave_frequencies((NaN, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0))
    @test_throws ArgumentError hybrid_wave_frequencies(
        (1.0, 0.0, 0.0),
        1.0,
        (1.0, 0.0, 0.0);
        Te = -0.1,
    )
    @test_throws ArgumentError hybrid_wave_frequencies(
        (1.0, 0.0, 0.0),
        1.0,
        (1.0, 0.0, 0.0);
        γe = 0.0,
    )
    @test_throws ArgumentError hybrid_wavenumbers(0.0, (1.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0))
    @test_throws ArgumentError hybrid_wavenumbers(0.5, (0.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0))
    @test_throws ArgumentError wave_group_velocity(
        (1.0, 0.0, 0.0),
        1.0,
        (1.0, 0.0, 0.0);
        branch = :bogus,
    )
    # slow branch is degenerate (ω = 0) for cold perpendicular propagation
    @test_throws ArgumentError wave_group_velocity(
        (1.0, 0.0, 0.0),
        1.0,
        (0.0, 0.0, 1.0);
        branch = :slow,
    )
    @test_throws ArgumentError AnalyticRayMedium((x, y, z) -> 1.0, (x, y, z) -> (1.0, 0, 0); h = 0)

    g = HybridPlasmaPIC.SpectralOperators.FourierGrid((16,), (8.0,))
    ones16 = fill(1.0, 16)
    @test_throws DimensionMismatch GridRayMedium(g, fill(1.0, 8), (ones16, ones16, ones16))
    @test_throws ArgumentError GridRayMedium(g, fill(0.0, 16), (ones16, ones16, ones16))
    badB = fill(NaN, 16)
    @test_throws ErrorException GridRayMedium(g, ones16, (badB, ones16, ones16))

    med = AnalyticRayMedium((x, y, z) -> 1.0, (x, y, z) -> (1.0, 0.0, 0.0))
    @test_throws ArgumentError trace_ray(
        med,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0);
        dt = 0.1,
        nsteps = 5,
    )
    @test_throws ArgumentError trace_ray(
        med,
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0);
        dt = 0.0,
        nsteps = 5,
    )
    @test_throws ArgumentError trace_ray(
        med,
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0);
        dt = 0.1,
        nsteps = 0,
    )
    @test_throws ArgumentError trace_ray(
        med,
        (NaN, 0.0, 0.0),
        (1.0, 0.0, 0.0);
        dt = 0.1,
        nsteps = 5,
    )
    @test_throws ArgumentError trace_ray(
        med,
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0);
        branch = :bogus,
        dt = 0.1,
        nsteps = 5,
    )
    # perpendicular cold slow branch: ω = 0 at launch
    medperp = AnalyticRayMedium((x, y, z) -> 1.0, (x, y, z) -> (0.0, 0.0, 1.0))
    @test_throws ArgumentError trace_ray(
        medperp,
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0);
        branch = :slow,
        dt = 0.1,
        nsteps = 5,
    )

    # early termination: the medium turns invalid (NaN density) downstream
    medbad = AnalyticRayMedium((x, y, z) -> x < 3.0 ? 1.0 : NaN, (x, y, z) -> (1.0, 0.0, 0.0))
    rbad =
        trace_ray(medbad, (0.0, 0.0, 0.0), (0.8, 0.0, 0.0); branch = :fast, dt = 0.1, nsteps = 200)
    @test rbad.status === :invalid_medium
    @test length(rbad.t) < 201
    @test size(rbad.x, 2) == length(rbad.t)

    # early termination: residual drift beyond residual_max (grossly unstable dt)
    medvar = AnalyticRayMedium((x, y, z) -> 1.0 + 0.3 * sin(x), (x, y, z) -> (0.0, 0.0, 1.0))
    rres =
        trace_ray(medvar, (0.0, 0.0, 0.0), (0.7, 0.0, 0.0); branch = :fast, dt = 8.0, nsteps = 60)
    @test rres.status === :residual
    @test length(rres.t) < 61
    @test rres.residual[end] > 1e-2 || isnan(rres.residual[end])
end

@testset "RAY-009 dimension parametricity: 1D ≡ y-invariant 2D ≡ yz-invariant 3D" begin
    L = 18.0
    nx = 96
    xg = [(i - 1) * L / nx for i = 1:nx]
    nprof = [1.0 + 0.25 * sin(2π * x / L) for x in xg]
    Bz = fill(1.0, nx)
    zx = zeros(nx)
    g1 = HybridPlasmaPIC.SpectralOperators.FourierGrid((nx,), (L,))
    m1 = GridRayMedium(g1, nprof, (zx, zx, Bz); Te = 0.04)

    ny = 6
    g2 = HybridPlasmaPIC.SpectralOperators.FourierGrid((nx, ny), (L, 6.0))
    rep2(v) = repeat(v, 1, ny)
    m2 = GridRayMedium(g2, rep2(nprof), (rep2(zx), rep2(zx), rep2(Bz)); Te = 0.04)

    nz = 4
    g3 = HybridPlasmaPIC.SpectralOperators.FourierGrid((nx, ny, nz), (L, 6.0, 4.0))
    rep3(v) = repeat(v, 1, ny, nz)
    m3 = GridRayMedium(g3, rep3(nprof), (rep3(zx), rep3(zx), rep3(Bz)); Te = 0.04)

    x0 = (2.0, 1.0, 1.0)
    k0 = (0.5, 0.2, 0.1)
    r1 = trace_ray(m1, x0, k0; branch = :fast, dt = 0.02, nsteps = 400)
    r2 = trace_ray(m2, x0, k0; branch = :fast, dt = 0.02, nsteps = 400)
    r3 = trace_ray(m3, x0, k0; branch = :fast, dt = 0.02, nsteps = 400)
    @test r1.status === r2.status === r3.status === :ok
    @test isapprox(r1.ω, r2.ω; rtol = 1e-12)
    @test isapprox(r1.ω, r3.ω; rtol = 1e-12)
    @test maximum(abs.(r1.x .- r2.x)) <= 1e-9
    @test maximum(abs.(r1.k .- r2.k)) <= 1e-9
    @test maximum(abs.(r1.x .- r3.x)) <= 1e-9
    @test maximum(abs.(r1.k .- r3.k)) <= 1e-9
end
