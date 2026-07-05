# Dimension-parametric electromagnetic PIC. These tests cover the 2D/3D paths
# that are intentionally separate from the legacy, more specialized EMPIC1D
# regression suite.

using HybridPlasmaPIC, FFTW, Test, Random

function _empic_peakfreq(series, dt)
    Nt = length(series)
    S = fft(series)
    mag = abs.(S)
    best = 2
    bm = -1.0
    for j = 2:Nt÷2
        if mag[j+1] > bm
            bm = mag[j+1]
            best = j
        end
    end
    a = mag[best]
    b = mag[best+1]
    c = mag[best+2]
    d = a - 2b + c
    δ = d != 0 ? 0.5 * (a - c) / d : 0.0
    return 2π * (best + δ) / (Nt * dt)
end

function _empic_gauss_residual(es, rho)
    g = es.g
    rhohat = complex.(rho)
    fft!(rhohat)
    for I in CartesianIndices(rhohat)
        k2 = zero(eltype(rho))
        for d = 1:length(g.n)
            k2 += g.kvec[d][I[d]]^2
        end
        k2 == 0 && (rhohat[I] = zero(eltype(rhohat)))
    end
    ifft!(rhohat)
    rhorepr = real.(rhohat)
    divE = similar(rho)
    divergence!(divE, es.E, g)
    scale = max(maximum(abs, rhorepr), 1.0)
    return maximum(abs, divE .- rhorepr) / scale
end

@testset "EMPIC constructor validation" begin
    g2 = FourierGrid((8, 6), (2π, 3π))
    es2 = EMPIC(g2, 5; n0 = 0.75, c = 4.0, mi = 100.0)
    @test es2.n0 == 0.75
    @test es2.c == 4.0
    @test es2.mi == 100.0
    @test size(es2.E[1]) == (8, 6)
    @test all(size(es2.B[c]) == (8, 6) for c = 1:3)
    @test all(length(es2.Ep[c]) == 5 for c = 1:3)
    @test all(length(es2.mide[d]) == 5 for d = 1:2)

    g3 = FourierGrid((6, 5, 4), (2π, 3π, 4π))
    es3 = EMPIC(g3, 7; n0 = 1.25)
    @test size(es3.rho_n) == (6, 5, 4)
    @test all(size(es3.J[c]) == (6, 5, 4) for c = 1:3)
    @test all(length(es3.mide[d]) == 7 for d = 1:3)

    @test_throws ArgumentError EMPIC(g2, 5; n0 = NaN)
    @test_throws ArgumentError EMPIC(g2, 5; n0 = -1.0)
    @test_throws ArgumentError EMPIC(g2, 5; c = NaN)
    @test_throws ArgumentError EMPIC(g2, 5; c = 0.0)
    @test_throws ArgumentError EMPIC(g2, 5; mi = 0.0)
    @test_throws ArgumentError EMPIC(g2, -1)

    g4 = FourierGrid((2, 2, 2, 2), (1.0, 1.0, 1.0, 1.0))
    @test_throws ArgumentError EMPIC(g4, 1)
end

@testset "EMPIC validates species and timestep before mutation" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 2π))
    e = ParticleSet{2,T}(48; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (8, 6))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, 48; n0 = 1.0, c = 5.0)
    init_empic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 2)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    E0 = ntuple(c -> copy(es.E[c]), 3)
    B0 = ntuple(c -> copy(es.B[c]), 3)
    rho_n0 = copy(es.rho_n)
    rho_np10 = copy(es.rho_np1)
    time0 = es.time[]
    step0 = es.step[]

    @test_throws ArgumentError step_empic!(es, e, NaN)
    @test_throws ArgumentError step_empic!(es, e, -0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.rho_n == rho_n0
    @test es.rho_np1 == rho_np10
    @test es.time[] == time0
    @test es.step[] == step0

    e.q = 0.0
    @test_throws ArgumentError step_empic!(es, e, 0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.time[] == time0
    @test es.step[] == step0
end

@testset "EMPIC init preserves transverse field and enforces Gauss law" begin
    T = Float64
    nx, ny = 16, 8
    L = 2π
    g = FourierGrid((nx, ny), (L, L))
    N = nx * ny
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (nx, ny))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, N; n0 = 1.0, c = 5.0)
    seed = zeros(T, nx, ny)
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        seed[i, j] = 1e-3 * cos(x)
        es.E[2][i, j] = seed[i, j]
    end

    init_empic!(es, e)
    @test maximum(abs, es.E[2] .- seed) < 1e-12
    @test _empic_gauss_residual(es, es.rho_n) < 1e-12
end

@testset "EMPIC 2D transverse dispersion ω²=ωpe²+c²k²" begin
    T = Float64
    nx, ny = 32, 4
    L = 2π
    cc = 5.0
    n0 = 1.0
    g = FourierGrid((nx, ny), (L, L))
    nppc = 80
    N = nppc * nx * ny
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (nx * nppc, ny))
    set_density_weight!(e, n0, g)
    es = EMPIC(g, N; n0 = n0, c = cc, shape = CIC())
    m = 1
    k = 2π * m / L
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        es.E[2][i, j] = 1e-3 * cos(k * x)
    end
    init_empic!(es, e)
    dt = 0.01
    series = ComplexF64[]
    for _ = 1:1000
        step_empic!(es, e, dt)
        push!(series, mode_amplitude(es.E[2], g, (m, 0)))
    end
    ω = _empic_peakfreq(series, dt)
    ωth = sqrt(n0 + cc^2 * k^2)
    @test abs(ω - ωth) / ωth < 0.03
    @test charge_conservation_residual(es, dt) < 1e-8
    @test _empic_gauss_residual(es, es.rho_np1) < 1e-8
end

@testset "EMPIC 2D/3D spectral current correction conserves charge" begin
    T = Float64

    rng2 = MersenneTwister(11)
    g2 = FourierGrid((16, 12), (2π, 2π))
    e2 = ParticleSet{2,T}(1200; q = -1.0, m = 1.0)
    load_uniform!(e2, rng2, (0.0, 0.0), g2.L)
    set_density_weight!(e2, 1.0, g2)
    for p in eachindex(e2.weight)
        e2.v[1][p] = 0.2 * randn(rng2)
        e2.v[2][p] = 0.2 * randn(rng2)
        e2.v[3][p] = 0.1 * randn(rng2)
    end
    es2 = EMPIC(g2, nparticles(e2); n0 = 1.0, c = 5.0)
    init_empic!(es2, e2)
    maxres2 = 0.0
    for _ = 1:20
        step_empic!(es2, e2, 0.005)
        maxres2 = max(maxres2, charge_conservation_residual(es2, 0.005))
    end
    @test maxres2 < 1e-10
    @test _empic_gauss_residual(es2, es2.rho_np1) < 1e-10

    rng3 = MersenneTwister(13)
    g3 = FourierGrid((8, 6, 4), (2π, 2π, 2π))
    e3 = ParticleSet{3,T}(384; q = -1.0, m = 1.0)
    load_uniform!(e3, rng3, (0.0, 0.0, 0.0), g3.L)
    set_density_weight!(e3, 1.0, g3)
    for p in eachindex(e3.weight)
        e3.v[1][p] = 0.15 * randn(rng3)
        e3.v[2][p] = 0.15 * randn(rng3)
        e3.v[3][p] = 0.15 * randn(rng3)
    end
    es3 = EMPIC(g3, nparticles(e3); n0 = 1.0, c = 4.0)
    init_empic!(es3, e3)
    maxres3 = 0.0
    for _ = 1:10
        step_empic!(es3, e3, 0.003)
        maxres3 = max(maxres3, charge_conservation_residual(es3, 0.003))
    end
    @test maxres3 < 1e-10
    @test _empic_gauss_residual(es3, es3.rho_np1) < 1e-10
end

@testset "EMPIC 2D mobile-ion path" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 2π))
    N = prod(g.n)
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0, 0.0), g.L, g.n)
    load_lattice!(ions, (0.0, 0.0), g.L, g.n)
    set_density_weight!(e, 1.0, g)
    set_density_weight!(ions, 1.0, g)
    es = EMPIC(g, N; mobile = true, mi = 25.0, c = 4.0)
    @test_throws ArgumentError init_empic!(es, e)

    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.01)
    @test charge_conservation_residual(es, 0.01) < 1e-12
    @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
    @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))

    badions = ParticleSet{2,T}(N; q = 0.0, m = 25.0)
    load_lattice!(badions, (0.0, 0.0), g.L, g.n)
    set_density_weight!(badions, 1.0, g)
    @test_throws ArgumentError step_empic!(es, e, badions, 0.01)
end

@testset "EMPIC resizes electron gather buffers when the electron count grows" begin
    T = Float64
    g = FourierGrid((16, 16), (2π, 2π))
    N = 64
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0, 0.0), g.L, (8, 8))
    load_lattice!(ions, (0.0, 0.0), g.L, (8, 8))
    es = EMPIC(g, N; mobile = true, mi = 25.0, c = 8.0)
    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.005)
    # grow the electron population (as ionize_mcc! secondaries would); the next step must not
    # DimensionMismatch on stale-sized electron gather buffers (they are ion-symmetric now).
    xe = ParticleSet{2,T}(5; q = -1.0, m = 1.0)
    for d = 1:2, k = 1:5
        xe.x[d][k] = 0.1k
    end
    HybridPlasmaPIC.append_particles!(e, xe)
    @test HybridPlasmaPIC.nparticles(e) == N + 5
    step_empic!(es, e, ions, 0.005)
    @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))
end

@testset "EMPIC-order: step_empic! is 2nd-order in dt (leapfrog v + B priming)" begin
    # The multi-dim EM PIC primes both species' velocities and the seeded B on the first step
    # (same fix as EMPIC1D). Cold 2D longitudinal Langmuir self-convergence on Ex — rate ≈ 2.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.02
    function ex(nsteps, Tf)
        n = 16
        pd = 80
        N = pd * pd
        g = FourierGrid((n, n), (T(L), T(L)))
        e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
        ions = ParticleSet{2,T}(N; q = 1.0, m = 100.0)
        load_lattice!(e, (0.0, 0.0), g.L, (pd, pd))
        load_lattice!(ions, (0.0, 0.0), g.L, (pd, pd))
        set_density_weight!(e, 1.0, g)
        set_density_weight!(ions, 1.0, g)
        for p = 1:N
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
        end
        es = EMPIC(g, N; mobile = true, mi = 100.0, c = 8.0)
        init_empic!(es, e, ions)
        for _ = 1:nsteps
            step_empic!(es, e, ions, Tf / nsteps)
        end
        abs(mode_amplitude(es.E[1], g, (1, 0)))
    end
    Tf = 1.0
    seq = (40, 80, 160, 320)
    v = [ex(ns, Tf) for ns in seq]
    r1 = log2(abs(v[1] - v[2]) / abs(v[2] - v[3]))
    r2 = log2(abs(v[2] - v[3]) / abs(v[3] - v[4]))
    @test r1 > 1.6
    @test r2 > 1.6
end

@testset "EMPIC ↔ EMPIC1D spectral equivalence incl. Nyquist current exclusion" begin
    # The continuity correction zeroes the pure-Nyquist current modes (as
    # _esirkepov_Jx! / _deposit_Jy_at_midpoints! do in EMPIC1D), so the two 1D
    # solvers must agree to roundoff on EVERY Fourier mode and Ex must never
    # accumulate a grid-Nyquist sawtooth. Pre-fix: the raw |Ĵx_Nyq| survived
    # only in EMPIC (1.5e-4 after one step here), Ex's Nyquist mode
    # random-walked to a few % of max|Ex|, and the solvers drifted apart by
    # ~2.5% within 200 steps.
    T = Float64
    n = 32
    L = 2π
    N = 1280
    dt = 0.01
    g1 = FourierGrid((n,), (L,))
    g2 = FourierGrid((n,), (L,))
    rng = MersenneTwister(1234)
    e1 = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_uniform!(e1, rng, (0.0,), (L,))
    set_density_weight!(e1, 1.0, g1)
    for p = 1:N
        e1.v[1][p] = 0.1 * randn(rng)
        e1.v[2][p] = 0.1 * randn(rng)     # exercise the transverse Jy channel
    end
    e2 = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    e2.x[1] .= e1.x[1]
    for c = 1:3
        e2.v[c] .= e1.v[c]
    end
    e2.weight .= e1.weight
    es1 = EMPIC1D(g1, N; n0 = 1.0, c = 5.0, shape = CIC())
    es2 = EMPIC(g2, N; n0 = 1.0, c = 5.0, shape = CIC())
    init_empic!(es1, e1)
    init_empic!(es2, e2)

    nyq(f) = abs(fft(complex.(f))[length(f)÷2+1]) / length(f)

    # one step: current spectra agree on ALL modes; Nyquist current is zero
    step_empic!(es1, e1, dt)
    step_empic!(es2, e2, dt)
    @test maximum(abs.(fft(complex.(es1.Jx)) .- fft(complex.(es2.J[1])))) / n < 1e-14
    @test maximum(abs.(fft(complex.(es1.Jy)) .- fft(complex.(es2.J[2])))) / n < 1e-14
    @test nyq(es2.J[1]) < 1e-15
    @test nyq(es2.J[2]) < 1e-15

    # 300 steps: solvers stay equivalent to roundoff and the Nyquist mode of
    # Ex stays at machine zero (measured 9e-15 / 5.5e-16 with the fix; the
    # pre-fix values 2.5e-2 / ≳1e-3 fail these bounds by many orders)
    nyqmax = 0.0
    for _ = 2:300
        step_empic!(es1, e1, dt)
        step_empic!(es2, e2, dt)
        nyqmax = max(nyqmax, nyq(es2.E[1]) / maximum(abs.(es2.E[1])))
    end
    @test maximum(abs.(es1.Ex .- es2.E[1])) / maximum(abs.(es1.Ex)) < 1e-12
    @test maximum(abs.(es1.Ey .- es2.E[2])) / maximum(abs.(es1.Ey)) < 1e-12
    @test nyqmax < 1e-13
    @test charge_conservation_residual(es2, dt) < 1e-10
end
