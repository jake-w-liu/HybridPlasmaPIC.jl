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
