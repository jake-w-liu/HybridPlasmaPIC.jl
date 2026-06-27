# Phase-6 diagnostics: conserved totals, J·E work, anisotropic temperatures,
# distributions, spectra, pressure–strain, shock front — vs analytic/loaded.

using HybridPlasmaPIC, Test, Random, Statistics, LinearAlgebra

@testset "total momentum" begin
    T = Float64
    N = 10_000
    ps = ParticleSet{1,T}(N; m = 2.0)
    rng = MersenneTwister(1)
    load_maxwellian!(ps, rng, (0.3, -0.1, 0.2), (0.5, 0.5, 0.5))
    p = total_momentum(ps)
    @test p[1] ≈ 2.0 * sum(ps.weight .* ps.v[1])
    @test isapprox(p[1] / (N * 2.0), 0.3; atol = 0.03)   # ≈ m·N·u0_x
end

@testset "electric work ∫J·E" begin
    T = Float64
    n = (8, 8)
    g = FourierGrid(n, (1.0, 1.0))
    J = (fill(2.0, n), fill(-1.0, n), fill(0.5, n))
    E = (fill(1.0, n), fill(3.0, n), fill(-2.0, n))
    @test electric_work(J, E, g) ≈ (2 * 1 + (-1) * 3 + 0.5 * (-2)) * prod(g.dx) * prod(n)
end

@testset "diagnostic grid shape checks reject equal-length wrong geometry" begin
    T = Float64
    g1 = FourierGrid((8,), (1.0,))
    good1 = ntuple(_ -> ones(T, 8), 3)
    bad_matrix = ntuple(_ -> ones(T, 4, 2), 3)

    @test_throws DimensionMismatch electric_work(good1, bad_matrix, g1)
    @test_throws DimensionMismatch resistive_dissipation(bad_matrix, 0.1, g1)
    @test_throws DimensionMismatch jdotE_density(good1, bad_matrix)
    @test_throws DimensionMismatch power_spectrum(ones(T, 4, 2), g1)

    g2 = FourierGrid((4, 2), (1.0, 1.0))
    wrong2 = ntuple(_ -> ones(T, 2, 4), 3)
    @test_throws DimensionMismatch magnetic_energy(wrong2, g2)
    @test_throws DimensionMismatch electron_internal_energy(
        ones(T, 2, 4),
        PolytropicElectrons(1.0, 1.0, 5 / 3),
        g2,
    )
    @test_throws DimensionMismatch mode_amplitude(ones(T, 2, 4), g2, (1, 0))

    P = ntuple(_ -> ones(T, 4, 2), 6)
    B = ntuple(_ -> ones(T, 4, 2), 3)
    @test_throws DimensionMismatch temperatures_par_perp(P, ones(T, 8), B)
    @test_throws DimensionMismatch pressure_strain(ntuple(_ -> ones(T, 2, 4), 6), B, g2)
end

@testset "parallel/perp temperatures (bi-Maxwellian, B=ẑ)" begin
    T = Float64
    n = 16
    L = 2π
    g = FourierGrid((n,), (L,))
    nppc = 3000
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_uniform!(ps, MersenneTwister(1), (0.0,), (L,))
    set_density_weight!(ps, 1.0, g)
    Tx, Ty, Tz = 0.8, 0.4, 1.5
    load_maxwellian!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (sqrt(Tx), sqrt(Ty), sqrt(Tz)))
    P = ntuple(_ -> zeros(T, n), 6)
    pressure_tensor!(P, ps, g, CIC())
    nb = zeros(T, n)
    density!(nb, ps, g, CIC())
    B = (zeros(T, n), zeros(T, n), ones(T, n))      # B along z
    Tpar, Tperp = temperatures_par_perp(P, nb, B)
    @test isapprox(mean(Tpar), Tz; rtol = 0.06)            # T∥ = T_z
    @test isapprox(mean(Tperp), (Tx + Ty) / 2; rtol = 0.06) # T⊥ = (T_x+T_y)/2
end

@testset "velocity + phase-space histograms" begin
    T = Float64
    N = 50_000
    ps = ParticleSet{1,T}(N)
    load_uniform!(ps, MersenneTwister(1), (0.0,), (1.0,))
    load_maxwellian!(ps, MersenneTwister(2), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0))
    ps.weight .= 1.0
    c, counts = velocity_histogram(ps, 1; nbins = 80, vmin = -6, vmax = 7)
    @test isapprox(sum(counts), N; rtol = 0.02)            # captures ~all weight
    @test abs(c[argmax(counts)] - 0.5) < 0.3               # peak near drift
    xc, vc, h =
        phase_space_histogram(ps, 1, 1; nx = 32, nv = 32, xmin = 0, xmax = 1, vmin = -6, vmax = 7)
    @test isapprox(sum(h), N; rtol = 0.02)
end

@testset "power spectrum (single mode)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    m = 5
    k = 2π * m / L
    x = [(i - 1) * g.dx[1] for i = 1:n]
    f = cos.(k .* x)
    kk, P = power_spectrum(f, g)
    @test argmax(P) == m + 1                                # peak at mode m
    @test sum(P[setdiff(1:length(P), m + 1)]) / P[m+1] < 1e-10
end

@testset "pressure–strain (analytic)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    p0 = 0.7
    A = 0.1
    k = 2π / L
    x = [(i - 1) * g.dx[1] for i = 1:n]
    P = (fill(p0, n), fill(p0, n), fill(p0, n), zeros(T, n), zeros(T, n), zeros(T, n))
    u = (A .* sin.(k .* x), zeros(T, n), zeros(T, n))       # u_x = A sin(kx)
    pid = pressure_strain(P, u, g)
    exact = @. -p0 * A * k * cos(k * x)                     # −Pxx ∂x u_x
    @test maximum(abs, pid .- exact) < 1e-10
end

@testset "shock front position and width (tanh ramp)" begin
    T = Float64
    n = 400
    L = 100.0
    x = collect(range(0, L; length = n))
    x0 = 30.0
    w = 2.0
    bz_down = 3.0
    bz_up = 1.0
    Bz = @. bz_down + (bz_up - bz_down) * 0.5 * (1 + tanh((x - x0) / w))
    xs, width = shock_front(Bz, x)
    @test abs(xs - x0) < 2 * (L / n)                        # within a cell of x0
    @test isapprox(width, 2w; rtol = 0.05)                  # tanh full width = 2w

    @test_throws ArgumentError shock_front(T[], T[])
    @test_throws ArgumentError shock_front(T[1.0], T[2.0])
    @test_throws DimensionMismatch shock_front([1.0, 2.0], [0.0])
    @test_throws ArgumentError shock_front([1.0, 2.0], [0.0, 0.0])
    @test_throws ArgumentError shock_front([1.0, 2.0], [1.0, 0.0])
end
