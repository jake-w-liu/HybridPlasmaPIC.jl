# Phase-6 pressure-tensor / temperature diagnostics, and a 3D smoke test that
# exercises the 3D operator and integrator code paths (the ultimate target
# dimension) for the first time.

using HybridPlasmaPIC, Test, Random, Statistics

@testset "pressure tensor / temperatures (bi-Maxwellian)" begin
    T = Float64
    n = 16
    L = 2π
    g = FourierGrid((n,), (L,))
    nppc = 2000
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_uniform!(ps, MersenneTwister(1), (0.0,), (L,))
    set_density_weight!(ps, 1.0, g)
    Tx, Ty, Tz = 1.0, 0.5, 0.3                       # m=1 ⇒ vth_c = √T_c
    load_maxwellian!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (sqrt(Tx), sqrt(Ty), sqrt(Tz)))
    P = ntuple(_ -> zeros(T, n), 6)
    pressure_tensor!(P, ps, g, CIC())
    nb = zeros(T, n)
    density!(nb, ps, g, CIC())
    Tc = temperature_components(P, nb)
    @test isapprox(mean(Tc[1]), Tx; rtol = 0.06)     # T_x = P_xx/n
    @test isapprox(mean(Tc[2]), Ty; rtol = 0.06)
    @test isapprox(mean(Tc[3]), Tz; rtol = 0.06)
    @test abs(mean(P[4])) < 0.05                     # off-diagonal P_xy ≈ 0
    @test abs(mean(P[5])) < 0.05
    @test abs(mean(P[6])) < 0.05
end

@testset "pressure tensor exact centered moment and scratch reuse" begin
    T = Float64
    g = FourierGrid((1,), (1.0,))
    ps = ParticleSet{1,T}(2; m = 2.0)
    ps.x[1] .= (0.0, 0.0)
    ps.weight .= 1.0
    ps.v[1] .= (-1.0, 3.0)
    ps.v[2] .= (2.0, 4.0)
    ps.v[3] .= (5.0, 1.0)

    P = ntuple(_ -> zeros(T, 1), 6)
    work = Vector{T}(undef, nparticles(ps))
    nbuf = zeros(T, 1)
    mom = ntuple(_ -> zeros(T, 1), 3)
    pressure_tensor!(P, ps, g, NGP(); work, nbuf, mom)
    expected = (16.0, 4.0, 16.0, 8.0, -16.0, -8.0)
    @test map(p -> p[1], P) == expected
    @test nbuf[1] == 2.0
    @test map(m -> m[1], mom) == (2.0, 6.0, 6.0)

    pressure_tensor!(P, ps, g, NGP(); work, nbuf, mom)
    @test (@allocated pressure_tensor!(P, ps, g, NGP(); work, nbuf, mom)) <= 128
    @test_throws DimensionMismatch pressure_tensor!(
        P,
        ps,
        g,
        NGP();
        work = Vector{T}(undef, 1),
        nbuf,
        mom,
    )
    @test_throws DimensionMismatch pressure_tensor!(P, ps, g, NGP(); work, nbuf = zeros(T, 2), mom)
    @test_throws DimensionMismatch pressure_tensor!(
        P,
        ps,
        g,
        NGP();
        work,
        nbuf,
        mom = (zeros(T, 2), zeros(T, 1), zeros(T, 1)),
    )
    @test_throws DimensionMismatch pressure_tensor!(
        (zeros(T, 2), P[2], P[3], P[4], P[5], P[6]),
        ps,
        g,
        NGP();
        work,
        nbuf,
        mom,
    )
end

@testset "3D smoke: operators + integrator" begin
    T = Float64
    nc = (8, 8, 8)
    L = (2π, 2π, 2π)
    g = FourierGrid(nc, L)

    # z-invariant 3D field ⇒ ∂z = 0, ∂x/∂y match the 2D gradient
    f2 = randn(MersenneTwister(1), 8, 8)
    f3 = repeat(reshape(f2, 8, 8, 1), 1, 1, 8)
    gr = ntuple(_ -> zeros(T, nc...), 3)
    gradient!(gr, f3, g)
    @test maximum(abs, gr[3]) < 1e-12

    # integrator runs a few 3D steps and stays finite
    counts = (16, 16, 16)                            # 4096 = 8 ppc per cell
    N = prod(counts)
    ps = ParticleSet{3,T}(N)
    load_lattice!(ps, (0.0, 0.0, 0.0), L, counts)
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    fill!(st.fields.B[3], 1.0)
    init!(st, ps)
    for _ = 1:10
        step!(st, ps, 0.02; NB = 2)
    end
    @test all(isfinite, st.fields.B[3])
    @test abs(mean(st.fields.B[3]) - 1.0) < 1e-10    # mean B_z conserved (curl has no k=0)
end

@testset "histograms handle degenerate ranges without crashing" begin
    # Regression: a single particle (min==max) gave dv=0 ⇒ floor(Int,NaN) crash.
    ps = ParticleSet{1,Float64}(1)
    ps.x[1][1] = 3.0
    ps.v[1][1] = 0.5
    _, h = velocity_histogram(ps, 1; nbins = 8)
    @test sum(h) ≈ ps.weight[1]
    _, _, ph = phase_space_histogram(ps, 1, 1; nx = 4, nv = 4)
    @test sum(ph) ≈ ps.weight[1]
end

@testset "histograms reject non-finite data and bounds" begin
    ps = ParticleSet{1,Float64}(2)
    ps.x[1] .= [0.0, 1.0]
    ps.v[1] .= [0.5, 1.0]
    ps.v[2] .= 0.0
    ps.v[3] .= 0.0

    ps_bad_v = deepcopy(ps)
    ps_bad_v.v[1][1] = NaN
    @test_throws ArgumentError velocity_histogram(ps_bad_v, 1; nbins = 8)
    @test_throws ArgumentError velocity_histogram(ps, 1; nbins = 8, vmin = NaN)

    ps_bad_w = deepcopy(ps)
    ps_bad_w.weight[1] = NaN
    @test_throws ArgumentError velocity_histogram(ps_bad_w, 1; nbins = 8)

    ps_bad_x = deepcopy(ps)
    ps_bad_x.x[1][1] = NaN
    @test_throws ArgumentError phase_space_histogram(ps_bad_x, 1, 1; nx = 4, nv = 4)
    @test_throws ArgumentError phase_space_histogram(ps, 1, 1; nx = 4, nv = 4, xmin = NaN, xmax = 1.0)
    @test_throws ArgumentError phase_space_histogram(ps_bad_w, 1, 1; nx = 4, nv = 4)
end
