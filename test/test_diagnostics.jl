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
    @test_throws ArgumentError temperature_components(P, nb; nfloor = 0.0)
    @test_throws ArgumentError temperature_components(P, nb; nfloor = -1.0)
    @test_throws ArgumentError temperature_components(P, nb; nfloor = NaN)
    @test_throws ArgumentError temperature_components(P, nb; nfloor = Inf)
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

    extreme = ParticleSet{1,T}(2)
    extreme.x[1] .= 0.0
    offset = 1.0e154
    extreme.v[1] .= (offset, nextfloat(offset))
    Pextreme = ntuple(_ -> zeros(T, 1), 6)
    pressure_tensor!(Pextreme, extreme, g, NGP())
    mean_big = (BigFloat(extreme.v[1][1]) + BigFloat(extreme.v[1][2])) / 2
    expected_big = sum((BigFloat(v) - mean_big)^2 for v in extreme.v[1])
    @test Pextreme[1][1] == Float64(expected_big)
    @test all(isfinite, Iterators.flatten(Pextreme))

    asymmetric = ParticleSet{1,T}(2)
    asymmetric.x[1] .= 0.0
    tiny_weight = nextfloat(0.0)
    asymmetric.weight .= tiny_weight
    asymmetric.v[1] .= (0.0, 1.0e-200)
    asymmetric.v[2] .= (0.0, 1.0e300)
    Pasymmetric = ntuple(_ -> zeros(T, 1), 6)
    pressure_tensor!(Pasymmetric, asymmetric, g, NGP(); nfloor = tiny_weight)
    expected_xy = Float64(BigFloat(tiny_weight) / 2 * BigFloat(1.0e-200) * BigFloat(1.0e300))
    expected_yy = Float64(BigFloat(tiny_weight) / 2 * BigFloat(1.0e300)^2)
    @test Pasymmetric[4][1] == expected_xy
    @test Pasymmetric[2][1] == expected_yy
    @test all(isfinite, Iterators.flatten(Pasymmetric))

    floored = ParticleSet{1,T}(2; m = 2.0)
    floored.x[1] .= 0.0
    floored.weight .= tiny_weight
    floored.v[1] .= (0.0, 1.0e300)
    Pfloored = ntuple(_ -> zeros(T, 1), 6)
    density_floor = 4tiny_weight
    pressure_tensor!(Pfloored, floored, g, NGP(); nfloor = density_floor)
    expected_floored = Float64(
        BigFloat(floored.m) * (
            BigFloat(tiny_weight) * BigFloat(1.0e300)^2 -
            (BigFloat(tiny_weight) * BigFloat(1.0e300))^2 / BigFloat(density_floor)
        ),
    )
    @test Pfloored[1][1] == expected_floored

    wide = ParticleSet{1,T}(2)
    wide.x[1] .= 0.0
    wide.weight .= (1.0e183, 1.0e25)
    wide.v[1] .= (1.0e-285, -1.0e-282)
    wide.v[2] .= (1.0e-289, -1.0e31)
    wide.v[3] .= (1.0e-235, 1.0e293)
    Pwide = ntuple(_ -> zeros(T, 1), 6)
    pressure_tensor!(Pwide, wide, g, NGP(); nfloor = nextfloat(0.0))
    expected_wide = setprecision(BigFloat, 1024) do
        total_weight = sum(BigFloat, wide.weight)
        mean_velocity = ntuple(
            c ->
                sum(BigFloat(wide.weight[p]) * BigFloat(wide.v[c][p]) for p = 1:2) / total_weight,
            3,
        )
        map(((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))) do (i, j)
            Float64(
                sum(
                    BigFloat(wide.weight[p]) *
                    (BigFloat(wide.v[i][p]) - mean_velocity[i]) *
                    (BigFloat(wide.v[j][p]) - mean_velocity[j]) for p = 1:2
                ),
            )
        end
    end
    @test map(first, Pwide) == expected_wide
    @test !any(isnan, Iterators.flatten(Pwide))

    cancelling = ParticleSet{1,T}(4)
    cancelling.x[1] .= 0.0
    cancelling.weight .= 1.0
    amplitude = 1.0e154
    cancelling.v[1] .= (amplitude, -amplitude, amplitude, -amplitude)
    cancelling.v[2] .= (amplitude, -amplitude, -amplitude, amplitude)
    Pcancelling = ntuple(_ -> zeros(T, 1), 6)
    pressure_tensor!(Pcancelling, cancelling, g, NGP(); nfloor = nextfloat(0.0))
    @test map(first, Pcancelling) == (Inf, Inf, 0.0, 0.0, 0.0, 0.0)

    separated = ParticleSet{1,T}(6)
    separated.x[1] .= 0.0
    separated.weight .= 1.0
    separated.v[1] .= (amplitude, 1.0, -amplitude, -1.0, amplitude, -amplitude)
    separated.v[2] .= (amplitude, 1.0, amplitude, -1.0, -amplitude, -amplitude)
    Pseparated = ntuple(_ -> zeros(T, 1), 6)
    pressure_tensor!(Pseparated, separated, g, NGP(); nfloor = nextfloat(0.0))
    @test map(first, Pseparated) == (Inf, Inf, 0.0, 2.0, 0.0, 0.0)

    for bad_floor in (0.0, -1.0, NaN, Inf)
        P0 = map(copy, P)
        @test_throws ArgumentError pressure_tensor!(
            P,
            ps,
            g,
            NGP();
            nfloor = bad_floor,
            work,
            nbuf,
            mom,
        )
        @test all(isequal(P[c], P0[c]) for c = 1:6)
    end

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

    weights0 = copy(ps.weight)
    P0 = map(copy, P)
    @test_throws ArgumentError pressure_tensor!(P, ps, g, NGP(); work = ps.weight, nbuf, mom)
    @test ps.weight == weights0
    @test all(P[c] == P0[c] for c = 1:6)
    @test_throws ArgumentError pressure_tensor!(P, ps, g, NGP(); work, nbuf = P[1], mom)
    @test all(P[c] == P0[c] for c = 1:6)
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

@testset "histograms handle empty particle sets with explicit bounds" begin
    ps = ParticleSet{1,Float64}(0)
    @test_throws ArgumentError velocity_histogram(ps, 1; nbins = 8)
    @test_throws ArgumentError phase_space_histogram(ps, 1, 1; nx = 4, nv = 4)
    _, h = velocity_histogram(ps, 1; nbins = 8, vmin = -1.0, vmax = 1.0)
    @test length(h) == 8
    @test sum(h) == 0.0
    _, _, ph = phase_space_histogram(
        ps,
        1,
        1;
        nx = 4,
        nv = 4,
        xmin = 0.0,
        xmax = 1.0,
        vmin = -1.0,
        vmax = 1.0,
    )
    @test size(ph) == (4, 4)
    @test sum(ph) == 0.0
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
    @test_throws ArgumentError phase_space_histogram(
        ps,
        1,
        1;
        nx = 4,
        nv = 4,
        xmin = NaN,
        xmax = 1.0,
    )
    @test_throws ArgumentError phase_space_histogram(ps_bad_w, 1, 1; nx = 4, nv = 4)
end

@testset "histograms reject invalid component selectors" begin
    ps = ParticleSet{1,Float64}(1)
    ps.x[1][1] = 0.0
    ps.v[1][1] = 0.0

    @test_throws ArgumentError velocity_histogram(ps, 0)
    @test_throws ArgumentError velocity_histogram(ps, 4)
    @test_throws ArgumentError phase_space_histogram(ps, 0, 1)
    @test_throws ArgumentError phase_space_histogram(ps, 2, 1)
    @test_throws ArgumentError phase_space_histogram(ps, 1, 0)
    @test_throws ArgumentError phase_space_histogram(ps, 1, 4)
end
