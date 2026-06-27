# test_phase_final.jl — remaining non-GPU/MPI Phase-11/13 items:
# semi-implicit integrator, integrator comparison (whistler + shock), four-
# spacecraft timing, load-imbalance metric, reference reproduction, upstream
# turbulence, 3-D restart, 3-D campaign + 1D/3D comparison. Measured tolerances.

using Random, LinearAlgebra

@testset "Semi-implicit CN integrator (conservative + unconditionally stable)" begin
    # CN multiplier modulus is exactly 1 (energy-conserving) for any ω, dt.
    for (ω, dt) in ((0.3, 0.05), (5.0, 0.7), (40.0, 1.3))
        @test abs(abs(cn_multiplier(ω, dt)) - 1) < 1e-12
    end
    # phase is the (2,2)-Padé of e^{-iωdt}: error O((ωdt)³).
    ω, dt = 0.7, 0.1
    @test abs(angle(cn_multiplier(ω, dt)) - (-ω * dt)) < (ω * dt)^3

    # full-spectrum stiff whistler: CN conserves energy, explicit Euler blows up.
    rc = run_whistler(; method = :cn, n = 128, dt = 0.6, nsteps = 200)
    @test abs(rc.energy_ratio - 1) < 1e-8
    re = run_whistler(; method = :euler, n = 128, dt = 0.6, nsteps = 200)
    @test !(isfinite(re.energy_ratio) && re.energy_ratio < 1.5)   # unstable

    # integrator comparison: on a resolved low-k mode CN and Euler agree.
    c = compare_integrators_whistler(; n = 128, nsteps = 200, dt_resolved = 0.01, dt_stiff = 0.6)
    @test c.agree_resolved < 0.03
    @test abs(c.cn_ratio_stiff - 1) < 1e-6
    @test !(isfinite(c.euler_ratio_stiff) && c.euler_ratio_stiff < 1.5)

    @test_throws ArgumentError run_whistler(; nsteps = -1)
    @test_throws ArgumentError run_whistler(; n = 0)
    @test_throws ArgumentError run_whistler(; L = 0.0)
    @test_throws ArgumentError run_whistler(; L = Inf)
    @test_throws ArgumentError run_whistler(; n0 = 0.0)
    @test_throws ArgumentError run_whistler(; n0 = Inf)
    @test_throws ArgumentError run_whistler(; dt = NaN)
    @test_throws ArgumentError run_whistler(; band = -1)
    @test_throws ArgumentError compare_integrators_whistler(; nsteps = -1)
    @test_throws ArgumentError compare_integrators_whistler(; band_resolved = -1)
end

@testset "Four-spacecraft timing" begin
    # exact recovery of a known planar boundary normal & speed
    n̂0 = (0.6, 0.48, 0.64)
    nn = sqrt(sum(abs2, n̂0))
    n̂ = n̂0 ./ nn
    V = 2.3
    pos = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    times = ntuple(i -> 5.0 + sum(n̂ .* pos[i]) / V, 4)
    r = four_spacecraft_timing(pos, times)
    @test abs(r.speed - V) < 1e-9
    @test norm(collect(r.normal) .- collect(n̂)) < 1e-9

    # Degenerate zero-slowness timing: identical crossing times imply undefined
    # normal and phase speed, so the API should return invalid markers rather
    # than a misleading infinite speed.
    rdeg = four_spacecraft_timing(pos, (1.0, 1.0, 1.0, 1.0))
    @test isnan(rdeg.speed)
    @test all(isnan, Tuple(rdeg.normal))
    @test all(isnan, Tuple(rdeg.slowness))

    # crossing_time linear interpolation: level 0.5 between (t=0,v=0) and (t=1,v=1)
    @test crossing_time([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], 0.5) ≈ 0.5 atol = 1e-12
    # level 1.0 between (t=1,v=0.5) and (t=2,v=2.0): t = 1 + (0.5/1.5) = 1.3333…
    @test crossing_time([0.0, 1.0, 2.0], [0.0, 0.5, 2.0], 1.0) ≈ 4 / 3 atol = 1e-12
    # exact level hits at the final sample are valid crossings, not "not found"
    @test crossing_time([0.0, 1.0], [0.2, 1.0], 1.0) ≈ 1.0 atol = 1e-12
    @test crossing_time([0.0, 1.0, 2.0], [0.0, 0.5, 1.0], 1.0) ≈ 2.0 atol = 1e-12
    @test isnan(crossing_time([0.0, 1.0], [0.0, 0.5], 2.0))
    @test_throws DimensionMismatch crossing_time([0.0], [0.0, 1.0], 0.5)
    @test_throws DimensionMismatch crossing_time([0.0, 1.0], [0.0], 0.5)
    @test_throws ArgumentError four_spacecraft_traces(nsteps = 0)
    @test_throws ArgumentError four_spacecraft_traces(nsteps = -1)

    # synthetic traces from a real 3-D shock recover normal ≈ x̂
    fs = four_spacecraft_traces(; MA = 3.0)
    @test all(isfinite, fs.crossings)
    @test abs(fs.normal[1]) > 0.85          # dominant x normal (measured 0.99)
    @test isfinite(fs.speed) && fs.speed > 0
end

@testset "Load-imbalance metric" begin
    @test load_imbalance(fill(10, 8)) == 1.0
    @test load_imbalance([100, 1, 1, 1, 1, 1, 1, 1]) > 5
    @test load_imbalance(Int[]) == 1.0
    @test load_imbalance(zeros(Int, 5)) == 1.0
    @test tile_loads([1, 2, 3, 4], 2) == [3, 7]
    @test sum(tile_loads(collect(1:10), 3)) == 55

    # clustered particles ⇒ imbalance > 1; uniform ⇒ near 1
    rng = MersenneTwister(3)
    g = FourierGrid((10,), (10.0,))
    psu = ParticleSet{1,Float64}(2000)
    load_uniform!(psu, rng, (0.0,), (10.0,))
    @test particle_load_imbalance(psu, g; ntiles = 5).imbalance < 1.5
    psc = ParticleSet{1,Float64}(2000)
    psc.x[1] .= 2.0 .* rand(rng, 2000)       # all in first 20% of the box
    @test particle_load_imbalance(psc, g; ntiles = 5).imbalance > 2.0
end

@testset "Reference comparison + established-shock reproduction" begin
    @test compare_to_reference((; a = 1.0, b = 2.05), (; a = 1.0, b = 2.0); rtol = 0.05).pass
    @test !compare_to_reference((; a = 1.0, b = 2.5), (; a = 1.0, b = 2.0); rtol = 0.05).pass
    @test !compare_to_reference((; a = 1.0), (; b = 2.0); rtol = 0.05).pass    # missing field
    re = reproduce_established_shock(; MA = 3.0, N = 256, nsteps = 500)
    @test re.pass
    @test abs(re.measured.frozen_ratio - 1) < 0.06
end

@testset "Shock integrator comparison (:rk4 vs :cn) + upstream turbulence" begin
    common = (;
        MA = 3.0,
        nx = 40,
        ny = 8,
        nz = 8,
        Lx = 70.0,
        Ly = 10.0,
        Lz = 10.0,
        nppc = 8,
        nsteps = 350,
        dt = 0.03,
    )
    a = run_perp_shock3d(; common..., field_method = :rk4)
    b = run_perp_shock3d(; common..., field_method = :cn)
    @test abs(a.n2 - b.n2) / a.n2 < 0.02                 # integrators agree on the shock
    @test 0.9 < b.frozen_ratio < 1.1

    tb = run_perp_shock3d(; common..., db_turb = 0.2)    # upstream turbulence: stable
    @test all(isfinite, tb.sh.B[3])
    @test tb.n2 > 1.5
end

@testset "3-D shock restart bitmatch covers full state" begin
    sh = PerpShock3D(8, 4, 4, 10.0, 4.0, 4.0)
    ps = ParticleSet{3,Float64}(3)
    load_uniform!(ps, MersenneTwister(11), (0.0, 0.0, 0.0), (10.0, 4.0, 4.0))
    load_maxwellian!(ps, MersenneTwister(12), (0.1, 0.2, 0.3), (0.01, 0.02, 0.03))
    sh2 = deepcopy(sh)
    ps2 = deepcopy(ps)
    @test HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)

    ps2.x[2][1] = nextfloat(ps2.x[2][1])
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
    ps2 = deepcopy(ps)
    ps2.v[2][1] = nextfloat(ps2.v[2][1])
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
    ps2 = deepcopy(ps)
    ps2.weight[1] = nextfloat(ps2.weight[1])
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
    ps2 = deepcopy(ps)
    ps2.id[1] += UInt64(1)
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
    ps2 = deepcopy(ps)
    ps2.tag[1] += UInt32(1)
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)

    ps2 = deepcopy(ps)
    sh2 = deepcopy(sh)
    sh2.sbp.H[1] = nextfloat(sh2.sbp.H[1])
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
    sh2 = deepcopy(sh)
    sh2.u[2][1, 1, 1] = one(eltype(sh2.u[2]))
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
    sh2 = deepcopy(sh)
    sh2.n[1, 1, 1] = one(eltype(sh2.n))
    @test !HybridPlasmaPIC._shock3d_restart_bitmatch(sh, ps, sh2, ps2)
end

@testset "3-D shock restart + campaign + 1D/3D comparison" begin
    @test_throws ArgumentError production_3d_case(; nsteps_pre = -1)
    @test_throws ArgumentError production_3d_case(; nsteps_post = -1)
    @test_throws ArgumentError shock_campaign_3d(; MAs = (4.0,), seeds = ())

    pc = production_3d_case(;
        MA = 4.0,
        nx = 40,
        ny = 8,
        nz = 8,
        nsteps_pre = 120,
        nsteps_post = 120,
        dt = 0.03,
    )
    @test pc.restart_bitmatch                            # checkpoint → identical restart
    @test pc.pass

    camp = shock_campaign_3d(;
        MAs = (4.0, 6.0),
        seeds = (1, 2),
        nx = 36,
        ny = 8,
        nz = 8,
        Lx = 70.0,
        Ly = 10.0,
        Lz = 10.0,
        nppc = 6,
        nsteps = 320,
        dt = 0.03,
    )
    @test length(camp) == 2
    for cse in camp
        @test cse.n2_mean > 1.5
        @test cse.robust                                 # cross-seed: not a single noisy run
    end

    cd = compare_dims_shock(; MA = 3.0)
    @test cd.frozen_consistent                           # matched physics, both flux-frozen
    @test cd.both_compress
end
