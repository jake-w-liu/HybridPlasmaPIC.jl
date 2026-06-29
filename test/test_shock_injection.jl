# Sustained-shock fixes: upstream injection (the reservoir no longer depletes),
# the flux-based reflected fraction α (Leroy 1982 definition), and a
# ripple/boundary-robust shock_front locator.

using HybridPlasmaPIC, Test, Random, Statistics

@testset "upstream injection sustains the inflow reservoir" begin
    T = Float64
    N, Lx, nppc, MA, vthi = 256, 80.0, 64, 4.0, 0.5
    U0 = MA
    wp = shock_density_weight(1.0, Lx, nppc * N)
    function run(; inject)
        sh = PerpShock(N, T(Lx); Te = 0.25, γe = 5 / 3, η = 0.02, τ = U0, B0 = 1.0)
        Np = nppc * N
        ps = ParticleSet{1,T}(Np)
        rng = MersenneTwister(1)
        for p = 1:Np
            ps.x[1][p] = Lx * rand(rng)
            ps.v[1][p] = -U0 + vthi * randn(rng)
            ps.v[2][p] = vthi * randn(rng)
            ps.v[3][p] = vthi * randn(rng)
        end
        ps.weight .= wp
        init_shock!(sh, ps)
        injector =
            inject ?
            ShockInjector(
                MersenneTwister(2);
                n0 = 1.0,
                drift = U0,
                vthi = vthi,
                weight = wp,
                first_id = Np + 1,
            ) : nothing
        for _ = 1:700
            step_shock!(sh, ps, 0.02; NB = 2, injector = injector)
        end
        return mean(sh.n[sh.x.>Lx-15])           # inflow-side density
    end
    up_noinj = run(; inject = false)
    up_inj = run(; inject = true)
    @test up_noinj < 0.1                          # reservoir drains without injection
    @test 0.8 < up_inj < 1.25                     # injection holds n₁ ≈ 1 at the inflow
end

@testset "step_shock! default (no injector) is unchanged" begin
    T = Float64
    r1 = run_perp_shock(; MA = 4.0, N = 256, Lx = 80.0, nppc = 32, nsteps = 400, seed = 3)
    r2 = run_perp_shock(;
        MA = 4.0,
        N = 256,
        Lx = 80.0,
        nppc = 32,
        nsteps = 400,
        seed = 3,
        inject = false,
    )
    @test r1.n2 == r2.n2 && r1.reflected_fraction == r2.reflected_fraction && r1.xf == r2.xf
end

@testset "reflected_flux_fraction α (flux at front / upstream flux)" begin
    T = Float64
    Vs, V1, xf, win = 2.0, 6.0, 10.0, 3.0       # n1·V1 = 6
    # build a particle set: upstream incident ions (vx−Vs<0) + a reflected beam
    # (x in (xf,xf+win), vx−Vs>0). α should be the reflected flux / (n1·V1).
    function αset(wrefl, vrefl)
        ps = ParticleSet{1,T}(0)
        # reflected beam at x=11.5, one node, weight wrefl, vx = Vs+vrefl
        push!(ps.x[1], 11.5)
        push!(ps.v[1], Vs + vrefl)
        for c = 2:3
            push!(ps.v[c], 0.0)
        end
        push!(ps.weight, wrefl)
        push!(ps.id, UInt64(1))
        push!(ps.tag, UInt32(0))
        return ps
    end
    α0 = reflected_flux_fraction(αset(0.6, 1.0), xf, Vs, V1; window = win)
    @test α0 ≈ (0.6 * 1.0 / win) / (1.0 * V1)    # exact flux ratio
    # zero reflected ions ⇒ α = 0
    psup = ParticleSet{1,T}(0)
    push!(psup.x[1], 11.5)
    push!(psup.v[1], Vs - 1.0)                    # incident (vx−Vs<0)
    for c = 2:3
        push!(psup.v[c], 0.0)
    end
    push!(psup.weight, 1.0)
    push!(psup.id, UInt64(1))
    push!(psup.tag, UInt32(0))
    @test reflected_flux_fraction(psup, xf, Vs, V1; window = win) == 0.0
    # linear in reflected weight, inverse in V1
    @test reflected_flux_fraction(αset(1.2, 1.0), xf, Vs, V1; window = win) ≈ 2 * α0
    @test reflected_flux_fraction(αset(0.6, 1.0), xf, Vs, 2V1; window = win) ≈ α0 / 2
    @test_throws ArgumentError reflected_flux_fraction(αset(0.6, 1.0), xf, Vs, V1; window = 0.0)
end

@testset "run_perp_shock_rh: sustained two-state-RH shock reproduces Leroy structure" begin
    # Two-state RH initialization (Leroy 1982 setup) → a sustained perpendicular
    # shock whose downstream holds the fluid-RH compression, with an emergent
    # magnetic overshoot and a reflected fraction that rises with M_A (Leroy trend).
    cfg = (β = 1.0, N = 256, nppc = 64, nsteps = 600, Lx = 160.0, t_avg_start = 6.0)
    r4 = run_perp_shock_rh(; MA = 4.0, seed = 1, cfg...)
    r6 = run_perp_shock_rh(; MA = 6.0, seed = 1, cfg...)
    r8 = run_perp_shock_rh(; MA = 8.0, seed = 1, cfg...)
    for (r, MA) in ((r4, 4.0), (r6, 6.0), (r8, 8.0))
        @test r.nsamples > 0
        @test isapprox(r.M_real, MA; rtol = 0.02)            # frame setup is consistent
        @test isapprox(r.compression, r.X_rh; rtol = 0.05)   # downstream holds the RH state
        @test 1.1 < r.overshoot < 1.9                        # a real magnetic overshoot
    end
    @test r4.reflected_flux < r6.reflected_flux < r8.reflected_flux  # α rises with M_A (Leroy)
    # deterministic for a fixed seed
    @test run_perp_shock_rh(; MA = 6.0, seed = 1, cfg...).compression == r6.compression
    # a non-compressive Mach is rejected (the RH solver finds no shock)
    @test_throws ArgumentError run_perp_shock_rh(; MA = 0.5, cfg...)
end

@testset "run_perp_shock_leroy: §11.3 wall-less two-ended-flux shock reaches Leroy α" begin
    # Leroy 1982's actual setup: no wall, shock-REST frame, upstream inflow at x=Lx /
    # downstream thermal-reservoir outflow at x=0. Unlike the reflecting-wall model, the
    # downstream reservoir (not a specular wall) lets the energetic foot develop, so the
    # reflected fraction α reaches Leroy's published 10–23% band at high Mach.
    cfg = (β = 1.0, N = 256, nppc = 64, nsteps = 600, Lx = 200.0, t_avg_start = 6.0)
    r4 = run_perp_shock_leroy(; MA = 4.0, seed = 1, cfg...)
    r6 = run_perp_shock_leroy(; MA = 6.0, seed = 1, cfg...)
    r8 = run_perp_shock_leroy(; MA = 8.0, seed = 1, cfg...)
    for (r, MA) in ((r4, 4.0), (r6, 6.0), (r8, 8.0))
        @test r.nsamples > 0
        @test r.M_real == MA                                  # inflow Mach, exact by construction
        @test isapprox(r.compression, r.X_rh; rtol = 0.05)    # downstream holds the fluid RH state
        @test 1.2 < r.overshoot < 2.0                         # a real magnetic overshoot
    end
    @test r4.reflected_flux < r6.reflected_flux < r8.reflected_flux   # α rises with M_A (Leroy)
    @test r8.reflected_flux > 0.08            # reaches Leroy's reflected-fraction regime (wall model < 0.02)
    # deterministic for a fixed seed; non-compressive Mach rejected
    @test run_perp_shock_leroy(; MA = 6.0, seed = 1, cfg...).reflected_flux == r6.reflected_flux
    @test_throws ArgumentError run_perp_shock_leroy(; MA = 0.5, cfg...)
    @test_throws ArgumentError run_perp_shock_leroy(; MA = 6.0, window = 0.0, cfg...)
end

@testset "shock_front is robust to ripples and boundary artifacts" begin
    T = Float64
    N = 200
    x = collect(range(0.0, 100.0, length = N))
    dx = x[2] - x[1]
    # clean step front at x≈60 (Bz 3→1 going upstream), + small 2-cell ripple, +
    # a sharp boundary spike near the wall (the kind that fooled the raw locator).
    Bz = [xi < 60 ? 3.0 : 1.0 for xi in x]
    for i = 2:N-1
        Bz[i] += 0.15 * sin(2π * i / 3)           # 3-cell ripple everywhere
    end
    Bz[5] = 6.0                                    # boundary spike near the wall
    xf, w = shock_front(Bz, x)
    @test 57 < xf < 63                             # locates the macroscopic ramp…
    @test xf > 20                                  # …not the near-wall spike
    @test w > 0
    # smoothing matters: the raw single-cell gradient is fooled by the spike
    xf_raw, _ = shock_front(Bz, x; smooth = 0)
    @test xf_raw < 20                              # raw locator picks the spike
    # a clean, ripple-free front is located identically with or without smoothing
    Bzc = [xi < 40 ? 4.0 : 1.0 for xi in x]
    @test shock_front(Bzc, x)[1] ≈ shock_front(Bzc, x; smooth = 0)[1] atol = 2dx
end
