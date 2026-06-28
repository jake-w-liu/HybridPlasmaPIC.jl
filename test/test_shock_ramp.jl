# SHK-RAMP — §11.3 / Phase-10: initial tanh-ramp generator + ramp-width and
# box-length sensitivity scans.
#
# Checks:
#   • initial_ramp! sets sh.Bz to the analytic tanh profile (machine precision),
#   • the CIC-deposited density reproduces the upstream/downstream endpoints
#     (n1 far upstream, n2 far downstream) within particle-statistics noise,
#   • ramp_width_scan returns finite measured widths,
#   • box_length_scan compressions agree across the two largest boxes (<10%).

using HybridPlasmaPIC, Test, Statistics, Random

@testset "SHK-RAMP initial ramp + sensitivity scans" begin

    @testset "finite-state runtime guard" begin
        @test HybridPlasmaPIC._require_all_finite("Bz", [1.0, 2.0], "guard test") === nothing
        @test_throws ErrorException HybridPlasmaPIC._require_all_finite(
            "Bz",
            [1.0, Inf],
            "guard test",
        )
        @test_throws ErrorException HybridPlasmaPIC._require_all_finite(
            "Bz",
            [1.0, NaN],
            "guard test",
        )
    end

    @testset "initial_ramp! field + density profile" begin
        T = Float64
        N = 512
        Lx = 120.0
        x_ramp = 60.0
        width = 4.0
        n1, n2 = 1.0, 3.0
        B1, B2 = 1.0, 3.0
        nppc = 256
        Np = nppc * N

        sh = PerpShock(N, Lx; Te = 0.125, γe = 5 / 3, η = 0.02, τ = 0.0, B0 = B1)
        ps = ParticleSet{1,T}(Np)
        rng = MersenneTwister(7)
        initial_ramp!(sh, ps, x_ramp, width, n1, n2, B1, B2; rng = rng)

        # 1. Bz matches the analytic tanh to ~1e-12
        Bz_exact = [B1 + (B2 - B1) * 0.5 * (1 - tanh((xi - x_ramp) / width)) for xi in sh.x]
        maxerr = maximum(abs.(sh.Bz .- Bz_exact))
        @info "ramp Bz error" maxerr
        @test maxerr < 1e-12

        # 2. density endpoints via mode-resolved (slab) mean. sh.n was deposited
        # by init_shock! inside initial_ramp!. Far upstream (x ≫ x_ramp) → n1,
        # far downstream (x ≪ x_ramp) → n2. Use slabs well clear of the ramp and
        # the SBP boundary nodes.
        up_mask = (sh.x .> x_ramp + 6 * width) .& (sh.x .< Lx - 3.0)
        dn_mask = (sh.x .> 3.0) .& (sh.x .< x_ramp - 6 * width)
        @test any(up_mask) && any(dn_mask)
        n_up = mean(sh.n[up_mask])
        n_dn = mean(sh.n[dn_mask])
        @info "ramp density endpoints" n_up n_dn n1 n2
        @test abs(n_up - n1) / n1 < 0.05    # upstream endpoint
        @test abs(n_dn - n2) / n2 < 0.05    # downstream endpoint
        # and density is monotone-ish: downstream mean > upstream mean
        @test n_dn > n_up

        # all finite
        @test all(isfinite, sh.Bz)
        @test all(isfinite, sh.n)
    end

    @testset "ramp_width_scan returns finite widths" begin
        @test_throws ArgumentError ramp_width_scan(; widths = (2.0,), nsteps = -1)

        res = ramp_width_scan(;
            widths = (2.0, 4.0, 8.0),
            N = 256,
            Lx = 120.0,
            nsteps = 40,
            nppc = 48,
            seed = 1,
        )
        @test length(res) == 3
        for r in res
            @info "ramp_width_scan" width0 = r.width0 xf = round(r.xf, digits = 2) wmeas =
                round(r.width_measured, digits = 3) n2 = round(r.n2_meas, digits = 3)
            @test isfinite(r.width_measured)
            @test r.width_measured > 0
            @test isfinite(r.xf)
            @test isfinite(r.n2_meas)
        end
    end

    @testset "box_length_scan compression insensitivity" begin
        @test_throws ArgumentError box_length_scan(; Lxs = (80.0,), N0 = 0)
        @test_throws ArgumentError box_length_scan(; Lxs = (80.0,), N0 = -1)
        @test_throws ArgumentError box_length_scan(; Lxs = (80.0,), Lx0 = 0.0)
        @test_throws ArgumentError box_length_scan(; Lxs = (80.0,), Lx0 = Inf)
        @test_throws ArgumentError box_length_scan(; Lxs = (80.0,), Lx0 = NaN)
        @test_throws ArgumentError box_length_scan(; Lxs = (NaN,), N0 = 64, nsteps = 0, nppc = 1)
        @test_throws ArgumentError box_length_scan(; Lxs = (Inf,), N0 = 64, nsteps = 0, nppc = 1)

        res = box_length_scan(;
            Lxs = (80.0, 120.0, 160.0),
            MA = 4.0,
            N0 = 256,
            Lx0 = 120.0,
            nsteps = 450,
            nppc = 48,
            seed = 1,
        )
        @test length(res) == 3
        for r in res
            @info "box_length_scan" Lx = r.Lx N = r.N n2 = round(r.n2, digits = 3) Bz2 =
                round(r.Bz2, digits = 3) Vs = round(r.Vs, digits = 4)
            @test isfinite(r.n2) && r.n2 > 1
        end
        # the two largest boxes should agree in compression within 10%
        n2_big = res[end].n2
        n2_mid = res[end-1].n2
        rel = abs(n2_big - n2_mid) / n2_mid
        @info "box-length compression insensitivity" n2_mid n2_big rel
        @test rel < 0.10
    end
end
