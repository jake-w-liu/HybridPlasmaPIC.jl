# SHK-CONVERGE — Phase 11 (1-D) perpendicular-shock Mach sweep + convergence.
#
# Exercises the two thin research drivers on top of the verified PerpShock model:
#
#   • mach_sweep over the full M_A = 1.2, 2, 4, 6 set returns 4 results, and every
#     supercritical (M_A ≥ 2) case forms a real, flux-frozen shock:
#       frozen-in ratio (Bz2/B0)/n2 ≈ 1  (|frozen − 1| < 6%), 2 < n2 < 4,
#       and reported shock speed consistent with Vs ≈ U0/(n2−1).
#
#   • convergence_study runs the same MA=4 shock at two grid resolutions
#     (N=256, 512) and two ppc (32, 64); the downstream compression n2 agrees
#     across both refinements within ~12% (converged in Δx and in ppc).

using HybridPlasmaPIC, Test

@testset "SHK-CONVERGE Mach sweep + convergence" begin
    @testset "mach_sweep: full M_A = 1.2,2,4,6 set" begin
        res = mach_sweep(; MAs = [1.2, 2.0, 4.0, 6.0], nsteps = 600, seed = 1)
        @test length(res) == 4
        @test [r.MA for r in res] == [1.2, 2.0, 4.0, 6.0]

        for r in res
            @info "mach_sweep" MA = r.MA n2 = round(r.n2, digits = 3) frozen =
                round(r.frozen_ratio, digits = 4) Vs = round(r.Vs, digits = 4) reflected =
                round(r.reflected_fraction, digits = 4)
            @test isfinite(r.n2) && isfinite(r.Bz2) && isfinite(r.frozen_ratio)
            @test 0.0 <= r.reflected_fraction <= 1.0
        end

        # supercritical cases (M_A ≥ 2): real, flux-frozen, 2 < n2 < 4
        for r in res
            if r.MA >= 2.0
                @test 2.0 < r.n2 < 4.0
                @test abs(r.frozen_ratio - 1) < 0.06
                Vs_mass = r.MA / (r.n2 - 1)
                @test isfinite(r.Vs)
                @test abs(r.Vs - Vs_mass) / Vs_mass < 0.08
            end
        end
    end

    @testset "convergence_study: Δx and ppc convergence at M_A=4" begin
        cv = convergence_study(; MA = 4.0, Ns = (256, 512), ppcs = (32, 64), nsteps = 500)
        @info "convergence_study" n2_base = round(cv.n2_base, digits = 3) n2_fine_N =
            round(cv.n2_fine_N, digits = 3) n2_fine_ppc = round(cv.n2_fine_ppc, digits = 3) rel_grid =
            round(cv.rel_grid, digits = 4) rel_ppc = round(cv.rel_ppc, digits = 4)

        @test isfinite(cv.n2_base) && isfinite(cv.n2_fine_N) && isfinite(cv.n2_fine_ppc)
        # each run is itself a real shock (sanity)
        @test 2.0 < cv.n2_base < 4.0
        @test 2.0 < cv.n2_fine_N < 4.0
        @test 2.0 < cv.n2_fine_ppc < 4.0
        # converged: compression insensitive to Δx and to ppc within ~12%
        @test cv.rel_grid < 0.12
        @test cv.rel_ppc < 0.12
        @test cv.rel_max < 0.12
        @test cv.converged
    end
end
