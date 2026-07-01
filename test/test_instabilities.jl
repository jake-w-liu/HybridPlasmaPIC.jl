# Phase-1 kinetic instabilities driven on the hybrid engine (instability_sweep.jl).
# The physics is engine-reuse (load an anisotropic distribution, run, measure the
# fluctuation growth), so the checks are: the theory-threshold flag is right, the
# unstable case grows well above the particle-noise floor, the sub-threshold case
# does not, and the runner is deterministic + input-validated.

using HybridPlasmaPIC, Test

@testset "firehose_growth: parallel firehose instability" begin
    # Fast config; strong anisotropy above threshold vs a clearly sub-threshold case.
    cfg = (N = 96, Lx = 19.2, nppc = 100, nsteps = 700, dt = 0.02, seed = 1)
    u = firehose_growth(; vth_par = 1.5, vth_perp = 0.3, cfg...)   # β_∥−β_⊥ = 4.32 > 2
    s = firehose_growth(; vth_par = 1.0, vth_perp = 0.8, cfg...)   # β_∥−β_⊥ = 0.72 < 2

    # threshold flag matches the analytic firehose condition vth_∥²−vth_⊥² > 1
    @test u.unstable_theory && !s.unstable_theory
    @test u.anisotropy > 1 && s.anisotropy < 1
    # the unstable case grows the transverse field well above the noise floor;
    # the sub-threshold case stays at noise — with a clear separation between them
    @test u.ratio_max > 0.02
    @test s.ratio_max < 0.01
    @test u.ratio_max > 8 * s.ratio_max
    @test u.nsamples == 700 && s.nsamples == 700        # neither run blew up
    # deterministic for a fixed seed
    @test firehose_growth(; vth_par = 1.5, vth_perp = 0.3, cfg...).ratio_max == u.ratio_max
    # input validation + whistler-CFL guard
    @test_throws ArgumentError firehose_growth(; vth_par = 1.5, vth_perp = 0.3, dt = 0.1)
    @test_throws ArgumentError firehose_growth(; vth_par = -1.0, vth_perp = 0.3)
    @test_throws ArgumentError firehose_growth(; vth_par = 1.5, vth_perp = 0.3, N = 4)
end
