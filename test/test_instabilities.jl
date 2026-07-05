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

@testset "ion_cyclotron_growth: EMIC (T_⊥>T_∥) anisotropy instability" begin
    cfg = (N = 96, Lx = 19.2, nppc = 100, nsteps = 700, dt = 0.02, seed = 1)
    u = ion_cyclotron_growth(; vth_par = 0.4, vth_perp = 1.3, cfg...)   # A=9.6 ≫ threshold
    s = ion_cyclotron_growth(; vth_par = 0.8, vth_perp = 0.9, cfg...)   # A=0.27 < threshold

    @test u.unstable_theory && !s.unstable_theory       # EMIC threshold flag (Gary 1993)
    @test u.T_anisotropy > 1 && s.T_anisotropy > 1      # both T_⊥ > T_∥, differ by margin
    @test u.ratio_max > 0.02                            # transverse δB grows above threshold
    @test s.ratio_max < 0.01                            # sub-threshold stays at noise
    @test u.ratio_max > 8 * s.ratio_max
    @test u.nsamples == 700
    # deterministic; reuses firehose_growth's guards
    @test ion_cyclotron_growth(; vth_par = 0.4, vth_perp = 1.3, cfg...).ratio_max == u.ratio_max
    @test_throws ArgumentError ion_cyclotron_growth(; vth_par = 0.4, vth_perp = 1.3, dt = 0.1)
end

@testset "weibel_growth: current-filamentation instability (full EM-PIC)" begin
    # Counter-streaming electron beams (full PIC): B_z grows from shot noise; the
    # hybrid model cannot do this (B=0 fixed point), so this uses EMPIC. Fast config.
    cfg = (N = (8, 64), L = (4π, 8π), nppc = 40, nsteps = 300, c = 3.0, dt = 0.05, seed = 1)
    u = weibel_growth(; u0 = 0.6, vth = 0.1, cfg...)   # A = (6)² = 36 ≫ A_c (unstable)
    s = weibel_growth(; u0 = 0.0, vth = 0.1, cfg...)   # A = 0 (single Maxwellian, stable)

    @test u.unstable_theory && !s.unstable_theory       # box filamentation threshold
    @test u.anisotropy > 1 && s.anisotropy == 0
    @test u.wBz_max > 1e-2                              # B_z grows to a macroscopic level
    @test s.wBz_max < 1e-3                              # no streaming ⇒ shot-noise floor
    @test u.wBz_max > 50 * s.wBz_max                    # large, unambiguous separation
    @test u.nsamples == 300 && s.nsamples == 300        # neither run blew up
    # deterministic for a fixed seed
    @test weibel_growth(; u0 = 0.6, vth = 0.1, cfg...).wBz_max == u.wBz_max
    # EM-Courant guard + input validation (both throw before running the sim)
    @test_throws ArgumentError weibel_growth(; u0 = 0.6, dt = 1.0)
    # 2-D corner CFL: an isotropic grid with c·dt·|k|max>1.9 must be rejected (a 1-D
    # min(dx) guard would wrongly pass this and silently blow up)
    @test_throws ArgumentError weibel_growth(;
        u0 = 0.6,
        N = (64, 64),
        L = (2π, 2π),
        c = 3.0,
        dt = 0.0173,
        nsteps = 1,
    )
    @test_throws ArgumentError weibel_growth(; u0 = 0.6, vth = -1.0)
    @test_throws ArgumentError weibel_growth(; u0 = -1.0)
end

@testset "weibel_growth: box-quantized filamentation threshold" begin
    # A_c = (c·k_min/ω_pe)² with k_min = 2π/L_y and ω_pe = √n₀ = 1: the default
    # c = 3, L_y = 12π give A_c = 0.25 — NOT the old A > 1 heuristic (Weibel/Fried
    # predict growth for ANY A > 0 in an infinite domain; the box quantizes k).
    # Classifier-only check: a minimal 1-step run keeps the same c, L (hence A_c).
    tiny = (N = (4, 4), L = (4π, 12π), nppc = 1, nsteps = 1, c = 3.0, dt = 0.05, seed = 1)
    bu = weibel_growth(; u0 = 0.08, vth = 0.1, tiny...)          # A = 0.64 > 0.25
    bs = weibel_growth(; u0 = 0.1 * sqrt(0.1), vth = 0.1, tiny...) # A = 0.10 < 0.25
    @test bu.anisotropy ≈ 0.64
    @test bs.anisotropy ≈ 0.1
    @test bu.unstable_theory        # was false under the A > 1 heuristic (0.25 < A < 1 band)
    @test !bs.unstable_theory       # below the box cutoff: k_min lies outside the unstable band
end

@testset "reconnection_growth: Harris-sheet tearing (2D hybrid)" begin
    # The clean ~12× sheet-vs-uniform separation needs the full 64×128 config (too slow
    # for CI) — that is validation case 32. Here a small smoke config: the runner builds
    # the pressure-balanced Harris equilibrium and tears the seeded m=1 mode.
    cfg = (Nx = 32, Ny = 64, Lx = 25.6, Ly = 25.6, nppc = 40, nsteps = 150, dt = 0.03, NB = 4)
    u = reconnection_growth(; sheet = true, cfg...)
    c = reconnection_growth(; sheet = false, cfg...)

    @test u.tearing_theory && !c.tearing_theory     # kx·λ<1 sheet unstable; no sheet stable
    @test u.nsamples == 150 && c.nsamples == 150     # neither run blew up (2-D whistler CFL ok)
    @test u.growth > 1.5                             # the sheet amplifies the m=1 tearing mode
    @test isfinite(u.m1_max) && isfinite(c.m1_max)
    @test reconnection_growth(; sheet = true, cfg...).growth == u.growth   # deterministic
    # 2-D whistler-CFL guard + equilibrium/input validation (all throw before the sim)
    @test_throws ArgumentError reconnection_growth(; dt = 0.1)
    @test_throws ArgumentError reconnection_growth(; λ = 10.0)   # Ly > 4λ violated
    @test_throws ArgumentError reconnection_growth(; Ti = -1.0)
    @test_throws ArgumentError reconnection_growth(; δ = 0.0)    # no seed ⇒ growth=Inf
end
