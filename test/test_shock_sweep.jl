# SHK-SWEEP — Phase 11 (1-D) collisionless perpendicular-shock Mach sweep.
#
# Drives the verified reflecting-wall PerpShock model over M_A ∈ {1.2, 2.0, 4.0}
# and checks the physics-grade, EOS-independent diagnostics on the supercritical
# (M_A ≥ 2) runs:
#   • a shock forms with frozen-in ratio (Bz2/B0)/n2 ≈ 1 (<5%),
#   • compression 1 < n2 < 4 (above unity, below the strong-shock fluid limit),
#   • reported shock speed is mass-flux consistent, Vs ≈ U0/(n2−1) (<8%).
# Reflected fraction and compression are reported per M_A. A 3-seed M_A=4 run
# bounds the kinetic-noise sensitivity of the measured compression (std/mean<15%).

using HybridPlasmaPIC, Test, Statistics

@testset "SHK-SWEEP perpendicular-shock Mach sweep" begin
    MAs = (1.2, 2.0, 4.0)
    results = Dict{Float64,NamedTuple}()
    for MA in MAs
        r = run_perp_shock(; MA = MA, nsteps = 700, seed = 1)
        results[MA] = r
        @info "M_A sweep" MA n2 = round(r.n2, digits = 3) Bz2 = round(r.Bz2, digits = 3) Vs =
            round(r.Vs, digits = 4) X_rh = round(r.X_rh, digits = 3) frozen =
            round(r.frozen_ratio, digits = 4) reflected = round(r.reflected_fraction, digits = 4) M_real =
            round(r.M_real, digits = 3)

        @test isfinite(r.n2) && isfinite(r.Bz2) && isfinite(r.Vs)
    end

    # supercritical cases (M_A ≥ 2) must form a real, flux-frozen shock
    for MA in (2.0, 4.0)
        r = results[MA]
        U0 = MA                                   # v_A = 1
        @test 1.0 < r.n2 < 4.0                     # compressed, below strong-shock limit
        @test abs(r.frozen_ratio - 1) < 0.05       # flux frozen to the flow
        Vs_mass = U0 / (r.n2 - 1)                  # mass-conservation prediction
        @test abs(r.Vs - Vs_mass) / Vs_mass < 0.08 # reported speed ⇔ mass conservation
        # fluid RH compresses MORE than the kinetic shock (sanity on the oracle)
        @test r.X_rh > r.n2
        @test 0.0 <= r.reflected_fraction <= 1.0
    end

    # seed sensitivity of the compression at M_A = 4 (kinetic-noise bound)
    n2s = Float64[]
    for s = 1:3
        r = run_perp_shock(; MA = 4.0, nsteps = 700, seed = s)
        push!(n2s, r.n2)
    end
    @info "M_A=4 seed spread" n2s = round.(n2s, digits = 3) mean = round(mean(n2s), digits = 3) std =
        round(std(n2s), digits = 4)
    @test std(n2s) / mean(n2s) < 0.15
end
