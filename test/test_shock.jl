# SHK-001 — Rankine–Hugoniot solver. Residuals of all six conservation laws
# must vanish; the hydrodynamic/parallel compression is cross-checked against
# the closed-form gas-dynamic jump  X = (γ+1)M²/((γ−1)M²+2).

using HybridPlasmaPIC, Test, Random

gas_compression(M, γ) = (γ + 1) * M^2 / ((γ - 1) * M^2 + 2)

# stub RNG whose uniform draw is exactly 0.0 — probes the flux_speed a=0
# Rayleigh branch at its rand()==0 corner (reachable with probability 2⁻⁵³)
struct FluxSpeedZeroRNG <: Random.AbstractRNG end
Base.rand(::FluxSpeedZeroRNG, ::Type{Float64}) = 0.0

@testset "SHK-001 Rankine–Hugoniot" begin
    γ = 5 / 3
    up_invalid = MHDState(1.0, 2.0, 0.0, 0.2, 0.0, 1.0)

    @testset "invalid γ is rejected" begin
        for badγ in (1.0, 0.9, NaN, Inf)
            @test_throws ArgumentError rankine_hugoniot(up_invalid, badγ)
            @test_throws ArgumentError rh_branches(up_invalid, badγ)
        end
    end

    @testset "invalid μ0 and upstream states are rejected" begin
        for badμ0 in (0.0, -1.0, NaN, Inf)
            @test_throws ArgumentError rankine_hugoniot(up_invalid, γ; μ0 = badμ0)
            @test_throws ArgumentError rh_branches(up_invalid, γ; μ0 = badμ0)
        end

        for badnscan in (0, -5)
            @test_throws ArgumentError rh_branches(up_invalid, γ; nscan = badnscan)
        end

        for up_bad in (
            MHDState(0.0, 2.0, 0.0, 0.2, 0.0, 1.0),
            MHDState(-1.0, 2.0, 0.0, 0.2, 0.0, 1.0),
            MHDState(1.0, 2.0, 0.0, -0.2, 0.0, 1.0),
            MHDState(1.0, NaN, 0.0, 0.2, 0.0, 1.0),
            MHDState(1.0, 2.0, NaN, 0.2, 0.0, 1.0),
            MHDState(1.0, 2.0, 0.0, 0.2, NaN, 1.0),
            MHDState(1.0, 2.0, 0.0, 0.2, 0.0, Inf),
        )
            @test_throws ArgumentError rankine_hugoniot(up_bad, γ)
            @test_throws ArgumentError rh_branches(up_bad, γ)
        end
    end

    @testset "hydrodynamic limit" begin
        ρ, p = 1.0, 1.0
        cs = sqrt(γ * p / ρ)
        for M in (1.5, 2.0, 4.0, 6.0)
            up = MHDState(ρ, M * cs, 0.0, p, 0.0, 0.0)
            sol = rankine_hugoniot(up, γ)
            @test isapprox(sol.X, gas_compression(M, γ); rtol = 1e-8)
            @test maximum(values(sol.residuals)) < 1e-10
        end
    end

    @testset "parallel shock (B ∥ normal decouples → gas compression)" begin
        ρ, p = 1.0, 1.0
        cs = sqrt(γ * p / ρ)
        M = 3.0
        up = MHDState(ρ, M * cs, 0.0, p, 1.0, 0.0)
        sol = rankine_hugoniot(up, γ)
        @test isapprox(sol.X, gas_compression(M, γ); rtol = 1e-8)
        @test sol.down.Bt == 0.0                 # tangential field stays zero
        @test maximum(values(sol.residuals)) < 1e-10
    end

    @testset "perpendicular shock (B ⟂ normal compresses with density)" begin
        up = MHDState(1.0, 4.0, 0.0, 1.0, 0.0, 0.5)
        sol = rankine_hugoniot(up, γ)
        @test sol.X > 1
        @test isapprox(sol.down.Bt, sol.X * up.Bt; rtol = 1e-10)   # B_t ∝ ρ for Bn=0
        @test maximum(values(sol.residuals)) < 1e-10
    end

    @testset "oblique fast shock (B amplified)" begin
        up = MHDState(1.0, 4.0, 0.0, 1.0, 0.5, 0.5)
        sol = rankine_hugoniot(up, γ)
        @test sol.X > 1
        @test sol.down.Bt > up.Bt                # fast shock amplifies tangential field
        @test sol.down.Bn == up.Bn               # normal B continuous
        @test maximum(values(sol.residuals)) < 1e-10
    end

    @testset "downstream is subsonic-compressed and physical" begin
        up = MHDState(1.0, 4.0, 0.0, 1.0, 0.3, 0.6)
        sol = rankine_hugoniot(up, γ)
        @test sol.down.ρ > up.ρ                  # compression
        @test sol.down.ux < up.ux                # deceleration
        @test sol.down.p > up.p                  # heating
        @test sol.down.ρ ≈ sol.X * up.ρ
        @test sol.down.ux ≈ up.ux / sol.X        # mass conservation
    end

    @testset "branch tracking (rh_branches)" begin
        γ = 5 / 3
        # perpendicular super-fast: a single (fast) branch, residuals ~0
        up = MHDState(1.0, 4.0, 0.0, 1.0, 0.0, 0.5)
        brs = rh_branches(up, γ)
        @test length(brs) >= 1
        for b in brs
            @test maximum(values(b.residuals)) < 1e-10   # every found branch satisfies the jumps
            @test b.X > 1
        end
        # the fast (largest-X) branch matches rankine_hugoniot
        sol = rankine_hugoniot(up, γ)
        @test isapprox(maximum(b.X for b in brs), sol.X; rtol = 1e-6)
        # oblique case: branches are tracked and all satisfy the jump conditions
        up2 = MHDState(1.0, 3.0, 0.0, 1.0, 0.6, 0.6)
        brs2 = rh_branches(up2, γ)
        @test length(brs2) >= 1
        for b in brs2
            @test maximum(values(b.residuals)) < 1e-10
        end
    end

    @testset "switch-on window (Bt₁ = 0 bifurcation)" begin
        # field-aligned low-β upstream INSIDE the switch-on window: M_An² = 2.25 ∈
        # (1, (γ+1)/(γ−1) = 4), β = 0.02. The physical downstream is the switch-on
        # branch X = M_An², Bt₂ = √(17/12) from the energy jump; the gasdynamic
        # root here crosses the Alfvén point (1→4) and is non-evolutionary.
        up = MHDState(1.0, 1.5, 0.0, 0.01, 1.0, 0.0)
        Bt2_ref = sqrt(17 / 12)                     # 1.190238…
        brs = rh_branches(up, γ)
        @test length(brs) >= 2                      # gasdynamic root AND switch-on
        for b in brs
            @test maximum(values(b.residuals)) < 1e-10
            @test b.X > 1
        end
        so = only(b for b in brs if b.down.Bt > 0.1)
        @test isapprox(so.X, 2.25; rtol = 1e-6)
        @test isapprox(so.down.Bt, Bt2_ref; rtol = 1e-6)
        @test isapprox(so.down.ux, 2 / 3; rtol = 1e-6)
        @test isapprox(so.down.uy, 2 / 3 * Bt2_ref; rtol = 1e-6)
        @test isapprox(so.down.p, 0.5516666666666666; rtol = 1e-6)
        # the single-solution API selects the evolutionary (switch-on) branch
        sol = rankine_hugoniot(up, γ)
        @test isapprox(sol.X, 2.25; rtol = 1e-6)
        @test isapprox(sol.down.Bt, Bt2_ref; rtol = 1e-6)
        @test maximum(values(sol.residuals)) < 1e-10
        # outside the window (M_An² = 9 > 4) the gasdynamic branch stays the
        # unique, evolutionary solution — no spurious switch-on branch appears
        brs_out = rh_branches(MHDState(1.0, 3.0, 0.0, 0.01, 1.0, 0.0), γ)
        @test all(b.down.Bt == 0 for b in brs_out)
    end
end

@testset "flux_speed stays finite when rand() returns exactly 0" begin
    # Xoshiro Float64 draws hit exactly 0.0 with probability 2⁻⁵³; the a=0
    # Rayleigh inverse CDF must invert on 1−U ∈ (0,1] (log1p), not log(U) = −Inf.
    s0 = flux_speed(FluxSpeedZeroRNG(), 0.0, 1.0)
    @test isfinite(s0)
    @test s0 >= 0
    # the drifting (a ≠ 0) bisection branch is finite at U = 0 too
    sa = flux_speed(FluxSpeedZeroRNG(), 2.0, 1.0)
    @test isfinite(sa)
    @test sa >= 0
end
