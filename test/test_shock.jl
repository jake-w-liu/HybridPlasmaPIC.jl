# SHK-001 — Rankine–Hugoniot solver. Residuals of all six conservation laws
# must vanish; the hydrodynamic/parallel compression is cross-checked against
# the closed-form gas-dynamic jump  X = (γ+1)M²/((γ−1)M²+2).

using HybridPlasmaPIC, Test

gas_compression(M, γ) = (γ + 1) * M^2 / ((γ - 1) * M^2 + 2)

@testset "SHK-001 Rankine–Hugoniot" begin
    γ = 5 / 3
    up_invalid = MHDState(1.0, 2.0, 0.0, 0.2, 0.0, 1.0)

    @testset "invalid γ is rejected" begin
        for badγ in (1.0, 0.9, NaN, Inf)
            @test_throws ArgumentError rankine_hugoniot(up_invalid, badγ)
            @test_throws ArgumentError rh_branches(up_invalid, badγ)
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
end
