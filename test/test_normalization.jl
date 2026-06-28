using HybridPlasmaPIC, Test

const KINDS = (:length, :time, :velocity, :magnetic, :electric, :density, :current, :pressure)

@testset "normalization SI <-> normalized" begin
    u = PlasmaUnits(n0 = 1e6, B0 = 5e-9, mi = 1.6726e-27)  # solar-wind-ish

    # Derived scales (recomputed independently from the definitions).
    vA = u.B0 / sqrt(u.mu0 * u.n0 * u.mi)
    Ωci = u.e * u.B0 / u.mi
    di = vA / Ωci

    @testset "derived accessors" begin
        @test alfven_speed(u) == vA
        @test gyrofrequency(u) == Ωci
        @test inertial_length(u) == di
        @test vA ≈ 1.087e5 rtol = 1e-2   # sanity: ~100 km/s
    end

    @testset "to_SI of unit normalized value == reference scale" begin
        @test to_SI(1.0, :velocity, u) == alfven_speed(u)
        @test to_SI(1.0, :length, u) == inertial_length(u)
        @test to_SI(1.0, :magnetic, u) == u.B0
        @test to_SI(1.0, :time, u) == 1 / gyrofrequency(u)
        @test to_SI(1.0, :density, u) == u.n0
        @test to_SI(1.0, :electric, u) == alfven_speed(u) * u.B0
        @test to_SI(1.0, :current, u) == u.B0 / (u.mu0 * inertial_length(u))
        @test to_SI(1.0, :pressure, u) == u.B0^2 / u.mu0
    end

    @testset "round-trip (scalars)" begin
        for kind in KINDS
            for x in (0.0, 1.0, -3.7, 2.5e3, 1e-9, 4.2e12)
                @test to_normalized(to_SI(x, kind, u), kind, u) ≈ x rtol = 1e-12
                @test to_SI(to_normalized(x, kind, u), kind, u) ≈ x rtol = 1e-12
            end
        end
    end

    @testset "to_SI and to_normalized are mutual inverses (scale check)" begin
        for kind in KINDS
            s_to_si = to_SI(1.0, kind, u)
            s_to_norm = to_normalized(1.0, kind, u)
            @test s_to_si * s_to_norm ≈ 1.0 rtol = 1e-14
        end
    end

    @testset "array broadcasting" begin
        xs = [0.0, 1.0, -2.0, 3.5]
        for kind in KINDS
            si = to_SI(xs, kind, u)
            @test si isa AbstractArray
            @test to_normalized(si, kind, u) ≈ xs rtol = 1e-12
        end
    end

    @testset "unknown kind throws" begin
        @test_throws ArgumentError to_SI(1.0, :bogus, u)
        @test_throws ArgumentError to_normalized(1.0, :foobar, u)
    end

    @testset "reference scales must be finite and positive" begin
        @test_throws ArgumentError PlasmaUnits(n0 = 0.0, B0 = 5e-9, mi = 1.6726e-27)
        @test_throws ArgumentError PlasmaUnits(n0 = -1.0, B0 = 5e-9, mi = 1.6726e-27)
        @test_throws ArgumentError PlasmaUnits(n0 = NaN, B0 = 5e-9, mi = 1.6726e-27)
        @test_throws ArgumentError PlasmaUnits(n0 = 1e6, B0 = 0.0, mi = 1.6726e-27)
        @test_throws ArgumentError PlasmaUnits(n0 = 1e6, B0 = -5e-9, mi = 1.6726e-27)
        @test_throws ArgumentError PlasmaUnits(n0 = 1e6, B0 = 5e-9, mi = 0.0)
        @test_throws ArgumentError PlasmaUnits(n0 = 1e6, B0 = 5e-9, mi = -1.0)
    end

    @testset "type promotion in constructor" begin
        ui = PlasmaUnits(n0 = 1, B0 = 5e-9, mi = 1.6726e-27)  # mixed Int/Float
        @test ui isa PlasmaUnits{Float64}
        @test ui.n0 == 1.0
    end
end
