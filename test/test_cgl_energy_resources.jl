using HybridPlasmaPIC, Test

@testset "CGL ordinary zero cells retain output and zero allocations" begin
    O = CGLEnergyRangeOracle
    for T in (Float32, Float64),
        C in (Float32, Float64),
        kind in (:perpendicular_null, :zero_density)

        g = FourierGrid((256,), (one(T),))
        n = ones(T, 256)
        B = (ones(T, 256), zeros(T, 256), zeros(T, 256))
        B[1][1] = zero(T)
        c = CGLElectrons(one(C), kind === :perpendicular_null ? zero(C) : one(C), one(C), one(C))
        kind === :zero_density && (n[1] = zero(T))
        before = map(copy, (n, B...))
        @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
        @test O.allocated(n, B, c, g) == 0
        @test all(isequal.(before, (n, B...)))
    end
    for T in (Float32, Float64), density in (T(Inf), T(-Inf), T(NaN))
        g = FourierGrid((16,), (one(T),))
        n = ones(T, 16)
        n[1] = density
        B = (ones(T, 16), zeros(T, 16), zeros(T, 16))
        B[1][1] = zero(T)
        c = CGLElectrons(one(T), zero(T), one(T), one(T))
        # The existing zero-field perpendicular convention ignores this cell.
        @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
    end
end

@testset "CGL positive bound retains large counts without overflow" begin
    Q = Rational{BigInt}
    for T in (Float32, Float64), count in (1, 2, 16, 1 << 22, 1 << 24, typemax(Int))
        spacing = one(T) / T(count)
        @test HybridPlasmaPIC._cgl_energy_positive_certified(one(T), (spacing,), count)
        # Two unit positive components in each of count cells give this exact
        # total. No huge array is represented by this count-arithmetic test.
        integral = 2Q(count) * Q(spacing)
        @test Q(floatmin(T)) <= integral < Q(floatmax(T))
    end
    for T in (Float32, Float64)
        @test !HybridPlasmaPIC._cgl_energy_positive_certified(zero(T), (one(T),), 1)
        @test !HybridPlasmaPIC._cgl_energy_positive_certified(one(T), (zero(T),), 1)
        @test !HybridPlasmaPIC._cgl_energy_positive_certified(one(T), (one(T),), 0)
        @test !HybridPlasmaPIC._cgl_energy_positive_certified(floatmax(T), (one(T),), 1)
        @test !HybridPlasmaPIC._cgl_energy_positive_certified(floatmin(T), (one(T),), 1)
    end
end

@testset "CGL positive grid above the signed-certificate count limit" begin
    T = Float32
    count = 1 << 22
    g = FourierGrid((count,), (one(T),))
    n = ones(T, count)
    B = (ones(T, count), zeros(T, count), zeros(T, count))
    c = CGLElectrons(one(T), one(T), one(T), one(T))
    value, absolute, normal, active, largest, nonnegative =
        HybridPlasmaPIC._cgl_energy_probe(n, B, c, g.dx, Val(true))
    @test normal && active && nonnegative
    @test !HybridPlasmaPIC._cgl_energy_normal_certified(value, absolute, count, T)
    @test HybridPlasmaPIC._cgl_energy_positive_certified(largest, g.dx, count)
    expected = CGLEnergyRangeOracle.ordinary(n, B, c, g.dx)
    @test expected === T(1.5)
    @test electron_internal_energy(n, B, c, g) === expected
    @test CGLEnergyRangeOracle.allocated(n, B, c, g) == 0
end
