using HybridPlasmaPIC, Test

# Uses the independent CGLEnergyRangeOracle from test_cgl_energy_range.jl.
@testset "CGL finite-looking overflow midpoint and decisive tails" begin
    O = CGLEnergyRangeOracle
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(4),))
        n = zeros(T, 4)
        B = (ones(T, 4), zeros(T, 4), zeros(T, 4))
        c = CGLElectrons(one(T), zero(T), one(T), one(T))
        quarter = (floatmax(T) - prevfloat(floatmax(T))) / T(4)
        for tail in (prevfloat(quarter), quarter, nextfloat(quarter))
            n .= (floatmax(T), quarter, tail, zero(T))
            @test isfinite(O.ordinary(n, B, c, g.dx))
            @test O.rounded(electron_internal_energy(n, B, c, g), O.interval(n, B, c, g.dx)...)
        end
        q = nextfloat(zero(T))
        n .= (one(T), ldexp(one(T), -precision(T)), q, zero(T))
        fill!(B[1], one(T))
        c = CGLElectrons(T(0.5), zero(T), one(T), one(T))
        @test O.ordinary(n, B, c, g.dx) === T(0.5)
        @test electron_internal_energy(n, B, c, g) === nextfloat(T(0.5))
        @test O.rounded(electron_internal_energy(n, B, c, g), O.interval(n, B, c, g.dx)...)
    end
end

@testset "CGL complete signed populations across pressure and conversion overflow" begin
    O = CGLEnergyRangeOracle
    for T in (Float32, Float64), C in (Float32, Float64), sign in (-one(C), one(C))
        g = FourierGrid((4,), (T(4),))
        large = ldexp(one(T), T === Float32 ? 100 : 600)
        n = T[large, -large, 1, 0]
        B = (ones(T, 4), zeros(T, 4), zeros(T, 4))
        c = CGLElectrons{C}(zero(C), C(2) * sign)
        @test isnan(O.ordinary(n, B, c, g.dx))
        @test electron_internal_energy(n, B, c, g) === T(sign)
        @test O.rounded(electron_internal_energy(n, B, c, g), O.interval(n, B, c, g.dx)...)
    end
end

@testset "CGL complete subnormal pressures with overflowing volume" begin
    O = CGLEnergyRangeOracle
    for T in (Float32, Float64)
        large = ldexp(one(T), exponent(floatmax(T)))
        g = FourierGrid((2, 2, 2), (large, large, large))
        n = fill(nextfloat(zero(T)), g.n)
        B = (ones(T, g.n), zeros(T, g.n), zeros(T, g.n))
        c = CGLElectrons(zero(T), large, one(T), one(T))
        @test isinf(prod(g.dx))
        @test isnan(O.ordinary(n, B, c, g.dx))
        value = electron_internal_energy(n, B, c, g)
        @test isfinite(value) && value > 0
        @test O.rounded(value, O.interval(n, B, c, g.dx)...)
        fill!(n, one(T))
        foreach(a -> fill!(a, zero(T)), B)
        cold = CGLElectrons(one(T), zero(T), one(T), one(T))
        @test isnan(O.ordinary(n, B, cold, g.dx))
        @test electron_internal_energy(n, B, cold, g) === zero(T)
        fill!(n, zero(T))
        @test electron_internal_energy(n, B, c, g) === zero(T)
    end
end
