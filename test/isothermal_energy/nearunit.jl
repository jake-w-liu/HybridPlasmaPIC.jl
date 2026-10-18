module NearUnitChecks
using HybridPlasmaPIC, Test, Random
import ..Oracle
using ..Oracle: Q, oracle_log, oracle_integral
const H = HybridPlasmaPIC
const O = Oracle
@testset "Fixed-point logarithm coefficients and signed products enclose exact rationals" begin
    for bits in (40, 69, 79), exponent in (-8, -24, -53), k = 2:12
        lo, hi = H._isothermal_coefficient_bounds(k, exponent, bits)
        exact = (iseven(k) ? 1 : -1) * Q(2)^(bits + exponent * (k - 2)) / (k * (k - 1))
        @test lo <= exact <= hi
        @test hi - lo <= 1
    end
    rng = MersenneTwister(0x6912)
    for bits in (40, 69, 79), q in (0, 15, 32, 45), repetition = 1:32
        bound = Int128(1) << bits
        lo, hi = minmax(rand(rng, -bound:bound), rand(rng, -bound:bound))
        scale = Int128(1) << q
        d = rand(rng, -scale:scale)
        lower, upper = H._isothermal_multiply_y(lo, hi, d, q, scale - 1)
        exactlo, exacthi = minmax(Q(lo) * Q(d, scale), Q(hi) * Q(d, scale))
        @test lower <= exactlo <= exacthi <= upper
        @test lower > exactlo - 1
        @test upper < exacthi + 1
    end
end

@testset "Near-unit complete integrals against independent alternating logarithms" begin
    rng = MersenneTwister(0x8713)
    for T in (Float32, Float64), C in (Float32, Float64)
        small = ldexp(one(T), -(precision(T) - 1))
        edge = ldexp(one(T), -8)
        cases = [
            T[],
            T[0, -1, 1],
            T[1+small, 1-small],
            T[1+edge, 1-edge],
            T[1+edge/2],
            T[1-edge/2],
            [one(T) + ldexp(T(rand(rng, -1024:1024)), -20) for i = 1:8],
        ]
        for n in cases, temperature in (C(-0.7), C(0), C(0.7))
            dx = (T(0.5), T(2), T(0.25))
            before = copy(n)
            bounds = O.interval(n, temperature, dx)
            @test bounds !== nothing
            reflo, refhi = oracle_integral(n, temperature, dx)
            @test bounds[1] <= reflo <= refhi <= bounds[2]
            value = O.evaluate(n, temperature, dx)
            @test value !== nothing
            @test O.rounded(value, reflo, refhi)
            @test isequal(n, before)
            @test isequal(O.evaluate(view(n, :), temperature, dx), value)
        end
        for density in (one(T) + edge, one(T) - edge)
            for (temperature, dx) in
                ((C(1), (floatmax(T), floatmax(T))), (nextfloat(zero(C)), (one(T),)))
                n = T[density]
                lo, hi = oracle_integral(n, temperature, dx)
                value = O.evaluate(n, temperature, dx)
                @test value !== nothing
                @test O.rounded(value, lo, hi)
            end
        end
        @test O.evaluate(T[2], C(1), (one(T),)) === nothing
        @test O.evaluate(T[Inf], C(1), (one(T),)) === nothing
        @test O.evaluate(T[1], C(Inf), (one(T),)) === nothing
        @test O.evaluate(T[1], C(1), (T(Inf),)) === nothing
    end
end
println("ISOTHERMAL_NEARUNIT_FOCUSED_COMPLETE")

end
