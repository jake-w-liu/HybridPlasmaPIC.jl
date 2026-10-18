module ConversionChecks
using Test, Random
using HybridPlasmaPIC
include("../polytropic_pressure/binomial_oracle.jl")
const P = (convert_scaled = HybridPlasmaPIC._from_exact,)
const B = BinomialOracle
const Q = Rational{BigInt}
function expected(::Type{T}, n::Int128, d::UInt64, e::Int) where {T}
    x = Q(n, d) * Q(2)^e
    return x < 0 ? -B.independent_round(T, -x) : B.independent_round(T, x)
end
println((VERSION,))
flush(stdout)
@testset "Scaled rational conversion against independent ordered float search" begin
    rng = MersenneTwister(823329)
    for T in (Float32, Float64)
        for i = 1:3000
            n = rand(rng, Int128)
            d = max(UInt64(1), rand(rng, UInt64))
            e = rand(rng, -1200:1200)
            @test isequal(P.convert_scaled(T, n, d, e), expected(T, n, d, e))
        end
        for n in (
                typemin(Int128),
                typemin(Int128) + 1,
                -Int128(1),
                Int128(0),
                Int128(1),
                typemax(Int128),
            ),
            d in (UInt64(1), UInt64(2), UInt64(3), UInt64(1) << 53, typemax(UInt64)),
            e in (
                -1200,
                -1075,
                -1074,
                -1022,
                -150,
                -149,
                -126,
                -127,
                -53,
                -1,
                0,
                1,
                53,
                127,
                128,
                1023,
                1024,
                1200,
            )

            @test isequal(P.convert_scaled(T, n, d, e), expected(T, n, d, e))
        end
    end
end
flush(stdout)
@testset "Exact boundaries, ties and extreme exponent guards" begin
    for T in (Float32, Float64)
        p = precision(T)
        minq = exponent(floatmin(T)) - p + 1
        for d in (UInt64(1), UInt64(3), UInt64(7), (UInt64(1) << 53) - 1),
            q in (Int128(0), Int128(1), Int128(2), (Int128(1) << p) - 2, (Int128(1) << p) - 1),
            e in (minq - 1, -p - 1, exponent(floatmax(T)) - p),
            sign in (-1, 1)
            # n/d = 2q+1 constructs exact half-quantum boundaries.
            n = Int128(sign) * (2q + 1) * Int128(d)
            @test isequal(P.convert_scaled(T, n, d, e), expected(T, n, d, e))
            @test isequal(P.convert_scaled(T, n + sign, d, e), expected(T, n + sign, d, e))
            @test isequal(P.convert_scaled(T, n - sign, d, e), expected(T, n - sign, d, e))
        end
        for n in (typemin(Int128), -Int128(1), Int128(1), typemax(Int128)),
            d in (UInt64(1), typemax(UInt64))

            @test P.convert_scaled(T, n, d, typemax(Int)) === copysign(T(Inf), n)
            @test P.convert_scaled(T, n, d, typemin(Int)) === copysign(zero(T), n)
        end
        @test P.convert_scaled(T, typemin(Int128), UInt64(1), minq - 128) === -zero(T)
        @test_throws ArgumentError P.convert_scaled(T, Int128(1), UInt64(0), 0)
    end
end
println("SCALED_RATIONAL_CONVERSION_NUMERICAL_COMPLETE")
flush(stdout)

end
