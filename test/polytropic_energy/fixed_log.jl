module FixedLogChecks
using Test, Random, HybridPlasmaPIC
const F = (
    wide_product = HybridPlasmaPIC._isothermal_fixed_wide_product,
    multiply_bounds = HybridPlasmaPIC._isothermal_fixed_multiply_bounds,
    fraction_bounds = HybridPlasmaPIC._isothermal_fixed_fraction_bounds,
    logarithm_ratio = HybridPlasmaPIC._isothermal_fixed_logarithm_ratio,
    logarithm_bounds = HybridPlasmaPIC._isothermal_fixed_logarithm_bounds,
    probe = HybridPlasmaPIC._isothermal_fixed_probe,
    FRACTION_BITS = HybridPlasmaPIC._ISOTHERMAL_FIXED_FRACTION_BITS,
    SCALE = HybridPlasmaPIC._ISOTHERMAL_FIXED_SCALE,
    MASK = HybridPlasmaPIC._ISOTHERMAL_FIXED_MASK,
    LIMB_MASK = HybridPlasmaPIC._ISOTHERMAL_FIXED_LIMB_MASK,
)
const Q = Rational{BigInt}
println((VERSION, root = pkgdir(HybridPlasmaPIC)))
flush(stdout)
@testset "Independent exact wide products and outward fixed operations" begin
    rng = MersenneTwister(1729)
    integers = vcat(
        UInt128[
            0,
            1,
            typemax(UInt64),
            UInt128(1)<<64,
            (UInt128(1)<<127)-1,
            UInt128(1)<<127,
            typemax(UInt128),
        ],
        rand(rng, UInt128, 1000),
    )
    for a in integers,
        b in (UInt128(0), UInt128(1), UInt128(1) << 64, typemax(UInt128), rand(rng, UInt128))

        low, high = F.wide_product(a, b)
        @test BigInt(low) + (BigInt(high) << 128) == BigInt(a) * BigInt(b)
    end
    for _ = 1:2000
        a = rand(rng, UInt128) % (F.SCALE + 1)
        b = rand(rng, UInt128) % (F.SCALE + 1)
        lo, hi = F.multiply_bounds(a, b)
        v = BigInt(a) * BigInt(b) // BigInt(F.SCALE)
        @test BigInt(lo) == fld(numerator(v), denominator(v))
        @test BigInt(hi) == cld(numerator(v), denominator(v))
        d = rand(rng, UInt64) | UInt64(1)
        n = rand(rng, UInt64) % d
        lo, hi = F.fraction_bounds(n, d)
        v = (BigInt(n) << 112) // BigInt(d)
        @test BigInt(lo) == fld(numerator(v), denominator(v))
        @test BigInt(hi) == cld(numerator(v), denominator(v))
    end
    @test F.multiply_bounds(F.SCALE, F.SCALE) == (F.SCALE, F.SCALE)
    @test F.fraction_bounds(typemax(UInt64), typemax(UInt64)) == (F.SCALE, F.SCALE)
    @test_throws ArgumentError F.fraction_bounds(UInt64(0), UInt64(0))
    @test_throws ArgumentError F.logarithm_ratio(UInt64(0), UInt64(0))
    @test_throws ArgumentError F.logarithm_ratio(UInt64(1), UInt64(2))
    @test 9big(2)^112 < 4 * 69 * big(3)^69
end
function mpfr_bounds(x)
    setprecision(BigFloat, 512) do
        lower = setrounding(BigFloat, RoundDown) do
            Q(log(BigFloat(x)))
        end
        upper = setrounding(BigFloat, RoundUp) do
            Q(log(BigFloat(x)))
        end
        return lower, upper
    end
end
@testset "Logarithm enclosures across exact binary input range" begin
    rng = MersenneTwister(2718)
    logtwo = F.logarithm_ratio(UInt64(1), UInt64(3))
    before = (precision(BigFloat), rounding(BigFloat))
    for T in (Float32, Float64)
        cases = T[
            nextfloat(zero(T)),
            floatmin(T),
            floatmax(T),
            prevfloat(one(T)),
            one(T),
            nextfloat(one(T)),
            T(0.5),
            T(0.75),
            T(1.5),
            T(2),
        ]
        for _ = 1:500
            bits =
                T === Float32 ? rand(rng, UInt32) & UInt32(0x7fffffff) :
                rand(rng, UInt64) & UInt64(0x7fffffffffffffff)
            x = reinterpret(T, bits)
            isfinite(x) && x > 0 && push!(cases, x)
        end
        for x in cases
            lo, hi, e = F.logarithm_bounds(x, logtwo)
            lower, upper = mpfr_bounds(x)
            @test Q(lo) * Q(2)^e <= lower <= upper <= Q(hi) * Q(2)^e
            @test hi >= lo
        end
        @test F.logarithm_bounds(one(T), logtwo) == (Int128(0), Int128(0), -112)
        for x in (zero(T), -one(T), T(Inf), T(NaN))
            @test_throws ArgumentError F.logarithm_bounds(x, logtwo)
        end
    end
    @test (precision(BigFloat), rounding(BigFloat)) == before
end
println("FIXED_LOG_COMPONENT_COMPLETE")

end
