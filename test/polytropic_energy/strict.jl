module StrictChecks
using Test, Random, HybridPlasmaPIC

const S = (
    bounds = HybridPlasmaPIC._polytropic_scaled_sum_state,
    dyadic_equal = HybridPlasmaPIC._polytropic_dyadic_equal,
    rounded = HybridPlasmaPIC._polytropic_round_sum,
)
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
println((VERSION,))
flush(stdout)
@testset "Strict flags enclose the entire interval with unchanged numeric bounds" begin
    rng = MersenneTwister(773492)
    for count = 0:8, bits in (1, 2, 24, 53, 69, 96), trial = 1:8
        intervals = Tuple{Q,Q,BigInt}[]
        for _ = 1:count
            lo = Q(rand(rng, 0:20), rand(rng, 1:9))
            hi = lo + Q(rand(rng, 0:7), rand(rng, 1:7))
            push!(intervals, (lo, hi, BigInt(rand(rng, -25:25))))
        end
        original = H._polytropic_scaled_sum_bounds(intervals, bits)
        lo, hi, e, ls, us = S.bounds(intervals, bits)
        @test (lo, hi, e) == original
        lower = Q(lo) * Q(2)^Int(e)
        upper = Q(hi) * Q(2)^Int(e)
        exact_lower = sum(x -> x[1] * Q(2)^Int(x[3]), intervals; init = Q(0))
        exact_upper = sum(x -> x[2] * Q(2)^Int(x[3]), intervals; init = Q(0))
        @test ls ? lower < exact_lower : lower <= exact_lower
        @test us ? upper > exact_upper : upper >= exact_upper
    end
    @test_throws ArgumentError S.bounds([(Q(2), Q(1), BigInt(0))], 10)
    @test_throws ArgumentError S.bounds(Tuple{Q,Q,BigInt}[], 0)
end
@testset "Scaled equality never constructs a large power of two" begin
    for c in BigInt[0, 1, 2, 3, 17, 1<<60], e = -100:17:100
        value = Q(c) * Q(2)^e
        @test S.dyadic_equal(c, BigInt(e), value)
        @test !S.dyadic_equal(c, BigInt(e), value + 1)
    end
    huge = BigInt(10)^300
    @test !S.dyadic_equal(BigInt(1), huge, Q(1))
    @test !S.dyadic_equal(BigInt(1), -huge, Q(1))
    @test !S.dyadic_equal(BigInt(1), BigInt(0), Q(1, 3))
end
flush(stdout)
@testset "Strict endpoints decide tiny positive terms and downward bounds" begin
    huge = BigInt(10)^300
    for T in (Float32, Float64),
        left in T[0, floatmin(T), 0.5, 1, prevfloat(floatmax(T)), floatmax(T)]

        right = nextfloat(left)
        boundary = H._polytropic_midpoint_after(left)
        _, mantissa, e = H._polytropic_normalized_factor(boundary)
        positive = ((mantissa, mantissa, e), (mantissa, mantissa, e - huge))
        state = S.bounds(positive, precision(T) + 16)
        @test state[4] && state[5]
        @test S.rounded(T, state) === right
        # A dyadic point strictly below the midpoint, at a much finer scale
        # than the summation precision, has a strict upper quantization bound.
        below = mantissa - Q(1, BigInt(1) << 200)
        state = S.bounds(((below, below, e),), precision(T) + 16)
        @test state[5]
        @test S.rounded(T, state) === left
        @test S.rounded(T, S.bounds(((mantissa, mantissa, e),), precision(T) + 16)) ===
              H._from_exact(T, boundary)
    end
    intervals = ((Q(1), Q(2), BigInt(0)),)
    @test S.rounded(Float64, S.bounds(intervals, 69)) === nothing
end
@testset "Tiny-term resources remain measured at extreme stored exponents" begin
    for magnitude in (BigInt(2)^52, BigInt(2)^54, BigInt(10)^300)
        a = Q(1) + Q(1, BigInt(1) << 26)
        b = Q(1) + Q(1, BigInt(1) << 27)
        boundary = a * b
        lower = 1.0 + 3 * 2.0^-27
        expected = nextfloat(lower)
        inputs = ((boundary, boundary, BigInt(0)), (boundary, boundary, -magnitude))
        S.rounded(Float64, S.bounds(inputs, 69))
        trial = @timed S.rounded(Float64, S.bounds(inputs, 69))
        @test trial.value === expected
        println((
            exponent_bits = ndigits(magnitude; base = 2),
            bytes = trial.bytes,
            seconds = trial.time,
        ))
        flush(stdout)
    end
end
println("STRICT_SCALED_SUM_CHECKS_COMPLETE")
flush(stdout)

end
