module DyadicChecks
using HybridPlasmaPIC, Test, Random
import ..Oracle
using ..Oracle: Q, oracle_log, oracle_integral
const H = HybridPlasmaPIC
const O = Oracle
exact(a) = (Q(a[1]) * Q(2)^a[3], Q(a[2]) * Q(2)^a[3])
@testset "Bounded dyadic interval operations against exact rationals" begin
    rng = MersenneTwister(0x6138)
    endpoints = Int128[typemin(Int128), -123, -1, 0, 1, 123, typemax(Int128)]
    for x in endpoints, shift in (0, 1, 63, 126, 127, 128, 1024)
        lower, upper = H._isothermal_outward_shift(x, shift)
        value = Q(x) / Q(2)^shift
        @test lower == floor(BigInt, value)
        @test upper == ceil(BigInt, value)
    end
    for repetition = 1:256
        lo, hi = minmax(rand(rng, Int128) >> 3, rand(rng, Int128) >> 3)
        a = (lo, hi, rand(rng, -2200:2200))
        lo, hi = minmax(rand(rng, Int128) >> 3, rand(rng, Int128) >> 3)
        b = (lo, hi, rand(rng, -2200:2200))
        alo, ahi = exact(a)
        blo, bhi = exact(b)
        lower, upper = exact(H._isothermal_add_dyadic(a, b))
        @test lower <= alo + blo <= ahi + bhi <= upper
        for T in (Float32, Float64)
            factor = reinterpret(T, rand(rng, T === Float32 ? UInt32 : UInt64))
            isfinite(factor) || continue
            lower, upper = exact(H._isothermal_multiply_dyadic(a, factor))
            targetlo, targethi = minmax(alo * Q(factor), ahi * Q(factor))
            @test lower <= targetlo <= targethi <= upper
        end
    end
    for T in (Float32, Float64), q in (0, 1, precision(T) - 8), active in (1, 256, 4194304)
        cap = min(124 - q, 124 - H._isothermal_ceil_log2(active))
        for bits in (precision(T) + 16, cap)
            K = T === Float32 ? 8 : 12
            r = precision(T) - q
            exact_tail = Q(active) * Q(2)^(bits - r * (K - 2)) / (K * (K + 1) * (BigInt(2)^r - 1))
            @test H._isothermal_tail_units(active, r, bits, Val(K)) == ceil(BigInt, exact_tail)
        end
    end
end

@testset "Direct dyadic rounding including all signs and range boundaries" begin
    rng = MersenneTwister(0x5817)
    for T in (Float32, Float64)
        exponents = (
            -2200,
            -1076,
            -1075,
            -1074,
            -1022,
            -151,
            -150,
            -149,
            -126,
            -53,
            -24,
            -1,
            0,
            1,
            100,
            127,
            128,
            1023,
            1024,
            2200,
        )
        for v in Int128[typemin(Int128), -1, 0, 1, typemax(Int128)], e in exponents
            @test isequal(H._from_exact(T, v, e), HybridPlasmaPIC._from_exact(T, Q(v) * Q(2)^e))
        end
        # Half-way points around successive significands, smallest subnormal and overflow.
        for q in (
            Base.exponent(floatmin(T)) - precision(T) + 1,
            -precision(T),
            -1,
            Base.exponent(floatmax(T)) - precision(T) + 1,
        )
            for k in (
                (Int128(1) << precision(T)) - 2,
                (Int128(1) << precision(T)) - 1,
                Int128(0),
                Int128(1),
            )
                for offset in (-1, 0, 1), sign in (-1, 1)
                    v = sign * (4k + 2 + offset)
                    @test isequal(
                        H._from_exact(T, v, q - 2),
                        HybridPlasmaPIC._from_exact(T, Q(v) * Q(2)^(q - 2)),
                    )
                end
            end
        end
        for repetition = 1:1024
            v, e = rand(rng, Int128), rand(rng, -2200:2200)
            @test isequal(H._from_exact(T, v, e), HybridPlasmaPIC._from_exact(T, Q(v) * Q(2)^e))
        end
        @test isequal(H._from_exact(T, Int128(-1), typemin(Int)), -zero(T))
        @test isequal(H._from_exact(T, Int128(1), typemax(Int)), T(Inf))
    end
end

@testset "Complete dyadic intervals enclose rational parent with full mantissas" begin
    rng = MersenneTwister(0x8781)
    for T in (Float32, Float64), C in (Float32, Float64), repetition = 1:48
        radius = Int128(1) << (precision(T) - 9)
        n = [one(T) + ldexp(T(rand(rng, -radius:radius)), -precision(T)) for i = 1:rand(rng, 1:32)]
        temperature = ldexp(C(rand(rng) - 0.5), rand(rng, -120:120))
        dx = ntuple(i -> ldexp(T(rand(rng) - 0.5), rand(rng, -120:120)), 3)
        state = H._isothermal_statistics(n)
        bits = precision(T) + 16
        degree = Val(T === Float32 ? 8 : 12)
        lower, upper = exact(H._isothermal_complete_dyadic(n, temperature, dx, state, bits, degree))
        reflo, refhi = O.complete_bounds(n, temperature, dx, state, bits, degree)
        @test lower <= reflo <= refhi <= upper
        expected_lo = HybridPlasmaPIC._from_exact(T, reflo)
        expected_hi = HybridPlasmaPIC._from_exact(T, refhi)
        value = O.evaluate(n, temperature, dx)
        @test value !== nothing
        if isequal(expected_lo, expected_hi)
            @test isequal(value, expected_lo)
        end
    end
end
println("ISOTHERMAL_DYADIC_COMPLETE")

@testset "Narrow polynomial path agrees bit for bit at its arithmetic bound" begin
    rng = MersenneTwister(0x7315)
    for T in (Float32, Float64), q in (0, 1, 8, 16), count in (1, 17, 256)
        radius = Int128(1) << q
        n = [one(T) + ldexp(T(rand(rng, -radius:radius)), -precision(T)) for i = 1:count]
        # Include the negative boundary exactly in either input precision.
        n[1] = one(T) - ldexp(one(T), q - precision(T))
        state = H._isothermal_statistics(n)
        degree = Val(T === Float32 ? 8 : 12)
        for bits in (40, 60 - state.q)
            @test bits + state.q <= 60
            narrow = H._isothermal_scaled_polynomial_storage(n, state, bits, degree, Int64)
            wide = H._isothermal_scaled_polynomial_storage(n, state, bits, degree, Int128)
            @test narrow == wide
            @test H._isothermal_scaled_polynomial(n, state, bits, degree) == wide
        end
        bits = 61 - state.q
        @test H._isothermal_scaled_polynomial(n, state, bits, degree) ==
              H._isothermal_scaled_polynomial_storage(n, state, bits, degree, Int128)
    end
end
println("ISOTHERMAL_NARROW_KERNEL_COMPLETE")

end
