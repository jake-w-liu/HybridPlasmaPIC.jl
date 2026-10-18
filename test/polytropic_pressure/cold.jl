module ColdChecks
using HybridPlasmaPIC, Test, Random, Statistics
import ..Oracle
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function public_pressure(::Type{S}, density::T, closure) where {S,T}
    source = T[density]
    destination = Vector{S}(undef, 1)
    electron_pressure!(destination, source, closure)
    return only(destination)
end
public_candidate(density::T, closure) where {T} = public_pressure(T, density, closure)
@testset "Single rational power equality has bounded intermediate integers" begin
    for base in BigInt[1, 2, 3, 7, 17], power in BigInt[1, 2, 17, 65]
        exact = base^Int(power)
        @test H._polytropic_integer_power_equal(base, power, exact)
        @test !H._polytropic_integer_power_equal(base, power, exact + 1)
        exact > 1 && @test !H._polytropic_integer_power_equal(base, power, exact - 1)
    end
    for power in (BigInt(10)^300, BigInt(1) << 1023)
        @test !H._polytropic_integer_power_equal(BigInt(2), power, BigInt(1) << 4096)
        @test H._polytropic_integer_power_equal(BigInt(1), power, BigInt(1))
    end
    for q in (1, 2, 4, 8), p in (1, 3, 17)
        base = Q(7, 5)^q
        target = Q(7, 5)^p
        @test H._polytropic_power_equal(base, Q(p, q), target)
        @test !H._polytropic_power_equal(base, Q(p, q), target + 1)
        @test !H._polytropic_power_equal(Q(2), Q(1, 2), target)
    end
end

@testset "Complete weighted powers recover exact rational values" begin
    for T in (Float32, Float64), a in (1, 2, 7), b in (1, 3, 5), q in (1, 2, 4), p in (1, 3, 17)
        coefficient = Q(T(0.7))
        base = Q(a, b)^q
        exact = coefficient * Q(a, b)^p
        value = H._polytropic_complete_power(T, base, Q(p, q), coefficient)
        @test Oracle.rounded(value, exact, exact)
    end
end

@testset "Exact midpoint and both sides including zero and overflow" begin
    for T in (Float32, Float64)
        U = T === Float32 ? UInt32 : UInt64
        for left in (zero(T), nextfloat(zero(T)), one(T), nextfloat(one(T)), floatmax(T)),
            q in (1, 2, 4, 8)

            right = nextfloat(left)
            boundary = H._polytropic_midpoint_after(left)
            base = Q(7, 5)^q
            gamma = Q(17, q)
            coefficient = boundary / Q(7, 5)^17
            expected = iseven(reinterpret(U, left)) ? left : right
            @test H._polytropic_power_equal(base, gamma, boundary / coefficient)
            @test isequal(H._polytropic_complete_power(T, base, gamma, coefficient), expected)
            delta = BigInt(1) // (BigInt(1) << (2precision(T)))
            @test isequal(
                H._polytropic_complete_power(T, base, gamma, coefficient * (1 - delta)),
                left,
            )
            @test isequal(
                H._polytropic_complete_power(T, base, gamma, coefficient * (1 + delta)),
                right,
            )
        end
    end
end

@testset "Irrational roots use independent integer root enclosures" begin
    for T in (Float32, Float64), m in (2, 3, 5), q in (2, 4, 8), p in (1, 3)
        scale = BigInt(1) << 160
        target = BigInt(m) * scale^q
        root = target
        for _ = 1:trailing_zeros(q)
            root = isqrt(root)
        end
        @test root^q <= target < (root + 1)^q
        coefficient = Q(T(0.7))
        lo, hi = coefficient * (root // scale)^p, coefficient * ((root + 1) // scale)^p
        @test Oracle.rounded(H._polytropic_complete_power(T, Q(m), Q(p, q), coefficient), lo, hi)
    end
end

@testset "Complete exponent range and minute departures from a midpoint" begin
    for T in (Float32, Float64)
        huge = Q(floatmax(Float64))
        @test isinf(H._polytropic_complete_power(T, Q(2), huge, Q(1)))
        @test isequal(H._polytropic_complete_power(T, Q(1, 2), huge, Q(1)), zero(T))
        @test H._polytropic_complete_power(T, Q(1), huge, Q(T(0.7))) == T(0.7)
        boundary = H._polytropic_midpoint_after(one(T))
        tiny = Q(nextfloat(0.0))
        @test H._polytropic_complete_power(T, Q(2), tiny, boundary) == nextfloat(one(T))
        @test H._polytropic_complete_power(T, Q(1, 2), tiny, boundary) == one(T)
        @test H._polytropic_complete_power(T, Q(0), Q(1, 2), Q(1)) == 0
        @test H._polytropic_complete_power(T, Q(2), Q(1, 2), Q(0)) == 0
    end
    @test_throws ArgumentError H._polytropic_complete_power(Float64, Q(-1), Q(1), Q(1))
    @test_throws ArgumentError H._polytropic_complete_power(Float64, Q(1), Q(0), Q(1))
    @test_throws ArgumentError H._polytropic_complete_power(Float64, Q(1), Q(1, 3), Q(1))
    @test_throws ArgumentError H._polytropic_complete_power(Float64, Q(1), Q(1), Q(-1))
end
println("COMPLETE_POINTWISE_POWER_ROUNDING_CHECKS_COMPLETE")

end
