module ColdChecks
using HybridPlasmaPIC, Test, Random
import ..Oracle
using ..Oracle: Q, oracle_log, oracle_integral
const H = HybridPlasmaPIC
const O = Oracle
@testset "Independent rational logarithm enclosures" begin
    logtwo = H._isothermal_atanh_log_bounds(Q(1, 3), 80)
    for base in BigInt.([2:64; 2^40 - 1; 2^40 + 1])
        lo, hi = H._isothermal_integer_log_bounds(base, 80, logtwo)
        reflo, refhi = oracle_log(Q(base); bits = 160)
        @test lo <= reflo <= refhi <= hi
        @test hi - lo <= (ndigits(base; base = 2) + 1) // (BigInt(1) << 80)
    end
    @test_throws ArgumentError H._isothermal_atanh_log_bounds(Q(-1, 3), 80)
    @test_throws ArgumentError H._isothermal_atanh_log_bounds(Q(1, 2), 80)
    @test_throws ArgumentError H._isothermal_atanh_log_bounds(Q(0), 0)
    @test_throws ArgumentError H._isothermal_integer_log_bounds(BigInt(0), 80, logtwo)
end

@testset "Isothermal exact zero and logarithm-map contract" begin
    for T in (Float32, Float64), C in (Float32, Float64)
        for n in (T[], T[0, -1, 1], T[2, 0.5, 0.5, 0.5, 0.5], [T(4); fill(T(0.5), 16)])
            @test H._isothermal_exact_energy(n, C(3), (T(2),)) === zero(T)
            @test H._isothermal_exact_energy(n, -C(3), (T(2),)) === zero(T)
        end
        @test H._isothermal_exact_energy(T[3, 4], zero(C), (floatmax(T),)) === zero(T)
        @test_throws ArgumentError H._isothermal_exact_energy(T[Inf], C(1), (one(T),))
        @test_throws ArgumentError H._isothermal_exact_energy(T[1], C(NaN), (one(T),))
        @test_throws ArgumentError H._isothermal_exact_energy(T[1], C(1), (T(Inf),))
    end
    @test_throws ArgumentError H._isothermal_finish_log_energy(
        Float64,
        Dict(UInt64(3) => Q(1), UInt64(9) => Q(-1)),
    )
    @test_throws ArgumentError H._isothermal_finish_log_energy(Float64, Dict(UInt64(4) => Q(1)))
    @test_throws ArgumentError H._isothermal_finish_log_energy(Float64, Dict(UInt64(1) => Q(1)))
end

@testset "Complete isothermal pressure, volume and result ranges" begin
    for T in (Float32, Float64)
        cases = [
            (T[2], floatmax(T), (T(0.25),)),
            (T[0.5], nextfloat(zero(T)), (floatmax(T),)),
            (T[0.5], -nextfloat(zero(T)), (floatmax(T),)),
            (T[2], floatmax(T), (floatmax(T),)),
            (T[0.5], floatmax(T), (floatmax(T),)),
            (T[2], nextfloat(zero(T)), (T(0.5),)),
            (T[0.5], nextfloat(zero(T)), (one(T),)),
            (T[prevfloat(one(T)), nextfloat(one(T)), T(1.5)], one(T), (one(T),)),
            (T[0.75, 3, 9, 0, -1], T(0.125), (T(0.5), T(2), T(0.25))),
        ]
        for (n, temperature, dx) in cases
            before = copy(n)
            bounds = oracle_integral(n, temperature, dx)
            result = H._isothermal_exact_energy(n, temperature, dx)
            @test O.rounded(result, bounds...)
            @test isequal(n, before)
            @test isequal(H._isothermal_exact_energy(view(n, :), temperature, dx), result)
        end
        for C in (Float32, Float64)
            n = T[0.75, 1.5, 3]
            temperature = C(0.7)
            @test O.rounded(
                H._isothermal_exact_energy(n, temperature, (T(0.125),)),
                oracle_integral(n, temperature, (T(0.125),))...,
            )
        end
    end
end
println("ISOTHERMAL_EXACT_PROTOTYPE_COMPLETE")

end
