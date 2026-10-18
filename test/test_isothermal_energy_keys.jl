module IsothermalKeyChecks
using HybridPlasmaPIC, Test, Random
include("isothermal_energy/oracles.jl")
using .Oracle: Q, oracle_integral
@testset "Finite binary density logarithm keys fit UInt64 exactly" begin
    rng = MersenneTwister(0x6917)
    checked = 0
    for T in (Float32, Float64)
        U = T === Float32 ? UInt32 : UInt64
        values = T[
            nextfloat(zero(T)),
            prevfloat(floatmin(T)),
            floatmin(T),
            nextfloat(floatmin(T)),
            prevfloat(one(T)),
            one(T),
            nextfloat(one(T)),
            floatmax(T),
        ]
        append!(values, [reinterpret(T, rand(rng, U) & (typemax(U) >> 1)) for i = 1:512])
        filter!(x -> isfinite(x) && x > zero(T), values)
        for value in values
            density = Q(value)
            odd = numerator(density)
            exponent = -(ndigits(denominator(density); base = 2) - 1)
            while iseven(odd)
                odd >>= 1
                exponent += 1
            end
            key = UInt64(odd)
            fraction, e = frexp(value)
            significand = UInt64(ldexp(fraction, precision(T)))
            shift = trailing_zeros(significand)
            @test key == significand >> shift
            @test exponent == e - precision(T) + shift
            @test Q(key) * Q(2)^exponent == density
            @test ndigits(odd; base = 2) <= precision(T)
            if key > 1
                terms = Dict{UInt64,Q}()
                HybridPlasmaPIC._isothermal_add_log!(terms, odd, density)
                @test terms[key] == density
                HybridPlasmaPIC._isothermal_add_log!(terms, odd, -density)
                @test isempty(terms)
            end
            checked += 1
        end
        # Include the full exponent extremes and fresh full-mantissa random values.
        for value in [values[1:8]; values[end-7:end]]
            spacing = ldexp(
                one(T),
                clamp(
                    -Base.exponent(value),
                    -Base.exponent(floatmax(T)),
                    Base.exponent(floatmax(T)),
                ),
            )
            n = T[value]
            bounds = oracle_integral(n, one(T), (spacing,))
            @test Oracle.rounded(
                HybridPlasmaPIC._isothermal_exact_energy(n, one(T), (spacing,)),
                bounds...,
            )
        end
    end
    @test checked > 1000
    println((checked_positive_finite_densities = checked,))
end
println("ISOTHERMAL_KEY_BOUND_COMPLETE")

end
