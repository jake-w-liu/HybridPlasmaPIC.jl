module ReuseChecks
using HybridPlasmaPIC, Test, Random, Statistics
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
apply!(pe, n, c, policy::Val{F}) where {F} =
    F ? electron_pressure!(pe, n, c) : H._polytropic_pressure_array!(pe, n, c, policy)
function reference(::Type{S}, n::Array{T}, c, ::Val{F}) where {S,T,F}
    return map(x -> H._polytropic_pressure(S, F ? max(x, zero(T)) : x, c), n)
end
function captured(f)
    try
        return (:value, f())
    catch error
        # Compare the unchanged scalar owner's actual failure type and message.
        return (:error, typeof(error), sprint(showerror, error))
    end
end
println((VERSION,))
flush(stdout)
@testset "Local pressure reuse preserves complete scalar outcomes" begin
    for S in (Float32, Float64),
        T in (Float32, Float64),
        C in (Float32, Float64),
        F in (false, true)

        for values in (
                T[],
                T[0, -0.0, 0, -0.0],
                T[0.5, 0.5, 1, 1, 1.5, 1.5, 2, 2],
                T[Inf, Inf, NaN, NaN, -Inf, -Inf],
            ),
            gamma in C[1, 1.5, 2, 5/3]

            c = PolytropicElectrons(C(0.7), C(1.25), gamma)
            n = copy(values)
            pe = similar(n, S)
            before = copy(n)
            expected = captured(() -> reference(S, n, c, Val(F)))
            actual = captured(() -> apply!(pe, n, c, Val(F)))
            @test actual[1] === expected[1]
            if expected[1] === :value
                @test actual[2] === pe
                @test isequal(pe, expected[2])
            else
                @test actual[2] === expected[2]
                @test actual[3] == expected[3]
            end
            @test isequal(n, before)
        end
        c = PolytropicElectrons(C(0.7), C(1.25), C(2))
        n = T[-2, -2, -1, -1, 0, 0, 1, 1]
        @test isequal(apply!(zeros(S, length(n)), n, c, Val(F)), reference(S, n, c, Val(F)))
        scalar = fill(T(1.5))
        out = fill(zero(S))
        @test apply!(out, scalar, c, Val(F)) === out
        @test isequal(out, reference(S, scalar, c, Val(F)))
        for coefficient in C[0.2, 0.7, 1.5]
            c = PolytropicElectrons(coefficient, C(1.25), C(2))
            @test isequal(apply!(zeros(S, length(n)), n, c, Val(F)), reference(S, n, c, Val(F)))
        end
    end
    for T in (Float32, Float64), F in (false, true)
        c = PolytropicElectrons(T(0.7), one(T), T(2))
        n = T[-1, -1, 0, 0, 0.5, 0.5, 1.5, 1.5]
        expected = reference(T, n, c, Val(F))
        @test apply!(n, n, c, Val(F)) === n
        @test isequal(n, expected)
        k = T === Float32 ? 100 : 600
        for sign in (-1, 1)
            density = ldexp(one(T), sign * k)
            n = fill(density, 256)
            c = PolytropicElectrons(inv(density), one(T), T(2))
            @test all(==(density), apply!(similar(n), n, c, Val(F)))
        end
        invalid = fill(T(-9), 2, 2)
        @test_throws DimensionMismatch apply!(invalid, ones(T, 4), c, Val(F))
        @test all(==(T(-9)), invalid)
    end
end
@testset "Raw failures and signed bits remain observable" begin
    for T in (Float32, Float64)
        c = PolytropicElectrons(T(0.7), one(T), T(1.5))
        n = T[-1, -1]
        @test_throws DomainError apply!(similar(n), n, c, Val(false))
        @test all(iszero, apply!(similar(n), n, c, Val(true)))
        U = T === Float32 ? UInt32 : UInt64
        # Adjacent NaN payloads must each pass through the existing scalar owner.
        bits = reinterpret(U, T(NaN))
        values = reinterpret.(T, U[bits, bits+1, bits+2])
        actual = apply!(similar(values), values, c, Val(false))
        expected = reference(T, values, c, Val(false))
        @test reinterpret.(U, actual) == reinterpret.(U, expected)
        zerosigned = T[-0.0, -0.0, 0.0, 0.0]
        c = PolytropicElectrons(-zero(T), one(T), one(T))
        actual = apply!(similar(zerosigned), zerosigned, c, Val(false))
        expected = reference(T, zerosigned, c, Val(false))
        @test reinterpret.(U, actual) == reinterpret.(U, expected)
    end
end
println("POLYTROPIC_ARRAY_REUSE_FOCUSED_COMPLETE")

end
