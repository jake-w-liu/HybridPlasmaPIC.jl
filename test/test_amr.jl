using HybridPlasmaPIC
using Random
using Test

@testset "AMR refinement curvature" begin
    exact_curvature(a, b, c) =
        typeof(a)(Rational{BigInt}(a) - 2Rational{BigInt}(b) + Rational{BigInt}(c))

    function check_curvature(a::T, b::T, c::T) where {T}
        u = T[a, b, c]
        expected = exact_curvature(a, b, c)
        @test HybridPlasmaPIC._amr_second_difference(u, 2) == expected
        @test refine_flags(u, zero(T)) == Bool[false, !iszero(expected), false]
        if isfinite(expected)
            @test !any(refine_flags(u, abs(expected)))
            if !iszero(expected)
                @test refine_flags(u, prevfloat(abs(expected))) == Bool[false, true, false]
            end
        end
    end

    rng = MersenneTwister(912)
    for (T, U) in ((Float16, UInt16), (Float32, UInt32), (Float64, UInt64))
        @testset "$T" begin
            for start in (-100, -1, 0, 1, 100), slope in (-20, -1, 0, 1, 20)
                a, b, c = T(start), T(start + slope), T(start + 2slope)
                check_curvature(a, b, c)
                @test !any(refine_flags(T[a, b, c], zero(T)))
            end

            M, tiny = floatmax(T), nextfloat(zero(T))
            for values in (
                (M, M, M),
                (M, -tiny, -M),
                (tiny, zero(T), tiny),
                (tiny, tiny, zero(T)),
                (M, -M, M),
                (one(T), nextfloat(one(T)), one(T)),
                (prevfloat(one(T)), one(T), nextfloat(one(T))),
            )
                check_curvature(values...)
            end

            for _ = 1:500
                values = ntuple(3) do _
                    while true
                        value = reinterpret(T, rand(rng, U))
                        isfinite(value) && return value
                    end
                end
                check_curvature(values...)
            end

            # These values exercise ordinary rounding without the allocating
            # exact fallback; public refine_flags allocates its output BitVector.
            for u in (T[1, 2, 3], T[0.25, -0.5, 0.125])
                HybridPlasmaPIC._amr_second_difference(u, 2)
                @test (@allocated HybridPlasmaPIC._amr_second_difference(u, 2)) == 0
            end
            @test !any(refine_flags(T[0, T(NaN), 0], zero(T)))
            @test refine_flags(T[0, T(Inf), 0], zero(T)) == Bool[false, true, false]
        end
    end

    @testset "generic precision and unchanged criterion boundaries" begin
        # Store a difference below the later working precision. Build the exact
        # oracle before lowering precision: Julia's Rational(BigFloat) conversion
        # otherwise rounds its intermediate scaling at the current precision.
        precise_values, precise_curvature = setprecision(256) do
            values = BigFloat[1, 1, 1+BigFloat(2)^(-100)]
            exact =
                Rational{BigInt}(values[1]) - 2Rational{BigInt}(values[2]) +
                Rational{BigInt}(values[3])
            values, exact
        end
        setprecision(24) do
            expected = BigFloat(precise_curvature)
            @test HybridPlasmaPIC._amr_second_difference(precise_values, 2) == expected
            @test refine_flags(precise_values, 0) == Bool[false, true, false]
            @test !any(refine_flags(precise_values, expected))
            @test precision(BigFloat) == 24
        end

        setprecision(24) do
            one_big = BigFloat(1)
            residuals = (eps(one_big) / 4, eps(one_big) / 2, BigFloat(2)^(-40))
            values = [
                (a, b, c) for residual in residuals for (a, b, c) in (
                    (one_big, zero(BigFloat), residual),
                    (-one_big, zero(BigFloat), residual),
                    (one_big, zero(BigFloat), -residual),
                    (one_big, one_big / 2, residual),
                )
            ]
            exacts = [
                Rational{BigInt}(a) - 2Rational{BigInt}(b) + Rational{BigInt}(c) for
                (a, b, c) in values
            ]
            for mode in (RoundDown, RoundUp, RoundToZero)
                setrounding(BigFloat, mode) do
                    for (abc, exact) in zip(values, exacts)
                        @test HybridPlasmaPIC._amr_second_difference(collect(abc), 2) ==
                              BigFloat(exact)
                    end
                    @test rounding(BigFloat) == mode
                end
            end
        end

        for bits in (24, 128)
            setprecision(bits) do
                check_curvature(BigFloat(1), BigFloat(2), BigFloat(3))
                for _ = 1:100
                    values = ntuple(_ -> ldexp(BigFloat(randn(rng)), rand(rng, -100:100)), 3)
                    check_curvature(values...)
                end
            end
        end
        for n = 0:2
            @test refine_flags(zeros(n), 0) == falses(n)
        end
        for threshold in (-1, Inf, NaN)
            @test_throws ArgumentError refine_flags([1.0, 2.0, 3.0], threshold)
        end
        for T in (Int8, UInt8, Int64, UInt64, Int128, UInt128)
            u = T[typemax(T), typemin(T), typemax(T)]
            expected = BigInt(u[1]) - 2BigInt(u[2]) + BigInt(u[3])
            @test HybridPlasmaPIC._amr_second_difference(u, 2) == expected
            @test refine_flags(u, 0) == Bool[false, true, false]
        end
        @test refine_flags([1 // 3, 2 // 3, 1 // 1], 0) == falses(3)
    end
end
