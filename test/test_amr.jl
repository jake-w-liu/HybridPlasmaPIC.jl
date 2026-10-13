using HybridPlasmaPIC
using Random
using Test

# AMR oracle expressions are dyadic. Represent the complete exact value at a
# precision large enough for both integers, then round directly to the target.
# Avoid Base's ambient-precision Rational-to-hardware-float conversion.
function amr_oracle_round(::Type{T}, exact::Rational{BigInt}) where {T}
    n, d = numerator(exact), denominator(exact)
    @assert ispow2(d)
    wide = setprecision(BigFloat, max(2, ndigits(n; base = 2), ndigits(d; base = 2))) do
        BigFloat(n) / BigFloat(d)
    end
    return T === BigFloat ? BigFloat(wide; precision = precision(BigFloat)) : T(wide, RoundNearest)
end

@testset "AMR exact outputs ignore ambient BigFloat settings" begin
    for T in (Float16, Float32, Float64)
        tiny, unit = nextfloat(zero(T)), eps(one(T))
        for bits in (16, 53, 256, 2048), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
            setprecision(bits) do
                setrounding(BigFloat, mode) do
                    for sign in (-one(T), one(T)), direction in (-1, 0, 1)
                        expected = sign * (direction > 0 ? nextfloat(one(T)) : one(T))
                        coarse = AMRGrid(sign .* T[8direction*tiny, 1+unit, 4unit], one(T))
                        @test refine(coarse).u[3] == expected
                        values = sign .* T[1, -unit/4, direction*tiny]
                        @test HybridPlasmaPIC._amr_second_difference(values, 2) == expected
                    end
                    # This index takes the exact complete-coordinate fallback.
                    # The oracle retains its half-integer before final rounding.
                    grid = AMRGrid(T[0], tiny; x0 = tiny)
                    i = typemax(Int)
                    exact = Rational{BigInt}(tiny) * (2BigInt(i) + 1) / 2
                    @test cell_center(grid, i) == amr_oracle_round(T, exact)
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                end
            end
        end
    end
end

@testset "AMR restriction preserves stored-precision cancellation" begin
    rng = MersenneTwister(318)
    for stored_bits in (128, 256, 512)
        inputs, exact = setprecision(stored_bits) do
            delta = BigFloat(2)^(-(stored_bits - 8))
            pairs = [
                (BigFloat(1) + delta, BigFloat(-1)),
                (BigFloat(-1), BigFloat(1) + delta),
                (-BigFloat(1) - delta, BigFloat(1)),
                (BigFloat(1), -BigFloat(1) - delta),
                (BigFloat(0), BigFloat(1) + delta),
            ]
            append!(
                pairs,
                [ntuple(_ -> ldexp(BigFloat(randn(rng)), rand(rng, -100:100)), 2) for _ = 1:100],
            )
            values = BigFloat[value for pair in pairs for value in pair]
            averages = [(Rational{BigInt}(a) + Rational{BigInt}(b)) / 2 for (a, b) in pairs]
            values, averages
        end
        for bits in (24, 64, 128), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
            setprecision(bits) do
                setrounding(BigFloat, mode) do
                    fine = AMRGrid(inputs, BigFloat(0.5))
                    coarse = AMRGrid(zeros(BigFloat, length(exact)), BigFloat(1))
                    output = coarse.u
                    @test restrict!(coarse, fine) === coarse
                    @test coarse.u === output
                    @test coarse.u == [amr_oracle_round(BigFloat, value) for value in exact]
                    @test all(x -> precision(x) == bits, coarse.u)
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                    @test !iszero(coarse.u[1])
                end
            end
        end
    end

    for T in (Float16, Float32, Float64)
        M, tiny = floatmax(T), nextfloat(zero(T))
        pairs = [(M, M), (M, -M), (tiny, tiny), (tiny, zero(T)), (-tiny, zero(T))]
        append!(pairs, [Tuple(rand(rng, T, 2)) for _ = 1:100])
        values = T[value for pair in pairs for value in pair]
        expected = [
            amr_oracle_round(T, (Rational{BigInt}(a) + Rational{BigInt}(b)) / 2) for (a, b) in pairs
        ]
        fine = AMRGrid(values, T(0.5))
        coarse = AMRGrid(zeros(T, length(pairs)), one(T))
        @test restrict!(coarse, fine) === coarse
        @test coarse.u == expected
        @test (@allocated restrict!(coarse, fine)) == 0
        for (a, b) in ((T(Inf), T(Inf)), (T(Inf), T(-Inf)), (T(NaN), one(T)))
            fine = AMRGrid(T[a, b], T(0.5))
            coarse = AMRGrid(zeros(T, 1), one(T))
            restrict!(coarse, fine)
            @test isequal(coarse.u[1], (a + b) / T(2))
        end
    end
end

@testset "AMR refinement curvature" begin
    exact_curvature(a, b, c) = amr_oracle_round(
        typeof(a),
        Rational{BigInt}(a) - 2Rational{BigInt}(b) + Rational{BigInt}(c),
    )

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
            expected = amr_oracle_round(BigFloat, precise_curvature)
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
                              amr_oracle_round(BigFloat, exact)
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

@testset "AMR prolongation rounds complete child values" begin
    function exact_children(values)
        T = eltype(values)
        exact = Rational{BigInt}.(values)
        n = length(values)
        result = T[]
        for i = 1:n
            adjustment =
                n == 1 ? zero(exact[i]) :
                i == 1 ? (exact[2] - exact[1]) / 4 :
                i == n ? (exact[n] - exact[n-1]) / 4 : (exact[i+1] - exact[i-1]) / 8
            push!(
                result,
                amr_oracle_round(T, exact[i] - adjustment),
                amr_oracle_round(T, exact[i] + adjustment),
            )
        end
        return result
    end

    rng = MersenneTwister(621)
    for (T, U) in ((Float16, UInt16), (Float32, UInt32), (Float64, UInt64))
        M, tiny = floatmax(T), nextfloat(zero(T))
        # Rounding the interior slope before subtracting it from M/8 loses
        # the representable minimum-subnormal child.
        coarse = AMRGrid(T[8tiny, M/8, M], one(T))
        @test refine(coarse).u[3] == tiny
        # A rounded normalized value can rescale to the smallest normal even
        # when the exact child rounds to its subnormal neighbor. Check both
        # signs, both child orientations, and nearby subnormal lattice values.
        for sign in (-one(T), one(T)), neighbor in (-8:8), small in (0:16)
            values = sign .* T[0, floatmin(T)+T(neighbor)*tiny, T(small)*tiny]
            @test refine(AMRGrid(values, one(T))).u == exact_children(values)
            @test refine(AMRGrid(reverse(values), one(T))).u == exact_children(reverse(values))
        end
        for values in (
            T[8tiny, M/8, M],
            T[M, M/8, 8tiny],
            T[-8tiny, -M/8, -M],
            T[M, M],
            T[-M, M],
            T[tiny, tiny],
            T[tiny, -tiny, tiny],
            T[1, nextfloat(one(T)), 1],
            T[prevfloat(one(T)), 1, nextfloat(one(T))],
        )
            @test refine(AMRGrid(values, one(T))).u == exact_children(values)
        end

        for n in (2, 3, 5), _ = 1:150
            values = map(1:n) do _
                while true
                    value = reinterpret(T, rand(rng, U))
                    isfinite(value) && return value
                end
            end
            @test refine(AMRGrid(values, one(T))).u == exact_children(values)
        end

        # Affine dyadic profiles reproduce exact fine-cell averages, including
        # both one-sided boundaries, and conserve their finite-volume integral.
        for slope in (-2, 0, 2)
            coarse = AMRGrid(T[3 + slope * (i - 0.5) for i = 1:8], one(T))
            fine = refine(coarse)
            @test fine.u == T[3 + slope * (i / 2 - 0.25) for i = 1:16]
            @test sum(Float64, fine.u) / 2 == sum(Float64, coarse.u)
        end
        for value in (zero(T), -zero(T), tiny, M, T(Inf), T(NaN))
            expected = T[value-zero(T), value+zero(T)]
            @test isequal(refine(AMRGrid(T[value], one(T))).u, expected)
        end
        zeros_grid = AMRGrid(fill(-zero(T), 3), one(T))
        @test isequal(refine(zeros_grid).u, repeat(T[-zero(T), zero(T)], 3))

        for exceptional in (T(Inf), T(-Inf), T(NaN))
            values = T[1, exceptional, 3]
            expected = T[]
            for i = 1:3
                adjustment =
                    i == 1 ? (values[2] - values[1]) / T(4) :
                    i == 3 ? (values[3] - values[2]) / T(4) : (values[3] - values[1]) / T(8)
                push!(expected, values[i] - adjustment, values[i] + adjustment)
            end
            @test isequal(refine(AMRGrid(values, one(T))).u, expected)
        end

        # Include rounded ordinary slopes and ties: exact reconstruction must
        # not require heap allocations in these common finite-volume profiles.
        for values in (
            T.(sin.(range(0, 2pi; length = 512))),
            rand(rng, T, 512),
            repeat(T[1, nextfloat(one(T)), 1], 100),
        )
            coarse = AMRGrid(values, one(T))
            fine = refine(coarse)
            prolong!(fine, coarse)
            @test (@allocated prolong!(fine, coarse)) == 0
            @test fine.u == exact_children(values)
        end
    end

    # Preserve stored input bits if a caller lowers its BigFloat precision.
    # Build independent exact values before entering the lower-precision scope.
    stored = setprecision(256) do
        delta = BigFloat(2)^(-100)
        inputs = (BigFloat[8delta, 1, 8], BigFloat[1+delta, 5], BigFloat[5, 1+delta])
        [(values, Rational{BigInt}.(values)) for values in inputs]
    end
    for bits in (24, 64), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
        setprecision(bits) do
            setrounding(BigFloat, mode) do
                for (values, exact) in stored
                    n = length(values)
                    expected = BigFloat[]
                    for i = 1:n
                        adjustment =
                            i == 1 ? (exact[2] - exact[1]) / 4 :
                            i == n ? (exact[n] - exact[n-1]) / 4 : (exact[i+1] - exact[i-1]) / 8
                        push!(
                            expected,
                            amr_oracle_round(BigFloat, exact[i] - adjustment),
                            amr_oracle_round(BigFloat, exact[i] + adjustment),
                        )
                    end
                    fine = refine(AMRGrid(values, BigFloat(1)))
                    @test fine.u == expected
                    @test all(x -> precision(x) == bits, fine.u)
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                end
                for _ = 1:100
                    values = BigFloat[ldexp(BigFloat(randn(rng)), rand(rng, -100:100)) for _ = 1:3]
                    @test refine(AMRGrid(values, BigFloat(1))).u == exact_children(values)
                end
            end
        end
    end
end

@testset "AMR cell centers retain the complete coordinate" begin
    exact_center(g, i) = amr_oracle_round(
        promote_type(eltype(g.u), typeof(i)),
        Rational{BigInt}(g.x0) + (2BigInt(i) - 1) * Rational{BigInt}(g.dx) / 2,
    )

    rng = MersenneTwister(725)
    for (T, U) in ((Float16, UInt16), (Float32, UInt32), (Float64, UInt64))
        M = floatmax(T)
        extreme = AMRGrid(T[0, 0], M; x0 = -M)
        @test cell_center(extreme, 2) == M / 2
        for i in (0, 1, 2, 1024, 1025, 2049, 10000)
            grid = AMRGrid(T[0], one(T))
            @test cell_center(grid, i) == exact_center(grid, i)
            cell_center(grid, i)
            @test (@allocated cell_center(grid, i)) == 0
        end

        for _ = 1:500
            x0, dx = ntuple(2) do _
                while true
                    value = reinterpret(T, rand(rng, U))
                    isfinite(value) && !iszero(value) && return value
                end
            end
            grid = AMRGrid(T[0], abs(dx); x0 = x0)
            i = rand(rng, -50000:50000)
            @test cell_center(grid, i) == exact_center(grid, i)
        end

        boundary = Int64(1) << (precision(T) - 1)
        for i in (1 - boundary, -boundary, boundary, boundary + 1, typemin(Int), typemax(Int))
            grid = AMRGrid(T[0], T(0.25); x0 = -one(T))
            @test cell_center(grid, i) == exact_center(grid, i)
        end
        for I in (Int8, UInt8, Int128, UInt128, BigInt)
            grid = AMRGrid(T[0, 0], one(T))
            @test cell_center(grid, I(2)) == exact_center(grid, I(2))
            @test typeof(cell_center(grid, I(2))) === promote_type(T, I)
        end
    end

    @test cell_center(AMRGrid(zeros(Float16, 2050), Float16(1); x0 = Float16(-2048)), 2049) ==
          Float16(0.5)
    half_grid = AMRGrid(Float16[0, 0], nextfloat(Float16(1)); x0 = -nextfloat(Float16(0)))
    @test cell_center(half_grid, 2) == exact_center(half_grid, 2)
    for i in (-1, 2), x0 in (-nextfloat(Float16(0)), Float16(0), nextfloat(Float16(0)))
        grid = AMRGrid(Float16[0, 0], Float16(43680); x0 = x0)
        # 1.5 * 43680 is the Float16 overflow midpoint, 65520.
        @test cell_center(grid, i) == exact_center(grid, i)
    end

    # Exercise large indices, including those near the largest half-integer
    # factor that the widened Float64 computation can represent exactly.
    wide_boundary = Int64(1) << (precision(Float64) - 1)
    for _ = 1:500
        grid = AMRGrid(Float16[0], Float16(0.001); x0 = rand(rng, Float16))
        i = rand(rng, (1-wide_boundary):wide_boundary)
        @test cell_center(grid, i) == exact_center(grid, i)
    end

    values = setprecision(256) do
        dx = BigFloat(1) + BigFloat(2)^(-100)
        x0 = -BigFloat(1.5)
        grid = AMRGrid(BigFloat[0, 0], dx; x0 = x0)
        grid, Rational{BigInt}(x0) + (3 // 2) * Rational{BigInt}(dx)
    end
    setprecision(24) do
        grid, exact = values
        @test cell_center(grid, 2) == amr_oracle_round(BigFloat, exact)
        @test precision(BigFloat) == 24
    end
end
