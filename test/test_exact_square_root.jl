module ExactSquareRootTests
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

dyadic(e::Integer) = e >= 0 ? (BigInt(1) << e) // BigInt(1) : BigInt(1) // (BigInt(1) << -e)

# Read the stored value without floating-point scaling or dependence on the
# caller's precision. MPFR's exponent range also applies to ldexp and division.
function stored_rational(value::BigFloat)
    significand, power, sign = Base.decompose(value)
    return (BigInt(sign) * significand) * dyadic(power)
end

function check_hardware_cell(::Type{T}, exact) where {T}
    actual = H._sqrt_from_exact(T, exact)
    @test actual isa T
    @test !signbit(actual)
    @test !isnan(actual)
    p = precision(T)
    tiny = dyadic(exponent(floatmin(T)) - p + 1)
    maximum = Rational{BigInt}(floatmax(T))
    overflow = maximum + dyadic(exponent(floatmax(T)) - p)
    if isinf(actual)
        @test exact >= overflow^2
    elseif iszero(actual)
        @test exact <= (tiny / 2)^2
    else
        value = Rational{BigInt}(actual)
        previous = Rational{BigInt}(prevfloat(actual))
        following = nextfloat(actual)
        upper_neighbor =
            isfinite(following) ? Rational{BigInt}(following) : dyadic(exponent(floatmax(T)) + 1)
        lower, upper = (previous + value) / 2, (value + upper_neighbor) / 2
        U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
        even = iseven(reinterpret(U, actual))
        # Squaring is monotone on these nonnegative rational cell boundaries.
        # This oracle neither takes an integer square root nor rounds an input.
        @test lower^2 < exact < upper^2 || (even && (exact == lower^2 || exact == upper^2))
    end
    return actual
end

function hardware_cases(::Type{T}, rng) where {T}
    p = precision(T)
    emin, emax = exponent(floatmin(T)), exponent(floatmax(T))
    tiny = dyadic(emin - p + 1)
    cases = Rational{BigInt}[0//1, 1//3, 2//1, 4//1]
    # Probe both sides of binade transitions, including the subnormal boundary.
    for power in (emin - p + 1, emin - 1, emin, -1, 0, 1, emax - 1, emax)
        value = ldexp(one(T), power)
        exact = Rational{BigInt}(value)
        previous, following = Rational{BigInt}(prevfloat(value)), Rational{BigInt}(nextfloat(value))
        push!(cases, exact^2)
        for (a, b) in ((previous, exact), (exact, following)),
            tail in (4, 1200),
            direction in (-1, 0, 1)

            midpoint, gap = (a + b) / 2, b - a
            push!(cases, (midpoint + direction * gap * dyadic(-tail))^2)
            # Nondyadic radicands near the same squared rounding boundary.
            push!(cases, midpoint^2 + direction * gap^2 * dyadic(-tail) / 3)
        end
    end
    maximum = Rational{BigInt}(floatmax(T))
    overflow = maximum + dyadic(emax - p)
    for boundary in (tiny / 2, dyadic(emin) - tiny / 2, overflow),
        tail in (4, 1200),
        direction in (-1, 0, 1)

        push!(cases, (boundary + direction * tiny * dyadic(-tail))^2)
    end
    append!(cases, (maximum^2, dyadic(2 * (emax + 100)), dyadic(2 * (emin - p - 100))))
    for _ = 1:100
        numerator = BigInt(rand(rng, UInt128)) + 1
        denominator = BigInt(rand(rng, UInt128)) + 1
        push!(cases, (numerator // denominator) * dyadic(rand(rng, (2*(emin-p-20)):(2*(emax+20)))))
    end
    return cases
end

function check_bigfloat_cell(exact)
    bits = precision(BigFloat)
    actual = H._sqrt_from_exact(BigFloat, exact)
    @test actual isa BigFloat
    @test precision(actual) == bits
    @test isfinite(actual)
    @test !signbit(actual)
    if iszero(exact)
        @test iszero(actual)
    else
        @test actual > 0
        value = stored_rational(actual)
        previous, following = stored_rational(prevfloat(actual)), stored_rational(nextfloat(actual))
        lower, upper = (previous + value) / 2, (value + following) / 2
        quantum = dyadic(exponent(actual) - bits + 1)
        significand = value / quantum
        @test denominator(significand) == 1
        even = iseven(numerator(significand))
        @test lower^2 < exact < upper^2 || (even && (exact == lower^2 || exact == upper^2))
    end
    return actual
end

@testset "Exact square roots round to nearest, ties to even" begin
    @testset "Every positive finite Float16 square" begin
        # Positive finite IEEE binary16 encodings run from 0x0001 to 0x7bff.
        for encoding = UInt16(0x0001):UInt16(0x7bff)
            value = reinterpret(Float16, encoding)
            @test H._sqrt_from_exact(Float16, Rational{BigInt}(value)^2) === value
        end
    end

    rng = MersenneTwister(0x5a17)
    @testset "IEEE boundaries and independent rational cells" begin
        for T in (Float16, Float32, Float64)
            cases = hardware_cases(T, rng)
            for bits in (24, 96, 256), mode in (RoundNearest, RoundUp, RoundDown)
                setprecision(BigFloat, bits) do
                    setrounding(BigFloat, mode) do
                        for exact in cases
                            check_hardware_cell(T, exact)
                        end
                        @test_throws DomainError H._sqrt_from_exact(T, BigInt(-1) // BigInt(3))
                        @test precision(BigFloat) == bits
                        @test rounding(BigFloat) == mode
                    end
                end
            end
        end
    end

    stored = setprecision(BigFloat, 512) do
        BigFloat[1+ldexp(BigFloat(1), -400), ldexp(BigFloat(1), -600), ldexp(BigFloat(1), 400)+1]
    end
    exact_stored = stored_rational.(stored)
    @testset "BigFloat working precision and caller rounding" begin
        for bits in (96, 128, 256)
            cases = Rational{BigInt}[0//1, 1//3, 2//1, 4//1]
            append!(cases, exact_stored)
            for power in (-1200, -1, 0, 1, 1200), offset in (0, 1)
                # Consecutive significands exercise both even and odd ties.
                quantum = dyadic(power - bits + 1)
                lower = dyadic(power) + offset * quantum
                midpoint = lower + quantum / 2
                for tail in (4, 1200), direction in (-1, 0, 1)
                    push!(cases, (midpoint + direction * quantum * dyadic(-tail))^2)
                    push!(cases, midpoint^2 + direction * quantum^2 * dyadic(-tail) / 3)
                end
            end
            for _ = 1:100
                push!(
                    cases,
                    ((BigInt(rand(rng, UInt128)) + 1) // (BigInt(rand(rng, UInt128)) + 1)) *
                    dyadic(rand(rng, -2400:2400)),
                )
            end
            setprecision(BigFloat, bits) do
                for mode in (RoundNearest, RoundUp, RoundDown)
                    setrounding(BigFloat, mode) do
                        @test H._exact_rational.(stored) == exact_stored
                        for exact in cases
                            check_bigfloat_cell(exact)
                        end
                        @test_throws DomainError H._sqrt_from_exact(
                            BigFloat,
                            BigInt(-1) // BigInt(3),
                        )
                        @test precision(BigFloat) == bits
                        @test rounding(BigFloat) == mode
                    end
                end
            end
        end
    end

    @testset "Concurrent calls preserve numerical context" begin
        for T in (Float16, Float32, Float64, BigFloat)
            cases =
                T === BigFloat ?
                Rational{BigInt}[(BigInt(i) // BigInt(i + 1)) * dyadic(i - 128) for i = 1:256] :
                hardware_cases(T, rng)
            expected =
                [T === BigFloat ? check_bigfloat_cell(x) : check_hardware_cell(T, x) for x in cases]
            actual = Vector{T}(undef, length(cases))
            preserved = Vector{Bool}(undef, length(cases))
            Threads.@threads for i in eachindex(cases)
                bits, mode = precision(BigFloat), rounding(BigFloat)
                actual[i] = H._sqrt_from_exact(T, cases[i])
                preserved[i] = precision(BigFloat) == bits && rounding(BigFloat) == mode
            end
            @test isequal(actual, expected)
            @test all(preserved)
        end
    end

    @testset "Stored BigFloat conversion preserves bits and context" begin
        cases = Tuple{BigFloat,Rational{BigInt}}[]
        for bits in (24, 96, 128, 512), power in (-600, 0, 600), sign in (-1, 1)
            value = setprecision(BigFloat, bits) do
                sign * ldexp(BigFloat(1) + ldexp(BigFloat(1), 2 - bits), power)
            end
            expected = sign * (1 + dyadic(2 - bits)) * dyadic(power)
            push!(cases, (value, expected))
        end
        for bits in (24, 128), mode in (RoundNearest, RoundUp, RoundDown)
            setprecision(BigFloat, bits) do
                setrounding(BigFloat, mode) do
                    actual = Vector{Rational{BigInt}}(undef, length(cases))
                    preserved = Vector{Bool}(undef, length(cases))
                    Threads.@threads for i in eachindex(cases)
                        before = precision(BigFloat), rounding(BigFloat)
                        actual[i] = H._exact_rational(cases[i][1])
                        preserved[i] = (precision(BigFloat), rounding(BigFloat)) == before
                    end
                    @test actual == last.(cases)
                    @test all(preserved)
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                    @test H._exact_rational(BigFloat(Inf)) == BigInt(1) // BigInt(0)
                    @test H._exact_rational(BigFloat(-Inf)) == BigInt(-1) // BigInt(0)
                    @test_throws ArgumentError H._exact_rational(BigFloat(NaN))
                end
            end
        end
    end

    @testset "MPFR exponent endpoints in an isolated process" begin
        project = Base.active_project()
        @test project !== nothing
        fixture = joinpath(@__DIR__, "ionization_bigfloat_endpoint_fixture.jl")
        package = dirname(@__DIR__)
        command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(project)) $fixture $package`
        log = IOBuffer()
        process = run(pipeline(ignorestatus(command); stdout = log, stderr = log))
        transcript = String(take!(log))
        success(process) || @error "BigFloat endpoint fixture failed" transcript
        @test success(process)
        @test occursin("BigFloat endpoint settings restored", transcript)
    end
end
end
