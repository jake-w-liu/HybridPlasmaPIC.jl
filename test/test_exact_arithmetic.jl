using HybridPlasmaPIC
using Random
using Test

@testset "Exact fallback conversions" begin
    dyadic(e) = e >= 0 ? (BigInt(1) << e) // BigInt(1) : BigInt(1) // (BigInt(1) << -e)

    # Certify the output's rounding cell with exact inequalities. This oracle
    # uses neither the production integer quotient nor Rational-to-float casts.
    function check_hardware(::Type{T}, exact) where {T}
        actual = HybridPlasmaPIC._from_exact(T, exact)
        U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
        overflow = Rational{BigInt}(floatmax(T)) + dyadic(exponent(floatmax(T)) - precision(T))
        if isinf(actual)
            @test abs(exact) >= overflow
        else
            value = Rational{BigInt}(actual)
            previous, following = prevfloat(actual), nextfloat(actual)
            lower_neighbor =
                isfinite(previous) ? Rational{BigInt}(previous) : -dyadic(exponent(floatmax(T)) + 1)
            upper_neighbor =
                isfinite(following) ? Rational{BigInt}(following) :
                dyadic(exponent(floatmax(T)) + 1)
            lower, upper = (lower_neighbor + value) / 2, (value + upper_neighbor) / 2
            even = iseven(reinterpret(U, actual))
            @test lower < exact < upper || (even && (exact == lower || exact == upper))
        end
        @test signbit(actual) == (exact < 0)
    end

    rng = MersenneTwister(821)
    for T in (Float16, Float32, Float64)
        emin, emax, p = exponent(floatmin(T)), exponent(floatmax(T)), precision(T)
        tiny = dyadic(emin - p + 1)
        cases = Rational{BigInt}[0//1]
        for e in (emin - p + 1, emin - 1, emin, -1, 0, 1, emax - 1, emax)
            value = ldexp(one(T), e)
            previous, following =
                Rational{BigInt}(prevfloat(value)), Rational{BigInt}(nextfloat(value))
            exact = Rational{BigInt}(value)
            for (a, b) in ((previous, exact), (exact, following))
                midpoint, gap = (a + b) / 2, b - a
                for bits in (2, 1200), direction in (-1, 0, 1), sign in (-1, 1)
                    delta = gap / (BigInt(1) << bits)
                    push!(cases, sign * (midpoint + direction * delta))
                    push!(cases, sign * (midpoint + direction * delta / 3))
                end
            end
        end
        overflow = Rational{BigInt}(floatmax(T)) + dyadic(emax - p)
        for midpoint in (tiny / 2, dyadic(emin) - tiny / 2, overflow),
            delta in (tiny / 8, tiny / (BigInt(1) << 1200)),
            direction in (-1, 0, 1),
            sign in (-1, 1)

            push!(cases, sign * (midpoint + direction * delta))
        end
        midpoint = 1 + Rational{BigInt}(eps(one(T))) / 2
        for sign in (-1, 1), direction in (-1, 0, 1)
            push!(cases, sign * (midpoint + direction * tiny))
        end
        for _ = 1:50
            n, d = BigInt(rand(rng, Int128)), BigInt(rand(rng, UInt128)) + 1
            push!(cases, (n // d) * dyadic(rand(rng, (emin-p-50):(emax+50))))
        end
        for bits in (16, 53, 256, 2048), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
            setprecision(bits) do
                setrounding(BigFloat, mode) do
                    for exact in cases
                        check_hardware(T, exact)
                    end
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                end
            end
        end
    end

    stored, expected = setprecision(512) do
        values = BigFloat[1+BigFloat(2)^(-400), -1-BigFloat(2)^(-400), BigFloat(2)^(-500)]
        values, Rational{BigInt}.(values)
    end
    for bits in (24, 128), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
        setprecision(bits) do
            setrounding(BigFloat, mode) do
                @test HybridPlasmaPIC._exact_rational.(stored) == expected
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
            end
        end
    end

    # Non-dyadic quotients must honor the caller's BigFloat rounding mode too.
    # Bound exponents so neighbor-to-Rational conversion is a practical oracle.
    cases = Rational{BigInt}[0//1, 1//3, -1//3]
    for sign in (-1, 1)
        push!(cases, sign * dyadic(-200) / (1 + dyadic(-200)))
    end
    for _ = 1:100
        push!(
            cases,
            (BigInt(rand(rng, Int128)) // (BigInt(rand(rng, UInt128)) + 1)) *
            dyadic(rand(rng, -100:100)),
        )
    end
    for bits in (16, 64, 128), mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)
        setprecision(bits) do
            setrounding(BigFloat, mode) do
                for exact in cases
                    actual = HybridPlasmaPIC._from_exact(BigFloat, exact)
                    @test precision(actual) == bits
                    if iszero(exact)
                        @test iszero(actual)
                        continue
                    end
                    value = Rational{BigInt}(actual)
                    previous, following =
                        Rational{BigInt}(prevfloat(actual)), Rational{BigInt}(nextfloat(actual))
                    if mode == RoundNearest
                        distance = abs(exact - value)
                        @test distance <= abs(exact - previous)
                        @test distance <= abs(exact - following)
                        even = iseven(numerator(value / dyadic(exponent(abs(actual)) - bits + 1)))
                        @test (
                            distance != abs(exact - previous) && distance != abs(exact - following)
                        ) || even
                    elseif mode == RoundDown || (mode == RoundToZero && exact > 0)
                        @test value <= exact < following
                    else
                        @test previous < exact <= value
                    end
                end
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
            end
        end
    end
end

@testset "Complete BigFloat conversion in an isolated process" begin
    project = Base.active_project()
    @test project !== nothing
    fixture = joinpath(@__DIR__, "exact_bigfloat_conversion_fixture.jl")
    package = dirname(@__DIR__)
    command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(project)) $fixture $package`
    log = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = log, stderr = log))
    transcript = String(take!(log))
    success(process) || @error "BigFloat conversion fixture failed" transcript
    @test success(process)
    @test occursin("BigFloat conversion settings restored", transcript)
end
