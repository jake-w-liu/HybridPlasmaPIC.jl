using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
@assert Threads.nthreads() == 1
@assert realpath(pkgdir(H)) == realpath(only(ARGS))
dyadic(e) = e >= 0 ? (big(1) << e) // big(1) : big(1) // (big(1) << -e)
function stored(value)
    n, e, d = Base.decompose(value)
    return (BigInt(n) // BigInt(d)) * dyadic(e)
end
get_emin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
get_emax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
set_emin(x) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), x)
set_emax(x) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), x)

# A component's exact rational must belong to the returned float's rounding
# cell. Construct neighbors as integer dyadics, including MPFR endpoints.
function cell(exact, result, p, emin, emax, mode)
    precision(result) == p && !isnan(result) || return false
    iszero(exact) && return iszero(result)
    signbit(result) == (exact < 0) || return false
    x = abs(exact)
    minimum = dyadic(emin - 1)
    maximum = ((big(1) << p) - 1) * dyadic(emax - p)
    away = mode === RoundFromZero || (exact < 0 ? mode === RoundDown : mode === RoundUp)
    if iszero(result)
        return mode === RoundNearest ? x <= minimum / 2 : !away && x < minimum
    elseif isinf(result)
        return mode === RoundNearest ? x >= maximum + dyadic(emax - p - 1) : away && x > maximum
    end
    y = abs(stored(result))
    power = exponent(result)
    step = dyadic(power - p + 1)
    previous = y == minimum ? zero(y) : y - (y == dyadic(power) ? step / 2 : step)
    following = y + step
    if mode === RoundNearest
        even = iseven(numerator(y / step))
        low, high = (previous + y) / 2, (y + following) / 2
        return (low < x || (even && y != minimum && low == x)) && (x < high || (even && x == high))
    elseif away
        return previous < x <= y
    else
        return y <= x && (y == maximum || x < following)
    end
end

saved_emin, saved_emax = get_emin(), get_emax()
saved_precision, saved_rounding = precision(BigFloat), rounding(BigFloat)
@testset "BigFloat exact intermediates and bounded-path transitions" begin
    rng = MersenneTwister(72053)
    for (pa, pb) in ((2, 3), (24, 53), (53, 128), (256, 256))
        # Cross the admissibility boundary with unrelated stored significands,
        # result precision and signs; integer rationals supply the cell oracle.
        width = pa + pb
        for power in (-width - 1, -width, -width + 1, -1, 0, 1, width - 1, width, width + 1)
            a = setprecision(BigFloat, pa) do
                ldexp(1 + rand(rng, BigFloat), fld(power, 2))
            end
            b = setprecision(BigFloat, pb) do
                ldexp(1 + rand(rng, BigFloat), power - fld(power, 2))
            end
            for p in (2, 24, 53, 256),
                mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)

                setprecision(BigFloat, p) do
                    setrounding(BigFloat, mode) do
                        for (x, y) in ((a, b), (b, a), (-a, b))
                            h = stored(x) * stored(y) / 2
                            expected = ((1 - h^2) / (1 + h^2), -2h / (1 + h^2))
                            before = (stored(x), stored(y), precision(x), precision(y))
                            actual = cn_multiplier(x, y)
                            # These inputs/results are strictly inside this
                            # bounded oracle range; no default MPFR endpoint
                            # needs expansion into a billion-bit integer.
                            @test cell(expected[1], real(actual), p, -4096, 4096, mode)
                            @test cell(expected[2], imag(actual), p, -4096, 4096, mode)
                            @test (stored(x), stored(y), precision(x), precision(y)) == before
                        end
                    end
                end
            end
        end
    end
end
@testset "CN exact component rounding in compact exponent ranges" begin
    try
        for p in (2, 3, 24, 53, 128), (emin, emax) in ((-64, 64), (4, 16), (-32, -4), (-1, 1))
            @assert set_emin(saved_emin) == 0
            @assert set_emax(saved_emax) == 0
            # Inputs retain more precision than the result, and are constructed
            # before changing the exponent interval, as allowed by scalar APIs.
            inputs = setprecision(BigFloat, 256) do
                v = BigFloat[
                    0,
                    ldexp(BigFloat(1), -70),
                    ldexp(BigFloat(1), -32),
                    0.125,
                    0.5,
                    1,
                    2,
                    4,
                    ldexp(BigFloat(1), 32),
                    ldexp(BigFloat(1), 63),
                ]
                push!(v, BigFloat(1) + ldexp(BigFloat(1), -80))
                push!(v, BigFloat(2) - ldexp(BigFloat(1), -79))
                vcat(v, -v)
            end
            @assert set_emin(emin) == 0
            @assert set_emax(emax) == 0
            setprecision(BigFloat, p) do
                for mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)
                    setrounding(BigFloat, mode) do
                        for a in inputs, b in inputs
                            h = stored(a) * stored(b) / 2
                            expected = ((1 - h^2) / (1 + h^2), -2h / (1 + h^2))
                            actual = cn_multiplier(a, b)
                            @test cell(expected[1], real(actual), p, emin, emax, mode)
                            @test cell(expected[2], imag(actual), p, emin, emax, mode)
                        end
                        @test (get_emin(), get_emax(), precision(BigFloat), rounding(BigFloat)) ==
                              (emin, emax, p, mode)
                    end
                end
            end
        end
    finally
        @assert set_emin(saved_emin) == 0
        @assert set_emax(saved_emax) == 0
    end
end

@testset "Exact-intermediate admission at compact exponent boundaries" begin
    try
        for (pa, pb) in ((2, 3), (3, 5)), (emin, emax) in ((-64, 64), (-16, 16), (-8, 8))
            @assert set_emin(saved_emin) == 0
            @assert set_emax(saved_emax) == 0
            inputs = [
                (
                    BigFloat(sign * 1.5; precision = pa) * BigFloat(2)^fld(power, 2),
                    BigFloat(1.75; precision = pb) * BigFloat(2)^(power - fld(power, 2)),
                ) for power = -12:12, sign in (-1, 1)
            ]
            # Scaling retains the intended stored input precisions, independent
            # of the eventual output precision and compact exponent interval.
            inputs =
                [(BigFloat(a; precision = pa), BigFloat(b; precision = pb)) for (a, b) in inputs]
            @assert set_emin(emin) == 0
            @assert set_emax(emax) == 0
            for p in (2, 24, 53),
                mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)

                setprecision(BigFloat, p) do
                    setrounding(BigFloat, mode) do
                        for (a, b) in inputs
                            h = stored(a) * stored(b) / 2
                            expected = ((1 - h^2) / (1 + h^2), -2h / (1 + h^2))
                            actual = cn_multiplier(a, b)
                            @test cell(expected[1], real(actual), p, emin, emax, mode)
                            @test cell(expected[2], imag(actual), p, emin, emax, mode)
                        end
                        @test (get_emin(), get_emax(), precision(BigFloat), rounding(BigFloat)) ==
                              (emin, emax, p, mode)
                    end
                end
            end
        end
    finally
        @assert set_emin(saved_emin) == 0
        @assert set_emax(saved_emax) == 0
    end
end

@testset "Huge binary exponents remain separate from significands" begin
    power = max(fld(saved_emax, 2) + 2, fld(-saved_emin, 2) + 3)
    for p in (2, 24, 53, 256)
        setprecision(BigFloat, p) do
            for mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)
                setrounding(BigFloat, mode) do
                    for sign in (-1, 1), direction in (-1, 1)
                        a = ldexp(BigFloat(sign), direction * power)
                        b = ldexp(BigFloat(1), direction * power)
                        @test isfinite(a) && isfinite(b) && !iszero(a) && !iszero(b)
                        actual = cn_multiplier(a, b)
                        real_expected = if direction > 0
                            mode in (RoundUp, RoundToZero) ? nextfloat(-one(BigFloat)) : -one(BigFloat)
                        else
                            mode in (RoundDown, RoundToZero) ? prevfloat(one(BigFloat)) : one(BigFloat)
                        end
                        away =
                            mode === RoundFromZero ||
                            (sign > 0 ? mode === RoundDown : mode === RoundUp)
                        imaginary_expected = away ? nextfloat(zero(BigFloat)) : zero(BigFloat)
                        sign > 0 && (imaginary_expected = -imaginary_expected)
                        @test isequal(real(actual), real_expected)
                        @test isequal(imag(actual), imaginary_expected)
                    end
                    a = ldexp(BigFloat(1), power)
                    b = ldexp(BigFloat(1), -power)
                    @test cn_multiplier(a, b) == complex(BigFloat(3) / 5, BigFloat(-4) / 5)
                    @test cn_multiplier(b, a) == cn_multiplier(a, b)
                    @test (get_emin(), get_emax(), precision(BigFloat), rounding(BigFloat)) ==
                          (saved_emin, saved_emax, p, mode)
                end
            end
        end
    end
end
@test (get_emin(), get_emax(), precision(BigFloat), rounding(BigFloat)) ==
      (saved_emin, saved_emax, saved_precision, saved_rounding)
println("BigFloat CN settings restored")
