using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

# MPFR exponent bounds are mutable native state (thread-local in TLS builds).
# Exercise them only in this dedicated single-thread process, and restore them
# even when a test fails.
@assert Threads.nthreads() == 1
@assert length(ARGS) == 1
@assert realpath(pkgdir(H)) == realpath(only(ARGS))

dyadic(e) = e >= 0 ? (BigInt(1) << e) // BigInt(1) : BigInt(1) // (BigInt(1) << -e)
function stored_rational(value::BigFloat)
    significand, power, sign = Base.decompose(value)
    return BigInt(sign) * significand * dyadic(power)
end
get_emin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
get_emax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
set_emin(value) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), value)
set_emax(value) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), value)

# An independent rounding-cell certificate: no integer quotient or rational to
# floating conversion is used. Adjacent representable values and their exact
# midpoints determine the answer, including the zero and infinity endpoints.
function in_rounding_cell(exact, actual, bits, emin, emax, mode)
    precision(actual) == bits || return false
    isnan(actual) && return false
    iszero(exact) && return iszero(actual) && !signbit(actual)
    signbit(actual) == (exact < 0) || return false
    magnitude = abs(exact)
    minimum = dyadic(emin - 1)
    largest = ((BigInt(1) << bits) - 1) * dyadic(emax - bits)
    away =
        mode == RoundFromZero || (mode == RoundUp && exact > 0) || (mode == RoundDown && exact < 0)
    if iszero(actual)
        return mode == RoundNearest ? magnitude <= minimum / 2 : !away && magnitude < minimum
    elseif isinf(actual)
        return mode == RoundNearest ? magnitude >= largest + dyadic(emax - bits - 1) :
               away && magnitude > largest
    end
    value = abs(stored_rational(actual))
    minimum <= value <= largest || return false
    power = exponent(actual)
    step = dyadic(power - bits + 1)
    previous = value == minimum ? zero(value) : value - (value == dyadic(power) ? step / 2 : step)
    following = value + step
    if mode == RoundNearest
        lower, upper = (previous + value) / 2, (value + following) / 2
        even = iseven(numerator(value / step))
        # At half the minimum MPFR rounds to zero. At the largest finite
        # value the conceptual next binade gives the overflow midpoint.
        lower_ok = lower < magnitude || (even && value != minimum && lower == magnitude)
        upper_ok = magnitude < upper || (even && magnitude == upper)
        return lower_ok && upper_ok
    elseif away
        return previous < magnitude <= value
    else
        return value <= magnitude && (value == largest || magnitude < following)
    end
end

function exact_cases(bits, emin, emax)
    minimum = dyadic(emin - 1)
    largest = ((BigInt(1) << bits) - 1) * dyadic(emax - bits)
    # Far-side perturbations retain a tie's direction even when the integer
    # factors cannot individually fit in the active exponent range.
    delta = dyadic(emin - bits - 260) / 3
    cases = Rational{BigInt}[0//1, 1//3, 1//1, minimum/8, 4*largest]
    boundaries = Rational{BigInt}[minimum/2, minimum, largest, largest+dyadic(emax - bits - 1)]
    for power in unique((emin - 1, clamp(0, emin - 1, emax - 1), emax - 1))
        step = dyadic(power - bits + 1)
        for significand in unique((
            (BigInt(1) << (bits - 1)),
            (BigInt(1) << (bits - 1)) + 1,
            (BigInt(1) << bits) - 2,
        ))
            push!(boundaries, (significand + 1 // 2) * step)
        end
    end
    for boundary in boundaries, side in (-1, 0, 1)
        push!(cases, boundary + side * delta)
    end
    for width in (200, 600)
        large = BigInt(1) << width
        push!(cases, (large + 1) // (large - 1))
        push!(cases, (large - 1) // (large + 1))
        push!(cases, minimum * ((large + 1) // (large - 1)))
    end
    rng = MersenneTwister(1129 + bits)
    for _ = 1:20
        numerator = (BigInt(rand(rng, UInt128)) << 128) + rand(rng, UInt128) + 1
        denominator = (BigInt(rand(rng, UInt128)) << 128) + rand(rng, UInt128) + 1
        push!(cases, (numerator // denominator) * dyadic(rand(rng, (emin-4):(emax+4))))
    end
    return unique(vcat(cases, -cases))
end

@testset "BigFloat complete quotient rounding cells and compact endpoints" begin
    saved_emin, saved_emax = get_emin(), get_emax()
    saved_precision, saved_rounding = precision(BigFloat), rounding(BigFloat)
    try
        for bits in (2, 3, 24, 53, 128, 256),
            (emin, emax) in ((-64, 64), (4, 16), (-32, -4), (-1, 1))

            cases = exact_cases(bits, emin, emax)
            @test set_emin(saved_emin) == 0
            @test set_emax(saved_emax) == 0
            @test set_emin(emin) == 0
            @test set_emax(emax) == 0
            setprecision(BigFloat, bits) do
                for mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)
                    setrounding(BigFloat, mode) do
                        for exact in cases
                            actual = H._from_exact(BigFloat, exact)
                            @test in_rounding_cell(exact, actual, bits, emin, emax, mode)
                        end
                        @test (get_emin(), get_emax()) == (emin, emax)
                        @test precision(BigFloat) == bits
                        @test rounding(BigFloat) == mode
                    end
                end
            end
        end
        @test set_emin(saved_emin) == 0
        @test set_emax(saved_emax) == 0
        dx, origin, samples = setprecision(BigFloat, 256) do
            ldexp(BigFloat(1), -50),
            BigFloat(1) + ldexp(BigFloat(1), -100),
            [BigFloat(1) + ldexp(BigFloat(k), -100) for k = 1:3]
        end
        center_exact = stored_rational(origin) + 3 * stored_rational(dx) / 2
        mean_exact = sum(stored_rational, samples) / 3
        grid = H.AMRGrid{BigFloat}(copy(samples), dx, origin)
        @test set_emin(-64) == 0
        @test set_emax(64) == 0
        for bits in (24, 128, 256),
            mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero)

            setprecision(BigFloat, bits) do
                setrounding(BigFloat, mode) do
                    center = H.cell_center(grid, 2)
                    @test in_rounding_cell(center_exact, center, bits, -64, 64, mode)
                    calls = Ref(0)
                    sample(_) = samples[(calls[]+=1)]
                    mean = H.gyroaverage(
                        sample,
                        (BigFloat(0), BigFloat(0), BigFloat(0)),
                        BigFloat(0),
                        (BigFloat(0), BigFloat(0), BigFloat(1));
                        n = 3,
                    )
                    @test calls[] == 3
                    @test in_rounding_cell(mean_exact, mean, bits, -64, 64, mode)
                    @test stored_rational(grid.x0) == stored_rational(origin)
                    @test (get_emin(), get_emax()) == (-64, 64)
                end
            end
        end
    finally
        @test set_emin(saved_emin) == 0
        @test set_emax(saved_emax) == 0
    end
    @test (get_emin(), get_emax()) == (saved_emin, saved_emax)
    @test precision(BigFloat) == saved_precision
    @test rounding(BigFloat) == saved_rounding
end
println("BigFloat conversion settings restored")
