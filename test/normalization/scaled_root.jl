module NormalizationScaledRootTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const MODES = (RoundNearest, RoundUp, RoundDown, RoundToZero, RoundFromZero)
dyadic(e) = e >= 0 ? (big(1) << Int(e)) // big(1) : big(1) // (big(1) << Int(-e))
function exact_stored(x)
    n, e, d = Base.decompose(x)
    return (BigInt(n) // BigInt(d)) * dyadic(e)
end
function check_cell(T, q, power, mode, negative)
    actual = H._sqrt_from_exact(T, q, Int128(power), mode, negative)
    @test actual isa T
    @test signbit(actual) == negative
    @test !isnan(actual)
    T === BigFloat && @test precision(actual) == precision(BigFloat)
    exact = q * dyadic(power)
    magnitude = abs(actual)
    nearest = mode === RoundNearest
    away = mode === RoundFromZero || (negative ? mode === RoundDown : mode === RoundUp)
    tiny = exact_stored(nextfloat(zero(T)))
    maximum = exact_stored(floatmax(T))
    if iszero(exact)
        @test iszero(actual)
    elseif isinf(magnitude)
        @test nearest || away
        threshold = maximum + dyadic(exponent(floatmax(T)) - precision(T))
        @test nearest ? exact >= threshold^2 : exact > maximum^2
    elseif iszero(magnitude)
        @test nearest || !away
        @test nearest ? exact <= (tiny / 2)^2 : exact < tiny^2
    else
        value = exact_stored(magnitude)
        previous = exact_stored(prevfloat(magnitude))
        following = nextfloat(magnitude)
        after = isfinite(following) ? exact_stored(following) : dyadic(exponent(floatmax(T)) + 1)
        if nearest
            quantum =
                T === BigFloat ? dyadic(exponent(magnitude) - precision(T) + 1) :
                dyadic(
                    max(
                        exponent(magnitude) - precision(T) + 1,
                        exponent(floatmin(T)) - precision(T) + 1,
                    ),
                )
            even = iseven(numerator(value / quantum))
            low, high = ((previous + value) / 2)^2, ((value + after) / 2)^2
            @test low < exact < high || (even && (low == exact || high == exact))
        elseif away
            @test previous^2 < exact <= value^2
        else
            @test value^2 <= exact && (isinf(following) || exact < after^2)
        end
    end
    return actual
end
function cases(bits, minimum, maximum, subnormal)
    tiny = dyadic(minimum)
    values = Rational{BigInt}[0//1, 1//3, 2//1, 4//1, tiny^2/4]
    for e in unique([minimum, minimum + 1, -1, 0, 1, maximum - 1, maximum]), offset in (0, 1)
        quantum = dyadic(subnormal ? max(e - bits + 1, minimum) : e - bits + 1)
        a = dyadic(e) + offset * quantum
        for boundary in (a, a + quantum / 2), direction in (-1, 0, 1)
            push!(values, boundary^2 + direction * quantum^2 * dyadic(-300) / 3)
        end
    end
    overflow = (dyadic(maximum + 1) - dyadic(maximum - bits))^2
    for direction in (-1, 0, 1)
        push!(values, overflow + direction * dyadic(2maximum - 2bits - 300))
        push!(values, (tiny / 2)^2 + direction * tiny^2 * dyadic(-300))
    end
    push!(values, dyadic(2maximum + 100), dyadic(2minimum - 100))
    return values
end
getemin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
getemax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
setemin(v) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), v)
setemax(v) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), v)
const saved = (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
try
    @testset "Scaled signed square-root exact cells" begin
        @testset "Hardware all rounding modes and split exponents" begin
            for T in (Float16, Float32, Float64)
                qs = cases(precision(T), exponent(nextfloat(zero(T))), exponent(floatmax(T)), true)
                for q in qs, mode in MODES, negative in (false, true), split in (-117, 0, 119)
                    check_cell(T, q / dyadic(split), split, mode, negative)
                end
            end
        end
        @testset "Compact MPFR endpoints and explicit rounding" begin
            for (emin, emax) in ((-64, 64), (4, 16), (-32, -4), (-1, 1))
                @test setemin(emin) == 0
                @test setemax(emax) == 0
                for bits in (2, 3, 24, 53, 128, 256)
                    setprecision(BigFloat, bits) do
                        qs = cases(bits, emin - 1, emax - 1, false)
                        for q in qs,
                            mode in MODES,
                            negative in (false, true),
                            split in (-131, 0, 137)
                            # Caller mode intentionally differs from the explicit requested mode.
                            setrounding(BigFloat, mode === RoundUp ? RoundDown : RoundUp) do
                                before =
                                    (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
                                check_cell(BigFloat, q / dyadic(split), split, mode, negative)
                                @test before == (
                                    getemin(),
                                    getemax(),
                                    precision(BigFloat),
                                    rounding(BigFloat),
                                )
                            end
                        end
                    end
                end
            end
        end
        @test setemin(saved[1]) == 0
        @test setemax(saved[2]) == 0
        @testset "Huge scale endpoints without huge integer allocation" begin
            huge = Int128(typemax(Int64)) * 4
            q = big(2) // big(3)
            for T in (Float16, Float32, Float64, BigFloat), mode in MODES, negative in (false, true)
                away = mode === RoundFromZero || (negative ? mode === RoundDown : mode === RoundUp)
                high = H._sqrt_from_exact(T, q, huge, mode, negative)
                low = H._sqrt_from_exact(T, q, -huge, mode, negative)
                @test signbit(high) == negative && signbit(low) == negative
                @test abs(high) == (mode === RoundNearest || away ? T(Inf) : floatmax(T))
                @test abs(low) == (away ? nextfloat(zero(T)) : zero(T))
            end
            # A huge yet representable exponent never participates in integer shifts.
            e = fld(min(getemax(), -getemin()), 2)
            for bits in (24, 128, 256), direction in (-1, 1)
                setprecision(BigFloat, bits) do
                    actual = H._sqrt_from_exact(
                        BigFloat,
                        big(1) // big(1),
                        Int128(2) * direction * e,
                        RoundNearest,
                        false,
                    )
                    @test actual == ldexp(BigFloat(1), direction * e)
                    @test precision(actual) == bits
                end
            end
        end
    end
finally
    @assert setemin(saved[1]) == 0
    @assert setemax(saved[2]) == 0
end
@test saved == (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
println("Scaled root numerical settings restored")

end
