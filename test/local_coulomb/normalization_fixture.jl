using HybridPlasmaPIC, Test
@assert realpath(pkgdir(HybridPlasmaPIC)) == realpath(only(ARGS))
const MODES = (RoundNearest, RoundUp, RoundDown, RoundToZero, RoundFromZero)
dyadic(e) = e >= 0 ? (big(1) << e) // big(1) : big(1) // (big(1) << -e)
function stored(x)
    n, e, d = Base.decompose(x)
    (BigInt(n) // BigInt(d)) * dyadic(e)
end
function physical_squares(u)
    n, B, m, e, mu, c = map(stored, (u.n0, u.B0, u.mi, u.e, u.mu0, u.c))
    length2 = m / (mu * n * e^2)
    time2 = m^2 / (e^2 * B^2)
    permittivity2 = 1 / (mu^2 * c^4)
    (permittivity2 * m^2 * length2^3 / (e^4 * time2^2), n^2 * length2^3)
end
function check_cell(y, q, mode)
    T = typeof(y)
    @test isfinite(y) && y > 0
    value = stored(y)
    previous = stored(prevfloat(y))
    following = nextfloat(y)
    after = isfinite(following) ? stored(following) : dyadic(exponent(floatmax(T)) + 1)
    if mode === RoundNearest
        quantum = dyadic(
            T === BigFloat ? exponent(y) - precision(y) + 1 :
            max(exponent(y) - precision(T) + 1, exponent(floatmin(T)) - precision(T) + 1),
        )
        even = iseven(numerator(value / quantum))
        low, high = ((previous + value) / 2)^2, ((value + after) / 2)^2
        @test low < q < high || (even && (q == low || q == high))
    elseif mode === RoundUp || mode === RoundFromZero
        @test previous^2 < q <= value^2
    else
        @test value^2 <= q < after^2
    end
end
function profile(u, mode)
    T = typeof(u.n0)
    before = (precision(BigFloat), rounding(BigFloat))
    q = physical_squares(u)
    invalid = if T === BigFloat
        false
    else
        tiny = stored(nextfloat(zero(T)))
        biggest = stored(floatmax(T))
        overflow = biggest + dyadic(exponent(floatmax(T)) - precision(T))
        any(x -> x <= (tiny / 2)^2 || x >= overflow^2, q)
    end
    if invalid
        @test_throws ArgumentError coulomb_normalization(u)
    else
        result = coulomb_normalization(u)
        @test keys(result) == (:epsilon0, :weight_to_number)
        for (y, square) in zip(values(result), q)
            @test y isa T
            T === BigFloat && @test precision(y) == precision(BigFloat)
            check_cell(y, square, mode)
        end
    end
    @test (precision(BigFloat), rounding(BigFloat)) == before
end
@testset "Coulomb helper independently derived physical rounding cells" begin
    for T in (Float16, Float32, Float64)
        span = T === Float16 ? 12 : T === Float32 ? 100 : 900
        for field = 1:6,
            power in (-span, -1, 0, 1, span),
            fraction in (one(T), T(1.5), nextfloat(one(T)))

            refs = ntuple(i -> i == field ? ldexp(fraction, power) : one(T), 6)
            profile(PlasmaUnits{T}(refs...), RoundNearest)
        end
    end
    for bits in (24, 53, 128, 256), mode in MODES, field = 1:6, power in (-128, -1, 0, 1, 128)
        u = setprecision(BigFloat, 512) do
            refs = ntuple(
                i -> i == field ? ldexp(BigFloat(1) + BigFloat(1) / 7, power) : BigFloat(1),
                6,
            )
            PlasmaUnits{BigFloat}(refs...)
        end
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                profile(u, mode)
            end
        end
    end
end
getemin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
getemax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
setemin(x) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), x)
setemax(x) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), x)
saved = (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
try
    @testset "Coulomb helper compact MPFR boundaries" begin
        for bits in (24, 256), mode in MODES, power in (-10000, 10000)
            u = setprecision(BigFloat, bits) do
                PlasmaUnits{BigFloat}(
                    ldexp(BigFloat(1), 6power),
                    ldexp(BigFloat(1), 2power),
                    1,
                    1,
                    ldexp(BigFloat(1), -2power),
                    1,
                )
            end
            try
                @test setemin(-32) == 0
                @test setemax(32) == 0
                setprecision(BigFloat, bits) do
                    setrounding(BigFloat, mode) do
                        result = coulomb_normalization(u)
                        @test result.epsilon0 == 1 && result.weight_to_number == 1
                        @test precision(result.epsilon0) == bits &&
                              precision(result.weight_to_number) == bits
                        @test (getemin(), getemax(), rounding(BigFloat)) == (-32, 32, mode)
                    end
                end
            finally
                setemin(saved[1])
                setemax(saved[2])
            end
        end
    end
finally
    setemin(saved[1])
    setemax(saved[2])
end
@test (getemin(), getemax(), precision(BigFloat), rounding(BigFloat)) == saved

println("Local Coulomb normalization settings restored")
