module NormalizationWorkingPrecisionTests
using HybridPlasmaPIC, Test
getemin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
getemax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
setemin(v) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), v)
setemax(v) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), v)
saved = (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
values, units = setprecision(BigFloat, 512) do
    (BigFloat[0, -0.0, Inf, -Inf, NaN], PlasmaUnits{BigFloat}(1, ldexp(BigFloat(1), 40), 1, 1, 1, 1))
end
try
    @assert setemin(-64) == 0
    @assert setemax(64) == 0
    @testset "Complete conversion result uses working BigFloat precision" begin
        for bits in (24, 53, 128, 256),
            mode in (RoundNearest, RoundUp, RoundDown, RoundToZero, RoundFromZero)

            setprecision(BigFloat, bits) do
                setrounding(BigFloat, mode) do
                    for value in values, inverse in (false, true)
                        actual =
                            inverse ? to_normalized(value, :pressure, units) :
                            to_SI(value, :pressure, units)
                        # Finite positive physical scaling preserves special values and signs;
                        # arithmetic creates its result at the current working precision.
                        @test isequal(actual, BigFloat(value))
                        @test precision(actual) == bits
                        @test precision(value) == 512
                    end
                    for value in values
                        z = complex(value, -value)
                        actual = to_SI(z, :pressure, units)
                        @test isequal(actual, complex(BigFloat(real(z)), BigFloat(imag(z))))
                        @test precision(real(actual)) == bits
                        @test precision(imag(actual)) == bits
                        ordinary = to_SI(value, :magnetic, PlasmaUnits{BigFloat}(1, 1, 1, 1, 1, 1))
                        @test isequal(ordinary, BigFloat(value))
                        @test precision(ordinary) == bits
                    end
                end
            end
        end
    end
finally
    @assert setemin(saved[1]) == 0
    @assert setemax(saved[2]) == 0
end
@test saved == (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))

end

module NormalizationDirectedEndpointTests
using HybridPlasmaPIC, Test
getemin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
getemax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
setemin(v) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), v)
setemax(v) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), v)
saved = (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
try
    @assert setemin(4) == 0
    @assert setemax(16) == 0
    @testset "Directed finite endpoints must not certify intermediate range" begin
        setprecision(BigFloat, 24) do
            setrounding(BigFloat, RoundUp) do
                u = PlasmaUnits{BigFloat}(8, 8, 8, 8, 64, 8)
                # The pressure scale is exactly8^2/64=1, below MPFR minimum8.
                # Complete conversion64*1 is nevertheless exactly representable.
                @test to_SI(BigFloat(64), :pressure, u) == BigFloat(64)
                @test to_normalized(BigFloat(64), :pressure, u) == BigFloat(64)
            end
            setrounding(BigFloat, RoundDown) do
                u = PlasmaUnits{BigFloat}(8, 8192, 8, 8, 64, 8)
                # The squared field overflows and rounds to maximum finite. The complete
                # result8*8192^2/64=2^23 must also round down to maximum finite.
                @test to_SI(BigFloat(8), :pressure, u) == floatmax(BigFloat)
            end
        end
    end
finally
    @assert setemin(saved[1]) == 0
    @assert setemax(saved[2]) == 0
end
@test saved == (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))

end

module NormalizationGuardTests
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
getemin() = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
getemax() = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
setemin(v) = ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), v)
setemax(v) = ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), v)
stored(x) = begin
    n, e, d = Base.decompose(x)
    q = BigInt(n) // BigInt(d)
    e >= 0 ? q * (big(1) << e) : q / (big(1) << -e)
end
saved = (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))
wide = setprecision(BigFloat, 512) do
    [
        ldexp(BigFloat(1) + ldexp(BigFloat(1), -400), e) for
        e in (-200, -65, -33, -5, -2, -1, 0, 1, 3, 4, 15, 63, 200)
    ]
end
try
    @testset "BigFloat guard retains strict endpoint semantics" begin
        for (emin, emax) in ((-64, 64), (4, 16), (-32, -4), (-1, 1), (1, 1), (-1, -1))
            @assert setemin(emin) == 0
            @assert setemax(emax) == 0
            for bits in (2, 3, 24, 53, 128, 256),
                mode in (RoundNearest, RoundUp, RoundDown, RoundToZero, RoundFromZero)

                setprecision(BigFloat, bits) do
                    setrounding(BigFloat, mode) do
                        lo, hi = floatmin(BigFloat), floatmax(BigFloat)
                        xs = BigFloat[
                            zero(BigFloat),
                            -zero(BigFloat),
                            BigFloat(Inf),
                            BigFloat(-Inf),
                            BigFloat(NaN),
                        ]
                        append!(xs, wide)
                        append!(
                            xs,
                            (prevfloat(lo), lo, nextfloat(lo), prevfloat(hi), hi, nextfloat(hi)),
                        )
                        for e in (emin - 1, emin, emax - 2, emax - 1)
                            push!(xs, ldexp(BigFloat(1), e))
                        end
                        for x in xs, value in (x, -x)
                            expected = isfinite(value) && lo < value < hi
                            @test H._ordinary_unit_scale(value) == expected
                            expected_complete =
                                !(isfinite(value) && stored(lo) < abs(stored(value)) < stored(hi))
                            @test H._needs_complete_normalization(2, value) == expected_complete
                            @test !H._needs_complete_normalization(0, value)
                            @test !H._needs_complete_normalization(Inf, value)
                        end
                        @test (getemin(), getemax(), precision(BigFloat), rounding(BigFloat)) ==
                              (emin, emax, bits, mode)
                    end
                end
            end
        end
    end
finally
    @assert setemin(saved[1]) == 0
    @assert setemax(saved[2]) == 0
end
@test saved == (getemin(), getemax(), precision(BigFloat), rounding(BigFloat))

end
