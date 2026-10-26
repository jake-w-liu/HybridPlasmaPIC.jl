using HybridPlasmaPIC, Test

# MPFR exponent limits are separate from task-scoped precision/rounding. Keep
# the test's temporary range changes inside this dedicated Julia process.
@testset "BigFloat narrowed exponent-range admission and restoration" begin
    H = HybridPlasmaPIC
    B = H._TimestepBounds
    g = FourierGrid((4,), (4π,))
    tiny_grid = FourierGrid((4,), (4e-100,))
    eta = BigFloat(2)
    outside = ldexp(BigFloat(1), -100)
    full = recommended_dt(g; η = eta, safety = 1.0)
    original_emin, original_emax = Base.MPFR.get_emin(), Base.MPFR.get_emax()
    original_precision, original_rounding = precision(BigFloat), rounding(BigFloat)
    try
        Base.MPFR.set_emin!(-3)
        Base.MPFR.set_emax!(128)
        step = recommended_dt(g; η = eta, safety = 1.0)
        @test isfinite(step) && step > 0
        @test step <= full
        @test_throws ArgumentError recommended_dt(g; η = outside)
        @test_throws ArgumentError B(outside)
        @test_throws ArgumentError B(1e-100)
        setprecision(BigFloat, 11) do
            setrounding(BigFloat, RoundUp) do
                # All keyword values pass admission; this error occurs inside
                # the bounded core when native grid spacing converts outside
                # the active MPFR range.
                @test_throws ArgumentError recommended_dt(tiny_grid; η = eta)
                @test precision(BigFloat) == 11
                @test rounding(BigFloat) == RoundUp
            end
        end
        @test Base.MPFR.get_emin() == -3
        @test Base.MPFR.get_emax() == 128
        Base.MPFR.set_emin!(2)
        @test_throws ArgumentError recommended_dt(g; η = eta)
    finally
        Base.MPFR.set_emin!(original_emin)
        Base.MPFR.set_emax!(original_emax)
    end
    @test Base.MPFR.get_emin() == original_emin
    @test Base.MPFR.get_emax() == original_emax
    @test precision(BigFloat) == original_precision
    @test rounding(BigFloat) == original_rounding
end
