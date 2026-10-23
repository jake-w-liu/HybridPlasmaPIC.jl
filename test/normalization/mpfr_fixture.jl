using HybridPlasmaPIC, Test
# Exponent limits are mutable native state. Contain them in this child and
# require the active resolved package, including when invoked by Pkg.test.
@assert Threads.nthreads() == 1
@assert length(ARGS) == 1
@assert realpath(pkgdir(HybridPlasmaPIC)) == realpath(only(ARGS))
include("scaled_root.jl")
include("mpfr_values.jl")
@testset "Huge complete conversion cancels exponents before rounding" begin
    maximum = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
    power = fld(maximum, 2) + 16
    for bits in (24, 53, 128, 256), direction in (-1, 1)
        setprecision(BigFloat, bits) do
            value = ldexp(BigFloat(1), direction * power)
            reciprocal = ldexp(BigFloat(1), -direction * power)
            u = PlasmaUnits(n0 = value, B0 = value, mi = value)
            @test to_SI(reciprocal, :pressure, u) ≈ ldexp(inv(u.mu0), direction * power) rtol =
                32eps(BigFloat)
            @test to_normalized(value, :pressure, u) ≈ ldexp(u.mu0, -direction * power) rtol =
                32eps(BigFloat)
        end
    end
end
println("Normalization numerical settings restored")
