module NormalizationRangeTests
using HybridPlasmaPIC, Random, Test
const KINDS = (:length, :time, :velocity, :magnetic, :electric, :density, :current, :pressure)
function reference(x, kind, u, inverse)
    T = promote_type(typeof(x), typeof(u.n0))
    setprecision(BigFloat, 1024) do
        n, B, m, e, mu = map(BigFloat, (u.n0, u.B0, u.mi, u.e, u.mu0))
        scale =
            kind === :length ? sqrt(m / (e^2 * mu * n)) :
            kind === :time ? m / (e * B) :
            kind === :velocity ? B / sqrt(mu * n * m) :
            kind === :magnetic ? B :
            kind === :electric ? B^2 / sqrt(mu * n * m) :
            kind === :density ? n : kind === :current ? e * B * sqrt(n / (mu * m)) : B^2 / mu
        return T(inverse ? BigFloat(x) / scale : BigFloat(x) * scale)
    end
end
@testset "Normalization complete expression and public contracts" begin
    @testset "Stored constructor invariant" begin
        for T in (Float16, Float32, Float64, BigFloat),
            i = 1:6,
            bad in (T(0), -zero(T), -one(T), T(Inf), T(-Inf), T(NaN))

            values = ntuple(j -> j == i ? bad : T(1), 6)
            @test_throws ArgumentError PlasmaUnits{T}(values...)
            @test_throws ArgumentError PlasmaUnits(values...)
        end
        @test_throws MethodError PlasmaUnits(1, 1.0, 1.0, 1.0, 1.0, 1.0)
        for T in (Float16, Float32, Float64, BigFloat)
            u = PlasmaUnits{T}(1, 2, 3, 4, 5, 6)
            @test (u.n0, u.B0, u.mi, u.e, u.mu0, u.c) == map(T, (1, 2, 3, 4, 5, 6))
        end
    end
    rng = MersenneTwister(0x43023)
    @testset "Independent whole formulas across binary exponents" begin
        for T in (Float16, Float32, Float64)
            lo = exponent(floatmin(T)) - precision(T) + 1
            hi = exponent(floatmax(T))
            for iteration = 1:96
                refs = ntuple(_ -> ldexp(T(1), rand(rng, lo:hi)), 6)
                u = PlasmaUnits{T}(refs...)
                for kind in KINDS, inverse in (false, true), sign in (-1, 1)
                    x = sign * ldexp(T(1), rand(rng, lo:hi))
                    actual = inverse ? to_normalized(x, kind, u) : to_SI(x, kind, u)
                    expected = reference(x, kind, u, inverse)
                    @test actual isa T
                    @test signbit(actual) == signbit(expected)
                    @test isapprox(actual, expected; rtol = 16eps(T), atol = zero(T))
                end
            end
        end
    end
    @testset "Exceptional containers, promotion, zeros and nonfinite values" begin
        u = PlasmaUnits{Float64}(1.0, 1e200, 1.0, 1.0, 1.0, 1.0)
        for x in (
            true,
            3,
            Int128(3),
            big(3),
            1 // 3,
            big(1) // big(3),
            Float16(3),
            Float32(3),
            3.0,
            BigFloat(3),
            big(1) << 5000,
            pi,
            MathConstants.e,
        )
            for inverse in (false, true)
                @testset "$(typeof(x)) inverse=$inverse" begin
                    expected = reference(x, :pressure, u, inverse)
                    actual = inverse ? to_normalized(x, :pressure, u) : to_SI(x, :pressure, u)
                    @test typeof(actual) === typeof(expected)
                    @test isapprox(
                        actual,
                        expected;
                        rtol = 32eps(typeof(expected)),
                        atol = zero(expected),
                    )
                end
            end
        end
        for inverse in (false, true), x in (0.0, -0.0, Inf, -Inf, NaN, missing)
            actual = inverse ? to_normalized(x, :pressure, u) : to_SI(x, :pressure, u)
            @test isequal(actual, x)
        end
        for x in (ComplexF64(1e-200, -1e-200), ComplexF64(0.0, -0.0), ComplexF64(Inf, -Inf)),
            inverse in (false, true)

            actual = inverse ? to_normalized(x, :pressure, u) : to_SI(x, :pressure, u)
            expected = complex(
                reference(real(x), :pressure, u, inverse),
                reference(imag(x), :pressure, u, inverse),
            )
            @test isequal(actual, expected)
        end
        values = [1e-200, -1e-200, 0.0, -0.0]
        for x in (
            values,
            reshape(values, 2, 2),
            view(values, 1:2),
            Tuple(values),
            fill(values[1]),
            Float64[],
            Union{Missing,Float64}[missing, values[1]],
        )
            expected = broadcast(y -> ismissing(y) ? missing : reference(y, :pressure, u, false), x)
            actual = to_SI(x, :pressure, u)
            @test typeof(actual) === typeof(expected)
            @test isequal(actual, expected)
        end
        @test_throws ArgumentError to_SI(Float64[], :bad, u)
        @test_throws ArgumentError to_normalized(Float64[], :bad, u)
    end
end

@testset "Stored IEEE finite/nonzero classification" begin
    for encoding = typemin(UInt16):typemax(UInt16)
        x = reinterpret(Float16, encoding)
        @test HybridPlasmaPIC._normalization_finite_nonzero(x) == (isfinite(x) && !iszero(x))
    end
    for T in (Float32, Float64),
        x in (
            zero(T),
            -zero(T),
            nextfloat(zero(T)),
            -nextfloat(zero(T)),
            prevfloat(floatmin(T)),
            floatmin(T),
            floatmax(T),
            T(Inf),
            T(-Inf),
            T(NaN),
        )

        @test HybridPlasmaPIC._normalization_finite_nonzero(x) == (isfinite(x) && !iszero(x))
    end
end
@testset "Finite complete conversions through unrepresentable scales" begin
    for T in (Float32, Float64), sign in (-1, 1), high in (false, true)
        power = T === Float32 ? 30 : 300
        large, small = T(10)^power, T(10)^(-power)
        u =
            high ? PlasmaUnits{T}(large, large, small, 1, 1, 1) :
            PlasmaUnits{T}(small, small, large, 1, 1, 1)
        for kind in (:pressure, :electric, :current, :time), inverse in (false, true)
            scale_large = (high && kind !== :time) || (!high && kind === :time)
            value = sign * (inverse == scale_large ? large : small)
            values = T[value, -value, 0, -0.0]
            expected = reference.(values, kind, Ref(u), inverse)
            actual = inverse ? to_normalized(values, kind, u) : to_SI(values, kind, u)
            @test typeof(actual) === typeof(values)
            @test all(isfinite, expected)
            @test all(isapprox.(actual, expected; rtol = 32eps(T), atol = zero(T)))
            @test signbit.(actual) == signbit.(expected)
        end
    end
end
@testset "Final value operation may overflow after an ordinary scale" begin
    for encoding = UInt16(0x3c01):UInt16(0x3fff)
        B = reinterpret(Float16, encoding)
        u = PlasmaUnits{Float16}(1, B, 1, 1, 1, 1)
        scale = B^2
        isfinite(scale) || continue
        x = Float16(Float64(floatmax(Float16)) / Float64(scale))
        for value in (prevfloat(x), x, nextfloat(x))
            # Three binary16 significands need at most33 bits: this Float64 product
            # is exact, so its one Float16 conversion is an independent reference.
            expected = Float16(Float64(value) * Float64(B)^2)
            actual = to_SI(value, :pressure, u)
            if isfinite(expected)
                @test isfinite(actual)
                @test isapprox(actual, expected; rtol = 2eps(Float16), atol = zero(Float16))
            end
        end
    end
end

end
