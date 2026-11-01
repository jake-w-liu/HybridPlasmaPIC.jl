module NativeDeviceFactorTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC

factor(ps) = H._validated_native_boris_qm(ps)

@testset "Native device quotient admission preserves original zero" begin
    for T in (Float16, Float32, Float64)
        sub = nextfloat(zero(T))
        for (q, m, expected) in (
            (zero(T), sub, zero(T)),
            (-zero(T), sub, -zero(T)),
            (sub, sub, one(T)),
            (-sub, sub, -one(T)),
            (floatmin(T), one(T), floatmin(T)),
            (T(3), T(2), T(1.5)),
        )
            @test isequal(factor(ParticleSet{1,T}(0; q, m)), expected)
        end
        for (q, m) in (
            (sub, T(2)), (-sub, T(2)), # A nonzero charge divided to signed zero.
            (sub, one(T)), (-sub, one(T)), # A nonzero subnormal quotient.
            (one(T), sub), (-one(T), sub), # An overflowing quotient.
            (T(Inf), one(T)), (T(NaN), one(T)),
            (one(T), zero(T)), (one(T), -one(T)),
            (zero(T), T(Inf)), (zero(T), T(NaN)),
        )
            @test_throws ArgumentError factor(ParticleSet{1,T}(0; q, m))
        end
    end
end
end
