using Test
using HybridPlasmaPIC

@testset "Implicit hybrid completed clock admission" begin
    for T in (Float32, Float64)
        report =
            (time, count, h) -> ImplicitHybridSolveReport(
                count,
                time,
                h,
                zero(T),
                one(T),
                zero(T),
                zero(T),
                one(T),
                one(T),
                0,
                0,
                2,
                0,
                zero(T),
                zero(T),
                zero(T),
                zero(T),
                zero(T),
                :represented_reservoir,
            )
        phase = HybridPlasmaPIC._implicit_hybrid_phase_values
        @test phase(true, one(T), 1, one(T), report(one(T), 1, one(T))) === nothing
        # Varying prior durations need not give time = step * last_dt.
        @test phase(true, T(3 / 2), 2, one(T), report(T(3 / 2), 2, one(T))) === nothing
        @test_throws ArgumentError phase(true, one(T), 1, T(2), report(one(T), 1, T(2)))
    end
end
