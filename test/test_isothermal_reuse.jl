module IsothermalReuseChecks
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
uncached(n, te, dx) = H._isothermal_general_probe_with_reuse(n, te, dx, Val(false))
@testset "Repeated active runs require at least one hit per miss" begin
    for T in (Float32, Float64)
        for n in (T[1.5, 1.5], T[0.5, 0, 1, -1, 0.5], T[0.5, 0.5, 2, 2])
            @test H._isothermal_repeated_active_runs(n)
        end
        for n in (T[], T[0, 1, -1], T[0.5], T[0.5, 0.5, 2], T[0.5, 2, 2], T[0.5, 0.5, 2, 2, 0.75])
            @test !H._isothermal_repeated_active_runs(n)
        end
    end
end
@testset "Per-density reuse retains complete probe states bit for bit" begin
    rng = MersenneTwister(91731)
    for T in (Float32, Float64), C in (Float32, Float64), D = 1:3
        arrays = [
            fill(T(1.5), 64),
            T[0.5, 0, 1, -1, 0.5, 0.75, 0.75, 1, 0.5],
            repeat(T[0.75, 0.75, 1.5, 1.5, 2, 2], 8),
            T(0.5) .+ T(2) .* rand(rng, T, 64),
            T[NaN, 0.5, 0.5, -Inf, 1.5],
            T[Inf, Inf, 0.75, 0.75],
            T[floatmax(T), floatmax(T), nextfloat(zero(T)), nextfloat(zero(T))],
            T[],
        ]
        spacings = [
            ntuple(_ -> one(T), D),
            ntuple(i -> isodd(i) ? floatmax(T) : floatmin(T), D),
            ntuple(_ -> nextfloat(zero(T)), D),
            ntuple(_ -> T(Inf), D),
        ]
        for n in arrays,
            te in
            (zero(C), C(0.7), -C(0.7), floatmax(C), nextfloat(zero(C)), C(NaN), C(Inf), -C(Inf)),
            dx in spacings

            before = copy(n)
            @test isequal(H._isothermal_general_probe(n, te, dx), uncached(n, te, dx))
            @test isequal(
                H._isothermal_general_probe_with_reuse(n, te, dx, Val(true)),
                uncached(n, te, dx),
            )
            @test isequal(
                H._isothermal_general_probe(view(n, :), te, dx),
                uncached(view(n, :), te, dx),
            )
            @test isequal(n, before)
        end
    end
end
end
