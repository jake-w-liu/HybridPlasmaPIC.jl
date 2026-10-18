using HybridPlasmaPIC, Test

# Independent rational oracles are loaded by the preceding diagnostic tests.
@testset "Complete CGL, particle and magnetic reservoirs in one budget" begin
    Q = Rational{BigInt}
    for T in (Float32, Float64), D = 1:3
        length_power = fld(T === Float32 ? -50 : -500, D)
        amplitude_power = T === Float32 ? 70 : 600
        amplitude = ldexp(one(T), amplitude_power)
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> ldexp(one(T), length_power), D))
        volume = prod(Q, g.L)
        ps = ParticleSet{D,T}(1; m = T(volume))
        ps.v[1][1] = amplitude
        B = (fill(amplitude, g.n), zeros(T, g.n), zeros(T, g.n))
        n = fill(amplitude, g.n)
        exact_kinetic = KineticRangeOracle.exact((ps,))
        exact_magnetic = sum(x -> Q(x)^2, B[1]) * prod(Q, g.dx) / 2
        @test isinf(amplitude^2)
        for (perpendicular, parallel) in
            ((one(T), zero(T)), (zero(T), amplitude), (one(T), amplitude))
            closure = CGLElectrons(perpendicular, parallel, one(T), one(T))
            lo, hi = CGLEnergyRangeOracle.interval(n, B, closure, g.dx)
            @test lo == hi
            exact_total = exact_kinetic + exact_magnetic + lo
            @test 0 < exact_total < Q(floatmax(T))
            ordinary = CGLEnergyRangeOracle.ordinary(n, B, closure, g.dx)
            @test !isfinite(ordinary) || iszero(ordinary)
            before = map(copy, (ps.x..., ps.v..., ps.weight, B..., n))
            for particles in (ps, [ps], view([ps], :))
                budget = energy_budget(particles, B, n, closure, g)
                @test isfinite(budget.total)
                @test KineticRangeOracle.rounded(budget.kinetic, exact_kinetic)
                @test KineticRangeOracle.rounded(budget.magnetic, exact_magnetic)
                @test CGLEnergyRangeOracle.rounded(budget.electron_internal, lo, hi)
                @test budget.electron_inertia === zero(T)
                # All terms are exact powers of two or their small integer sums.
                @test KineticRangeOracle.rounded(budget.total, exact_total)
            end
            @test all(isequal.(before, (ps.x..., ps.v..., ps.weight, B..., n)))
        end
    end
end
