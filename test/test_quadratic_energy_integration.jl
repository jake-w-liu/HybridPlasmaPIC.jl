using HybridPlasmaPIC, Test

# KineticRangeOracle supplies independent rational rounding-cell checks.
@testset "Complete particle, magnetic and stored-pressure energy in one budget" begin
    Q = Rational{BigInt}
    for T in (Float32, Float64), D = 1:3
        # The individual squares and pressure sum overflow, while all complete
        # reservoirs and their total are finite. Powers of two isolate range
        # handling from the ordinary summation error tested elsewhere.
        length_power = fld(T === Float32 ? -50 : -500, D)
        speed_power = T === Float32 ? 70 : 600
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> ldexp(one(T), length_power), D))
        volume = prod(Q, g.L)
        ps = ParticleSet{D,T}(1; m = T(volume))
        ps.v[1][1] = ldexp(one(T), speed_power)
        B = ntuple(_ -> zeros(T, g.n), 3)
        fill!(B[1], ps.v[1][1])
        n = ones(T, g.n)
        pe = fill(ldexp(one(T), exponent(floatmax(T))), g.n)
        closure = EvolvingElectrons(T(2))
        exact_kinetic = Q(ps.m) * Q(ps.v[1][1])^2 / 2
        exact_magnetic = sum(x -> Q(x)^2, B[1]) * prod(Q, g.dx) / 2
        exact_electron = sum(Q, pe) * prod(Q, g.dx)
        exact_total = exact_kinetic + exact_magnetic + exact_electron
        before = map(copy, (ps.x..., ps.v..., ps.weight, B..., n, pe))
        @test isinf(ps.v[1][1]^2)
        @test isinf(sum(pe))
        @test 0 < exact_total < Q(floatmax(T))
        for particles in (ps, [ps], view([ps], :))
            budget = energy_budget(particles, B, n, closure, g; pe)
            @test isfinite(budget.total)
            @test KineticRangeOracle.rounded(budget.kinetic, exact_kinetic)
            @test KineticRangeOracle.rounded(budget.magnetic, exact_magnetic)
            @test KineticRangeOracle.rounded(budget.electron_internal, exact_electron)
            @test budget.electron_inertia === zero(T)
            @test KineticRangeOracle.rounded(budget.total, exact_total)
        end
        @test all(isequal.(before, (ps.x..., ps.v..., ps.weight, B..., n, pe)))
    end
end
