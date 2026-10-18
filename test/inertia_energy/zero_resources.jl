module ZeroQuadraticEnergyTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
budget_cost(ps, B, n, c, g) = @timed energy_budget(ps, B, n, c, g; de2 = eltype(n)(0.25))
field_cost(es) = @allocated em_field_energy(es)
magnetic_cost(B, g) = @allocated magnetic_energy(B, g)
@testset "Zero current budgets retain ordinary allocation at overflowing volume" begin
    for T in (Float32, Float64), count in (4, 16, 64)
        side = ldexp(one(T), T === Float32 ? 50 : 400)
        large = FourierGrid(ntuple(_ -> count, 3), ntuple(_ -> side, 3))
        ordinary = FourierGrid(large.n, ntuple(_ -> T(2pi), 3))
        B = ntuple(_ -> zeros(T, large.n), 3)
        n = ones(T, large.n)
        ps = ParticleSet{3,T}(0)
        c = IsothermalElectrons(zero(T))
        @test isinf(prod(large.dx))
        budget_cost(ps, B, n, c, ordinary)
        budget_cost(ps, B, n, c, large)
        reference = [budget_cost(ps, B, n, c, ordinary) for _ = 1:3]
        samples = [budget_cost(ps, B, n, c, large) for _ = 1:3]
        for (a, b) in zip(reference, samples)
            @test all(iszero, b.value)
            @test isequal(a.value, b.value)
        end
        count >= 16 &&
            (@test maximum(s.bytes for s in samples) <= 1.1minimum(s.bytes for s in reference))
        magnetic_cost(B, large)
        @test magnetic_cost(B, large) == 0
    end
end
@testset "Zero-field certificate distinguishes lost nonzero squares and invalid factors" begin
    for T in (Float32, Float64), D = 1:3
        e = T === Float32 ? 80 : 600
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
        for legacy in (false, true)
            legacy && D != 1 && continue
            es = legacy ? EMPIC1D(g, 0; c = ldexp(one(T), e)) : EMPIC(g, 0; c = ldexp(one(T), e))
            @test em_field_energy(es) === zero(T)
            field_cost(es)
            @test field_cost(es) == 0
            other =
                legacy ? EMPIC1D(g, 0; c = ldexp(one(T), -e)) : EMPIC(g, 0; c = ldexp(one(T), -e))
            fields = legacy ? (other.Bz,) : other.B
            fill!(first(fields), ldexp(one(T), e))
            # c² rounds to zero and B² overflows, but their complete product is one.
            value = em_field_energy(other)
            @test isfinite(value) && value > 0
            @test value ≈ T(0.5) * prod(g.L)
        end
        z = ntuple(_ -> zeros(T, g.n), 3)
        @test isnan(H._electron_inertia_energy(z, T(Inf), g.dx))
        @test isnan(H._em_field_energy_finish(z, z, T(Inf), g.dx, T(NaN), T(Inf), T(Inf)))
        @test isnan(H._em_field_energy_finish(z, z, one(T), (T(Inf),), zero(T), T(Inf), T(Inf)))
        @test isequal(
            H._em_field_energy_finish(z, z, one(T), g.dx, -zero(T), T(Inf), T(Inf)),
            -zero(T),
        )
    end
end
end
