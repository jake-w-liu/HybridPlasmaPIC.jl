using HybridPlasmaPIC, Test, Random

# Uses KineticRangeOracle from test_kinetic_energy_range.jl.
@testset "Kinetic small normalized weights keep ordinary outputs and allocations" begin
    rng = MersenneTwister(249852)
    for T in (Float32, Float64),
        count in (16, 256),
        power in (-5, -10, -20),
        mass in (one(T), T(1 / 1836), -one(T)),
        mixed in (false, true)

        ps = ParticleSet{1,T}(count; m = mass)
        ps.weight .= T(10)^power
        mixed && (ps.weight[1:2:end] .*= -one(T))
        for v in ps.v
            v .= T(0.2) .+ T(0.2) .* rand(rng, T, count)
        end
        expected = KineticRangeOracle.ordinary(ps)
        @test isfinite(expected) && abs(expected) >= floatmin(T)
        view_ps = ParticleSet{1,T}(
            map(a -> view(a, :), ps.x),
            map(a -> view(a, :), ps.v),
            view(ps.weight, :),
            view(ps.id, :),
            view(ps.tag, :),
            ps.q,
            ps.m,
        )
        before = map(copy, (ps.v..., ps.weight))
        for arg in (ps, [ps], view_ps)
            @test isequal(kinetic_energy(arg), expected)
            @test kinetic_range_alloc(arg) == 0
        end
        @test all(isequal.(before, (ps.v..., ps.weight)))
    end
end

@testset "Kinetic direct operation enclosures vs complete rational sums" begin
    Q = Rational{BigInt}
    rng = MersenneTwister(128572)
    accepted = 0
    checked = 0
    for T in (Float16, Float32, Float64), count in (1, 3, 16, 129), trial = 1:20
        wp =
            T === Float16 ? (-5, -2, 0, 2) :
            T === Float32 ? (-100, -65, -10, 0, 20, 60) : (-500, -200, -20, 0, 200, 500)
        mp = T === Float16 ? (-2, 0, 2) : T === Float32 ? (-20, 0, 20) : (-200, 0, 200)
        vp = T === Float16 ? (-1, 0, 1) : T === Float32 ? (-10, -1, 2) : (-50, -1, 50)
        populations = [ParticleSet{1,T}(count ÷ 2), ParticleSet{1,T}(count - count ÷ 2)]
        observed = zero(T)
        absolute_exact = Q(0)
        for (j, ps) in enumerate(populations)
            ps.m = ldexp(one(T), rand(rng, mp)) * (isodd(trial + j) ? -one(T) : one(T))
            for i in eachindex(ps.weight)
                ps.weight[i] = ldexp(T(0.75) + rand(rng, T) / T(2), rand(rng, wp))
                isodd(trial) && rand(rng, Bool) && (ps.weight[i] = -ps.weight[i])
                for v in ps.v
                    v[i] = ldexp(T(0.75) + rand(rng, T) / T(2), rand(rng, vp))
                    absolute_exact += abs(Q(ps.m) * Q(ps.weight[i]) * Q(v[i])^2 / 2)
                end
            end
            observed += KineticRangeOracle.ordinary(ps)
        end
        absolute, normal = HybridPlasmaPIC._kinetic_energy_absolute_probe(populations, T, Val(true))
        exact = KineticRangeOracle.exact(populations)
        certified = HybridPlasmaPIC._kinetic_energy_operations_certified(
            populations,
            observed,
            Int128(count),
        )
        @test !certified || normal
        if normal
            checked += 1
            u = Q(1) / (BigInt(1) << precision(T))
            delta = (count + 4) * u
            gamma = delta / (1 - delta)
            theoretical = delta / (1 - 2delta) * Q(absolute) + count * Q(nextfloat(zero(T)))
            @test abs(Q(absolute) - absolute_exact) <= gamma * absolute_exact
            @test abs(Q(observed) - exact) <= theoretical
        end
        if certified
            accepted += 1
            @test Q(floatmin(T)) <= abs(exact) <= Q(floatmax(T))
            @test isequal(kinetic_energy(populations), observed)
        end
    end
    @test checked > 0
    @test accepted > 0
end

@testset "Kinetic signed population underflow stays inside the operation enclosure" begin
    Q = Rational{BigInt}
    for T in (Float16, Float32, Float64)
        small = ParticleSet{1,T}(2; m = T(6) * floatmin(T))
        small.v[1] .= one(T)
        small.weight .= (one(T), -prevfloat(one(T)))
        ordinary_small = KineticRangeOracle.ordinary(small)
        @test zero(T) < ordinary_small < floatmin(T)
        @test !isequal(Q(ordinary_small), KineticRangeOracle.exact((small,)))
        normal = ParticleSet{1,T}(1; m = T(2))
        normal.v[1][1] = one(T)
        populations = [small, normal]
        observed = ordinary_small + KineticRangeOracle.ordinary(normal)
        @test HybridPlasmaPIC._kinetic_energy_operations_certified(populations, observed, Int128(3))
        @test isequal(kinetic_energy(populations), observed)
    end
end
