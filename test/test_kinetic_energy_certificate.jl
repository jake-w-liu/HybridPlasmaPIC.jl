using HybridPlasmaPIC, Test, Random

@testset "Kinetic forward-error enclosures vs exact expressions" begin
    Q = Rational{BigInt}
    rng = MersenneTwister(182619)
    for T in (Float16, Float32, Float64), count in (1, 3, 16, 129), trial = 1:16
        bounds = HybridPlasmaPIC._kinetic_energy_bounds(T, Int128(count))
        lo, hi = bounds[1:2]
        pool = T[lo, hi, prevfloat(hi), one(T), clamp(T(0.5) + rand(rng, T), lo, hi)]
        filter!(x -> lo <= x <= hi, pool)
        particles =
            [ParticleSet{1,T}(count ÷ 2), ParticleSet{1,T}(count - count ÷ 2), ParticleSet{1,T}(0)]
        for (j, ps) in enumerate(particles)
            ps.m = rand(rng, pool) * (isodd(j) ? one(T) : -one(T))
            for i in eachindex(ps.weight)
                ps.weight[i] = rand(rng, pool) * (rand(rng, Bool) ? one(T) : -one(T))
                for v in ps.v
                    v[i] = rand(rng, pool)
                end
            end
        end
        observed = zero(T)
        exact = Q(0)
        absolute_exact = Q(0)
        for ps in particles
            observed += KineticRangeOracle.ordinary(ps)
            for i in eachindex(ps.weight), v in ps.v
                term = Q(ps.m) * Q(ps.weight[i]) * Q(v[i])^2 / 2
                exact += term
                absolute_exact += abs(term)
            end
        end
        absolute = HybridPlasmaPIC._kinetic_energy_absolute(particles, T)
        error = HybridPlasmaPIC._kinetic_energy_error_bound(absolute, Int128(count))
        @test isfinite(absolute) && absolute >= 0
        @test abs(Q(observed) - exact) <= Q(error)
        # Verify the absolute computation itself against an independent exact
        # sum, rather than assume it bounds the signed total without rounding.
        u = Q(1) / (BigInt(1) << precision(T))
        k = count + 4
        gamma = k * u / (1 - k * u)
        @test abs(Q(absolute) - absolute_exact) <= gamma * absolute_exact
        certified =
            HybridPlasmaPIC._kinetic_energy_normal_certified(particles, observed, Int128(count))
        if certified
            @test abs(exact) >= Q(floatmin(T))
            @test abs(exact) < Q(floatmax(T))
            @test isequal(kinetic_energy(particles), observed)
        else
            probes =
                [HybridPlasmaPIC._kinetic_energy_probe(ps, bounds, Val(true)) for ps in particles]
            existing_certificate =
                all(r -> r[3], probes) ||
                (all(r -> r[2], probes) && !(any(r -> r[4], probes) && any(r -> r[5], probes)))
            if existing_certificate
                @test isequal(kinetic_energy(particles), observed)
            else
                @test KineticRangeOracle.rounded(kinetic_energy(particles), exact)
            end
        end
    end
end

@testset "Kinetic enclosure arithmetic and count thresholds" begin
    Q = Rational{BigInt}
    for T in (Float16, Float32, Float64)
        limit = Int128(1) << (precision(T) - 2)
        u = Q(1) / (BigInt(1) << precision(T))
        for count in (Int128(1), limit - 5, limit - 4, limit - 3, limit)
            bound = HybridPlasmaPIC._kinetic_energy_error_bound(one(T), count)
            if count + 4 <= limit
                delta = ldexp(T(count + 4), -precision(T))
                @test Q(delta) == (count + 4) * u
                @test Q(one(T) - T(2) * delta) == 1 - 2 * (count + 4) * u
                @test Q(bound) >= (count + 4) * u / (1 - 2 * (count + 4) * u)
            else
                @test bound == T(Inf)
            end
        end
    end
end
