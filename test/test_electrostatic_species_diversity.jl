module ElectrostaticSpeciesDiversityTests
using HybridPlasmaPIC, Test, Random

@testset "ES three active charges and unequal population sizes" begin
    for T in (Float32, Float64), D = 1:3
        rng = MersenneTwister(82469)
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D))
        counts = (19, 7, 29, 5)
        charges = (-3, 1, 2, 0)
        masses = (1, 5, 0.7, 7)
        populations = [ParticleSet{D,T}(counts[k]; q = charges[k], m = masses[k]) for k = 1:4]
        for ps in populations
            for a in ps.x
                a .= rand(rng, T, nparticles(ps))
            end
            for a in ps.v
                a .= T(0.1) .* randn(rng, T, nparticles(ps))
            end
            ps.weight .= one(T) / T(nparticles(ps))
        end
        original = deepcopy(populations)
        reversed = reverse(deepcopy(populations))
        state = ElectrostaticPIC(g, populations)
        reverse_state = ElectrostaticPIC(FourierGrid(g.n, g.L), reversed)
        init_espic!(state, populations)
        init_espic!(reverse_state, reversed)
        expected, scratch = zeros(T, g.n), zeros(T, g.n)
        for ps in populations
            density!(scratch, ps, g, CIC())
            expected .-= ps.q .* scratch
        end
        @test state.ne == expected
        @test all(
            isapprox(a, b; atol = 64eps(T), rtol = 64eps(T)) for
            (a, b) in zip(state.E, reverse_state.E)
        )
        # An imposed uniform old E isolates each q/m kick from the new Poisson
        # solve and detects any field update between population pushes.
        force = ntuple(c -> c <= D ? T(0.03c) : zero(T), 3)
        for c = 1:3
            fill!(state.E[c], force[c])
            fill!(reverse_state.E[c], force[c])
        end
        dt = T(0.02)
        step_espic!(state, populations, dt)
        step_espic!(reverse_state, reversed, dt)
        for (ps, before, permuted) in zip(populations, original, reverse(reversed))
            @test ps.x == permuted.x
            @test ps.v == permuted.v
            for c = 1:3
                expected_v = before.v[c] .+ (dt / 2 * ps.q / ps.m) * force[c]
                @test isapprox(ps.v[c], expected_v; atol = 32eps(T), rtol = 32eps(T))
            end
            for d = 1:D
                @test ps.x[d] ≈ mod.(before.x[d] .+ dt .* ps.v[d], one(T))
            end
        end
        @test all(
            isapprox(a, b; atol = 128eps(T), rtol = 128eps(T)) for
            (a, b) in zip(state.E, reverse_state.E)
        )
        diagnostic = electrostatic_diagnostics(state, populations)
        @test diagnostic.field == field_energy(state)
        @test diagnostic.kinetic == kinetic_energy(populations)
        @test diagnostic.momentum == total_momentum(populations)
        # Staged logical lengths match current population counts; the model
        # carries one shared field, not one mesh field per population.
        @test length(state.collection.staged) == length(populations)
        @test map(nparticles, state.collection.staged) == collect(counts)
        @test all(size(a) == g.n for a in state.collection.E)
    end
end
@testset "ES particle growth shares the existing half level and workspaces" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        species = [ParticleSet{D,T}(2; q = q) for q in (-1, 1)]
        for ps in species
            foreach(a -> (a .= T[0.2, 0.6]), ps.x)
            foreach(a -> fill!(a, T(0.02)), ps.v)
        end
        es = ElectrostaticPIC(g, species)
        init_espic!(es, species)
        step_espic!(es, species, T(0.01))
        mesh = es.collection.E
        stages = copy(es.collection.staged)
        stage_x = stages[1].x[1]
        for ps in species
            for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
                resize!(a, 7)
            end
            foreach(a -> (a[3:7] .= T(0.4)), ps.x)
            # New velocities are supplied at the current shared half level.
            foreach(a -> (a[3:7] .= T(-0.03)), ps.v)
            ps.weight[3:7] .= one(T)
            ps.id[3:7] .= UInt64.(3:7)
            ps.tag[3:7] .= UInt32(9)
        end
        positions = deepcopy(species[1].x)
        velocities = deepcopy(species[1].v)
        step_espic!(es, species, T(0.04))
        @test es.last_dt[] == T(0.04)
        @test all(ps.v == velocities for ps in species)
        for ps in species, d = 1:D
            @test ps.x[d] == positions[d] .+ T(0.04) .* velocities[d]
        end
        @test all(es.collection.staged[k] === stages[k] for k = 1:2)
        @test es.collection.staged[1].x[1] === stage_x
        @test all(a === b for (a, b) in zip(es.collection.E, mesh))
        @test all(nparticles(ps) == 7 for ps in es.collection.staged)
        @test all(ps.id[3:7] == UInt64.(3:7) && all(==(UInt32(9)), ps.tag[3:7]) for ps in species)
        empty = species[1:0]
        init_espic!(es, empty)
        step_espic!(es, empty, T(0.01))
        @test isempty(es.collection.staged)
        @test iszero(field_energy(es))
    end
end

end
