module PeriodicMigrationTests
using HybridPlasmaPIC, Test

@testset "Periodic rank lookup and migration use the half-open seam" begin
    for T in (Float32, Float64), D = 1:3, axis = 1:D, mask = 0:(2^D-1)
        periodic = ntuple(d -> !iszero(mask & (1 << (d - 1))), D)
        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        layout = LogicalRankLayout(ntuple(_ -> 2, D); periodic)
        for (x, slab) in (
            (-eps(T) / T(4), 1),
            (-nextfloat(zero(T)), 1),
            (-one(T), 1),
            (zero(T), 1),
            (one(T), 1),
            (nextfloat(one(T)), 1),
            (prevfloat(one(T)), 2),
            (T(1.25), 1),
            (-T(0.25), 2),
        )
            escaped = !periodic[axis] && !(zero(T) <= x < one(T))
            expected = escaped ? nothing : 1 + (slab - 1) * 2^(axis - 1)
            original = ntuple(d -> d == axis ? x : T(0.25), D)
            @test rank_of_position(original, grid, layout) == expected
            ps = ParticleSet{D,T}(1; q = T(2), m = T(3))
            for d = 1:D
                ps.x[d][1] = original[d]
            end
            for c = 1:3
                ps.v[c][1] = T(c) / T(8)
            end
            ps.weight[1] = T(1.25)
            ps.id[1] = 123
            ps.tag[1] = 456
            reference = deepcopy(ps)
            periodic[axis] && apply_periodic!(reference, ntuple(_ -> zero(T), D), grid.L)
            populations = [ParticleSet{D,T}(0; q = ps.q, m = ps.m) for _ = 1:nranks(layout)]
            populations[1] = ps
            result = migrate_particles!(populations, grid, layout)
            @test populations[1] === ps
            @test result == (moved = !escaped && expected != 1 ? 1 : 0, lost = escaped ? 1 : 0)
            @test sum(nparticles, populations) == (escaped ? 0 : 1)
            if !escaped
                destination = populations[expected]
                @test nparticles(destination) == 1
                @test destination.x == reference.x
                @test destination.v == reference.v
                @test destination.weight == reference.weight
                @test destination.id == reference.id
                @test destination.tag == reference.tag
                @test destination.q === reference.q && destination.m === reference.m
                @test all(d -> zero(T) <= destination.x[d][1] < grid.L[d], 1:D)
                @test rank_of_position(ntuple(d -> destination.x[d][1], D), grid, layout) ==
                      expected
            end
        end
    end
end

@testset "Periodic translation and nonunit seam" begin
    for T in (Float32, Float64), D = 1:3
        # Binary endpoints put the predecessor of lo halfway between hi and
        # its predecessor after translation; ties round to hi, which folds to lo.
        lo, hi = ntuple(_ -> one(T), D), ntuple(_ -> T(1.125), D)
        ps = ParticleSet{D,T}(2)
        for d = 1:D
            ps.x[d] .= (prevfloat(one(T)), T(1.0625))
        end
        apply_periodic!(ps, lo, hi)
        @test all(d -> ps.x[d] == T[1, 1.0625], 1:D)

        lengths = ntuple(d -> T(2)^d, D)
        grid = FourierGrid(ntuple(_ -> 4, D), lengths)
        layout = LogicalRankLayout(ntuple(_ -> 2, D))
        position = ntuple(d -> -eps(T) * lengths[d] / T(4), D)
        @test rank_of_position(position, grid, layout) == 1
        for d = 1:D
            ps.x[d] .= (position[d], lengths[d] / T(4))
        end
        apply_periodic!(ps, ntuple(_ -> zero(T), D), lengths)
        @test all(d -> ps.x[d] == T[0, lengths[d]/T(4)], 1:D)
    end
end
end
