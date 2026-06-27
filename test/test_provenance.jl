using HybridPlasmaPIC, Random, Test

function _sorted_ids(rank_particles)
    ids = UInt64[]
    for ps in rank_particles
        append!(ids, ps.id)
    end
    return sort(ids)
end

@testset "rank-local RNG streams" begin
    layout = LogicalRankLayout((2, 3); periodic = (true, false))

    seed = rank_seed(1234, layout, 4)
    @test seed == rank_seed(1234, layout, 4)
    @test seed != rank_seed(1234, layout, 5)
    @test seed != rank_seed(1234, layout, 4; stream = 1)

    @test rand(rank_rng(1234, layout, 4), 8) == rand(rank_rng(1234, layout, 4), 8)
    @test rand(rank_rng(1234, layout, 4; stream = 2), 8) ==
          rand(MersenneTwister(rank_seed(1234, layout, 4; stream = 2)), 8)

    @test_throws ArgumentError rank_seed(-1, layout, 1)
    @test_throws ArgumentError rank_seed(1234, layout, 0)
    @test_throws ArgumentError rank_seed(1234, layout, 1; stream = -1)
end

@testset "rank-independent particle IDs" begin
    @test global_particle_id(1) == 0x0000000000000001
    @test global_particle_id(42; species = 3) == (UInt64(3) << 48) | UInt64(42)
    @test_throws ArgumentError global_particle_id(0)
    @test_throws ArgumentError global_particle_id(UInt64(1) << 48)
    @test_throws ArgumentError global_particle_id(1; species = UInt64(1) << 16)

    serial = ParticleSet{2,Float64}(8)
    assign_global_particle_ids!(serial, 1; species = 7)

    left = ParticleSet{2,Float64}(3)
    right = ParticleSet{2,Float64}(5)
    assign_global_particle_ids!(left, 1; species = 7)
    assign_global_particle_ids!(right, 4; species = 7)
    @test vcat(left.id, right.id) == serial.id
    @test length(unique(serial.id)) == nparticles(serial)

    subset = ParticleSet{1,Float64}(3)
    assign_global_particle_ids!(subset, [2, 5, 8]; species = 1)
    @test subset.id == UInt64[
        global_particle_id(2; species = 1),
        global_particle_id(5; species = 1),
        global_particle_id(8; species = 1),
    ]
    @test_throws DimensionMismatch assign_global_particle_ids!(subset, [1, 2])
    @test_throws ArgumentError assign_global_particle_ids!(subset, [1, 0, 2])
end

@testset "particle IDs survive migration" begin
    g = FourierGrid((10,), (10.0,))
    layout = LogicalRankLayout((2,); periodic = (true,))
    ranks = [ParticleSet{1,Float64}(2), ParticleSet{1,Float64}(1)]

    ranks[1].x[1] .= [5.2, -0.2]
    ranks[2].x[1] .= [0.1]
    assign_global_particle_ids!(ranks[1], [1, 2]; species = 2)
    assign_global_particle_ids!(ranks[2], [3]; species = 2)

    before = _sorted_ids(ranks)
    stats = migrate_particles!(ranks, g, layout)
    @test stats == (moved = 3, lost = 0)
    @test _sorted_ids(ranks) == before
    @test collect(ranks[1].id) == UInt64[global_particle_id(3; species = 2)]
    @test sort(ranks[2].id) ==
          UInt64[global_particle_id(1; species = 2), global_particle_id(2; species = 2)]
end
