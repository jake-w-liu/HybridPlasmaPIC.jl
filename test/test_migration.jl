using HybridPlasmaPIC, Test

function _ps1(xs, ids; q = 1.0, m = 1.0)
    ps = ParticleSet{1,Float64}(length(xs); q, m)
    ps.x[1] .= xs
    ps.v[1] .= 10.0 .+ ids
    ps.v[2] .= 20.0 .+ ids
    ps.v[3] .= 30.0 .+ ids
    ps.weight .= 0.5 .+ ids
    ps.id .= UInt64.(ids)
    ps.tag .= UInt32.(ids .+ 100)
    return ps
end

function _ps2(xy, ids; q = 1.0, m = 1.0)
    ps = ParticleSet{2,Float64}(length(ids); q, m)
    ps.x[1] .= first.(xy)
    ps.x[2] .= last.(xy)
    ps.v[1] .= 10.0 .+ ids
    ps.v[2] .= 20.0 .+ ids
    ps.v[3] .= 30.0 .+ ids
    ps.weight .= 0.5 .+ ids
    ps.id .= UInt64.(ids)
    ps.tag .= UInt32.(ids .+ 100)
    return ps
end

@testset "logical rank layout" begin
    g = FourierGrid((8, 6), (8.0, 6.0))
    layout = LogicalRankLayout((2, 3); periodic = (true, false))

    @test nranks(layout) == 6
    @test rank_coords(layout, 1) == (1, 1)
    @test rank_coords(layout, 4) == (2, 2)
    @test rank_index(layout, (3, 1)) == 1       # x wraps
    @test rank_index(layout, (1, 4)) === nothing # y does not

    b = rank_bounds(g, layout, 4)
    @test b.lo == (4.0, 2.0)
    @test b.hi == (8.0, 4.0)

    @test rank_of_position((0.1, 0.1), g, layout) == 1
    @test rank_of_position((4.1, 0.1), g, layout) == 2
    @test rank_of_position((8.2, 0.1), g, layout) == 1
    @test rank_of_position((0.1, 6.1), g, layout) === nothing

    @test_throws ArgumentError LogicalRankLayout((2, 0))
    @test_throws ArgumentError LogicalRankLayout((2, 2); periodic = (true,))
    @test_throws ArgumentError rank_coords(layout, 0)
    @test_throws ArgumentError rank_of_position((NaN, 0.1), g, layout)
    @test_throws ArgumentError rank_of_position((0.1, Inf), g, layout)

    huge = big(typemax(Int)) + 1
    @test_throws ArgumentError LogicalRankLayout((huge,))
    @test_throws ArgumentError LogicalRankLayout((typemax(Int), 2))
    @test_throws ArgumentError rank_coords(LogicalRankLayout((2,)), huge)
    @test_throws ArgumentError rank_index(LogicalRankLayout((2,)), (huge,))
end

@testset "multi-dimensional migration and transverse wrapping" begin
    g = FourierGrid((10, 12), (10.0, 12.0))
    layout = LogicalRankLayout((2, 2); periodic = (true, true))
    ranks = [
        _ps2([(5.1, 2.0), (-0.2, 3.0), (1.0, 12.2), (1.0, -0.1)], [11, 12, 13, 14]),
        _ps2([(4.9, 2.0), (10.1, 3.0), (7.0, 12.1), (7.0, -0.2)], [21, 22, 23, 24]),
        _ps2([(5.2, 8.0), (-0.1, 8.0), (1.0, 5.9), (1.0, 12.1)], [31, 32, 33, 34]),
        _ps2([(4.8, 8.0), (10.2, 8.0), (7.0, 5.8), (7.0, 12.2)], [41, 42, 43, 44]),
    ]

    result = migrate_particles!(ranks, g, layout)
    @test result == (moved = 14, lost = 0)
    @test sort(vcat((collect(r.id) for r in ranks)...)) ==
          UInt64[11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44]

    expected = Dict(
        1 => UInt64[13, 21, 22, 33, 34],
        2 => UInt64[11, 12, 23, 43, 44],
        3 => UInt64[14, 41, 42],
        4 => UInt64[24, 31, 32],
    )
    for r = 1:nranks(layout)
        @test sort(ranks[r].id) == expected[r]
        for p = 1:nparticles(ranks[r])
            pos = (ranks[r].x[1][p], ranks[r].x[2][p])
            @test rank_of_position(pos, g, layout) == r
            @test 0.0 <= pos[1] < 10.0
            @test 0.0 <= pos[2] < 12.0
        end
    end
end

@testset "nonperiodic transverse migration drops escaped particles" begin
    g = FourierGrid((10, 12), (10.0, 12.0))
    layout = LogicalRankLayout((1, 2); periodic = (true, false))
    ranks = [
        _ps2([(10.2, 2.0), (1.0, -0.1), (1.0, 6.2)], [1, 2, 3]),
        _ps2([(-0.2, 8.0), (1.0, 5.8), (1.0, 12.0)], [4, 5, 6]),
    ]

    result = migrate_particles!(ranks, g, layout)
    @test result == (moved = 2, lost = 2)
    @test sort(ranks[1].id) == UInt64[1, 5]
    @test sort(ranks[2].id) == UInt64[3, 4]
    for r = 1:nranks(layout), p = 1:nparticles(ranks[r])
        pos = (ranks[r].x[1][p], ranks[r].x[2][p])
        @test rank_of_position(pos, g, layout) == r
        @test 0.0 <= pos[1] < 10.0
        @test 0.0 <= pos[2] < 12.0
    end
end

@testset "particle migration across logical ranks" begin
    g = FourierGrid((10,), (10.0,))
    layout = LogicalRankLayout((4,); periodic = (true,))
    ranks = [_ps1([2.6, -0.2], [11, 12]), _ps1([5.1], [21]), _ps1([7.8], [31]), _ps1([10.2], [41])]

    result = migrate_particles!(ranks, g, layout)
    @test result.moved == 5
    @test result.lost == 0
    @test sort(vcat((collect(r.id) for r in ranks)...)) == UInt64[11, 12, 21, 31, 41]

    @test sort(ranks[1].id) == UInt64[41]
    @test sort(ranks[2].id) == UInt64[11]
    @test sort(ranks[3].id) == UInt64[21]
    @test sort(ranks[4].id) == UInt64[12, 31]

    for r = 1:nranks(layout), p = 1:nparticles(ranks[r])
        @test rank_of_position((ranks[r].x[1][p],), g, layout) == r
        @test 0.0 <= ranks[r].x[1][p] < 10.0
    end

    byid = Dict{UInt64,Tuple{Float64,Float64,Float64,Float64,UInt32}}()
    for ps in ranks, p = 1:nparticles(ps)
        byid[ps.id[p]] = (ps.v[1][p], ps.v[2][p], ps.v[3][p], ps.weight[p], ps.tag[p])
    end
    @test byid[0x000000000000000b] == (21.0, 31.0, 41.0, 11.5, 111)
    @test byid[0x000000000000000c] == (22.0, 32.0, 42.0, 12.5, 112)
    @test byid[0x0000000000000029] == (51.0, 61.0, 71.0, 41.5, 141)
end

@testset "nonperiodic migration drops escaped particles" begin
    g = FourierGrid((8,), (8.0,))
    layout = LogicalRankLayout((2,); periodic = (false,))
    ranks = [_ps1([-0.1, 1.0], [1, 2]), _ps1([7.9, 8.0], [3, 4])]

    result = migrate_particles!(ranks, g, layout)
    @test result.moved == 0
    @test result.lost == 2
    @test collect(ranks[1].id) == UInt64[2]
    @test collect(ranks[2].id) == UInt64[3]
end

@testset "migration rejects non-finite particle positions" begin
    g = FourierGrid((8,), (8.0,))
    periodic = LogicalRankLayout((2,); periodic = (true,))
    nonperiodic = LogicalRankLayout((2,); periodic = (false,))

    @test_throws ArgumentError migrate_particles!(
        [_ps1([NaN], [1]), ParticleSet{1,Float64}(0)],
        g,
        periodic,
    )
    @test_throws ArgumentError migrate_particles!(
        [_ps1([Inf], [1]), ParticleSet{1,Float64}(0)],
        g,
        nonperiodic,
    )
end

@testset "migration validates before mutation" begin
    g = FourierGrid((8,), (8.0,))
    periodic = LogicalRankLayout((2,); periodic = (true,))

    species = [_ps1([5.0], [1]; q = 1.0), _ps1(Float64[], Int[]; q = 2.0)]
    @test_throws ArgumentError migrate_particles!(species, g, periodic)
    @test nparticles(species[1]) == 1
    @test species[1].x[1] == [5.0]
    @test species[1].id == UInt64[1]
    @test nparticles(species[2]) == 0

    nonfinite = [_ps1([-0.25, NaN], [1, 2]), ParticleSet{1,Float64}(0)]
    @test_throws ArgumentError migrate_particles!(nonfinite, g, periodic)
    @test nonfinite[1].x[1][1] == -0.25
    @test isnan(nonfinite[1].x[1][2])
    @test nonfinite[1].id == UInt64[1, 2]
    @test nparticles(nonfinite[2]) == 0
end

@testset "append_particles! rejects species mismatch" begin
    dest = _ps1([0.1], [1]; q = 1.0, m = 1.0)
    src_q = _ps1([0.2], [2]; q = 2.0, m = 1.0)
    src_m = _ps1([0.3], [3]; q = 1.0, m = 2.0)
    @test_throws ArgumentError append_particles!(dest, src_q)
    @test_throws ArgumentError append_particles!(dest, src_m)
end

@testset "slab field halo exchange" begin
    layout = LogicalRankLayout((3,); periodic = (true,))
    a = [
        [100.0, 10.0, 11.0, 12.0, 101.0],
        [200.0, 20.0, 21.0, 22.0, 201.0],
        [300.0, 30.0, 31.0, 32.0, 301.0],
    ]

    stats = exchange_field_halos!(a, layout; halo = 1)
    @test stats == (exchanged = 6, filled = 0)
    @test a[1] == [32.0, 10.0, 11.0, 12.0, 20.0]
    @test a[2] == [12.0, 20.0, 21.0, 22.0, 30.0]
    @test a[3] == [22.0, 30.0, 31.0, 32.0, 10.0]

    nonperiodic = LogicalRankLayout((2,); periodic = (false,))
    b = [[-1.0, 10.0, 11.0, 12.0, -1.0], [-2.0, 20.0, 21.0, 22.0, -2.0]]
    stats = exchange_field_halos!(b, nonperiodic; halo = 1, fill_value = -99.0)
    @test stats == (exchanged = 2, filled = 2)
    @test b[1] == [-99.0, 10.0, 11.0, 12.0, 20.0]
    @test b[2] == [12.0, 20.0, 21.0, 22.0, -99.0]

    huge = big(typemax(Int)) + 1
    @test_throws ArgumentError exchange_field_halos!(
        [copy(a[1])],
        LogicalRankLayout((1,));
        halo = huge,
    )
end

@testset "tuple field halo exchange and validation" begin
    layout = LogicalRankLayout((2,); periodic = (false,))
    fields = [
        ([1.0, 10.0, 11.0, 12.0, 2.0], [3.0, 20.0, 21.0, 22.0, 4.0]),
        ([5.0, 30.0, 31.0, 32.0, 6.0], [7.0, 40.0, 41.0, 42.0, 8.0]),
    ]

    stats = exchange_field_halos!(fields, layout; halo = 1, fill_value = -1.0)
    @test stats == (exchanged = 4, filled = 4)
    @test fields[1][1] == [-1.0, 10.0, 11.0, 12.0, 30.0]
    @test fields[1][2] == [-1.0, 20.0, 21.0, 22.0, 40.0]
    @test fields[2][1] == [12.0, 30.0, 31.0, 32.0, -1.0]
    @test fields[2][2] == [22.0, 40.0, 41.0, 42.0, -1.0]

    pencil = LogicalRankLayout((2, 2))
    mats = [zeros(4, 4) for _ = 1:4]
    @test_throws ArgumentError exchange_field_halos!(mats, pencil; halo = 1)
    @test_throws DimensionMismatch exchange_field_halos!([zeros(2), zeros(4)], layout; halo = 1)
end

@testset "slab ghost moment exchange" begin
    layout = LogicalRankLayout((3,); periodic = (true,))
    a = [
        [100.0, 10.0, 11.0, 12.0, 101.0],
        [200.0, 20.0, 21.0, 22.0, 201.0],
        [300.0, 30.0, 31.0, 32.0, 301.0],
    ]

    stats = exchange_ghost_moments!(a, layout; halo = 1)
    @test stats == (exchanged = 6, dropped = 0)
    @test a[1] == [0.0, 311.0, 11.0, 212.0, 0.0]
    @test a[2] == [0.0, 121.0, 21.0, 322.0, 0.0]
    @test a[3] == [0.0, 231.0, 31.0, 132.0, 0.0]

    nonperiodic = LogicalRankLayout((2,); periodic = (false,))
    b = [[5.0, 10.0, 11.0, 12.0, 6.0], [7.0, 20.0, 21.0, 22.0, 8.0]]
    stats = exchange_ghost_moments!(b, nonperiodic; halo = 1)
    @test stats == (exchanged = 2, dropped = 2)
    @test b[1] == [0.0, 10.0, 11.0, 19.0, 0.0]
    @test b[2] == [0.0, 26.0, 21.0, 22.0, 0.0]
end

@testset "tuple ghost moment exchange and validation" begin
    layout = LogicalRankLayout((2,); periodic = (true,))
    moments = [
        ([1.0, 10.0, 11.0, 12.0, 2.0], [3.0, 20.0, 21.0, 22.0, 4.0]),
        ([5.0, 30.0, 31.0, 32.0, 6.0], [7.0, 40.0, 41.0, 42.0, 8.0]),
    ]

    stats = exchange_ghost_moments!(moments, layout; halo = 1)
    @test stats == (exchanged = 8, dropped = 0)
    @test moments[1][1] == [0.0, 16.0, 11.0, 17.0, 0.0]
    @test moments[1][2] == [0.0, 28.0, 21.0, 29.0, 0.0]
    @test moments[2][1] == [0.0, 32.0, 31.0, 33.0, 0.0]
    @test moments[2][2] == [0.0, 44.0, 41.0, 45.0, 0.0]

    pencil = LogicalRankLayout((2, 2))
    mats = [zeros(4, 4) for _ = 1:4]
    @test_throws ArgumentError exchange_ghost_moments!(mats, pencil; halo = 1)
    @test_throws DimensionMismatch exchange_ghost_moments!([zeros(2), zeros(4)], layout; halo = 1)
end
