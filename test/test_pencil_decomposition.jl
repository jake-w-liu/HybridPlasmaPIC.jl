using HybridPlasmaPIC
using Test

@testset "fully periodic 3D pencil decomposition" begin
    dec = PencilDecomposition3D((17, 19, 23), (3, 4))
    @test pencil_nranks(dec) == 12
    @test pencil_rank_coords(dec, 1) == (1, 1)
    @test pencil_rank_coords(dec, 5) == (2, 2)
    @test pencil_rank_index(dec, (3, 4)) == 12
    @test_throws ArgumentError pencil_rank_coords(dec, 0)
    @test_throws ArgumentError pencil_rank_index(dec, (4, 1))
    @test_throws ArgumentError pencil_rank_index(dec, (1,))
    @test_throws ArgumentError pencil_owner(dec, (1, 1))
    @test_throws ArgumentError pencil_bounds(dec, 1, :bad)

    for orientation in (:x, :y, :z)
        axes = pencil_orientation_axes(orientation)
        seen = falses(dec.n)
        total = 0

        for rank = 1:pencil_nranks(dec)
            ranges = pencil_bounds(dec, rank, orientation)
            @test ranges[axes.full_axis] == 1:dec.n[axes.full_axis]
            @test pencil_local_size(dec, rank, orientation) == map(length, ranges)

            for I in CartesianIndices(ranges)
                idx = Tuple(I)
                @test !seen[idx...]
                seen[idx...] = true
                total += 1
                @test pencil_owner(dec, idx, orientation) == rank
            end
        end

        @test total == prod(dec.n)
        @test all(seen)
        @test pencil_owner(dec, (dec.n[1] + 1, 1, 1), orientation) ==
              pencil_owner(dec, (1, 1, 1), orientation)
        @test pencil_owner(dec, (1, dec.n[2] + 2, 1), orientation) ==
              pencil_owner(dec, (1, 2, 1), orientation)
        @test pencil_owner(dec, (1, 1, -1), orientation) ==
              pencil_owner(dec, (1, 1, dec.n[3] - 1), orientation)
    end
end

@testset "pencil decomposition validation" begin
    huge = big(typemax(Int)) + 1
    @test_throws ArgumentError PencilDecomposition3D((8, 8), (2, 2))
    @test_throws ArgumentError PencilDecomposition3D((8, 8, 8), (0, 2))
    @test_throws ArgumentError PencilDecomposition3D((8, 8, 8), (2, 0))
    @test_throws ArgumentError PencilDecomposition3D((8, 8, 8), (9, 2))
    @test_throws ArgumentError PencilDecomposition3D((8, 8, 8), (2, 9))
    @test_throws ArgumentError PencilDecomposition3D((huge, 8, 8), (2, 2))
    @test_throws ArgumentError PencilDecomposition3D((8, 8, 8), (huge, 2))
    @test_throws ArgumentError pencil_rank_coords(PencilDecomposition3D((8, 8, 8), (2, 2)), huge)
    @test_throws ArgumentError pencil_rank_index(
        PencilDecomposition3D((8, 8, 8), (2, 2)),
        (huge, 1),
    )
    @test_throws ArgumentError pencil_owner(PencilDecomposition3D((8, 8, 8), (2, 2)), (huge, 1, 1))
    @test_throws ArgumentError PencilDecomposition3D(
        (8, 8, 8),
        (2, 2);
        periodic = (true, false, true),
    )
end
