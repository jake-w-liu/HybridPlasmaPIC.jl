module RankClassificationTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
function rank_oracle(layout, coords)
    result = big(1)
    stride = big(1)
    for d in eachindex(coords)
        c = BigInt(coords[d])
        n = BigInt(layout.ranks[d])
        if layout.periodic[d]
            c = mod(c - 1, n) + 1
        elseif !(1 <= c <= n)
            return nothing
        end
        result += (c - 1) * stride
        stride *= n
    end
    return Int(result)
end
function captured_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end
function check_range_error(f, label)
    err = captured_error(f)
    @test err isa ArgumentError
    @test err isa ArgumentError && err.msg == "$label must fit in Int"
end
function classification_bytes(layout, coords, g, x)
    rank_index(layout, coords)
    rank_of_position(x, g, layout)
    return (@allocated(rank_index(layout, coords)), @allocated(rank_of_position(x, g, layout)))
end
function pencil_bytes(dec, orientation)
    pencil_rank_index(dec, (1, 1))
    pencil_owner(dec, (1, 1, 1), orientation)
    return (
        @allocated(pencil_rank_index(dec, (1, 1))),
        @allocated(pencil_owner(dec, (1, 1, 1), orientation))
    )
end
@testset "Rank lookup retains integer and periodic boundary semantics" begin
    for D = 1:3, mask = 0:(2^D-1)
        periodic = ntuple(d -> !iszero(mask & (1 << (d - 1))), D)
        layout = LogicalRankLayout(ntuple(d -> d + 1, D); periodic)
        for dimension = 1:D, value in (typemin(Int), -1, 0, 1, 2, typemax(Int))
            coords = ntuple(d -> d == dimension ? value : 1, D)
            expected = rank_oracle(layout, coords)
            @test rank_index(layout, coords) == expected
            @test rank_index(layout, map(BigInt, coords)) == expected
        end
        for T in (Float32, Float64), dimension = 1:D, numerator in (-5, -1, 0, 1, 2, 3, 4, 5, 9)
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
            x = ntuple(d -> T(d == dimension ? numerator : 1) / T(4), D)
            rational = ntuple(d -> (d == dimension ? numerator : 1) // big(4), D)
            escaped = any(!periodic[d] && !(0 <= rational[d] < 1) for d = 1:D)
            expected = if escaped
                nothing
            else
                coords = ntuple(d -> Int(floor(mod(rational[d], 1) * layout.ranks[d])) + 1, D)
                rank_oracle(layout, coords)
            end
            @test rank_of_position(x, g, layout) == expected
        end
        for dimension = 1:D,
            bad in (big(typemax(Int)) + 1, big(typemin(Int)) - 1, BigInt(typemax(UInt)))
            # Even an earlier nonperiodic coordinate outside the domain must
            # not bypass conversion validation of a later coordinate.
            coords = ntuple(d -> d == dimension ? bad : big(0), D)
            check_range_error(() -> rank_index(layout, coords), "coords[$dimension]")
        end
    end
end
@testset "Indexed labels retain pencil and MPI-neighbor validation text" begin
    dec = PencilDecomposition3D((8, 10, 12), (2, 2))
    huge = big(typemax(Int)) + 1
    for k = 1:2
        coords = ntuple(d -> d == k ? huge : big(1), 2)
        check_range_error(() -> pencil_rank_index(dec, coords), "coords[$k]")
        check_range_error(
            () -> invoke(pencil_rank_index, Tuple{typeof(dec),Tuple}, dec, coords),
            "coords[$k]",
        )
    end
    for k = 1:3, orientation in (:x, :y, :z)
        index = ntuple(d -> d == k ? huge : big(1), 3)
        check_range_error(() -> pencil_owner(dec, index, orientation), "index[$k]")
        check_range_error(
            () -> invoke(pencil_owner, Tuple{typeof(dec),Tuple,Any}, dec, index, orientation),
            "index[$k]",
        )
    end
    for label in ("offset", "neighbor coordinate", "index")
        @test H._rank_coordinate(label, typemin(Int)) == typemin(Int)
        @test H._rank_coordinate(label, typemax(Int)) == typemax(Int)
        check_range_error(() -> H._rank_coordinate(label, huge), label)
    end
    @test invoke(pencil_rank_index, Tuple{typeof(dec),Tuple}, dec, (Int8(1), UInt16(1))) == 1
    for orientation in (:x, :y, :z)
        @test invoke(
            pencil_owner,
            Tuple{typeof(dec),Tuple,Any},
            dec,
            (Int8(1), UInt16(1), big(1)),
            orientation,
        ) == 1
    end
end
@testset "Ordinary classification does not allocate indexed error labels" begin
    for T in (Float32, Float64), D = 1:3, periodic in (false, true)
        layout = LogicalRankLayout(ntuple(_ -> 2, D); periodic = ntuple(_ -> periodic, D))
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        coords = ntuple(_ -> 1, D)
        x = ntuple(_ -> T(0.25), D)
        classification_bytes(layout, coords, g, x)
        @test classification_bytes(layout, coords, g, x) == (0, 0)
    end
    dec = PencilDecomposition3D((8, 10, 12), (2, 2))
    for orientation in (:x, :y, :z)
        pencil_bytes(dec, orientation)
        @test pencil_bytes(dec, orientation) == (0, 0)
    end
end

@testset "Pencil orientation normalization preserves values and rejection order" begin
    dec = PencilDecomposition3D((8, 10, 12), (2, 2))
    for orientation in (:x, :y, :z)
        value = Val(orientation)
        @test pencil_orientation_axes(orientation) == pencil_orientation_axes(value)
        for rank = 1:pencil_nranks(dec)
            @test pencil_bounds(dec, rank, orientation) == pencil_bounds(dec, rank, value)
            @test pencil_local_size(dec, rank, orientation) == pencil_local_size(dec, rank, value)
        end
        for index in ((1, 1, 1), (0, -1, 13), (8, 10, 12))
            @test pencil_owner(dec, index, orientation) == pencil_owner(dec, index, value)
        end
    end
    for orientation in (:bad, Val(:bad))
        for f in (
            () -> pencil_orientation_axes(orientation),
            () -> pencil_bounds(dec, 0, orientation),
            () -> pencil_owner(dec, (1, 1, 1), orientation),
        )
            err = captured_error(f)
            @test err isa ArgumentError
            @test err isa ArgumentError &&
                  err.msg == "pencil orientation must be :x, :y, or :z, got bad"
        end
        huge = big(typemax(Int)) + 1
        check_range_error(() -> pencil_owner(dec, (huge, big(1), big(1)), orientation), "index[1]")
    end
end
end
