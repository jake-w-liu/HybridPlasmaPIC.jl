# pencil_decomposition.jl — 3-D pencil decomposition planner for fully periodic domains.

"""
    PencilDecomposition3D(n, ranks; periodic=(true, true, true))

Logical 3-D pencil decomposition for fully periodic grids. `n` is the global
cell count and `ranks=(p, q)` is the two-dimensional process grid. For an
`x`-pencil, each rank owns all `x` cells and a partition of `y,z`; for `y`- and
`z`-pencils the full axis changes accordingly.

This is a transport-independent planner. It provides exact local index ranges
and owner lookups for future distributed FFT and MPI transpose layers, but does
not itself perform communication.
"""
struct PencilDecomposition3D
    n::NTuple{3,Int}
    ranks::NTuple{2,Int}
    periodic::NTuple{3,Bool}

    function PencilDecomposition3D(n::NTuple{3,Int}, ranks::NTuple{2,Int}, periodic::NTuple{3,Bool})
        return new(n, ranks, periodic)
    end
end

function PencilDecomposition3D(
    n::NTuple{3,<:Integer},
    ranks::NTuple{2,<:Integer};
    periodic::NTuple{3,Bool} = (true, true, true),
)
    nn = ntuple(d -> _require_positive_intlike("n[$d]", n[d]), 3)
    rr = ntuple(d -> _require_positive_intlike("ranks[$d]", ranks[d]), 2)
    return _pencil_decomposition3d_checked(nn, rr, periodic)
end

function PencilDecomposition3D(
    n::Tuple{A,B,C},
    ranks::Tuple{R,S};
    periodic::NTuple{3,Bool} = (true, true, true),
) where {A,B,C,R,S}
    all(x -> x isa Integer, n) ||
        throw(ArgumentError("PencilDecomposition3D grid sizes must be integers, got $n"))
    all(x -> x isa Integer, ranks) ||
        throw(ArgumentError("PencilDecomposition3D rank dimensions must be integers, got $ranks"))
    nn = ntuple(d -> _require_positive_intlike("n[$d]", n[d]), 3)
    rr = ntuple(d -> _require_positive_intlike("ranks[$d]", ranks[d]), 2)
    return _pencil_decomposition3d_checked(nn, rr, periodic)
end

function PencilDecomposition3D(
    n::Tuple,
    ranks::Tuple;
    periodic::NTuple{3,Bool} = (true, true, true),
)
    length(n) == 3 || throw(
        ArgumentError("PencilDecomposition3D requires exactly three grid sizes, got $(length(n))"),
    )
    length(ranks) == 2 || throw(
        ArgumentError(
            "PencilDecomposition3D requires exactly two pencil rank dimensions, got $(length(ranks))",
        ),
    )
    throw(
        ArgumentError(
            "PencilDecomposition3D expects three integer grid sizes and two integer ranks",
        ),
    )
end

function _pencil_decomposition3d_checked(nn::NTuple{3,Int}, rr::NTuple{2,Int}, pp::NTuple{3,Bool})
    all(>(0), nn) || throw(ArgumentError("global pencil grid sizes must be positive, got $nn"))
    all(>(0), rr) || throw(ArgumentError("pencil rank grid sizes must be positive, got $rr"))
    all(pp) || throw(
        ArgumentError(
            "PencilDecomposition3D requires a fully periodic 3-D domain, got periodic=$pp",
        ),
    )
    rr[1] <= min(nn[1], nn[2]) || throw(
        ArgumentError(
            "first pencil rank dimension $(rr[1]) exceeds available x/y cells $(nn[1:2])",
        ),
    )
    rr[2] <= min(nn[2], nn[3]) || throw(
        ArgumentError(
            "second pencil rank dimension $(rr[2]) exceeds available y/z cells $(nn[2:3])",
        ),
    )
    return PencilDecomposition3D(nn, rr, pp)
end

pencil_decomposition(n, ranks; kwargs...) = PencilDecomposition3D(n, ranks; kwargs...)
pencil_nranks(dec::PencilDecomposition3D) = prod(dec.ranks)

function _check_pencil_rank(dec::PencilDecomposition3D, rank::Integer)
    1 <= rank <= pencil_nranks(dec) ||
        throw(ArgumentError("pencil rank must be in 1:$(pencil_nranks(dec)), got $rank"))
    return Int(rank)
end

"""
    pencil_rank_coords(dec, rank) -> NTuple{2,Int}

Return the 1-based coordinates of `rank` in the pencil process grid.
"""
function pencil_rank_coords(dec::PencilDecomposition3D, rank::Integer)
    r = _check_pencil_rank(dec, rank)
    return (mod(r - 1, dec.ranks[1]) + 1, fld(r - 1, dec.ranks[1]) + 1)
end

"""
    pencil_rank_index(dec, coords) -> Int

Return the 1-based linear rank for 1-based pencil process coordinates.
"""
function pencil_rank_index(dec::PencilDecomposition3D, coords::NTuple{2,<:Integer})
    cc = ntuple(d -> _rank_coordinate("coords[$d]", coords[d]), 2)
    return _pencil_rank_index_checked(dec, cc)
end

function pencil_rank_index(dec::PencilDecomposition3D, coords::Tuple)
    length(coords) == 2 ||
        throw(ArgumentError("pencil rank coordinates must have length 2, got $(length(coords))"))
    all(x -> x isa Integer, coords) ||
        throw(ArgumentError("pencil rank coordinates must be integers, got $coords"))
    cc = ntuple(d -> _rank_coordinate("coords[$d]", coords[d]), 2)
    return _pencil_rank_index_checked(dec, cc)
end

function _pencil_rank_index_checked(dec::PencilDecomposition3D, cc::NTuple{2,Int})
    all(d -> 1 <= cc[d] <= dec.ranks[d], 1:2) ||
        throw(ArgumentError("pencil rank coordinates $cc are outside 1:$(dec.ranks)"))
    return (cc[2] - 1) * dec.ranks[1] + cc[1]
end

function _pencil_axes(::Val{:x})
    return 1, (2, 3)
end

function _pencil_axes(::Val{:y})
    return 2, (1, 3)
end

function _pencil_axes(::Val{:z})
    return 3, (1, 2)
end

function _pencil_axes(::Val{name}) where {name}
    throw(ArgumentError("pencil orientation must be :x, :y, or :z, got $(name)"))
end

_pencil_val(orientation::Symbol) = Val(orientation)
_pencil_val(orientation::Val) = orientation

function _partition_range(n::Int, parts::Int, coord::Int)
    1 <= coord <= parts || throw(ArgumentError("partition coordinate $coord is outside 1:$parts"))
    lo = fld((coord - 1) * n, parts) + 1
    hi = fld(coord * n, parts)
    return lo:hi
end

function _partition_owner(index::Integer, n::Int, parts::Int)
    i = mod1(_rank_coordinate("index", index), n)
    return cld(i * parts, n)
end

"""
    pencil_bounds(dec, rank, orientation=:x) -> NTuple{3,UnitRange{Int}}

Return global 1-based half-open-equivalent owned ranges as inclusive Julia
`UnitRange`s for `rank` in the requested pencil orientation. The orientation's
full axis is `1:n[axis]`; the other two axes are partitioned across `dec.ranks`.
"""
function pencil_bounds(dec::PencilDecomposition3D, rank::Integer, orientation = :x)
    full_axis, split_axes = _pencil_axes(_pencil_val(orientation))
    coords = pencil_rank_coords(dec, rank)
    ranges = ntuple(d -> 1:dec.n[d], 3)
    ranges = Base.setindex(
        ranges,
        _partition_range(dec.n[split_axes[1]], dec.ranks[1], coords[1]),
        split_axes[1],
    )
    ranges = Base.setindex(
        ranges,
        _partition_range(dec.n[split_axes[2]], dec.ranks[2], coords[2]),
        split_axes[2],
    )
    ranges = Base.setindex(ranges, 1:dec.n[full_axis], full_axis)
    return ranges
end

"""
    pencil_local_size(dec, rank, orientation=:x) -> NTuple{3,Int}

Local array shape for a rank-local pencil in the requested orientation.
"""
function pencil_local_size(dec::PencilDecomposition3D, rank::Integer, orientation = :x)
    return map(length, pencil_bounds(dec, rank, orientation))
end

"""
    pencil_owner(dec, index, orientation=:x) -> Int

Return the rank that owns global cell `index` in the requested pencil
orientation. Since pencil decompositions here are fully periodic, out-of-range
indices wrap onto the global grid before owner lookup.
"""
function pencil_owner(dec::PencilDecomposition3D, index::NTuple{3,<:Integer}, orientation = :x)
    idx = ntuple(d -> _rank_coordinate("index[$d]", index[d]), 3)
    return _pencil_owner_checked(dec, idx, orientation)
end

function pencil_owner(dec::PencilDecomposition3D, index::Tuple, orientation = :x)
    length(index) == 3 ||
        throw(ArgumentError("pencil owner index must have length 3, got $(length(index))"))
    all(x -> x isa Integer, index) ||
        throw(ArgumentError("pencil owner index entries must be integers, got $index"))
    idx = ntuple(d -> _rank_coordinate("index[$d]", index[d]), 3)
    return _pencil_owner_checked(dec, idx, orientation)
end

function _pencil_owner_checked(dec::PencilDecomposition3D, index::NTuple{3,Int}, orientation)
    _, split_axes = _pencil_axes(_pencil_val(orientation))
    coords = (
        _partition_owner(index[split_axes[1]], dec.n[split_axes[1]], dec.ranks[1]),
        _partition_owner(index[split_axes[2]], dec.n[split_axes[2]], dec.ranks[2]),
    )
    return pencil_rank_index(dec, coords)
end

"""
    pencil_orientation_axes(orientation) -> (; full_axis, split_axes)

Return the full and distributed axes for `orientation`.
"""
function pencil_orientation_axes(orientation = :x)
    full_axis, split_axes = _pencil_axes(_pencil_val(orientation))
    return (; full_axis, split_axes)
end
