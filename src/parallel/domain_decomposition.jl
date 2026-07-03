# domain_decomposition.jl (§5.3 parallel/)
#
# Logical Cartesian rank decomposition used by the particle-migration core. This
# is transport-agnostic: the same rank/bounds/destination rules are exercised in
# serial tests and can be called from an MPI layer without changing particle
# semantics.

"""
    LogicalRankLayout(ranks; periodic=ntuple(_ -> true, length(ranks)))

Cartesian logical-rank layout for a `D`-dimensional physical domain. `ranks[d]`
is the number of subdomains along axis `d`; ranks are 1-based and stored in
column-major order, matching Julia's `LinearIndices`.
"""
struct LogicalRankLayout{D}
    ranks::NTuple{D,Int}
    periodic::NTuple{D,Bool}
end

function LogicalRankLayout(ranks::NTuple{D,<:Integer}; periodic = ntuple(_ -> true, D)) where {D}
    D >= 1 || throw(ArgumentError("rank layout dimension must be >= 1"))
    rr = ntuple(d -> _require_positive_intlike("ranks[$d]", ranks[d]), D)
    all(>(0), rr) || throw(ArgumentError("all rank counts must be positive, got $rr"))
    length(periodic) == D || throw(ArgumentError("periodic length must equal rank dimension $D"))
    pp = ntuple(d -> Bool(periodic[d]), D)
    return LogicalRankLayout{D}(rr, pp)
end

nranks(layout::LogicalRankLayout) = prod(layout.ranks)

function _check_rank(layout::LogicalRankLayout, rank::Integer)
    1 <= rank <= nranks(layout) ||
        throw(ArgumentError("rank must be in 1:$(nranks(layout)), got $rank"))
    r = Int(rank)
    return r
end

function _rank_coordinate(name::AbstractString, coord::Integer)
    typemin(Int) <= coord <= typemax(Int) || throw(ArgumentError("$name must fit in Int"))
    return Int(coord)
end

"""
    rank_coords(layout, rank) -> NTuple{D,Int}

Return the 1-based Cartesian coordinate of `rank` in `layout`.
"""
function rank_coords(layout::LogicalRankLayout{D}, rank::Integer) where {D}
    r = _check_rank(layout, rank)
    return Tuple(CartesianIndices(layout.ranks)[r])
end

"""
    rank_index(layout, coords) -> Union{Int,Nothing}

Return the 1-based linear rank for 1-based Cartesian `coords`. Periodic axes
wrap. Nonperiodic out-of-range coordinates return `nothing`.
"""
function rank_index(layout::LogicalRankLayout{D}, coords::NTuple{D,<:Integer}) where {D}
    cc = ntuple(d -> _rank_coordinate("coords[$d]", coords[d]), D)
    wrapped = ntuple(d -> begin
        c = cc[d]
        if layout.periodic[d]
            mod1(c, layout.ranks[d])
        elseif 1 <= c <= layout.ranks[d]
            c
        else
            0
        end
    end, D)
    any(==(0), wrapped) && return nothing
    return LinearIndices(layout.ranks)[wrapped...]
end

"""
    rank_bounds(g, layout, rank) -> (; lo, hi)

Physical half-open bounds `[lo[d], hi[d])` for `rank` on grid `g`. The last
subdomain on each axis ends exactly at `g.L[d]`.
"""
function rank_bounds(g::FourierGrid{D,T}, layout::LogicalRankLayout{D}, rank::Integer) where {D,T}
    coords = rank_coords(layout, rank)
    lo = ntuple(d -> T((coords[d] - 1) * g.L[d] / layout.ranks[d]), D)
    hi = ntuple(d -> T(coords[d] * g.L[d] / layout.ranks[d]), D)
    return (; lo, hi)
end

"""
    rank_of_position(x, g, layout) -> Union{Int,Nothing}

Destination rank for physical position tuple `x`. Periodic axes are classified
after wrapping into `[0, L)`. Nonperiodic out-of-domain positions return
`nothing`.
"""
function rank_of_position(
    x::NTuple{D,<:Real},
    g::FourierGrid{D,T},
    layout::LogicalRankLayout{D},
) where {D,T}
    coords = ntuple(d -> begin
        xd = T(x[d])
        isfinite(xd) || throw(ArgumentError("x[$d] must be finite"))
        Ld = g.L[d]
        if layout.periodic[d]
            xd = mod(xd, Ld)
            c = floor(Int, xd * layout.ranks[d] / Ld) + 1
            clamp(c, 1, layout.ranks[d])
        elseif xd < zero(T) || xd >= Ld
            0
        else
            c = floor(Int, xd * layout.ranks[d] / Ld) + 1
            clamp(c, 1, layout.ranks[d])
        end
    end, D)
    any(==(0), coords) && return nothing
    return rank_index(layout, coords)
end

function _slab_axis(layout::LogicalRankLayout{D}) where {D}
    axes = findall(>(1), layout.ranks)
    length(axes) <= 1 || throw(
        ArgumentError(
            "slab ghost exchange supports one decomposed axis; got ranks=$(layout.ranks)",
        ),
    )
    return isempty(axes) ? 0 : axes[1]
end

function _face_ranges(A::AbstractArray{T,D}, axis::Integer, r::UnitRange{Int}) where {T,D}
    return ntuple(d -> d == axis ? r : axes(A, d), D)
end

function _validate_halo_arrays!(
    rank_arrays::AbstractVector{<:AbstractArray{T,D}},
    layout::LogicalRankLayout{D},
    halo::Integer,
) where {T,D}
    halo >= 1 || throw(ArgumentError("halo must be >= 1, got $halo"))
    length(rank_arrays) == nranks(layout) ||
        throw(ArgumentError("rank_arrays length must equal nranks(layout)"))
    isempty(rank_arrays) && return nothing
    sz = size(first(rank_arrays))
    for (r, A) in pairs(rank_arrays)
        size(A) == sz || throw(DimensionMismatch("rank $r has size $(size(A)); expected $sz"))
    end
    axis = _slab_axis(layout)
    axis == 0 && return axis
    sz[axis] >= 3halo ||
        throw(DimensionMismatch("axis $axis size $(sz[axis]) is too small for halo=$halo"))
    return axis
end

"""
    exchange_field_halos!(rank_arrays, layout; halo=1, fill_value=zero(T))

Copy slab field halo values between neighboring logical ranks. `rank_arrays[r]`
is a rank-local field array with `halo` ghost cells on both sides of the single
decomposed axis. Lower ghost cells receive the lower neighbor's upper owned
boundary cells; upper ghost cells receive the upper neighbor's lower owned
boundary cells. Periodic rank layouts wrap through [`rank_index`](@ref).
Nonperiodic exterior ghost cells are filled with `fill_value`.

This is a deterministic local reference for MPI field-halo exchange. It does
not implement message passing and deliberately rejects pencil decompositions.
Returns `(; exchanged, filled)`, the number of scalar ghost entries copied from
neighbors or filled at nonperiodic exterior boundaries.
"""
function exchange_field_halos!(
    rank_arrays::AbstractVector{<:AbstractArray{T,D}},
    layout::LogicalRankLayout{D};
    halo::Integer = 1,
    fill_value = zero(T),
) where {T,D}
    h = _require_positive_intlike("halo", halo)
    axis = _validate_halo_arrays!(rank_arrays, layout, h)
    axis == 0 && return (; exchanged = 0, filled = 0)

    transfers = Tuple{Int,NTuple{D,Any},Array{T,D}}[]
    exchanged = 0
    filled = 0

    for r = 1:length(rank_arrays)
        A = rank_arrays[r]
        coords = rank_coords(layout, r)
        nloc = size(A, axis)

        lower_dst = _face_ranges(A, axis, 1:h)
        upper_dst = _face_ranges(A, axis, (nloc-h+1):nloc)
        lower_src_rank = rank_index(layout, ntuple(d -> d == axis ? coords[d] - 1 : coords[d], D))
        upper_src_rank = rank_index(layout, ntuple(d -> d == axis ? coords[d] + 1 : coords[d], D))

        if lower_src_rank === nothing
            A[lower_dst...] .= fill_value
            filled += length(view(A, lower_dst...))
        else
            B = rank_arrays[lower_src_rank]
            src = _face_ranges(B, axis, (size(B, axis)-2h+1):(size(B, axis)-h))
            push!(transfers, (r, lower_dst, copy(view(B, src...))))
            exchanged += length(view(B, src...))
        end

        if upper_src_rank === nothing
            A[upper_dst...] .= fill_value
            filled += length(view(A, upper_dst...))
        else
            B = rank_arrays[upper_src_rank]
            src = _face_ranges(B, axis, (h+1):(2h))
            push!(transfers, (r, upper_dst, copy(view(B, src...))))
            exchanged += length(view(B, src...))
        end
    end

    for (dest, dst, payload) in transfers
        rank_arrays[dest][dst...] .= payload
    end

    return (; exchanged, filled)
end

function exchange_field_halos!(
    rank_fields::AbstractVector{<:NTuple{N,<:AbstractArray{T,D}}},
    layout::LogicalRankLayout{D};
    halo::Integer = 1,
    fill_value = zero(T),
) where {N,T,D}
    exchanged = 0
    filled = 0
    for c = 1:N
        arrays = [rank_fields[r][c] for r = 1:length(rank_fields)]
        stats = exchange_field_halos!(arrays, layout; halo, fill_value)
        exchanged += stats.exchanged
        filled += stats.filled
    end
    return (; exchanged, filled)
end

"""
    exchange_ghost_moments!(rank_arrays, layout; halo=1, clear_ghosts=true)

Accumulate slab ghost-zone moment contributions into neighboring rank interiors.
`rank_arrays[r]` is a rank-local moment array with `halo` ghost cells on both
sides of the single decomposed axis. Lower ghost cells are added to the lower
neighbor's upper owned boundary cells; upper ghost cells are added to the upper
neighbor's lower owned boundary cells. Periodic rank layouts wrap through
[`rank_index`](@ref); nonperiodic exterior ghost contributions are dropped.

This is a deterministic local reference for MPI ghost-moment exchange. It does
not implement message passing and deliberately rejects pencil decompositions.
Returns `(; exchanged, dropped)`, the number of scalar ghost entries added or
dropped.
"""
function exchange_ghost_moments!(
    rank_arrays::AbstractVector{<:AbstractArray{T,D}},
    layout::LogicalRankLayout{D};
    halo::Integer = 1,
    clear_ghosts::Bool = true,
) where {T,D}
    h = Int(halo)
    axis = _validate_halo_arrays!(rank_arrays, layout, h)
    axis == 0 && return (; exchanged = 0, dropped = 0)

    transfers = Tuple{Int,NTuple{D,Any},Array{T,D},NTuple{D,Any}}[]
    dropped = 0
    exchanged = 0

    for r = 1:length(rank_arrays)
        A = rank_arrays[r]
        coords = rank_coords(layout, r)
        nloc = size(A, axis)
        lower_src = _face_ranges(A, axis, 1:h)
        upper_src = _face_ranges(A, axis, (nloc-h+1):nloc)

        lower_dest = rank_index(layout, ntuple(d -> d == axis ? coords[d] - 1 : coords[d], D))
        upper_dest = rank_index(layout, ntuple(d -> d == axis ? coords[d] + 1 : coords[d], D))

        if lower_dest === nothing
            dropped += length(view(A, lower_src...))
        else
            B = rank_arrays[lower_dest]
            target = _face_ranges(B, axis, (size(B, axis)-2h+1):(size(B, axis)-h))
            push!(transfers, (lower_dest, target, copy(view(A, lower_src...)), lower_src))
            exchanged += length(view(A, lower_src...))
        end

        if upper_dest === nothing
            dropped += length(view(A, upper_src...))
        else
            B = rank_arrays[upper_dest]
            target = _face_ranges(B, axis, (h+1):(2h))
            push!(transfers, (upper_dest, target, copy(view(A, upper_src...)), upper_src))
            exchanged += length(view(A, upper_src...))
        end
    end

    for (dest, target, payload, _) in transfers
        rank_arrays[dest][target...] .+= payload
    end

    if clear_ghosts
        for A in rank_arrays
            nloc = size(A, axis)
            A[_face_ranges(A, axis, 1:h)...] .= zero(T)
            A[_face_ranges(A, axis, (nloc-h+1):nloc)...] .= zero(T)
        end
    end

    return (; exchanged, dropped)
end

function exchange_ghost_moments!(
    rank_moments::AbstractVector{<:NTuple{N,<:AbstractArray{T,D}}},
    layout::LogicalRankLayout{D};
    halo::Integer = 1,
    clear_ghosts::Bool = true,
) where {N,T,D}
    exchanged = 0
    dropped = 0
    for c = 1:N
        arrays = [rank_moments[r][c] for r = 1:length(rank_moments)]
        stats = exchange_ghost_moments!(arrays, layout; halo, clear_ghosts)
        exchanged += stats.exchanged
        dropped += stats.dropped
    end
    return (; exchanged, dropped)
end

function _diagnostic_reduction_op(op::Symbol)
    op === :sum && return +
    op === :min && return min
    op === :max && return max
    throw(ArgumentError("diagnostic reduction op must be :sum, :min, or :max, got $op"))
end

function _reduce_diagnostic_values(values::AbstractVector, op)
    isempty(values) && throw(ArgumentError("cannot reduce an empty diagnostic collection"))
    first_value = first(values)

    if first_value isa NamedTuple
        key_tuple = keys(first_value)
        for value in values
            value isa NamedTuple && keys(value) == key_tuple ||
                throw(ArgumentError("all diagnostic NamedTuples must have the same keys"))
        end
        reduced = ntuple(
            i -> _reduce_diagnostic_values([getfield(value, key_tuple[i]) for value in values], op),
            length(key_tuple),
        )
        return NamedTuple{key_tuple}(reduced)
    elseif first_value isa Tuple
        n = length(first_value)
        for value in values
            value isa Tuple && length(value) == n ||
                throw(ArgumentError("all diagnostic Tuples must have the same length"))
        end
        return ntuple(i -> _reduce_diagnostic_values([value[i] for value in values], op), n)
    elseif first_value isa AbstractArray
        ax = axes(first_value)
        out = copy(first_value)
        for value in Iterators.drop(values, 1)
            value isa AbstractArray ||
                throw(ArgumentError("all diagnostic array leaves must be arrays"))
            axes(value) == ax ||
                throw(DimensionMismatch("diagnostic array axes $(axes(value)) do not match $ax"))
            out .= op.(out, value)
        end
        return out
    elseif first_value isa Number
        for value in values
            value isa Number || throw(ArgumentError("all diagnostic scalar leaves must be numbers"))
        end
        return reduce(op, values)
    else
        throw(ArgumentError("unsupported diagnostic reduction leaf type $(typeof(first_value))"))
    end
end

"""
    reduce_diagnostics(local_values; op=:sum)

Reduce rank-local diagnostic values with deterministic local semantics. Supports
numeric scalars, arrays with identical axes, tuples, and named tuples composed
of those leaves. `op` may be `:sum`, `:min`, or `:max`.

This is the serial reference for MPI diagnostic allreduces: an MPI backend can
replace the outer collection with communicator-local values and preserve these
structure and error semantics.
"""
function reduce_diagnostics(local_values::AbstractVector; op::Symbol = :sum)
    return _reduce_diagnostic_values(local_values, _diagnostic_reduction_op(op))
end

sum_diagnostics(local_values::AbstractVector) = reduce_diagnostics(local_values; op = :sum)
min_diagnostics(local_values::AbstractVector) = reduce_diagnostics(local_values; op = :min)
max_diagnostics(local_values::AbstractVector) = reduce_diagnostics(local_values; op = :max)
