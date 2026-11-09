# Serial Cartesian dyadic forest. Only valid leaves are physical integration cells.
struct _CartesianAMRPatch{D}
    parent::Int
    level::Int
    cells::Vector{Tuple{Int,NTuple{D,Int}}}
    valid::Vector{Bool}
    covered::Vector{Bool}
    leaves::Vector{Int}
end
struct _CartesianAMRFace{D}
    axis::Int
    left::Int
    right::Int
    lower::NTuple{D,_CompositeQ}
    upper::NTuple{D,_CompositeQ}
    area::_CompositeQ
end
struct CartesianAMRHierarchy{D,T<:Union{Float32,Float64},R<:Tuple}
    nbase::NTuple{D,Int}
    L::NTuple{D,T}
    refinements::R
    patches::Vector{_CartesianAMRPatch{D}}
    addresses::Vector{Tuple{Int,NTuple{D,Int}}}
    lower::Vector{NTuple{D,_CompositeQ}}
    upper::Vector{NTuple{D,_CompositeQ}}
    volumes::Vector{_CompositeQ}
    vertices::Vector{NTuple{D,_CompositeQ}}
    corners::Vector{Vector{Int}}
    rows::Vector{Vector{Pair{Int,_CompositeQ}}}
    masters::Vector{Int}
    free::Vector{Int}
    boundary::Vector{Int}
    free_index::Vector{Int}
    faces::Vector{_CartesianAMRFace{D}}
    refined::Dict{Tuple{Int,NTuple{D,Int}},Int}
    leaf_index::Dict{Tuple{Int,NTuple{D,Int}},Int}
end
_camr_bits(::Val{D}) where {D} = [ntuple(d -> (k >> (d - 1)) & 1, D) for k = 0:((1 << D)-1)]
function _camr_indices(n::NTuple{D,Int}) where {D}
    count = 1
    for value in n
        count = Base.checked_mul(count, value)
    end
    result = NTuple{D,Int}[]
    sizehint!(result, count)
    for I in CartesianIndices(n)
        push!(result, Tuple(I))
    end
    return result
end
_camr_parent(a) = (a[1] - 1, map(i -> i ÷ 2 + Int(isodd(i)), a[2]))
function _camr_children(a::Tuple{Int,NTuple{D,Int}}) where {D}
    level = Base.checked_add(a[1], 1)
    return [(level, ntuple(d -> Base.checked_add(Base.checked_mul(2, a[2][d]) - 1, b[d]), D)) for b in _camr_bits(Val(D))]
end
function _camr_specs(n::NTuple{D,Int}, refinements) where {D}
    refinements isa Union{Tuple,AbstractVector} || throw(ArgumentError("Cartesian refinements require a tuple or vector"))
    refinements isa AbstractVector && Base.require_one_based_indexing(refinements)
    result = Tuple{Int,NTuple{D,Int}}[]
    for record in refinements
        record isa Tuple && length(record) == 2 && record[1] isa Integer &&
        record[2] isa Tuple && length(record[2]) == D && all(i -> i isa Integer, record[2]) ||
            throw(ArgumentError("Cartesian refinement is (level, D-tuple of logical cell indices)"))
        0 <= record[1] <= typemax(Int) && all(i -> 0 < i <= typemax(Int), record[2]) ||
            throw(ArgumentError("Cartesian refinement indices must fit Int"))
        push!(result, (Int(record[1]), map(Int, record[2])))
    end
    sort!(result)
    length(unique(result)) == length(result) || throw(ArgumentError("duplicate Cartesian refinement address"))
    maximum_level = isempty(result) ? 0 : result[end][1]
    denominators = NTuple{D,Int}[n]
    for level = 1:(isempty(result) ? 0 : Base.checked_add(maximum_level, 1))
        push!(denominators, map(i -> Base.checked_mul(2, i), denominators[end]))
    end
    seen = Set{Tuple{Int,NTuple{D,Int}}}()
    for a in result
        all(d -> a[2][d] <= denominators[a[1]+1][d], 1:D) || throw(ArgumentError("Cartesian refinement outside its level"))
        (a[1] == 0 || _camr_parent(a) in seen) || throw(ArgumentError("Cartesian refinement is missing an ancestor"))
        push!(seen, a)
    end
    return Tuple(result), denominators
end
function _camr_bounds(a::Tuple{Int,NTuple{D,Int}}, L, denominators) where {D}
    h = ntuple(d -> _exact_rational(L[d]) / denominators[a[1]+1][d], D)
    lower = ntuple(d -> h[d] * (a[2][d] - 1), D)
    return lower, ntuple(d -> lower[d] + h[d], D)
end
function _camr_query_data(n::NTuple{D,Int}, L, refined, leaves, lo, hi; closed = false, roots = nothing, segment = nothing) where {D}
    any(d -> hi[d] < 0 || lo[d] > _exact_rational(L[d]) || lo[d] > hi[d], 1:D) && return Int[]
    ranges = ntuple(D) do d
        h = _exact_rational(L[d]) / n[d]
        low = max(zero(_CompositeQ), lo[d]) / h
        high = min(_exact_rational(L[d]), hi[d]) / h
        first = closed ? max(1, ceil(Int, low)) : min(n[d], floor(Int, low) + 1)
        last = closed ? min(n[d], floor(Int, high) + 1) : min(n[d], ceil(Int, high))
        first:last
    end
    out = Int[]
    function visit(a, lower, upper)
        intersects = all(d -> closed ? lo[d] <= upper[d] && hi[d] >= lower[d] : lo[d] < upper[d] && hi[d] > lower[d], 1:D)
        intersects || return nothing
        if segment !== nothing
            start, finish = segment
            first, last = zero(_CompositeQ), one(_CompositeQ)
            for d = 1:D
                delta = finish[d] - start[d]
                if iszero(delta)
                    lower[d] <= start[d] <= upper[d] || return nothing
                else
                    a0, a1 = (lower[d] - start[d]) / delta, (upper[d] - start[d]) / delta
                    first, last = max(first, min(a0, a1)), min(last, max(a0, a1))
                    first <= last || return nothing
                end
            end
        end
        if haskey(refined, a)
            middle = ntuple(d -> (lower[d] + upper[d]) / 2, D)
            for (child, bit) in zip(_camr_children(a), _camr_bits(Val(D)))
                l = ntuple(d -> iszero(bit[d]) ? lower[d] : middle[d], D)
                u = ntuple(d -> iszero(bit[d]) ? middle[d] : upper[d], D)
                visit(child, l, u)
            end
        else
            push!(out, leaves[a])
        end
        return nothing
    end
    for I in (roots === nothing ? Iterators.product(ranges...) : roots)
        lower = ntuple(d -> _exact_rational(L[d]) * (I[d] - 1) / n[d], D)
        upper = ntuple(d -> _exact_rational(L[d]) * I[d] / n[d], D)
        visit((0, Tuple(I)), lower, upper)
    end
    sort!(out)
    return out
end
_camr_query(mesh, lo, hi; closed = false) = _camr_query_data(mesh.nbase, mesh.L, mesh.refined, mesh.leaf_index, lo, hi; closed)
function _camr_locate(mesh::CartesianAMRHierarchy{D}, x; trace = false) where {D}
    all(d -> trace ? 0 <= x[d] <= _exact_rational(mesh.L[d]) : 0 < x[d] < _exact_rational(mesh.L[d]), 1:D) ||
        throw(ArgumentError("Cartesian point is outside its admitted physical box"))
    index = ntuple(d -> min(mesh.nbase[d], floor(Int, x[d] * mesh.nbase[d] / _exact_rational(mesh.L[d])) + 1), D)
    a = (0, index)
    while haskey(mesh.refined, a)
        level = Base.checked_add(a[1], 1)
        a = (level, ntuple(d -> begin
            base = Base.checked_mul(2, a[2][d])
            midpoint = _exact_rational(mesh.L[d]) * (base - 1) / (BigInt(mesh.nbase[d]) << level)
            x[d] < midpoint ? base - 1 : base
        end, D))
    end
    return mesh.leaf_index[a]
end
function _camr_row(corners, rows, lower, upper, x)
    D = length(lower)
    result = Dict{Int,_CompositeQ}()
    for (raw, bit) in zip(corners, _camr_bits(Val(D)))
        value = prod(d -> iszero(bit[d]) ? (upper[d] - x[d]) / (upper[d] - lower[d]) : (x[d] - lower[d]) / (upper[d] - lower[d]), 1:D)
        iszero(value) && continue
        for (master, coefficient) in rows[raw]
            result[master] = get(result, master, zero(_CompositeQ)) + value * coefficient
        end
    end
    return sort!([i => value for (i, value) in result if !iszero(value)]; by = first)
end
function _camr_layout(n::NTuple{D,Int}, L::NTuple{D,T}, specs, denominators) where {D,T}
    refined = Dict{Tuple{Int,NTuple{D,Int}},Int}(a => k + 1 for (k, a) in enumerate(specs))
    patchcells = Vector{Vector{Tuple{Int,NTuple{D,Int}}}}()
    push!(patchcells, [(0, I) for I in _camr_indices(n)])
    for a in specs
        push!(patchcells, _camr_children(a))
    end
    total = 0
    addresses = Tuple{Int,NTuple{D,Int}}[]
    for cells in patchcells
        total = Base.checked_add(total, length(cells))
        for a in cells
            haskey(refined, a) || push!(addresses, a)
        end
    end
    Base.checked_add(total, 1)
    sort!(addresses; by = a -> (_camr_bounds(a, L, denominators)[1], a[1]))
    leaf_index = Dict(a => k for (k, a) in enumerate(addresses))
    lower, upper = NTuple{D,_CompositeQ}[], NTuple{D,_CompositeQ}[]
    for a in addresses
        l, u = _camr_bounds(a, L, denominators)
        push!(lower, l); push!(upper, u)
    end
    volumes = [prod(d -> upper[k][d] - lower[k][d], 1:D) for k in eachindex(addresses)]
    sum(volumes; init = zero(_CompositeQ)) == prod(_exact_rational, L) || throw(ArgumentError("Cartesian valid volume does not cover its box"))
    patches = _CartesianAMRPatch{D}[]
    for (k, cells) in enumerate(patchcells)
        parent = k == 1 ? 0 : (specs[k-1][1] == 0 ? 1 : refined[_camr_parent(specs[k-1])])
        covered = [haskey(refined, a) for a in cells]
        push!(patches, _CartesianAMRPatch(parent, cells[1][1], cells, .!covered, covered, [get(leaf_index, a, 0) for a in cells]))
    end
    query(lo, hi) = _camr_query_data(n, L, refined, leaf_index, lo, hi; closed = true)
    faces = _CartesianAMRFace{D}[]
    for k in eachindex(addresses)
        for j in query(lower[k], upper[k])
            abs(addresses[k][1] - addresses[j][1]) <= 1 || throw(ArgumentError("Cartesian face/edge/corner neighbors violate 2:1 balance"))
        end
        for d = 1:D
            area = prod(j -> j == d ? one(_CompositeQ) : upper[k][j] - lower[k][j], 1:D)
            if iszero(lower[k][d])
                lo = lower[k]; hi = ntuple(j -> j == d ? lo[j] : upper[k][j], D)
                push!(faces, _CartesianAMRFace(d, 0, k, lo, hi, area))
            end
            if upper[k][d] == _exact_rational(L[d])
                hi = upper[k]; lo = ntuple(j -> j == d ? hi[j] : lower[k][j], D)
                push!(faces, _CartesianAMRFace(d, k, 0, lo, hi, area))
            else
                lo = ntuple(j -> j == d ? upper[k][j] : lower[k][j], D)
                for j in query(lo, upper[k])
                    lower[j][d] == upper[k][d] || continue
                    l = ntuple(c -> max(lo[c], lower[j][c]), D)
                    u = ntuple(c -> min(upper[k][c], upper[j][c]), D)
                    all(c -> c == d || l[c] < u[c], 1:D) || continue
                    a = prod(c -> c == d ? one(_CompositeQ) : u[c] - l[c], 1:D)
                    push!(faces, _CartesianAMRFace(d, k, j, l, u, a))
                end
            end
        end
    end
    points = NTuple{D,_CompositeQ}[]
    for k in eachindex(addresses), bit in _camr_bits(Val(D))
        push!(points, ntuple(d -> iszero(bit[d]) ? lower[k][d] : upper[k][d], D))
    end
    vertices = sort!(unique(points))
    Base.checked_add(length(vertices), 1)
    rawindex = Dict(v => k for (k, v) in enumerate(vertices))
    corners = [[rawindex[ntuple(d -> iszero(bit[d]) ? lower[k][d] : upper[k][d], D)] for bit in _camr_bits(Val(D))] for k in eachindex(addresses)]
    owners = zeros(Int, length(vertices))
    for (v, x) in enumerate(vertices)
        candidates = [k for k in query(x, x) if any(d -> x[d] != lower[k][d] && x[d] != upper[k][d], 1:D)]
        isempty(candidates) || (owners[v] = first(sort!(candidates; by = k -> addresses[k])))
    end
    masters = findall(iszero, owners)
    masterindex = zeros(Int, length(vertices))
    for (i, raw) in enumerate(masters)
        masterindex[raw] = i
    end
    births = [first(level for level = 0:length(denominators)-1 if all(d -> denominator(x[d] * denominators[level+1][d] / _exact_rational(L[d])) == 1, 1:D)) for x in vertices]
    rows = [Pair{Int,_CompositeQ}[] for _ in vertices]
    visiting = falses(length(vertices))
    function expand(v)
        isempty(rows[v]) || return rows[v]
        visiting[v] && throw(ArgumentError("cyclic Cartesian hanging constraint"))
        visiting[v] = true
        if iszero(owners[v])
            rows[v] = [masterindex[v] => one(_CompositeQ)]
        else
            k = owners[v]
            for raw in corners[k]
                births[raw] < births[v] || throw(ArgumentError("Cartesian hanging dependency does not lower vertex birth level"))
                expand(raw)
            end
            rows[v] = _camr_row(corners[k], rows, lower[k], upper[k], vertices[v])
        end
        visiting[v] = false
        all(pair -> pair.second >= 0, rows[v]) && sum(last, rows[v]; init = zero(_CompositeQ)) == 1 || throw(ArgumentError("Cartesian constraint lost positive partition of unity"))
        return rows[v]
    end
    foreach(expand, eachindex(vertices))
    for (v, x) in enumerate(vertices), k in query(x, x)
        isequal(rows[v], _camr_row(corners[k], rows, lower[k], upper[k], x)) || throw(ArgumentError("Cartesian edge/corner constraint trace mismatch"))
    end
    for face in faces
        if face.left > 0 && face.right > 0
            for bit in _camr_bits(Val(D))
                x = ntuple(d -> iszero(bit[d]) ? face.lower[d] : face.upper[d], D)
                left = _camr_row(corners[face.left], rows, lower[face.left], upper[face.left], x)
                right = _camr_row(corners[face.right], rows, lower[face.right], upper[face.right], x)
                isequal(left, right) || throw(ArgumentError("Cartesian coarse/fine face trace mismatch"))
            end
        end
    end
    boundary = [i for (i, raw) in enumerate(masters) if any(d -> iszero(vertices[raw][d]) || vertices[raw][d] == _exact_rational(L[d]), 1:D)]
    boundaryset = Set(boundary)
    free = [i for i in eachindex(masters) if !(i in boundaryset)]
    isempty(free) && throw(ArgumentError("Cartesian grounded hierarchy has no free potential"))
    free_index = zeros(Int, length(masters))
    for (i, master) in enumerate(free)
        free_index[master] = i
    end
    for (v, x) in enumerate(vertices)
        if any(d -> iszero(x[d]) || x[d] == _exact_rational(L[d]), 1:D)
            all(pair -> pair.first in boundaryset, rows[v]) || throw(ArgumentError("Cartesian boundary constraint uses an interior master"))
        end
    end
    return CartesianAMRHierarchy(n, L, specs, patches, addresses, lower, upper, volumes, vertices, corners, rows, masters, free, boundary, free_index, faces, refined, leaf_index)
end
"Owned balanced 2D/3D dyadic hierarchy; each refinement is (level, logical-cell tuple)."
function CartesianAMRHierarchy(nbase::NTuple{D,Integer}, L::NTuple{D,T}, refinements = ()) where {D,T<:Union{Float32,Float64}}
    D in (2, 3) || throw(ArgumentError("Cartesian AMR requires D=2 or D=3"))
    all(i -> 2 <= i < typemax(Int), nbase) || throw(ArgumentError("Cartesian base counts must be >=2 and fit node indices"))
    all(x -> isfinite(x) && _stored_positive(x), L) || throw(ArgumentError("Cartesian lengths must be finite positive original IEEE values"))
    n = map(Int, nbase)
    specs, denominators = _camr_specs(n, refinements)
    return _camr_layout(n, L, specs, denominators)
end
function _camr_record(mesh)
    patch(p) = (p.parent, p.level, p.cells, p.valid, p.covered, p.leaves)
    face(f) = (f.axis, f.left, f.right, f.lower, f.upper, f.area)
    return (mesh.nbase, mesh.L, mesh.refinements, map(patch, mesh.patches), mesh.addresses, mesh.lower, mesh.upper, mesh.volumes, mesh.vertices, mesh.corners, mesh.rows, mesh.masters, mesh.free, mesh.boundary, mesh.free_index, map(face, mesh.faces), mesh.refined, mesh.leaf_index)
end
function _camr_validate(mesh::CartesianAMRHierarchy)
    expected = CartesianAMRHierarchy(mesh.nbase, mesh.L, mesh.refinements)
    isequal(_camr_record(mesh), _camr_record(expected)) || throw(ArgumentError("Cartesian AMR derived geometry/constraint cache changed"))
    return nothing
end
function _camr_arrays(mesh)
    arrays = Any[mesh.patches, mesh.addresses, mesh.lower, mesh.upper, mesh.volumes, mesh.vertices, mesh.corners, mesh.rows, mesh.masters, mesh.free, mesh.boundary, mesh.free_index, mesh.faces]
    append!(arrays, mesh.corners); append!(arrays, mesh.rows)
    for p in mesh.patches
        append!(arrays, (p.cells, p.valid, p.covered, p.leaves))
    end
    return arrays
end

# Root-cell ray traversal plus exact slab pruning at every refined tree node.
function _camr_segment(mesh::CartesianAMRHierarchy{D}, x, y) where {D}
    parameters = _CompositeQ[0, 1]
    for d = 1:D
        delta = y[d] - x[d]
        iszero(delta) && continue
        scale = mesh.nbase[d] / _exact_rational(mesh.L[d])
        first = max(1, ceil(Int, min(x[d], y[d]) * scale))
        last = min(mesh.nbase[d] - 1, floor(Int, max(x[d], y[d]) * scale))
        for j = first:last
            push!(parameters, (j / scale - x[d]) / delta)
        end
    end
    sort!(unique!(parameters))
    samples = copy(parameters)
    append!(samples, [(parameters[i] + parameters[i+1]) / 2 for i = 1:length(parameters)-1])
    roots = Set{NTuple{D,Int}}()
    for parameter in samples
        point = ntuple(d -> x[d] + parameter * (y[d] - x[d]), D)
        ranges = ntuple(d -> begin
            coordinate = point[d] * mesh.nbase[d] / _exact_rational(mesh.L[d])
            max(1, ceil(Int, coordinate)):min(mesh.nbase[d], floor(Int, coordinate) + 1)
        end, D)
        for I in Iterators.product(ranges...)
            push!(roots, Tuple(I))
        end
    end
    lo, hi = ntuple(d -> min(x[d], y[d]), D), ntuple(d -> max(x[d], y[d]), D)
    return _camr_query_data(mesh.nbase, mesh.L, mesh.refined, mesh.leaf_index, lo, hi; closed = true, roots = sort!(collect(roots)), segment = (x, y))
end
