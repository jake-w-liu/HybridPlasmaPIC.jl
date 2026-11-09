# Static 1D interval hierarchy. Only unique leaves are physical field degrees of
# freedom; covered patch cells are derived views and this direct solver has no halo.
const _CompositeQ = Rational{BigInt}

struct _StaticAMRPatch1D
    parent::Int
    level::Int
    first::Int
    ncell::Int
    valid::Vector{Bool}
    covered::Vector{Bool}
    leaf_index::Vector{Int}
    leaf_first::Vector{Int}
    leaf_last::Vector{Int}
    nodes::Vector{Int}

    function _StaticAMRPatch1D(
        parent::Int,
        level::Int,
        first::Int,
        ncell::Int,
        valid::AbstractVector{Bool},
        covered::AbstractVector{Bool},
        leaf_index::Vector{Int},
        leaf_first::Vector{Int},
        leaf_last::Vector{Int},
        nodes::Vector{Int},
    )
        return new(
            parent,
            level,
            first,
            ncell,
            Vector{Bool}(valid),
            Vector{Bool}(covered),
            leaf_index,
            leaf_first,
            leaf_last,
            nodes,
        )
    end
end

struct StaticAMRHierarchy1D{T<:Union{Float32,Float64},S<:Tuple}
    nbase::Int
    L::T
    refinements::S
    patches::Vector{_StaticAMRPatch1D}
    vertices::Vector{_CompositeQ}
    widths::Vector{_CompositeQ}
    dual_volumes::Vector{_CompositeQ}
    owner_patch::Vector{Int}
    owner_cell::Vector{Int}
    interfaces::Vector{Int}
end

function _static_amr_specs(refinements)
    refinements isa Union{Tuple,AbstractVector} ||
        throw(ArgumentError("static refinements must be an ordered tuple or vector"))
    refinements isa AbstractVector && Base.require_one_based_indexing(refinements)
    return Tuple(
        map(refinements) do record
            record isa Tuple && length(record) == 3 && all(x -> x isa Integer, record) || throw(
                ArgumentError("each static refinement must be (parent, first_cell, last_cell)"),
            )
            all(x -> 0 < x <= typemax(Int), record) || throw(
                ArgumentError("static refinement indices must be positive Int-sized integers"),
            )
            (Int(record[1]), Int(record[2]), Int(record[3]))
        end,
    )
end

function _static_amr_layout(nbase::Int, L::T, specs::Tuple) where {T}
    nbase >= 2 || throw(ArgumentError("static hierarchy requires at least two base cells"))
    isfinite(L) && _stored_positive(L) ||
        throw(ArgumentError("static hierarchy length must be finite and positive"))
    Base.checked_add(nbase, 1)
    parents, levels, firsts, counts, denominators = [0], [0], [0], [nbase], [nbase]
    covered = [fill(false, nbase)]
    total = nbase
    for (parent, first, last) in specs
        1 <= parent <= length(counts) ||
            throw(ArgumentError("static refinement parent must precede its child"))
        1 <= first <= last <= counts[parent] ||
            throw(ArgumentError("static refinement must cover existing parent cells"))
        any(view(covered[parent], first:last)) &&
            throw(ArgumentError("static sibling refinements must not overlap"))
        count = Base.checked_mul(2, Base.checked_add(last - first, 1))
        origin = Base.checked_mul(2, Base.checked_add(firsts[parent], first - 1))
        denominator = Base.checked_mul(2, denominators[parent])
        level = Base.checked_add(levels[parent], 1)
        total = Base.checked_add(total, count)
        Base.checked_add(total, 1)
        push!(parents, parent)
        push!(levels, level)
        push!(firsts, origin)
        push!(counts, count)
        push!(denominators, denominator)
        push!(covered, fill(false, count))
        fill!(view(covered[parent], first:last), true)
    end
    length0 = _exact_rational(L)
    leaves = Tuple{_CompositeQ,_CompositeQ,Int,Int}[]
    for patch in eachindex(counts), cell = 1:counts[patch]
        covered[patch][cell] && continue
        j = Base.checked_add(firsts[patch], cell - 1)
        lo = length0 * j / denominators[patch]
        hi = length0 * Base.checked_add(j, 1) / denominators[patch]
        push!(leaves, (lo, hi, patch, cell))
    end
    sort!(leaves; by = x -> x[1])
    isempty(leaves) && throw(ArgumentError("static hierarchy has no valid leaves"))
    vertices = _CompositeQ[zero(_CompositeQ)]
    widths, owner_patch, owner_cell = _CompositeQ[], Int[], Int[]
    previous_level = levels[leaves[1][3]]
    for (lo, hi, patch, cell) in leaves
        lo == vertices[end] && hi > lo ||
            throw(ArgumentError("static leaves must cover the domain once without holes"))
        abs(levels[patch] - previous_level) <= 1 ||
            throw(ArgumentError("static adjacent leaves must satisfy 2:1 balance"))
        push!(vertices, hi)
        push!(widths, hi - lo)
        push!(owner_patch, patch)
        push!(owner_cell, cell)
        previous_level = levels[patch]
    end
    vertices[end] == length0 || throw(ArgumentError("static leaves do not cover the full domain"))
    patches = _StaticAMRPatch1D[]
    for patch in eachindex(counts)
        nodes = Int[]
        for j = 0:counts[patch]
            coordinate = length0 * Base.checked_add(firsts[patch], j) / denominators[patch]
            node = searchsortedfirst(vertices, coordinate)
            node <= length(vertices) && vertices[node] == coordinate ||
                throw(ArgumentError("static patch node is absent from the composite"))
            push!(nodes, node)
        end
        lo, hi, indices = Int[], Int[], Int[]
        for cell = 1:counts[patch]
            push!(lo, nodes[cell])
            push!(hi, nodes[cell+1] - 1)
            nodes[cell] < nodes[cell+1] || throw(ArgumentError("empty static patch cell"))
            valid = !covered[patch][cell]
            valid &&
                nodes[cell+1] != nodes[cell] + 1 &&
                throw(ArgumentError("valid static cell must own exactly one leaf"))
            push!(indices, valid ? nodes[cell] : 0)
        end
        push!(
            patches,
            _StaticAMRPatch1D(
                parents[patch],
                levels[patch],
                firsts[patch],
                counts[patch],
                .!covered[patch],
                covered[patch],
                indices,
                lo,
                hi,
                nodes,
            ),
        )
    end
    dual = _CompositeQ[widths[1]/2]
    for i = 2:length(widths)
        push!(dual, (widths[i-1] + widths[i]) / 2)
    end
    push!(dual, widths[end] / 2)
    interfaces = Int[i + 1 for i = 1:length(widths)-1 if owner_patch[i] != owner_patch[i+1]]
    return (; patches, vertices, widths, dual_volumes = dual, owner_patch, owner_cell, interfaces)
end

"""
    StaticAMRHierarchy1D(nbase, L, refinements=())

Owned static nested 2:1 interval patches on `[0,L]`. Each refinement is the tuple
`(parent_patch, first_parent_cell, last_parent_cell)`; parents precede children.
Only valid leaves contribute to the composite field. Covered cells are derived,
shared nodes occur once, and the direct 1D operator requires zero ghost width.
Geometry is exact from the original Float32/64 length. No dynamic regridding,
multidimensional halos or automatic field/particle transfer is implied.
"""
function StaticAMRHierarchy1D(
    nbase::Integer,
    L::T,
    refinements = (),
) where {T<:Union{Float32,Float64}}
    2 <= nbase <= typemax(Int) ||
        throw(ArgumentError("static base cell count must fit Int and be >=2"))
    specs = _static_amr_specs(refinements)
    layout = _static_amr_layout(Int(nbase), L, specs)
    return StaticAMRHierarchy1D{T,typeof(specs)}(
        Int(nbase),
        L,
        specs,
        layout.patches,
        layout.vertices,
        layout.widths,
        layout.dual_volumes,
        layout.owner_patch,
        layout.owner_cell,
        layout.interfaces,
    )
end

function _static_amr_validate(mesh::StaticAMRHierarchy1D)
    expected = _static_amr_layout(mesh.nbase, mesh.L, mesh.refinements)
    for name in (:vertices, :widths, :dual_volumes, :owner_patch, :owner_cell, :interfaces)
        isequal(getfield(mesh, name), getproperty(expected, name)) ||
            throw(ArgumentError("static hierarchy $name cache changed"))
    end
    length(mesh.patches) == length(expected.patches) ||
        throw(ArgumentError("static hierarchy patch cache changed"))
    for (patch, wanted) in zip(mesh.patches, expected.patches),
        name in fieldnames(_StaticAMRPatch1D)

        isequal(getfield(patch, name), getfield(wanted, name)) ||
            throw(ArgumentError("static hierarchy patch $name cache changed"))
    end
    return nothing
end

function _static_amr_arrays(mesh::StaticAMRHierarchy1D)
    arrays = Any[
        mesh.vertices,
        mesh.widths,
        mesh.dual_volumes,
        mesh.owner_patch,
        mesh.owner_cell,
        mesh.interfaces,
    ]
    for patch in mesh.patches
        append!(
            arrays,
            (
                patch.valid,
                patch.covered,
                patch.leaf_index,
                patch.leaf_first,
                patch.leaf_last,
                patch.nodes,
            ),
        )
    end
    return arrays
end

function _static_amr_leaf(mesh::StaticAMRHierarchy1D, x::_CompositeQ)
    zero(x) < x < mesh.vertices[end] ||
        throw(ArgumentError("composite particle center must lie strictly inside (0,L)"))
    return searchsortedlast(mesh.vertices, x)
end
