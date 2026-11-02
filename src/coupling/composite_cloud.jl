# Exact physical cloud geometry on the unique unequal leaf intervals. The
# primitive and original macrocharge remain the existing closed-cloud owners.
struct _CompositeParticleBins
    offsets::Vector{Int}
    species::Vector{Int}
    particle::Vector{Int}
end

function _sce_cloud_span(mesh::StaticAMRHierarchy1D, x::_CompositeQ, y::_CompositeQ, a::_CompositeQ)
    lower = max(zero(a), min(x, y) - a)
    upper = min(mesh.vertices[end], max(x, y) + a)
    first = max(1, searchsortedlast(mesh.vertices, lower))
    last = min(length(mesh.widths), searchsortedfirst(mesh.vertices, upper) - 1)
    first <= last || throw(ArgumentError("composite cloud has no physical support interval"))
    return first:last
end

function _sce_cloud_charge!(out::Vector{_CompositeQ}, species, mesh, shape)
    fill!(out, zero(_CompositeQ))
    a = _exact_rational(shape.half_width)
    for ps in species, p in eachindex(ps.weight)
        charge = _yee_macrocharge(ps, p, Val(true))
        iszero(charge) && continue
        x = _exact_rational(ps.x[1][p])
        span = _sce_cloud_span(mesh, x, x, a)
        previous = one(_CompositeQ)
        for cell in span
            cumulative, _ = _yee_cloud_edge(x, mesh.vertices[cell], mesh.widths[cell], a)
            out[cell] += charge * (previous - cumulative)
            previous = cumulative
        end
        out[last(span)+1] += charge * previous
    end
    return out
end

function _sce_source_rhs!(out, species, mesh, shape, background, left, right)
    _sce_cloud_charge!(out, species, mesh, shape)
    b = _exact_rational(background)
    for i in eachindex(out)
        out[i] += b * mesh.dual_volumes[i]
    end
    out[1] += _exact_rational(left)
    out[end] += right
    sum(out; init = zero(_CompositeQ)) == 0 ||
        throw(ArgumentError("composite Neumann source must have exact zero total including plates"))
    return out
end

function _sce_cloud_gather(E, x::_CompositeQ, mesh, shape)
    a = _exact_rational(shape.half_width)
    total = zero(_CompositeQ)
    for cell in _sce_cloud_span(mesh, x, x, a)
        _, derivative = _yee_cloud_edge(x, mesh.vertices[cell], mesh.widths[cell], a)
        total += mesh.widths[cell] * _exact_rational(E[cell]) * derivative
    end
    return total
end

function _sce_cloud_impulse!(out, old, new, mesh, shape)
    fill!(out, zero(_CompositeQ))
    a = _exact_rational(shape.half_width)
    for (before, after) in zip(old, new), p in eachindex(before.weight)
        charge = _yee_macrocharge(before, p, Val(true))
        iszero(charge) && continue
        x, y = _exact_rational(before.x[1][p]), _exact_rational(after.x[1][p])
        for cell in _sce_cloud_span(mesh, x, y, a)
            left, _ = _yee_cloud_edge(x, mesh.vertices[cell], mesh.widths[cell], a)
            right, _ = _yee_cloud_edge(y, mesh.vertices[cell], mesh.widths[cell], a)
            out[cell] += charge * (right - left)
        end
    end
    return out
end

function _sce_bins(species, mesh)
    leaves = length(mesh.widths)
    counts = zeros(Int, leaves)
    local_leaf = Vector{Vector{Int}}()
    total = 0
    for ps in species
        total = Base.checked_add(total, length(ps.weight))
        Base.checked_add(total, 1)
        indices = Vector{Int}(undef, length(ps.weight))
        for p in eachindex(ps.weight)
            leaf = _static_amr_leaf(mesh, _exact_rational(ps.x[1][p]))
            counts[leaf] = Base.checked_add(counts[leaf], 1)
            indices[p] = leaf
        end
        push!(local_leaf, indices)
    end
    offsets = Vector{Int}(undef, Base.checked_add(leaves, 1))
    offsets[1] = 1
    for cell = 1:leaves
        offsets[cell+1] = Base.checked_add(offsets[cell], counts[cell])
    end
    slots = copy(offsets[1:leaves])
    owners, particles = Vector{Int}(undef, total), Vector{Int}(undef, total)
    for (s, indices) in enumerate(local_leaf), p in eachindex(indices)
        cell = indices[p]
        slot = slots[cell]
        owners[slot], particles[slot] = s, p
        slots[cell] += 1
    end
    return _CompositeParticleBins(offsets, owners, particles)
end

function _sce_crossings(old, new, mesh)
    out = zeros(Int, length(mesh.interfaces))
    boundaries = mesh.vertices[mesh.interfaces]
    for (before, after) in zip(old, new), p in eachindex(before.weight)
        x, y = _exact_rational(before.x[1][p]), _exact_rational(after.x[1][p])
        x == y && continue
        first = searchsortedlast(boundaries, min(x, y)) + 1
        last = searchsortedlast(boundaries, max(x, y))
        for i = first:last
            out[i] = Base.checked_add(out[i], 1)
        end
    end
    return out
end

function _sce_validate_bins(bins::_CompositeParticleBins, species, mesh)
    total = 0
    seen = Vector{Vector{Bool}}()
    for ps in species
        total = Base.checked_add(total, length(ps.weight))
        push!(seen, fill(false, length(ps.weight)))
    end
    tail = Base.checked_add(total, 1)
    length(bins.offsets) == length(mesh.widths) + 1 &&
        length(bins.species) == length(bins.particle) == total &&
        bins.offsets[1] == 1 &&
        bins.offsets[end] == tail ||
        throw(ArgumentError("composite ownership bins have inconsistent sizes"))
    for cell in eachindex(mesh.widths)
        lo, stop = bins.offsets[cell], bins.offsets[cell+1]
        1 <= lo <= stop <= tail || throw(ArgumentError("invalid composite ownership offsets"))
        for slot = lo:stop-1
            s, p = bins.species[slot], bins.particle[slot]
            1 <= s <= length(species) && 1 <= p <= length(species[s].weight) ||
                throw(ArgumentError("invalid composite ownership particle index"))
            !seen[s][p] && _static_amr_leaf(mesh, _exact_rational(species[s].x[1][p])) == cell ||
                throw(
                    ArgumentError(
                        "composite ownership must contain each particle once in its actual leaf",
                    ),
                )
            seen[s][p] = true
        end
    end
    all(a -> all(a), seen) || throw(ArgumentError("composite ownership omits particles"))
    return nothing
end
