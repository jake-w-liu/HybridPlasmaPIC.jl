# Private nodal companions of the existing closed-cloud edge primitives.
# Charge keeps endpoint tails; transverse PEC degrees of freedom exclude them.
@inline function _yee_cloud_node(x, node, h, a)
    left, _ = _yee_cloud_edge(x, (node - 1) * h, h, a)
    right, _ = _yee_cloud_edge(x, node * h, h, a)
    return left - right
end
@inline function _yee_cloud_node_path(x, y, node, h, a)
    x == y && return _yee_cloud_node(x, node, h, a)
    primitive(t) =
        _yee_cloud_integrated_edge(t, (node - 1) * h, h, a) -
        _yee_cloud_integrated_edge(t, node * h, h, a)
    return (primitive(y) - primitive(x)) / (y - x)
end
@inline function _yee_cloud_node_span(x, y, h, a, n)
    edges = _yee_cloud_edge_span(x, y, h, a, n)
    return max(1, first(edges)):min(n - 1, last(edges) + 1)
end

# Caller validates the complete field once. PEC endpoint E is stored zero;
# G0 needs no delta/h factor and no wall renormalization.
function _yee_cloud_nodal_rational(E, x, g::YeeGrid{1}, shape::YeeTriangularCloud)
    h = _exact_rational(g.L[1]) / g.n[1]
    a, position = _exact_rational(shape.half_width), _exact_rational(x)
    total = zero(h)
    for node in _yee_cloud_node_span(position, position, h, a, g.n[1])
        total += _exact_rational(E[node+1]) * _yee_cloud_node(position, node, h, a)
    end
    return total
end

# Complete original signed current entries. No per-particle grid matrices or
# physical partial currents are published; final PEC endpoints remain zero.
function _yee_transverse_impulse!(out, species, candidates, g::YeeGrid{1,T}, shape, dt::T) where {T}
    h, delta = _exact_rational(g.L[1]) / g.n[1], _exact_rational(g.dx[1])
    a, tau = _exact_rational(shape.half_width), _exact_rational(dt)
    sums = ntuple(_ -> zeros(Rational{BigInt}, g.n[1] + 1), 2)
    for (ps, candidate) in zip(species, candidates), p in eachindex(ps.weight)
        charge = _exact_rational(ps.q) * _exact_rational(ps.weight[p])
        iszero(charge) && continue
        x, y = _exact_rational(ps.x[1][p]), _exact_rational(candidate.x[1][p])
        factors = ntuple(c -> charge * tau * _exact_rational(candidate.v[c+1][p]) / delta, 2)
        for node in _yee_cloud_node_span(x, y, h, a, g.n[1])
            coefficient = _yee_cloud_node_path(x, y, node, h, a)
            for c = 1:2
                sums[c][node+1] += factors[c] * coefficient
            end
        end
    end
    for c = 1:2, i in eachindex(out[c+1])
        value = _from_exact(T, sums[c][i])
        isfinite(value) || throw(ArgumentError("Yee transverse impulse does not round finite"))
        out[c+1][i] = value
    end
    return nothing
end
