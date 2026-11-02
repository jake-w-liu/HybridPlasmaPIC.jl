# Periodization of the canonical free-space cloud, only for the compatible
# periodic implicit model. Closed PEC public coupling keeps its own policy.
function _pyee_span(x, y, h, a, L, n, nodal::Bool)
    radius = nodal ? a + h : a
    max(x, y) + radius <= 0 && return 1:0
    min(x, y) - a >= L && return 1:0
    lo = floor(Int, clamp(min(x, y) - a, zero(h), L) / h)
    hi = floor(Int, clamp(max(x, y) + a, zero(h), L) / h)
    nodal && hi < n - 1 && (hi += 1)
    return min(lo, n - 1):min(hi, n - 1)
end
function _pyee_segments(x, y, L)
    z = zero(L)
    y < z && return ((x, z), (L, y + L))
    y > L && return ((x, L), (z, y - L))
    return ((x, y), (z, z))
end
function _pyee_each_basis(f, x, y, g::YeeGrid{1}, a, ::Val{nodal}) where {nodal}
    L, h = _exact_rational(g.L[1]), _exact_rational(g.dx[1])
    d = y - x
    if iszero(d)
        for image = -1:1
            position = x + image * L
            for i in _pyee_span(position, position, h, a, L, g.n[1], nodal)
                coefficient =
                    nodal ? _yee_cloud_node(position, i, h, a) :
                    last(_yee_cloud_edge(position, i * h, h, a))
                iszero(coefficient) || f(i + 1, coefficient)
            end
        end
    else
        for (left, right) in _pyee_segments(x, y, L)
            left == right && continue
            for image = -1:1
                before, after = left + image * L, right + image * L
                for i in _pyee_span(before, after, h, a, L, g.n[1], nodal)
                    coefficient = if nodal
                        (right - left) * _yee_cloud_node_path(before, after, i, h, a) / d
                    else
                        (
                            first(_yee_cloud_edge(after, i * h, h, a)) -
                            first(_yee_cloud_edge(before, i * h, h, a))
                        ) / d
                    end
                    iszero(coefficient) || f(i + 1, coefficient)
                end
            end
        end
    end
    return nothing
end
function _pyee_charge!(out::Vector{Rational{BigInt}}, species, g, shape)
    fill!(out, zero(Rational{BigInt}))
    a = _exact_rational(shape.half_width)
    for ps in species, p in eachindex(ps.weight)
        charge = _exact_rational(ps.q) * _exact_rational(ps.weight[p])
        iszero(charge) && continue
        x = _exact_rational(ps.x[1][p])
        _pyee_each_basis(x, x, g, a, Val(true)) do i, coefficient
            out[i] += charge * coefficient
        end
    end
    return out
end

# Both deposit and gather consume each identical original-Q coefficient in this
# one pass; no rounded sample or per-particle/grid coefficient matrix intervenes.
function _pyee_couple!(w, species, g, shape)
    Q = _exact_rational
    h, a = Q(g.dx[1]), Q(shape.half_width)
    foreach(a -> fill!(a, zero(Rational{BigInt})), (w.exact_current..., w.Ep..., w.Bp...))
    fill!(w.kinematic_bound, zero(Rational{BigInt}))
    offset = 0
    for ps in species
        for p in eachindex(ps.weight)
            index = offset + p
            x, y = Q(ps.x[1][p]), w.lift[index]
            charge = Q(ps.q) * Q(ps.weight[p])
            velocity = ntuple(c -> w.average_velocity[c][index], 3)
            bound = abs(charge * w.path_error[index])
            _pyee_each_basis(x, y, g, a, Val(false)) do i, coefficient
                w.exact_current[1][i] += charge * velocity[1] * coefficient
                w.Ep[1][index] += h * w.mean_E[1][i] * coefficient
                w.Bp[2][index] += h * w.mean_B[2][i] * coefficient
                w.Bp[3][index] += h * w.mean_B[3][i] * coefficient
                w.kinematic_bound[i] += bound * coefficient
                w.kinematic_bound[i == g.n[1] ? 1 : i + 1] += bound * coefficient
            end
            _pyee_each_basis(x, y, g, a, Val(true)) do i, coefficient
                for c = 2:3
                    w.exact_current[c][i] += charge * velocity[c] * coefficient / h
                    w.Ep[c][index] += w.mean_E[c][i] * coefficient
                end
                w.Bp[1][index] += w.mean_B[1][i] * coefficient
            end
        end
        offset += nparticles(ps)
    end
    return nothing
end
