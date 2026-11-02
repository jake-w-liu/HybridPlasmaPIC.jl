# Gather.jl — field gather (from deposit.jl)

"""
    gather_scalar!(out, field, ps, g, shape)

Interpolate `out[p] = Σ_g field_g · S_g(x_p)` — the transpose of
[`deposit_scalar!`](@ref) with the same shape.
"""
function gather_scalar!(
    out::AbstractVector{T},
    field::Array{T,D},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    length(out) == nparticles(ps) || throw(DimensionMismatch("out length ≠ particle count"))
    size(field) == g.n ||
        throw(DimensionMismatch("field size $(size(field)) does not match grid size $(g.n)"))
    _validate_particle_positions(ps)
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for (p, oi) in enumerate(eachindex(out))
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        acc = zero(T)
        for c in stamp
            o = Tuple(c)
            w = one(T)
            for d = 1:D
                w *= st[d][2][o[d]]
            end
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
            acc += field[idx...] * w
        end
        out[oi] = acc
    end
    return out
end

# One canonical vector gather row; distributed ownership supplies only the lookup.
@inline function _gather_vector_stencil(
    field_at,
    st,
    n,
    stamp,
    ::Type{T},
    boundary = Val(:periodic),
    complete = Val(false),
    ops::_ScalarOps = _DefaultScalarOps(),
) where {T}
    complete isa Val{true} &&
        return _gather_vector_stencil_exact(field_at, st, n, stamp, T, boundary)
    ax = ay = az = zero(T)
    for c in stamp
        o = Tuple(c)
        w = one(T)
        for d = 1:length(n)
            w = _scalar_mul(ops, w, st[d][2][o[d]])
        end
        idx = _stencil_node_index(st, n, o, boundary)
        values = field_at(idx)
        ax = _scalar_add(ops, ax, _scalar_mul(ops, values[1], w))
        ay = _scalar_add(ops, ay, _scalar_mul(ops, values[2], w))
        az = _scalar_add(ops, az, _scalar_mul(ops, values[3], w))
    end
    return ax, ay, az
end

# Complete the same stored per-axis stencil, with one publication per component.
function _gather_vector_stencil_exact(field_at, st, n, stamp, ::Type{T}, boundary) where {T}
    Q = _exact_rational
    sums = (zero(Rational{BigInt}), zero(Rational{BigInt}), zero(Rational{BigInt}))
    for c in stamp
        o = Tuple(c)
        weight = one(Rational{BigInt})
        for d = 1:length(n)
            weight *= Q(st[d][2][o[d]])
        end
        values = field_at(_stencil_node_index(st, n, o, boundary))
        sums = ntuple(i -> sums[i] + Q(values[i]) * weight, 3)
    end
    out = map(x -> _from_exact(T, x), sums)
    all(isfinite, out) ||
        throw(ArgumentError("complete vector gather exceeds its materialized field range"))
    return out
end

"Gather a 3-component grid field to particles using a shared interpolation stencil."
function gather_vector!(
    out::NTuple{3,<:AbstractVector{T}},
    field::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    for c = 1:3
        length(out[c]) == nparticles(ps) ||
            throw(DimensionMismatch("out[$c] length ≠ particle count"))
        size(field[c]) == g.n || throw(
            DimensionMismatch("field[$c] size $(size(field[c])) does not match grid size $(g.n)"),
        )
    end
    _validate_particle_positions(ps)
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    output_indices = zip(eachindex(out[1]), eachindex(out[2]), eachindex(out[3]))
    @inbounds for (p, oi) in enumerate(output_indices)
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        ax, ay, az = _gather_vector_stencil(st, n, stamp, T) do idx
            (field[1][idx...], field[2][idx...], field[3][idx...])
        end
        out[1][oi[1]] = ax
        out[2][oi[2]] = ay
        out[3][oi[3]] = az
    end
    return out
end

"Number density n_g = (Σ_p w_p S_g)/ΔV."
