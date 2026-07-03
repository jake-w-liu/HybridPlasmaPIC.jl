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

"Gather a 3-component grid field to particles, component by component."
function gather_vector!(
    out::NTuple{3,<:AbstractVector{T}},
    field::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    for c = 1:3
        gather_scalar!(out[c], field[c], ps, g, shape)
    end
    return out
end

"Number density n_g = (Σ_p w_p S_g)/ΔV."
