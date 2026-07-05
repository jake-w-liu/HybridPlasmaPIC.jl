# Deposit.jl — scalar deposition (from deposit.jl)

"""
    deposit_scalar!(out, ps, vals, g, shape)

Accumulate `out_g = Σ_p vals[p] · S_g(x_p)` for the given shape (periodic wrap).
`out` is zeroed first. Partition of unity ⇒ `Σ_g out_g = Σ_p vals[p]`.
"""
function deposit_scalar!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    vals::AbstractVector,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    length(vals) == nparticles(ps) || throw(DimensionMismatch("vals length ≠ particle count"))
    size(out) == g.n ||
        throw(DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"))
    _validate_particle_positions(ps)
    fill!(out, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for (p, vi) in enumerate(eachindex(vals))
        st = ntuple(d -> _stencil1d(shape, _particle_cell_position(ps, g, d, p)), D)
        val = vals[vi]
        for c in stamp
            o = Tuple(c)
            w = val
            for d = 1:D
                w *= st[d][2][o[d]]
            end
            idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
            out[idx...] += w
        end
    end
    return out
end
