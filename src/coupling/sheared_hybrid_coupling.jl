# Exact collection moments and complete gather use the existing canonical TSC
# and Hall iterator. Positions here are logical; velocity/current are Cartesian.
function _sheared_moments!(st, species)
    Q = _exact_rational
    sums, f = st.work.moments, st.work.fields
    foreach(a -> fill!(a, zero(Rational{BigInt})), sums)
    for p in species
        charge = Q(p.q)
        _foreach_exact_particle_contribution(p, (p.weight,), st.g, TSC(), Val(:periodic)) do index, value, marker
            contribution = charge * value
            sums[1][index...] += contribution
            for c in 1:3
                sums[c + 1][index...] += contribution * Q(p.v[c][marker])
            end
        end
    end
    T = eltype(f.n)
    volume = prod(Q, st.g.dx)
    for I in eachindex(f.n)
        f.n[I] = _shear_round(T, sums[1][I] / volume, "sheared charge density"; preserve_nonzero = true)
        f.n[I] > 0 || throw(ArgumentError("sheared actual density must be strictly positive; no floor is supplied"))
        f.ninv[I] = _shear_round(T, inv(Q(f.n[I])), "sheared reciprocal density"; preserve_nonzero = true)
        for c in 1:3
            f.Ji[c][I] = _shear_round(T, sums[c + 1][I] / volume, "sheared ion current")
            f.ui[c][I] = _shear_round(T, Q(f.Ji[c][I]) / Q(f.n[I]), "sheared ion velocity")
        end
    end
    return nothing
end

function _sheared_gather(st, p::ParticleSet{3,T}, marker) where {T}
    stencil = ntuple(d -> _stencil1d(TSC(), _particle_cell_position(p, st.g, d, marker)), 3)
    stamp = CartesianIndices((3, 3, 3))
    E = _gather_vector_stencil(stencil, st.g.n, stamp, T, Val(:periodic), Val(true)) do index
        ntuple(c -> st.work.fields.E[c][index...], 3)
    end
    B = _gather_vector_stencil(stencil, st.g.n, stamp, T, Val(:periodic), Val(true)) do index
        ntuple(c -> st.work.fields.B[c][index...], 3)
    end
    return E, B
end
