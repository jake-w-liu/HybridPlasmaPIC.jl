_tah_particle_arrays(p) = (p.x..., p.v..., p.weight, p.id, p.tag)

function _tah_particle_storage(p, ::Type{T}) where {T}
    p isa _ElectrostaticHostParticle{3,T} ||
        throw(ArgumentError("annular ions require grid-precision ordinary 3D host particle leaves"))
    _check_particle_storage_axes(p)
    _require_particle_ordinary_stores(p, "toroidal annulus")
    return nothing
end

function _tah_independent(arrays)
    for i in eachindex(arrays), j = 1:(i-1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("annular public/private/source/geometry roles overlap"))
    end
    return nothing
end

function _tah_clone_species(species, ::Type{T}) where {T}
    copies = _ElectrostaticHostParticle{3,T}[]
    _pic_prepare_species_stages!(copies, species)
    return copies
end

function _tah_particle_values(p, mesh)
    isfinite(p.q) &&
        p.q > 0 &&
        isfinite(p.m) &&
        p.m > 0 &&
        all(w -> isfinite(w) && w >= 0, p.weight) &&
        all(a -> all(isfinite, a), (p.x..., p.v...)) || throw(
        ArgumentError(
            "annular particles require finite phase space, positive q/m and nonnegative weights",
        ),
    )
    all(x -> 0 <= x <= mesh.L[1], p.x[1]) ||
        throw(ArgumentError("annular particle stage leaves its admitted radial interval"))
    return nothing
end

function _particle_cell_position(
    p::ParticleSet{3,T},
    mesh::ToroidalAnnulusGrid{T},
    d::Int,
    marker::Int,
) where {T}
    x = p.x[d][marker]
    isfinite(x) && (d != 1 || 0 <= x <= mesh.L[1]) ||
        throw(ArgumentError("annular CIC coordinate is outside its logical domain"))
    coordinate = d == 1 ? x / mesh.dx[d] : _periodic_cell_position(x, mesh.L[d], mesh.dx[d])
    isfinite(coordinate) &&
        _exact_rational(typemin(Int)) <=
        _exact_rational(coordinate) <=
        _exact_rational(typemax(Int) - 2) ||
        throw(ArgumentError("annular CIC coordinate is not safely Int-convertible"))
    return coordinate
end

function _tah_moments!(st, species)
    sums, f, mesh = st.work.moments, st.work.fields, st.mesh
    foreach(a -> fill!(a, zero(_CompositeQ)), sums)
    Q = _exact_rational
    for p in species
        _foreach_exact_particle_contribution(
            p,
            (p.weight,),
            mesh,
            CIC(),
            Val(:fold_x),
        ) do index, value, marker
            charge = Q(p.q) * value
            sums[1][index...] += charge
            for c = 1:3
                sums[c+1][index...] += charge * Q(p.v[c][marker])
            end
        end
    end
    T = eltype(f.n)
    for I in CartesianIndices(f.n)
        volume = _tah_volume(mesh, I)
        f.n[I] =
            _tah_round(T, sums[1][I] / volume, "annular actual density"; preserve_nonzero = true)
        f.n[I] > 0 ||
            throw(ArgumentError("annular actual density must be positive; no floor is supplied"))
        f.ninv[I] =
            _tah_round(T, inv(Q(f.n[I])), "annular reciprocal density"; preserve_nonzero = true)
        for c = 1:3
            f.Ji[c][I] = _tah_round(T, sums[c+1][I] / volume, "annular ion current")
            f.ui[c][I] = _tah_round(T, Q(f.Ji[c][I]) / Q(f.n[I]), "annular ion velocity")
        end
    end
    return nothing
end

function _tah_gather(st, p::ParticleSet{3,T}, marker) where {T}
    mesh, f = st.mesh, st.work.fields
    stencil = ntuple(d -> _stencil1d(CIC(), _particle_cell_position(p, mesh, d, marker)), 3)
    stamp = CartesianIndices((2, 2, 2))
    E = _gather_vector_stencil(stencil, mesh.n, stamp, T, Val(:fold_x), Val(true)) do index
        ntuple(c -> f.E[c][index...], 3)
    end
    B = _gather_vector_stencil(stencil, mesh.n, stamp, T, Val(:fold_x), Val(true)) do index
        ntuple(c -> f.B[c][index...], 3)
    end
    return E, B
end
