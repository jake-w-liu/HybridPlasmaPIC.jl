# Moments.jl — density/momentum/current/pressure moments (from deposit.jl)

@inline function _deposit_weighted_scalar!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    weight::AbstractVector{T},
    value::AbstractVector{T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    np = nparticles(ps)
    length(weight) == np || throw(
        DimensionMismatch("weight length $(length(weight)) must equal particle count $np"),
    )
    length(value) == np || throw(
        DimensionMismatch("value length $(length(value)) must equal particle count $np"),
    )
    size(out) == g.n || throw(
        DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"),
    )
    fill!(out, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for (p, wi) in enumerate(eachindex(weight))
        st = ntuple(d -> _stencil1d(shape, ps.x[d][p] / g.dx[d]), D)
        val = weight[wi] * value[wi]
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

@inline function _deposit_weighted_product!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    weight::AbstractVector{T},
    left::AbstractVector{T},
    right::AbstractVector{T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    np = nparticles(ps)
    length(weight) == np || throw(
        DimensionMismatch("weight length $(length(weight)) must equal particle count $np"),
    )
    length(left) == np || throw(
        DimensionMismatch("left length $(length(left)) must equal particle count $np"),
    )
    length(right) == np || throw(
        DimensionMismatch("right length $(length(right)) must equal particle count $np"),
    )
    size(out) == g.n || throw(
        DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"),
    )
    fill!(out, zero(T))
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))
    @inbounds for (p, wi) in enumerate(eachindex(weight))
        st = ntuple(d -> _stencil1d(shape, ps.x[d][p] / g.dx[d]), D)
        val = weight[wi] * left[wi] * right[wi]
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

function density!(
    nout::Array{T,D},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    deposit_scalar!(nout, ps, ps.weight, g, shape)
    nout ./= prod(g.dx)
    return nout
end

"Momentum density (n u)_c = (Σ_p w_p v_{c,p} S_g(x_p)) / ΔV for c = 1,2,3."
function momentum!(
    mom::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    work::Union{Nothing,AbstractVector{T}} = nothing,
) where {D,T}
    ΔV = prod(g.dx)
    if work === nothing
        for c = 1:3
            _deposit_weighted_scalar!(mom[c], ps, ps.weight, ps.v[c], g, shape)
            mom[c] ./= ΔV
        end
    else
        length(work) == nparticles(ps) || throw(
            DimensionMismatch(
                "work length $(length(work)) must equal particle count $(nparticles(ps))",
            ),
        )
        for c = 1:3
            @. work = ps.weight * ps.v[c]
            deposit_scalar!(mom[c], ps, work, g, shape)
            mom[c] ./= ΔV
        end
    end
    return mom
end

"Ion current density J_c = q · (Σ_p w_p v_{c,p} S_g(x_p))/ΔV."
function current!(
    J::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    work::Union{Nothing,AbstractVector{T}} = nothing,
) where {D,T}
    momentum!(J, ps, g, shape; work)
    for c = 1:3
        J[c] .*= ps.q
    end
    return J
end

# (i,j) index pairs of the 6 independent symmetric-tensor components, in order
# (xx, yy, zz, xy, xz, yz).
const _PT_PAIRS = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))

"""
    pressure_tensor!(P, ps, g, shape; nfloor=1e-6, work=..., nbuf=..., mom=...)

Deposit the ion pressure tensor `P_ij = m·Σ_p w_p (v_i−U_i)(v_j−U_j) S / ΔV`
(`U` the local bulk velocity), as 6 components in the order (xx,yy,zz,xy,xz,yz).
Computed as the centered second moment: second-moment density − ρ U_i U_j.
"""
function pressure_tensor!(
    P::NTuple{6,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    nfloor = 1e-6,
    work::Union{Nothing,AbstractVector{T}} = nothing,
    nbuf::Union{Nothing,Array{T,D}} = nothing,
    mom::Union{Nothing,NTuple{3,<:Array{T,D}}} = nothing,
) where {D,T}
    Np = nparticles(ps)
    nbuf === nothing && (nbuf = similar(P[1]))
    mom === nothing && (mom = ntuple(_ -> similar(P[1]), 3))
    for c = 1:6
        size(P[c]) == g.n || throw(
            DimensionMismatch("P[$c] size $(size(P[c])) does not match grid size $(g.n)"),
        )
    end

    size(nbuf) == g.n || throw(
        DimensionMismatch("nbuf size $(size(nbuf)) does not match grid size $(g.n)"),
    )
    for c = 1:3
        size(mom[c]) == g.n || throw(
            DimensionMismatch(
                "mom[$c] size $(size(mom[c])) does not match grid size $(g.n)",
            ),
        )
    end

    if work === nothing
        density!(nbuf, ps, g, shape)
        momentum!(mom, ps, g, shape)         # mom = n·U
    else
        length(work) == Np ||
            throw(
                DimensionMismatch("work length $(length(work)) must equal particle count $Np"),
            )
        density!(nbuf, ps, g, shape)
        momentum!(mom, ps, g, shape; work)
    end

    ΔV = prod(g.dx)
    nf = T(nfloor)
    mq = ps.m
    for (idx, (i, j)) in enumerate(_PT_PAIRS)
        vi = ps.v[i]
        vj = ps.v[j]
        if work === nothing
            _deposit_weighted_product!(P[idx], ps, ps.weight, vi, vj, g, shape)
        else
            @inbounds @. work = ps.weight * vi * vj
            deposit_scalar!(P[idx], ps, work, g, shape)
        end
        Pij = P[idx]
        mi = mom[i]
        mj = mom[j]
        @inbounds for I in eachindex(Pij)
            nv = max(nbuf[I], nf)
            Pij[I] = mq * (Pij[I] / ΔV - mi[I] * mj[I] / nv)
        end
    end
    return P
end

"""
    temperature_components(P, n; nfloor=1e-6)

Per-cell directional temperatures `(T_x, T_y, T_z) = (P_xx, P_yy, P_zz)/n` from a
pressure tensor produced by [`pressure_tensor!`](@ref).
"""
function temperature_components(P::NTuple{6,<:Array{T,D}}, n::Array{T,D}; nfloor = 1e-6) where {T,D}
    nf = T(nfloor)
    return ntuple(c -> P[c] ./ max.(n, nf), 3)
end
