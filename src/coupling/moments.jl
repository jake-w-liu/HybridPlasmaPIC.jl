# Moments.jl — density/momentum/current/pressure moments (from deposit.jl)

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

"Momentum density (n u)_c = (Σ_p w_p v_{c,p} S_g)/ΔV for c = 1,2,3."
function momentum!(
    mom::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    work::AbstractVector{T} = Vector{T}(undef, nparticles(ps)),
) where {D,T}
    length(work) == nparticles(ps) || throw(
        DimensionMismatch(
            "work length $(length(work)) must equal particle count $(nparticles(ps))",
        ),
    )
    ΔV = prod(g.dx)
    for c = 1:3
        @. work = ps.weight * ps.v[c]
        deposit_scalar!(mom[c], ps, work, g, shape)
        mom[c] ./= ΔV
    end
    return mom
end

"Ion current density J_c = q · (Σ_p w_p v_{c,p} S_g)/ΔV."
function current!(
    J::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction;
    work::AbstractVector{T} = Vector{T}(undef, nparticles(ps)),
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
    work::AbstractVector{T} = Vector{T}(undef, nparticles(ps)),
    nbuf::Array{T,D} = Array{T,D}(undef, g.n),
    mom::NTuple{3,<:Array{T,D}} = ntuple(_ -> Array{T,D}(undef, g.n), 3),
) where {D,T}
    Np = nparticles(ps)
    length(work) == Np ||
        throw(DimensionMismatch("work length $(length(work)) must equal particle count $Np"))
    size(nbuf) == g.n ||
        throw(DimensionMismatch("nbuf size $(size(nbuf)) does not match grid size $(g.n)"))
    for c = 1:3
        size(mom[c]) == g.n ||
            throw(DimensionMismatch("mom[$c] size $(size(mom[c])) does not match grid size $(g.n)"))
    end
    for c = 1:6
        size(P[c]) == g.n ||
            throw(DimensionMismatch("P[$c] size $(size(P[c])) does not match grid size $(g.n)"))
    end
    density!(nbuf, ps, g, shape)
    momentum!(mom, ps, g, shape; work)         # mom = n·U
    ΔV = prod(g.dx)
    nf = T(nfloor)
    mq = ps.m
    for (idx, (i, j)) in enumerate(_PT_PAIRS)
        vi = ps.v[i]
        vj = ps.v[j]
        @inbounds @. work = ps.weight * vi * vj
        deposit_scalar!(P[idx], ps, work, g, shape)
        Pij = P[idx]
        mi = mom[i]
        mj = mom[j]
        @inbounds for I in eachindex(Pij)
            nv = max(nbuf[I], nf)
            Pij[I] = mq * (Pij[I] / ΔV - mi[I] * mj[I] / nv)   # Π_ij/ΔV − ρ U_i U_j
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
