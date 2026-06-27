# migration.jl (§5.3 particles/)
#
# Transport-agnostic particle migration across logical ranks. The MPI layer only
# needs to ship the outbound subsets produced by the same destination rules; this
# serial implementation is the deterministic reference used by tests.

function _subset_particles(ps::ParticleSet{D,T}, idx::AbstractVector{<:Integer}) where {D,T}
    out = ParticleSet{D,T}(length(idx); q = ps.q, m = ps.m)
    for d = 1:D
        out.x[d] .= ps.x[d][idx]
    end
    for c = 1:3
        out.v[c] .= ps.v[c][idx]
    end
    out.weight .= ps.weight[idx]
    out.id .= ps.id[idx]
    out.tag .= ps.tag[idx]
    return out
end

function _replace_particles!(ps::ParticleSet{D,T}, idx::AbstractVector{<:Integer}) where {D,T}
    ps.x = ntuple(d -> ps.x[d][idx], D)
    ps.v = ntuple(c -> ps.v[c][idx], 3)
    ps.weight = ps.weight[idx]
    ps.id = ps.id[idx]
    ps.tag = ps.tag[idx]
    return ps
end

"""
    append_particles!(dest, src) -> dest

Append every particle in `src` to `dest`, preserving all per-particle fields.
Species charge and mass must match exactly.
"""
function append_particles!(dest::ParticleSet{D,T}, src::ParticleSet{D,T}) where {D,T}
    dest.q == src.q || throw(ArgumentError("cannot append particles with different charge"))
    dest.m == src.m || throw(ArgumentError("cannot append particles with different mass"))
    nparticles(src) == 0 && return dest
    dest.x = ntuple(d -> vcat(dest.x[d], src.x[d]), D)
    dest.v = ntuple(c -> vcat(dest.v[c], src.v[c]), 3)
    dest.weight = vcat(dest.weight, src.weight)
    dest.id = vcat(dest.id, src.id)
    dest.tag = vcat(dest.tag, src.tag)
    return dest
end

function _wrap_for_layout!(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    layout::LogicalRankLayout{D},
) where {D,T}
    @inbounds for d = 1:D
        if layout.periodic[d]
            Ld = g.L[d]
            for p = 1:nparticles(ps)
                ps.x[d][p] = mod(ps.x[d][p], Ld)
            end
        end
    end
    return ps
end

"""
    migrate_particles!(rank_particles, g, layout) -> (; moved, lost)

Move particles between a vector of rank-local [`ParticleSet`](@ref)s according to
`layout`. Periodic axes are wrapped into the global domain before destination
classification. Particles outside a nonperiodic axis are removed and counted as
`lost`.

This routine is deliberately local and deterministic; an MPI implementation can
use the same classification and subset rules before replacing the in-process
append with send/receive.
"""
function migrate_particles!(
    rank_particles::AbstractVector{<:ParticleSet{D,T}},
    g::FourierGrid{D,T},
    layout::LogicalRankLayout{D},
) where {D,T}
    length(rank_particles) == nranks(layout) ||
        throw(ArgumentError("rank_particles length must equal nranks(layout)"))

    moves = Tuple{Int,ParticleSet{D,T}}[]
    moved = 0
    lost = 0

    for r = 1:length(rank_particles)
        ps = rank_particles[r]
        _wrap_for_layout!(ps, g, layout)
        keep = Int[]
        lost_idx = Int[]
        by_dest = Dict{Int,Vector{Int}}()

        @inbounds for p = 1:nparticles(ps)
            pos = ntuple(d -> ps.x[d][p], D)
            dest = rank_of_position(pos, g, layout)
            if dest === nothing
                push!(lost_idx, p)
            elseif dest == r
                push!(keep, p)
            else
                push!(get!(by_dest, dest, Int[]), p)
            end
        end

        for (dest, idx) in by_dest
            push!(moves, (dest, _subset_particles(ps, idx)))
            moved += length(idx)
        end
        lost += length(lost_idx)
        _replace_particles!(ps, keep)
    end

    for (dest, payload) in moves
        append_particles!(rank_particles[dest], payload)
    end

    return (; moved, lost)
end
