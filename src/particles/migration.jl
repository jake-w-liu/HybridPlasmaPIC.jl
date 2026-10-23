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

# Construct the exact destination type before publishing any field. Indexing and
# field conversion can fail for custom storage, including a late id/tag field.
function _prepare_particle_replacement(
    ps::P,
    idx::AbstractVector{<:Integer},
) where {D,T,P<:ParticleSet{D,T}}
    prepared = P(
        ntuple(d -> ps.x[d][idx], D),
        ntuple(c -> ps.v[c][idx], 3),
        ps.weight[idx],
        ps.id[idx],
        ps.tag[idx],
        ps.q,
        ps.m,
    )
    _check_particle_storage_axes(prepared)
    _particle_sets_mightalias(prepared, ps, _particle_backing_mightalias) &&
        throw(ArgumentError("migration indexed storage must be independent of the source"))
    _particle_components_mightalias(prepared, _particle_backing_mightalias) &&
        throw(ArgumentError("migration indexed components must have independent storage"))
    return prepared
end

function _publish_particle_replacement!(ps::P, prepared::P) where {P<:ParticleSet}
    ps.x = prepared.x
    ps.v = prepared.v
    ps.weight = prepared.weight
    ps.id = prepared.id
    ps.tag = prepared.tag
    return ps
end

"""
    append_particles!(dest, src) -> dest

Append every particle in `src` to `dest`, preserving all per-particle fields.
Species charge and mass must match exactly. For a nonempty source, destination
components must have independent resizable storage. Distinct populations must
not share backing storage; self-append is supported. Detectable storage input
errors are rejected before either population changes.
"""
function append_particles!(dest::ParticleSet{D,T}, src::ParticleSet{D,T}) where {D,T}
    dest.q == src.q || throw(ArgumentError("cannot append particles with different charge"))
    dest.m == src.m || throw(ArgumentError("cannot append particles with different mass"))
    nparticles(src) == 0 && return dest
    count = _check_particle_storage_axes(src)
    _validate_particle_growth(dest, count, "append_particles!: destination")
    dest === src ||
        !_particle_sets_mightalias(dest, src, _particle_backing_mightalias) ||
        throw(
            ArgumentError("append_particles!: distinct populations must not share backing storage"),
        )
    # In-place append! (not vcat): vcat reallocates and copies every destination array on
    # each call, so repeated appends (ionization sources, per-step MPI receives) cost O(N²)
    # over a run; append! amortizes to O(total). Requires resizable (Vector) fields, which
    # the ParticleSet constructor produces; it preserves the one-based-axes contract.
    for d = 1:D
        append!(dest.x[d], src.x[d])
    end
    for c = 1:3
        append!(dest.v[c], src.v[c])
    end
    append!(dest.weight, src.weight)
    append!(dest.id, src.id)
    append!(dest.tag, src.tag)
    return dest
end

function _validate_rank_particle_species!(
    rank_particles::AbstractVector{<:ParticleSet{D,T}},
) where {D,T}
    isempty(rank_particles) && return nothing
    ref_q = first(rank_particles).q
    ref_m = first(rank_particles).m
    for (r, ps) in pairs(rank_particles)
        ps.q == ref_q || throw(
            ArgumentError("rank $r particle charge $(ps.q) differs from rank 1 charge $(ref_q)"),
        )
        ps.m == ref_m ||
            throw(ArgumentError("rank $r particle mass $(ps.m) differs from rank 1 mass $(ref_m)"))
    end
    return nothing
end

function _has_nonfinite_particle_positions(ps::ParticleSet{D}) where {D}
    @inbounds for d = 1:D, p = 1:nparticles(ps)
        isfinite(ps.x[d][p]) || return true
    end
    return false
end

function _validate_migration_positions!(ps::ParticleSet{D}) where {D}
    @inbounds for d = 1:D, p = 1:nparticles(ps)
        isfinite(ps.x[d][p]) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
    end
    return nothing
end

function _validate_particle_migration_inputs!(
    rank_particles::AbstractVector{<:ParticleSet{D,T}},
) where {D,T}
    Base.require_one_based_indexing(rank_particles)
    owners = Base.IdSet{Any}()
    for ps in rank_particles
        ps in owners && throw(ArgumentError("migration ranks must have distinct particle owners"))
        push!(owners, ps)
        _check_particle_storage_axes(ps)
    end
    _validate_rank_particle_species!(rank_particles)
    for ps in rank_particles
        _validate_migration_positions!(ps)
    end
    return nothing
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
                ps.x[d][p] = _periodic_coordinate(ps.x[d][p], zero(T), Ld, Ld)
            end
        end
    end
    return ps
end

function _partition_particle_migration(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    layout::LogicalRankLayout{D},
    rank::Int,
) where {D,T}
    keep = Int[]
    by_dest = Dict{Int,Vector{Int}}()
    lost = 0
    @inbounds for p = 1:nparticles(ps)
        # Rank lookup owns periodic wrapping; classification does not write input storage.
        pos = ntuple(d -> ps.x[d][p], D)
        dest = rank_of_position(pos, g, layout)
        if dest === nothing
            lost += 1
        elseif dest == rank
            push!(keep, p)
        else
            push!(get!(by_dest, dest, Int[]), p)
        end
    end
    return (; keep, by_dest, lost)
end

"""
    migrate_particles!(rank_particles, g, layout) -> (; moved, lost)

Move particles between a vector of rank-local [`ParticleSet`](@ref)s according to
`layout`. Periodic axes are wrapped into the global domain before destination
classification. Particles outside a nonperiodic axis are removed and counted as
`lost`.

This routine is deliberately local and deterministic; an MPI implementation can
use the same classification and subset rules before replacing the in-process
append with send/receive. Distinct mutable population owners are required. All
results are prepared before publication, so a storage preparation or assembly
failure leaves every input population unchanged.
"""
function migrate_particles!(
    rank_particles::AbstractVector{<:ParticleSet{D,T}},
    g::FourierGrid{D,T},
    layout::LogicalRankLayout{D},
) where {D,T}
    length(rank_particles) == nranks(layout) ||
        throw(ArgumentError("rank_particles length must equal nranks(layout)"))
    _validate_particle_migration_inputs!(rank_particles)

    moves = Tuple{Int,ParticleSet{D,T}}[]
    prepared = Vector{eltype(rank_particles)}(undef, length(rank_particles))
    moved = 0
    lost = 0

    for r = 1:length(rank_particles)
        ps = rank_particles[r]
        partition = _partition_particle_migration(ps, g, layout, r)
        for (dest, idx) in partition.by_dest
            payload = _subset_particles(ps, idx)
            _wrap_for_layout!(payload, g, layout)
            push!(moves, (dest, payload))
            moved += length(idx)
        end
        lost += partition.lost
        prepared[r] = _prepare_particle_replacement(ps, partition.keep)
        _wrap_for_layout!(prepared[r], g, layout)
    end

    for (dest, payload) in moves
        append_particles!(prepared[dest], payload)
    end
    for r = 1:length(rank_particles)
        _publish_particle_replacement!(rank_particles[r], prepared[r])
    end

    return (; moved, lost)
end
