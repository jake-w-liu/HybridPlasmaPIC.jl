# provenance.jl (§5.3 particles/)
#
# Deterministic per-rank random streams and rank-independent particle IDs.
# Run-level provenance/repro is in io/restart.jl (archive_run/operators_match).

const _PARTICLE_ID_SPECIES_BITS = 16
const _PARTICLE_ID_INDEX_BITS = 64 - _PARTICLE_ID_SPECIES_BITS
const _PARTICLE_ID_MAX_SPECIES = (UInt64(1) << _PARTICLE_ID_SPECIES_BITS) - UInt64(1)
const _PARTICLE_ID_MAX_INDEX = (UInt64(1) << _PARTICLE_ID_INDEX_BITS) - UInt64(1)

@inline function _nonnegative_uint64(name::AbstractString, x::Integer)
    x >= 0 || throw(ArgumentError("$name must be nonnegative, got $x"))
    x <= typemax(UInt64) || throw(ArgumentError("$name must fit in UInt64, got $x"))
    return UInt64(x)
end

@inline function _positive_particle_index(x::Integer)
    x >= 1 || throw(ArgumentError("global particle index must be >= 1, got $x"))
    x <= _PARTICLE_ID_MAX_INDEX ||
        throw(ArgumentError("global particle index must be <= $_PARTICLE_ID_MAX_INDEX, got $x"))
    return UInt64(x)
end

@inline function _particle_species(x::Integer)
    x >= 0 || throw(ArgumentError("species must be nonnegative, got $x"))
    x <= _PARTICLE_ID_MAX_SPECIES ||
        throw(ArgumentError("species must be <= $_PARTICLE_ID_MAX_SPECIES, got $x"))
    return UInt64(x)
end

@inline function _mix64(x::UInt64)
    x = xor(x, x >> 30) * 0xbf58476d1ce4e5b9
    x = xor(x, x >> 27) * 0x94d049bb133111eb
    return xor(x, x >> 31)
end

@inline _combine_seed(seed::UInt64, value::UInt64) = _mix64(xor(seed + 0x9e3779b97f4a7c15, value))

"""
    rank_seed(base_seed, layout, rank; stream=0) -> UInt64

Return a deterministic seed for `rank` in a [`LogicalRankLayout`](@ref). The
seed is derived from the base seed, layout dimensions, rank coordinates,
periodicity flags, and optional nonnegative stream number.
"""
function rank_seed(
    base_seed::Integer,
    layout::LogicalRankLayout{D},
    rank::Integer;
    stream::Integer = 0,
) where {D}
    coords = rank_coords(layout, rank)
    seed = _mix64(_nonnegative_uint64("base_seed", base_seed))
    seed = _combine_seed(seed, UInt64(D))
    seed = _combine_seed(seed, _nonnegative_uint64("stream", stream))
    for d = 1:D
        seed = _combine_seed(seed, UInt64(d))
        seed = _combine_seed(seed, UInt64(layout.ranks[d]))
        seed = _combine_seed(seed, UInt64(coords[d]))
        seed = _combine_seed(seed, layout.periodic[d] ? UInt64(1) : UInt64(0))
    end
    return seed
end

"""
    rank_rng(base_seed, layout, rank; stream=0) -> MersenneTwister

Create a reproducible rank-local RNG initialized by [`rank_seed`](@ref).
"""
function rank_rng(base_seed::Integer, layout::LogicalRankLayout, rank::Integer; stream::Integer = 0)
    return MersenneTwister(rank_seed(base_seed, layout, rank; stream))
end

@inline function _global_particle_id(index::UInt64, species::UInt64)
    return (species << _PARTICLE_ID_INDEX_BITS) | index
end

"""
    global_particle_id(global_index; species=0) -> UInt64

Return a rank-independent 64-bit particle ID. IDs are encoded as a 16-bit
species namespace and a 48-bit one-based global particle index, so the same
physical particle receives the same ID regardless of rank decomposition.
"""
function global_particle_id(global_index::Integer; species::Integer = 0)
    return _global_particle_id(_positive_particle_index(global_index), _particle_species(species))
end

"""
    assign_global_particle_ids!(ps, first_global_index=1; species=0) -> ps

Assign consecutive rank-independent IDs to every particle in `ps`, starting at
`first_global_index`.
"""
function assign_global_particle_ids!(
    ps::ParticleSet,
    first_global_index::Integer = 1;
    species::Integer = 0,
)
    start = _positive_particle_index(first_global_index)
    sp = _particle_species(species)
    n = nparticles(ps)
    if n > 0
        nminus1 = UInt64(n - 1)
        nminus1 <= _PARTICLE_ID_MAX_INDEX - start ||
            throw(ArgumentError("particle ID range exceeds $_PARTICLE_ID_MAX_INDEX"))
        @inbounds for p = 1:n
            ps.id[p] = _global_particle_id(start + UInt64(p - 1), sp)
        end
    end
    return ps
end

"""
    assign_global_particle_ids!(ps, global_indices; species=0) -> ps

Assign IDs from explicit one-based global particle indices. This form is useful
when local particles are a non-contiguous subset of a global deterministic load.
"""
function assign_global_particle_ids!(
    ps::ParticleSet,
    global_indices::AbstractVector{<:Integer};
    species::Integer = 0,
)
    length(global_indices) == nparticles(ps) ||
        throw(DimensionMismatch("global_indices length must equal nparticles(ps)"))
    sp = _particle_species(species)
    @inbounds for p = 1:nparticles(ps)
        ps.id[p] = _global_particle_id(_positive_particle_index(global_indices[p]), sp)
    end
    return ps
end
