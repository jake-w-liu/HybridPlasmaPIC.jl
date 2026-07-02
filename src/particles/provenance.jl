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
    # Validate every index BEFORE writing any id, so an invalid entry leaves ps.id unchanged
    # (strong exception guarantee), matching the range form above rather than corrupting the
    # id array half-way through.
    @inbounds for p = 1:nparticles(ps)
        _positive_particle_index(global_indices[p])
    end
    @inbounds for p = 1:nparticles(ps)
        ps.id[p] = _global_particle_id(_positive_particle_index(global_indices[p]), sp)
    end
    return ps
end

# §9.4 particle provenance — an ID-KEYED EVENT LOG (the checklist explicitly allows
# "event logs rather than per-particle dense arrays"). Keying by the stable global
# particle ID (not the local SoA index) makes the log invariant under sorting,
# compaction, and MPI migration, which a parallel dense array would not be. The
# species is already carried in the high bits of the ID, so the log + ID together
# supply the full §9.4 retain-list: id, species, source region, injection batch,
# injection time, first/last shock crossing, crossing count, reflection flag, and
# maximum kinetic energy.

@enum ProvenanceKind PROV_INJECTION = 1 PROV_CROSSING = 2 PROV_REFLECTION = 3 PROV_MAXKE = 4

struct ProvenanceEvent
    id::UInt64
    kind::ProvenanceKind
    time::Float64
    a::Float64            # injection: source_region; maxke: kinetic energy; else unused
    b::Float64            # injection: batch; else unused
end

"Append-only per-particle provenance event log keyed by global particle ID (§9.4)."
struct ParticleProvenanceLog
    events::Vector{ProvenanceEvent}
end
ParticleProvenanceLog() = ParticleProvenanceLog(ProvenanceEvent[])

"Record a particle injection: its `time`, `source_region`, and `batch`."
function record_injection!(
    log::ParticleProvenanceLog,
    id::Integer,
    time::Real;
    source_region::Integer = 0,
    batch::Integer = 0,
)
    push!(
        log.events,
        ProvenanceEvent(
            UInt64(id),
            PROV_INJECTION,
            Float64(time),
            Float64(source_region),
            Float64(batch),
        ),
    )
    return log
end

"Record a shock-front crossing for particle `id` at `time`."
function record_crossing!(log::ParticleProvenanceLog, id::Integer, time::Real)
    push!(log.events, ProvenanceEvent(UInt64(id), PROV_CROSSING, Float64(time), 0.0, 0.0))
    return log
end

"Record that particle `id` was reflected at `time`."
function record_reflection!(log::ParticleProvenanceLog, id::Integer, time::Real)
    push!(log.events, ProvenanceEvent(UInt64(id), PROV_REFLECTION, Float64(time), 0.0, 0.0))
    return log
end

"""
    record_max_kinetic_energy!(log, ps, time) -> log

Log each particle's current kinetic energy ½m|v|² at `time`. The per-particle
maximum over the run is recovered by [`provenance_summary`](@ref).
"""
function record_max_kinetic_energy!(
    log::ParticleProvenanceLog,
    ps::ParticleSet{D,T},
    time::Real,
) where {D,T}
    vx, vy, vz = ps.v
    m = ps.m
    @inbounds for p = 1:nparticles(ps)
        ke = 0.5 * m * (vx[p]^2 + vy[p]^2 + vz[p]^2)
        push!(log.events, ProvenanceEvent(ps.id[p], PROV_MAXKE, Float64(time), Float64(ke), 0.0))
    end
    return log
end

"""
    provenance_summary(log, id) -> NamedTuple

Reduce the event log for one particle `id` to the full §9.4 provenance record:
`(; id, species, injection_time, source_region, injection_batch, crossing_count,
first_crossing_time, last_crossing_time, reflection_flag, max_kinetic_energy)`.
Times default to `NaN` and `max_kinetic_energy` to `0.0` when no event is logged.
"""
# ponytail: O(n_events) linear scan per id; build the Dict form below for a full
# reduction, or add an id→indices index if per-id queries become hot.
function provenance_summary(log::ParticleProvenanceLog, id::Integer)
    idU = UInt64(id)
    inj_t = NaN
    src = NaN
    batch = NaN
    nc = 0
    first_c = NaN
    last_c = NaN
    refl = false
    maxke = 0.0
    @inbounds for e in log.events
        e.id == idU || continue
        if e.kind == PROV_INJECTION
            inj_t = e.time
            src = e.a
            batch = e.b
        elseif e.kind == PROV_CROSSING
            nc += 1
            first_c = isnan(first_c) ? e.time : min(first_c, e.time)
            last_c = isnan(last_c) ? e.time : max(last_c, e.time)
        elseif e.kind == PROV_REFLECTION
            refl = true
        elseif e.kind == PROV_MAXKE
            maxke = max(maxke, e.a)
        end
    end
    species = idU >> _PARTICLE_ID_INDEX_BITS
    return (;
        id = idU,
        species,
        injection_time = inj_t,
        source_region = src,
        injection_batch = batch,
        crossing_count = nc,
        first_crossing_time = first_c,
        last_crossing_time = last_c,
        reflection_flag = refl,
        max_kinetic_energy = maxke,
    )
end

"Reduce the whole log to a `Dict{UInt64,NamedTuple}` of per-ID provenance summaries."
function provenance_summary(log::ParticleProvenanceLog)
    ids = Set(e.id for e in log.events)
    return Dict(id => provenance_summary(log, id) for id in ids)
end
