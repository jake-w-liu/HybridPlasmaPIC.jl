# particle_sort.jl — spatial (cell) sorting of the structure-of-arrays particle
# storage, a per-cell occupancy histogram, and a static memory estimator.
#
# Sorting particles by their column-major cell index improves cache locality for
# deposition/gather and makes per-cell operations (binning, diagnostics)
# contiguous. The grid geometry comes from the FourierGrid (n, dx); cells are
# indexed exactly like the deposition mesh: cell_d = clamp(floor(x_d/dx_d), 0, n_d-1).

"""
    cell_index(ps::ParticleSet{D,T}, g::FourierGrid{D,T}) -> Vector{Int}

Column-major linear cell index of every particle. For particle `p` the per-axis
cell is `cell_d = clamp(floor(x_d / dx_d), 0, n_d-1)` and the returned 1-based
linear index is `1 + Σ_d cell_d * stride_d` with `stride_1 = 1` and
`stride_d = prod(n[1:d-1])`. Indices lie in `1:prod(g.n)`.
"""
function cell_index(ps::ParticleSet{D,T}, g::FourierGrid{D,T}) where {D,T}
    N = nparticles(ps)
    n = g.n
    dx = g.dx
    # column-major strides: stride[1]=1, stride[d]=prod(n[1:d-1])
    strides = ntuple(d -> d == 1 ? 1 : prod(ntuple(k -> n[k], d - 1)), D)
    out = Vector{Int}(undef, N)
    @inbounds for p = 1:N
        lin = 1
        for d = 1:D
            c = floor(Int, ps.x[d][p] / dx[d])
            if c < 0
                c = 0
            elseif c > n[d] - 1
                c = n[d] - 1
            end
            lin += c * strides[d]
        end
        out[p] = lin
    end
    return out
end

"""
    sort_particles!(ps::ParticleSet{D,T}, g::FourierGrid{D,T}) -> ps

Stably reorder every particle array (`x`, `v`, `weight`, `id`, `tag`) in place by
ascending [`cell_index`](@ref). The same stable permutation is applied to all
arrays, so per-particle data stays consistent and the original order is preserved
within a cell. Returns `ps`.
"""
function sort_particles!(ps::ParticleSet{D,T}, g::FourierGrid{D,T}) where {D,T}
    N = nparticles(ps)
    N <= 1 && return ps                       # nothing to reorder
    ci = cell_index(ps, g)
    perm = sortperm(ci; alg = MergeSort)      # stable sort
    for d = 1:D
        ps.x[d] .= ps.x[d][perm]
    end
    for c = 1:3
        ps.v[c] .= ps.v[c][perm]
    end
    ps.weight .= ps.weight[perm]
    ps.id .= ps.id[perm]
    ps.tag .= ps.tag[perm]
    return ps
end

"""
    particles_per_cell(ps::ParticleSet{D,T}, g::FourierGrid{D,T}) -> Vector{Int}

Histogram of particle occupancy over the `prod(g.n)` column-major cells; entry
`i` counts the particles whose [`cell_index`](@ref) equals `i`. The result sums
to `nparticles(ps)`.
"""
function particles_per_cell(ps::ParticleSet{D,T}, g::FourierGrid{D,T}) where {D,T}
    ncells = prod(g.n)
    counts = zeros(Int, ncells)
    ci = cell_index(ps, g)
    @inbounds for c in ci
        counts[c] += 1
    end
    return counts
end

"""
    load_imbalance(counts) -> Float64

Load-imbalance factor `max(counts) / mean(counts)` of a per-partition work
histogram: `1.0` is perfect balance, `P` is the worst case (all work on one of
`P` partitions). Returns `1.0` for an empty or all-zero histogram. This is the
CPU-side metric a domain decomposition would minimize (item: load-imbalance
metrics) — no MPI required to compute or to study a proposed partition.
"""
function load_imbalance(counts::AbstractArray{<:Real})
    n = length(counts)
    n == 0 && return 1.0
    tot = sum(counts)
    tot == 0 && return 1.0
    return Float64(maximum(counts)) / (Float64(tot) / n)
end

"""
    tile_loads(percell::AbstractVector{<:Integer}, ntiles::Integer) -> Vector{Int}

Sum a per-cell occupancy histogram into `ntiles` contiguous (column-major) tiles
of as-equal-as-possible cell count — a 1-D slab decomposition of the work. Used
with [`load_imbalance`](@ref) to evaluate a candidate partition.
"""
function tile_loads(percell::AbstractVector{<:Integer}, ntiles::Integer)
    ntiles >= 1 || throw(ArgumentError("ntiles must be ≥ 1"))
    ncells = length(percell)
    loads = zeros(Int, ntiles)
    @inbounds for c = 1:ncells
        t = min(ntiles, (c - 1) * ntiles ÷ max(ncells, 1) + 1)
        loads[t] += percell[c]
    end
    return loads
end

"""
    particle_load_imbalance(ps, g; ntiles) -> (; per_tile, imbalance)

Particle load per contiguous tile and the [`load_imbalance`](@ref) factor for a
`ntiles`-way slab decomposition of grid `g`. A clustered (shock-compressed) load
gives `imbalance > 1`; a uniform load gives `imbalance ≈ 1`.
"""
function particle_load_imbalance(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T};
    ntiles::Integer,
) where {D,T}
    per_tile = tile_loads(particles_per_cell(ps, g), ntiles)
    return (; per_tile, imbalance = load_imbalance(per_tile))
end

function _validated_loads(percell::AbstractVector{<:Integer})
    loads = Vector{Int}(undef, length(percell))
    for (j, i) in enumerate(eachindex(percell))
        w = percell[i]
        w >= 0 || throw(ArgumentError("per-cell loads must be nonnegative"))
        loads[j] = Int(w)
    end
    return loads
end

function _equal_cell_ranges(ncells::Integer, ntiles::Integer)
    ranges = Vector{UnitRange{Int}}(undef, ntiles)
    base = ncells ÷ ntiles
    extra = ncells % ntiles
    start = 1
    for t = 1:ntiles
        width = base + (t <= extra ? 1 : 0)
        ranges[t] = start:(start+width-1)
        start += width
    end
    return ranges
end

function _can_partition_with_cap(percell::AbstractVector{Int}, nparts::Integer, cap::Integer)
    nparts >= 1 || return isempty(percell)
    isempty(percell) && return true
    parts = 1
    acc = 0
    @inbounds for w in percell
        w > cap && return false
        if acc + w <= cap
            acc += w
        else
            parts += 1
            parts > nparts && return false
            acc = w
        end
    end
    return true
end

function _minimax_load_cap(percell::AbstractVector{Int}, nparts::Integer)
    isempty(percell) && return 0
    lo = maximum(percell)
    hi = sum(percell)
    while lo < hi
        mid = (lo + hi) >>> 1
        if _can_partition_with_cap(percell, nparts, mid)
            hi = mid
        else
            lo = mid + 1
        end
    end
    return lo
end

function _ranges_for_cap(percell::AbstractVector{Int}, nparts::Integer, cap::Integer)
    ncells = length(percell)
    nparts == 0 && return UnitRange{Int}[]
    ranges = UnitRange{Int}[]
    start = 1
    acc = 0
    @inbounds for c = 1:ncells
        w = percell[c]
        if c > start && acc + w > cap && length(ranges) < nparts - 1
            push!(ranges, start:(c-1))
            start = c
            acc = 0
        end
        acc += w

        remaining_cells = ncells - c
        remaining_parts = nparts - length(ranges) - 1
        if remaining_parts > 0 && remaining_cells == remaining_parts
            push!(ranges, start:c)
            start = c + 1
            acc = 0
        end
    end
    push!(ranges, start:ncells)
    return ranges
end

"""
    balanced_tile_ranges(percell, ntiles) -> Vector{UnitRange{Int}}

Partition a column-major per-cell particle histogram into exactly `ntiles`
contiguous ranges. For positive total load, the nonempty ranges minimize the
maximum tile load, so the returned partition is the deterministic 1-D
load-balancing plan a slab decomposition would use before moving rank
boundaries. If `ntiles > length(percell)`, the extra trailing ranges are empty.

Cells with zero load are still assigned to preserve a complete domain covering.
For all-zero histograms, cells are split evenly by count because every load
partition is equivalent.
"""
function balanced_tile_ranges(percell::AbstractVector{<:Integer}, ntiles::Integer)
    ntiles >= 1 || throw(ArgumentError("ntiles must be >= 1"))
    loads = _validated_loads(percell)
    ncells = length(loads)
    ncells == 0 && return [1:0 for _ = 1:ntiles]
    sum(loads) == 0 && return _equal_cell_ranges(ncells, Int(ntiles))

    nonempty = min(Int(ntiles), ncells)
    cap = _minimax_load_cap(loads, nonempty)
    ranges = _ranges_for_cap(loads, nonempty, cap)
    while length(ranges) < ntiles
        push!(ranges, (ncells+1):ncells)
    end
    return ranges
end

"""
    balanced_tile_loads(percell, ranges) -> Vector{Int}

Load of each tile range returned by [`balanced_tile_ranges`](@ref). Empty
ranges contribute zero.
"""
function balanced_tile_loads(
    percell::AbstractVector{<:Integer},
    ranges::AbstractVector{<:UnitRange{Int}},
)
    loads = _validated_loads(percell)
    ncells = length(loads)
    out = Vector{Int}(undef, length(ranges))
    @inbounds for (i, r) in pairs(ranges)
        if isempty(r)
            out[i] = 0
        else
            first(r) >= 1 && last(r) <= ncells ||
                throw(ArgumentError("tile range $r lies outside 1:$ncells"))
            out[i] = sum(view(loads, r))
        end
    end
    return out
end

"""
    particle_load_balance(ps, g; ntiles) -> (; ranges, per_tile, imbalance)

Compute a deterministic contiguous-cell load-balancing plan for particles on
grid `g`. `ranges` are column-major cell ranges, `per_tile` is the particle
count assigned to each range, and `imbalance` is [`load_imbalance`](@ref) of the
balanced plan.
"""
function particle_load_balance(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T};
    ntiles::Integer,
) where {D,T}
    percell = particles_per_cell(ps, g)
    ranges = balanced_tile_ranges(percell, ntiles)
    per_tile = balanced_tile_loads(percell, ranges)
    return (; ranges, per_tile, imbalance = load_imbalance(per_tile))
end

"""
    memory_bytes(; ncells, nppc, nspecies=1, D, Tbytes=4) -> Int

Estimated particle-storage footprint in bytes for `nspecies` species, each with
`ncells * nppc` particles. Per particle: `D` positions and `3` velocities of
`Tbytes` each, an 8-byte `UInt64` id, and a 4-byte `UInt32` tag:

    nspecies * ncells * nppc * ((D + 3) * Tbytes + 8 + 4)

This counts only the per-particle SoA payload (not the field grids or scratch).
"""
@inline function _checked_memory_bytes_per_particle(D::Int, Tbytes::Int)
    scalars = Base.checked_mul(Base.checked_add(D, 3), Tbytes)
    return Base.checked_add(Base.checked_add(scalars, 8), 4)
end

@inline function _checked_memory_product(nspecies::Int, ncells::Int, nppc::Int, per_particle::Int)
    return Base.checked_mul(
        Base.checked_mul(Base.checked_mul(nspecies, ncells), nppc),
        per_particle,
    )
end

function memory_bytes(; ncells::Int, nppc::Int, nspecies::Int = 1, D::Int, Tbytes::Int = 4)
    ncells >= 0 || throw(ArgumentError("ncells must be ≥ 0"))
    nppc >= 0 || throw(ArgumentError("nppc must be ≥ 0"))
    nspecies >= 0 || throw(ArgumentError("nspecies must be ≥ 0"))
    D >= 1 || throw(ArgumentError("D must be ≥ 1"))
    Tbytes >= 1 || throw(ArgumentError("Tbytes must be ≥ 1"))
    per_particle = _checked_memory_bytes_per_particle(D, Tbytes)
    return _checked_memory_product(nspecies, ncells, nppc, per_particle)
end
