using HybridPlasmaPIC, Test, Random, LinearAlgebra

# Independent (brute-force) column-major linear index oracle, computed directly
# from positions, dx and n — no reuse of the production strides code.
function ref_cell_index(ps::ParticleSet{D,T}, g) where {D,T}
    N = nparticles(ps)
    out = Vector{Int}(undef, N)
    for p = 1:N
        sub = ntuple(d -> clamp(floor(Int, ps.x[d][p] / g.dx[d]), 0, g.n[d] - 1), D)
        # LinearIndices is column-major; +1 because sub is 0-based
        out[p] = LinearIndices(g.n)[(sub .+ 1)...]
    end
    return out
end

function brute_minimax_load(percell, ntiles)
    loads = Int.(percell)
    n = length(loads)
    n == 0 && return 0
    k = min(Int(ntiles), n)
    best = typemax(Int)

    function search(start, parts_left, current_max)
        if parts_left == 1
            tail = sum(view(loads, start:n))
            best = min(best, max(current_max, tail))
            return nothing
        end

        acc = 0
        last_stop = n - parts_left + 1
        for stop = start:last_stop
            acc += loads[stop]
            next_max = max(current_max, acc)
            next_max >= best && continue
            search(stop + 1, parts_left - 1, next_max)
        end
        return nothing
    end

    search(1, k, 0)
    return best
end

function ordered_cell_cover(ranges, ncells)
    covered = Int[]
    for r in ranges
        append!(covered, collect(r))
    end
    return covered == collect(1:ncells)
end

@testset "particle_sort" begin

    @testset "cell_index matches oracle (D=$D)" for D in (1, 2, 3)
        rng = MersenneTwister(100 + D)
        n = ntuple(d -> 4 + d, D)            # distinct sizes -> catches stride bugs
        L = ntuple(d -> Float64(n[d]), D)    # dx = 1 here, but code uses dx generally
        g = FourierGrid(n, L)
        N = 500
        ps = ParticleSet{D,Float64}(N)
        load_uniform!(ps, rng, ntuple(_ -> 0.0, D), L)
        ci = cell_index(ps, g)
        @test ci == ref_cell_index(ps, g)
        @test all(1 .<= ci .<= prod(n))
    end

    @testset "cell_index clamps out-of-box positions" begin
        n = (5, 3)
        L = (5.0, 3.0)
        g = FourierGrid(n, L)
        ps = ParticleSet{2,Float64}(4)
        # below lo, above hi, exactly on edges
        ps.x[1] .= [-2.0, 100.0, 0.0, 4.9999]
        ps.x[2] .= [-0.5, 50.0, 0.0, 2.9999]
        ci = cell_index(ps, g)
        @test ci == ref_cell_index(ps, g)
        @test all(1 .<= ci .<= 15)
        # particle 1 clamps to cell (0,0)->1 ; particle 2 clamps to (4,2)->15
        @test ci[1] == 1
        @test ci[2] == 15
    end

    @testset "sort_particles! sorts + preserves data (D=$D)" for D in (1, 2, 3)
        rng = MersenneTwister(7 + D)
        n = ntuple(d -> 3 + 2d, D)
        L = ntuple(d -> 2.0 * n[d], D)       # dx = 2 -> exercises division
        g = FourierGrid(n, L)
        N = 2000
        ps = ParticleSet{D,Float64}(N)
        load_uniform!(ps, rng, ntuple(_ -> 0.0, D), L)
        load_maxwellian!(ps, rng, (0.1, -0.2, 0.3), (1.0, 1.0, 1.0))
        ps.weight .= rand(rng, N)
        ps.id .= UInt64.(randperm(rng, N))       # scrambled unique ids
        ps.tag .= UInt32.(rand(rng, 0:9, N))

        # snapshots keyed by id to verify per-particle data travels together
        ids0 = copy(ps.id)
        x0 = ntuple(d -> copy(ps.x[d]), D)
        v0 = ntuple(c -> copy(ps.v[c]), 3)
        w0 = copy(ps.weight)
        tag0 = copy(ps.tag)
        byid = Dict(ids0[p] => p for p = 1:N)

        ret = sort_particles!(ps, g)
        @test ret === ps

        # 1) cell index is non-decreasing
        ci = cell_index(ps, g)
        @test issorted(ci)

        # 2) multiset of ids unchanged
        @test sort(ps.id) == sort(ids0)

        # 3) per-particle data still consistent (look up each id's original row)
        for p = 1:N
            o = byid[ps.id[p]]
            for d = 1:D
                @test ps.x[d][p] == x0[d][o]
            end
            for c = 1:3
                @test ps.v[c][p] == v0[c][o]
            end
            @test ps.weight[p] == w0[o]
            @test ps.tag[p] == tag0[o]
        end

        # 4) stable within a cell: original id order preserved among equal cells
        for c in unique(ci)
            rows = findall(==(c), ci)
            orig = [byid[ps.id[r]] for r in rows]
            @test issorted(orig)
        end
    end

    @testset "particles_per_cell histogram (D=$D)" for D in (1, 2, 3)
        rng = MersenneTwister(55 + D)
        n = ntuple(d -> 4 + d, D)
        L = ntuple(d -> 1.5 * n[d], D)
        g = FourierGrid(n, L)
        N = 1234
        ps = ParticleSet{D,Float64}(N)
        load_uniform!(ps, rng, ntuple(_ -> 0.0, D), L)

        ppc = particles_per_cell(ps, g)
        @test length(ppc) == prod(n)
        @test sum(ppc) == N

        # direct histogram from the oracle indices
        ref = zeros(Int, prod(n))
        for c in ref_cell_index(ps, g)
            ref[c] += 1
        end
        @test ppc == ref

        # invariant under sorting (same multiset of cells)
        sort_particles!(ps, g)
        @test particles_per_cell(ps, g) == ppc
    end

    @testset "balanced particle load partition" begin
        cases = (
            ([9, 1, 1, 1, 1, 9], 3),
            ([4, 4, 4, 4, 4], 2),
            ([0, 7, 0, 3, 10, 0, 2], 4),
            ([80, 10, 0, 0, 0, 0, 0, 0, 10, 0], 3),
        )

        for (percell, ntiles) in cases
            ranges = balanced_tile_ranges(percell, ntiles)
            loads = balanced_tile_loads(percell, ranges)
            @test length(ranges) == ntiles
            @test length(loads) == ntiles
            @test ordered_cell_cover(ranges, length(percell))
            @test sum(loads) == sum(percell)
            @test maximum(loads) == brute_minimax_load(percell, ntiles)
            @test load_imbalance(loads) <= load_imbalance(tile_loads(percell, ntiles))
        end

        zero_ranges = balanced_tile_ranges(zeros(Int, 3), 5)
        zero_loads = balanced_tile_loads(zeros(Int, 3), zero_ranges)
        @test length(zero_ranges) == 5
        @test ordered_cell_cover(zero_ranges, 3)
        @test zero_loads == zeros(Int, 5)

        empty_ranges = balanced_tile_ranges(Int[], 3)
        @test length(empty_ranges) == 3
        @test all(isempty, empty_ranges)
        @test balanced_tile_loads(Int[], empty_ranges) == zeros(Int, 3)

        @test_throws ArgumentError balanced_tile_ranges([1, -1, 2], 2)
        @test_throws ArgumentError balanced_tile_ranges([1, 2, 3], 0)
        @test_throws ArgumentError balanced_tile_loads([1, 2], [1:3])
    end

    @testset "particle_load_balance wrapper" begin
        g = FourierGrid((10,), (10.0,))
        ps = ParticleSet{1,Float64}(100)
        ps.x[1][1:80] .= 0.25
        ps.x[1][81:90] .= 1.25
        ps.x[1][91:100] .= 8.25

        percell = particles_per_cell(ps, g)
        balanced = particle_load_balance(ps, g; ntiles = 3)
        equal = particle_load_imbalance(ps, g; ntiles = 3)

        @test ordered_cell_cover(balanced.ranges, length(percell))
        @test sum(balanced.per_tile) == nparticles(ps)
        @test maximum(balanced.per_tile) == brute_minimax_load(percell, 3)
        @test balanced.imbalance < equal.imbalance
        @test balanced.per_tile == balanced_tile_loads(percell, balanced.ranges)
    end

    @testset "sorting/load helpers reject invalid inputs" begin
        g = FourierGrid((4,), (4.0,))

        ps1 = ParticleSet{1,Float64}(1)
        ps1.x[1][1] = NaN
        @test_throws ArgumentError cell_index(ps1, g)
        @test_throws ArgumentError particles_per_cell(ps1, g)
        @test_throws ArgumentError particle_load_imbalance(ps1, g; ntiles = 2)
        @test_throws ArgumentError particle_load_balance(ps1, g; ntiles = 2)
        @test_throws ArgumentError sort_particles!(ps1, g)

        ps2 = ParticleSet{1,Float64}(2)
        ps2.x[1] .= [NaN, 1.0]
        @test_throws ArgumentError sort_particles!(ps2, g)

        @test_throws ArgumentError load_imbalance([1.0, NaN])
        @test_throws ArgumentError load_imbalance([1.0, Inf])
        @test_throws ArgumentError load_imbalance([1.0, -1.0])
        @test_throws ArgumentError tile_loads([1, -1, 2], 2)
    end

    @testset "memory_bytes formula" begin
        # hand-computed case: D=3, Tbytes=4 -> (3+3)*4+8+4 = 36 per particle
        @test memory_bytes(ncells = 1000, nppc = 50, D = 3) == 1000 * 50 * ((3 + 3) * 4 + 8 + 4)
        @test memory_bytes(ncells = 1000, nppc = 50, D = 3) == 1000 * 50 * 36

        # multi-species scales linearly
        @test memory_bytes(ncells = 10, nppc = 4, nspecies = 3, D = 2) ==
              3 * 10 * 4 * ((2 + 3) * 4 + 8 + 4)

        # Float64 payload, D=1
        @test memory_bytes(ncells = 7, nppc = 9, D = 1, Tbytes = 8) == 7 * 9 * ((1 + 3) * 8 + 8 + 4)

        # zero particles -> zero bytes
        @test memory_bytes(ncells = 0, nppc = 50, D = 3) == 0
        @test memory_bytes(ncells = 100, nppc = 0, D = 3) == 0

        # bad args rejected
        @test_throws ArgumentError memory_bytes(ncells = -1, nppc = 1, D = 3)
        @test_throws ArgumentError memory_bytes(ncells = 1, nppc = 1, D = 0)
        @test_throws ArgumentError memory_bytes(ncells = 1, nppc = 1, D = 3, Tbytes = 0)
        @test_throws OverflowError memory_bytes(ncells = typemax(Int), nppc = 2, D = 3)
    end

    @testset "edge cases: empty / single particle" begin
        n = (4, 4)
        L = (4.0, 4.0)
        g = FourierGrid(n, L)
        ps0 = ParticleSet{2,Float64}(0)
        @test cell_index(ps0, g) == Int[]
        @test sort_particles!(ps0, g) === ps0
        @test sum(particles_per_cell(ps0, g)) == 0
        @test length(particles_per_cell(ps0, g)) == 16

        ps1 = ParticleSet{2,Float64}(1)
        ps1.x[1] .= [2.5]
        ps1.x[2] .= [1.5]
        ci = cell_index(ps1, g)
        @test ci == ref_cell_index(ps1, g)
        sort_particles!(ps1, g)
        @test cell_index(ps1, g) == ci
    end
end
