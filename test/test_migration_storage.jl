module MigrationStorageTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC

components(ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
state(ps) = (map(copy, components(ps)), ps.q, ps.m)
function population(D, T, xs, firstid = 0)
    ps = ParticleSet{D,T}(length(xs); q = T(2), m = T(3))
    for d = 1:D
        ps.x[d] .= d == 1 ? T.(xs) : T(1.25)
    end
    for c = 1:3
        ps.v[c] .= T.(firstid .+ (1:length(xs)) .+ 10c)
    end
    ps.weight .= T.(firstid .+ (1:length(xs)) .+ 0.5)
    ps.id .= UInt64.(firstid .+ (1:length(xs)))
    ps.tag .= UInt32.(firstid .+ (1:length(xs)) .+ 100)
    return ps
end
function rebuild(ps::ParticleSet{D,T}, a; abstract_storage = false) where {D,T}
    x, v = ntuple(d -> a[d], D), ntuple(c -> a[D+c], 3)
    if abstract_storage
        return ParticleSet{
            D,
            T,
            AbstractVector{T},
            AbstractVector{T},
            AbstractVector{T},
            AbstractVector{UInt64},
            AbstractVector{UInt32},
        }(
            x,
            v,
            a[D+4],
            a[D+5],
            a[D+6],
            ps.q,
            ps.m,
        )
    end
    return ParticleSet{D,T}(x, v, a[D+4], a[D+5], a[D+6], ps.q, ps.m)
end
function payloads(ranks)
    Dict(
        ps.id[p] => (ntuple(c -> ps.v[c][p], 3), ps.weight[p], ps.tag[p], ps.q, ps.m) for
        ps in ranks for p = 1:nparticles(ps)
    )
end

# A compatible concrete vector: indexed copies retain their destination type.
struct MigrationVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(a::MigrationVector) = size(a.data)
Base.getindex(a::MigrationVector, i::Int) = a.data[i]
Base.getindex(a::MigrationVector, i::AbstractVector{<:Integer}) = MigrationVector(a.data[i])
Base.setindex!(a::MigrationVector, v, i::Int) = (a.data[i] = v)
Base.parent(a::MigrationVector) = a.data
Base.resize!(a::MigrationVector, n::Integer) = (resize!(a.data, n); a)
Base.append!(a::MigrationVector, b::AbstractVector) = (append!(a.data, b); a)

struct SharedIndexVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(a::SharedIndexVector) = size(a.data)
Base.getindex(a::SharedIndexVector, i::Int) = a.data[i]
Base.getindex(a::SharedIndexVector, ::AbstractVector{<:Integer}) = a
Base.parent(a::SharedIndexVector) = a.data

struct ShiftedRanks{P} <: AbstractVector{P}
    data::Vector{P}
end
Base.size(a::ShiftedRanks) = size(a.data)
Base.axes(a::ShiftedRanks) = (0:(length(a.data)-1),)
Base.getindex(a::ShiftedRanks, i::Int) = a.data[i+1]

@testset "migration prepares all fields and ranks before publication" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D))
        layout = LogicalRankLayout(ntuple(d -> d == 1 ? 2 : 1, D))
        for late_rank in (1, 2), field in (:x, :v, :weight, :id, :tag)
            ranks = ParticleSet{D,T}[population(D, T, [1.75], 0), population(D, T, [0.25], 10)]
            original = ranks[late_rank]
            a = components(original)
            selected =
                field == :x ? (1:D) :
                field == :v ? ((D+1):(D+3)) :
                field == :weight ? (D+4:D+4) : field == :id ? (D+5:D+5) : (D+6:D+6)
            ranks[late_rank] =
                rebuild(original, ntuple(k -> k in selected ? view(a[k], :) : a[k], length(a)))
            before, owners = map(state, ranks), copy(ranks)
            refs = map(components, ranks)
            @test_throws MethodError migrate_particles!(ranks, g, layout)
            @test map(state, ranks) == before
            @test all(ranks[r] === owners[r] for r = 1:2)
            @test all(components(ranks[r])[k] === refs[r][k] for r = 1:2 for k = 1:(D+6))
        end
        a = population(D, T, [0.25, 0.75])
        before = state(a)
        @test_throws ArgumentError migrate_particles!([a, a], g, layout)
        @test state(a) == before
        ranks = [population(D, T, [1.25]), population(D, T, [0.75], 10)]
        pop!(ranks[2].tag)
        before = map(state, ranks)
        @test_throws ArgumentError migrate_particles!(ranks, g, layout)
        @test map(state, ranks) == before
    end
end

@testset "migration preserves aliased input values and backing storage" begin
    for T in (Float32, Float64), D = 1:3, cross_population in (false, true), k = 1:(D+4)
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D))
        layout = LogicalRankLayout(ntuple(d -> d == 1 ? 2 : 1, D))
        a, b = population(D, T, [1.75]), population(D, T, [1.25], 10)
        target = cross_population ? b : a
        fields = components(target)
        changed = rebuild(target, ntuple(j -> j == k ? a.x[1] : fields[j], length(fields)))
        if cross_population
            b = changed
        else
            a = changed
        end
        ranks = [a, b]
        expected = payloads(ranks)
        originals = map(components, ranks)
        before = map(state, ranks)
        stats = migrate_particles!(ranks, g, layout)
        @test stats.lost == 0
        @test payloads(ranks) == expected
        @test all(map(copy, originals[r]) == before[r][1] for r = 1:2)
        @test ranks[1] === a && ranks[2] === b
        @test all(!H._particle_components_mightalias(ps) for ps in ranks)
        @test !H._particle_sets_mightalias(a, b)
        @test all(
            rank_of_position(ntuple(d -> ps.x[d][p], D), g, layout) == r for (r, ps) in
                                                                             pairs(ranks) for
            p = 1:nparticles(ps)
        )
    end
end

@testset "migration compatible storage, empty ranks, loss and ordering" begin
    for T in (Float32, Float64), D = 1:3, storage in (:vector, :custom, :abstract_views)
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D))
        for periodic in (false, true), xs in ([1.75, -0.25], [0.25, 0.75], [0.25, 0.125], Float64[])
            layout = LogicalRankLayout(
                ntuple(d -> d == 1 ? 2 : 1, D);
                periodic = ntuple(d -> d == 1 ? periodic : true, D),
            )
            a, b = population(D, T, xs), population(D, T, [0.25, 0.75], 10)
            if storage == :custom
                a, b = map(
                    ps -> rebuild(ps, map(v -> MigrationVector(copy(v)), components(ps))),
                    (a, b),
                )
            elseif storage == :abstract_views
                a, b = map(
                    ps -> rebuild(
                        ps,
                        map(v -> view(copy(v), :), components(ps));
                        abstract_storage = true,
                    ),
                    (a, b),
                )
            end
            ranks = [a, b]
            expected_payload = payloads(ranks)
            expected_ids = [UInt64[] for _ = 1:2]
            # Retained particles precede arrivals; arrivals follow source-rank order.
            for retained in (true, false), (r, ps) in pairs(ranks), p = 1:nparticles(ps)
                dest = rank_of_position(ntuple(d -> ps.x[d][p], D), g, layout)
                dest === nothing && continue
                (dest == r) == retained && push!(expected_ids[dest], ps.id[p])
            end
            before = map(components, ranks)
            snapshots = map(state, ranks)
            stats = migrate_particles!(ranks, g, layout)
            @test map(ps -> collect(ps.id), ranks) == expected_ids
            @test stats.lost == length(xs) + 2 - sum(length, expected_ids)
            @test payloads(ranks) ==
                  Dict(id => expected_payload[id] for id in vcat(expected_ids...))
            @test all(map(copy, before[r]) == snapshots[r][1] for r = 1:2)
            @test all(zero(T) <= x < one(T) for ps in ranks for xd in ps.x for x in xd)
        end
    end
    ps = population(1, Float64, [1.25])
    a = components(ps)
    ps = rebuild(ps, ntuple(k -> k == 1 ? SharedIndexVector(a[k]) : a[k], length(a)))
    before = state(ps)
    @test_throws ArgumentError migrate_particles!(
        [ps],
        FourierGrid((8,), (1.0,)),
        LogicalRankLayout((1,)),
    )
    @test state(ps) == before
end

@testset "migration rank axes and disjoint source views" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D))
        layout = LogicalRankLayout(ntuple(d -> d == 1 ? 2 : 1, D))
        ranks = [population(D, T, [1.75]), population(D, T, [1.25], 10)]
        before = map(state, ranks)
        @test_throws ArgumentError migrate_particles!(ShiftedRanks(ranks), g, layout)
        @test map(state, ranks) == before
        parents = Vector{T}[]
        viewed = map(ranks) do ps
            a = components(ps)
            parent_storage = vcat(a[1:(D+4)]...)
            push!(parents, parent_storage)
            arrays = ntuple(k -> k <= D + 4 ? view(parent_storage, k:k) : a[k], D + 6)
            rebuild(ps, arrays; abstract_storage = true)
        end
        original_parents = map(copy, parents)
        expected = payloads(viewed)
        stats = migrate_particles!(viewed, g, layout)
        @test stats == (moved = 2, lost = 0)
        @test payloads(viewed) == expected
        @test parents == original_parents
        @test viewed[1].id == UInt64[11] && viewed[2].id == UInt64[1]
    end
end
end
