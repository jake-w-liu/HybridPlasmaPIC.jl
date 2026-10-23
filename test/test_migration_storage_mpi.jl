module MigrationStorageMPITests
using HybridPlasmaPIC, Test
import MPI
const H = HybridPlasmaPIC
H.ensure_mpi_initialized!()
const COMM = MPI.COMM_WORLD
const SIZE = MPI.Comm_size(COMM)
const RANK = MPI.Comm_rank(COMM)

struct MetadataBytes <: AbstractVector{UInt8}
    n::Int
    fail::Bool
end
Base.size(a::MetadataBytes) = (a.n,)
Base.getindex(a::MetadataBytes, i::Int) =
    a.fail ? throw(ArgumentError("injected byte packing failure")) : UInt8(0)

const POSITION_READ_FAILURE = ArgumentError("injected position read failure")
struct FailingMigrationPosition{T} <: AbstractVector{T}
    data::Vector{T}
    fail::Base.RefValue{Bool}
end
Base.size(a::FailingMigrationPosition) = size(a.data)
Base.parent(a::FailingMigrationPosition) = a.data
Base.getindex(a::FailingMigrationPosition, i::Int) =
    a.fail[] ? throw(POSITION_READ_FAILURE) : a.data[i]

# Indexed copies are compatible, but this storage cannot accept arrivals.
struct FixedMigrationTags <: AbstractVector{UInt32}
    data::Vector{UInt32}
end
Base.size(a::FixedMigrationTags) = size(a.data)
Base.getindex(a::FixedMigrationTags, i::Int) = a.data[i]
Base.getindex(a::FixedMigrationTags, i::AbstractVector{<:Integer}) = FixedMigrationTags(a.data[i])
Base.parent(a::FixedMigrationTags) = a.data

state(ps) = (map(copy, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)), ps.q, ps.m)
function rejection(f, local_type; peer_text, local_exception = nothing)
    err = try
        f()
        nothing
    catch caught
        caught
    end
    @test err isa (RANK == 0 ? local_type : ErrorException)
    if RANK == 0 && local_exception !== nothing
        @test err === local_exception
    end
    if RANK != 0 && err isa Exception
        @test occursin(peer_text, sprint(showerror, err))
    end
    @test MPI.Allreduce(err === nothing ? 0 : 1, +, COMM) == SIZE
end
function byte_recovery()
    send = [UInt8[RANK, dest-1, 0xff] for dest = 1:SIZE]
    recv = H._mpi_alltoallv_bytes(send, COMM)
    @test recv == [UInt8[src-1, RANK, 0xff] for src = 1:SIZE]
end

@testset "MPI byte displacement metadata boundaries" begin
    for counts in (
        Int[],
        [0],
        [1, 0, 2],
        [typemax(Cint)],
        [typemax(Cint), 1],
        [typemax(Cint), typemax(Cint), 0],
    )
        offsets, total = H._mpi_byte_displacements(counts)
        expected = Cint[
            counts[i] == 0 ? 0 : sum(Int, counts[1:(i-1)]; init = 0) for i in eachindex(counts)
        ]
        @test offsets == expected
        @test total == sum(Int, counts; init = 0)
    end
    @test_throws ArgumentError H._mpi_byte_displacements([-1])
    @test_throws ArgumentError H._mpi_byte_displacements([Int(typemax(Cint)) + 1])
    @test_throws ArgumentError H._mpi_byte_displacements([typemax(Cint), 1, 1])
    @test_throws ArgumentError H._mpi_byte_displacements([typemax(Cint), 1, 0, 1])
end

@testset "MPI byte preparation rejects collectively and recovers" begin
    for failure in (:length, :count, :packing)
        send = AbstractVector{UInt8}[UInt8[] for _ = 1:SIZE]
        if RANK == 0
            if failure == :length
                pop!(send)
            elseif failure == :count
                send[1] = MetadataBytes(Int(typemax(Cint)) + 1, false)
            else
                send[1] = MetadataBytes(1, true)
            end
        end
        rejection(ArgumentError; peer_text = "MPI byte transport send preparation") do
            H._mpi_alltoallv_bytes(send, COMM)
        end
        byte_recovery()
    end
    if SIZE >= 3
        send = AbstractVector{UInt8}[UInt8[] for _ = 1:SIZE]
        if RANK == 0
            send[1] = MetadataBytes(Int(typemax(Cint)), false)
            send[2] = MetadataBytes(1, false)
            send[3] = MetadataBytes(1, false)
        end
        rejection(ArgumentError; peer_text = "MPI byte transport send preparation") do
            H._mpi_alltoallv_bytes(send, COMM)
        end
        byte_recovery()
    end
    # Checkpoint callers delegate to the same boundary, retaining their peer text.
    ctx = (; comm = COMM)
    rejection(ArgumentError; peer_text = "MPI checkpoint test rejection") do
        H._mpi_checkpoint_stage(ctx, "test rejection") do
            RANK == 0 && throw(ArgumentError("checkpoint local exception"))
            nothing
        end
    end
    @test H._mpi_checkpoint_stage(() -> 17, ctx, "recovery") == 17
    @test H._mpi_alltoallv_bytes([UInt8[] for _ = 1:SIZE], COMM) == [UInt8[] for _ = 1:SIZE]
end

function population(D, T, ctx)
    ps = ParticleSet{D,T}(1; q = T(2), m = T(3))
    destination = mod(ctx.logical_rank, SIZE) + 1
    bounds =
        rank_bounds(FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D)), ctx.layout, destination)
    for d = 1:D
        ps.x[d][1] = (bounds.lo[d] + bounds.hi[d]) / T(2) + one(T)
    end
    for c = 1:3
        ps.v[c][1] = T(10c + ctx.logical_rank)
    end
    ps.weight[1] = T(ctx.logical_rank) + T(0.5)
    ps.id[1] = UInt64(ctx.logical_rank)
    ps.tag[1] = UInt32(100 + ctx.logical_rank)
    return ps
end
function migrate_recovery(D, T, ctx, g)
    ps = population(D, T, ctx)
    original = state(ps)
    arrays = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
    stats = mpi_migrate_particles!(ps, g, ctx)
    source = mod(ctx.logical_rank - 2, SIZE) + 1
    @test stats == (
        moved = SIZE == 1 ? 0 : SIZE,
        lost = 0,
        sent = SIZE == 1 ? 0 : 1,
        received = SIZE == 1 ? 0 : 1,
    )
    @test ps.id == UInt64[source]
    @test ps.tag == UInt32[100+source]
    @test ps.weight == T[T(source)+T(0.5)]
    @test all(ps.v[c] == T[10c+source] for c = 1:3)
    @test rank_of_position(ntuple(d -> ps.x[d][1], D), g, ctx.layout) == ctx.logical_rank
    @test all(zero(T) <= ps.x[d][1] < one(T) for d = 1:D)
    @test map(copy, arrays) == original[1]
end

@testset "MPI migration rejection preserves every rank and recovers" begin
    D = SIZE == 8 ? 3 : SIZE == 4 ? 2 : 1
    shape = SIZE in (4, 8) ? ntuple(_ -> 2, D) : (SIZE,)
    layout = LogicalRankLayout(shape)
    ctx = create_cartesian_communicator(layout; comm = COMM)
    for T in (Float32, Float64)
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> one(T), D))
        for failure in (:axes, :position_read, :conversion, :assembly)
            SIZE == 1 && failure == :assembly && continue
            ps = population(D, T, ctx)
            fail_read = Ref(false)
            if RANK == 0
                if failure == :axes
                    pop!(ps.tag)
                elseif failure == :position_read
                    positions = map(a -> FailingMigrationPosition(a, fail_read), ps.x)
                    ps = ParticleSet{D,T}(positions, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
                else
                    tag = failure == :conversion ? view(ps.tag, :) : FixedMigrationTags(ps.tag)
                    ps = ParticleSet{D,T}(ps.x, ps.v, ps.weight, ps.id, tag, ps.q, ps.m)
                end
            end
            before = state(ps)
            refs = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
            local_type = failure == :conversion ? MethodError : ArgumentError
            fail_read[] = failure == :position_read
            local_exception = failure == :position_read ? POSITION_READ_FAILURE : nothing
            rejection(local_type; peer_text = "MPI particle migration", local_exception) do
                mpi_migrate_particles!(ps, g, ctx)
            end
            fail_read[] = false
            @test state(ps) == before
            @test all(
                (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)[k] === refs[k] for k in eachindex(refs)
            )
            migrate_recovery(D, T, ctx, g)
        end
        # Positions may share read-only inputs with weights; wrapping must not change weights.
        ps = population(D, T, ctx)
        ps.weight = ps.x[1]
        original_weight = copy(ps.weight)
        expected_source = mod(ctx.logical_rank - 2, SIZE) + 1
        all_weights = MPI.Allgather(original_weight[1], ctx.comm)
        mpi_migrate_particles!(ps, g, ctx)
        source_mpi = H._mpi_rank_for_logical(ctx, expected_source)
        @test ps.weight == T[all_weights[source_mpi+1]]
        @test !H._particle_components_mightalias(ps)
    end
end
end
