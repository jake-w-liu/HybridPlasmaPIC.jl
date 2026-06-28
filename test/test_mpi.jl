using HybridPlasmaPIC, Test
import MPI
import Serialization

struct WrappedMPIArray{T,N,A<:Array{T,N}} <: AbstractArray{T,N}
    data::A
end

Base.size(A::WrappedMPIArray) = size(A.data)
Base.axes(A::WrappedMPIArray) = axes(A.data)
Base.IndexStyle(::Type{<:WrappedMPIArray}) = IndexLinear()
Base.getindex(A::WrappedMPIArray, i::Int) = A.data[i]
Base.getindex(A::WrappedMPIArray, I::Vararg{Int,N}) where {N} = A.data[I...]
Base.setindex!(A::WrappedMPIArray, v, i::Int) = (A.data[i] = v)
Base.setindex!(A::WrappedMPIArray, v, I::Vararg{Int,N}) where {N} = (A.data[I...] = v)

function _mpi_parity_particles(xs, vx; L = 10.0)
    ps = ParticleSet{1,Float64}(length(xs))
    ps.x[1] .= xs
    ps.v[1] .= vx
    ps.v[2] .= 2 .* vx
    ps.v[3] .= -vx
    ps.weight .= 1.0 .+ 0.1 .* eachindex(xs)
    ps.id .= UInt64.(100 .+ eachindex(xs))
    ps.tag .= UInt32.(200 .+ eachindex(xs))
    return ps
end

@testset "MPI initialization and dimensions" begin
    ensure_mpi_initialized!()
    @test mpi_initialized()
    @test mpi_comm_size() >= 1
    @test mpi_comm_rank() >= 0
    @test prod(mpi_dims_create(1, 3)) == 1
    @test_throws ArgumentError mpi_dims_create(0, 2)
    @test_throws ArgumentError mpi_dims_create(1, 0)
end

@testset "MPI Cartesian communicator over logical layout" begin
    ensure_mpi_initialized!()
    layout = LogicalRankLayout((1, 1); periodic = (true, false))
    ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
    try
        @test ctx.layout == layout
        @test ctx.mpi_rank == 0
        @test ctx.mpi_size == 1
        @test ctx.coords == (1, 1)
        @test ctx.logical_rank == 1
        @test mpi_cartesian_neighbor(ctx, 1, 1) == 1
        @test mpi_cartesian_neighbor(ctx, 1, -1) == 1
        @test mpi_cartesian_neighbor(ctx, 2, 1) === nothing
        @test mpi_cartesian_neighbor(ctx, 2, -1) === nothing
        @test_throws ArgumentError mpi_cartesian_neighbor(ctx, 0, 1)

        desc = mpi_rank_layout_description(ctx)
        @test occursin("mpi;ranks=1", desc)
        @test occursin("dims=(1,1)", desc)
        @test occursin("coords=(0,0)", desc)
        @test occursin("periodic=(true,false)", desc)
        @test occursin("mpi_rank=0", desc)

        meta = capture_metadata(; rng_seed = 7, rank_layout = desc)
        @test meta.rank_layout == desc
    finally
        free_mpi_communicator!(ctx)
    end

    mismatch = LogicalRankLayout((2,); periodic = (true,))
    @test_throws ArgumentError create_cartesian_communicator(mismatch; comm = MPI.COMM_SELF)
end

@testset "one-rank MPI checkpoint restart" begin
    ensure_mpi_initialized!()
    @test MPI_CHECKPOINT_SCHEMA_VERSION == 1
    g = FourierGrid((8,), (2π,))
    layout = LogicalRankLayout((1,); periodic = (true,))
    ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
    try
        ps = ParticleSet{1,Float64}(8)
        load_lattice_1d!(ps, 0.0, 2π)
        set_density_weight!(ps, 1.0, g)
        st = HybridStepper(g, HybridModel(IsothermalElectrons(0.2)), CIC(), nparticles(ps))
        fill!(st.fields.B[3], 1.0)
        mpi_init!(st, ps, ctx)
        mpi_step!(st, ps, ctx, 0.01; NB = 2)

        mktempdir() do dir
            manifest_path = save_mpi_checkpoint(dir, st, ps, ctx)
            @test isfile(manifest_path)
            restored = HybridStepper(g, HybridModel(IsothermalElectrons(0.2)), CIC(), 1)
            ps_restored = ParticleSet{1,Float64}(1)
            load_mpi_checkpoint!(restored, ps_restored, dir, ctx)
            @test restored.step[] == st.step[]
            @test restored.time[] == st.time[]
            @test nparticles(ps_restored) == nparticles(ps)
            @test length(restored.work) == nparticles(ps)
            @test ps_restored.x[1] == ps.x[1]
            for c = 1:3
                @test ps_restored.v[c] == ps.v[c]
                @test restored.fields.B[c] == st.fields.B[c]
                @test restored.fields.E[c] == st.fields.E[c]
            end
            @test ps_restored.weight == ps.weight
            @test ps_restored.id == ps.id
            @test ps_restored.tag == ps.tag

            manifest = Serialization.deserialize(manifest_path)
            bad = merge(manifest, (schema = MPI_CHECKPOINT_SCHEMA_VERSION + 1,))
            Serialization.serialize(manifest_path, bad)
            poisoned = HybridStepper(g, HybridModel(IsothermalElectrons(0.2)), CIC(), 1)
            ps_poisoned = ParticleSet{1,Float64}(1)
            fill!(poisoned.fields.B[1], 7.0)
            @test_throws ArgumentError load_mpi_checkpoint!(poisoned, ps_poisoned, dir, ctx)
            @test nparticles(ps_poisoned) == 1
            @test all(==(7.0), poisoned.fields.B[1])
        end
    finally
        free_mpi_communicator!(ctx)
    end
end

@testset "one-rank MPI layout agrees with serial reference operations" begin
    ensure_mpi_initialized!()
    g = FourierGrid((8,), (10.0,))
    layout = LogicalRankLayout((1,); periodic = (true,))
    ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
    try
        @test ctx.logical_rank == 1
        @test rank_bounds(g, ctx.layout, ctx.logical_rank) == (lo = (0.0,), hi = (10.0,))

        field = [[-1.0, 10.0, 11.0, 12.0, 13.0, -2.0]]
        field_stats = exchange_field_halos!(field, ctx.layout; halo = 1)
        @test field_stats == (exchanged = 0, filled = 0)
        @test field[1] == [-1.0, 10.0, 11.0, 12.0, 13.0, -2.0]
        mpi_field = copy(field[1])
        mpi_field_stats = mpi_exchange_field_halos!(mpi_field, ctx; halo = 1, fill_value = -99.0)
        @test mpi_field_stats == (exchanged = 0, filled = 0)
        @test mpi_field == field[1]

        moments = [[2.0, 10.0, 11.0, 12.0, 3.0]]
        moment_stats = exchange_ghost_moments!(moments, ctx.layout; halo = 1)
        @test moment_stats == (exchanged = 0, dropped = 0)
        @test moments[1] == [2.0, 10.0, 11.0, 12.0, 3.0]
        mpi_moments = copy(moments[1])
        mpi_moment_stats = mpi_exchange_ghost_moments!(mpi_moments, ctx; halo = 1)
        @test mpi_moment_stats == (exchanged = 0, dropped = 0)
        @test mpi_moments == moments[1]

        serial = _mpi_parity_particles([9.9, 0.2, 5.0, 0.1], [1.0, -2.0, 0.5, 3.0])
        rank_particles = [_mpi_parity_particles([-0.1, 0.2, 5.0, 10.1], [1.0, -2.0, 0.5, 3.0])]
        migration = migrate_particles!(rank_particles, g, ctx.layout)
        @test migration == (moved = 0, lost = 0)
        @test rank_particles[1].x[1] ≈ serial.x[1]
        @test rank_particles[1].v == serial.v
        @test rank_particles[1].weight == serial.weight
        @test rank_particles[1].id == serial.id
        @test rank_particles[1].tag == serial.tag

        local_value =
            (number = sum(rank_particles[1].weight), momentum = total_momentum(rank_particles[1]))
        mpi_reduced = mpi_allreduce_diagnostics(local_value, ctx; op = :sum)
        serial_reduced = reduce_diagnostics([local_value]; op = :sum)
        @test mpi_reduced.number ≈ serial_reduced.number
        @test mpi_reduced.momentum[1] ≈ total_momentum(serial)[1]
        @test mpi_reduced.momentum[2] ≈ total_momentum(serial)[2]
        @test mpi_reduced.momentum[3] ≈ total_momentum(serial)[3]

        empty = ParticleSet{1,Float64}(0)
        f = HybridFields{1,Float64}((8,))
        fill!(f.n, 1.0)
        fill!(f.ui[1], 2.0)
        fill!(f.ui[2], 3.0)
        fill!(f.ui[3], 4.0)
        status = GPUAwareMPIStatus(false, false, false, :test, "host-only one-rank test")

        @test_throws ArgumentError mpi_compute_moments!(
            f,
            empty,
            g,
            NGP(),
            0.0,
            ctx;
            gpu_status = status,
        )
        @test all(==(1.0), f.n)
        @test all(==(2.0), f.ui[1])
        @test all(==(3.0), f.ui[2])
        @test all(==(4.0), f.ui[3])

        @test_throws ArgumentError mpi_compute_moments!(
            f,
            empty,
            g,
            NGP(),
            NaN,
            ctx;
            gpu_status = status,
        )
        @test all(==(1.0), f.n)
        @test all(==(2.0), f.ui[1])
        @test all(==(3.0), f.ui[2])
        @test all(==(4.0), f.ui[3])
    finally
        free_mpi_communicator!(ctx)
    end
end

@testset "GPU-aware MPI status and host staging fallback" begin
    ensure_mpi_initialized!()
    status = gpu_aware_mpi_status()
    @test status isa GPUAwareMPIStatus
    @test status.enabled == (status.cuda || status.rocm)
    @test status.source === :mpi
    @test occursin("cuda=", status.reason)
    @test occursin("rocm=", status.reason)

    withenv("JULIA_MPI_HAS_CUDA" => "true", "JULIA_MPI_HAS_ROCM" => "false") do
        forced = gpu_aware_mpi_status(; initialize = false)
        @test forced.enabled
        @test forced.cuda
        @test !forced.rocm
    end

    host = [1.0, 2.0, 3.0]
    @test !mpi_buffer_uses_host_staging(host; status)
    @test host_staging_buffer(host) === host

    wrapped = WrappedMPIArray(copy(host))
    @test mpi_buffer_uses_host_staging(wrapped; status)

    send_plan = prepare_mpi_buffer(wrapped; status, intent = :send)
    @test send_plan.used_host_staging
    @test !send_plan.copy_back
    @test send_plan.buffer == host
    send_plan.buffer[1] = 9.0
    finish_mpi_buffer!(send_plan)
    @test wrapped[1] == 1.0

    recv_plan = prepare_mpi_buffer(wrapped; status, intent = :recv)
    @test recv_plan.used_host_staging
    @test recv_plan.copy_back
    recv_plan.buffer .= [4.0, 5.0, 6.0]
    @test finish_mpi_buffer!(recv_plan) === wrapped
    @test collect(wrapped) == [4.0, 5.0, 6.0]

    @test_throws ArgumentError prepare_mpi_buffer(wrapped; status, intent = :bad)
    @test_throws DimensionMismatch copy_from_host_staging!(wrapped, ones(2))
end

@testset "one-rank MPI diagnostics agree with serial reference" begin
    ensure_mpi_initialized!()
    layout = LogicalRankLayout((1,); periodic = (true,))
    ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
    try
        local_value = (energy = [1.0, 2.0, 3.0], extrema = (lo = -2.0, hi = 5.0), count = 4)

        @test mpi_allreduce_diagnostics(local_value, ctx; op = :sum) ==
              reduce_diagnostics([local_value]; op = :sum)
        @test mpi_allreduce_diagnostics(local_value, ctx; op = :min) ==
              reduce_diagnostics([local_value]; op = :min)
        @test mpi_allreduce_diagnostics(local_value, ctx; op = :max) ==
              reduce_diagnostics([local_value]; op = :max)
        @test_throws ArgumentError mpi_allreduce_diagnostics("bad leaf", ctx; op = :sum)
        @test_throws ArgumentError mpi_allreduce_diagnostics(1.0, ctx; op = :prod)

        wrapped = WrappedMPIArray([2.0, 3.0])
        reduced = mpi_allreduce_diagnostics((field = wrapped,), ctx; op = :sum)
        @test reduced.field isa Vector{Float64}
        @test reduced.field == [2.0, 3.0]
    finally
        free_mpi_communicator!(ctx)
    end
end

@testset "mpi_step! validates timestep and subcycles before mutation" begin
    ensure_mpi_initialized!()
    g = FourierGrid((4,), (2π,))
    layout = LogicalRankLayout((1,); periodic = (true,))
    ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
    try
        T = Float64
        ps = ParticleSet{1,T}(4)
        load_lattice_1d!(ps, 0.0, 2π)
        set_density_weight!(ps, 1.0, g)
        st = HybridStepper(g, HybridModel(IsothermalElectrons(0.0)), NGP(), nparticles(ps))
        fill!(st.fields.B[3], 1.0)
        init!(st, ps)
        status = gpu_aware_mpi_status()

        x0 = ntuple(d -> copy(ps.x[d]), 1)
        v0 = ntuple(c -> copy(ps.v[c]), 3)
        B0 = ntuple(c -> copy(st.fields.B[c]), 3)
        time0 = st.time[]
        step0 = st.step[]
        work_len0 = length(st.work)

        @test_throws ArgumentError mpi_step!(st, ps, ctx, 0.1; NB = 0, gpu_status = status)
        @test_throws ArgumentError mpi_step!(st, ps, ctx, NaN; NB = 1, gpu_status = status)
        @test_throws ArgumentError mpi_step!(st, ps, ctx, -0.1; NB = 1, gpu_status = status)
        @test st.time[] == time0
        @test st.step[] == step0
        @test length(st.work) == work_len0
        @test all(ps.x[d] == x0[d] for d = 1:1)
        @test all(ps.v[c] == v0[c] for c = 1:3)
        @test all(st.fields.B[c] == B0[c] for c = 1:3)
    finally
        free_mpi_communicator!(ctx)
    end
end
