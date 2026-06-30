#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

import MPI
import Serialization

struct ValidationWrappedMPIArray{T,N,A<:Array{T,N}} <: AbstractArray{T,N}
    data::A
end

Base.size(A::ValidationWrappedMPIArray) = size(A.data)
Base.axes(A::ValidationWrappedMPIArray) = axes(A.data)
Base.IndexStyle(::Type{<:ValidationWrappedMPIArray}) = IndexLinear()
Base.getindex(A::ValidationWrappedMPIArray, i::Int) = A.data[i]
Base.getindex(A::ValidationWrappedMPIArray, I::Vararg{Int,N}) where {N} = A.data[I...]
Base.setindex!(A::ValidationWrappedMPIArray, v, i::Int) = (A.data[i] = v)
Base.setindex!(A::ValidationWrappedMPIArray, v, I::Vararg{Int,N}) where {N} = (A.data[I...] = v)

function _max_field_error(a::HybridFields{D,T}, b::HybridFields{D,T}) where {D,T}
    err = maximum(abs, a.n .- b.n)
    err = max(err, maximum(abs, a.pe .- b.pe))
    err = max(err, maximum(abs, a.ninv .- b.ninv))
    for d = 1:D
        err = max(err, maximum(abs, a.gradp[d] .- b.gradp[d]))
    end
    for c = 1:3
        err = max(err, maximum(abs, a.ui[c] .- b.ui[c]))
        err = max(err, maximum(abs, a.B[c] .- b.B[c]))
        err = max(err, maximum(abs, a.E[c] .- b.E[c]))
        err = max(err, maximum(abs, a.J[c] .- b.J[c]))
    end
    return err
end

function _max_checkpoint_field_error(a::HybridFields, b::HybridFields)
    err = 0.0
    for c = 1:3
        err = max(err, maximum(abs, a.B[c] .- b.B[c]))
        err = max(err, maximum(abs, a.E[c] .- b.E[c]))
    end
    return err
end

function _max_particle_error(a::ParticleSet{D,T}, b::ParticleSet{D,T}) where {D,T}
    err = Float64(nparticles(a) == nparticles(b) && a.q == b.q && a.m == b.m ? 0.0 : 1.0)
    for d = 1:D
        err = max(err, maximum(abs, a.x[d] .- b.x[d]))
    end
    for c = 1:3
        err = max(err, maximum(abs, a.v[c] .- b.v[c]))
    end
    err = max(err, maximum(abs, a.weight .- b.weight))
    err = max(err, a.id == b.id ? 0.0 : 1.0)
    err = max(err, a.tag == b.tag ? 0.0 : 1.0)
    return err
end

function case_15_mpi_single_rank_validation(artifact_dir::AbstractString)
    id = "15_mpi_single_rank_validation"
    ensure_mpi_initialized!()
    gpu_status =
        GPUAwareMPIStatus(false, false, false, :validation, "host-only single-rank validation")
    queried_gpu_status = gpu_aware_mpi_status(; initialize = false)

    init_error =
        mpi_initialized() &&
        mpi_comm_size(MPI.COMM_WORLD) >= 1 &&
        mpi_comm_rank(MPI.COMM_WORLD) >= 0 &&
        prod(mpi_dims_create(1, 3)) == 1 ? 0.0 : 1.0

    layout = LogicalRankLayout((1,); periodic = (true,))
    ctx = create_cartesian_communicator(layout; comm = MPI.COMM_SELF)
    try
        rank_layout_error =
            ctx isa MPICartesianCommunicator &&
            ctx.mpi_rank == 0 &&
            ctx.mpi_size == 1 &&
            ctx.coords == (1,) &&
            ctx.logical_rank == 1 &&
            mpi_cartesian_neighbor(ctx, 1, 1) == 1 &&
            occursin("mpi;ranks=1", mpi_rank_layout_description(ctx)) ? 0.0 : 1.0

        local_value = (energy = [1.0, 2.0, 3.0], extrema = (lo = -2.0, hi = 5.0), count = 4)
        sum_error = maximum(
            abs,
            mpi_allreduce_diagnostics(local_value, ctx; op = :sum, gpu_status).energy .-
            local_value.energy,
        )
        min_value = mpi_allreduce_diagnostics(local_value, ctx; op = :min, gpu_status)
        max_value = mpi_allreduce_diagnostics(local_value, ctx; op = :max, gpu_status)
        allreduce_error =
            max(sum_error, min_value.extrema.lo == -2.0 && max_value.extrema.hi == 5.0 ? 0.0 : 1.0)

        wrapped = ValidationWrappedMPIArray([1.0, 2.0, 3.0])
        plan = prepare_mpi_buffer(wrapped; status = gpu_status, intent = :recv)
        plan.buffer .= [4.0, 5.0, 6.0]
        finish_mpi_buffer!(plan)
        direct_copy = ValidationWrappedMPIArray([0.0, 0.0, 0.0])
        copy_from_host_staging!(direct_copy, [7.0, 8.0, 9.0])
        gpu_status_query_error =
            queried_gpu_status isa GPUAwareMPIStatus &&
            queried_gpu_status.enabled == (queried_gpu_status.cuda || queried_gpu_status.rocm) &&
            queried_gpu_status.source === :mpi ? 0.0 : 1.0
        direct_copy_error = direct_copy.data == [7.0, 8.0, 9.0] ? 0.0 : 1.0
        staging_error = max(
            plan isa MPIBufferPlan &&
            plan.used_host_staging &&
            plan.copy_back &&
            wrapped.data == [4.0, 5.0, 6.0] &&
            !mpi_buffer_uses_host_staging([1.0]; status = gpu_status) &&
            host_staging_buffer([1.0, 2.0]) == [1.0, 2.0] ? 0.0 : 1.0,
            gpu_status_query_error,
            direct_copy_error,
        )

        g = FourierGrid((8,), (2π,))
        ps = ParticleSet{1,Float64}(8)
        load_lattice_1d!(ps, 0.0, 2π)
        set_density_weight!(ps, 1.0, g)
        for p = 1:nparticles(ps)
            ps.v[1][p] = 0.1 * sin(ps.x[1][p])
            ps.v[2][p] = 0.2 * cos(ps.x[1][p])
            ps.v[3][p] = -0.05 * sin(2 * ps.x[1][p])
        end
        serial_fields = HybridFields{1,Float64}(g.n)
        mpi_fields = HybridFields{1,Float64}(g.n)
        compute_moments!(serial_fields, ps, g, CIC(), 1e-6)
        mpi_compute_moments!(mpi_fields, ps, g, CIC(), 1e-6, ctx; gpu_status)
        moments_error = _max_field_error(serial_fields, mpi_fields)

        rank_particles = ParticleSet{1,Float64}(4)
        rank_particles.x[1] .= [-0.1, 0.2, 5.0, 10.1]
        rank_particles.v[1] .= [1.0, -2.0, 0.5, 3.0]
        rank_particles.v[2] .= 2 .* rank_particles.v[1]
        rank_particles.v[3] .= -rank_particles.v[1]
        rank_particles.weight .= [1.1, 1.2, 1.3, 1.4]
        rank_particles.id .= UInt64[101, 102, 103, 104]
        rank_particles.tag .= UInt32[201, 202, 203, 204]
        migration_stats = mpi_migrate_particles!(rank_particles, FourierGrid((8,), (10.0,)), ctx)
        migration_error =
            migration_stats == (moved = 0, lost = 0, sent = 0, received = 0) &&
            maximum(abs, rank_particles.x[1] .- [9.9, 0.2, 5.0, 0.1]) < 1e-12 ? 0.0 : 1.0

        field = [-1.0, 10.0, 11.0, 12.0, 13.0, -2.0]
        field_stats = mpi_exchange_field_halos!(field, ctx; halo = 1, fill_value = -99.0)
        moment = [2.0, 10.0, 11.0, 12.0, 3.0]
        moment_stats = mpi_exchange_ghost_moments!(moment, ctx; halo = 1)
        halo_error =
            field_stats == (exchanged = 0, filled = 0) &&
            moment_stats == (exchanged = 0, dropped = 0) &&
            field == [-1.0, 10.0, 11.0, 12.0, 13.0, -2.0] &&
            moment == [2.0, 10.0, 11.0, 12.0, 3.0] ? 0.0 : 1.0

        st = HybridStepper(g, HybridModel(IsothermalElectrons(0.2)), CIC(), nparticles(ps))
        fill!(st.fields.B[3], 1.0)
        mpi_init!(st, ps, ctx; gpu_status)
        mpi_step!(st, ps, ctx, 0.01; NB = 2, gpu_status, migrate_particles = true)
        restored = HybridStepper(g, HybridModel(IsothermalElectrons(0.2)), CIC(), 1)
        ps_restored = ParticleSet{1,Float64}(1)
        checkpoint_error = 0.0
        mktempdir() do dir
            manifest_path = save_mpi_checkpoint(dir, st, ps, ctx)
            isfile(manifest_path) || (checkpoint_error = 1.0)
            manifest = Serialization.deserialize(manifest_path)
            manifest.schema == MPI_CHECKPOINT_SCHEMA_VERSION || (checkpoint_error = 1.0)
            load_mpi_checkpoint!(restored, ps_restored, dir, ctx)
        end
        checkpoint_error = max(
            checkpoint_error,
            restored.step[] == st.step[] && restored.time[] == st.time[] ? 0.0 : 1.0,
            _max_checkpoint_field_error(restored.fields, st.fields),
            _max_particle_error(ps_restored, ps),
        )

        artifact = joinpath(artifact_dir, "15_mpi_single_rank_validation.csv")
        rows = (
            ("mpi_initialization_contract_error", init_error, 0.0, "absolute", init_error, 0.0),
            (
                "mpi_cartesian_rank_layout_error",
                rank_layout_error,
                0.0,
                "absolute",
                rank_layout_error,
                0.0,
            ),
            (
                "mpi_allreduce_nested_contract_error",
                allreduce_error,
                0.0,
                "absolute",
                allreduce_error,
                0.0,
            ),
            ("mpi_host_staging_contract_error", staging_error, 0.0, "absolute", staging_error, 0.0),
            (
                "mpi_compute_moments_serial_max_abs_error",
                moments_error,
                0.0,
                "absolute",
                moments_error,
                1e-12,
            ),
            (
                "mpi_single_rank_migration_error",
                migration_error,
                0.0,
                "absolute",
                migration_error,
                0.0,
            ),
            ("mpi_single_rank_halo_exchange_error", halo_error, 0.0, "absolute", halo_error, 0.0),
            (
                "mpi_checkpoint_roundtrip_error",
                checkpoint_error,
                0.0,
                "absolute",
                checkpoint_error,
                1e-12,
            ),
        )
        _write_metric_csv(artifact, rows)
        return _metric_rows_to_results(
            id = id,
            category = "mpi_single_rank",
            reference_kind = "external_library",
            reference = "MPI.jl COMM_SELF transport and serial reference invariants",
            rows = rows,
            artifact = artifact,
        )
    finally
        free_mpi_communicator!(ctx)
    end
end

VALIDATION_CASE = ValidationCase(
    id = "15_mpi_single_rank_validation",
    default = true,
    description = "MPI.jl single-rank transport, reductions, staging, moments, migration, halos, and checkpoint roundtrip.",
    runner = case_15_mpi_single_rank_validation,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_case_main(
            VALIDATION_CASE,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
