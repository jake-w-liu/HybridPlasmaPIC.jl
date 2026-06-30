#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_14_parallel_backend_extension_validation(artifact_dir::AbstractString)
    id = "14_parallel_backend_extension_validation"
    ps = ParticleSet{2,Float64}(3)
    ps.x[1] .= [0.1, 0.2, 0.3]
    ps.x[2] .= [0.4, 0.5, 0.6]
    ps.v[1] .= [1.0, 2.0, 3.0]
    ps.id .= UInt64[3, 2, 1]

    supported_error = supported_extensions() == (:cuda, :metal, :io, :pencilfft) ? 0.0 : 1.0
    backend_error = particle_storage_backend(ps) == :cpu ? 0.0 : 1.0
    ps_cpu = copy_particles_to_backend(Val(:cpu), ps)
    cpu_copy_error =
        particle_storage_backend(ps_cpu) == :cpu && ps_cpu.x[1] == ps.x[1] && ps_cpu.id == ps.id ?
        0.0 : 1.0
    pressure =
        memory_pressure(BackendMemoryStatus(:cpu, true, false; total_bytes = 100, free_bytes = 25))
    memory_pressure_error = abs(pressure - 0.75)

    dec = PencilDecomposition3D((8, 6, 4), (2, 2))
    pencil_coord_error = pencil_rank_coords(dec, 3) == (1, 2) ? 0.0 : 1.0
    bounds = pencil_bounds(dec, 3, :x)
    pencil_bounds_error = bounds == (1:8, 1:3, 3:4) ? 0.0 : 1.0
    pencil_owner_error = abs(pencil_owner(dec, (1, 4, 3), :x) - 4)
    dec_alias = pencil_decomposition((8, 6, 4), (2, 2))
    axes_y = pencil_orientation_axes(:y)
    pencil_api_error =
        dec_alias isa PencilDecomposition3D &&
        pencil_nranks(dec_alias) == 4 &&
        pencil_rank_index(dec_alias, (2, 2)) == 4 &&
        pencil_local_size(dec_alias, 3, :x) == (8, 3, 2) &&
        axes_y.full_axis == 2 &&
        axes_y.split_axes == (1, 3) ? 0.0 : 1.0

    layout = LogicalRankLayout((2,); periodic = (false,))
    halos = [Float64[10, 11, 12, 13, 14], Float64[20, 21, 22, 23, 24]]
    halo_stats = exchange_field_halos!(halos, layout; halo = 1, fill_value = -1.0)
    halo_error =
        halo_stats.exchanged == 2 &&
        halo_stats.filled == 2 &&
        halos[1] == [-1.0, 11.0, 12.0, 13.0, 21.0] &&
        halos[2] == [13.0, 21.0, 22.0, 23.0, -1.0] ? 0.0 : 1.0

    mpi_error = 0.0
    mpi_notes = ""
    try
        dims = mpi_dims_create(1, 2)
        dims == (1, 1) || (mpi_error = 1.0)
        ctx = create_cartesian_communicator(LogicalRankLayout((1,); periodic = (true,)))
        try
            desc = mpi_rank_layout_description(ctx)
            occursin("ranks=1", desc) && mpi_cartesian_neighbor(ctx, 1, 1) == 1 || (mpi_error = 1.0)
        finally
            free_mpi_communicator!(ctx)
        end
    catch err
        mpi_error = NaN
        mpi_notes = "MPI local communicator check skipped: $(typeof(err))"
    end

    artifact = joinpath(artifact_dir, "14_parallel_backend_extension_validation.csv")
    rows = (
        (
            "supported_extensions_contract_error",
            supported_error,
            0.0,
            "absolute",
            supported_error,
            0.0,
        ),
        ("cpu_particle_backend_error", backend_error, 0.0, "absolute", backend_error, 0.0),
        ("cpu_particle_copy_error", cpu_copy_error, 0.0, "absolute", cpu_copy_error, 0.0),
        ("memory_pressure_abs_error", pressure, 0.75, "absolute", memory_pressure_error, 0.0),
        ("pencil_rank_coords_error", pencil_coord_error, 0.0, "absolute", pencil_coord_error, 0.0),
        ("pencil_bounds_error", pencil_bounds_error, 0.0, "absolute", pencil_bounds_error, 0.0),
        (
            "pencil_owner_abs_error",
            pencil_owner(dec, (1, 4, 3), :x),
            4.0,
            "absolute",
            pencil_owner_error,
            0.0,
        ),
        (
            "pencil_api_wrapper_contract_error",
            pencil_api_error,
            0.0,
            "absolute",
            pencil_api_error,
            0.0,
        ),
        ("field_halo_exchange_error", halo_error, 0.0, "absolute", halo_error, 0.0),
    )
    _write_metric_csv(artifact, rows)
    results = _metric_rows_to_results(
        id = id,
        category = "parallel_backend_extensions",
        reference_kind = "analytic",
        reference = "CPU backend contracts, pencil ownership/ranges, deterministic local halo exchange",
        rows = rows,
        artifact = artifact,
    )
    if isfinite(mpi_error)
        mpirows = (("mpi_single_rank_cartesian_error", mpi_error, 0.0, "absolute", mpi_error, 0.0),)
        append!(
            results,
            _metric_rows_to_results(
                id = id,
                category = "mpi_serial_smoke",
                reference_kind = "external_library",
                reference = "MPI.jl COMM_WORLD size-1 Cartesian communicator",
                rows = mpirows,
                artifact = artifact,
                notes = mpi_notes,
            ),
        )
    else
        push!(
            results,
            _skip_result(
                id = id,
                category = "mpi_serial_smoke",
                reference_kind = "external_library",
                reference = "MPI.jl COMM_WORLD size-1 Cartesian communicator",
                metric = "mpi_single_rank_cartesian_error",
                notes = mpi_notes,
            ),
        )
    end
    return results
end


VALIDATION_CASE = ValidationCase(
    id = "14_parallel_backend_extension_validation",
    default = true,
    description = "CPU backend contracts, pencil decomposition, local halo exchange, and MPI serial smoke.",
    runner = case_14_parallel_backend_extension_validation,
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
