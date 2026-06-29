#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_15_mpi_single_rank_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "15_mpi_single_rank_validation.csv",
        "15_mpi_single_rank_validation.pdf";
        title = "MPI single-rank validation",
    )
end

VALIDATION_PLOT = plot_15_mpi_single_rank_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
