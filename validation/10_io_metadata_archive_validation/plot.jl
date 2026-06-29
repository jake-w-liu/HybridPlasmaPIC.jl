#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_10_io_metadata_archive_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "10_io_metadata_archive_validation.csv",
        "10_io_metadata_archive_validation.pdf";
        title = "Metadata, archive, and IO validation",
    )
end

VALIDATION_PLOT = plot_10_io_metadata_archive_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
