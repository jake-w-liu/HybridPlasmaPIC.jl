#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_32_magnetic_reconnection(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "32_magnetic_reconnection.csv",
        "32_magnetic_reconnection.pdf";
        title = "Harris-sheet reconnection validation",
    )
end

VALIDATION_PLOT = plot_32_magnetic_reconnection

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
