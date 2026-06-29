#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_06_boundary_loading_kdv_smoothing(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "06_boundary_loading_kdv_smoothing.csv",
        "06_boundary_loading_kdv_smoothing.pdf";
        title = "Boundary, loading, smoothing, and KdV validation",
    )
end

VALIDATION_PLOT = plot_06_boundary_loading_kdv_smoothing

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
