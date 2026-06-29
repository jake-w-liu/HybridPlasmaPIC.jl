#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_13_threaded_backend_api_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "13_threaded_backend_api_validation.csv",
        "13_threaded_backend_api_validation.pdf";
        title = "Threaded backend API validation",
    )
end

VALIDATION_PLOT = plot_13_threaded_backend_api_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
