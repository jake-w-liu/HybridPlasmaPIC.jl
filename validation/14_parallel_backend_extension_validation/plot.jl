#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_14_parallel_backend_extension_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "14_parallel_backend_extension_validation.csv",
        "14_parallel_backend_extension_validation.pdf";
        title = "Parallel, backend, and extension validation",
    )
end

VALIDATION_PLOT = plot_14_parallel_backend_extension_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
