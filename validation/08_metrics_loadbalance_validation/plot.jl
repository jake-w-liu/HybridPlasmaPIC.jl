#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_08_metrics_loadbalance_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "08_metrics_loadbalance_validation.csv",
        "08_metrics_loadbalance_validation.pdf";
        title = "Closure and load-balance validation",
    )
end

VALIDATION_PLOT = plot_08_metrics_loadbalance_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
