#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_12_model_runtime_smoke_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "12_model_runtime_smoke_validation.csv",
        "12_model_runtime_smoke_validation.pdf";
        title = "Model runtime smoke validation",
    )
end

VALIDATION_PLOT = plot_12_model_runtime_smoke_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
