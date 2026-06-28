#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_07_closure_budget_filter_diagnostics(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "07_closure_budget_filter_diagnostics.csv",
        "07_closure_budget_filter_diagnostics.pdf";
        title = "Closure, budget, and diagnostic validation",
    )
end

VALIDATION_PLOT = plot_07_closure_budget_filter_diagnostics

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
