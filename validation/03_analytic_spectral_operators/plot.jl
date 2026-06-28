#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_03_analytic_spectral_operators(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "03_analytic_spectral_operators.csv",
        "03_analytic_spectral_operators.pdf";
        title = "Spectral operator analytic validation",
    )
end

VALIDATION_PLOT = plot_03_analytic_spectral_operators

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
