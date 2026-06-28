#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_22_shock_spacecraft_diagnostics(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "22_shock_spacecraft_diagnostics.csv",
        "22_shock_spacecraft_diagnostics.pdf";
        title = "Shock and spacecraft diagnostics validation",
    )
end

VALIDATION_PLOT = plot_22_shock_spacecraft_diagnostics

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
