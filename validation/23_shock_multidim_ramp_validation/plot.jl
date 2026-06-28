#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_23_shock_multidim_ramp_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "23_shock_multidim_ramp_validation.csv",
        "23_shock_multidim_ramp_validation.pdf";
        title = "Multidimensional shock and ramp validation",
    )
end

VALIDATION_PLOT = plot_23_shock_multidim_ramp_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
