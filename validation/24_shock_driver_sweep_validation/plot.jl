#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_24_shock_driver_sweep_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "24_shock_driver_sweep_validation.csv",
        "24_shock_driver_sweep_validation.pdf";
        title = "Shock driver sweep validation",
    )
end

VALIDATION_PLOT = plot_24_shock_driver_sweep_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
