#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_19_integrator_camcl_semiimplicit(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "19_integrator_camcl_semiimplicit.csv",
        "19_integrator_camcl_semiimplicit.pdf";
        title = "CAM-CL and semi-implicit integrator validation",
    )
end

VALIDATION_PLOT = plot_19_integrator_camcl_semiimplicit

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
