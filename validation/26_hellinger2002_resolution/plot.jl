#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_26_hellinger2002_resolution(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "26_hellinger2002_resolution.csv",
        "26_hellinger2002_resolution.pdf";
        title = "Hellinger 2002 shock-resolution validation",
    )
end

VALIDATION_PLOT = plot_26_hellinger2002_resolution

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
