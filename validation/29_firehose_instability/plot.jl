#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_29_firehose_instability(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "29_firehose_instability.csv",
        "29_firehose_instability.pdf";
        title = "Parallel firehose instability validation",
    )
end

VALIDATION_PLOT = plot_29_firehose_instability

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
