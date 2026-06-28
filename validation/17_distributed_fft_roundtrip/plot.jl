#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_17_distributed_fft_roundtrip(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "17_distributed_fft_roundtrip.csv",
        "17_distributed_fft_roundtrip.pdf";
        title = "Distributed FFT extension validation",
    )
end

VALIDATION_PLOT = plot_17_distributed_fft_roundtrip

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
