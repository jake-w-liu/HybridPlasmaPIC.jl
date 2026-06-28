#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_16_pic1d_hdf5_extension_validation(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "16_pic1d_hdf5_extension_validation.csv",
        "16_pic1d_hdf5_extension_validation.pdf";
        title = "1D PIC and HDF5 validation",
    )
end

VALIDATION_PLOT = plot_16_pic1d_hdf5_extension_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
