#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_18_hybrid_ion_acoustic_dispersion(artifact_dir::AbstractString)
    return _time_series_plot(
        artifact_dir,
        "18_hybrid_ion_acoustic_dispersion.csv",
        "18_hybrid_ion_acoustic_dispersion.pdf";
        title = "Hybrid PIC ion-acoustic mode history",
        ycols = ("density_mode_real",),
        yaxis_title = "Re density mode",
    )
end

VALIDATION_PLOT = plot_18_hybrid_ion_acoustic_dispersion

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
