#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_20_empic_transverse_dispersion_2d(artifact_dir::AbstractString)
    return _time_series_plot(
        artifact_dir,
        "20_empic_transverse_dispersion_2d.csv",
        "20_empic_transverse_dispersion_2d.pdf";
        title = "2D EM PIC transverse mode history",
        ycols = ("mode_real", "mode_imag", "mode_abs"),
        yaxis_title = "mode amplitude",
    )
end

VALIDATION_PLOT = plot_20_empic_transverse_dispersion_2d

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
