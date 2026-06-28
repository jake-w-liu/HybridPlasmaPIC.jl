#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_01_analytic_poisson_2d(artifact_dir::AbstractString)
    return _line_compare_plot(
        artifact_dir,
        "01_analytic_poisson_2d_slice.csv",
        "01_analytic_poisson_2d_slice.pdf";
        title = "Electrostatic Poisson single-mode field",
        xcol = "x",
        measured_expected_pairs = (("measured_Ex", "expected_Ex", "Ex"), ("measured_Ey", "expected_Ey", "Ey")),
        yaxis_title = "electric field",
    )
end

VALIDATION_PLOT = plot_01_analytic_poisson_2d

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
