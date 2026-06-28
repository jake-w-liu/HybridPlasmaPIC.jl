#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_02_analytic_hall_mhd_continuity(artifact_dir::AbstractString)
    return _line_compare_plot(
        artifact_dir,
        "02_analytic_hall_mhd_continuity.csv",
        "02_analytic_hall_mhd_continuity.pdf";
        title = "Hall-MHD continuity RHS",
        xcol = "x",
        measured_expected_pairs = (("measured_dn_dt", "expected_dn_dt", "dn/dt"),),
        yaxis_title = "dn/dt",
    )
end

VALIDATION_PLOT = plot_02_analytic_hall_mhd_continuity

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
