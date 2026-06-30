#!/usr/bin/env julia
#
# Overlay figure for case 27: our Hall-MHD whistler ω(k) (solid) vs the external
# NHDS kinetic solver (dashed). "Lines near overlay" = the dispersions agree.
# Reads the overlay CSV written by run.jl; PDF is gitignored.

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_27_nhds_dispersion_comparison(artifact_dir::AbstractString)
    return _line_compare_plot(
        artifact_dir,
        "27_nhds_dispersion_comparison_overlay.csv",
        "27_nhds_dispersion_comparison.pdf";
        title = "Whistler ω(k): our Hall-MHD oracle vs NHDS (kinetic, Vlasov-Maxwell)",
        xcol = "k_dp",
        measured_expected_pairs = [("omega_oracle_hallmhd", "omega_NHDS_kinetic", "whistler ω(k)")],
        yaxis_title = "ω / Ω_p",
    )
end

VALIDATION_PLOT = plot_27_nhds_dispersion_comparison

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
