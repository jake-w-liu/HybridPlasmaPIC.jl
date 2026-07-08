#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_25_leroy1982_perp_shock(artifact_dir::AbstractString)
    return _metric_plot(
        artifact_dir,
        "25_leroy1982_perp_shock.csv",
        "25_leroy1982_perp_shock.pdf";
        title = "Leroy 1982 perpendicular-shock validation",
    )
end

VALIDATION_PLOT = plot_25_leroy1982_perp_shock

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
