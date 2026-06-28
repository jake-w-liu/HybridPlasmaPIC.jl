#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_21_published_preisser2020_summary(artifact_dir::AbstractString)
    rows = _read_csv(joinpath(artifact_dir, "21_published_preisser2020_summary.csv"))
    rows = [row for row in rows if startswith(row["metric"], "Bavg_y_")]
    isempty(rows) && return nothing
    x = [replace(row["metric"], "Bavg_y_" => "") for row in rows]
    measured = [_num(row["measured"]) for row in rows]
    expected = [_num(row["expected"]) for row in rows]
    fig = PlotlySupply.plot_bar(
        x,
        [measured, expected];
        title = "Preisser 2020 published hybrid-code summary",
        xlabel = "Bavg_y statistic",
        ylabel = "value",
        legend = ["measured", "published summary"],
        show = false,
    )
    return _save_pdf(joinpath(artifact_dir, "21_published_preisser2020_summary.pdf"), fig)
end

VALIDATION_PLOT = plot_21_published_preisser2020_summary

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_plot_main(VALIDATION_PLOT, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
