#!/usr/bin/env julia
#
# Plot the Hybrid-VPIC perpendicular-shock Bz(x) profile (the external code's shock
# structure). Reads the profile CSV written by run.jl; PDF is gitignored.

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_28_hybridvpic_perp_shock(artifact_dir::AbstractString)
    rows = _read_csv(joinpath(artifact_dir, "28_hybridvpic_perp_shock_vpic_profile.csv"))
    isempty(rows) && return nothing
    x = [_num(r["x_di"]) for r in rows]
    bz = [_num(r["Bz_vpic"]) for r in rows]
    fig = PlotlySupply.plot_scatter(
        x,
        [bz];
        title = "Hybrid-VPIC perpendicular shock Bz(x)  (β=1, M_A≈6) — compression+overshoot vs ours",
        xlabel = "x / d_i  (wall/downstream at left, upstream at right)",
        ylabel = "Bz / B0",
        legend = ["Hybrid-VPIC"],
        show = false,
    )
    return _save_pdf(joinpath(artifact_dir, "28_hybridvpic_perp_shock.pdf"), fig)
end

VALIDATION_PLOT = plot_28_hybridvpic_perp_shock

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
