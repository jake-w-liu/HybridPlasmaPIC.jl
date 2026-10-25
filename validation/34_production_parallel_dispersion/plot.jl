#!/usr/bin/env julia

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_34_production_parallel_dispersion(artifact_dir::AbstractString)
    id = "34_production_parallel_dispersion"
    fits = _read_csv(joinpath(artifact_dir, "$(id)_fits.csv"))
    isempty(fits) && return nothing
    outputs = String[]
    indices = collect(eachindex(fits))
    _write_label_key(
        joinpath(artifact_dir, "$(id)_plot_labels.csv"),
        indices,
        [row["label"] for row in fits],
    )
    for (metric, title, ylabel) in (
        ("omega", "Signed frequency from production wave histories", "omega"),
        ("gamma", "Growth and damping from production wave histories", "gamma"),
    )
        fig = PlotlySupply.plot_scatter(
            indices,
            [[_num(row[metric]) for row in fits], [_num(row[metric*"_reference"]) for row in fits]];
            title,
            xlabel = "profile index (see plot-label CSV)",
            ylabel,
            legend = ["measured", "continuum reference"],
            dash = ["", "dash"],
            show = false,
        )
        push!(outputs, _save_pdf(joinpath(artifact_dir, "$(id)_$(metric).pdf"), fig))
    end
    rows = _read_csv(joinpath(artifact_dir, "$(id)_orders.csv"))
    if !isempty(rows)
        groups = [[row for row in rows if row["solver"] == solver] for solver in ("hall", "hybrid")]
        x = [log10(_num(row["dt"])) for row in groups[1]]
        values = [[log10(_num(row["error"])) for row in group] for group in groups]
        fig = PlotlySupply.plot_scatter(
            x,
            values;
            title = "Temporal error against independent linear evolution",
            xlabel = "log10(timestep)",
            ylabel = "log10(relative maximum history error)",
            legend = ["Hall RK4", "hybrid"],
            show = false,
        )
        push!(outputs, _save_pdf(joinpath(artifact_dir, "$(id)_orders.pdf"), fig))
    end
    history = _time_series_plot(
        artifact_dir,
        "$(id)_hybrid_eta0.0_ioncyclotron_m2_baseline.csv",
        "$(id)_history.pdf";
        title = "Hybrid ion-cyclotron mode 2: measured baseline history",
        ycols = ("mode_real", "mode_imag", "mode_abs"),
        yaxis_title = "complex Fourier coefficient",
    )
    history === nothing || push!(outputs, history)
    return outputs
end

VALIDATION_PLOT = plot_34_production_parallel_dispersion

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
