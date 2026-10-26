#!/usr/bin/env julia
if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function plot_35_production_oblique_hall(artifact_dir::AbstractString)
    id = "35_production_oblique_hall"
    fits = _read_csv(joinpath(artifact_dir, "$(id)_fits.csv"))
    isempty(fits) && throw(ArgumentError("oblique fit artifact is missing or empty"))
    indices = collect(eachindex(fits))
    outputs = String[]
    _write_label_key(
        joinpath(artifact_dir, "$(id)_plot_labels.csv"),
        indices,
        [r["label"] for r in fits],
    )
    for (metric, title, ylabel) in (
        ("omega", "Signed oblique Hall frequency", "omega"),
        ("gamma", "Numerical oblique Hall growth/damping", "gamma"),
    )
        fig = PlotlySupply.plot_scatter(
            indices,
            [[_num(r[metric]) for r in fits], [_num(r[metric*"_reference"]) for r in fits]];
            title,
            xlabel = "profile index (see plot-label CSV)",
            ylabel,
            legend = ["measured", "independent ideal-fluid reference"],
            dash = ["", "dash"],
            show = false,
        )
        push!(outputs, _save_pdf(joinpath(artifact_dir, "$(id)_$(metric).pdf"), fig))
    end
    orders = _read_csv(joinpath(artifact_dir, "$(id)_orders.csv"))
    isempty(orders) && throw(ArgumentError("oblique order artifact is missing or empty"))
    groups = [
        [r for r in orders if _num(r["Ti"]) == Ti && _num(r["direction"]) == direction] for
        Ti in (0.0, 0.3) for direction in (-1, 1)
    ]
    x = log10.([_num(r["dt"]) for r in first(groups)])
    values = [log10.([_num(r["history_error"]) for r in group]) for group in groups]
    fig = PlotlySupply.plot_scatter(
        x,
        values;
        title = "Oblique RK4 full-vector history error",
        xlabel = "log10(timestep)",
        ylabel = "log10(relative maximum error)",
        legend = ["cold negative", "cold positive", "warm negative", "warm positive"],
        show = false,
    )
    push!(outputs, _save_pdf(joinpath(artifact_dir, "$(id)_orders.pdf"), fig))
    selected = only(
        filter(
            r ->
                r["branch"] == "acoustic" &&
                    r["refinement"] == "baseline" &&
                    _num(r["angle"]) == 60 &&
                    _num(r["direction"]) == 1,
            fits,
        ),
    )
    component = selected["component"]
    history = _time_series_plot(
        artifact_dir,
        "$(id)_$(selected["label"]).csv",
        "$(id)_history.pdf";
        title = "Measured warm acoustic continuation at 60 degrees",
        ycols = (component * "_real", component * "_imag"),
        yaxis_title = "complex Fourier coefficient",
    )
    history === nothing && throw(ArgumentError("selected oblique raw history is unavailable"))
    push!(outputs, history)
    return outputs
end

VALIDATION_PLOT = plot_35_production_oblique_hall
if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
