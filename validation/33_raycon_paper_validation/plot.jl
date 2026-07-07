#!/usr/bin/env julia
#
# Raycon paper-validation overlays. The PDFs compare computed trajectories and
# conversion markers against figure/text targets recorded by the runner.

if !isdefined(@__MODULE__, :_save_pdf)
    include(joinpath(@__DIR__, "..", "plot_common.jl"))
end

function _group_rows(rows, key)
    groups = Dict{String,Vector{Dict{String,String}}}()
    order = String[]
    for row in rows
        name = row[key]
        if !haskey(groups, name)
            groups[name] = Dict{String,String}[]
            push!(order, name)
        end
        push!(groups[name], row)
    end
    return order, groups
end

function _plot_phase_space_overlay(artifact_dir::AbstractString)
    rows = _read_csv(joinpath(artifact_dir, "33_raycon_paper_validation_phase_space.csv"))
    isempty(rows) && return nothing

    order, groups = _group_rows(rows, "series")
    xseries = Vector{Float64}[]
    yseries = Vector{Float64}[]
    modes = String[]
    legends = String[]
    marker_symbols = String[]
    marker_sizes = Int[]
    linewidths = Float64[]
    for name in order
        grp = sort(groups[name]; by = row -> _num(row["point_index"]))
        kind = grp[1]["kind"]
        push!(xseries, [_num(row["R_m"]) for row in grp])
        push!(yseries, [_num(row["kR_inv_m"]) for row in grp])
        push!(modes, kind == "computed_ray" ? "lines" : "markers")
        push!(legends, name)
        push!(
            marker_symbols,
            kind == "paper_target" ? "x" :
            kind == "saddle" ? "diamond" :
            kind == "transmitted" ? "square" : "circle",
        )
        push!(marker_sizes, kind == "computed_ray" ? 0 : 9)
        push!(linewidths, kind == "computed_ray" ? 2.0 : 0.0)
    end

    fig = PlotlySupply.plot_scatter(
        xseries,
        yseries;
        title = "Raycon Jaun Fig. 3 phase-space comparison",
        xlabel = "major radius R [m]",
        ylabel = "radial wave vector kR [1/m]",
        xrange = [1.9, 3.75],
        yrange = [-45, 45],
        mode = modes,
        legend = legends,
        marker_symbol = marker_symbols,
        marker_size = marker_sizes,
        linewidth = linewidths,
    )
    return _save_pdf(joinpath(artifact_dir, "33_raycon_paper_validation_phase_space.pdf"), fig)
end

function _plot_rz_trace_overlay(
    artifact_dir::AbstractString,
    csv_name::AbstractString,
    pdf_name::AbstractString;
    title,
    xrange,
    yrange,
)
    rows = _read_csv(joinpath(artifact_dir, csv_name))
    isempty(rows) && return nothing

    order, groups = _group_rows(rows, "series")
    xseries = Vector{Float64}[]
    yseries = Vector{Float64}[]
    modes = String[]
    legends = String[]
    marker_symbols = String[]
    marker_sizes = Int[]
    linewidths = Float64[]
    for name in order
        grp = sort(groups[name]; by = row -> _num(row["point_index"]))
        kind = grp[1]["kind"]
        push!(xseries, [_num(row["R_m"]) for row in grp])
        push!(yseries, [_num(row["Z_m"]) for row in grp])
        push!(modes, kind == "computed_ray" ? "lines" : "markers")
        push!(legends, name)
        push!(
            marker_symbols,
            kind == "paper_target" ? "x" :
            kind == "saddle" ? "diamond" :
            kind == "transmitted" ? "square" : "circle",
        )
        push!(marker_sizes, kind == "computed_ray" ? 0 : 9)
        push!(linewidths, kind == "computed_ray" ? 2.0 : 0.0)
    end

    fig = PlotlySupply.plot_scatter(
        xseries,
        yseries;
        title = title,
        xlabel = "major radius R [m]",
        ylabel = "vertical coordinate Z [m]",
        xrange = xrange,
        yrange = yrange,
        mode = modes,
        legend = legends,
        marker_symbol = marker_symbols,
        marker_size = marker_sizes,
        linewidth = linewidths,
    )
    return _save_pdf(joinpath(artifact_dir, pdf_name), fig)
end

function plot_33_raycon_paper_validation(artifact_dir::AbstractString)
    outputs = Any[
        _plot_phase_space_overlay(artifact_dir),
        _plot_rz_trace_overlay(
            artifact_dir,
            "33_raycon_paper_validation_fig6_7_rz.csv",
            "33_raycon_paper_validation_fig6_7_rz.pdf";
            title = "Raycon Jaun Fig. 6/7 off-mid-plane R-Z comparison",
            xrange = [1.8, 3.8],
            yrange = [-0.65, 0.15],
        ),
        _plot_rz_trace_overlay(
            artifact_dir,
            "33_raycon_paper_validation_fig9_10_dh_rz.csv",
            "33_raycon_paper_validation_fig9_10_dh_rz.pdf";
            title = "Raycon Jaun Fig. 9/10 D-H strong-coupling R-Z comparison",
            xrange = [1.8, 3.8],
            yrange = [-0.9, 0.3],
        ),
    ]
    return filter(!isnothing, outputs)
end

VALIDATION_PLOT = plot_33_raycon_paper_validation

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_plot_main(
            VALIDATION_PLOT,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
