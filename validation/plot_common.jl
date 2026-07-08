#!/usr/bin/env julia

using DelimitedFiles
using TOML
import PlotlySupply

const REQUIRED_PLOTLYSUPPLY_VERSION = v"1.8.0"

function _parse_plot_args(args; default_artifact_dir::AbstractString, allow_cases::Bool = false)
    artifact_dir = default_artifact_dir
    selected = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--artifact-dir"
            i == length(args) && throw(ArgumentError("--artifact-dir requires a directory"))
            artifact_dir = args[i+1]
            i += 1
        elseif allow_cases && arg == "--case"
            i == length(args) && throw(ArgumentError("--case requires a case id"))
            push!(selected, args[i+1])
            i += 1
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
        i += 1
    end
    return (; artifact_dir = abspath(artifact_dir), selected)
end

function _plotlysupply_version()
    pkgdir = dirname(dirname(pathof(PlotlySupply)))
    project = TOML.parsefile(joinpath(pkgdir, "Project.toml"))
    return VersionNumber(project["version"]), pathof(PlotlySupply)
end

function _require_plotlysupply_180()
    version, path = _plotlysupply_version()
    version == REQUIRED_PLOTLYSUPPLY_VERSION || error(
        "global PlotlySupply must resolve to $REQUIRED_PLOTLYSUPPLY_VERSION; got $version at $path",
    )
    return version, path
end

function _read_csv(path::AbstractString)
    isfile(path) || return Dict{String,String}[]
    table = readdlm(path, ',', String; quotes = true)
    ndims(table) == 1 && (table = reshape(table, 1, :))
    size(table, 1) >= 1 || return Dict{String,String}[]
    header = vec(table[1, :])
    rows = Dict{String,String}[]
    for i = 2:size(table, 1)
        row = Dict{String,String}()
        for j in eachindex(header)
            row[header[j]] = table[i, j]
        end
        push!(rows, row)
    end
    return rows
end

_num(s) = tryparse(Float64, s) === nothing ? NaN : parse(Float64, s)

const _CASE_LABELS = Dict(
    "01_analytic_poisson_2d" => "01_",
    "02_analytic_hall_mhd_continuity" => "02_",
    "03_analytic_spectral_operators" => "03_",
    "04_particle_push_collision_boundaries" => "04_",
    "05_particle_coupling_diagnostics" => "05_",
    "06_boundary_loading_kdv_smoothing" => "06_",
    "07_closure_budget_filter_diagnostics" => "07_",
    "08_metrics_loadbalance_validation" => "08_",
    "09_normalization_migration_io" => "09_",
    "10_io_metadata_archive_validation" => "10_",
    "11_api_contract_regression_validation" => "11_",
    "12_model_runtime_smoke_validation" => "12_",
    "13_threaded_backend_api_validation" => "13_",
    "14_parallel_backend_extension_validation" => "14_",
    "15_mpi_single_rank_validation" => "15_",
    "16_pic1d_hdf5_extension_validation" => "16_",
    "17_distributed_fft_roundtrip" => "17_",
    "18_hybrid_ion_acoustic_dispersion" => "18_",
    "19_integrator_camcl_semiimplicit" => "19_",
    "20_empic_transverse_dispersion_2d" => "20_",
    "21_published_preisser2020_summary" => "21_",
    "22_shock_spacecraft_diagnostics" => "22_",
    "23_shock_multidim_ramp_validation" => "23_",
    "24_shock_driver_sweep_validation" => "24_",
    "25_leroy1982_perp_shock" => "25_",
    "26_hellinger2002_resolution" => "26_",
    "27_nhds_dispersion_comparison" => "27_",
    "28_hybridvpic_perp_shock" => "28_",
    "29_firehose_instability" => "29_",
    "30_ion_cyclotron_instability" => "30_",
    "31_weibel_instability" => "31_",
    "32_magnetic_reconnection" => "32_",
    "33_raycon_paper_validation" => "33_",
)

function _case_label(id::AbstractString)
    return get(_CASE_LABELS, id, replace(id, "_" => "<br>"))
end

function _metric_label(metric::AbstractString, index::Integer)
    return "m" * lpad(string(index), 2, "0")
end

function _plot_csv_escape(x)
    s = string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function _write_label_key(path::AbstractString, labels, full_names)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "plot_label,full_name")
        for i in eachindex(labels)
            println(io, _plot_csv_escape(labels[i]), ",", _plot_csv_escape(full_names[i]))
        end
    end
    return path
end

function _paired_dash(n::Integer)
    return [iseven(i) ? "dash" : "" for i = 1:n]
end

function _save_pdf(path::AbstractString, fig)
    mkpath(dirname(path))
    PlotlySupply.savefig(path, fig)
    return path
end

function _plot_output_paths(output)
    output === nothing && return String[]
    if output isa AbstractVector
        paths = String[]
        for item in output
            append!(paths, _plot_output_paths(item))
        end
        return paths
    end
    return [String(output)]
end

function _summary_plot(artifact_dir::AbstractString, selected::Vector{String} = String[])
    rows = _read_csv(joinpath(artifact_dir, "validation_summary.csv"))
    selected_ids = Set(selected)
    comparable = [
        row for row in rows if row["status"] != "skip" &&
        (isempty(selected_ids) || row["id"] in selected_ids) &&
        isfinite(_num(row["error"])) &&
        isfinite(_num(row["tolerance"]))
    ]
    isempty(comparable) && return nothing
    by_case = Dict{String,Float64}()
    for row in comparable
        tol = _num(row["tolerance"])
        err = _num(row["error"])
        ratio = tol == 0 ? (err == 0 ? 0.0 : Inf) : err / tol
        by_case[row["id"]] = max(get(by_case, row["id"], 0.0), ratio)
    end
    ids = sort(collect(keys(by_case)))
    labels = [_case_label(id) for id in ids]
    ratios = [by_case[id] for id in ids]
    _write_label_key(joinpath(artifact_dir, "validation_summary_plot_labels.csv"), labels, ids)
    fig = PlotlySupply.plot_bar(
        labels,
        ratios;
        title = "HybridPlasmaPIC validation max error budget by case",
        xlabel = "case",
        ylabel = "max error / tolerance",
        show = false,
    )
    return _save_pdf(joinpath(artifact_dir, "validation_summary.pdf"), fig)
end

function _metric_plot(
    artifact_dir::AbstractString,
    csv_name::AbstractString,
    pdf_name::AbstractString;
    title,
)
    all_rows = _read_csv(joinpath(artifact_dir, csv_name))
    rows = [
        row for row in all_rows if haskey(row, "metric") &&
        haskey(row, "error") &&
        haskey(row, "tolerance") &&
        haskey(row, "status") &&
        row["status"] != "skip"
    ]
    labels = String[]
    values = Float64[]
    if isempty(rows)
        skipped = [
            row for row in all_rows if
            haskey(row, "metric") && haskey(row, "status") && row["status"] == "skip"
        ]
        isempty(skipped) && return nothing
        labels = [row["metric"] * " (skip)" for row in skipped]
        values = fill(0.0, length(labels))
    else
        labels = [row["metric"] for row in rows]
        for row in rows
            err = _num(row["error"])
            tol = _num(row["tolerance"])
            push!(values, tol == 0 ? err : err / tol)
        end
    end
    full_labels = copy(labels)
    labels = [_metric_label(full_labels[i], i) for i in eachindex(full_labels)]
    _write_label_key(
        joinpath(artifact_dir, splitext(pdf_name)[1] * "_plot_labels.csv"),
        labels,
        full_labels,
    )
    fig = PlotlySupply.plot_bar(
        labels,
        values;
        title = title,
        xlabel = "metric",
        ylabel = "error / tolerance",
        show = false,
    )
    return _save_pdf(joinpath(artifact_dir, pdf_name), fig)
end

function _line_compare_plot(
    artifact_dir::AbstractString,
    csv_name::AbstractString,
    pdf_name::AbstractString;
    title,
    xcol,
    measured_expected_pairs,
    yaxis_title,
)
    rows = _read_csv(joinpath(artifact_dir, csv_name))
    isempty(rows) && return nothing
    x = [_num(row[xcol]) for row in rows]
    series = Vector{Float64}[]
    legends = String[]
    for pair in measured_expected_pairs
        measured_col, expected_col, label = pair
        push!(series, [_num(row[measured_col]) for row in rows])
        push!(legends, "measured " * label)
        push!(series, [_num(row[expected_col]) for row in rows])
        push!(legends, "expected " * label)
    end
    fig = PlotlySupply.plot_scatter(
        x,
        series;
        title = title,
        xlabel = xcol,
        ylabel = yaxis_title,
        legend = legends,
        dash = _paired_dash(length(series)),
        show = false,
    )
    return _save_pdf(joinpath(artifact_dir, pdf_name), fig)
end

function _time_series_plot(
    artifact_dir::AbstractString,
    csv_name::AbstractString,
    pdf_name::AbstractString;
    title,
    ycols,
    yaxis_title,
)
    rows = _read_csv(joinpath(artifact_dir, csv_name))
    isempty(rows) && return nothing
    x = [_num(row["time"]) for row in rows]
    series = [[_num(row[col]) for row in rows] for col in ycols]
    fig = PlotlySupply.plot_scatter(
        x,
        series;
        title = title,
        xlabel = "time",
        ylabel = yaxis_title,
        legend = collect(ycols),
        dash = _paired_dash(length(series)),
        show = false,
    )
    return _save_pdf(joinpath(artifact_dir, pdf_name), fig)
end

function _write_plot_metadata(
    artifact_dir::AbstractString,
    version::VersionNumber,
    path::AbstractString,
    outputs,
)
    output_paths = String[]
    for output in outputs
        append!(output_paths, _plot_output_paths(output))
    end
    metadata_path = joinpath(artifact_dir, "plot_metadata.csv")
    mkpath(artifact_dir)
    open(metadata_path, "w") do io
        println(io, "key,value")
        println(io, "plotlysupply_version,", version)
        println(io, "plotlysupply_path,\"", replace(path, "\"" => "\"\""), "\"")
        println(io, "template,", PlotlySupply.get_default_template())
        println(io, "pdf_count,", length(output_paths))
    end
    return metadata_path
end

function _run_single_plot_main(plotter::Function, args; default_artifact_dir::AbstractString)
    options = _parse_plot_args(args; default_artifact_dir)
    version, path = _require_plotlysupply_180()
    output = plotter(options.artifact_dir)
    println("PlotlySupply ", version, " at ", path)
    for path in _plot_output_paths(output)
        println("PDF written: ", path)
    end
    return 0
end
