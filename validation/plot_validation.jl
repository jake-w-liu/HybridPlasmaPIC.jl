#!/usr/bin/env julia

include(joinpath(@__DIR__, "plot_common.jl"))

const DEFAULT_ARTIFACT_DIR = joinpath(@__DIR__, "artifacts")

function _case_artifact_root(summary_artifact_dir::AbstractString)
    return basename(summary_artifact_dir) == "artifacts" ? dirname(summary_artifact_dir) : summary_artifact_dir
end

function _case_artifact_dir(case_id::AbstractString, summary_artifact_dir::AbstractString)
    return joinpath(_case_artifact_root(summary_artifact_dir), case_id, "artifacts")
end

function _case_plot_dirs()
    dirs = [path for path in readdir(@__DIR__; join = true) if isdir(path) && isfile(joinpath(path, "plot.jl"))]
    sort!(dirs; by = basename)
    return dirs
end

function _load_plots()
    plots = Pair{String,Function}[]
    for dir in _case_plot_dirs()
        script = joinpath(dir, "plot.jl")
        global VALIDATION_PLOT = nothing
        Base.include(@__MODULE__, script)
        plotter = getfield(@__MODULE__, :VALIDATION_PLOT)
        plotter isa Function || error("$script did not define VALIDATION_PLOT")
        push!(plots, basename(dir) => plotter)
    end
    isempty(plots) && error("no validation plot scripts found under $(@__DIR__)")
    return plots
end

function main(args = ARGS)
    options = _parse_plot_args(args; default_artifact_dir = DEFAULT_ARTIFACT_DIR)
    version, path = _require_plotlysupply_180()
    outputs = Any[_summary_plot(options.artifact_dir)]
    for (case_id, plotter) in _load_plots()
        push!(outputs, Base.invokelatest(plotter, _case_artifact_dir(case_id, options.artifact_dir)))
    end
    metadata = _write_plot_metadata(options.artifact_dir, version, path, outputs)
    println("PlotlySupply ", version, " at ", path)
    println("Plot metadata written to ", metadata)
    for output in outputs
        output === nothing || println("PDF written: ", output)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
