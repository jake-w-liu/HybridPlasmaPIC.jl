#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))

function _case_artifact_root(summary_artifact_dir::AbstractString)
    return basename(summary_artifact_dir) == "artifacts" ? dirname(summary_artifact_dir) : summary_artifact_dir
end

function _case_artifact_dir(case::ValidationCase, summary_artifact_dir::AbstractString)
    return joinpath(_case_artifact_root(summary_artifact_dir), case.id, "artifacts")
end

function _case_script_dirs()
    dirs = [path for path in readdir(@__DIR__; join = true) if isdir(path) && isfile(joinpath(path, "run.jl"))]
    sort!(dirs; by = basename)
    return dirs
end

function _load_cases()
    cases = ValidationCase[]
    for dir in _case_script_dirs()
        script = joinpath(dir, "run.jl")
        global VALIDATION_CASE = nothing
        Base.include(@__MODULE__, script)
        case = getfield(@__MODULE__, :VALIDATION_CASE)
        case isa ValidationCase || error("$script did not define VALIDATION_CASE")
        case.id == basename(dir) || error("$script defines case id $(case.id), expected $(basename(dir))")
        push!(cases, case)
    end
    isempty(cases) && error("no validation cases found under $(@__DIR__)")
    return cases
end

function _parse_args(args)
    selected = String[]
    artifact_dir = DEFAULT_ARTIFACT_DIR
    all_cases = false
    quick = false
    list = false
    plots = true
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--all"
            all_cases = true
        elseif arg == "--quick"
            quick = true
        elseif arg == "--list"
            list = true
        elseif arg == "--no-plots"
            plots = false
        elseif arg == "--plots"
            plots = true
        elseif arg == "--case"
            i == length(args) && throw(ArgumentError("--case requires a case id"))
            push!(selected, args[i + 1])
            i += 1
        elseif arg == "--artifact-dir"
            i == length(args) && throw(ArgumentError("--artifact-dir requires a directory"))
            artifact_dir = abspath(args[i + 1])
            i += 1
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
        i += 1
    end
    quick && all_cases && throw(ArgumentError("choose only one of --quick and --all"))
    return (; selected, artifact_dir = abspath(artifact_dir), all_cases, quick, list, plots)
end

function _selected_cases(cases::Vector{ValidationCase}, options)
    by_id = Dict(case.id => case for case in cases)
    if !isempty(options.selected)
        repeated = String[]
        seen = Set{String}()
        for id in options.selected
            id in seen && push!(repeated, id)
            push!(seen, id)
        end
        isempty(repeated) || throw(ArgumentError("duplicate validation case(s): $(join(unique(repeated), ", "))"))
        unknown = setdiff(options.selected, keys(by_id))
        isempty(unknown) || throw(ArgumentError("unknown validation case(s): $(join(unknown, ", "))"))
        return [by_id[id] for id in options.selected]
    end
    options.all_cases && return cases
    return [case for case in cases if case.default]
end

function _print_cases(cases::Vector{ValidationCase})
    width = maximum(length(case.id) for case in cases) + 2
    for case in cases
        marker = case.default ? "quick" : "all"
        println(rpad(case.id, width), marker, "  ", case.description)
    end
    return nothing
end

function _run_plotter(artifact_dir::AbstractString, case_ids::Vector{String})
    script = joinpath(@__DIR__, "plot_validation.jl")
    args = ["--project=@v#.#", script, "--artifact-dir", artifact_dir]
    for case_id in case_ids
        push!(args, "--case", case_id)
    end
    cmd = Cmd(vcat(Base.julia_cmd().exec, args))
    println("Plotting with global PlotlySupply: ", cmd)
    run(cmd)
    return nothing
end

function main(args = ARGS)
    options = _parse_args(args)
    all_cases = _load_cases()
    if options.list
        _print_cases(all_cases)
        return 0
    end

    cases = _selected_cases(all_cases, options)
    mkpath(options.artifact_dir)
    _clean_generated_artifacts!(options.artifact_dir)
    for case in all_cases
        _clean_generated_artifacts!(_case_artifact_dir(case, options.artifact_dir))
    end

    results = ValidationResult[]
    for case in cases
        case_artifact_dir = _case_artifact_dir(case, options.artifact_dir)
        mkpath(case_artifact_dir)
        print("running ", case.id, " ... ")
        started = time()
        case_results = Base.invokelatest(case.runner, case_artifact_dir)
        append!(results, case_results)
        elapsed = time() - started
        statuses = unique(r.status for r in case_results)
        println(join(statuses, ","), " (", @sprintf("%.2fs", elapsed), ")")
    end

    summary_path = joinpath(options.artifact_dir, "validation_summary.csv")
    _write_summary_csv(summary_path, results)
    metadata_path = joinpath(options.artifact_dir, "validation_metadata.csv")
    _write_csv(
        metadata_path,
        ("key", "value"),
        (
            ("timestamp_utc", Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ")),
            ("julia_version", string(VERSION)),
            ("artifact_dir", options.artifact_dir),
            ("case_artifact_root", _case_artifact_root(options.artifact_dir)),
            ("case_count", length(cases)),
            ("result_count", length(results)),
        ),
    )

    println("\nSummary written to ", summary_path)
    _print_results(results)

    if options.plots
        _run_plotter(options.artifact_dir, [case.id for case in cases])
    end

    any(r.status == "fail" for r in results) && return 1
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
