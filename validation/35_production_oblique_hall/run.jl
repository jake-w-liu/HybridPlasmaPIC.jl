#!/usr/bin/env julia
if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using SHA, TOML
if !isdefined(@__MODULE__, :_oblique_hall_history)
    include(joinpath(@__DIR__, "..", "oblique_hall_history.jl"))
end

function _oblique_hall_profiles()
    profiles = NamedTuple[]
    for angle in (30.0, 60.0),
        (Ti, Te) in ((0.0, 0.0), (0.3, 0.2)),
        branch in
        (Ti + Te == 0 ? (:whistler, :ioncyclotron) : (:whistler, :ioncyclotron, :acoustic)),
        direction in (-1, 1),
        refinement in (:baseline, :grid, :time)

        label = "angle$(angle)_Ti$(Ti)_Te$(Te)_$(branch)_sign$(direction)_$(refinement)"
        push!(
            profiles,
            (;
                label,
                angle,
                Ti,
                Te,
                branch,
                direction,
                refinement,
                mode = 1,
                n = refinement == :grid ? 32 : 16,
                dt = refinement == :time ? 0.005 : 0.01,
                duration = 80.0,
                amplitude = 1e-5,
            ),
        )
    end
    return profiles
end

function _oblique_order_profiles()
    return [
        (;
            label = "order_Ti$(Ti)_Te$(Te)_sign$(direction)_dt$(dt)",
            angle = 30.0,
            Ti,
            Te,
            branch = :whistler,
            direction,
            mode = 2,
            n = 16,
            dt,
            duration = 1.0,
            amplitude = 1e-6,
        ) for (Ti, Te) in ((0.0, 0.0), (0.3, 0.2)) for direction in (-1, 1) for
        dt in (0.04, 0.02, 0.01)
    ]
end

function _oblique_source_hashes()
    root = dirname(dirname(pathof(HybridPlasmaPIC)))
    paths = String[]
    for (dir, _, files) in walkdir(joinpath(root, "src")), file in files
        endswith(file, ".jl") && push!(paths, joinpath(dir, file))
    end
    append!(
        paths,
        [
            joinpath(root, "Project.toml"),
            @__FILE__,
            joinpath(@__DIR__, "..", "common.jl"),
            joinpath(@__DIR__, "..", "mode_history.jl"),
            joinpath(@__DIR__, "..", "oblique_hall_history.jl"),
            joinpath(@__DIR__, "..", "..", "test", "oracles", "hybrid_dispersion_oracle.jl"),
            joinpath(@__DIR__, "..", "..", "docs", "OBLIQUE_HALL_VALIDATION_PLAN.md"),
        ],
    )
    return Dict(abspath(p) => bytes2hex(sha256(read(p))) for p in paths)
end

_oblique_profile_dict(p) = Dict(string(k) => (v isa Symbol ? string(v) : v) for (k, v) in pairs(p))

function _oblique_receipt(profiles, orders)
    project = Base.active_project()
    environment = Dict{String,String}()
    if project !== nothing
        environment[project] = bytes2hex(sha256(read(project)))
        manifest = joinpath(dirname(project), "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
        isfile(manifest) || (manifest = joinpath(dirname(project), "Manifest.toml"))
        isfile(manifest) && (environment[manifest] = bytes2hex(sha256(read(manifest))))
    end
    return Dict{String,Any}(
        "status" => "running",
        "julia_version" => string(VERSION),
        "machine" => Sys.MACHINE,
        "julia_threads" => Threads.nthreads(),
        "fftw_threads" => FFTW.get_num_threads(),
        "source_sha256" => _oblique_source_hashes(),
        "environment_sha256" => environment,
        "profiles" => _oblique_profile_dict.(profiles),
        "order_profiles" => _oblique_profile_dict.(orders),
        "normalization" => "2*positive-mode Fourier coefficient/grid cells for each real field",
        "fit_window" => "final three quarters including both endpoints",
        "uncertainty" => "conditional OLS slope errors and separate half-window sensitivity; not independent-noise confidence intervals",
        "completed_profiles" => String[],
        "artifact_sha256" => Dict{String,String}(),
        "initializations" => Dict{String,Any}(),
    )
end

function _write_oblique_receipt(path, receipt, artifacts)
    receipt["artifact_sha256"] =
        Dict(basename(p) => bytes2hex(sha256(read(p))) for p in artifacts if isfile(p))
    open(io -> TOML.print(io, receipt; sorted = true), path, "w")
    return nothing
end

function _oblique_run_profile(p, path)
    return _oblique_hall_history(;
        p.branch,
        p.angle,
        p.Ti,
        p.Te,
        p.direction,
        p.mode,
        p.n,
        p.dt,
        p.duration,
        p.amplitude,
        output = path,
    )
end

function _oblique_metric(
    id,
    label,
    metric,
    measured,
    expected,
    error,
    tolerance,
    artifact;
    kind = "absolute",
)
    return _result(;
        id,
        category = "production_dispersion",
        reference_kind = "independent_linear_equations",
        reference = "Isothermal oblique Hall-MHD continuum eigenmode with tracked signed branch",
        metric = "$(label)_$(metric)",
        measured,
        expected,
        error_kind = kind,
        error,
        tolerance,
        artifact = basename(artifact),
        notes = "Actual production Fourier history; full-vector linearization error includes finite-amplitude effects",
    )
end

function _run_oblique_hall_case(artifact_dir, profiles, orders)
    id = "35_production_oblique_hall"
    mkpath(artifact_dir)
    receipt_path = joinpath(artifact_dir, "$(id)_provenance.toml")
    ispath(receipt_path) &&
        throw(ArgumentError("oblique run receipt already exists; use a fresh artifact directory"))
    receipt = _oblique_receipt(profiles, orders)
    artifacts = String[]
    _write_oblique_receipt(receipt_path, receipt, artifacts)
    results = ValidationResult[]
    fits = NamedTuple[]
    order_rows = NamedTuple[]
    try
        for p in profiles
            receipt["active_profile"] = p.label
            _write_oblique_receipt(receipt_path, receipt, artifacts)
            @info "Evolving oblique Hall mode" profile = p.label
            path = joinpath(artifact_dir, "$(id)_$(p.label).csv")
            push!(artifacts, path)
            history = _oblique_run_profile(p, path)
            first_sample = history.steps ÷ 4 + 1
            fit = _fit_mode_history(
                @view(history.samples[history.reference.component, first_sample:end]),
                p.dt,
            )
            errors = _oblique_history_errors(history)
            trace_path = joinpath(artifact_dir, "$(id)_$(p.label)_branches.csv")
            _write_csv(
                trace_path,
                keys(first(history.reference.trace)),
                (values(r) for r in history.reference.trace),
            )
            push!(artifacts, trace_path)
            ref = history.reference
            receipt["initializations"][p.label] = Dict(
                "omega" => ref.omega,
                "gamma" => ref.gamma,
                "component" => ref.component,
                "vector_real" => real.(ref.vector),
                "vector_imag" => imag.(ref.vector),
            )
            push!(
                fits,
                merge(
                    p,
                    (;
                        omega_reference = ref.omega,
                        gamma_reference = ref.gamma,
                        component = _OBLIQUE_COMPONENTS[ref.component],
                        fit_first_sample = first_sample,
                        fit_last_sample = history.steps + 1,
                    ),
                    fit,
                    (; duration = p.duration, fit_duration = fit.duration),
                    errors,
                    (;
                        history.mass_error,
                        history.longitudinal_error,
                        history.density_min,
                        history.density_max,
                        history.floor_count,
                    ),
                ),
            )
            metrics = (
                (
                    "frequency",
                    fit.omega,
                    ref.omega,
                    abs(fit.omega - ref.omega) / abs(ref.omega),
                    1e-3,
                ),
                ("growth_rate", fit.gamma, ref.gamma, abs(fit.gamma - ref.gamma), 1e-4),
                ("phase_residual", fit.phase_rms, 0.0, fit.phase_rms, 1e-3),
                ("log_amplitude_residual", fit.logamp_rms, 0.0, fit.logamp_rms, 1e-3),
                (
                    "frequency_window_difference",
                    fit.omega_window_difference,
                    0.0,
                    fit.omega_window_difference,
                    1e-4,
                ),
                (
                    "growth_window_difference",
                    fit.gamma_window_difference,
                    0.0,
                    fit.gamma_window_difference,
                    1e-4,
                ),
                ("minimum_periods", fit.periods, 2.0, max(0.0, 2 - fit.periods), 0.0),
                ("full_vector_history", errors.history_error, 0.0, errors.history_error, 1e-3),
                ("polarization_leakage", errors.leakage, 0.0, errors.leakage, 1e-3),
                ("periodic_mass", history.mass_error, 0.0, history.mass_error, 1e-12),
                (
                    "longitudinal_magnetic_field",
                    history.longitudinal_error,
                    0.0,
                    history.longitudinal_error,
                    1e-12,
                ),
            )
            for (metric, measured, expected, error, tolerance) in metrics
                kind =
                    metric == "frequency" ? "relative" :
                    metric == "minimum_periods" ? "margin" : "absolute"
                push!(
                    results,
                    _oblique_metric(
                        id,
                        p.label,
                        metric,
                        measured,
                        expected,
                        error,
                        tolerance,
                        path;
                        kind,
                    ),
                )
            end
            push!(receipt["completed_profiles"], p.label)
            _write_oblique_receipt(receipt_path, receipt, artifacts)
        end
        fit_path = joinpath(artifact_dir, "$(id)_fits.csv")
        isempty(fits) && throw(ArgumentError("oblique case requires fit profiles"))
        _write_csv(fit_path, keys(first(fits)), (values(f) for f in fits))
        push!(artifacts, fit_path)
        for baseline in filter(f -> f.refinement == :baseline, fits), refinement in (:grid, :time)
            matching = filter(
                f ->
                    f.refinement == refinement &&
                        (f.angle, f.Ti, f.Te, f.branch, f.direction) == (
                            baseline.angle,
                            baseline.Ti,
                            baseline.Te,
                            baseline.branch,
                            baseline.direction,
                        ),
                fits,
            )
            length(matching) == 1 ||
                throw(ArgumentError("oblique refinement matrix is incomplete or duplicated"))
            other = only(matching)
            for (metric, measured, expected, tolerance) in (
                ("refined_frequency", other.omega, baseline.omega, 1e-4),
                ("refined_growth", other.gamma, baseline.gamma, 1e-5),
            )
                push!(
                    results,
                    _oblique_metric(
                        id,
                        other.label,
                        metric,
                        measured,
                        expected,
                        abs(measured - expected),
                        tolerance,
                        fit_path,
                    ),
                )
            end
        end
        for p in orders
            receipt["active_profile"] = p.label
            _write_oblique_receipt(receipt_path, receipt, artifacts)
            path = joinpath(artifact_dir, "$(id)_$(p.label).csv")
            push!(artifacts, path)
            history = _oblique_run_profile(p, path)
            errors = _oblique_history_errors(history)
            push!(order_rows, merge(p, errors))
            push!(receipt["completed_profiles"], p.label)
            _write_oblique_receipt(receipt_path, receipt, artifacts)
        end
        order_path = joinpath(artifact_dir, "$(id)_orders.csv")
        isempty(order_rows) && throw(ArgumentError("oblique case requires timestep-order profiles"))
        _write_csv(order_path, keys(first(order_rows)), (values(r) for r in order_rows))
        push!(artifacts, order_path)
        for (Ti, Te) in ((0.0, 0.0), (0.3, 0.2)), direction in (-1, 1)
            rows = sort(
                filter(p -> (p.Ti, p.Te, p.direction) == (Ti, Te, direction), order_rows);
                by = p -> -p.dt,
            )
            length(rows) == 3 ||
                throw(ArgumentError("oblique order matrix is incomplete or duplicated"))
            for pair = 1:2
                order = log2(rows[pair].history_error / rows[pair+1].history_error)
                push!(
                    results,
                    _oblique_metric(
                        id,
                        "Ti$(Ti)_Te$(Te)_sign$(direction)",
                        "time_order_$(pair)",
                        order,
                        3.5,
                        isfinite(order) ? max(0.0, 3.5 - order) : Inf,
                        0.0,
                        order_path;
                        kind = "margin",
                    ),
                )
            end
        end
        _oblique_source_hashes() == receipt["source_sha256"] ||
            error("oblique validation source changed during run")
        receipt["status"] = all(r -> r.status == "pass", results) ? "pass" : "fail"
        receipt["active_profile"] = "complete"
        receipt["metric_count"] = length(results)
        _write_oblique_receipt(receipt_path, receipt, artifacts)
    catch err
        # Preserve partial raw data and the cause; this boundary never hides failure.
        receipt["status"] = "error"
        receipt["failure_type"] = string(typeof(err))
        receipt["failure_cause"] = sprint(showerror, err)
        _write_oblique_receipt(receipt_path, receipt, artifacts)
        rethrow()
    end
    return results
end

case_35_production_oblique_hall(artifact_dir::AbstractString) =
    _run_oblique_hall_case(artifact_dir, _oblique_hall_profiles(), _oblique_order_profiles())

VALIDATION_CASE = ValidationCase(
    id = "35_production_oblique_hall",
    default = false,
    description = "Actual cold/warm oblique Hall histories, tracked signed branches and independent refinement.",
    runner = case_35_production_oblique_hall,
)
if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_case_main(
            VALIDATION_CASE,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
