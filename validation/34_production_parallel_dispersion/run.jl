#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using SHA, TOML
if !isdefined(@__MODULE__, :_parallel_wave_history)
    include(joinpath(@__DIR__, "..", "parallel_wave_history.jl"))
end

function _parallel_wave_profiles()
    profiles = NamedTuple[]
    for solver in (:hall, :hybrid),
        eta in (solver == :hall ? (0.0, 0.03) : (0.0,)),
        branch in (:whistler, :ioncyclotron),
        mode in (1, 2)

        refinements =
            solver == :hall ? (:baseline, :grid, :time) : (:baseline, :grid, :time, :particles)
        for refinement in refinements
            label = "$(solver)_eta$(eta)_$(branch)_m$(mode)_$(refinement)"
            n = refinement == :grid ? 32 : 16
            dt = refinement == :time ? 0.005 : 0.01
            nppc = refinement == :particles ? 32 : 16
            push!(
                profiles,
                (;
                    label,
                    solver,
                    eta,
                    branch,
                    mode,
                    refinement,
                    n,
                    dt,
                    nppc,
                    duration = 40.0,
                    amplitude = 1e-4,
                    subcycles = 8,
                ),
            )
        end
    end
    return profiles
end

function _parallel_wave_receipt(profiles)
    source_root = dirname(dirname(pathof(HybridPlasmaPIC)))
    files = Dict{String,String}()
    for (dir, _, names) in walkdir(joinpath(source_root, "src")), name in names
        endswith(name, ".jl") || continue
        path = joinpath(dir, name)
        files[path] = bytes2hex(sha256(read(path)))
    end
    for path in (
        joinpath(source_root, "Project.toml"),
        @__FILE__,
        joinpath(@__DIR__, "..", "common.jl"),
        joinpath(@__DIR__, "..", "mode_history.jl"),
        joinpath(@__DIR__, "..", "parallel_wave_history.jl"),
    )
        files[abspath(path)] = bytes2hex(sha256(read(path)))
    end
    environment = Dict{String,String}()
    project = Base.active_project()
    if project !== nothing
        environment[project] = bytes2hex(sha256(read(project)))
        # Prefer the version-specific manifest if that is the active environment's layout.
        manifest = joinpath(dirname(project), "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
        isfile(manifest) || (manifest = joinpath(dirname(project), "Manifest.toml"))
        isfile(manifest) && (environment[manifest] = bytes2hex(sha256(read(manifest))))
    end
    return Dict(
        "status" => "running",
        "julia_version" => string(VERSION),
        "machine" => Sys.MACHINE,
        "julia_threads" => Threads.nthreads(),
        "fftw_threads" => FFTW.get_num_threads(),
        "package_root" => source_root,
        "source_sha256" => files,
        "environment_sha256" => environment,
        "fit_window" => "last three quarters, including endpoints",
        "uncertainty" => "OLS slope standard errors are conditional; window differences are separate sensitivity diagnostics",
        "profiles" => [
            Dict(string(k) => (v isa Symbol ? string(v) : v) for (k, v) in pairs(p)) for
            p in profiles
        ],
    )
end

function _write_parallel_wave_history(path, history)
    return _write_csv(
        path,
        ("time", "mode_real", "mode_imag", "mode_abs"),
        ((t, real(z), imag(z), abs(z)) for (t, z) in zip(history.times, history.samples)),
    )
end

function case_34_production_parallel_dispersion(artifact_dir::AbstractString)
    id = "34_production_parallel_dispersion"
    mkpath(artifact_dir)
    profiles = _parallel_wave_profiles()
    receipt = _parallel_wave_receipt(profiles)
    receipt_path = joinpath(artifact_dir, "$(id)_provenance.toml")
    open(io -> TOML.print(io, receipt; sorted = true), receipt_path, "w")
    results = ValidationResult[]
    fits = NamedTuple[]
    artifacts = String[]
    for p in profiles
        @info "Evolving parallel mode" profile = p.label
        history = _parallel_wave_history(;
            p.solver,
            p.branch,
            p.mode,
            p.n,
            p.dt,
            p.nppc,
            p.eta,
            p.duration,
            p.amplitude,
            p.subcycles,
        )
        filename = "$(id)_$(p.label).csv"
        path = joinpath(artifact_dir, filename)
        _write_parallel_wave_history(path, history)
        push!(artifacts, path)
        first_sample = history.steps ÷ 4 + 1
        fit = _fit_mode_history(@view(history.samples[first_sample:end]), p.dt)
        omega_reference = -imag(history.reference.rate)
        gamma_reference = real(history.reference.rate)
        push!(
            fits,
            merge(
                p,
                (;
                    omega_reference,
                    gamma_reference,
                    coupling = history.discrete_reference.coupling,
                ),
                fit,
                (; duration = p.duration, fit_duration = fit.duration),
            ),
        )
        rows = (
            (
                "frequency",
                fit.omega,
                omega_reference,
                "relative",
                abs(fit.omega - omega_reference) / abs(omega_reference),
                0.03,
            ),
            (
                "growth_rate",
                fit.gamma,
                gamma_reference,
                "absolute",
                abs(fit.gamma - gamma_reference),
                0.01,
            ),
            ("phase_residual", fit.phase_rms, 0.0, "absolute", fit.phase_rms, 0.03),
            ("log_amplitude_residual", fit.logamp_rms, 0.0, "absolute", fit.logamp_rms, 0.03),
            (
                "frequency_window_difference",
                fit.omega_window_difference,
                0.0,
                "absolute",
                fit.omega_window_difference,
                0.01,
            ),
            (
                "growth_window_difference",
                fit.gamma_window_difference,
                0.0,
                "absolute",
                fit.gamma_window_difference,
                0.01,
            ),
            ("minimum_periods", fit.periods, 2.0, "margin", max(0.0, 2 - fit.periods), 0.0),
            (
                "uniform_density",
                history.density_error,
                0.0,
                "absolute",
                history.density_error,
                1e-6,
            ),
        )
        for (metric, measured, expected, error_kind, error, tolerance) in rows
            push!(
                results,
                _result(;
                    id,
                    category = "production_dispersion",
                    reference_kind = "independent_linear_equations",
                    reference = "Cold parallel continuum Hall/ion equations; single discrete CIC eigenmode initializes hybrid",
                    metric = "$(p.label)_$(metric)",
                    measured,
                    expected,
                    error_kind,
                    error,
                    tolerance,
                    artifact = filename,
                    notes = "Signed omega in exp((gamma-i*omega)t); OLS errors and window sensitivity retained separately in fit CSV",
                ),
            )
        end
    end
    fit_path = joinpath(artifact_dir, "$(id)_fits.csv")
    _write_csv(fit_path, keys(first(fits)), (values(f) for f in fits))
    push!(artifacts, fit_path)
    order_rows = NamedTuple[]
    for solver in (:hall, :hybrid)
        errors = Float64[]
        for dt in (0.04, 0.02, 0.01)
            history = _parallel_wave_history(;
                solver,
                branch = :whistler,
                mode = 2,
                n = 16,
                nppc = 16,
                dt,
                duration = 1.0,
                amplitude = 1e-6,
            )
            filename = "$(id)_order_$(solver)_dt$(dt).csv"
            path = joinpath(artifact_dir, filename)
            _write_parallel_wave_history(path, history)
            push!(artifacts, path)
            rate = history.discrete_reference.rate
            error =
                maximum(
                    abs(history.samples[j] - history.samples[1] * exp(rate * history.times[j])) for
                    j in eachindex(history.times)
                ) / abs(history.samples[1])
            push!(errors, error)
            push!(order_rows, (; solver, dt, error))
        end
        minimum_order = solver == :hall ? 3.5 : 1.5
        for pair = 1:2
            order = log2(errors[pair] / errors[pair+1])
            push!(
                results,
                _result(;
                    id,
                    category = "production_dispersion",
                    reference_kind = "independent_linear_equations",
                    reference = "Complete complex history against fixed continuum/discrete linear eigenmode",
                    metric = "$(solver)_time_order_$(pair)",
                    measured = order,
                    expected = minimum_order,
                    error_kind = "margin",
                    error = isfinite(order) ? max(0.0, minimum_order - order) : Inf,
                    tolerance = 0.0,
                    artifact = "$(id)_orders.csv",
                ),
            )
        end
    end
    order_path = joinpath(artifact_dir, "$(id)_orders.csv")
    _write_csv(order_path, keys(first(order_rows)), (values(r) for r in order_rows))
    push!(artifacts, order_path)
    receipt["status"] = all(r -> r.status == "pass", results) ? "pass" : "fail"
    receipt["artifact_sha256"] = Dict(basename(p) => bytes2hex(sha256(read(p))) for p in artifacts)
    open(io -> TOML.print(io, receipt; sorted = true), receipt_path, "w")
    return results
end

VALIDATION_CASE = ValidationCase(
    id = "34_production_parallel_dispersion",
    default = false,
    description = "Production hybrid/Hall wave frequencies, damping, fit stability and independent refinement.",
    runner = case_34_production_parallel_dispersion,
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
