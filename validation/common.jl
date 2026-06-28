#!/usr/bin/env julia

using Dates
using FFTW
using HybridPlasmaPIC
using LinearAlgebra
using Printf
using Random
using Serialization

const DEFAULT_ARTIFACT_DIR = joinpath(@__DIR__, "artifacts")

Base.@kwdef struct ValidationCase
    id::String
    default::Bool
    description::String
    runner::Function
end

Base.@kwdef struct ValidationResult
    id::String
    category::String
    reference_kind::String
    reference::String
    metric::String
    measured::Float64
    expected::Float64
    error_kind::String
    error::Float64
    tolerance::Float64
    status::String
    artifact::String
    notes::String
end

_finite(x) = isfinite(x) ? x : NaN

function _result(;
    id,
    category,
    reference_kind,
    reference,
    metric,
    measured,
    expected,
    error_kind,
    error,
    tolerance,
    artifact = "",
    notes = "",
)
    status = isfinite(error) && error <= tolerance ? "pass" : "fail"
    return ValidationResult(
        String(id),
        String(category),
        String(reference_kind),
        String(reference),
        String(metric),
        Float64(measured),
        Float64(expected),
        String(error_kind),
        Float64(error),
        Float64(tolerance),
        status,
        String(artifact),
        String(notes),
    )
end

function _skip_result(; id, category, reference_kind, reference, metric, artifact = "", notes)
    return ValidationResult(
        String(id),
        String(category),
        String(reference_kind),
        String(reference),
        String(metric),
        NaN,
        NaN,
        "skip",
        NaN,
        NaN,
        "skip",
        String(artifact),
        String(notes),
    )
end

function _csv_escape(x)
    s = string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function _write_csv(path::AbstractString, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(_csv_escape.(header), ","))
        for row in rows
            println(io, join(_csv_escape.(row), ","))
        end
    end
    return path
end

function _write_summary_csv(path::AbstractString, results::Vector{ValidationResult})
    header = (
        "id",
        "category",
        "reference_kind",
        "reference",
        "metric",
        "measured",
        "expected",
        "error_kind",
        "error",
        "tolerance",
        "status",
        "artifact",
        "notes",
    )
    rows = map(results) do r
        (
            r.id,
            r.category,
            r.reference_kind,
            r.reference,
            r.metric,
            r.measured,
            r.expected,
            r.error_kind,
            r.error,
            r.tolerance,
            r.status,
            r.artifact,
            r.notes,
        )
    end
    return _write_csv(path, header, rows)
end

function _clean_generated_artifacts!(artifact_dir::AbstractString)
    isdir(artifact_dir) || return nothing
    for name in readdir(artifact_dir)
        path = joinpath(artifact_dir, name)
        isfile(path) || continue
        if endswith(name, ".csv") || endswith(name, ".pdf")
            rm(path; force = true)
        end
    end
    return nothing
end

function _metric_rows_to_results(;
    id,
    category,
    reference_kind,
    reference,
    rows,
    artifact,
    notes = "",
)
    return [
        _result(
            id = id,
            category = category,
            reference_kind = reference_kind,
            reference = reference,
            metric = row[1],
            measured = row[2],
            expected = row[3],
            error_kind = row[4],
            error = row[5],
            tolerance = row[6],
            artifact = basename(artifact),
            notes = notes,
        ) for row in rows
    ]
end

function _write_metric_csv(path::AbstractString, rows)
    outrows = map(rows) do row
        metric, measured, expected, error_kind, error, tolerance = row
        status = isfinite(error) && error <= tolerance ? "pass" : "fail"
        (metric, measured, expected, error_kind, error, tolerance, status)
    end
    return _write_csv(
        path,
        ("metric", "measured", "expected", "error_kind", "error", "tolerance", "status"),
        outrows,
    )
end

function _write_skip_metric_csv(path::AbstractString, metric::AbstractString)
    return _write_csv(
        path,
        ("metric", "measured", "expected", "error_kind", "error", "tolerance", "status"),
        ((metric, NaN, NaN, "skip", NaN, NaN, "skip"),),
    )
end

function _peak_frequency(series::AbstractVector{<:Complex}, dt::Real)
    nt = length(series)
    nt >= 5 || return NaN
    spectrum = fft(series)
    mag = abs.(spectrum)
    best = 2
    bestmag = -Inf
    for j = 2:(nt ÷ 2)
        if mag[j + 1] > bestmag
            bestmag = mag[j + 1]
            best = j
        end
    end
    best <= 1 && return NaN
    a = mag[best]
    b = mag[best + 1]
    c = mag[best + 2]
    denom = a - 2b + c
    offset = denom != 0 ? 0.5 * (a - c) / denom : 0.0
    return 2π * (best + offset) / (nt * dt)
end

function _zero_cross_frequency(series::AbstractVector{<:Real}, dt::Real)
    crossings = Int[]
    for i = 2:length(series)
        if series[i - 1] < 0 && series[i] >= 0
            push!(crossings, i)
        end
    end
    length(crossings) >= 2 || return NaN
    period = (crossings[end] - crossings[1]) / (length(crossings) - 1) * dt
    return 2π / period
end


function _weighted_totals(ps)
    px = sum(ps.weight .* ps.v[1])
    py = sum(ps.weight .* ps.v[2])
    pz = sum(ps.weight .* ps.v[3])
    energy = sum(@. ps.weight * (ps.v[1]^2 + ps.v[2]^2 + ps.v[3]^2))
    return (px, py, pz), energy
end


function _checkpoint_validation_run(seed)
    g = FourierGrid((8,), (2π,))
    np = 20 * g.n[1]
    ps = ParticleSet{1,Float64}(np)
    load_lattice_1d!(ps, 0.0, g.L[1])
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), (0.05, 0.05, 0.05))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.25)), CIC(), np)
    fill!(st.fields.B[1], 1.0)
    xs = [(i - 1) * g.dx[1] for i = 1:g.n[1]]
    st.fields.B[2] .= 0.01 .* cos.(xs)
    init!(st, ps)
    return ps, st
end


function _setup_hybrid_1d(n, l, seed; nppc = 400, vth = (0.0, 0.0, 0.0), te = 0.0, b0 = (0.0, 0.0, 0.0))
    g = FourierGrid((n,), (Float64(l),))
    np = nppc * n
    ps = ParticleSet{1,Float64}(np)
    load_lattice_1d!(ps, 0.0, Float64(l))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), vth)
    st = HybridStepper(g, HybridModel(IsothermalElectrons(Float64(te))), CIC(), np)
    for c = 1:3
        fill!(st.fields.B[c], Float64(b0[c]))
    end
    return g, ps, st
end


function _ion_acoustic_frequency(integrator::Symbol)
    te = 1.0
    mode = 1
    l = 2π
    k = 2π * mode / l
    g, ps, hybrid = _setup_hybrid_1d(48, l, 3; nppc = 400, te = te)
    for p = 1:nparticles(ps)
        ps.v[1][p] += 0.005 * sin(k * ps.x[1][p])
    end

    dt = 0.02
    nt = 700
    series = Float64[]
    if integrator === :hybrid
        init!(hybrid, ps)
        for _ = 1:nt
            step!(hybrid, ps, dt)
            push!(series, real(mode_amplitude(hybrid.fields.n, g, (mode,))))
        end
    elseif integrator === :camcl
        camcl = CAMCLStepper(g, HybridModel(IsothermalElectrons(te)), CIC(), nparticles(ps))
        init_camcl!(camcl, ps)
        for _ = 1:nt
            step_camcl!(camcl, ps, dt; NB = 2)
            push!(series, real(mode_amplitude(camcl.fields.n, g, (mode,))))
        end
    else
        throw(ArgumentError("unknown ion-acoustic integrator: $integrator"))
    end
    return _zero_cross_frequency(series, dt), k * sqrt(te)
end

function _parse_single_artifact_args(args; default_artifact_dir::AbstractString)
    artifact_dir = default_artifact_dir
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--artifact-dir"
            i == length(args) && throw(ArgumentError("--artifact-dir requires a directory"))
            artifact_dir = args[i + 1]
            i += 1
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
        i += 1
    end
    return (; artifact_dir = abspath(artifact_dir))
end

function _print_results(results::Vector{ValidationResult})
    for r in results
        println(
            r.status == "pass" ? "[PASS] " : r.status == "skip" ? "[SKIP] " : "[FAIL] ",
            r.id,
            " / ",
            r.metric,
            " error=",
            r.error,
            " tolerance=",
            r.tolerance,
        )
    end
    return nothing
end

function _run_single_case_main(case::ValidationCase, args; default_artifact_dir::AbstractString)
    options = _parse_single_artifact_args(args; default_artifact_dir)
    mkpath(options.artifact_dir)
    _clean_generated_artifacts!(options.artifact_dir)

    print("running ", case.id, " ... ")
    started = time()
    results = case.runner(options.artifact_dir)
    elapsed = time() - started
    statuses = unique(r.status for r in results)
    println(join(statuses, ","), " (", @sprintf("%.2fs", elapsed), ")")

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
            ("case_id", case.id),
            ("result_count", length(results)),
        ),
    )

    println("Summary written to ", summary_path)
    _print_results(results)
    any(r.status == "fail" for r in results) && return 1
    return 0
end
