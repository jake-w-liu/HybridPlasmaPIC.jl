#!/usr/bin/env julia
#
# Phase-1 kinetic instability: the Weibel / current-filamentation instability, driven
# in the full electromagnetic PIC model (EMPIC). Two counter-streaming cold electron
# beams ±u₀ x̂ over an immobile neutralizing ion background carry a velocity-space
# anisotropy ⟨vₓ²⟩ = u₀²+vth² > vth²; the unstable k ∥ ŷ grows the out-of-plane field
# B_z(y) from shot noise. Unstable when A = (u₀/vth)² > 1 (Weibel 1959, Fried 1959).
# Unlike the anisotropy instabilities (cases 29–30, hybrid), the Weibel needs kinetic
# electrons — in the quasineutral hybrid model B=0 is an exact fixed point.

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using HybridPlasmaPIC

function case_31_weibel_instability(artifact_dir::AbstractString)
    id = "31_weibel_instability"
    cfg = (N = (8, 96), L = (4π, 12π), nppc = 50, c = 3.0, dt = 0.05, nsteps = 600, seed = 1)
    u = weibel_growth(; u0 = 0.6, vth = 0.1, cfg...)   # A = 36 (unstable)
    s = weibel_growth(; u0 = 0.0, vth = 0.1, cfg...)   # A = 0 (single Maxwellian, stable)

    unstable_grows = max(0.0, 0.02 - u.wBz_max)
    stable_quiet = max(0.0, s.wBz_max - 0.005)
    threshold_ok = (u.unstable_theory && !s.unstable_theory) ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "$(id).csv")
    rows = (
        ("weibel_unstable_Bz_growth", u.wBz_max, 0.02, "margin", unstable_grows, 0.0),
        ("weibel_no_stream_stays_at_noise", s.wBz_max, 0.0, "margin", stable_quiet, 0.0),
        ("weibel_threshold_matches_theory", threshold_ok, 0.0, "absolute", threshold_ok, 0.0),
    )
    _write_metric_csv(artifact, rows)
    gated = _metric_rows_to_results(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "Weibel current-filamentation threshold A=(u₀/vth)²>1 (Weibel 1959; Fried 1959)",
        rows = rows,
        artifact = artifact,
    )
    skip = _skip_result(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "Weibel B_z growth vs threshold (full EM-PIC)",
        metric = "weibel_Bz_energy",
        artifact = basename(artifact),
        notes = "counter-streaming (A=$(round(Int, u.anisotropy))): peak B_z energy → " *
                "$(round(u.wBz_max, digits = 4)); no streaming (A=0): → " *
                "$(round(s.wBz_max, digits = 5)) (shot noise); separation " *
                "$(round(u.wBz_max / max(s.wBz_max, 1e-12), digits = 0))×.",
    )
    return vcat(gated, [skip])
end

VALIDATION_CASE = ValidationCase(
    id = "31_weibel_instability",
    default = false,
    description = "Weibel current-filamentation instability on the full EM-PIC engine (B_z grows from counter-streaming, quiet without).",
    runner = case_31_weibel_instability,
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
