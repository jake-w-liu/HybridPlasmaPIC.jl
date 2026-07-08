#!/usr/bin/env julia
#
# Phase-1 kinetic instability: the parallel FIREHOSE, driven on the hybrid engine.
# An anisotropic bi-Maxwellian ion distribution (T_∥ > T_⊥) in a uniform B₀ = x̂ is
# firehose-unstable when β_∥ − β_⊥ > 2 (ions only; the fluid electrons are isotropic),
# i.e. vth_∥² − vth_⊥² > 1 in n=1, B₀=1, μ₀=1 units (Gary 1993, kinetic plasma theory).
# Above threshold the transverse field δB_⊥ grows; below it stays at the noise floor.

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using HybridPlasmaPIC

function case_29_firehose_instability(artifact_dir::AbstractString)
    id = "29_firehose_instability"
    cfg = (N = 128, Lx = 25.6, nppc = 200, nsteps = 1600, dt = 0.02, seed = 1)
    u = firehose_growth(; vth_par = 1.5, vth_perp = 0.3, cfg...)   # β_∥−β_⊥ = 4.32 (unstable)
    s = firehose_growth(; vth_par = 1.0, vth_perp = 0.8, cfg...)   # β_∥−β_⊥ = 0.72 (stable)

    unstable_grows = max(0.0, 0.10 - u.ratio_max)                  # δB_⊥ reaches ≥10% of B₀ energy
    stable_quiet = max(0.0, s.ratio_max - 0.02)                    # sub-threshold stays ≤2% (noise)
    threshold_ok = (u.unstable_theory && !s.unstable_theory) ? 0.0 : 1.0
    separation = u.ratio_max / max(s.ratio_max, eps(Float64))
    separation_error = max(0.0, 100.0 - separation)

    artifact = joinpath(artifact_dir, "$(id).csv")
    rows = (
        ("firehose_unstable_transverse_growth", u.ratio_max, 0.10, "margin", unstable_grows, 0.0),
        ("firehose_sub_threshold_stays_at_noise", s.ratio_max, 0.0, "margin", stable_quiet, 0.0),
        ("firehose_threshold_matches_theory", threshold_ok, 0.0, "absolute", threshold_ok, 0.0),
        (
            "firehose_unstable_vs_stable_energy_separation",
            separation,
            100.0,
            "margin",
            separation_error,
            0.0,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "Parallel firehose threshold β_∥−β_⊥ > 2 (Gary 1993, Theory of Space Plasma Microinstabilities)",
        rows = rows,
        artifact = artifact,
        notes = "The unstable case must grow well above noise, the sub-threshold control must stay quiet, and the measured energy separation must exceed 100x.",
    )
end

VALIDATION_CASE = ValidationCase(
    id = "29_firehose_instability",
    default = false,
    description = "Parallel firehose instability on the hybrid engine (grows above β_∥−β_⊥>2, quiet below).",
    runner = case_29_firehose_instability,
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
