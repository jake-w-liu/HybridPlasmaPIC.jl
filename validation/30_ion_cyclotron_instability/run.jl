#!/usr/bin/env julia
#
# Phase-1 kinetic instability: the electromagnetic ion-cyclotron (EMIC) anisotropy
# instability — the T_⊥ > T_∥ counterpart of the firehose (case 29), driven on the
# hybrid engine. Anisotropic bi-Maxwellian ions in a uniform B₀ = x̂. With
# A = T_⊥/T_∥ − 1 and β_∥ = 2 T_∥, the plasma is EMIC-unstable when A > 0.43/β_∥^0.43
# (Gary 1993). Above threshold the transverse δB_⊥ grows; below it stays at noise.

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using HybridPlasmaPIC

function case_30_ion_cyclotron_instability(artifact_dir::AbstractString)
    id = "30_ion_cyclotron_instability"
    cfg = (N = 128, Lx = 25.6, nppc = 200, nsteps = 1600, dt = 0.02, seed = 1)
    u = ion_cyclotron_growth(; vth_par = 0.4, vth_perp = 1.3, cfg...)   # A=9.6 (unstable)
    s = ion_cyclotron_growth(; vth_par = 0.8, vth_perp = 0.9, cfg...)   # A=0.27 (sub-threshold)

    unstable_grows = max(0.0, 0.05 - u.ratio_max)
    stable_quiet = max(0.0, s.ratio_max - 0.02)
    threshold_ok = (u.unstable_theory && !s.unstable_theory) ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "$(id).csv")
    rows = (
        ("emic_unstable_transverse_growth", u.ratio_max, 0.05, "margin", unstable_grows, 0.0),
        ("emic_sub_threshold_stays_at_noise", s.ratio_max, 0.0, "margin", stable_quiet, 0.0),
        ("emic_threshold_matches_theory", threshold_ok, 0.0, "absolute", threshold_ok, 0.0),
    )
    _write_metric_csv(artifact, rows)
    gated = _metric_rows_to_results(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "EMIC threshold A > 0.43/β_∥^0.43 (Gary 1993, Theory of Space Plasma Microinstabilities)",
        rows = rows,
        artifact = artifact,
    )
    skip = _skip_result(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "EMIC growth vs threshold",
        metric = "emic_transverse_energy",
        artifact = basename(artifact),
        notes = "unstable (A=$(round(u.T_anisotropy-1,digits=2)),β∥=$(round(u.beta_par,digits=2))): " *
                "δB_⊥/B₀ energy → $(round(u.ratio_max,digits=3)); sub-threshold " *
                "(A=$(round(s.T_anisotropy-1,digits=2))): → $(round(s.ratio_max,digits=4)) (noise); " *
                "separation $(round(u.ratio_max/max(s.ratio_max,1e-12),digits=0))×.",
    )
    return vcat(gated, [skip])
end

VALIDATION_CASE = ValidationCase(
    id = "30_ion_cyclotron_instability",
    default = false,
    description = "EMIC (T_⊥>T_∥) ion-cyclotron instability on the hybrid engine (grows above Gary threshold, quiet below).",
    runner = case_30_ion_cyclotron_instability,
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
