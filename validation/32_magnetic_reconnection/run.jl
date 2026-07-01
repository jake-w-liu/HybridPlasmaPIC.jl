#!/usr/bin/env julia
#
# Phase-1 kinetic phenomenon: collisionless magnetic reconnection — the tearing mode
# of a Harris current sheet, driven in the 2-D hybrid model. A periodic double-Harris
# equilibrium (B_x(y)=B0[tanh((y-y1)/λ)-tanh((y-y2)/λ)-1], pressure-balanced Harris
# density carried by per-particle weights, sheet current carried by the electron fluid)
# is perturbed with a divergence-free flux-function seed. The reconnected flux is the
# coherent m=1 power of B_y (rfft along x), which isolates the growing tearing island
# from broadband particle noise. Tearing-unstable when kx·λ<1 (Furth-Killeen-Rosenbluth
# 1963); with sheet=false (uniform B_x) there is no free energy and the mode stays flat.

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end
using HybridPlasmaPIC

function case_32_magnetic_reconnection(artifact_dir::AbstractString)
    id = "32_magnetic_reconnection"
    cfg = (
        Nx = 64,
        Ny = 128,
        Lx = 25.6,
        Ly = 25.6,
        nppc = 40,
        dt = 0.015,
        nsteps = 400,
        NB = 4,
        seed = 1,
    )
    u = reconnection_growth(; sheet = true, cfg...)    # tearing-unstable Harris sheet
    s = reconnection_growth(; sheet = false, cfg...)   # uniform B_x (no free energy)

    unstable_grows = max(0.0, 3.0 - u.growth)          # sheet m=1 grows ≥3× (measured ~13×)
    stable_quiet = max(0.0, s.growth - 2.0)            # uniform m=1 stays ≤2× (measured ~1.1×)
    threshold_ok = (u.tearing_theory && !s.tearing_theory) ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "$(id).csv")
    rows = (
        ("reconnection_tearing_mode_grows", u.growth, 3.0, "margin", unstable_grows, 0.0),
        ("reconnection_no_sheet_stays_flat", s.growth, 0.0, "margin", stable_quiet, 0.0),
        (
            "reconnection_tearing_threshold_matches_theory",
            threshold_ok,
            0.0,
            "absolute",
            threshold_ok,
            0.0,
        ),
    )
    _write_metric_csv(artifact, rows)
    gated = _metric_rows_to_results(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "Harris-sheet tearing instability, unstable ⇔ kx·λ<1 (Furth-Killeen-Rosenbluth 1963)",
        rows = rows,
        artifact = artifact,
    )
    skip = _skip_result(
        id = id,
        category = "kinetic_instability",
        reference_kind = "analytic",
        reference = "reconnection m=1 tearing growth vs uniform control",
        metric = "reconnected_flux_m1",
        artifact = basename(artifact),
        notes = "Harris sheet (kx·λ=$(round(2π / 25.6 * 0.5, digits = 3))<1): coherent m=1 " *
                "B_y power → $(round(u.growth, digits = 1))× seed; uniform B_x: → " *
                "$(round(s.growth, digits = 2))× (flat); separation " *
                "$(round(u.growth / max(s.growth, 1e-12), digits = 0))×.",
    )
    return vcat(gated, [skip])
end

VALIDATION_CASE = ValidationCase(
    id = "32_magnetic_reconnection",
    default = false,
    description = "Collisionless magnetic reconnection (Harris-sheet tearing) on the 2-D hybrid engine (m=1 island grows; uniform control flat).",
    runner = case_32_magnetic_reconnection,
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
