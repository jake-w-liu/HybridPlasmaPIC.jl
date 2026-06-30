#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_21_published_preisser2020_summary(artifact_dir::AbstractString)
    id = "21_published_preisser2020_summary"
    ref = published_hybrid_reference(; alpha_fraction = 0.05)

    # (1) Integrity / round-trip: the bundled Preisser scalar summary (real Zenodo data,
    #     with file checksum + provenance) is self-consistent and the comparison harness
    #     reproduces it to ~machine precision.
    rt = compare_to_published_hybrid_reference(ref; alpha_fraction = 0.05, rtol = 1e-9)
    roundtrip_err = rt.comparison.maxrelerr

    # (2) Discrimination: the comparator must REJECT a perturbed copy (+10% on one
    #     observable). This is what makes (1) meaningful — previously this case compared
    #     the reference to itself, which passes trivially regardless of the harness.
    perturbed = merge(ref, (; Bavg_y_max = ref.Bavg_y_max * 1.10))
    rejects_bad =
        compare_to_published_hybrid_reference(perturbed; alpha_fraction = 0.05, rtol = 0.01).pass ?
        1.0 : 0.0

    artifact = joinpath(artifact_dir, "21_published_preisser2020_summary.csv")
    rows = (
        (
            "preisser_summary_roundtrip_rel_error",
            roundtrip_err,
            0.0,
            "relative",
            roundtrip_err,
            1e-8,
        ),
        ("comparator_rejects_perturbed_input", rejects_bad, 0.0, "absolute", rejects_bad, 0.0),
    )
    _write_metric_csv(artifact, rows)
    gated = _metric_rows_to_results(
        id = id,
        category = "hybrid_pic",
        reference_kind = "published_external_summary",
        reference = "Preisser et al. 2020 Zenodo DOI 10.5281/zenodo.3697360, 65deg 5perc Bavg_y summary",
        rows = rows,
        artifact = artifact,
    )
    skip = _skip_result(
        id = id,
        category = "hybrid_pic",
        reference_kind = "published_external_summary",
        reference = "Preisser et al. 2020 (65deg oblique)",
        metric = "our_model_vs_preisser_physics",
        artifact = basename(artifact),
        notes = "BLOCKED (not faked): Preisser is a 65° OBLIQUE shock (Bavg_y); our shock models " *
                "are PERPENDICULAR (Bz), so a live our-sim-vs-Preisser overlay is geometry-incompatible " *
                "(would need an oblique driver). This case verifies the bundled Preisser data integrity " *
                "and that the comparison harness discriminates. For a LIVE external plasma-code physics " *
                "comparison, see case 27 (our whistler vs the NHDS kinetic solver).",
    )
    return vcat(gated, [skip])
end


VALIDATION_CASE = ValidationCase(
    id = "21_published_preisser2020_summary",
    default = true,
    description = "Published external hybrid shock scalar summary oracle.",
    runner = case_21_published_preisser2020_summary,
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
