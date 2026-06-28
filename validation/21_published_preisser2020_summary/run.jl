#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_21_published_preisser2020_summary(artifact_dir::AbstractString)
    id = "21_published_preisser2020_summary"
    reference = published_hybrid_reference(; alpha_fraction = 0.05)
    comparison = compare_to_published_hybrid_reference(reference; alpha_fraction = 0.05, rtol = 0.0)
    details = comparison.comparison.details
    maxerr = comparison.comparison.maxrelerr

    artifact = joinpath(artifact_dir, "21_published_preisser2020_summary.csv")
    rows = map(details) do detail
        key, measured, expected, relerr, ok = detail
        (String(key), measured, expected, relerr, ok)
    end
    _write_csv(artifact, ("metric", "measured", "expected", "relative_error", "pass"), rows)
    return [
        _result(
            id = id,
            category = "hybrid_pic",
            reference_kind = "published_external_summary",
            reference = "Preisser et al. 2020 Zenodo DOI 10.5281/zenodo.3697360, 65deg 5perc Bavg_y summary",
            metric = "max_relative_summary_error",
            measured = maxerr,
            expected = 0.0,
            error_kind = "relative",
            error = maxerr,
            tolerance = 0.0,
            artifact = basename(artifact),
            notes = "Compact bundled scalar summary; full HDF5 source is not vendored.",
        ),
    ]
end


VALIDATION_CASE = ValidationCase(
    id = "21_published_preisser2020_summary",
    default = true,
    description = "Published external hybrid shock scalar summary oracle.",
    runner = case_21_published_preisser2020_summary,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
