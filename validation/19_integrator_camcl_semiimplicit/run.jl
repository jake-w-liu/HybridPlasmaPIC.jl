#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_19_integrator_camcl_semiimplicit(artifact_dir::AbstractString)
    id = "19_integrator_camcl_semiimplicit"
    omega_camcl, omega_ref = _ion_acoustic_frequency(:camcl)
    omega_hybrid, _ = _ion_acoustic_frequency(:hybrid)
    camcl_relerr = abs(omega_camcl - omega_ref) / omega_ref
    camcl_hybrid_relerr = abs(omega_camcl - omega_hybrid) / omega_hybrid

    cn_modulus_error = abs(abs(cn_multiplier(7.0, 0.3)) - 1.0)
    cn = run_whistler(; method = :cn, n = 16, dt = 0.3, nsteps = 25, seed = 4)
    euler = run_whistler(; method = :euler, n = 16, dt = 0.05, nsteps = 3, seed = 4)
    cn_energy_error = abs(cn.energy_ratio - 1.0)
    euler_growth_error = max(0.0, 1.1 - euler.energy_ratio)
    compare = compare_integrators_whistler(;
        n = 16,
        nsteps = 20,
        dt_resolved = 0.005,
        dt_stiff = 0.3,
        seed = 4,
        band_resolved = 1,
    )
    compare_cn_error = abs(compare.cn_ratio_stiff - 1.0)
    compare_euler_growth_error = max(0.0, 10.0 - compare.euler_ratio_stiff)

    artifact = joinpath(artifact_dir, "19_integrator_camcl_semiimplicit.csv")
    rows = (
        ("camcl_ion_acoustic_relative_frequency_error", omega_camcl, omega_ref, "relative", camcl_relerr, 0.04),
        ("camcl_vs_hybrid_frequency_relative_difference", omega_camcl, omega_hybrid, "relative", camcl_hybrid_relerr, 0.05),
        ("cn_multiplier_modulus_error", abs(cn_multiplier(7.0, 0.3)), 1.0, "absolute", cn_modulus_error, 1e-14),
        ("cn_whistler_energy_ratio_error", cn.energy_ratio, 1.0, "absolute", cn_energy_error, 1e-12),
        ("euler_stiff_growth_margin_error", euler.energy_ratio, 1.1, "margin", euler_growth_error, 0.0),
        ("compare_integrators_resolved_agreement_error", compare.agree_resolved, 0.0, "absolute", compare.agree_resolved, 1e-3),
        ("compare_integrators_cn_stiff_energy_error", compare.cn_ratio_stiff, 1.0, "absolute", compare_cn_error, 1e-12),
        ("compare_integrators_euler_stiff_growth_margin_error", compare.euler_ratio_stiff, 10.0, "margin", compare_euler_growth_error, 0.0),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "integrators",
        reference_kind = "analytic",
        reference = "ion-acoustic dispersion and linear whistler amplification factors",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "19_integrator_camcl_semiimplicit",
    default = false,
    description = "CAM-CL ion-acoustic and semi-implicit whistler integrator checks.",
    runner = case_19_integrator_camcl_semiimplicit,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
