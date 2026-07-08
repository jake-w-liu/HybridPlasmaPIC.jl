#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function _max_sweep_difference(a, b)
    length(a) == length(b) || return Inf
    err = 0.0
    for i in eachindex(a)
        err = max(err, abs(Float64(a[i].MA) - Float64(b[i].MA)))
        for name in (:n2, :Bz2, :Vs, :X_rh, :frozen_ratio, :reflected_fraction, :M_real, :xf)
            err = max(err, abs(getproperty(a[i], name) - getproperty(b[i], name)))
        end
    end
    return err
end

function _compression_violation(results)
    err = 0.0
    for r in results
        err = max(err, max(0.0, 2.0 - r.n2), max(0.0, r.n2 - 4.0))
    end
    return err
end

function _mass_speed_error(results)
    err = 0.0
    for r in results
        expected = Float64(r.MA) / (r.n2 - 1.0)
        err = max(err, abs(r.Vs - expected) / expected)
    end
    return err
end

function case_24_shock_driver_sweep_validation(artifact_dir::AbstractString)
    id = "24_shock_driver_sweep_validation"

    direct = run_perp_shock(; MA = 2.0, N = 128, Lx = 60.0, nppc = 8, nsteps = 160, seed = 2)
    direct_frozen_error = abs(direct.frozen_ratio - 1.0)
    direct_contract_error =
        all(isfinite, (direct.n2, direct.Bz2, direct.Vs, direct.M_real, direct.xf)) &&
        1.0 < direct.n2 < 4.0 &&
        0.0 <= direct.reflected_fraction <= 1.0 ? 0.0 : 1.0

    established = reproduce_established_shock(; MA = 3.0, N = 256, nsteps = 500)
    established_error = established.pass ? 0.0 : 1.0

    ramp_scan = ramp_width_scan(;
        widths = (1.0, 2.0),
        N = 64,
        Lx = 20.0,
        x_ramp = 10.0,
        nppc = 4,
        nsteps = 0,
        seed = 1,
    )
    ramp_scan_error =
        length(ramp_scan) == 2 &&
        ramp_scan[2].width_measured > ramp_scan[1].width_measured &&
        all(r -> all(isfinite, (r.xf, r.width_measured, r.n2_meas, r.Bz_jump)), ramp_scan) &&
        maximum(abs(r.Bz_jump - 2.0) for r in ramp_scan) < 1e-2 ? 0.0 : 1.0

    box_scan = box_length_scan(;
        Lxs = (30.0, 40.0),
        MA = 2.0,
        N0 = 64,
        Lx0 = 30.0,
        nppc = 4,
        nsteps = 80,
        seed = 2,
    )
    box_scan_error =
        length(box_scan) == 2 &&
        box_scan[1].N == 64 &&
        box_scan[2].N == 85 &&
        all(r -> all(isfinite, (r.n2, r.Bz2, r.Vs)), box_scan) ? 0.0 : 1.0

    direct3 = run_perp_shock3d(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = 0,
        dt = 0.03,
        seed = 1,
    )
    direct3_error =
        all(
            isfinite,
            (direct3.n2, direct3.frozen_ratio, direct3.maxdivB, direct3.mean_xs, direct3.M_real),
        ) ? 0.0 : 1.0

    traces = four_spacecraft_traces(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = 1,
        dt = 0.03,
        seed = 1,
        level = 1.0,
    )
    trace_values = [traces.traces[q][1] for q = 1:4]
    trace_driver_error =
        length(traces.times) == 1 &&
        all(q -> length(traces.traces[q]) == 1, 1:4) &&
        all(isfinite, trace_values) &&
        isnan(traces.speed) &&
        all(isnan, traces.normal) ? 0.0 : 1.0

    sweep = perp_shock_sweep(
        (2.0, 4.0);
        N = 256,
        Lx = 100.0,
        Te = 0.125,
        vthi = 0.35,
        η = 0.02,
        nppc = 16,
        nsteps = 400,
        seed = 1,
    )
    mach = mach_sweep(;
        MAs = [2.0, 4.0],
        N = 256,
        Lx = 100.0,
        Te = 0.125,
        vthi = 0.35,
        η = 0.02,
        nppc = 16,
        nsteps = 400,
        seed = 1,
    )
    driver_equivalence_error = _max_sweep_difference(sweep, mach)
    frozen_error = maximum(abs(r.frozen_ratio - 1.0) for r in sweep)
    compression_error = _compression_violation(sweep)
    speed_error = _mass_speed_error(sweep)
    reflected_error =
        all(
            0.0 <= r.reflected_fraction <= 1.0 && all(isfinite, (r.n2, r.Bz2, r.Vs, r.M_real)) for
            r in sweep
        ) ? 0.0 : 1.0

    cv = convergence_study(;
        MA = 4.0,
        Ns = (192, 256),
        ppcs = (12, 16),
        Lx = 100.0,
        Te = 0.125,
        vthi = 0.35,
        η = 0.02,
        nsteps = 350,
        seed = 1,
        tol = 0.15,
    )
    convergence_error = cv.rel_max
    convergence_flag_error = cv.converged ? 0.0 : 1.0

    prod = production_3d_case(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps_pre = 0,
        nsteps_post = 0,
        dt = 0.03,
        seed = 1,
    )
    production_restart_error = prod.restart_bitmatch ? 0.0 : 1.0
    production_finite_error = all(isfinite, (prod.n2, prod.frozen_ratio, prod.maxdivB)) ? 0.0 : 1.0

    campaign = shock_campaign_3d(;
        MAs = (3.0,),
        seeds = (1,),
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = 0,
        dt = 0.03,
    )
    campaign_error =
        length(campaign) == 1 &&
        campaign[1].MA == 3.0 &&
        isfinite(campaign[1].n2_mean) &&
        isfinite(campaign[1].frozen_mean) ? 0.0 : 1.0

    dims = compare_dims_shock(;
        MA = 3.0,
        oned_kwargs = (; N = 64, Lx = 30.0, nppc = 4, nsteps = 80),
        threed_kwargs = (;
            nx = 8,
            ny = 4,
            nz = 4,
            Lx = 12.0,
            Ly = 4.0,
            Lz = 4.0,
            nppc = 1,
            nsteps = 0,
            dt = 0.03,
        ),
    )
    compare_dims_error =
        isfinite(dims.oned.n2) &&
        isfinite(dims.oned.frozen_ratio) &&
        isfinite(dims.threed.n2) &&
        isfinite(dims.threed.frozen_ratio) ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "24_shock_driver_sweep_validation.csv")
    rows = (
        (
            "run_perp_shock_direct_frozen_abs_error",
            direct.frozen_ratio,
            1.0,
            "absolute",
            direct_frozen_error,
            0.15,
        ),
        (
            "run_perp_shock_direct_contract_error",
            direct_contract_error,
            0.0,
            "absolute",
            direct_contract_error,
            0.0,
        ),
        (
            "reproduce_established_shock_oracle_error",
            established_error,
            0.0,
            "absolute",
            established_error,
            0.0,
        ),
        ("ramp_width_scan_contract_error", ramp_scan_error, 0.0, "absolute", ramp_scan_error, 0.0),
        ("box_length_scan_contract_error", box_scan_error, 0.0, "absolute", box_scan_error, 0.0),
        (
            "run_perp_shock3d_direct_finite_error",
            direct3_error,
            0.0,
            "absolute",
            direct3_error,
            0.0,
        ),
        (
            "four_spacecraft_traces_one_step_contract_error",
            trace_driver_error,
            0.0,
            "absolute",
            trace_driver_error,
            0.0,
        ),
        (
            "perp_vs_mach_sweep_max_abs_error",
            driver_equivalence_error,
            0.0,
            "absolute",
            driver_equivalence_error,
            0.0,
        ),
        (
            "supercritical_frozen_in_max_abs_error",
            maximum(abs(r.frozen_ratio - 1.0) for r in sweep),
            0.0,
            "absolute",
            frozen_error,
            0.08,
        ),
        (
            "supercritical_compression_band_violation",
            compression_error,
            0.0,
            "margin",
            compression_error,
            0.0,
        ),
        (
            "supercritical_mass_speed_relative_error",
            speed_error,
            0.0,
            "relative",
            speed_error,
            1e-12,
        ),
        (
            "supercritical_reflected_fraction_contract_error",
            reflected_error,
            0.0,
            "absolute",
            reflected_error,
            0.0,
        ),
        ("shock_convergence_rel_max", cv.rel_max, 0.0, "absolute", convergence_error, 0.15),
        (
            "shock_convergence_flag_error",
            convergence_flag_error,
            0.0,
            "absolute",
            convergence_flag_error,
            0.0,
        ),
        (
            "production_3d_restart_bitmatch_error",
            production_restart_error,
            0.0,
            "absolute",
            production_restart_error,
            0.0,
        ),
        (
            "production_3d_finite_diagnostics_error",
            production_finite_error,
            0.0,
            "absolute",
            production_finite_error,
            0.0,
        ),
        (
            "shock_campaign_3d_driver_contract_error",
            campaign_error,
            0.0,
            "absolute",
            campaign_error,
            0.0,
        ),
        (
            "compare_dims_shock_driver_contract_error",
            compare_dims_error,
            0.0,
            "absolute",
            compare_dims_error,
            0.0,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "shock_driver_sweeps",
        reference_kind = "analytic_or_contract",
        reference = "1D shock sweep flux-freezing/mass-conservation checks plus 3D campaign driver smoke contracts",
        rows = rows,
        artifact = artifact,
    )
end

VALIDATION_CASE = ValidationCase(
    id = "24_shock_driver_sweep_validation",
    default = false,
    description = "Perpendicular-shock sweep, convergence, 3D campaign restart, and dimension-comparison driver checks.",
    runner = case_24_shock_driver_sweep_validation,
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
