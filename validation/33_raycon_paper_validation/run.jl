#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

const RAYCON_TRACY_REFERENCE =
    "Tracy et al. 2001 Eq. (22); Tracy, Kaufman, and Jaun 2007 Eq. (58)"

const RAYCON_JAUN_REFERENCE =
    "Jaun et al. 2007 PPCF Figs. 3-4 and Sec. 3.1 D(20%He3) JET mid-plane case"

const RAYCON_JAUN_OFFMID_REFERENCE =
    "Jaun et al. 2007 PPCF Figs. 6-7 D(20%He3) off-mid-plane case"

const RAYCON_JAUN_DH_REFERENCE =
    "Jaun et al. 2007 PPCF Figs. 9-10 and Sec. 4 D(20%H) strong-coupling case"

function _jet_dhe3_paper_problem()
    # The Jaun Fig. 3 trajectory follows the bundled upstream jet06 preset:
    # iaspr=0.3 and ne=7e19. That historical figure-matching setup is not the
    # later quasi-neutral edit with ne=nD+2*nHe3.
    eq = Raycon.SolovevEquilibrium(; b0 = 3.35, r0 = 3.0, q0 = 1.5, iaspr = 0.3, elong = 1.4)
    nD = 5.6e19
    nHe3 = 1.4e19
    ne = 7.0e19
    return Raycon.RayconProblem(;
        eq,
        amass = [1 / 1836, 2.0, 3.0],
        acharge = [-1.0, 1.0, 2.0],
        n0 = [ne, nD, nHe3],
        na = [1.0, 1.0, 1.0],
        nb = [1.0, 1.0, 1.0],
        t0 = [5.0, 5.0, 5.0],
        ta = [0.6, 0.6, 0.6],
        tb = [1.0, 1.0, 1.0],
        freq = 37e6,
        kphi = -4.67,
        model = :cld2x2,
    )
end

function _jet_dh_strong_coupling_problem()
    # The published D-H text says to replace He3 by H and reduce B0. The
    # bundled upstream ray05 figure preset also uses a narrower H profile
    # (naH=0.5); using naH=1 produces a weak second event instead of Fig. 9/10's
    # negligibly small transmissions.
    eq = Raycon.SolovevEquilibrium(; b0 = 3.0, r0 = 3.0, q0 = 1.5, iaspr = 0.3, elong = 1.4)
    nD = 5.6e19
    nH = 1.4e19
    ne = 7.0e19
    return Raycon.RayconProblem(;
        eq,
        amass = [1 / 1836, 2.0, 1.0],
        acharge = [-1.0, 1.0, 1.0],
        n0 = [ne, nD, nH],
        na = [1.0, 1.0, 0.5],
        nb = [1.0, 1.0, 1.0],
        t0 = [5.0, 5.0, 5.0],
        ta = [0.6, 0.6, 0.6],
        tb = [1.0, 1.0, 1.0],
        freq = 37e6,
        kphi = -4.67,
        model = :cld2x2,
    )
end

function _flux_from_rz(eq, r::Real, z::Real)
    rho = sqrt((Float64(r) - eq.r0)^2 + Float64(z)^2)
    theta = atan(eq.r0 - Float64(r), -Float64(z)) + pi / 2
    s = Raycon.solovev_flux(eq, rho, theta).sflx
    return (; rho, theta, s)
end

function _push_abs_metric!(rows, metric::AbstractString, measured, expected, tolerance)
    m = Float64(measured)
    e = Float64(expected)
    err = isfinite(m) ? abs(m - e) : Inf
    push!(rows, (String(metric), m, e, "absolute", err, Float64(tolerance)))
    return rows
end

function _push_count_at_least_metric!(rows, metric::AbstractString, measured, expected_min)
    m = Float64(measured)
    e = Float64(expected_min)
    err = m >= e ? 0.0 : e - m
    push!(rows, (String(metric), m, e, "at_least", err, 0.0))
    return rows
end

function _push_at_most_metric!(rows, metric::AbstractString, measured, expected_max)
    m = Float64(measured)
    e = Float64(expected_max)
    err = isfinite(m) ? max(0.0, m - e) : Inf
    push!(rows, (String(metric), m, e, "at_most", err, 0.0))
    return rows
end

function _trace_rays_from_state(
    prob::Raycon.RayconProblem,
    y0::AbstractVector{<:Real};
    sigma_span::Real,
    max_conversions::Integer,
    rtol::Real = 1e-6,
    atol::Real = 1e-7,
)
    length(y0) == 4 || throw(ArgumentError("state must be (R, Z, kR, kZ)"))
    all(isfinite, y0) || throw(ArgumentError("state must be finite"))
    (isfinite(sigma_span) && sigma_span > 0) || throw(ArgumentError("sigma_span must be positive"))
    max_conversions >= 0 || throw(ArgumentError("max_conversions must be nonnegative"))

    span = Float64(sigma_span)
    queue = [(parent = 0, y = Float64.(collect(y0)), sigma0 = 0.0)]
    rays = @NamedTuple{parent::Int, sigma::Vector{Float64}, y::Matrix{Float64}, status::Symbol}[]
    conversions = @NamedTuple{ray::Int, sigma::Float64, conversion::Raycon.RayconConversion}[]

    while !isempty(queue)
        item = popfirst!(queue)
        segs_sigma = Float64[]
        segs_y = Matrix{Float64}(undef, 4, 0)
        y = collect(item.y)
        sigma = item.sigma0
        nconv = 0
        nseg = 0
        status = :end_of_span
        while sigma < span
            (nseg += 1) <= 50 || (status = :segment_limit; break)
            det = prob.model in (:cld2x2, :cld3x3) && nconv < max_conversions
            tr = Raycon.integrate_ray(prob, y, sigma, span; rtol, atol, detect_conversion = det)
            segs_sigma = vcat(segs_sigma, tr.sigma)
            segs_y = hcat(segs_y, tr.y)
            status = tr.status
            if tr.status !== :conversion_event
                break
            end
            sigma = tr.sigma[end]
            y = tr.y[:, end]
            conv = try
                Raycon.analyze_conversion(prob, y, tr.zdot, tr.zddot)
            catch e
                (e isa DomainError || e isa ArgumentError || e isa ErrorException) || rethrow()
                nothing
            end
            rayidx = length(rays) + 1
            if conv !== nothing && Raycon.is_valid(conv)
                nconv += 1
                push!(conversions, (; ray = rayidx, sigma, conversion = conv))
                push!(queue, (parent = rayidx, y = collect(conv.transmitted), sigma0 = sigma))
            end
            sigma >= span && break
        end
        push!(rays, (; parent = item.parent, sigma = segs_sigma, y = segs_y, status))
    end
    return (; rays, conversions)
end

function _write_phase_space_overlay!(
    path::AbstractString,
    result,
    conversions,
)
    rows = Tuple{String,Int,Float64,Float64,String}[]
    for (i, ray) in enumerate(result.rays)
        for j = 1:size(ray.y, 2)
            push!(rows, ("computed_ray_$i", j, ray.y[1, j], ray.y[3, j], "computed_ray"))
        end
    end
    for (i, c) in enumerate(conversions)
        push!(rows, ("computed_conversion_$i", 1, c.incoming[1], c.incoming[3], "incoming"))
        push!(rows, ("computed_conversion_$i", 2, c.saddle[1], c.saddle[3], "saddle"))
        push!(rows, ("computed_conversion_$i", 3, c.transmitted[1], c.transmitted[3], "transmitted"))
    end
    paper_points = (
        ("paper_label_1_launch", 3.50, -27.0),
        ("paper_label_2_incoming", 2.52, -11.0),
        ("paper_label_2prime_transmitted", 2.48, -27.0),
        ("paper_label_4_converted", 2.52, 11.0),
        ("paper_label_4prime_transmitted", 3.07, 27.0),
    )
    for (i, (label, r, kr)) in enumerate(paper_points)
        push!(rows, (label, i, r, kr, "paper_target"))
    end
    return _write_csv(
        path,
        ("series", "point_index", "R_m", "kR_inv_m", "kind"),
        rows,
    )
end

function _write_rz_trace_overlay!(path::AbstractString, result, paper_points)
    rows = Tuple{String,Int,Float64,Float64,String}[]
    for (i, ray) in enumerate(result.rays)
        for j = 1:size(ray.y, 2)
            push!(rows, ("computed_ray_$i", j, ray.y[1, j], ray.y[2, j], "computed_ray"))
        end
    end
    for entry in result.conversions
        c = entry.conversion
        i = length(rows) + 1
        push!(rows, ("conversion_$(entry.ray)_incoming", i, c.incoming[1], c.incoming[2], "incoming"))
        push!(rows, ("conversion_$(entry.ray)_saddle", i + 1, c.saddle[1], c.saddle[2], "saddle"))
        push!(
            rows,
            ("conversion_$(entry.ray)_transmitted", i + 2, c.transmitted[1], c.transmitted[2], "transmitted"),
        )
    end
    for (i, point) in enumerate(paper_points)
        label, r, z = point
        push!(rows, (String(label), i, Float64(r), Float64(z), "paper_target"))
    end
    return _write_csv(
        path,
        ("series", "point_index", "R_m", "Z_m", "kind"),
        rows,
    )
end

function _write_conversion_summary!(path::AbstractString, result)
    rows = map(enumerate(result.conversions)) do (i, entry)
        c = entry.conversion
        (
            i,
            entry.ray,
            entry.sigma,
            c.incoming[1],
            c.incoming[2],
            c.incoming[3],
            c.incoming[4],
            c.saddle[1],
            c.saddle[2],
            c.saddle[3],
            c.saddle[4],
            c.transmitted[1],
            c.transmitted[2],
            c.transmitted[3],
            c.transmitted[4],
            c.eta2,
            c.eta2_estimate,
            c.tau,
            abs(c.beta),
            c.converged,
            c.hyperbola_ok,
            c.transmitted_ok,
        )
    end
    return _write_csv(
        path,
        (
            "conversion_index",
            "ray_index",
            "sigma",
            "incoming_R_m",
            "incoming_Z_m",
            "incoming_kR_inv_m",
            "incoming_kZ_inv_m",
            "saddle_R_m",
            "saddle_Z_m",
            "saddle_kR_inv_m",
            "saddle_kZ_inv_m",
            "transmitted_R_m",
            "transmitted_Z_m",
            "transmitted_kR_inv_m",
            "transmitted_kZ_inv_m",
            "eta2",
            "eta2_estimate",
            "tau",
            "abs_beta",
            "converged",
            "hyperbola_ok",
            "transmitted_ok",
        ),
        rows,
    )
end

function _case_33_raycon_paper_validation(artifact_dir::AbstractString)
    id = "33_raycon_paper_validation"
    prob = _jet_dhe3_paper_problem()
    dh_prob = _jet_dh_strong_coupling_problem()

    launch = _flux_from_rz(prob.eq, 3.5, 0.0)
    launch_y = Raycon.launch_ray(prob; s = launch.s, theta = launch.theta, kr = -27.0, kz = 0.0)
    result = Raycon.trace_rays(
        prob;
        s = launch.s,
        theta = launch.theta,
        kr = -27.0,
        kz = 0.0,
        sigma_span = 3e-2,
        max_conversions = 2,
        rtol = 1e-7,
        atol = 1e-9,
    )

    offmid_y0 = [3.45, -0.31, -27.0, 0.0]
    offmid_result = _trace_rays_from_state(
        prob,
        offmid_y0;
        sigma_span = 4e-2,
        max_conversions = 2,
        rtol = 1e-7,
        atol = 1e-9,
    )

    dh_launch_y = Raycon.launch_ray(dh_prob; s = 0.6, theta = 0.4, kr = -7.5, kz = 0.0)
    dh_result = Raycon.trace_rays(
        dh_prob;
        s = 0.6,
        theta = 0.4,
        kr = -7.5,
        kz = 0.0,
        sigma_span = 8e-2,
        max_conversions = 4,
        rtol = 1e-7,
        atol = 1e-9,
    )

    conversions = [entry.conversion for entry in result.conversions]
    offmid_conversions = [entry.conversion for entry in offmid_result.conversions]
    dh_conversions = [entry.conversion for entry in dh_result.conversions]
    c1 = length(conversions) >= 1 ? conversions[1] : nothing
    c2 = length(conversions) >= 2 ? conversions[2] : nothing
    offmid_c1 = length(offmid_conversions) >= 1 ? offmid_conversions[1] : nothing
    offmid_c2 = length(offmid_conversions) >= 2 ? offmid_conversions[2] : nothing

    phase_space_artifact = joinpath(artifact_dir, "33_raycon_paper_validation_phase_space.csv")
    conversion_artifact = joinpath(artifact_dir, "33_raycon_paper_validation_conversions.csv")
    offmid_conversion_artifact =
        joinpath(artifact_dir, "33_raycon_paper_validation_fig6_7_conversions.csv")
    offmid_rz_artifact = joinpath(artifact_dir, "33_raycon_paper_validation_fig6_7_rz.csv")
    dh_conversion_artifact =
        joinpath(artifact_dir, "33_raycon_paper_validation_fig9_10_dh_conversions.csv")
    dh_rz_artifact = joinpath(artifact_dir, "33_raycon_paper_validation_fig9_10_dh_rz.csv")
    _write_phase_space_overlay!(phase_space_artifact, result, conversions)
    _write_conversion_summary!(conversion_artifact, result)
    _write_conversion_summary!(offmid_conversion_artifact, offmid_result)
    _write_rz_trace_overlay!(
        offmid_rz_artifact,
        offmid_result,
        (("paper_fig6_launch", offmid_y0[1], offmid_y0[2]),),
    )
    _write_conversion_summary!(dh_conversion_artifact, dh_result)
    _write_rz_trace_overlay!(
        dh_rz_artifact,
        dh_result,
        (("paper_fig9_launch", dh_launch_y[1], dh_launch_y[2]),),
    )

    rows = Tuple{String,Float64,Float64,String,Float64,Float64}[]
    _push_abs_metric!(rows, "jaun_fig3_projected_launch_R_m", launch_y[1], 3.5, 0.02)
    _push_abs_metric!(rows, "jaun_fig3_projected_launch_kR_inv_m", launch_y[3], -27.0, 6.0)
    _push_count_at_least_metric!(
        rows,
        "jaun_fig3_midplane_conversion_count",
        length(conversions),
        2,
    )

    if c1 === nothing
        for metric in (
            "tracy_tau_formula_max_abs_error",
            "tracy_beta_unitarity_max_abs_error",
            "tracy_transmitted_dispersion_max_abs",
            "jaun_fig3_first_incoming_R_m",
            "jaun_fig3_first_incoming_kR_inv_m",
            "jaun_fig3_first_transmitted_R_m",
            "jaun_fig3_first_transmitted_kR_inv_m",
            "jaun_fig4_first_tau",
            "jaun_fig4_first_eta2_over_B",
        )
            _push_abs_metric!(rows, metric, NaN, 0.0, 0.0)
        end
    else
        tau_formula_errors = [
            abs(c.tau - exp(-pi * c.eta2)) for c in conversions if isfinite(c.tau) && isfinite(c.eta2)
        ]
        beta_unitarity_errors = [
            abs(c.tau^2 + abs2(c.beta) - 1.0) for c in conversions if isfinite(c.tau) && isfinite(abs(c.beta))
        ]
        transmitted_U = [
            abs(Raycon.dispersion_U(prob, collect(c.transmitted))) for c in conversions
        ]

        _push_abs_metric!(
            rows,
            "tracy_tau_formula_max_abs_error",
            isempty(tau_formula_errors) ? NaN : maximum(tau_formula_errors),
            0.0,
            1e-12,
        )
        _push_abs_metric!(
            rows,
            "tracy_beta_unitarity_max_abs_error",
            isempty(beta_unitarity_errors) ? NaN : maximum(beta_unitarity_errors),
            0.0,
            1e-10,
        )
        _push_abs_metric!(
            rows,
            "tracy_transmitted_dispersion_max_abs",
            isempty(transmitted_U) ? NaN : maximum(transmitted_U),
            0.0,
            1e-8,
        )

        _push_abs_metric!(rows, "jaun_fig3_first_incoming_R_m", c1.incoming[1], 2.52, 0.12)
        _push_abs_metric!(rows, "jaun_fig3_first_incoming_kR_inv_m", c1.incoming[3], -11.0, 6.0)
        _push_abs_metric!(rows, "jaun_fig3_first_transmitted_R_m", c1.transmitted[1], 2.48, 0.12)
        _push_abs_metric!(
            rows,
            "jaun_fig3_first_transmitted_kR_inv_m",
            c1.transmitted[3],
            -27.0,
            6.0,
        )
        _push_abs_metric!(rows, "jaun_fig4_first_tau", c1.tau, 0.27, 0.08)
        _push_abs_metric!(rows, "jaun_fig4_first_eta2_over_B", c1.eta2, 0.41, 0.12)
    end

    if c1 !== nothing && c2 !== nothing
        _push_abs_metric!(rows, "jaun_fig4_midplane_tau_symmetry", abs(c1.tau - c2.tau), 0.0, 0.03)
        _push_abs_metric!(rows, "jaun_fig4_second_tau", c2.tau, 0.27, 0.08)
    else
        _push_abs_metric!(rows, "jaun_fig4_midplane_tau_symmetry", NaN, 0.0, 0.03)
        _push_abs_metric!(rows, "jaun_fig4_second_tau", NaN, 0.27, 0.08)
    end

    _push_abs_metric!(rows, "jaun_fig6_launch_R_m", offmid_y0[1], 3.45, 0.01)
    _push_abs_metric!(rows, "jaun_fig6_launch_Z_m", offmid_y0[2], -0.31, 0.01)
    _push_abs_metric!(rows, "jaun_fig6_launch_kR_inv_m", offmid_y0[3], -27.0, 0.01)
    _push_abs_metric!(rows, "jaun_fig6_launch_kZ_inv_m", offmid_y0[4], 0.0, 0.01)
    _push_count_at_least_metric!(
        rows,
        "jaun_fig6_7_off_midplane_conversion_count",
        length(offmid_conversions),
        2,
    )
    if offmid_c1 === nothing
        for metric in (
            "jaun_fig7_off_midplane_first_incoming_kZ_inv_m",
            "jaun_fig7_off_midplane_first_tau",
            "jaun_fig7_off_midplane_second_tau",
        )
            _push_abs_metric!(rows, metric, NaN, 0.0, 0.0)
        end
    else
        _push_abs_metric!(
            rows,
            "jaun_fig7_off_midplane_first_incoming_kZ_inv_m",
            offmid_c1.incoming[4],
            6.8,
            1.5,
        )
        _push_abs_metric!(rows, "jaun_fig7_off_midplane_first_tau", offmid_c1.tau, 0.25, 0.08)
        if offmid_c2 === nothing
            _push_abs_metric!(rows, "jaun_fig7_off_midplane_second_tau", NaN, 0.28, 0.08)
        else
            _push_abs_metric!(rows, "jaun_fig7_off_midplane_second_tau", offmid_c2.tau, 0.28, 0.08)
        end
    end

    _push_abs_metric!(rows, "jaun_fig9_dh_launch_R_m", dh_launch_y[1], 3.476, 0.03)
    _push_abs_metric!(rows, "jaun_fig9_dh_launch_Z_m", dh_launch_y[2], -0.201, 0.03)
    _push_count_at_least_metric!(
        rows,
        "jaun_fig9_10_dh_strong_coupling_conversion_count",
        length(dh_conversions),
        4,
    )
    if length(dh_conversions) >= 4
        first_four = dh_conversions[1:4]
        _push_abs_metric!(rows, "jaun_fig9_dh_first_eta2_over_B", first_four[1].eta2, 7.0, 0.8)
        _push_abs_metric!(rows, "jaun_fig9_dh_second_eta2_over_B", first_four[2].eta2, 5.2, 1.8)
        _push_count_at_least_metric!(
            rows,
            "jaun_fig9_10_dh_min_eta2_over_B_first4",
            minimum(c.eta2 for c in first_four),
            4.6,
        )
        _push_at_most_metric!(
            rows,
            "jaun_fig9_10_dh_max_tau_first4",
            maximum(c.tau for c in first_four),
            1e-7,
        )
    else
        for metric in (
            "jaun_fig9_dh_first_eta2_over_B",
            "jaun_fig9_dh_second_eta2_over_B",
            "jaun_fig9_10_dh_min_eta2_over_B_first4",
            "jaun_fig9_10_dh_max_tau_first4",
        )
            _push_abs_metric!(rows, metric, NaN, 0.0, 0.0)
        end
    end

    artifact = joinpath(artifact_dir, "33_raycon_paper_validation.csv")
    _write_metric_csv(artifact, rows)

    notes =
        "Jaun figure targets were visually inspected from the local PDF; Tracy metrics are exact " *
        "connection-coefficient identities. Jaun coordinates are approximate figure/text targets. " *
        "This case uses the historical upstream jet06 density/aspect-ratio preset that matches Fig. 3. " *
        "Figs. 6-7 use the published fixed lab-frame off-mid-plane launch, and Figs. 9-10 use the " *
        "D(20%H), B0=3.0 T strong-coupling setup with the upstream ray05 H-profile choice. " *
        "Additional phase-space, R-Z, and conversion-summary CSV artifacts are written for figure comparison."
    gated = _metric_rows_to_results(
        id = id,
        category = "raycon",
        reference_kind = "published_external_paper",
        reference =
            RAYCON_TRACY_REFERENCE * "; " *
            RAYCON_JAUN_REFERENCE * "; " *
            RAYCON_JAUN_OFFMID_REFERENCE * "; " *
            RAYCON_JAUN_DH_REFERENCE,
        rows = rows,
        artifact = artifact,
        notes = notes,
    )
    return gated
end


VALIDATION_CASE = ValidationCase(
    id = "33_raycon_paper_validation",
    default = false,
    description = "Raycon validation against Jaun/Tracy mode-conversion papers.",
    runner = _case_33_raycon_paper_validation,
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
