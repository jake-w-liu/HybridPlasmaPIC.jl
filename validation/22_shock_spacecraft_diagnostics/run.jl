#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_22_shock_spacecraft_diagnostics(artifact_dir::AbstractString)
    id = "22_shock_spacecraft_diagnostics"

    xf, width = shock_front([2.0, 2.0, 1.0, 1.0], [0.0, 1.0, 2.0, 3.0])
    shock_front_error = max(abs(xf - 2.0), abs(width - 1.0))

    ps = ParticleSet{1,Float64}(3)
    ps.id .= UInt64[10, 20, 30]
    ps.m = 1.0
    logger = CrossingLogger(Float64)
    ps.x[1] .= [4.0, 4.0, 6.0]
    ps.v[1] .= [1.0, 0.0, 0.0]
    ps.v[2] .= [0.0, 0.0, 0.0]
    ps.v[3] .= [0.0, 0.0, 0.0]
    log_crossings!(logger, ps, 5.0)
    ps.x[1] .= [6.0, 4.5, 4.0]
    ps.v[1] .= [2.0, 0.0, 0.0]
    ps.v[2] .= [0.0, 0.0, 3.0]
    ps.v[3] .= [0.0, 0.0, 0.0]
    log_crossings!(logger, ps, 5.0)
    ps.x[1] .= [4.0, 4.6, 3.5]
    ps.v[1] .= [0.0, 0.0, 0.0]
    ps.v[2] .= [0.0, 0.0, 0.0]
    ps.v[3] .= [1.0, 0.0, 0.0]
    log_crossings!(logger, ps, 5.0)
    crossing_error = abs(crossing_count(logger) - 3)
    gain_error = abs(energy_gain(logger) - 4.5)

    n_hat = (1.0, 0.0, 0.0)
    u = (-3.0, 1.7, -0.4)
    B = (0.0, 0.0, 1.0)
    shock_frame_error = abs(shock_frame(-3.0, 1.25) + 4.25)
    Vnif = normal_incidence_frame(u, B, n_hat)
    urel = (u[1] - Vnif[1], u[2] - Vnif[2], u[3] - Vnif[3])
    nif_tangent_error = hypot(urel[2], urel[3])

    sh = PerpShock(24, 12.0; B0 = 1.0)
    dx = sh.s.dx
    edge = sh.x[end] - 3dx
    pband = ParticleSet{1,Float64}(6)
    pband.x[1] .= [sh.x[end], sh.x[end] - dx, sh.x[end] - 2dx, edge - 0.5, 1.0, 5.0]
    pband.v[1] .= [1.0, 2.0, -1.0, 5.0, 1.0, 1.0]
    reflection_error = abs(boundary_reflection_fraction(sh, pband; ncells = 3) - 2 / 3)

    n = length(sh.x)
    sh.Bz[1] = 2.0
    sh.Bz[n] = 1.0
    sh.Ey[1] = 0.5
    sh.Ey[n] = -3.0
    sh.n[1] = 4.0
    sh.n[n] = 1.0
    sh.ux[1] = 0.0
    sh.ux[n] = -3.0
    sh.uy[1] = 0.1
    sh.uy[n] = 0.2
    sh.pe[1] = 3.0
    sh.pe[n] = 2.0
    flux = boundary_energy_flux(sh)
    gamma_factor = sh.γe / (sh.γe - 1)
    expected_total_inflow =
        -3.0 + 0.5 * 1.0 * ((-3.0)^2 + 0.2^2) * (-3.0) + gamma_factor * 2.0 * (-3.0)
    expected_total_wall = 1.0
    flux_error =
        max(abs(flux.total[1] - expected_total_inflow), abs(flux.total[2] - expected_total_wall))

    g = FourierGrid((8,), (8.0,))
    field = [(i - 1) * g.dx[1] for i = 1:g.n[1]]
    gather_error = abs(gather_at(field, g, 2.25) - 2.25)
    probe = SyntheticProbe(2.25)
    sample_error = abs(sample!(probe, field, g, 0.0) - 2.25)
    advance!(probe, 0.5, 2.0)
    probe_error = max(sample_error, abs(probe.x - 3.25))
    crossing_time_error = abs(crossing_time([0.0, 1.0], [0.0, 2.0], 1.0) - 0.5)

    normal = ntuple(_ -> 1 / sqrt(3), 3)
    speed = 2.0
    positions = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    times = ntuple(i -> dot(normal, positions[i]) / speed, 4)
    timing = four_spacecraft_timing(positions, times)
    timing_error =
        max(maximum(abs(timing.normal[i] - normal[i]) for i = 1:3), abs(timing.speed - speed))

    vht = dehoffmann_teller_velocity((2.0, 3.0, 4.0), (0.0, 0.0, 2.0))
    dht_residual_error = hypot(2.0 - vht[1], 3.0 - vht[2])

    pref = ParticleSet{1,Float64}(3)
    pref.x[1] .= [6.0, 6.0, 4.0]
    pref.v[1] .= [2.0, -1.0, 3.0]
    flags = classify_reflected(pref, 5.0, 0.0)
    classify_error = flags == [true, false, false] ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "22_shock_spacecraft_diagnostics.csv")
    rows = (
        (
            "shock_front_position_width_max_abs_error",
            max(xf, width),
            2.0,
            "absolute",
            shock_front_error,
            1e-12,
        ),
        (
            "crossing_logger_count_error",
            crossing_count(logger),
            3.0,
            "absolute",
            crossing_error,
            0.0,
        ),
        (
            "crossing_logger_energy_gain_error",
            energy_gain(logger),
            4.5,
            "absolute",
            gain_error,
            1e-12,
        ),
        (
            "shock_frame_velocity_abs_error",
            shock_frame(-3.0, 1.25),
            -4.25,
            "absolute",
            shock_frame_error,
            0.0,
        ),
        (
            "normal_incidence_tangent_residual",
            nif_tangent_error,
            0.0,
            "absolute",
            nif_tangent_error,
            1e-12,
        ),
        (
            "boundary_reflection_fraction_abs_error",
            boundary_reflection_fraction(sh, pband; ncells = 3),
            2 / 3,
            "absolute",
            reflection_error,
            1e-12,
        ),
        (
            "boundary_energy_flux_total_abs_error",
            maximum(abs, flux.total),
            maximum(abs, (expected_total_inflow, expected_total_wall)),
            "absolute",
            flux_error,
            1e-12,
        ),
        (
            "gather_at_linear_interpolation_abs_error",
            gather_at(field, g, 2.25),
            2.25,
            "absolute",
            gather_error,
            1e-12,
        ),
        ("synthetic_probe_sample_advance_abs_error", probe.x, 3.25, "absolute", probe_error, 1e-12),
        (
            "crossing_time_linear_abs_error",
            crossing_time([0.0, 1.0], [0.0, 2.0], 1.0),
            0.5,
            "absolute",
            crossing_time_error,
            1e-12,
        ),
        (
            "four_spacecraft_timing_max_abs_error",
            timing.speed,
            speed,
            "absolute",
            timing_error,
            1e-12,
        ),
        (
            "dehoffmann_teller_perp_residual",
            dht_residual_error,
            0.0,
            "absolute",
            dht_residual_error,
            1e-12,
        ),
        ("classify_reflected_boolean_error", sum(flags), 1.0, "absolute", classify_error, 0.0),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "shock_diagnostics",
        reference_kind = "analytic",
        reference = "constructed shock surfaces, particle crossings, boundary fluxes, and spacecraft timing geometry",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "22_shock_spacecraft_diagnostics",
    default = true,
    description = "Shock diagnostics and synthetic-spacecraft geometry against constructed references.",
    runner = case_22_shock_spacecraft_diagnostics,
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
