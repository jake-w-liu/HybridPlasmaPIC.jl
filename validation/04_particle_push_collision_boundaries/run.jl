#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_04_particle_push_collision_boundaries(artifact_dir::AbstractString)
    id = "04_particle_push_collision_boundaries"
    ps = ParticleSet{2,Float64}(1)
    ps.v[1][1] = 1.0
    speed0 = sqrt(ps.v[1][1]^2 + ps.v[2][1]^2 + ps.v[3][1]^2)
    dt = 0.025
    for _ = 1:800
        push_uniform!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), dt)
    end
    speed = sqrt(ps.v[1][1]^2 + ps.v[2][1]^2 + ps.v[3][1]^2)
    speed_relerr = abs(speed - speed0) / speed0

    wrapped = ParticleSet{2,Float64}(3)
    wrapped.x[1] .= (-0.25, 0.5, 2.25)
    wrapped.x[2] .= (0.25, -0.5, 2.5)
    apply_periodic!(wrapped, (0.0, 0.0), (2.0, 2.0))
    periodic_expected_x = [1.75, 0.5, 0.25]
    periodic_expected_y = [0.25, 1.5, 0.5]
    periodic_error = max(
        maximum(abs, wrapped.x[1] .- periodic_expected_x),
        maximum(abs, wrapped.x[2] .- periodic_expected_y),
    )

    reflected = ParticleSet{1,Float64}(2)
    reflected.x[1] .= (-0.2, 1.2)
    reflected.v[1] .= (-3.0, 4.0)
    apply_reflecting!(reflected, (0.0,), (1.0,))
    reflecting_error =
        max(maximum(abs, reflected.x[1] .- [0.2, 0.8]), maximum(abs, reflected.v[1] .- [3.0, -4.0]))

    rng = MersenneTwister(31)
    coll = ParticleSet{1,Float64}(4_000)
    load_maxwellian!(coll, rng, (0.4, -0.3, 0.2), (1.0, 1.6, 0.7))
    coll.weight .= 0.5 .+ rand(rng, nparticles(coll))
    P0, E0 = _weighted_totals(coll)
    collide_bgk!(coll, 5.0, 0.1; rng = MersenneTwister(99))
    P1, E1 = _weighted_totals(coll)
    pscale = sqrt(sum(abs2, P0)) + 1.0
    momentum_relerr = sqrt(sum(abs2, (P1[1] - P0[1], P1[2] - P0[2], P1[3] - P0[3]))) / pscale
    energy_relerr = abs(E1 - E0) / E0

    artifact = joinpath(artifact_dir, "04_particle_push_collision_boundaries.csv")
    rows = (
        (
            "boris_uniform_B_speed_relative_error",
            speed_relerr,
            0.0,
            "relative",
            speed_relerr,
            1e-12,
        ),
        ("periodic_boundary_max_abs_error", periodic_error, 0.0, "absolute", periodic_error, 1e-12),
        (
            "reflecting_boundary_max_abs_error",
            reflecting_error,
            0.0,
            "absolute",
            reflecting_error,
            1e-12,
        ),
        ("bgk_momentum_relative_error", momentum_relerr, 0.0, "relative", momentum_relerr, 1e-10),
        ("bgk_energy_relative_error", energy_relerr, 0.0, "relative", energy_relerr, 1e-10),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "particles_boundaries_collisions",
        reference_kind = "analytic",
        reference = "Boris magnetic speed conservation, exact boundary maps, BGK conservation laws",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "04_particle_push_collision_boundaries",
    default = true,
    description = "Boris pusher, particle boundaries, and BGK conservation laws.",
    runner = case_04_particle_push_collision_boundaries,
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
