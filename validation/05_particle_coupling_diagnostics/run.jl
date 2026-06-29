#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_05_particle_coupling_diagnostics(artifact_dir::AbstractString)
    id = "05_particle_coupling_diagnostics"
    g = FourierGrid((12, 10), (2π, 2π))
    ps = ParticleSet{2,Float64}(prod(g.n))
    load_lattice!(ps, (0.0, 0.0), g.L, g.n)
    set_density_weight!(ps, 1.0, g)
    density_field = zeros(Float64, g.n)
    density!(density_field, ps, g, CIC())
    particle_weight = sum(ps.weight)
    density_integral = sum(density_field) * prod(g.dx)
    density_relerr = abs(density_integral - particle_weight) / particle_weight

    field = zeros(Float64, g.n)
    a = 0.4
    bx, by = 1.2, -0.7
    probe = ParticleSet{2,Float64}(16)
    idx = 1
    for j = 1:4, i = 1:4
        probe.x[1][idx] = (2 + i) * g.dx[1]
        probe.x[2][idx] = (2 + j) * g.dx[2]
        idx += 1
    end
    for index in CartesianIndices(field)
        i, j = Tuple(index)
        x = (i - 1) * g.dx[1]
        y = (j - 1) * g.dx[2]
        field[index] = a + bx * x + by * y
    end
    gathered = zeros(Float64, nparticles(probe))
    gather_scalar!(gathered, field, probe, g, CIC())
    exact = [a + bx * probe.x[1][p] + by * probe.x[2][p] for p = 1:nparticles(probe)]
    gather_error = maximum(abs, gathered .- exact)

    rng = MersenneTwister(23)
    vals = randn(rng, nparticles(ps))
    random_field = randn(rng, g.n...)
    dep = zeros(Float64, g.n)
    gat = zeros(Float64, nparticles(ps))
    deposit_scalar!(dep, ps, vals, g, CIC())
    gather_scalar!(gat, random_field, ps, g, CIC())
    adjoint_scale = abs(dot(vec(dep), vec(random_field))) + abs(dot(vals, gat)) + eps(Float64)
    adjoint_relerr = abs(dot(vec(dep), vec(random_field)) - dot(vals, gat)) / adjoint_scale

    onecell = FourierGrid((1,), (1.0,))
    pair = ParticleSet{1,Float64}(2; m = 2.0)
    pair.x[1] .= (0.0, 0.0)
    pair.weight .= 1.0
    pair.v[1] .= (-1.0, 3.0)
    pair.v[2] .= (2.0, 4.0)
    pair.v[3] .= (5.0, 1.0)
    P = ntuple(_ -> zeros(Float64, 1), 6)
    pressure_tensor!(P, pair, onecell, NGP())
    expected_pressure = (16.0, 4.0, 16.0, 8.0, -16.0, -8.0)
    pressure_error = maximum(abs(P[c][1] - expected_pressure[c]) for c = 1:6)

    _, hist = velocity_histogram(pair, 1; nbins = 4)
    histogram_relerr = abs(sum(hist) - sum(pair.weight)) / sum(pair.weight)
    _, _, phase_hist = phase_space_histogram(pair, 1, 1; nx = 4, nv = 4)
    phase_histogram_relerr = abs(sum(phase_hist) - sum(pair.weight)) / sum(pair.weight)

    artifact = joinpath(artifact_dir, "05_particle_coupling_diagnostics.csv")
    rows = (
        ("density_integral_relative_error", density_relerr, 0.0, "relative", density_relerr, 1e-12),
        ("cic_linear_gather_max_abs_error", gather_error, 0.0, "absolute", gather_error, 1e-12),
        (
            "deposit_gather_adjoint_relative_error",
            adjoint_relerr,
            0.0,
            "relative",
            adjoint_relerr,
            1e-12,
        ),
        (
            "pressure_tensor_exact_max_abs_error",
            pressure_error,
            0.0,
            "absolute",
            pressure_error,
            1e-12,
        ),
        (
            "velocity_histogram_weight_relative_error",
            histogram_relerr,
            0.0,
            "relative",
            histogram_relerr,
            1e-12,
        ),
        (
            "phase_space_histogram_weight_relative_error",
            phase_histogram_relerr,
            0.0,
            "relative",
            phase_histogram_relerr,
            1e-12,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "particles_coupling_diagnostics",
        reference_kind = "analytic",
        reference = "partition of unity, CIC linear reproduction, adjoint identity, centered pressure tensor, histogram weight conservation",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "05_particle_coupling_diagnostics",
    default = true,
    description = "Particle deposition/gather/pressure/histogram diagnostics against exact invariants.",
    runner = case_05_particle_coupling_diagnostics,
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
