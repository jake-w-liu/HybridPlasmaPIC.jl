#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_06_boundary_loading_kdv_smoothing(artifact_dir::AbstractString)
    id = "06_boundary_loading_kdv_smoothing"

    ps_abs = ParticleSet{1,Float64}(4)
    ps_abs.x[1] .= [-0.1, 0.2, 0.8, 1.2]
    ps_abs.v[1] .= [1.0, 2.0, 3.0, 4.0]
    ps_abs.id .= UInt64[10, 20, 30, 40]
    removed = apply_absorbing!(ps_abs, (0.0,), (1.0,))
    absorbing_error = removed == 2 && ps_abs.id == UInt64[20, 30] ? 0.0 : 1.0

    gload = FourierGrid((10,), (2π,))
    ps_quiet = ParticleSet{1,Float64}(10)
    load_lattice_1d!(ps_quiet, 0.0, 2π)
    set_density_weight!(ps_quiet, 1.25, gload)
    load_quiet_velocities!(ps_quiet, MersenneTwister(91), (0.1, -0.2, 0.3), (1.0, 2.0, 3.0))
    quiet_error = maximum(
        abs(sum(ps_quiet.v[c]) / nparticles(ps_quiet) - (0.1, -0.2, 0.3)[c]) for c = 1:3
    )
    weight_error = abs(sum(ps_quiet.weight) - 1.25 * prod(gload.L))

    rng = MersenneTwister(5)
    ps_inj = ParticleSet{1,Float64}(0)
    acc = Ref(0.0)
    nextid = Ref(UInt64(1))
    ninj = inject_face_1d!(ps_inj, rng, 0.0, +1, 1.0, 2.0, 0.0, (0.0, 0.0), 0.3, 1.0, 0.5, acc, nextid)
    injection_error =
        ninj == 4 && nparticles(ps_inj) == 4 && all(==(2.0), ps_inj.v[1]) &&
        ps_inj.id == UInt64[1, 2, 3, 4] ? 0.0 : 1.0
    flux_error = abs(flux_per_density(2.0, 0.0) - 2.0)
    flux_speed_error = abs(flux_speed(MersenneTwister(6), 2.0, 0.0) - 2.0)

    gs = FourierGrid((32,), (2π,))
    mode = 3
    passes = 2
    xs = [(i - 1) * gs.dx[1] for i = 1:gs.n[1]]
    f = cos.(mode .* xs)
    f0 = copy(f)
    binomial_smooth!(f, gs; passes = passes)
    sfac = smoothing_transfer(mode, gs.dx[1]; passes = passes)
    smoothing_error = maximum(abs, f .- sfac .* f0)

    Ld = 40.0
    nkdv = 256
    xkdv = collect(range(0, Ld; length = nkdv + 1)[1:nkdv])
    c0, α, β, amp, x0 = 0.0, 6.0, 1.0, 1.0, 10.0
    tend = 1.0
    dt = 0.004
    u0 = kdv_soliton(amp, c0, α, β, xkdv, 0.0, x0, Ld)
    uf = kdv_solve(u0, Ld, c0, α, β, dt, round(Int, tend / dt))
    ua = kdv_soliton(amp, c0, α, β, xkdv, tend, x0, Ld)
    kdv_relerr = sqrt(sum(abs2, uf .- ua) / sum(abs2, ua))

    artifact = joinpath(artifact_dir, "06_boundary_loading_kdv_smoothing.csv")
    rows = (
        ("absorbing_boundary_compaction_error", absorbing_error, 0.0, "absolute", absorbing_error, 0.0),
        ("quiet_velocity_mean_max_abs_error", quiet_error, 0.0, "absolute", quiet_error, 1e-14),
        ("density_weight_total_abs_error", sum(ps_quiet.weight), 1.25 * prod(gload.L), "absolute", weight_error, 1e-12),
        ("cold_flux_per_density_abs_error", flux_per_density(2.0, 0.0), 2.0, "absolute", flux_error, 0.0),
        ("cold_flux_speed_abs_error", flux_speed(MersenneTwister(6), 2.0, 0.0), 2.0, "absolute", flux_speed_error, 0.0),
        ("cold_injection_batch_error", injection_error, 0.0, "absolute", injection_error, 0.0),
        ("binomial_smoothing_single_mode_error", smoothing_error, 0.0, "absolute", smoothing_error, 1e-12),
        ("kdv_soliton_relative_l2_error", kdv_relerr, 0.0, "relative", kdv_relerr, 2e-3),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "boundaries_loading_spectral_nonlinear",
        reference_kind = "analytic",
        reference = "absorbing compaction, cold flux injection, binomial transfer function, and KdV sech^2 soliton",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "06_boundary_loading_kdv_smoothing",
    default = true,
    description = "Open boundaries, loading, flux injection, binomial smoothing, and KdV soliton.",
    runner = case_06_boundary_loading_kdv_smoothing,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
