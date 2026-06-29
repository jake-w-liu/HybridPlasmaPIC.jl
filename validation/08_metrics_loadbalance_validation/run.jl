#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_08_metrics_loadbalance_validation(artifact_dir::AbstractString)
    id = "08_metrics_loadbalance_validation"
    g = FourierGrid((16,), (2π,))
    n = 128
    ps = ParticleSet{1,Float64}(n; q = 1.7, m = 1.0)
    rng = MersenneTwister(72)
    load_uniform!(ps, rng, (0.0,), g.L)
    for p = 1:n
        ps.weight[p] = 0.5 + rand(rng)
        ps.v[1][p] = randn(rng)
        ps.v[2][p] = randn(rng)
        ps.v[3][p] = randn(rng)
    end
    E0 = (0.3, -0.45, 0.8)
    Egrid = ntuple(c -> fill(Float64(E0[c]), g.n), 3)
    dt = 0.05
    work_error = 0.0
    for shape in (NGP(), CIC(), TSC())
        work = zeros(Float64, n)
        particle_work!(work, ps, Egrid, g, shape, dt)
        ref =
            [ps.q * (ps.v[1][p] * E0[1] + ps.v[2][p] * E0[2] + ps.v[3][p] * E0[3]) * dt for p = 1:n]
        work_error = max(work_error, maximum(abs, work .- ref))
    end

    residuals = [mixed_divcurl_residual(SBP1D(nx, 1.0), nx, 16, 2π) for nx in (33, 65, 129)]
    mixed_monotone_error = max(0.0, residuals[2] - residuals[1], residuals[3] - residuals[2])

    percell = [9, 1, 1, 1, 1, 9]
    ranges = balanced_tile_ranges(percell, 3)
    loads = balanced_tile_loads(percell, ranges)
    balance_cap_error = abs(maximum(loads) - 9)
    coverage = Int[]
    for r in ranges
        append!(coverage, collect(r))
    end
    coverage_error = coverage == collect(1:length(percell)) ? 0.0 : 1.0

    sort_grid = FourierGrid((4, 3), (4.0, 3.0))
    ps2 = ParticleSet{2,Float64}(5)
    ps2.x[1] .= [3.5, 0.1, 2.2, 1.1, 0.2]
    ps2.x[2] .= [2.5, 0.1, 1.4, 2.1, 1.2]
    ps2.id .= UInt64[50, 10, 30, 20, 40]
    ci = cell_index(ps2, sort_grid)
    expected_ci = [12, 1, 7, 10, 5]
    cell_index_error = maximum(abs, ci .- expected_ci)
    ids0 = sort(ps2.id)
    sort_particles!(ps2, sort_grid)
    sort_error = issorted(cell_index(ps2, sort_grid)) && sort(ps2.id) == ids0 ? 0.0 : 1.0

    mem = memory_bytes(; ncells = 10, nppc = 4, nspecies = 2, D = 3, Tbytes = 8)
    memory_error = abs(mem - 4800)

    artifact = joinpath(artifact_dir, "08_metrics_loadbalance_validation.csv")
    rows = (
        (
            "particle_work_uniform_field_max_abs_error",
            work_error,
            0.0,
            "absolute",
            work_error,
            1e-12,
        ),
        ("mixed_divcurl_finest_residual", residuals[end], 0.0, "absolute", residuals[end], 1e-2),
        (
            "mixed_divcurl_monotonicity_error",
            mixed_monotone_error,
            0.0,
            "absolute",
            mixed_monotone_error,
            0.0,
        ),
        (
            "balanced_tile_minimax_cap_error",
            maximum(loads),
            9.0,
            "absolute",
            balance_cap_error,
            0.0,
        ),
        ("balanced_tile_coverage_error", coverage_error, 0.0, "absolute", coverage_error, 0.0),
        (
            "cell_index_oracle_max_abs_error",
            maximum(ci),
            maximum(expected_ci),
            "absolute",
            cell_index_error,
            0.0,
        ),
        ("sort_particles_ordering_error", sort_error, 0.0, "absolute", sort_error, 0.0),
        ("memory_bytes_exact_error", mem, 4800.0, "absolute", memory_error, 0.0),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "closure_loadbalance",
        reference_kind = "analytic",
        reference = "uniform-field work identity, SBP/Fourier mixed residual convergence, brute cell/load/memory arithmetic",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "08_metrics_loadbalance_validation",
    default = true,
    description = "Particle work, mixed-grid residual, sorting, load-balance, and memory arithmetic.",
    runner = case_08_metrics_loadbalance_validation,
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
