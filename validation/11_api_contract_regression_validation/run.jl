#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_11_api_contract_regression_validation(artifact_dir::AbstractString)
    id = "11_api_contract_regression_validation"

    layout = LogicalRankLayout((2, 3); periodic = (true, false))
    rng_a = rand(rank_rng(1234, layout, 4; stream = 2), 6)
    rng_b = rand(MersenneTwister(rank_seed(1234, layout, 4; stream = 2)), 6)
    seed_rng_error =
        rank_seed(1234, layout, 4) == rank_seed(1234, layout, 4) &&
        rank_seed(1234, layout, 4) != rank_seed(1234, layout, 5) &&
        rng_a == rng_b ? 0.0 : 1.0

    ids = ParticleSet{1,Float64}(3)
    assign_global_particle_ids!(ids, [2, 5, 8]; species = 1)
    expected_ids = UInt64[
        global_particle_id(2; species = 1),
        global_particle_id(5; species = 1),
        global_particle_id(8; species = 1),
    ]
    id_error =
        ids.id == expected_ids &&
        global_particle_id(42; species = 3) == (UInt64(3) << 48) | UInt64(42) ? 0.0 : 1.0

    nfd = 33
    s = SBP1D(nfd, 1.0)
    xfd = collect(range(0.0, 1.0; length = nfd))
    u = sin.(2 .* xfd) .+ 0.5 .* xfd .^ 2
    v = cos.(3 .* xfd) .- 0.25 .* xfd
    Du = sbp_deriv(u, s)
    Dv = similar(v)
    sbp_deriv!(Dv, v, s)
    sbp_identity_error =
        abs(sum(s.H .* u .* Dv) + sum(s.H .* Du .* v) - (u[end] * v[end] - u[1] * v[1]))

    ly = 2π
    ny = 16
    ky = 3.0
    F = [sin(2 * xfd[i]) * cos(ky * (j - 1) * ly / ny) for i = 1:nfd, j = 1:ny]
    Fy = similar(F)
    ywork = FourierDerivYWorkspace(F, ly)
    fourier_deriv_y!(Fy, F, ywork)
    Fy_expected = [-ky * sin(2 * xfd[i]) * sin(ky * (j - 1) * ly / ny) for i = 1:nfd, j = 1:ny]
    mixed_derivative_error = maximum(abs, Fy .- Fy_expected)
    Fx_linear = [2.0 * xfd[i] + 0.25 * j for i = 1:nfd, j = 1:ny]
    Dx_linear = similar(Fx_linear)
    sbp_deriv_x!(Dx_linear, Fx_linear, s)
    sbp_x_error = maximum(abs, Dx_linear .- 2.0)

    gs = FourierGrid((48,), (2π,))
    xs = [(i - 1) * gs.dx[1] for i = 1:gs.n[1]]
    c = fill(3.7, gs.n)
    c0 = copy(c)
    exp_filter!(c, gs)
    filter_error = maximum(abs, c .- c0)
    high = cos.(20 .* xs)
    dealias_two_thirds!(high, gs)
    dealias_error = maximum(abs, high)
    friendly_error = fft_friendly_size(17) == 18 && fft_friendly_size(101) == 105 ? 0.0 : 1.0
    namespace_error = SpectralOperators.deriv! === HybridPlasmaPIC.deriv! ? 0.0 : 1.0
    wisdom_error = mktempdir() do dir
        path = joinpath(dir, "wisdom.dat")
        val = with_fftw_wisdom(path) do
            plan_fft!(zeros(ComplexF64, 8))
            42
        end
        val == 42 && isfile(path) ? 0.0 : 1.0
    end
    spectral_contract_error =
        maximum((filter_error, dealias_error, friendly_error, namespace_error, wisdom_error))

    hf = HybridFields{2,Float64}((3, 4))
    nvals = [1.0, 2.0, 3.0]
    pe = similar(nvals)
    closure = PolytropicElectrons(2.0, 1.0, 2.0)
    electron_pressure!(pe, nvals, closure)
    fields_closure_error = maximum((
        maximum(abs, pe .- 2.0 .* nvals .^ 2),
        abs(closure_gamma(closure) - 2.0),
        abs(closure_gamma(IsothermalElectrons(1.0)) - 1.0),
        maximum(abs, hf.n),
        Float64(hf.floor_count[]),
    ),)

    gvec = FourierGrid((4, 4), (1.0, 1.0))
    psg = ParticleSet{2,Float64}(2)
    psg.x[1] .= [0.1, 0.7]
    psg.x[2] .= [0.2, 0.8]
    field = (fill(1.0, gvec.n), fill(2.0, gvec.n), fill(3.0, gvec.n))
    gathered = (zeros(2), zeros(2), zeros(2))
    gather_vector!(gathered, field, psg, gvec, CIC())
    gather_vector_error = maximum((
        maximum(abs, gathered[1] .- 1.0),
        maximum(abs, gathered[2] .- 2.0),
        maximum(abs, gathered[3] .- 3.0),
    ),)

    gcur = FourierGrid((4,), (4.0,))
    pscur = ParticleSet{1,Float64}(2; q = 2.0, m = 3.0)
    pscur.x[1] .= [0.2, 1.2]
    pscur.weight .= [1.0, 2.0]
    pscur.v[1] .= [3.0, -1.0]
    pscur.v[2] .= [0.5, 2.0]
    pscur.v[3] .= [-2.0, 1.0]
    mom = ntuple(_ -> zeros(Float64, gcur.n), 3)
    cur = ntuple(_ -> zeros(Float64, gcur.n), 3)
    momentum!(mom, pscur, gcur, NGP())
    current!(cur, pscur, gcur, NGP())
    P = (
        fill(4.0, gcur.n),
        fill(6.0, gcur.n),
        fill(8.0, gcur.n),
        zeros(gcur.n),
        zeros(gcur.n),
        zeros(gcur.n),
    )
    tc = temperature_components(P, fill(2.0, gcur.n))
    moment_temperature_error = maximum((
        maximum(abs, cur[1] .- 2.0 .* mom[1]),
        maximum(abs, cur[2] .- 2.0 .* mom[2]),
        maximum(abs, cur[3] .- 2.0 .* mom[3]),
        maximum(abs, tc[1] .- 2.0),
        maximum(abs, tc[2] .- 3.0),
        maximum(abs, tc[3] .- 4.0),
    ),)

    up = MHDState(1.0, 4.0, 0.0, 0.1, 0.0, 1.0)
    rh = rankine_hugoniot(up, 5 / 3)
    branches = rh_branches(up, 5 / 3; nscan = 256)
    rh_error = maximum(abs(Float64(getfield(rh.residuals, k))) for k in keys(rh.residuals))
    rh_contract_error = max(rh_error, rh.X > 1.0 && length(branches) == 1 ? 0.0 : 1.0)

    local_values = [
        (energy = (kinetic = 1.0, magnetic = 2.0), momentum = (1.0, 2.0, 3.0), hist = [1, 0]),
        (energy = (kinetic = 4.0, magnetic = 5.0), momentum = (-1.0, 1.0, 0.0), hist = [0, 3]),
    ]
    reduced = sum_diagnostics(local_values)
    minv = min_diagnostics([3.0, -2.0, 5.0])
    maxv = max_diagnostics([3.0, -2.0, 5.0])
    reduction_error =
        reduced.energy == (kinetic = 5.0, magnetic = 7.0) &&
        reduced.momentum == (0.0, 3.0, 3.0) &&
        reduced.hist == [1, 3] &&
        minv == -2.0 &&
        maxv == 5.0 ? 0.0 : 1.0

    balance_ps = ParticleSet{1,Float64}(4)
    balance_ps.x[1] .= [0.2, 1.2, 2.2, 3.2]
    balance_g = FourierGrid((4,), (4.0,))
    percell = particles_per_cell(balance_ps, balance_g)
    ranges = balanced_tile_ranges(percell, 2)
    balance = particle_load_balance(balance_ps, balance_g; ntiles = 2)
    imbalance = particle_load_imbalance(balance_ps, balance_g; ntiles = 2)
    load_balance_error =
        percell == [1, 1, 1, 1] &&
        tile_loads(percell, 2) == [2, 2] &&
        balanced_tile_loads(percell, ranges) == [2, 2] &&
        balance.per_tile == [2, 2] &&
        abs(imbalance.imbalance - 1.0) <= eps() ? 0.0 : 1.0

    status = GPUAwareMPIStatus(false, false, false, :validation, "synthetic host-buffer status")
    host_buffer = [1.0, 2.0]
    mpi_plan = prepare_mpi_buffer(host_buffer; status, intent = :recv)
    mpi_buffer_error =
        !mpi_buffer_uses_host_staging(host_buffer; status) &&
        host_staging_buffer(host_buffer) === host_buffer &&
        !mpi_plan.used_host_staging &&
        finish_mpi_buffer!(mpi_plan) === host_buffer ? 0.0 : 1.0

    sh3 = PerpShock3D(6, 3, 2, 3.0, 2.0, 1.0)
    ps3 = ParticleSet{3,Float64}(1)
    ps3.x[1][1] = 0.25
    ps3.id[1] = 99
    checkpoint_error = mktempdir() do dir
        path = joinpath(dir, "shock3d.ser")
        checkpoint_shock3d(path, sh3, ps3)
        sh_loaded, ps_loaded = restore_shock3d(path)
        sh_loaded.nx == sh3.nx &&
            sh_loaded.B[3] == sh3.B[3] &&
            ps_loaded.id == ps3.id &&
            ps_loaded.x[1] == ps3.x[1] ? 0.0 : 1.0
    end

    reference_cmp =
        compare_to_reference((a = 1.01, b = 2.0, c = 99.0), (a = 1.0, b = 2.0); rtol = 0.02)
    published_contract =
        published_hybrid_reference_ids() == (:preisser2020_65deg_Bavg_y,) &&
        published_hybrid_reference_metadata().doi == "10.5281/zenodo.3697360" &&
        reference_cmp.pass &&
        reference_cmp.maxrelerr <= 0.02 ? 0.0 : 1.0

    artifact = joinpath(artifact_dir, "11_api_contract_regression_validation.csv")
    rows = (
        ("rank_seed_rng_contract_error", seed_rng_error, 0.0, "absolute", seed_rng_error, 0.0),
        ("global_particle_id_contract_error", id_error, 0.0, "absolute", id_error, 0.0),
        ("sbp_identity_abs_error", sbp_identity_error, 0.0, "absolute", sbp_identity_error, 1e-12),
        (
            "mixed_sbp_fourier_y_max_abs_error",
            mixed_derivative_error,
            0.0,
            "absolute",
            mixed_derivative_error,
            1e-12,
        ),
        ("sbp_deriv_x_linear_max_abs_error", sbp_x_error, 0.0, "absolute", sbp_x_error, 1e-12),
        (
            "spectral_utilities_contract_error",
            spectral_contract_error,
            0.0,
            "absolute",
            spectral_contract_error,
            1e-10,
        ),
        (
            "electron_closure_fields_contract_error",
            fields_closure_error,
            0.0,
            "absolute",
            fields_closure_error,
            1e-12,
        ),
        (
            "gather_vector_constant_field_error",
            gather_vector_error,
            0.0,
            "absolute",
            gather_vector_error,
            1e-12,
        ),
        (
            "current_temperature_contract_error",
            moment_temperature_error,
            0.0,
            "absolute",
            moment_temperature_error,
            1e-12,
        ),
        (
            "rankine_hugoniot_residual_contract_error",
            rh_contract_error,
            0.0,
            "absolute",
            rh_contract_error,
            1e-12,
        ),
        (
            "diagnostic_reduction_contract_error",
            reduction_error,
            0.0,
            "absolute",
            reduction_error,
            0.0,
        ),
        (
            "load_balance_contract_error",
            load_balance_error,
            0.0,
            "absolute",
            load_balance_error,
            0.0,
        ),
        (
            "mpi_host_buffer_contract_error",
            mpi_buffer_error,
            0.0,
            "absolute",
            mpi_buffer_error,
            0.0,
        ),
        (
            "shock3d_checkpoint_roundtrip_error",
            checkpoint_error,
            0.0,
            "absolute",
            checkpoint_error,
            0.0,
        ),
        (
            "reference_comparison_contract_error",
            published_contract,
            0.0,
            "absolute",
            published_contract,
            0.0,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "api_contracts",
        reference_kind = "analytic_or_contract",
        reference = "deterministic API contracts, exact SBP identity, Fourier-y derivative, RH residuals, and reference-comparison invariants",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "11_api_contract_regression_validation",
    default = true,
    description = "Deterministic exported API contracts, FD/spectral identities, RH residuals, and reference comparison.",
    runner = case_11_api_contract_regression_validation,
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
