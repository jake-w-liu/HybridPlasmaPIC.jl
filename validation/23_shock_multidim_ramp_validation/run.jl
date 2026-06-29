#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_23_shock_multidim_ramp_validation(artifact_dir::AbstractString)
    id = "23_shock_multidim_ramp_validation"

    sh2 = PerpShock2D(16, 32, 10.0, 8.0; B0 = 1.0)
    mseed = 3
    base = 8
    amp_nodes = 2
    expected_xs = zeros(Float64, sh2.ny)
    for j = 1:sh2.ny
        ip = base + round(Int, amp_nodes * cos(2π * mseed * (j - 1) / sh2.ny))
        expected_xs[j] = sh2.x[ip+1]
        for i = 1:sh2.nx
            sh2.Bz[i, j] = i <= ip ? 2.0 : 1.0
        end
    end
    xs2, mean2, sigma2 = shock_surface(sh2)
    spec = shock_surface_spectrum(sh2)
    pnon = copy(spec.Ps)
    pnon[1] = -Inf
    surface2_error = maximum(abs, xs2 .- expected_xs)
    spectrum_peak_error = abs((argmax(pnon) - 1) - mseed)
    coh = transverse_coherence(sh2)
    expected_fluct = expected_xs .- sum(expected_xs) / length(expected_xs)
    expected_var = sum(abs2, expected_fluct) / length(expected_fluct)
    expected_coh = zeros(Float64, sh2.ny)
    for lag = 0:sh2.ny-1
        acc = 0.0
        for j = 1:sh2.ny
            jj = mod(j - 1 + lag, sh2.ny) + 1
            acc += expected_fluct[j] * expected_fluct[jj]
        end
        expected_coh[lag+1] = (acc / sh2.ny) / expected_var
    end
    coherence_error = max(
        maximum(abs, coh.C .- expected_coh),
        maximum(abs, coh.dy .- [k * sh2.dy for k = 0:sh2.ny-1]),
    )

    sh3 = PerpShock3D(12, 8, 6, 6.0, 4.0, 3.0; B0 = 1.0)
    expected3 = zeros(Float64, sh3.ny, sh3.nz)
    for k3 = 1:sh3.nz, j = 1:sh3.ny
        ip = 5 + ((j + k3) % 2)
        expected3[j, k3] = sh3.x[ip+1]
        for i = 1:sh3.nx
            sh3.B[3][i, j, k3] = i <= ip ? 2.0 : 1.0
        end
    end
    xs3, mean3, sigma3 = shock_surface3d(sh3)
    surface3_error = maximum(abs, xs3 .- expected3)
    div3_error = maximum(abs, magnetic_divergence3d(PerpShock3D(8, 4, 4, 2.0, 2.0, 2.0)))

    N = 64
    Lx = 20.0
    xramp = 10.0
    width = 2.0
    ramp = PerpShock(N, Lx; B0 = 1.0, τ = 0.0)
    particles = ParticleSet{1,Float64}(4 * N)
    initial_ramp!(ramp, particles, xramp, width, 1.0, 3.0, 1.0, 3.0; rng = MersenneTwister(12))
    expected_B = [1.0 + (3.0 - 1.0) * 0.5 * (1 - tanh((xi - xramp) / width)) for xi in ramp.x]
    ramp_field_error = maximum(abs, ramp.Bz .- expected_B)

    artifact = joinpath(artifact_dir, "23_shock_multidim_ramp_validation.csv")
    rows = (
        (
            "shock2d_surface_position_max_abs_error",
            surface2_error,
            0.0,
            "absolute",
            surface2_error,
            1e-12,
        ),
        (
            "shock2d_surface_spectrum_peak_error",
            argmax(pnon) - 1,
            mseed,
            "absolute",
            spectrum_peak_error,
            0.0,
        ),
        (
            "shock2d_transverse_coherence_max_abs_error",
            coherence_error,
            0.0,
            "absolute",
            coherence_error,
            1e-12,
        ),
        (
            "shock3d_surface_position_max_abs_error",
            surface3_error,
            0.0,
            "absolute",
            surface3_error,
            1e-12,
        ),
        (
            "shock3d_uniform_divergence_max_abs_error",
            div3_error,
            0.0,
            "absolute",
            div3_error,
            1e-12,
        ),
        (
            "initial_ramp_field_max_abs_error",
            ramp_field_error,
            0.0,
            "absolute",
            ramp_field_error,
            1e-12,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "shock_multidim_initial_conditions",
        reference_kind = "analytic",
        reference = "constructed 2D/3D shock surfaces, uniform div-B identity, and tanh ramp field",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "23_shock_multidim_ramp_validation",
    default = true,
    description = "Constructed 2D/3D shock surfaces, 3D div-B, and initial tanh ramp.",
    runner = case_23_shock_multidim_ramp_validation,
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
