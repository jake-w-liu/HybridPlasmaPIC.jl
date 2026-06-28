#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_03_analytic_spectral_operators(artifact_dir::AbstractString)
    id = "03_analytic_spectral_operators"
    n = (16, 12, 10)
    modes = (2, 3, 1)
    g = FourierGrid(n, ntuple(_ -> 2π, 3))
    f = Array{Float64,3}(undef, n)
    dfdx = similar(f)
    lap = similar(f)
    k2 = sum(abs2, modes)
    for index in CartesianIndices(f)
        i, j, k = Tuple(index)
        x = (i - 1) * g.dx[1]
        y = (j - 1) * g.dx[2]
        z = (k - 1) * g.dx[3]
        base = sin(modes[1] * x) * cos(modes[2] * y) * cos(modes[3] * z)
        f[index] = base
        dfdx[index] = modes[1] * cos(modes[1] * x) * cos(modes[2] * y) * cos(modes[3] * z)
        lap[index] = -k2 * base
    end
    dfdx_measured = similar(f)
    deriv!(dfdx_measured, f, g, 1)
    derivative_relerr = norm(dfdx_measured .- dfdx) / norm(dfdx)

    lap_measured = similar(f)
    laplacian!(lap_measured, f, g)
    laplacian_relerr = norm(lap_measured .- lap) / norm(lap)

    rng = MersenneTwister(17)
    A = ntuple(_ -> randn(rng, n...), 3)
    B = ntuple(_ -> zeros(Float64, n), 3)
    curl!(B, A, g)
    divB = zeros(Float64, n)
    divergence!(divB, B, g)
    bnorm = sum(norm, B)
    kmax = maximum(maximum(abs, g.kvec[d]) for d = 1:3)
    divcurl_residual = norm(divB) / (kmax * bnorm + eps(Float64))

    target = ntuple(c -> copy(B[c]), 3)
    phi = randn(rng, n...)
    gradphi = ntuple(_ -> zeros(Float64, n), 3)
    gradient!(gradphi, phi, g)
    projected = ntuple(c -> target[c] .+ gradphi[c], 3)
    project_divfree!(projected, g)
    divergence!(divB, projected, g)
    projection_div_residual = norm(divB) / (kmax * sum(norm, projected) + eps(Float64))
    projection_recovery_relerr =
        sum(norm(projected[c] .- target[c]) for c = 1:3) / (sum(norm, target) + eps(Float64))

    artifact = joinpath(artifact_dir, "03_analytic_spectral_operators.csv")
    rows = (
        ("derivative_relative_l2_error", derivative_relerr, 0.0, "relative", derivative_relerr, 1e-10),
        ("laplacian_relative_l2_error", laplacian_relerr, 0.0, "relative", laplacian_relerr, 1e-10),
        ("divergence_of_curl_residual", divcurl_residual, 0.0, "relative", divcurl_residual, 1e-10),
        ("projection_divergence_residual", projection_div_residual, 0.0, "relative", projection_div_residual, 1e-10),
        ("projection_recovery_relative_error", projection_recovery_relerr, 0.0, "relative", projection_recovery_relerr, 1e-10),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "spectral_operators",
        reference_kind = "analytic",
        reference = "Fourier derivative/laplacian identities and vector-calculus invariants",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "03_analytic_spectral_operators",
    default = true,
    description = "3D Fourier derivative/laplacian/vector identities against analytic references.",
    runner = case_03_analytic_spectral_operators,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
