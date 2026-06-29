#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_01_analytic_poisson_2d(artifact_dir::AbstractString)
    id = "01_analytic_poisson_2d"
    nx, ny = 32, 24
    lx, ly = 2π, 2π
    g = FourierGrid((nx, ny), (lx, ly))
    es = ElectrostaticPIC(g, 0; n0 = 1.25)
    amp = 0.2
    mx, my = 2, 3
    kx = 2π * mx / lx
    ky = 2π * my / ly
    k2 = kx^2 + ky^2
    ex = zeros(Float64, nx, ny)
    ey = zeros(Float64, nx, ny)
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        y = (j - 1) * g.dx[2]
        rho = amp * cos(kx * x) * cos(ky * y)
        es.ne[i, j] = es.n0 - rho
        ex[i, j] = amp * kx / k2 * sin(kx * x) * cos(ky * y)
        ey[i, j] = amp * ky / k2 * cos(kx * x) * sin(ky * y)
    end
    poisson_E!(es)
    err = maximum((maximum(abs, es.E[1] .- ex), maximum(abs, es.E[2] .- ey), maximum(abs, es.E[3])))

    artifact = joinpath(artifact_dir, "01_analytic_poisson_2d_slice.csv")
    rows = [((i - 1) * g.dx[1], es.E[1][i, 1], ex[i, 1], es.E[2][i, 1], ey[i, 1]) for i = 1:nx]
    _write_csv(artifact, ("x", "measured_Ex", "expected_Ex", "measured_Ey", "expected_Ey"), rows)
    return [
        _result(
            id = id,
            category = "electrostatic_pic",
            reference_kind = "analytic",
            reference = "Fourier Poisson single-mode electric field",
            metric = "max_abs_field_error",
            measured = err,
            expected = 0.0,
            error_kind = "absolute",
            error = err,
            tolerance = 1e-11,
            artifact = basename(artifact),
        ),
    ]
end


VALIDATION_CASE = ValidationCase(
    id = "01_analytic_poisson_2d",
    default = true,
    description = "2D spectral Poisson solve against exact Fourier-mode electric field.",
    runner = case_01_analytic_poisson_2d,
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
