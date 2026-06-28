#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_20_empic_transverse_dispersion_2d(artifact_dir::AbstractString)
    id = "20_empic_transverse_dispersion_2d"
    nx, ny = 32, 4
    l = 2π
    c = 5.0
    n0 = 1.0
    g = FourierGrid((nx, ny), (l, l))
    nppc = 80
    np = nppc * nx * ny
    electrons = ParticleSet{2,Float64}(np; q = -1.0, m = 1.0)
    load_lattice!(electrons, (0.0, 0.0), g.L, (nx * nppc, ny))
    set_density_weight!(electrons, n0, g)
    em = EMPIC(g, np; n0 = n0, c = c, shape = CIC())
    m = 1
    k = 2π * m / l
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        em.E[2][i, j] = 1e-3 * cos(k * x)
    end
    init_empic!(em, electrons)
    dt = 0.01
    nt = 1000
    series = ComplexF64[]
    for _ = 1:nt
        step_empic!(em, electrons, dt)
        push!(series, mode_amplitude(em.E[2], g, (m, 0)))
    end
    omega = _peak_frequency(series, dt)
    omega_ref = sqrt(n0 + c^2 * k^2)
    relerr = abs(omega - omega_ref) / omega_ref
    charge_residual = charge_conservation_residual(em, dt)

    artifact = joinpath(artifact_dir, "20_empic_transverse_dispersion_2d.csv")
    rows = [(i * dt, real(series[i]), imag(series[i]), abs(series[i])) for i = 1:length(series)]
    _write_csv(artifact, ("time", "mode_real", "mode_imag", "mode_abs"), rows)
    return [
        _result(
            id = id,
            category = "em_pic",
            reference_kind = "analytic",
            reference = "cold transverse EM dispersion omega^2 = omega_pe^2 + c^2 k^2",
            metric = "relative_frequency_error",
            measured = omega,
            expected = omega_ref,
            error_kind = "relative",
            error = relerr,
            tolerance = 0.03,
            artifact = basename(artifact),
        ),
        _result(
            id = id,
            category = "em_pic",
            reference_kind = "analytic",
            reference = "discrete charge conservation residual",
            metric = "charge_conservation_residual",
            measured = charge_residual,
            expected = 0.0,
            error_kind = "absolute",
            error = charge_residual,
            tolerance = 1e-8,
            artifact = basename(artifact),
        ),
    ]
end


VALIDATION_CASE = ValidationCase(
    id = "20_empic_transverse_dispersion_2d",
    default = false,
    description = "2D EM PIC transverse wave frequency against cold EM dispersion.",
    runner = case_20_empic_transverse_dispersion_2d,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
