#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_02_analytic_hall_mhd_continuity(artifact_dir::AbstractString)
    id = "02_analytic_hall_mhd_continuity"
    n = 64
    l = 2π
    k = 2π / l
    g = FourierGrid((n,), (l,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.0); Ti = 0.0))
    amp = 0.2
    drift = 0.3
    st.fields.n .= @. 1 + amp * cos(k * x)
    fill!(st.fields.ui[1], drift)
    rhs = hall_mhd_rhs!(st)
    expected = @. drift * amp * k * sin(k * x)
    err = maximum(abs, rhs.dn .- expected)

    artifact = joinpath(artifact_dir, "02_analytic_hall_mhd_continuity.csv")
    rows = [(x[i], rhs.dn[i], expected[i]) for i = 1:n]
    _write_csv(artifact, ("x", "measured_dn_dt", "expected_dn_dt"), rows)
    return [
        _result(
            id = id,
            category = "hall_mhd",
            reference_kind = "analytic",
            reference = "-d(nu)/dx for n=1+A cos(kx), u=constant",
            metric = "max_abs_rhs_error",
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
    id = "02_analytic_hall_mhd_continuity",
    default = true,
    description = "Hall-MHD continuity RHS against exact advected density wave.",
    runner = case_02_analytic_hall_mhd_continuity,
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
