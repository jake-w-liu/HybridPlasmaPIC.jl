#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_18_hybrid_ion_acoustic_dispersion(artifact_dir::AbstractString)
    id = "18_hybrid_ion_acoustic_dispersion"
    te = 1.0
    m = 1
    l = 2π
    k = 2π * m / l
    g, ps, st = _setup_hybrid_1d(64, l, 1; nppc = 600, te = te)
    for p = 1:nparticles(ps)
        ps.v[1][p] += 0.005 * sin(k * ps.x[1][p])
    end
    init!(st, ps)
    dt = 0.02
    nt = 700
    series = Float64[]
    for _ = 1:nt
        step!(st, ps, dt)
        push!(series, real(mode_amplitude(st.fields.n, g, (m,))))
    end
    omega = _zero_cross_frequency(series, dt)
    omega_ref = k * sqrt(te)
    relerr = abs(omega - omega_ref) / omega_ref

    artifact = joinpath(artifact_dir, "18_hybrid_ion_acoustic_dispersion.csv")
    rows = [(i * dt, series[i]) for i = 1:length(series)]
    _write_csv(artifact, ("time", "density_mode_real"), rows)
    return [
        _result(
            id = id,
            category = "hybrid_pic",
            reference_kind = "analytic",
            reference = "ion-acoustic dispersion omega = k sqrt(Te)",
            metric = "relative_frequency_error",
            measured = omega,
            expected = omega_ref,
            error_kind = "relative",
            error = relerr,
            tolerance = 0.03,
            artifact = basename(artifact),
        ),
    ]
end


VALIDATION_CASE = ValidationCase(
    id = "18_hybrid_ion_acoustic_dispersion",
    default = false,
    description = "Hybrid PIC ion-acoustic frequency against omega = k sqrt(Te).",
    runner = case_18_hybrid_ion_acoustic_dispersion,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
