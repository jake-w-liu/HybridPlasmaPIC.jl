#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_07_closure_budget_filter_diagnostics(artifact_dir::AbstractString)
    id = "07_closure_budget_filter_diagnostics"

    n = 32
    g = FourierGrid((n,), (2π,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    pe = fill(2.0, n)
    ue = (0.3 .* sin.(x), zeros(Float64, n), zeros(Float64, n))
    expected_pe = @. 2.0 + 0.05 * (-(5 / 3) * 2.0 * 0.3 * cos(x))
    HybridPlasmaPIC.advance_electron_pressure!(pe, ue, 0.05, 5 / 3, g)
    pressure_evolution_error = maximum(abs, pe .- expected_pe)

    ps = ParticleSet{1,Float64}(2; m = 2.0)
    ps.weight .= [1.0, 2.0]
    ps.v[1] .= [1.0, -2.0]
    ps.v[2] .= [0.5, 1.0]
    ps.v[3] .= [0.0, -0.5]
    B = (fill(0.2, n), fill(-0.1, n), fill(0.3, n))
    density = fill(1.5, n)
    closure = PolytropicElectrons(0.75, 1.5, 5 / 3)
    budget = energy_budget(ps, B, density, closure, g)
    kinetic_ref =
        sum(ps.weight[p] * 0.5 * ps.m * (ps.v[1][p]^2 + ps.v[2][p]^2 + ps.v[3][p]^2) for p = 1:2)
    magnetic_ref = 0.5 * n * (0.2^2 + (-0.1)^2 + 0.3^2) * prod(g.dx)
    electron_ref = (0.75 / (5 / 3 - 1)) * prod(g.L)
    budget_error = maximum(
        abs,
        (
            budget.kinetic - kinetic_ref,
            budget.magnetic - magnetic_ref,
            budget.electron_internal - electron_ref,
        ),
    )
    isothermal_budget = energy_budget(ps, B, density, IsothermalElectrons(1.0), g)
    isothermal_error =
        isnan(isothermal_budget.electron_internal) && isnan(isothermal_budget.total) ? 0.0 : 1.0

    J = (fill(2.0, n), fill(0.0, n), fill(-1.0, n))
    E = (fill(3.0, n), fill(5.0, n), fill(4.0, n))
    jE = jdotE_density(J, E)
    jdot_error = max(maximum(abs, jE .- 2.0), abs(sum(jE) * prod(g.dx) - electric_work(J, E, g)))
    diss_error = abs(resistive_dissipation(J, 0.3, g) - 0.3 * (2.0^2 + (-1.0)^2) * prod(g.L))

    P = (fill(4.0, n), fill(2.0, n), fill(2.0, n), zeros(n), zeros(n), zeros(n))
    Tpar, Tperp = temperatures_par_perp(P, fill(2.0, n), (ones(n), zeros(n), zeros(n)))
    temp_error = max(maximum(abs, Tpar .- 2.0), maximum(abs, Tperp .- 1.0))

    strain = pressure_strain(P, (sin.(x), zeros(n), zeros(n)), g)
    strain_error = maximum(abs, strain .- (-4.0 .* cos.(x)))
    k, power = power_spectrum(cos.(3 .* x), g)
    spectrum_error = abs((argmax(power) - 1) - 3) + abs(k[4] - 3.0)

    artifact = joinpath(artifact_dir, "07_closure_budget_filter_diagnostics.csv")
    rows = (
        (
            "electron_pressure_evolution_max_abs_error",
            pressure_evolution_error,
            0.0,
            "absolute",
            pressure_evolution_error,
            1e-12,
        ),
        (
            "energy_budget_component_max_abs_error",
            budget.total,
            kinetic_ref + magnetic_ref + electron_ref,
            "absolute",
            budget_error,
            1e-12,
        ),
        (
            "isothermal_budget_nan_contract_error",
            isothermal_error,
            0.0,
            "absolute",
            isothermal_error,
            0.0,
        ),
        (
            "jdotE_density_integral_max_abs_error",
            sum(jE) * prod(g.dx),
            electric_work(J, E, g),
            "absolute",
            jdot_error,
            1e-12,
        ),
        (
            "resistive_dissipation_abs_error",
            resistive_dissipation(J, 0.3, g),
            0.3 * 5.0 * prod(g.L),
            "absolute",
            diss_error,
            1e-12,
        ),
        (
            "temperature_parallel_perp_max_abs_error",
            maximum(Tpar),
            2.0,
            "absolute",
            temp_error,
            1e-12,
        ),
        (
            "pressure_strain_single_mode_max_abs_error",
            maximum(abs, strain),
            4.0,
            "absolute",
            strain_error,
            1e-12,
        ),
        (
            "power_spectrum_peak_index_error",
            argmax(power) - 1,
            3.0,
            "absolute",
            spectrum_error,
            1e-12,
        ),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "closures_budgets_diagnostics",
        reference_kind = "analytic",
        reference = "adiabatic pressure update, exact energy/work/heating budgets, pressure tensor projections, and spectral peak",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "07_closure_budget_filter_diagnostics",
    default = true,
    description = "Electron pressure, energy/work budgets, pressure diagnostics, and spectra.",
    runner = case_07_closure_budget_filter_diagnostics,
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
