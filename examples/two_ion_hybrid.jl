# Uniform proton/alpha evolution in a shared magnetic field.
#   julia --project=. examples/two_ion_hybrid.jl
# Velocities are compared at their stored half-time level. The explicit
# predictor has a convergent momentum error, not exact discrete conservation.
using HybridPlasmaPIC, Printf

function two_ion_hybrid()
    g = FourierGrid((32,), (2π,))
    protons = ParticleSet{1,Float64}(32 * 64; q = 1.0, m = 1.0)
    alphas = ParticleSet{1,Float64}(32 * 32; q = 2.0, m = 4.0)
    species = [protons, alphas]
    densities = (0.8, 0.1) # Electron charge density: 1*0.8 + 2*0.1 = 1.
    initial = ((0.0, 0.12, -0.04), (0.0, -0.08, 0.09))
    for (s, ps) in enumerate(species)
        load_lattice_1d!(ps, 0.0, g.L[1])
        set_density_weight!(ps, densities[s], g)
        for c = 1:3
            fill!(ps.v[c], initial[s][c])
        end
    end
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.1)), CIC(), 0)
    fill!(st.fields.B[1], 1.0)
    init!(st, species)
    for _ = 1:200
        step!(st, species, 0.002; NB = 2)
    end

    # For uniform fields, relative velocity rotates at Ω while the mass
    # center stays fixed: Ω = B*(a₂*q₁/m₁ + a₁*q₂/m₂), a_s=q_s*n_s/n_e.
    charge = [ps.q * densities[s] for (s, ps) in enumerate(species)]
    mass = [ps.m * densities[s] for (s, ps) in enumerate(species)]
    frequency = (charge[2] * protons.q / protons.m + charge[1] * alphas.q / alphas.m) / sum(charge)
    center = [sum(mass[s] * initial[s][c] for s = 1:2) / sum(mass) for c = 1:3]
    relative = collect(initial[1]) - collect(initial[2])
    phase = frequency * (st.time[] - st.last_dt[] / 2)
    rotated = [
        relative[1],
        relative[2] * cos(phase) + relative[3] * sin(phase),
        relative[3] * cos(phase) - relative[2] * sin(phase),
    ]
    expected = [center + mass[2] / sum(mass) * rotated, center - mass[1] / sum(mass) * rotated]
    measured = [[sum(ps.v[c]) / nparticles(ps) for c = 1:3] for ps in species]
    velocity_error = maximum(abs(measured[s][c] - expected[s][c]) for s = 1:2, c = 1:3)
    momentum_error =
        maximum(abs(sum(mass[s] * (measured[s][c] - initial[s][c]) for s = 1:2)) for c = 1:3)
    all(isfinite, (velocity_error, momentum_error)) || error("Nonfinite two-ion result")
    # This fixed dt resolves the analytic orbit; the independent convergence
    # tests in test_multispecies_hybrid.jl establish the temporal order.
    velocity_error < 1e-6 || error("Two-ion trajectory disagrees with the analytic solution")
    @printf(
        "Two ions: %d protons, %d alpha particles; electron density %.3f\n",
        nparticles(protons),
        nparticles(alphas),
        sum(charge)
    )
    @printf("t = %.3f, relative gyrofrequency = %.3f\n", st.time[], frequency)
    @printf(
        "maximum velocity error = %.3e; momentum-density error = %.3e\n",
        velocity_error,
        momentum_error
    )
    budget = energy_budget(species, st.fields.B, st.fields.n, st.model.closure, g)
    println("Energy budget at stored time levels: ", budget)
    return st, species
end

two_ion_hybrid()
