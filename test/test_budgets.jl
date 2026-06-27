# Energy & momentum BUDGET diagnostics (§6/§29). Oracles:
#  - energy_budget.total drift on a periodic adiabatic hybrid wave is bounded
#    and decreases with dt (HYB-007 protocol: polytropic γ=5/3, η=0, B=x̂ + small
#    By perturbation, cold quiet-start protons, n0=1).
#  - momentum_budget.particle/.total are conserved to particle noise, and
#    total ≡ particle (field momentum omitted in the hybrid model).
#  - jdotE_density matches J·E elementwise; its integral matches electric_work.
#  - resistive_dissipation on uniform J,η matches η|J|²·V.

using HybridPlasmaPIC, Test, LinearAlgebra, Random, Statistics

# ---------------------------------------------------------------- helpers

function build_wave_stepper(; n = 32, L = 2π, npc = 300, seed = 2)
    T = Float64
    clo = PolytropicElectrons(0.5, 1.0, 5 / 3)        # adiabatic, η=0
    k = 2π / L
    g = FourierGrid((n,), (L,))
    N = npc * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, L)
    set_density_weight!(ps, 1.0, g)                   # n0 = 1
    load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(clo), CIC(), N)
    fill!(st.fields.B[1], 1.0)                        # B0 = x̂
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st.fields.B[2] .= 0.05 .* cos.(k .* x)            # small By perturbation
    init!(st, ps)
    return st, ps, g, clo
end

# ---------------------------------------------------------------- energy budget

@testset "energy_budget total drift bounded & decreases with dt" begin
    function max_drift(dt; nsteps = 1200, NB = 4)
        st, ps, g, clo = build_wave_stepper()
        Etot() = energy_budget(ps, st.fields.B, st.fields.n, clo, g).total
        E0 = Etot()
        d = 0.0
        for _ = 1:nsteps
            step!(st, ps, dt; NB)
            d = max(d, abs(Etot() - E0) / E0)
        end
        return d
    end
    d_coarse = max_drift(0.04)
    d_fine = max_drift(0.01)
    @info "energy_budget drift" d_coarse d_fine
    @test isfinite(d_coarse)
    @test d_coarse < 0.15                  # bounded (no secular blow-up)
    @test d_fine < d_coarse                # drift decreases with Δt
end

@testset "energy_budget component agreement with primitives" begin
    st, ps, g, clo = build_wave_stepper()
    b = energy_budget(ps, st.fields.B, st.fields.n, clo, g)
    @test b.kinetic == kinetic_energy(ps)
    @test b.magnetic == magnetic_energy(st.fields.B, g)
    @test b.electron_internal == electron_internal_energy(st.fields.n, clo, g)
    @test b.total == b.kinetic + b.magnetic + b.electron_internal
    @test b.kinetic > 0 && b.magnetic > 0 && b.electron_internal > 0
end

@testset "energy_budget isothermal electron_internal is NaN" begin
    st, ps, g, _ = build_wave_stepper()
    iso = IsothermalElectrons(0.5)
    b = energy_budget(ps, st.fields.B, st.fields.n, iso, g)
    @test isnan(b.electron_internal)
    @test isnan(b.total)
    @test isfinite(b.kinetic) && isfinite(b.magnetic)
end

# ---------------------------------------------------------------- momentum budget

@testset "momentum_budget conserved to particle noise; total ≡ particle" begin
    st, ps, g, _ = build_wave_stepper()
    mb0 = momentum_budget(ps, st.fields.B, g)
    # total ≡ particle (field momentum omitted)
    @test mb0.total == mb0.particle
    @test mb0.particle == total_momentum(ps)
    # scale for the noise tolerance: mean |momentum-per-particle contribution|
    pscale = ps.m * sum(abs, ps.weight) * 0.1   # ~ m * Σw * v_th
    dt = 0.02
    for _ = 1:600
        step!(st, ps, dt; NB = 4)
    end
    mb1 = momentum_budget(ps, st.fields.B, g)
    @test mb1.total == mb1.particle
    for c = 1:3
        drift = abs(mb1.particle[c] - mb0.particle[c])
        @info "momentum drift" c drift pscale
        @test drift < 0.05 * pscale          # conserved to particle noise
    end
end

# ---------------------------------------------------------------- jdotE_density

@testset "jdotE_density matches J·E elementwise and integral" begin
    T = Float64
    n = 24
    L = 2π
    g = FourierGrid((n,), (L,))
    rng = MersenneTwister(7)
    J = ntuple(_ -> randn(rng, n), 3)
    E = ntuple(_ -> randn(rng, n), 3)
    jE = jdotE_density(J, E)
    @test size(jE) == size(J[1])
    for I in eachindex(jE)
        @test jE[I] ≈ J[1][I] * E[1][I] + J[2][I] * E[2][I] + J[3][I] * E[3][I]
    end
    @test sum(jE) * prod(g.dx) ≈ electric_work(J, E, g)
end

@testset "jdotE_density on a known J,E" begin
    T = Float64
    n = 8
    J = (fill(2.0, n), fill(0.0, n), fill(-1.0, n))
    E = (fill(3.0, n), fill(5.0, n), fill(4.0, n))
    jE = jdotE_density(J, E)
    @test all(≈(2.0 * 3.0 + 0.0 * 5.0 + (-1.0) * 4.0), jE)   # = 2.0
end

# ---------------------------------------------------------------- resistive heating

@testset "resistive_dissipation on uniform J,η == η|J|²·V" begin
    T = Float64
    for D = 1:3
        nc = ntuple(_ -> 12, D)
        g = FourierGrid(nc, ntuple(_ -> 2π, D))
        Jc = (1.5, -0.5, 2.0)
        J = ntuple(c -> fill(Jc[c], nc), 3)
        η = 0.3
        V = prod(g.L)                               # full periodic volume
        absJ2 = Jc[1]^2 + Jc[2]^2 + Jc[3]^2
        @test resistive_dissipation(J, η, g) ≈ η * absJ2 * V
    end
end

@testset "resistive_dissipation η=0 gives zero" begin
    T = Float64
    n = 16
    g = FourierGrid((n,), (2π,))
    J = ntuple(_ -> randn(MersenneTwister(3), n), 3)
    @test resistive_dissipation(J, 0.0, g) == 0.0
    @test_throws ArgumentError resistive_dissipation(J, NaN, g)
    @test_throws ArgumentError resistive_dissipation(J, Inf, g)
    @test_throws ArgumentError resistive_dissipation(J, -0.1, g)
end
