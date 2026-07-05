# Energy & momentum BUDGET diagnostics (§6/§29). Oracles:
#  - energy_budget.total drift on a periodic adiabatic hybrid wave is bounded
#    and decreases with dt (HYB-007 protocol: polytropic γ=5/3, η=0, B=x̂ + small
#    By perturbation, cold quiet-start protons, n0=1) — and likewise for the
#    isothermal closure, whose free energy T_e ∫ n ln n dV closes the budget.
#  - the electron-inertia reservoir ∫ d_e²|J|²/2 dV (de2 keyword) conserves the
#    2D Hall+inertia field energy where ∫½|B|² alone wanders by O(fluc energy).
#  - kinetic_energy_relativistic matches Σ w m (γ−1)c² and its Newtonian limit.
#  - momentum_budget.particle/.total are conserved to particle noise, and
#    total ≡ particle (field momentum omitted in the hybrid model).
#  - jdotE_density matches J·E elementwise; its integral matches electric_work.
#  - resistive_dissipation on uniform J,η matches η|J|²·V.

using HybridPlasmaPIC, Test, LinearAlgebra, Random, Statistics

# ---------------------------------------------------------------- helpers

function build_wave_stepper(;
    n = 32,
    L = 2π,
    npc = 300,
    seed = 2,
    clo = PolytropicElectrons(0.5, 1.0, 5 / 3),       # adiabatic, η=0
)
    T = Float64
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
    @test b.electron_inertia == 0.0                    # de2 defaults to 0
    @test b.total == b.kinetic + b.magnetic + b.electron_internal
    @test b.kinetic > 0 && b.magnetic > 0 && b.electron_internal > 0
end

@testset "energy_budget isothermal free energy closes the budget" begin
    # F_e = T_e ∫ n ln n dV is the exact γ→1 limit of the polytropic invariant:
    # same HYB-007 wave protocol as above with the isothermal closure — the
    # budget drift is bounded and decreases with dt (measured 0.099 at dt=0.04,
    # 0.036 at dt=0.01 over 1200 steps).
    function max_drift(dt; nsteps = 1200, NB = 4)
        st, ps, g, clo = build_wave_stepper(clo = IsothermalElectrons(0.5))
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
    @info "isothermal energy_budget drift" d_coarse d_fine
    @test isfinite(d_coarse)
    @test d_coarse < 0.15                  # bounded (no secular blow-up)
    @test d_fine < d_coarse                # drift decreases with Δt
end

@testset "isothermal free energy value: T_e ∫ n ln n dV" begin
    g = FourierGrid((16,), (2π,))
    iso = IsothermalElectrons(0.5)
    @test electron_internal_energy(fill(2.0, 16), iso, g) ≈ 0.5 * 2 * log(2) * 2π
    @test electron_internal_energy(ones(16), iso, g) == 0.0    # n ≡ 1 ⇒ ln n = 0
    n0 = fill(2.0, 16)
    n0[3] = 0.0                                                # n=0 cell contributes 0 (n ln n → 0⁺)
    @test electron_internal_energy(n0, iso, g) ≈ 0.5 * 2 * log(2) * 2π * 15 / 16
    b = energy_budget(ParticleSet{1,Float64}(0), ntuple(_ -> zeros(16), 3), fill(2.0, 16), iso, g)
    @test isfinite(b.electron_internal) && isfinite(b.total)
end

# ---------------------------------------------------------------- electron inertia

@testset "energy_budget electron-inertia reservoir (de2 keyword)" begin
    # 2D field-only Hall+inertia run (u_i = 0, n ≡ 1, η = 0): the inertia-filtered
    # Hall channel conserves ∫(|B|² + d_e²|J|²)/2 dV, while ∫½|B|² alone wanders
    # by O(fluctuation energy) as B exchanges energy with the electron flow
    # (measured: E_B wander 25% of the fluctuation energy, corrected total < 1e-8;
    # the exchange is multi-D — 1D parallel runs are blind to it).
    T = Float64
    nx = 16
    de2 = 0.25
    g = FourierGrid((nx, nx), (2π, 2π))
    # analytic reservoir: B = cos(x) ẑ ⇒ J = ∇×B = sin(x) ŷ ⇒ ∫ de²|J|²/2 dV = de²π²
    ps0 = ParticleSet{2,T}(0)
    n1 = ones(T, nx, nx)
    iso = IsothermalElectrons(0.5)
    x = [(i - 1) * g.dx[1] for i = 1:nx]
    Bt = (zeros(T, nx, nx), zeros(T, nx, nx), [cos(x[i]) for i = 1:nx, j = 1:nx])
    bt = energy_budget(ps0, Bt, n1, iso, g; de2 = de2)
    @test bt.electron_inertia ≈ de2 * π^2
    @test bt.total ≈ energy_budget(ps0, Bt, n1, iso, g).total + de2 * π^2
    @test energy_budget(ps0, Bt, n1, iso, g).electron_inertia == 0.0   # default: no change
    @test_throws ArgumentError energy_budget(ps0, Bt, n1, iso, g; de2 = -0.1)
    @test_throws ArgumentError energy_budget(ps0, Bt, n1, iso, g; de2 = NaN)
    @test_throws ArgumentError energy_budget(ps0, Bt, n1, iso, g; de2 = Inf)
    # --- conservation oracle: RK4 on ∂B/∂t = −∇×E, E = (1+d_e²∇×∇×)⁻¹(J×B), J = ∇×B
    B = ntuple(_ -> zeros(T, nx, nx), 3)
    y = [(j - 1) * g.dx[2] for j = 1:nx]
    fill!(B[3], 1.0)
    modes =
        ((1, (1, 2, 0.3)), (2, (2, 1, 1.1)), (3, (3, 2, 2.0)), (1, (2, 3, 4.0)), (2, (1, 1, 0.7)))
    for (c, (mx, my, ph)) in modes, j = 1:nx, i = 1:nx
        B[c][i, j] += 0.3 * cos(mx * x[i] + my * y[j] + ph)
    end
    project_divfree!(B, g)
    J = ntuple(_ -> zeros(T, nx, nx), 3)
    E = ntuple(_ -> zeros(T, nx, nx), 3)
    function rhs!(dB, Bin)
        curl!(J, Bin, g)
        @inbounds for I in eachindex(E[1])
            jx, jy, jz = J[1][I], J[2][I], J[3][I]
            bx, by, bz = Bin[1][I], Bin[2][I], Bin[3][I]
            E[1][I] = jy * bz - jz * by
            E[2][I] = jz * bx - jx * bz
            E[3][I] = jx * by - jy * bx
        end
        HybridPlasmaPIC._apply_electron_inertia!(E, T(de2), g)
        curl!(dB, E, g)
        for c = 1:3
            @. dB[c] = -dB[c]
        end
        return dB
    end
    k1, k2, k3, k4, Btmp = ntuple(_ -> ntuple(_ -> zeros(T, nx, nx), 3), 5)
    budget(d) = energy_budget(ps0, B, n1, iso, g; de2 = d).total
    EBs = [budget(0.0)]
    TOTs = [budget(de2)]
    dt = 0.02
    for _ = 1:400
        rhs!(k1, B)
        for c = 1:3
            @. Btmp[c] = B[c] + 0.5dt * k1[c]
        end
        rhs!(k2, Btmp)
        for c = 1:3
            @. Btmp[c] = B[c] + 0.5dt * k2[c]
        end
        rhs!(k3, Btmp)
        for c = 1:3
            @. Btmp[c] = B[c] + dt * k3[c]
        end
        rhs!(k4, Btmp)
        for c = 1:3
            @. B[c] += (dt / 6) * (k1[c] + 2k2[c] + 2k3[c] + k4[c])
        end
        push!(EBs, budget(0.0))
        push!(TOTs, budget(de2))
    end
    Efluc = EBs[1] - 0.5 * prod(g.L)                   # subtract the uniform-B0 energy
    wEB = maximum(EBs) - minimum(EBs)
    wTOT = maximum(TOTs) - minimum(TOTs)
    @info "electron-inertia budget wander" Efluc wEB wTOT
    @test wEB > 0.05 * Efluc                # E_B alone is NOT conserved (measured 25%)
    @test wTOT < 1e-3 * wEB                 # corrected total is (measured < 1e-8 abs)
end

# ---------------------------------------------------------------- relativistic KE

@testset "kinetic_energy_relativistic: value, Newtonian limit, guards" begin
    T = Float64
    c = 5.0
    ps = ParticleSet{1,T}(1)
    ps.x[1] .= 0.0
    ps.v[1] .= 0.9 * c                                 # γ = 2.294…
    γ = 1 / sqrt(1 - 0.81)
    @test kinetic_energy_relativistic(ps, c) ≈ (γ - 1) * c^2 rtol = 1e-12   # 32.3539…
    @test kinetic_energy(ps) ≈ 10.125                  # Newtonian undercounts 3.2×
    # weights and mass: Σ_p w m (γ_p − 1) c², any velocity component
    ps2 = ParticleSet{1,T}(2; m = 2.0)
    ps2.x[1] .= 0.0
    ps2.weight .= [1.0, 3.0]
    ps2.v[1] .= [0.9 * c, 0.0]
    ps2.v[2] .= [0.0, 0.6 * c]                         # γ = 1.25
    @test kinetic_energy_relativistic(ps2, c) ≈ 2.0 * (1.0 * (γ - 1) + 3.0 * 0.25) * c^2 rtol =
        1e-12
    @test kinetic_energy_relativistic(ParticleSet{1,T}(0), c) == 0.0
    # Newtonian limit: rel − newt = (3/8) w m v⁴/c² (1 + O(β²)) — quartic in v
    dEs = map((0.1, 0.05)) do β
        ps.v[1] .= β * c
        d = kinetic_energy_relativistic(ps, c) - kinetic_energy(ps)
        @test d ≈ (3 / 8) * (β * c)^4 / c^2 rtol = 0.02
        d
    end
    @test dEs[1] / dEs[2] ≈ 16 rtol = 0.02             # v⁴ scaling, not v²
    # guards: |v| ≥ c is corrupted state; c must be finite and positive
    ps.v[1] .= c
    @test_throws ArgumentError kinetic_energy_relativistic(ps, c)
    ps.v[1] .= 1.1 * c
    @test_throws ArgumentError kinetic_energy_relativistic(ps, c)
    ps.v[1] .= 0.5
    @test_throws ArgumentError kinetic_energy_relativistic(ps, 0.0)
    @test_throws ArgumentError kinetic_energy_relativistic(ps, -1.0)
    @test_throws ArgumentError kinetic_energy_relativistic(ps, NaN)
    @test_throws ArgumentError kinetic_energy_relativistic(ps, Inf)
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
