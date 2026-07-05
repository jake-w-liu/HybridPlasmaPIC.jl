# §6.4 electron velocity u_e = u_i − J/n, and §10.3 whistler-CFL recommended_dt.

using HybridPlasmaPIC, Test, Random

@testset "§6.4 electron_velocity! u_e = u_i − J/n (array form, exact)" begin
    T = Float64
    nc = (6,)
    ui = ntuple(c -> T[c + 0.1i for i = 1:6], 3)
    J = ntuple(c -> T[0.5c - 0.2i for i = 1:6], 3)
    n = T[1.0, 2.0, 0.5, 4.0, 0.25, 8.0]
    ue = ntuple(_ -> zeros(T, nc), 3)
    electron_velocity!(ue, ui, J, n; nfloor = 1e-6)
    for c = 1:3, i = 1:6
        @test ue[c][i] ≈ ui[c][i] - J[c][i] / max(n[i], 1e-6)
    end
    # density floor: a below-floor cell uses the floor, not the raw (tiny) density
    n2 = copy(n)
    n2[3] = 1e-12
    electron_velocity!(ue, ui, J, n2; nfloor = 1e-6)
    @test ue[1][3] ≈ ui[1][3] - J[1][3] / 1e-6
    @test_throws ArgumentError electron_velocity!(ue, ui, J, n; nfloor = 0.0)
    # mismatched component sizes are rejected, not read out of bounds
    Jbad = (J[1], T[1.0, 2.0], J[3])             # J[2] too short
    @test_throws DimensionMismatch electron_velocity!(ue, ui, Jbad, n)
    uibad = (ui[1], ui[2], T[1.0])               # ui[3] too short
    @test_throws DimensionMismatch electron_velocity!(ue, uibad, J, n)
end

@testset "§6.4 electron_velocity! from HybridFields (analytic curl)" begin
    T = Float64
    n, L = 16, 2π
    g = FourierGrid((n,), (T(L),))
    k = 2π / L
    x = [(i - 1) * g.dx[1] for i = 1:n]
    f = HybridFields{1,T}(g.n)
    fill!(f.n, 2.0)                                   # n0 = 2
    f.B[2] .= sin.(k .* x)                            # B = (0, sin kx, 0)
    for c = 1:3
        fill!(f.ui[c], 0.3 * c)                       # uniform ion velocity
    end
    ue = ntuple(_ -> zeros(T, g.n), 3)
    electron_velocity!(ue, f, g)
    # J = ∇×B = (0, 0, ∂x B_y) = (0, 0, k cos kx); u_e = u_i − J/n
    for i = 1:n
        @test ue[1][i] ≈ 0.3 * 1
        @test ue[2][i] ≈ 0.3 * 2
        @test isapprox(ue[3][i], 0.3 * 3 - (k * cos(k * x[i])) / 2.0; atol = 1e-9)
    end
end

@testset "§10.3 recommended_dt whistler CFL" begin
    T = Float64
    # First-principles oracle, written independently of the implementation: the
    # frozen-moment B-subcycle is electron-MHD with whistler ω = (k·b̂)|k| ≤ |k|²
    # (Hall branch ω_W below, d_i = Ω_ci = 1), stiffest at the k-space CORNER
    # |k|² = Σ_d (π/Δx_d)² — the largest |k| the grid represents, reached when B0
    # is along the corner direction. In 1-D the corner IS the axis Nyquist π/Δx.
    # C are the classical imaginary-axis mover limits (RK4 2·√2 ≈ 2.8; leapfrog
    # family 2.0), not values read from the implementation.
    ωw(K) = 0.5 * (sqrt(K^4 + 4K^2) + K^2)
    g = FourierGrid((64,), (T(2π),))
    kmax = π / g.dx[1]                                 # 1-D: corner = axis Nyquist
    # whistler-limited (fine grid): ω_W(k_max)·Δt_B = ω_W·(dt/NB) ≈ safety·C
    for (integ, C) in ((:rk4, 2.8), (:leapfrog, 2.0), (:camcl, 2.0))
        dt = recommended_dt(g; NB = 1, integrator = integ, safety = 0.8)
        @test dt > 0
        @test isapprox(ωw(kmax) * dt, 0.8 * C; rtol = 1e-12)   # binding whistler limit
    end
    # 2-D anisotropic grid: the binding mode is the corner √((π/Δx)²+(π/Δy)²),
    # strictly stiffer than the single-axis Nyquist π/min(Δx) that a 1-D-minded
    # formula would use (a factor 2 in ω on a square grid, 3 in 3-D).
    g2d = FourierGrid((32, 48), (T(2π), T(2π)))
    kcorner = sqrt((π / g2d.dx[1])^2 + (π / g2d.dx[2])^2)
    for (integ, C) in ((:rk4, 2.8), (:camcl, 2.0))
        dt = recommended_dt(g2d; NB = 1, integrator = integ, safety = 0.8)
        @test isapprox(ωw(kcorner) * dt, 0.8 * C; rtol = 1e-12)
        @test dt < 0.8 * C / ωw(π / minimum(g2d.dx))   # tighter than the axis-Nyquist formula
    end
    # scaling: halving Δx (doubling n) shrinks dt as ω_W(k_max1)/ω_W(k_max2) ≈ 1/4
    g2 = FourierGrid((128,), (T(2π),))
    dt1 = recommended_dt(g; integrator = :rk4)
    dt2 = recommended_dt(g2; integrator = :rk4)
    @test dt2 < dt1
    @test isapprox(dt2 / dt1, ωw(π / g.dx[1]) / ωw(π / g2.dx[1]); rtol = 1e-12)
    @test dt2 / dt1 < 0.3                              # ≈1/4 in the asymptotic k² regime
    # subcycling relaxes the whistler limit linearly until the gyro cap
    @test recommended_dt(g; NB = 4, integrator = :rk4) ≈
          4 * recommended_dt(g; NB = 1, integrator = :rk4)
    # Crank–Nicolson is unconditional on the whistler ⇒ only the gyro limit applies
    @test recommended_dt(g; integrator = :cn) ≈ 0.8 * 0.3
    @test recommended_dt(g; integrator = :cn) > recommended_dt(g; integrator = :rk4)
    # resistive/hyper-resistive real-axis cap: the subcycle integrates the decay
    # rate ν(k) = η k² + ηH k⁴ explicitly, so dt must also satisfy
    # ν(k_max)·Δt_B ≤ safety·C_real with RK4's classical negative-real-axis
    # interval 2.785 (|1+z+z²/2+z³/6+z⁴/24| = 1 at z ≈ −2.785).
    η0, ηH0 = 1e-3, 1e-2                               # grid-scale hyper-resistivity: ν(k_max)·Δt ≫ C_real at the whistler-only dt
    ν(k) = η0 * k^2 + ηH0 * k^4
    dtR = recommended_dt(g; NB = 1, integrator = :rk4, η = η0, ηH = ηH0)
    @test isapprox(dtR, min(0.8 * 2.8 / ωw(kmax), 0.8 * 2.785 / ν(kmax)); rtol = 1e-12)
    @test dtR < recommended_dt(g; integrator = :rk4)   # ηH k⁴ binds on this grid
    # η = ηH = 0 is exactly the undamped recommendation (backward compatible)
    @test recommended_dt(g; integrator = :rk4, η = 0.0, ηH = 0.0) ==
          recommended_dt(g; integrator = :rk4)
    # CN is A-stable: real-axis decay never restricts it (gyro cap only)
    @test recommended_dt(g; integrator = :cn, η = 1.0, ηH = 1.0) ≈ 0.8 * 0.3
    # the leapfrog/CL family tolerates essentially no real-axis damping: with
    # η/ηH > 0 both a small ν·Δt_B budget and a reduced whistler constant apply,
    # so the recommendation must be far below the undamped one
    dtCL = recommended_dt(g; NB = 1, integrator = :camcl, η = η0, ηH = ηH0)
    @test dtCL < 0.2 * recommended_dt(g; NB = 1, integrator = :camcl)
    @test dtCL <= 0.8 * min(0.3 / ωw(kmax), 0.25 / ν(kmax)) * (1 + 1e-12)
    # validation
    @test_throws ArgumentError recommended_dt(g; NB = 0)
    @test_throws ArgumentError recommended_dt(g; safety = 1.5)
    @test_throws ArgumentError recommended_dt(g; integrator = :bogus)
    @test_throws ArgumentError recommended_dt(g; d_i = -1.0)
    @test_throws ArgumentError recommended_dt(g; η = -1.0)
    @test_throws ArgumentError recommended_dt(g; ηH = -1.0)
    @test_throws ArgumentError recommended_dt(g; η = NaN)
    @test_throws ArgumentError recommended_dt(g; ηH = Inf)
end

@testset "§10.3 a HybridStepper run at recommended_dt is stable" begin
    T = Float64
    n, L = 32, 2π
    g = FourierGrid((n,), (T(L),))
    N = 200 * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(3), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    fill!(st.fields.B[1], 1.0)
    k = 2π / L
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st.fields.B[2] .= 0.01 .* cos.(k .* x)
    init!(st, ps)
    NB = 4
    dt = recommended_dt(g; NB = NB, integrator = :rk4)
    E0 = magnetic_energy(st.fields.B, g) + kinetic_energy(ps)
    for _ = 1:400
        step!(st, ps, dt; NB = NB)
    end
    for c = 1:3
        @test all(isfinite, st.fields.B[c])
    end
    E1 = magnetic_energy(st.fields.B, g) + kinetic_energy(ps)
    @test isfinite(E1)
    @test E1 < 5 * E0                                  # bounded (no whistler blow-up)
end

@testset "§10.3 2-D HybridStepper at recommended_dt is stable (diagonal B0)" begin
    # Regression guard for the multi-D CFL: with B0 along the in-plane diagonal
    # the stiffest whistler mode is the k-space CORNER (|k| = √2·π/Δx here), a
    # geometry the 1-D oracle above cannot probe. With the axis-Nyquist kmax
    # formula this run went NON-FINITE within a few steps at its own
    # recommended_dt; the corner-mode formula must keep it bounded.
    T = Float64
    n, L = 32, 2π
    g = FourierGrid((n, n), (T(L), T(L)))
    N = 8 * n * n
    ps = ParticleSet{2,T}(N)
    load_uniform!(ps, MersenneTwister(5), (0.0, 0.0), (T(L), T(L)))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(6), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    b0 = 1 / sqrt(2)
    fill!(st.fields.B[1], b0)                          # B0 ∥ (1,1,0)/√2: corner-aligned
    fill!(st.fields.B[2], b0)
    k = 2π / L
    x = [(i - 1) * g.dx[1] for i = 1:n]
    for j = 1:n
        st.fields.B[3][:, j] .= 0.01 .* cos.(k .* x)   # small transverse seed
    end
    init!(st, ps)
    NB = 4
    dt = recommended_dt(g; NB = NB, integrator = :rk4)
    E0 = magnetic_energy(st.fields.B, g) + kinetic_energy(ps)
    for _ = 1:300
        step!(st, ps, dt; NB = NB)
    end
    for c = 1:3
        @test all(isfinite, st.fields.B[c])
    end
    E1 = magnetic_energy(st.fields.B, g) + kinetic_energy(ps)
    @test isfinite(E1)
    @test E1 < 5 * E0                                  # bounded (no corner-mode blow-up)
end

@testset "§10.3 CAM-CL at recommended_dt(:camcl) is stable (stiffest mode)" begin
    # Regression guard for the cyclic-leapfrog stability fix: the pre-fix
    # restart-every-step CL amplified the stiffest active whistler mode (k = 31,
    # ω·h ≈ 1.5 at this dt) by ≈1.02/step — from a 1e-4 seed that is non-finite
    # long before 1500 steps. The persistent-CL scheme must stay bounded.
    T = Float64
    n, L = 64, 2π
    g = FourierGrid((n,), (T(L),))
    N = 200 * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(7), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = CAMCLStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    fill!(st.fields.B[1], 1.0)
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st.fields.B[2] .= 1e-4 .* cos.(31 .* x)            # stiffest active whistler mode
    st.fields.B[3] .= 1e-4 .* sin.(31 .* x)
    init_camcl!(st, ps)
    NB = 4
    dt = recommended_dt(g; NB = NB, integrator = :camcl)
    E0 = magnetic_energy(st.fields.B, g) + kinetic_energy(ps)
    for _ = 1:1500
        step_camcl!(st, ps, dt; NB = NB)
    end
    for c = 1:3
        @test all(isfinite, st.fields.B[c])
    end
    E1 = magnetic_energy(st.fields.B, g) + kinetic_energy(ps)
    @test isfinite(E1)
    @test E1 < 5 * E0                                  # bounded (no CL restart pumping)
    @test maximum(abs, st.fields.B[3]) < 0.1           # seed stayed linear, no blow-up
end
