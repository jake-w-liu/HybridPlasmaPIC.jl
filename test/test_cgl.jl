# CGL anisotropic (double-adiabatic) electron closure: the invariant laws, the pressure-
# tensor divergence ∇·P_e vs analytic references, and end-to-end stability in the
# well-posed regime. The scalar Ohm's-law path is unchanged (compile-time dispatch on the
# closure type) — that no-regression is covered by test_hybrid.jl / test_instabilities.jl.

using HybridPlasmaPIC, Test, Random
using HybridPlasmaPIC:
    FourierGrid,
    deriv!,
    HybridModel,
    HybridStepper,
    CAMCLStepper,
    electron_internal_energy,
    ParticleSet,
    load_uniform!,
    load_maxwellian!,
    set_density_weight!,
    init!,
    step!

@testset "CGL-001 double-adiabatic invariant laws" begin
    c = CGLElectrons(0.6, 0.3, 1.0, 1.0)                 # C⊥ = 0.6, C∥ = 0.3
    for (n, B) in ((1.0, 1.0), (2.0, 1.0), (1.0, 2.0), (0.5, 1.5), (3.0, 0.7))
        @test cgl_pperp(c, n, B) / (n * B) ≈ 0.6 rtol = 1e-14      # μ invariant  p⊥/(nB)
        @test cgl_ppar(c, n, B) * B^2 / n^3 ≈ 0.3 rtol = 1e-14     # J invariant  p∥B²/n³
    end
    @test cgl_pperp(c, 2.0, 1.0) ≈ 2 * cgl_pperp(c, 1.0, 1.0)      # p⊥ ∝ nB
    @test cgl_ppar(c, 2.0, 1.0) ≈ 8 * cgl_ppar(c, 1.0, 1.0)        # p∥ ∝ n³
    @test cgl_pperp(c, 1.0, 2.0) ≈ 2 * cgl_pperp(c, 1.0, 1.0)      # p⊥ ∝ B
    @test cgl_ppar(c, 1.0, 2.0) ≈ cgl_ppar(c, 1.0, 1.0) / 4        # p∥ ∝ 1/B²
    ci = CGLElectrons(0.5, 0.5, 1.0, 1.0)
    @test cgl_pperp(ci, 1.0, 1.0) ≈ 0.5 && cgl_ppar(ci, 1.0, 1.0) ≈ 0.5   # isotropic reference
    @test is_anisotropic(ci) && !is_anisotropic(IsothermalElectrons(0.5))
    @test_throws ArgumentError CGLElectrons(0.6, 0.3, -1.0, 1.0)
    @test_throws ArgumentError CGLElectrons(-0.6, 0.3, 1.0, 1.0)

    extreme_perp = CGLElectrons(1.0e308, 0.0, 1.0e308, 2.0)
    @test extreme_perp.Cperp == 0.5
    extreme_par = CGLElectrons(0.0, 1.0e-300, 1.0e30, 1.0e200)
    @test extreme_par.Cpar ≈ 1.0e10 rtol = 2eps(Float64)
end

@testset "CGL-002 pressure-tensor divergence ∇·P_e vs analytic" begin
    T = Float64
    N = 32
    Lx = 2π
    g = FourierGrid((N,), (Lx,))
    xs = [(i - 1) * Lx / N for i = 1:N]
    ε = 0.3
    kx = 2π / Lx
    n = [1.0 * (1 + ε * sin(kx * x)) for x in xs]
    c = CGLElectrons(0.6, 0.3, 1.0, 1.0)
    Fp = (zeros(N), zeros(N), zeros(N))
    # B ∥ ẑ (uniform): b ⊥ ∇n ⇒ ∇·P_e = (∂ₓp⊥, 0, 0)
    anisotropic_pressure_force!(Fp, n, (zeros(N), zeros(N), fill(2.0, N)), c, g)
    ref = similar(n)
    deriv!(ref, [cgl_pperp(c, ni, 2.0) for ni in n], g, 1)
    @test maximum(abs.(Fp[1] .- ref)) < 1e-12
    @test maximum(abs.(Fp[2])) < 1e-14 && maximum(abs.(Fp[3])) < 1e-14
    # B ∥ x̂ (uniform): b ∥ ∇n ⇒ ∇·P_e = (∂ₓp∥, 0, 0) — exercises the field-aligned stress
    anisotropic_pressure_force!(Fp, n, (fill(2.0, N), zeros(N), zeros(N)), c, g)
    ref2 = similar(n)
    deriv!(ref2, [cgl_ppar(c, ni, 2.0) for ni in n], g, 1)
    @test maximum(abs.(Fp[1] .- ref2)) < 1e-12
    @test maximum(abs.(Fp[2])) < 1e-14 && maximum(abs.(Fp[3])) < 1e-14
    # uniform n, B ⇒ zero force
    anisotropic_pressure_force!(Fp, fill(1.3, N), (fill(0.5, N), fill(0.5, N), fill(0.7, N)), c, g)
    @test maximum(abs.(Fp[1]) .+ abs.(Fp[2]) .+ abs.(Fp[3])) < 1e-14

    # Squaring this finite uniform field overflows, and an unscaled FFT of the
    # finite pressure DC mode does too. Its analytic pressure force is exactly zero.
    huge_B = (fill(1.0e308, N), zeros(N), zeros(N))
    anisotropic_pressure_force!(Fp, ones(N), huge_B, CGLElectrons(0.5, 0.5, 1.0, 1.0), g)
    @test all(all(isfinite, Fp[c]) for c = 1:3)
    @test all(all(iszero, Fp[c]) for c = 1:3)

    # The Fourier coefficients remain finite here, but multiplying them by the
    # short-domain wavenumber overflows unless the rescaling bound includes the
    # spectral multiplier as well as the unnormalised transform.
    high_k = 100.0
    high_N = 16
    high_L = 2π / high_k
    high_g = FourierGrid((high_N,), (high_L,))
    high_x = (0:high_N-1) .* (high_L / high_N)
    high_n = 2.0e306 .+ 1.0e306 .* sin.(high_k .* high_x)
    high_Fp = ntuple(_ -> zeros(high_N), 3)
    anisotropic_pressure_force!(
        high_Fp,
        high_n,
        (zeros(high_N), zeros(high_N), ones(high_N)),
        CGLElectrons(1.0, 0.0, 1.0, 1.0),
        high_g,
    )
    high_gradient = 1.0e308 .* cos.(high_k .* high_x)
    @test all(isfinite, high_Fp[1])
    @test maximum(abs.(high_Fp[1] ./ 1.0e308 .- high_gradient ./ 1.0e308)) < 32eps()

    # Exercise the same multiplier-aware scaling through the field-aligned
    # stress divergence, independently of the perpendicular gradient.
    stress_N = 32
    stress_g = FourierGrid((stress_N,), (high_L,))
    stress_x = (0:stress_N-1) .* (high_L / stress_N)
    stress_n = 1 .+ 0.1 .* sin.(high_k .* stress_x)
    stress_Fp = ntuple(_ -> zeros(stress_N), 3)
    anisotropic_pressure_force!(
        stress_Fp,
        stress_n,
        (ones(stress_N), zeros(stress_N), zeros(stress_N)),
        CGLElectrons(0.0, 1.0e306, 1.0, 1.0),
        stress_g,
    )
    stress_reference = 3.0e307 .* stress_n .^ 2 .* cos.(high_k .* stress_x)
    @test all(isfinite, stress_Fp[1])
    @test maximum(abs.(stress_Fp[1] ./ 1.0e307 .- stress_reference ./ 1.0e307)) < 256eps()
    @test all(iszero, stress_Fp[2]) && all(iszero, stress_Fp[3])
end

@testset "CGL-003 well-posed hybrid run: stable, deterministic, distinct from scalar" begin
    T = Float64
    N = 64
    Lx = 12.8
    g = FourierGrid((N,), (Lx,))
    function runit(closure)
        m = HybridModel(closure; η = 0.001)
        st = HybridStepper(g, m, CIC(), 100 * N)
        ps = ParticleSet{1,T}(100 * N)
        rng = MersenneTwister(1)
        load_uniform!(ps, rng, (0.0,), (Lx,))
        load_maxwellian!(ps, rng, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        set_density_weight!(ps, 1.0, g)
        fill!(st.fields.B[1], 1.0)
        init!(st, ps)
        ok = 0
        for _ = 1:150
            step!(st, ps, 0.02; NB = 2)
            all(isfinite, st.fields.B[2]) || break
            ok += 1
        end
        return (ok, sum(abs2, st.fields.B[2]) + sum(abs2, st.fields.B[3]))
    end
    cgl = runit(CGLElectrons(0.55, 0.45, 1.0, 1.0))       # β⊥·A = 0.24 < 1: mirror-stable
    @test cgl[1] == 150 && isfinite(cgl[2])                # completes without blowing up
    @test runit(CGLElectrons(0.55, 0.45, 1.0, 1.0))[2] == cgl[2]   # deterministic
    iso = runit(IsothermalElectrons(0.5))
    @test iso[1] == 150
    @test cgl[2] != iso[2]                                 # the anisotropic closure changes B
end

@testset "CGL-004 robustness: size guards, diagnostics, integrator support" begin
    T = Float64
    N = 16
    g = FourierGrid((N,), (2π,))
    c = CGLElectrons(0.6, 0.3, 1.0, 1.0)
    # the anisotropic Ohm helper validates array sizes (mirrors its scalar twin)
    @test_throws DimensionMismatch anisotropic_pressure_force!(
        ntuple(_ -> zeros(N), 3),
        zeros(2N),
        ntuple(_ -> ones(N), 3),
        c,
        g,
    )
    # electron_internal_energy: the density-only method degrades to NaN for CGL (needs |B|;
    # kept for backward compatibility), while the (n, B, closure, g) method computes the
    # gyrotropic internal energy ∫ (p_⊥ + p_∥/2) dV — checked against the uniform-field value.
    @test isnan(electron_internal_energy(ones(N), c, g))
    Bu = (zeros(N), zeros(N), fill(2.0, N))
    @test electron_internal_energy(ones(N), Bu, c, g) ≈
          (cgl_pperp(c, 1.0, 2.0) + cgl_ppar(c, 1.0, 2.0) / 2) * 2π
    # CAM-CL rejects anisotropic closures at construction (frozen scalar ∇p_e is incompatible)
    @test_throws ArgumentError CAMCLStepper(g, HybridModel(c), CIC(), 100)
end

@testset "CGL-005 numeric-type generality: Float32 fields + Float64 closure" begin
    T = Float32
    g = FourierGrid((16,), (T(2π),))
    f = HybridPlasmaPIC.HybridFields{1,T}(g.n; anisotropic = true)
    fill!(f.n, one(T))
    fill!(f.B[3], one(T))
    # a CGL closure built from ordinary Float64 literals must drive Float32 fields (like the
    # scalar closures) — the anisotropic force's closure precision is independent of field T.
    model = HybridModel(CGLElectrons(0.3, 0.6, 1.0, 1.0))
    ohms_law!(f, model, g)                                    # must not MethodError
    @test eltype(f.E[1]) === T
    @test all(all(isfinite, f.E[c]) for c = 1:3)
end

@testset "CGL solver pressure workspace avoids per-call grid allocations" begin
    N = (64, 64)
    g = FourierGrid(N, (2π, 2π))
    f = HybridPlasmaPIC.HybridFields{2,Float64}(N; anisotropic = true)
    fill!(f.n, 1.0)
    fill!(f.B[1], 1.0)
    model = HybridModel(CGLElectrons(0.5, 0.4, 1.0, 1.0))
    ohms_law!(f, model, g)

    bytes = @allocated ohms_law!(f, model, g)

    @test bytes < sizeof(f.n)
    @test size(f.cgl_dp) == N
    scalar = HybridPlasmaPIC.HybridFields{2,Float64}(N)
    @test isempty(scalar.cgl_dp)
end

@testset "CGL-006 energy budget closes with the gyrotropic internal energy" begin
    # 1D compressive run ⊥ B (B0 = ẑ, bulk u_x = 0.2 sin x): KE + ∫½|B|² alone
    # swings by ~half the KE swing as compression does work against P_e; adding
    # ∫(p_⊥ + p_∥/2) dV (the electron_internal term energy_budget now computes
    # from B) closes the budget to a few per cent of the KE swing — exact modulo
    # the anisotropic battery term (measured: 0.505·sKE → 0.046·sKE, corr 0.998).
    T = Float64
    N = 32
    L = 2π
    npc = 300
    g = FourierGrid((N,), (L,))
    clo = CGLElectrons(0.5, 0.5, 1.0, 1.0)
    Np = npc * N
    ps = ParticleSet{1,T}(Np)
    load_lattice_1d!(ps, 0.0, L)
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(clo), CIC(), Np)
    fill!(st.fields.B[3], 1.0)
    ps.v[1] .+= 0.2 .* sin.((2π / L) .* ps.x[1])       # compressive bulk flow ⊥ B
    init!(st, ps)
    b0 = energy_budget(ps, st.fields.B, st.fields.n, clo, g)
    @test isfinite(b0.total) && b0.electron_internal > 0
    @test b0.electron_internal == electron_internal_energy(st.fields.n, st.fields.B, clo, g)
    KE = [kinetic_energy(ps)]
    EB = [magnetic_energy(st.fields.B, g)]
    EI = [electron_internal_energy(st.fields.n, st.fields.B, clo, g)]
    for _ = 1:300
        step!(st, ps, 0.02; NB = 4)
        push!(KE, kinetic_energy(ps))
        push!(EB, magnetic_energy(st.fields.B, g))
        push!(EI, electron_internal_energy(st.fields.n, st.fields.B, clo, g))
    end
    swing(v) = maximum(v) - minimum(v)
    sKE = swing(KE)
    sKEB = swing(KE .+ EB)
    sTOT = swing(KE .+ EB .+ EI)
    @info "CGL budget swings" sKE sKEB sTOT
    @test sKEB > 0.3 * sKE                  # without E_int the budget is open
    @test sTOT < 0.15 * sKE                 # with E_int it closes (measured 4.6%)
    @test sTOT < 0.3 * sKEB                 # the E_int term does the closing
end
