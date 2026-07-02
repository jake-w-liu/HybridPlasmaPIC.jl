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
    # electron_internal_energy degrades to NaN for CGL (no scalar pressure) rather than throwing
    @test isnan(electron_internal_energy(ones(N), c, g))
    # CAM-CL rejects anisotropic closures at construction (frozen scalar ∇p_e is incompatible)
    @test_throws ArgumentError CAMCLStepper(g, HybridModel(c), CIC(), 100)
end
