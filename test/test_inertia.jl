# Finite-electron-mass (electron inertia) in the hybrid Ohm's law: the term is
# (d_e²/n)∂_t J with J = ∇×B, which is divergence-free, so it filters ONLY the
# k-transverse projection of E: Ê_⊥ ← Ê_⊥/(1 + d_e² k²) with Ê_∥ untouched
# (regularizes the whistler / reconnection at the electron inertial scale d_e).
# de2=0 recovers the massless hybrid exactly (no-op). The decisive oracles here:
# a purely longitudinal E (B=0, ∇n≠0) must be independent of de2, and at an
# oblique k only the transverse projection may be divided by (1+d_e²k²).

using HybridPlasmaPIC, Test, Random
using HybridPlasmaPIC:
    FourierGrid,
    _inertia_filter!,
    _apply_electron_inertia!,
    HybridModel,
    HybridStepper,
    HybridFields,
    ParticleSet,
    load_uniform!,
    load_maxwellian!,
    set_density_weight!,
    ohms_law!,
    init!,
    step!

@testset "INERTIA-001 scalar transverse multiplier f ← f/(1+d_e²k²) vs analytic" begin
    T = Float64
    N = 32
    Lx = 2π
    g = FourierGrid((N,), (Lx,))
    xs = [(i - 1) * Lx / N for i = 1:N]
    for (m, de2) in ((2, 0.5), (3, 1.0), (5, 0.2))
        k = m * (2π / Lx)
        f = [sin(k * x) for x in xs]
        _inertia_filter!(f, de2, g)                        # single Fourier mode ⇒ exact factor
        @test maximum(abs.(f .- [sin(k * x) / (1 + de2 * k^2) for x in xs])) < 1e-13
    end
    f = [cos(2 * (2π / Lx) * x) + 0.3 for x in xs]
    f0 = copy(f)
    _inertia_filter!(f, 0.0, g)                            # de2=0 ⇒ identity
    @test maximum(abs.(f .- f0)) < 1e-13
    # pure Nyquist mode: the discrete ∇×∇× (Nyquist-zeroed first derivatives) vanishes
    # identically there, so the filter must pass it unchanged
    fny = [cos(π * (i - 1)) for i = 1:N]                   # (−1)^i, k = k_Nyq
    @test maximum(abs.(_inertia_filter!(copy(fny), 1.0, g) .- fny)) < 1e-13
end

@testset "INERTIA-003 longitudinal (electrostatic) E is de2-independent" begin
    # B = 0 ⇒ J = ∇×B = 0 and ∂_t J = 0: NO inertia correction can physically exist.
    # E = −Te ∇n/n is purely longitudinal (k ∥ x̂ ∥ E), so ohms_law! must return the
    # SAME field for any de2 — this is the oracle the transverse projection enforces.
    T = Float64
    N = 32
    Lx = 2π
    g = FourierGrid((N,), (Lx,))
    xs = [(i - 1) * Lx / N for i = 1:N]
    function runE(de2)
        f = HybridFields{1,T}((N,))
        f.n .= 1 .+ 0.3 .* sin.(xs)                        # u_i = 0, B = 0
        ohms_law!(f, HybridModel(IsothermalElectrons(1.0); de2 = de2), g)
        return deepcopy(f.E)
    end
    E0 = runE(0.0)
    Ean = @. -0.3 * cos(xs) / (1 + 0.3 * sin(xs))          # analytic −Te ∇n/n
    @test maximum(abs.(E0[1] .- Ean)) < 1e-13
    for de2 in (0.5, 1.0)
        E1 = runE(de2)
        @test maximum(abs.(E1[1] .- E0[1])) < 1e-13        # longitudinal: UNfiltered
        @test maximum(abs.(E1[2] .- E0[2])) < 1e-13
        @test maximum(abs.(E1[3] .- E0[3])) < 1e-13
    end
end

@testset "INERTIA-004 1D split: Ex (∥) passes, Ey/Ez (⊥) get exactly 1/(1+d_e²k²)" begin
    # In 1D every mode has k ∥ x̂, so Ex is the longitudinal projection (must pass
    # unfiltered) while Ey, Ez are transverse (must match the scalar multiplier
    # _inertia_filter! exactly, mode by mode).
    T = Float64
    N = 64
    Lx = 4π
    de2 = 0.7
    g = FourierGrid((N,), (Lx,))
    xs = [(i - 1) * Lx / N for i = 1:N]
    k1 = 2π / Lx
    function runE(d2)
        f = HybridFields{1,T}((N,))
        f.n .= 1 .+ 0.2 .* cos.(k1 .* xs)
        f.B[1] .= 1.0
        f.B[2] .= 0.1 .* cos.(3k1 .* xs) .+ 0.05 .* sin.(5k1 .* xs)
        f.B[3] .= 0.1 .* sin.(3k1 .* xs) .- 0.02 .* cos.(7k1 .* xs)
        f.ui[1] .= 0.05 .* sin.(k1 .* xs)
        f.ui[2] .= 0.03
        f.ui[3] .= -0.04 .* cos.(2k1 .* xs)
        ohms_law!(f, HybridModel(IsothermalElectrons(0.5); η = 0.01, de2 = d2), g)
        return deepcopy(f.E)
    end
    E0 = runE(0.0)
    E1 = runE(de2)
    @test maximum(abs.(E1[1] .- E0[1])) < 1e-13            # longitudinal component untouched
    for c = 2:3                                            # transverse = scalar multiplier
        @test maximum(abs.(E1[c] .- _inertia_filter!(copy(E0[c]), de2, g))) < 1e-13
    end
end

@testset "INERTIA-005 2D oblique k: only the transverse projection is filtered" begin
    # Single oblique mode k = (3,2)·2π/L with known longitudinal amplitude aL (along k̂)
    # and transverse amplitudes aT1 (in-plane ⊥ k̂) and aT2 (ẑ): after the filter the
    # longitudinal part must be untouched and both transverse parts divided by (1+d_e²k²).
    T = Float64
    N = 32
    L = 2π
    de2 = 0.8
    g = FourierGrid((N, N), (L, L))
    kx = 3 * (2π / L)
    ky = 2 * (2π / L)
    k2 = kx^2 + ky^2
    kh = (kx, ky) ./ sqrt(k2)                              # k̂ (longitudinal direction)
    th = (-kh[2], kh[1])                                   # in-plane transverse direction
    aL, aT1, aT2 = 0.9, 0.4, -0.6
    E = ntuple(_ -> zeros(T, N, N), 3)
    for j = 1:N, i = 1:N
        ph = kx * (i - 1) * L / N + ky * (j - 1) * L / N
        E[1][i, j] = (aL * kh[1] + aT1 * th[1]) * cos(ph)
        E[2][i, j] = (aL * kh[2] + aT1 * th[2]) * cos(ph)
        E[3][i, j] = aT2 * sin(ph)
    end
    _apply_electron_inertia!(E, de2, g)
    mult = 1 / (1 + de2 * k2)
    err = 0.0
    for j = 1:N, i = 1:N
        ph = kx * (i - 1) * L / N + ky * (j - 1) * L / N
        err = max(
            err,
            abs(E[1][i, j] - (aL * kh[1] + mult * aT1 * th[1]) * cos(ph)),
            abs(E[2][i, j] - (aL * kh[2] + mult * aT1 * th[2]) * cos(ph)),
            abs(E[3][i, j] - mult * aT2 * sin(ph)),
        )
    end
    @test err < 1e-13
    # de2=0 vector path is a bitwise no-op
    rng = MersenneTwister(7)
    E2 = ntuple(_ -> randn(rng, T, N, N), 3)
    E2c = deepcopy(E2)
    _apply_electron_inertia!(E2, 0.0, g)
    @test all(E2[c] == E2c[c] for c = 1:3)
end

@testset "INERTIA-002 massless no-op + de2>0 stable, deterministic, distinct" begin
    T = Float64
    N = 64
    Lx = 12.8
    g = FourierGrid((N,), (Lx,))
    function runit(de2)
        m = HybridModel(IsothermalElectrons(0.5); η = 0.001, de2 = de2)
        st = HybridStepper(g, m, CIC(), 100 * N)
        ps = ParticleSet{1,T}(100 * N)
        rng = MersenneTwister(1)
        load_uniform!(ps, rng, (0.0,), (Lx,))
        load_maxwellian!(ps, rng, (0.0, 0.0, 0.0), (1.4, 0.5, 0.5))    # anisotropic (firehose) ions
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
    massless = runit(0.0)
    @test massless[1] == 150
    inertial = runit(0.5)
    @test inertial[1] == 150 && isfinite(inertial[2])       # stable
    @test inertial[2] != massless[2]                        # electron inertia changes the field
    @test runit(0.5)[2] == inertial[2]                      # deterministic
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.5); de2 = -1.0)
end
