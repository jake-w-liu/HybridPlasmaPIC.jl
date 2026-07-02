# Finite-electron-mass (electron inertia) in the hybrid Ohm's law: the leading-order term
# is the spectral multiplier E ← E/(1 + d_e² k²) (regularizes the whistler / reconnection at
# the electron inertial scale d_e). de2=0 recovers the massless hybrid exactly (no-op).

using HybridPlasmaPIC, Test, Random
using HybridPlasmaPIC:
    FourierGrid,
    _inertia_filter!,
    HybridModel,
    HybridStepper,
    ParticleSet,
    load_uniform!,
    load_maxwellian!,
    set_density_weight!,
    init!,
    step!

@testset "INERTIA-001 electron-inertia filter E ← E/(1+d_e²k²) vs analytic" begin
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
