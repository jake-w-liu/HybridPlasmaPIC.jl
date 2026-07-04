# Hybrid field model: HYB-001 algebraic equilibrium plus each generalized
# Ohm's-law term, Faraday, divergence control, and the density floor — every
# term checked against a closed-form oracle (band-limited fields → spectrally
# exact). Full time-evolution HYB-001..008 come with the integrator (Phase 5).

using HybridPlasmaPIC, Test, LinearAlgebra, Random, Statistics

coords1d(g) = [(i - 1) * g.dx[1] for i = 1:g.n[1]]

@testset "electron closure parameter validation" begin
    @test IsothermalElectrons(0.0).Te == 0.0
    @test IsothermalElectrons(0.5).Te == 0.5
    @test_throws ArgumentError IsothermalElectrons(NaN)
    @test_throws ArgumentError IsothermalElectrons(Inf)
    @test_throws ArgumentError IsothermalElectrons(-1.0)

    clo = PolytropicElectrons(0.5, 1.0, 5 / 3)
    @test clo.pe0 == 0.5
    @test clo.n0 == 1.0
    @test clo.γ ≈ 5 / 3
    @test PolytropicElectrons(0.0, 1.0, 1.0).pe0 == 0.0
    @test_throws ArgumentError PolytropicElectrons(NaN, 1.0, 5 / 3)
    @test_throws ArgumentError PolytropicElectrons(-0.1, 1.0, 5 / 3)
    @test_throws ArgumentError PolytropicElectrons(0.5, 0.0, 5 / 3)
    @test_throws ArgumentError PolytropicElectrons(0.5, -1.0, 5 / 3)
    @test_throws ArgumentError PolytropicElectrons(0.5, 1.0, NaN)
    @test_throws ArgumentError PolytropicElectrons(0.5, 1.0, 0.0)
    @test_throws ArgumentError PolytropicElectrons(0.5, 1.0, -1.0)
end

@testset "HybridModel parameter validation" begin
    model = HybridModel(IsothermalElectrons(0.0); η = 0.1, ηH = 0.05, nfloor = 1e-3)
    @test model.η == 0.1
    @test model.ηH == 0.05
    @test model.nfloor == 1e-3
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.0); η = NaN)
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.0); η = -0.1)
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.0); ηH = NaN)
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.0); ηH = -0.1)
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.0); nfloor = 0.0)
    @test_throws ArgumentError HybridModel(IsothermalElectrons(0.0); nfloor = NaN)
end

@testset "spatial dimension validation" begin
    g4 = FourierGrid((2, 2, 2, 2), (1.0, 1.0, 1.0, 1.0))
    model = HybridModel(IsothermalElectrons(0.0))
    @test_throws ArgumentError HybridFields{4,Float64}(g4.n)
    @test_throws ArgumentError HybridFields{4,Float64}(g4.n; anisotropic = true)
    @test_throws ArgumentError HybridFields{true,Float64}((2,))
    @test_throws ArgumentError HybridFields{Int8(1),Float64}((2,))
    @test_throws ArgumentError HybridFields{UInt(2),Float64}((2, 2))
    @test_throws ArgumentError HybridStepper(g4, model, NGP(), 1)
    @test_throws ArgumentError CAMCLStepper(g4, model, NGP(), 1)
end

@testset "HYB-001 algebraic uniform equilibrium" begin
    for D = 1:3
        T = Float64
        nc = ntuple(_ -> 16, D)
        g = FourierGrid(nc, ntuple(_ -> 2π, D))
        f = HybridFields{D,T}(nc)
        fill!(f.n, 1.0)
        fill!(f.B[3], 1.0)                       # ui=0, B=ẑ
        ohms_law!(f, HybridModel(IsothermalElectrons(0.5)), g)
        for c = 1:3
            @test maximum(abs, f.J[c]) < 1e-12
            @test maximum(abs, f.E[c]) < 1e-12
        end
        @test f.floor_count[] == 0
    end
end

@testset "Ohm motional term E = −u_i×B" begin
    for D = 1:3
        T = Float64
        nc = ntuple(_ -> 16, D)
        g = FourierGrid(nc, ntuple(_ -> 2π, D))
        f = HybridFields{D,T}(nc)
        fill!(f.n, 1.0)
        u0 = (0.3, -0.1, 0.2)
        B0 = (0.2, 0.5, 1.0)
        for c = 1:3
            fill!(f.ui[c], u0[c])
            fill!(f.B[c], B0[c])
        end
        ohms_law!(f, HybridModel(IsothermalElectrons(0.0)), g)   # Te=0, uniform ⇒ only motional
        ex = -(u0[2] * B0[3] - u0[3] * B0[2])
        ey = -(u0[3] * B0[1] - u0[1] * B0[3])
        ez = -(u0[1] * B0[2] - u0[2] * B0[1])
        @test maximum(abs, f.E[1] .- ex) < 1e-12
        @test maximum(abs, f.E[2] .- ey) < 1e-12
        @test maximum(abs, f.E[3] .- ez) < 1e-12
    end
end

@testset "Ohm Hall term (1D analytic J×B/n)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    f = HybridFields{1,T}((n,))
    fill!(f.n, 1.0)
    x = coords1d(g)
    k = 3.0
    f.B[3] .= cos.(k .* x)                       # Bz=cos(kx) ⇒ J=(0, k sin kx, 0)
    ohms_law!(f, HybridModel(IsothermalElectrons(0.0)), g)
    @test maximum(abs, f.J[2] .- (k .* sin.(k .* x))) < 1e-10
    @test maximum(abs, f.E[1] .- (@. k * sin(k * x) * cos(k * x))) < 1e-10
    @test maximum(abs, f.E[2]) < 1e-10
    @test maximum(abs, f.E[3]) < 1e-10
end

@testset "Ohm electron-pressure term (1D analytic)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    f = HybridFields{1,T}((n,))
    x = coords1d(g)
    k = 2.0
    amp = 0.2
    Te = 0.5
    f.n .= 1 .+ amp .* cos.(k .* x)              # B=0, ui=0 ⇒ only −∇p_e/n
    ohms_law!(f, HybridModel(IsothermalElectrons(Te)), g)
    Ex = @. -(Te * (-amp * k * sin(k * x))) / (1 + amp * cos(k * x))
    @test maximum(abs, f.E[1] .- Ex) < 1e-10
    @test maximum(abs, f.E[2]) < 1e-10
    @test maximum(abs, f.E[3]) < 1e-10
end

@testset "Faraday dB/dt = −∇×E (1D analytic)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    x = coords1d(g)
    k = 3.0
    E = (zeros(T, n), zeros(T, n), cos.(k .* x))  # Ez=cos(kx) ⇒ ∇×E=(0, k sin kx, 0)
    dB = (zeros(T, n), zeros(T, n), zeros(T, n))
    faraday_rhs!(dB, E, g)
    @test maximum(abs, dB[1]) < 1e-10
    @test maximum(abs, dB[2] .- (.-k .* sin.(k .* x))) < 1e-10
    @test maximum(abs, dB[3]) < 1e-10
end

@testset "divergence control via project_b!" begin
    for D = 1:3
        T = Float64
        nc = ntuple(_ -> 16, D)
        g = FourierGrid(nc, ntuple(_ -> 2π, D))
        f = HybridFields{D,T}(nc)
        for c = 1:3
            f.B[c] .= randn(MersenneTwister(10c + D), nc...)
        end
        out = zeros(T, nc)
        d0 = magnetic_divergence!(out, f, g)
        project_b!(f, g)
        d1 = magnetic_divergence!(out, f, g)
        kmax = maximum(maximum(abs, g.kvec[d]) for d = 1:D)
        scale = sum(norm(f.B[c]) for c = 1:3)
        @test d1 / (kmax * scale) < 1e-10
        @test d1 < d0
    end
end

@testset "density floor activation + finiteness" begin
    T = Float64
    n = 16
    g = FourierGrid((n,), (2π,))
    f = HybridFields{1,T}((n,))
    f.n .= 1.0
    f.n[1] = 1e-12
    fill!(f.B[3], 1.0)
    fill!(f.ui[1], 0.5)
    model = HybridModel(IsothermalElectrons(0.5); nfloor = 1e-3)
    ohms_law!(f, model, g)
    @test f.floor_count[] >= 1
    @test all(isfinite, f.E[1]) && all(isfinite, f.E[2]) && all(isfinite, f.E[3])
end

@testset "Ohm's-law defensive shape checks" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    model = HybridModel(IsothermalElectrons(0.5))

    f = HybridFields{1,T}((8,))
    f.ninv = zeros(T, 7)
    @test_throws DimensionMismatch ohms_law!(f, model, g)

    f = HybridFields{1,T}((8,))
    f.ui = (zeros(T, 7), zeros(T, 8), zeros(T, 8))
    @test_throws DimensionMismatch ohms_law!(f, model, g)

    f = HybridFields{1,T}((8,))
    f.E = (zeros(T, 7), zeros(T, 8), zeros(T, 8))
    @test_throws DimensionMismatch ohms_law!(f, model, g)
end

@testset "step! validates timestep and subcycles before mutation" begin
    T = Float64
    g = FourierGrid((4,), (2π,))
    ps = ParticleSet{1,T}(4)
    load_lattice_1d!(ps, 0.0, 2π)
    set_density_weight!(ps, 1.0, g)
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.0)), NGP(), nparticles(ps))
    fill!(st.fields.B[3], 1.0)
    init!(st, ps)

    x0 = ntuple(d -> copy(ps.x[d]), 1)
    v0 = ntuple(c -> copy(ps.v[c]), 3)
    B0 = ntuple(c -> copy(st.fields.B[c]), 3)
    time0 = st.time[]
    step0 = st.step[]

    @test_throws ArgumentError step!(st, ps, 0.1; NB = 0)
    @test_throws ArgumentError step!(st, ps, NaN; NB = 1)
    @test_throws ArgumentError step!(st, ps, -0.1; NB = 1)
    @test st.time[] == time0
    @test st.step[] == step0
    @test all(ps.x[d] == x0[d] for d = 1:1)
    @test all(ps.v[c] == v0[c] for c = 1:3)
    @test all(st.fields.B[c] == B0[c] for c = 1:3)
end

@testset "step! resizes particle workspaces after particle count changes" begin
    T = Float64
    g = FourierGrid((8,), (1.0,))
    ps = ParticleSet{1,T}(4)
    ps.x[1] .= T[0.2, 0.4, 0.6, 0.8]
    ps.weight .= T(1 / 8)
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.0)), CIC(), nparticles(ps))
    st.fields.B[1] .= 1.0
    init!(st, ps)

    ps.x[1][1] = -0.1
    ps.x[1][2] = 1.1
    @test apply_absorbing!(ps, (0.0,), (1.0,)) == 2
    @test step!(st, ps, 0.01) === st
    @test length(st.Ep[1]) == nparticles(ps)
    @test length(st.work) == nparticles(ps)

    ps2 = ParticleSet{1,T}(0)
    st2 = HybridStepper(g, HybridModel(IsothermalElectrons(0.0)), CIC(), nparticles(ps2))
    st2.fields.B[1] .= 1.0
    init!(st2, ps2)
    acc = Ref(0.0)
    nextid = Ref(UInt64(1))
    @test inject_face_1d!(
        ps2,
        MersenneTwister(4),
        0.0,
        +1,
        1.0,
        2.0,
        0.0,
        (0.0, 0.0),
        0.3,
        1.0,
        0.5,
        acc,
        nextid,
    ) == 4
    @test step!(st2, ps2, 0.01) === st2
    @test length(st2.Ep[1]) == nparticles(ps2)
    @test length(st2.work) == nparticles(ps2)
end

@testset "compute_moments! from particles" begin
    T = Float64
    n = 32
    L = 2π
    g = FourierGrid((n,), (L,))
    nppc = 200
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_uniform!(ps, MersenneTwister(1), (0.0,), (L,))
    load_maxwellian!(ps, MersenneTwister(2), (0.1, 0.0, 0.0), (0.5, 0.5, 0.5))
    f = HybridFields{1,T}((n,))
    compute_moments!(f, ps, g, CIC(), 1e-6)
    @test isapprox(mean(f.n), nppc / g.dx[1]; rtol = 0.02)
    @test isapprox(mean(f.ui[1]), 0.1; atol = 0.02)
    @test isapprox(mean(f.ui[2]), 0.0; atol = 0.02)
end

@testset "compute_moments! validates density floor before mutation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    ps = ParticleSet{1,T}(0)
    f = HybridFields{1,T}((8,))
    fill!(f.n, 1.0)
    fill!(f.ui[1], 2.0)
    fill!(f.ui[2], 3.0)
    fill!(f.ui[3], 4.0)

    @test_throws ArgumentError compute_moments!(f, ps, g, NGP(), 0.0)
    @test all(==(1.0), f.n)
    @test all(==(2.0), f.ui[1])
    @test all(==(3.0), f.ui[2])
    @test all(==(4.0), f.ui[3])

    @test_throws ArgumentError compute_moments!(f, ps, g, NGP(), NaN)
    @test all(==(1.0), f.n)
    @test all(==(2.0), f.ui[1])
    @test all(==(3.0), f.ui[2])
    @test all(==(4.0), f.ui[3])
end

@testset "HYB-order: magnetized hybrid step! is 2nd-order in dt" begin
    # The carried-E u_i re-centering (predictor) + one-time leapfrog velocity priming make a
    # real init!+step! run 2nd-order accurate on a magnetized problem (previously 1st-order:
    # the carried E used the half-step-lagged u_i^{n+1/2}, and the loaded v^0 was used as
    # v^{-1/2} unprimed). Self-convergence on the seeded transverse magnetic energy: the
    # temporal order must sit near 2, well above the old 1.
    T = Float64
    L = 2π
    k = 2π / L
    function bE(nsteps, Tf)
        n = 16
        N = 1000 * n
        g = FourierGrid((n,), (T(L),))
        ps = ParticleSet{1,T}(N)
        load_lattice_1d!(ps, 0.0, T(L))
        set_density_weight!(ps, 1.0, g)
        load_quiet_velocities!(ps, MersenneTwister(7), (0.0, 0.0, 0.0), (0.02, 0.02, 0.02))
        st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
        x = [(i - 1) * g.dx[1] for i = 1:n]
        fill!(st.fields.B[1], 1.0)
        st.fields.B[2] .= 0.01 .* cos.(k .* x)
        st.fields.B[3] .= 0.01 .* sin.(k .* x)
        init!(st, ps)
        dt = Tf / nsteps
        for _ = 1:nsteps
            step!(st, ps, dt; NB = 4)
        end
        magnetic_energy(st.fields.B, g)
    end
    Tf = 0.2
    s = (10, 20, 40, 80)
    v = [bE(ns, Tf) for ns in s]
    r1 = log2(abs(v[1] - v[2]) / abs(v[2] - v[3]))
    r2 = log2(abs(v[2] - v[3]) / abs(v[3] - v[4]))
    @test r1 > 1.6        # ≈2.0 with the fix; ≈1.0 (would fail) without it
    @test r2 > 1.6
end

@testset "step! dt=0 is a true no-op (preserves state and the priming guard)" begin
    T = Float64
    L = 2π
    n = 16
    N = 2000
    g = FourierGrid((n,), (T(L),))
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(3), (0.0, 0.0, 0.0), (0.02, 0.02, 0.02))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    fill!(st.fields.B[1], 1.0)
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st.fields.B[2] .= 0.01 .* cos.(2π / L .* x)
    init!(st, ps)
    x0 = copy(ps.x[1])
    v0 = copy(ps.v[1])
    B0 = copy(st.fields.B[2])
    s0 = st.step[]
    step!(st, ps, 0.0; NB = 4)
    # a dt=0 step advances nothing AND must NOT consume the one-time leapfrog priming (step==0),
    # which would otherwise silently drop the run back to 1st-order.
    @test ps.x[1] == x0 && ps.v[1] == v0 && st.fields.B[2] == B0
    @test st.step[] == s0
end
