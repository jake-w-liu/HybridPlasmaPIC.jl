# Phase-4 forcing operators: Landau-Lifshitz radiation reaction (synchrotron cooling) and
# the external antenna/RF field source (divergence-free Faraday injection). Uniform applied
# fields + gravity + a uniform RF drive are already provided by push_uniform!.

using HybridPlasmaPIC, Test, Random
using HybridPlasmaPIC: FourierGrid, ParticleSet, divergence!, curl!

@testset "RR-001 radiation reaction: synchrotron cooling rate vs analytic" begin
    T = Float64
    γof(vx) = 1 / sqrt(1 - vx^2)
    K = 0.001
    B0 = 1.0
    dt = 0.05
    ps = ParticleSet{1,T}(1)
    ps.v[1][1] = 0.995                                  # ultrarelativistic, ⊥ B
    γ0 = γof(ps.v[1][1])
    inv = Float64[1/γ0]
    for _ = 1:60
        apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, B0), dt; K = K)
        push!(inv, 1 / γof(ps.v[1][1]))
    end
    @test γ0 > 9                                        # started ultrarelativistic
    @test all(diff(inv) .>= 0)                          # 1/γ increases ⇒ monotone cooling
    # early-time slope d(1/γ)/dt ≈ K·B²·v0⁴ (Landau-Lifshitz synchrotron)
    slope = (inv[21] - inv[1]) / (20 * dt)
    @test isapprox(slope, K * B0^2 * 0.995^4; rtol = 0.02)
end

@testset "RR-002 radiation reaction: no-op, guards, determinism" begin
    T = Float64
    ps = ParticleSet{1,T}(10)
    for i = 1:10
        ps.v[1][i] = 0.5
    end
    snap = copy(ps.v[1])
    apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.1; K = 0.0)   # K=0
    @test ps.v[1] == snap
    apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.0; K = 1.0)   # dt=0
    @test ps.v[1] == snap
    @test_throws ArgumentError apply_radiation_reaction!(
        ps,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = -1.0,
    )
    @test_throws ArgumentError apply_radiation_reaction!(
        ps,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = 1.0,
        c = 0.0,
    )
    a = ParticleSet{1,T}(5)
    b = ParticleSet{1,T}(5)
    for i = 1:5
        a.v[1][i] = 0.8
        b.v[1][i] = 0.8
    end
    apply_radiation_reaction!(a, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.05; K = 0.01)
    apply_radiation_reaction!(b, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.05; K = 0.01)
    @test a.v[1] == b.v[1]
    # exponential integrator is unconditionally stable: a large K·dt cools, never reverses/heats
    ps3 = ParticleSet{1,T}(1)
    ps3.v[1][1] = 0.99
    γ0 = 1 / sqrt(1 - 0.99^2)
    apply_radiation_reaction!(ps3, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.5; K = 1.0)   # overshoot regime
    @test ps3.v[1][1] > 0 && 1 / sqrt(1 - ps3.v[1][1]^2) < γ0
    @test_throws ArgumentError apply_radiation_reaction!(
        ps3,
        (NaN, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = 0.01,
    )
end

@testset "ANT-001 antenna injects −dt·∇×E_ant, divergence-free" begin
    T = Float64
    N = (16, 12)
    L = (2π, 2π)
    g = FourierGrid(N, L)
    xs = [(i - 1) * L[1] / N[1] for i = 1:N[1]]
    ys = [(j - 1) * L[2] / N[2] for j = 1:N[2]]
    # a smooth structured antenna E-field (all components), exercises every curl term
    Eant = (
        [sin(2xs[i]) * cos(ys[j]) for i = 1:N[1], j = 1:N[2]],
        [cos(xs[i]) * sin(3ys[j]) for i = 1:N[1], j = 1:N[2]],
        [sin(xs[i] + 2ys[j]) for i = 1:N[1], j = 1:N[2]],
    )
    B = ntuple(_ -> zeros(T, N), 3)
    dt = 0.1
    apply_antenna!(B, Eant, dt, g)
    # matches −dt·∇×E_ant
    ce = ntuple(_ -> zeros(T, N), 3)
    curl!(ce, Eant, g)
    for c = 1:3
        @test maximum(abs.(B[c] .- (-dt .* ce[c]))) < 1e-12
    end
    # divergence-free: ∇·B ≈ 0 (a curl has zero divergence)
    dv = zeros(T, N)
    divergence!(dv, B, g)
    @test maximum(abs.(dv)) < 1e-12
end

@testset "ANT-002 antenna: no-op (dt=0) & validation" begin
    T = Float64
    N = (8, 8)
    g = FourierGrid(N, (2π, 2π))
    Eant = ntuple(_ -> [sin(i + j) for i = 1:N[1], j = 1:N[2]] .* 1.0, 3)
    B = ntuple(_ -> ones(T, N), 3)
    snap = map(copy, B)
    apply_antenna!(B, Eant, 0.0, g)                     # dt=0 ⇒ untouched
    @test all(B[c] == snap[c] for c = 1:3)
    @test_throws ArgumentError apply_antenna!(B, Eant, -0.1, g)
    Bbad = ntuple(_ -> zeros(T, (N[1] + 1, N[2])), 3)   # B not matching the grid ⇒ guarded
    @test_throws DimensionMismatch apply_antenna!(Bbad, Eant, 0.1, g)
end
