# Particle benchmarks: PUSH-001..004 (exact single-particle orbits),
# LOAD-001..002 (loading moments / quiet start), and particle boundaries.
# Oracles are closed-form (analytic orbits, analytic drifts) — never the mover.

using HybridPlasmaPIC, Test, LinearAlgebra, Random, Statistics

struct OffsetParticleVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    first_index::Int
end

Base.size(v::OffsetParticleVector) = size(v.data)
Base.axes(v::OffsetParticleVector) = (v.first_index:(v.first_index+length(v.data)-1),)
Base.IndexStyle(::Type{<:OffsetParticleVector}) = IndexLinear()
Base.getindex(v::OffsetParticleVector, i::Int) = v.data[i-v.first_index+1]
Base.setindex!(v::OffsetParticleVector, x, i::Int) = (v.data[i-v.first_index+1] = x)

# Run a single gyrating particle; return (speed series, accumulated signed phase,
# x series, y series). q=m=1 so Ω_c = B0.
function run_gyro(::Type{T}, dt, nsteps; vperp = 1.0, B0 = 1.0) where {T}
    ps = ParticleSet{2,T}(1)
    ps.x[1][1] = 0
    ps.x[2][1] = 0
    ps.v[1][1] = T(vperp)
    ps.v[2][1] = 0
    ps.v[3][1] = 0
    E = (zero(T), zero(T), zero(T))
    B = (zero(T), zero(T), T(B0))
    xs = Vector{T}(undef, nsteps)
    ys = similar(xs)
    spd = similar(xs)
    phase = 0.0
    prev = atan(ps.v[2][1], ps.v[1][1])
    for n = 1:nsteps
        push_uniform!(ps, E, B, dt)
        vx = ps.v[1][1]
        vy = ps.v[2][1]
        vz = ps.v[3][1]
        cur = atan(vy, vx)
        dphi = cur - prev
        dphi -= 2π * round(dphi / (2π))          # unwrap to (−π,π]
        phase += dphi
        prev = cur
        xs[n] = ps.x[1][1]
        ys[n] = ps.x[2][1]
        spd[n] = sqrt(vx^2 + vy^2 + vz^2)
    end
    return spd, phase, xs, ys
end

@testset "PUSH-001 uniform B: gyration" begin
    T = Float64
    dt = 0.05
    B0 = 1.0
    vperp = 1.0
    Ωc = B0
    nsteps = 2000
    spd, phase, xs, ys = run_gyro(T, dt, nsteps; vperp, B0)
    # constant speed (pure rotation): roundoff
    @test maximum(abs.(spd .- vperp)) < 1e-12
    # gyrofrequency: Boris error is O((Ωdt)²)
    ωmeas = abs(phase) / (nsteps * dt)
    @test abs(ωmeas - Ωc) / Ωc < 1e-3
    # gyroradius from the orbit: ρ = vperp/Ωc, leapfrog radius ρ·√(1+(Ωdt/2)²)
    cx = mean(xs)
    cy = mean(ys)
    r = mean(@. sqrt((xs - cx)^2 + (ys - cy)^2))
    ρ = vperp / Ωc
    @test abs(r - ρ) / ρ < 5e-3
end

@testset "PUSH-002 E×B drift" begin
    T = Float64
    dt = 0.02
    B0 = 1.0
    E0 = 0.1
    ps = ParticleSet{3,T}(1)
    for d = 1:3
        ps.x[d][1] = 0
        ps.v[d][1] = 0
    end
    E = (T(E0), zero(T), zero(T))
    B = (zero(T), zero(T), T(B0))
    nsteps = 6000                                   # ~19 gyroperiods
    sx = 0.0
    sy = 0.0
    sz = 0.0
    for n = 1:nsteps
        push_uniform!(ps, E, B, dt)
        sx += ps.v[1][1]
        sy += ps.v[2][1]
        sz += ps.v[3][1]
    end
    vdx = sx / nsteps
    vdy = sy / nsteps
    vdz = sz / nsteps
    # v_d = E×B/B² = (0, −E0/B0, 0)
    @test abs(vdy + E0 / B0) / (E0 / B0) < 2e-2
    @test abs(vdx) < 2e-3
    @test abs(vdz) < 1e-12
end

@testset "PUSH-003 parallel E acceleration" begin
    T = Float64
    dt = 0.05
    B0 = 1.0
    E0 = 0.3
    v0 = 0.2
    ps = ParticleSet{3,T}(1)
    for d = 1:3
        ps.x[d][1] = 0
        ps.v[d][1] = 0
    end
    ps.v[3][1] = T(v0)                              # initial v∥
    E = (zero(T), zero(T), T(E0))
    B = (zero(T), zero(T), T(B0))
    nsteps = 100
    for _ = 1:nsteps
        push_uniform!(ps, E, B, dt)
    end
    qm = ps.q / ps.m
    expected = v0 + qm * E0 * dt * nsteps           # Boris is exact for E∥B
    @test abs(ps.v[3][1] - expected) < 1e-12
    @test abs(ps.v[1][1]) < 1e-12 && abs(ps.v[2][1]) < 1e-12
end

@testset "PUSH gathered fields match uniform fields" begin
    T = Float64
    N = 5
    dt = T(0.03)
    E0 = (T(0.2), T(-0.1), T(0.05))
    B0 = (T(0.03), T(0.04), T(0.2))
    ps_uniform = ParticleSet{2,T}(N; q = T(1.5), m = T(2.0))
    ps_gathered = ParticleSet{2,T}(N; q = ps_uniform.q, m = ps_uniform.m)
    ps_uniform.x[1] .= T[0.1, 0.2, 0.3, 0.4, 0.5]
    ps_uniform.x[2] .= T[0.6, 0.7, 0.8, 0.9, 1.0]
    ps_uniform.v[1] .= T[0.2, -0.1, 0.4, -0.3, 0.5]
    ps_uniform.v[2] .= T[-0.4, 0.3, -0.2, 0.1, 0.0]
    ps_uniform.v[3] .= T[0.7, -0.6, 0.5, -0.4, 0.3]
    for d = 1:2
        ps_gathered.x[d] .= ps_uniform.x[d]
    end
    for c = 1:3
        ps_gathered.v[c] .= ps_uniform.v[c]
    end
    E = ntuple(c -> fill(E0[c], N), 3)
    B = ntuple(c -> fill(B0[c], N), 3)
    xold = ntuple(d -> copy(ps_gathered.x[d]), 2)
    xmid = ntuple(_ -> zeros(T, N), 2)
    push_uniform!(ps_uniform, E0, B0, dt)
    push_gathered!(ps_gathered, E, B, dt; xmid)
    for c = 1:3
        @test ps_gathered.v[c] ≈ ps_uniform.v[c] rtol = 1e-14 atol = 1e-14
    end
    for d = 1:2
        @test ps_gathered.x[d] ≈ ps_uniform.x[d] rtol = 1e-14 atol = 1e-14
        @test xmid[d] ≈ xold[d] .+ (dt / 2) .* ps_gathered.v[d] rtol = 1e-14 atol = 1e-14
    end
    mixed = ParticleSet{1,T}(2)
    E_mixed = (T[0.1, 0.2], view(T[0.3, 0.4, 0.5], 1:2), T[0.0, 0.0])
    B_mixed = (T[0.0, 0.0], T[0.0, 0.0], T[0.2, 0.2])
    @test push_gathered!(mixed, E_mixed, B_mixed, dt) === mixed
    @test_throws DimensionMismatch push_gathered!(ps_gathered, (E[1][1:4], E[2], E[3]), B, dt)
end

@testset "pushers reject invalid inputs" begin
    T = Float64

    function seeded_particle(; q = 1.0, m = 1.0)
        ps = ParticleSet{1,T}(1; q, m)
        ps.x[1][1] = 0.0
        ps.v[1][1] = 1.0
        ps.v[2][1] = 0.0
        ps.v[3][1] = 0.0
        return ps
    end

    @test_throws ArgumentError push_uniform!(
        seeded_particle(),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        NaN,
    )
    @test_throws ArgumentError push_uniform!(
        seeded_particle(),
        (NaN, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1,
    )
    @test_throws ArgumentError push_uniform!(
        seeded_particle(),
        (0.0, 0.0, 0.0),
        (0.0, Inf, 1.0),
        0.1,
    )
    @test_throws ArgumentError push_uniform!(
        seeded_particle(; q = NaN),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1,
    )
    @test_throws ArgumentError push_uniform!(
        seeded_particle(; m = 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1,
    )

    @test_throws ArgumentError push_gathered!(
        seeded_particle(),
        ([NaN], [0.0], [0.0]),
        ([0.0], [0.0], [1.0]),
        0.1,
    )
    @test_throws ArgumentError push_gathered!(
        seeded_particle(),
        ([0.0], [0.0], [0.0]),
        ([0.0], [0.0], [1.0]),
        NaN,
    )
    @test_throws ArgumentError push_gathered!(
        seeded_particle(; q = NaN),
        ([0.0], [0.0], [0.0]),
        ([0.0], [0.0], [1.0]),
        0.1,
    )
    @test_throws ArgumentError push_gathered!(
        seeded_particle(; m = 0.0),
        ([0.0], [0.0], [0.0]),
        ([0.0], [0.0], [1.0]),
        0.1,
    )
end

@testset "PUSH-004 timestep convergence order" begin
    T = Float64
    Ωc = 1.0
    Tphys = 20.0
    dts = [0.2, 0.1, 0.05, 0.025]
    errs = Float64[]
    for dt in dts
        nsteps = round(Int, Tphys / dt)
        _, phase, _, _ = run_gyro(T, dt, nsteps; vperp = 1.0, B0 = Ωc)
        elapsed = nsteps * dt
        push!(errs, abs(abs(phase) - Ωc * elapsed))
    end
    # slope of log(err) vs log(dt)
    p = (log(errs[end]) - log(errs[1])) / (log(dts[end]) - log(dts[1]))
    @test 1.8 < p < 2.2
end

@testset "LOAD-001 Maxwellian moments" begin
    T = Float64
    N = 200_000
    rng = MersenneTwister(20240625)
    u0 = (0.3, -0.2, 0.1)
    vth = (1.0, 1.5, 0.7)
    ps = ParticleSet{1,T}(N)
    load_maxwellian!(ps, rng, u0, vth)
    z = 6.0                                         # 6σ CI → effectively never flaky
    for c = 1:3
        m̂ = mean(ps.v[c])
        v̂ = var(ps.v[c])
        @test abs(m̂ - u0[c]) < z * vth[c] / sqrt(N)          # mean CI
        @test abs(v̂ - vth[c]^2) < z * vth[c]^2 * sqrt(2 / N) # variance CI
    end
end

@testset "LOAD-002 quiet start" begin
    T = Float64
    N = 4096
    rng = MersenneTwister(7)
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, 1.0)
    load_quiet_velocities!(ps, rng, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    # exact thermal-momentum cancellation
    for c = 1:3
        @test abs(sum(ps.v[c])) < 1e-10
    end
    # low-k density noise: ρ̂(k₁) for the lattice is ~0; random load ~ √N
    k1 = 2π / 1.0
    ρ_quiet = abs(sum(p -> cis(-k1 * ps.x[1][p]), 1:N))
    psr = ParticleSet{1,T}(N)
    load_uniform!(psr, MersenneTwister(99), (0.0,), (1.0,))
    ρ_rand = abs(sum(p -> cis(-k1 * psr.x[1][p]), 1:N))
    @test ρ_quiet < ρ_rand / 50
    @test ρ_quiet < 1e-8
end

@testset "particle loaders reject invalid inputs" begin
    rng = MersenneTwister(11)
    g = FourierGrid((4,), (4.0,))
    ps = ParticleSet{1,Float64}(2)

    @test_throws ArgumentError load_uniform!(ps, rng, (0.0,), (NaN,))
    @test_throws ArgumentError load_uniform!(ps, rng, (1.0,), (0.0,))
    @test_throws ArgumentError load_lattice!(ps, (0.0,), (NaN,), (2,))
    @test_throws ArgumentError load_lattice!(ps, (1.0,), (0.0,), (2,))

    @test_throws ArgumentError set_density_weight!(ps, NaN, g)
    @test_throws ArgumentError set_density_weight!(ps, -1.0, g)

    @test_throws ArgumentError load_maxwellian!(ps, rng, (NaN, 0.0, 0.0), (1.0, 1.0, 1.0))
    @test_throws ArgumentError load_maxwellian!(ps, rng, (0.0, 0.0, 0.0), (-1.0, 1.0, 1.0))

    @test_throws ArgumentError load_quiet_velocities!(ps, rng, (NaN, 0.0, 0.0), (1.0, 1.0, 1.0))
    @test_throws ArgumentError load_quiet_velocities!(ps, rng, (0.0, 0.0, 0.0), (Inf, 1.0, 1.0))
end

@testset "particle boundaries" begin
    T = Float64
    # periodic wrap into [0,1)
    ps = ParticleSet{1,T}(4)
    ps.x[1] .= T[-0.3, 1.2, 0.5, 2.7]
    apply_periodic!(ps, (0.0,), (1.0,))
    @test all(0 .<= ps.x[1] .< 1)
    @test ps.x[1] ≈ [0.7, 0.2, 0.5, 0.7]
    @test nparticles(ps) == 4

    # reflecting wall flips normal velocity
    ps = ParticleSet{1,T}(2)
    ps.x[1] .= T[1.2, -0.3]
    ps.v[1] .= T[1.0, -1.0]
    apply_reflecting!(ps, (0.0,), (1.0,))
    @test ps.x[1] ≈ [0.8, 0.3]
    @test ps.v[1] ≈ [-1.0, 1.0]
    @test nparticles(ps) == 2

    # absorbing removes only out-of-box particles
    ps = ParticleSet{1,T}(4)
    ps.x[1] .= T[0.5, -0.1, 1.3, 0.9]
    ps.id .= UInt64[10, 11, 12, 13]
    nrem = apply_absorbing!(ps, (0.0,), (1.0,))
    @test nrem == 2
    @test nparticles(ps) == 2
    @test ps.x[1] ≈ [0.5, 0.9]
    @test ps.id == UInt64[10, 13]
end

@testset "apply_reflecting! keeps a particle off the exact upper wall" begin
    # Regression: x==hi mapped to 2h−h==h (outside the half-open box [lo,hi)).
    ps = ParticleSet{1,Float64}(1)
    ps.x[1][1] = 10.0
    ps.v[1][1] = 1.0
    apply_reflecting!(ps, (0.0,), (10.0,))
    @test ps.x[1][1] < 10.0          # strictly inside; was == 10.0 before the fix
    @test ps.v[1][1] == -1.0
end

@testset "particle boundaries reject invalid bounds" begin
    function boundary_seed()
        ps = ParticleSet{1,Float64}(2)
        ps.x[1] .= [-0.2, 1.2]
        ps.v[1] .= [1.0, -1.0]
        return ps
    end

    @test_throws ArgumentError apply_periodic!(boundary_seed(), (0.0,), (NaN,))
    @test_throws ArgumentError apply_periodic!(boundary_seed(), (1.0,), (0.0,))

    @test_throws ArgumentError apply_reflecting!(boundary_seed(), (0.0,), (NaN,))
    @test_throws ArgumentError apply_reflecting!(boundary_seed(), (1.0,), (0.0,))

    @test_throws ArgumentError apply_absorbing!(boundary_seed(), (0.0,), (NaN,))
    @test_throws ArgumentError apply_absorbing!(boundary_seed(), (1.0,), (0.0,))
end

@testset "ParticleSet custom arrays must use one-based axes" begin
    N = 3
    x_bad = (OffsetParticleVector(zeros(Float64, N), -2),)
    v = ntuple(_ -> zeros(Float64, N), 3)
    weight = ones(Float64, N)
    id = UInt64.(1:N)
    tag = zeros(UInt32, N)
    @test_throws ArgumentError ParticleSet{1,Float64}(x_bad, v, weight, id, tag, 1.0, 1.0)

    x = (zeros(Float64, N),)
    v_bad = ntuple(_ -> OffsetParticleVector(zeros(Float64, N), 0), 3)
    @test_throws ArgumentError ParticleSet{1,Float64}(x, v_bad, weight, id, tag, 1.0, 1.0)

    weight_bad = OffsetParticleVector(ones(Float64, N), 4)
    @test_throws ArgumentError ParticleSet{1,Float64}(x, v, weight_bad, id, tag, 1.0, 1.0)
end
