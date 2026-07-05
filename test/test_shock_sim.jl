# SHK-002 — reflecting-wall perpendicular collisionless shock. A supercritical
# (M_A=3) shock forms; the EOS-independent checks (magnetic-flux freezing and
# mass conservation) are tight; the fluid Rankine–Hugoniot compression is a
# looser comparison (kinetic shocks compress less than a γ=5/3 fluid — the
# downstream is hotter with reflected ions).

using HybridPlasmaPIC, Test, Random, Statistics

front(sh) = sh.x[argmax(abs.(diff(sh.Bz)))]   # steepest Bz gradient

@testset "step_shock! validates timestep and subcycles before mutation" begin
    T = Float64
    sh = PerpShock(8, T(1.0))
    ps = ParticleSet{1,T}(2)
    ps.x[1] .= (T(0.25), T(0.75))
    ps.v[1] .= (T(-0.1), T(-0.2))
    ps.v[2] .= (T(0.05), T(0.06))
    ps.v[3] .= (T(0.0), T(0.0))
    ps.weight .= shock_density_weight(T(1.0), T(1.0), nparticles(ps))
    init_shock!(sh, ps)

    x0 = copy(ps.x[1])
    v0 = ntuple(c -> copy(ps.v[c]), 3)
    B0 = copy(sh.Bz)
    n0 = copy(sh.n)

    @test_throws ArgumentError step_shock!(sh, ps, 0.1; NB = 0)
    @test_throws ArgumentError step_shock!(sh, ps, NaN; NB = 1)
    @test_throws ArgumentError step_shock!(sh, ps, -0.1; NB = 1)
    @test ps.x[1] == x0
    @test all(ps.v[c] == v0[c] for c = 1:3)
    @test sh.Bz == B0
    @test sh.n == n0
end

@testset "PerpShock validates physical parameters" begin
    @test_throws ArgumentError PerpShock(8, 1.0; Te = -0.1)
    @test_throws ArgumentError PerpShock(8, 1.0; γe = 1.0)
    @test_throws ArgumentError PerpShock(8, 1.0; η = NaN)
    @test_throws ArgumentError PerpShock(8, 1.0; η = -0.1)
    @test_throws ArgumentError PerpShock(8, 1.0; nfloor = 0.0)
    @test_throws ArgumentError shock_density_weight(1.0, 1.0, 0)
end

@testset "PerpShock moment velocity floor uses density units" begin
    T = Float64
    sh = PerpShock(3, T(1.0); nfloor = T(0.2))
    ps = ParticleSet{1,T}(1)
    ps.x[1][1] = T(0.5)
    ps.v[1][1] = T(1.0)
    ps.weight[1] = T(0.15)
    deposit_moments!(sh, ps)
    @test sh.n[2] > sh.nfloor
    @test sh.ux[2] ≈ 1.0
end

@testset "SHK-002 reflecting-wall perpendicular shock" begin
    T = Float64
    N = 512
    Lx = 120.0
    B0 = 1.0
    U0 = 3.0      # M_A = U0/v_A = 3 (supercritical)
    Te = 0.125
    γe = 5 / 3
    vthi = 0.35
    η = 0.02
    sh = PerpShock(N, Lx; Te, γe, η, τ = U0, B0)
    nppc = 64
    Np = nppc * N
    ps = ParticleSet{1,T}(Np)
    rng = MersenneTwister(1)
    for p = 1:Np
        ps.x[1][p] = Lx * rand(rng)
        ps.v[1][p] = -U0 + vthi * randn(rng)
        ps.v[2][p] = vthi * randn(rng)
        ps.v[3][p] = vthi * randn(rng)
    end
    ps.weight .= shock_density_weight(1.0, Lx, Np)
    init_shock!(sh, ps)

    dt = 0.02
    nsteps = 900
    pos = Float64[]
    tt = Float64[]
    for st = 1:nsteps
        step_shock!(sh, ps, dt; NB = 2)
        if st % 100 == 0
            push!(pos, front(sh))
            push!(tt, st * dt)
        end
    end
    @test all(isfinite, sh.Bz)

    xf = front(sh)
    dmask = (sh.x .> 5.0) .& (sh.x .< xf - 5.0)
    n2 = mean(sh.n[dmask])
    Bz2 = mean(sh.Bz[dmask])
    half = length(pos) ÷ 2
    Vs_front = (pos[end] - pos[half]) / (tt[end] - tt[half])
    Vs_mass = U0 / (n2 - 1)

    # a real, supercritical shock formed
    @test 2.0 < n2 < 4.0                          # supercritical, below strong-shock limit 4
    @test xf > 10.0                                # front propagated away from the wall

    # PRIMARY (EOS-independent): magnetic flux frozen to the flow ⇒ Bz2/Bz1 = n2/n1
    @test abs((Bz2 / B0) / n2 - 1) < 0.03
    # PRIMARY: mass conservation (independent front-tracking vs jump relation)
    @test abs(Vs_front - Vs_mass) / Vs_mass < 0.06

    # SECONDARY (kinetic band): compression vs fluid RH at the realized Mach
    M = (U0 + Vs_front) / 1.0
    p1 = vthi^2 + Te
    sol = rankine_hugoniot(MHDState(1.0, U0 + Vs_front, 0.0, p1, 0.0, B0), 5 / 3)
    @test sol.X > n2                               # kinetic compresses less than fluid
    @test abs(n2 - sol.X) / sol.X < 0.25           # within the kinetic overshoot band
end
