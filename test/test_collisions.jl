# BGK collision operator (§3.3 / §21.5): conservation, isotropization, and
# Maxwellianization. Oracles are the analytic conservation laws (Σ w v, Σ w |v|²)
# and the qualitative relaxation behavior (anisotropy → 0, kurtosis → Maxwellian).

using HybridPlasmaPIC, Test, Random, Statistics

# Whole-set weighted totals: momentum (3-vec) and kinetic energy Σ w |v|².
function set_totals(ps)
    w = ps.weight
    vx, vy, vz = ps.v
    px = sum(w .* vx)
    py = sum(w .* vy)
    pz = sum(w .* vz)
    E = sum(@. w * (vx^2 + vy^2 + vz^2))
    return (px, py, pz), E
end

# Per-component weighted variance about the weighted mean (∝ temperature).
function comp_temp(ps, c)
    w = ps.weight
    vc = ps.v[c]
    W = sum(w)
    u = sum(w .* vc) / W
    return sum(@. w * (vc - u)^2) / W
end

# Excess kurtosis of an (unweighted) velocity component: 0 for a Gaussian.
function excess_kurtosis(v)
    μ = mean(v)
    σ2 = mean((v .- μ) .^ 2)
    m4 = mean((v .- μ) .^ 4)
    return m4 / σ2^2 - 3
end

@testset "BGK-001 momentum & energy conservation" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(2024)
    ps = ParticleSet{1,T}(N)
    # Drifting bi-Maxwellian + nonuniform weights (conservation must be weighted).
    load_maxwellian!(ps, rng, (0.4, -0.3, 0.2), (1.0, 1.6, 0.7))
    for p = 1:N
        ps.weight[p] = 0.5 + rand(rng)            # heterogeneous weights
    end

    P0, E0 = set_totals(ps)
    collide_bgk!(ps, 5.0, 0.1; rng = MersenneTwister(99))   # νdt = 0.5 → ~39% scatter
    P1, E1 = set_totals(ps)

    scale = sqrt(P0[1]^2 + P0[2]^2 + P0[3]^2) + 1.0
    @test abs(P1[1] - P0[1]) / scale < 1e-10
    @test abs(P1[2] - P0[2]) / scale < 1e-10
    @test abs(P1[3] - P0[3]) / scale < 1e-10
    @test abs(E1 - E0) / E0 < 1e-10
end

@testset "BGK-002 conservation with uniform weights & large νdt" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(7)
    ps = ParticleSet{2,T}(N)
    load_maxwellian!(ps, rng, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    P0, E0 = set_totals(ps)
    # Large νdt → almost all particles scatter (Pcoll ≈ 0.993); still exact.
    collide_bgk!(ps, 100.0, 0.05; rng = MersenneTwister(123))
    P1, E1 = set_totals(ps)
    @test abs(P1[1] - P0[1]) < 1e-9 && abs(P1[2] - P0[2]) < 1e-9 && abs(P1[3] - P0[3]) < 1e-9
    @test abs(E1 - E0) / E0 < 1e-10
end

@testset "BGK conserves energy when selected subset is cold" begin
    T = Float64
    ps = ParticleSet{1,T}(3)
    ps.v[1] .= T[0.0, 0.0, 10.0]
    P0, E0 = set_totals(ps)
    collide_bgk!(ps, 1.0, 1.0; rng = MersenneTwister(4))
    P1, E1 = set_totals(ps)
    @test all(P1[c] ≈ P0[c] for c = 1:3)
    @test E1 ≈ E0 rtol = 1e-14 atol = 1e-14
end

@testset "BGK-003 isotropization of a bi-Maxwellian" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(31)
    ps = ParticleSet{1,T}(N)
    # Strongly anisotropic: Tx ≫ Ty,Tz.
    load_maxwellian!(ps, rng, (0.0, 0.0, 0.0), (2.0, 0.5, 0.5))

    function anisotropy(ps)
        Tx = comp_temp(ps, 1)
        Tperp = (comp_temp(ps, 2) + comp_temp(ps, 3)) / 2
        return (Tx - Tperp) / (Tx + Tperp)
    end

    A = Float64[anisotropy(ps)]
    crng = MersenneTwister(555)
    for _ = 1:60
        collide_bgk!(ps, 1.0, 0.2; rng = crng)    # νdt = 0.2 per step
        push!(A, anisotropy(ps))
    end
    # initial anisotropy clearly nonzero; decays monotonically toward 0
    @test A[1] > 0.5
    @test A[end] < 0.05
    @test A[end] < A[1]
    # monotone-ish decrease: every block-average shrinks
    @test mean(A[2:11]) > mean(A[21:30]) > mean(A[end-9:end])
end

@testset "BGK-004 two-beam relaxes toward Maxwellian (kurtosis → 0)" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(101)
    ps = ParticleSet{1,T}(N)
    # Two counter-streaming cold-ish beams in vx → strongly NON-Gaussian (bimodal,
    # negative excess kurtosis). vy,vz thermal.
    half = N ÷ 2
    for p = 1:N
        ps.v[1][p] = (p <= half ? -1.5 : 1.5) + 0.25 * randn(rng)
        ps.v[2][p] = randn(rng)
        ps.v[3][p] = randn(rng)
    end

    k0 = excess_kurtosis(ps.v[1])
    crng = MersenneTwister(777)
    ks = Float64[k0]
    for _ = 1:40
        collide_bgk!(ps, 1.0, 0.25; rng = crng)
        push!(ks, excess_kurtosis(ps.v[1]))
    end
    # bimodal beams start with large negative excess kurtosis; relaxes toward 0
    @test k0 < -0.8
    @test abs(ks[end]) < 0.1
    @test abs(ks[end]) < abs(k0)
end

@testset "BGK-005 no-op edge cases" begin
    T = Float64
    rng = MersenneTwister(1)
    ps = ParticleSet{1,T}(1000)
    load_maxwellian!(ps, rng, (0.1, 0.0, 0.0), (1.0, 1.0, 1.0))
    snap = (copy(ps.v[1]), copy(ps.v[2]), copy(ps.v[3]))
    collide_bgk!(ps, 0.0, 0.1)            # ν = 0 → untouched
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    collide_bgk!(ps, 5.0, 0.0)            # dt = 0 → untouched
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    @test_throws ArgumentError collide_bgk!(ps, -1.0, 0.1)
    @test_throws ArgumentError collide_bgk!(ps, 1.0, -0.1)
    # single particle: no fluctuation to relax → no-op, no error
    ps1 = ParticleSet{1,T}(1)
    ps1.v[1][1] = 0.7
    collide_bgk!(ps1, 10.0, 1.0)
    @test ps1.v[1][1] == 0.7
end

# ---- Takizuka-Abe binary Coulomb collisions (collide_coulomb!) ----------------
# Oracles: exact per-pair (hence whole-set) momentum + energy conservation for EQUAL
# weights (pure rotation about the pair midpoint); conservation in EXPECTATION for
# unequal weights (Higginson et al. 2020 rejection scheme); relaxation rate independent
# of the macro-weight distribution; isotropization monotone in collisionality (isotropic
# large-angle fallback); odd-N triplet so every particle collides at the nominal rate.

@testset "TA-001 Coulomb momentum & energy conservation (equal weights: exact)" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(2024)
    ps = ParticleSet{1,T}(N)
    load_maxwellian!(ps, rng, (0.4, -0.3, 0.2), (1.0, 1.6, 0.7))
    fill!(ps.weight, 1.0)                         # equal weights ⇒ both partners always
    # accept their update ⇒ each pair is a pure midpoint rotation ⇒ exact conservation
    P0, E0 = set_totals(ps)
    crng = MersenneTwister(99)
    for _ = 1:20
        collide_coulomb!(ps, 0.8, 0.05; rng = crng)
    end
    P1, E1 = set_totals(ps)

    scale = sqrt(P0[1]^2 + P0[2]^2 + P0[3]^2) + 1.0
    @test abs(P1[1] - P0[1]) / scale < 1e-10       # midpoint update ⇒ exact momentum
    @test abs(P1[2] - P0[2]) / scale < 1e-10
    @test abs(P1[3] - P0[3]) / scale < 1e-10
    @test abs(E1 - E0) / E0 < 1e-10                # pure rotation |g'|=|g| ⇒ exact energy
end

@testset "TA-001b Coulomb conservation in expectation (unequal weights)" begin
    # With unequal weights the rejection scheme (Higginson, Holod & Link 2020) conserves
    # Σwv and Σw|v|² in EXPECTATION only; each run's drift is a zero-mean random walk.
    # Measured over the 8 pinned collision seeds below (deterministic): per-run dE/E0
    # std = 4.4e-4, max |dE/E0| = 9.3e-4, max component |dP|/scale = 3.8e-3;
    # mean dE/E0 = -1.7e-4 (1.1 SE from zero), mean dPx/scale = +2.1e-4 (0.6 SE).
    # Bounds: per-run ≥ 3.2x the measured max; means ≈ 3.6-3.9x the measured SE.
    T = Float64
    N = 20_000
    dEs = Float64[]
    dPxs = Float64[]
    dPmaxs = Float64[]
    for seed = 1:8
        ps = ParticleSet{1,T}(N)
        rng = MersenneTwister(2024)
        load_maxwellian!(ps, rng, (0.4, -0.3, 0.2), (1.0, 1.6, 0.7))
        for p = 1:N
            ps.weight[p] = 0.5 + rand(rng)        # heterogeneous weights
        end
        P0, E0 = set_totals(ps)
        crng = MersenneTwister(1000 + seed)
        for _ = 1:20
            collide_coulomb!(ps, 0.8, 0.05; rng = crng)
        end
        P1, E1 = set_totals(ps)
        scale = sqrt(P0[1]^2 + P0[2]^2 + P0[3]^2) + 1.0
        push!(dEs, (E1 - E0) / E0)
        push!(dPxs, (P1[1] - P0[1]) / scale)
        push!(dPmaxs, maximum(abs.(P1 .- P0)) / scale)
    end
    @test maximum(abs.(dEs)) < 3e-3                # per-run energy drift bounded
    @test maximum(dPmaxs) < 1.2e-2                 # per-run momentum drift bounded
    @test abs(mean(dEs)) < 6e-4                    # conserved in expectation (energy)
    @test abs(mean(dPxs)) < 1.2e-3                 # conserved in expectation (momentum)
end

@testset "TA-002 Coulomb isotropization of a bi-Maxwellian" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(31)
    ps = ParticleSet{1,T}(N)
    load_maxwellian!(ps, rng, (0.0, 0.0, 0.0), (2.0, 0.5, 0.5))   # Tx ≫ Ty,Tz
    fill!(ps.weight, 1.0)

    aniso(ps) = (Tx = comp_temp(ps, 1);
    Tperp = (comp_temp(ps, 2) + comp_temp(ps, 3)) / 2;
    (Tx - Tperp) / (Tx + Tperp))
    _, E0 = set_totals(ps)
    A = Float64[aniso(ps)]
    crng = MersenneTwister(555)
    for _ = 1:150
        collide_coulomb!(ps, 10.0, 0.05; rng = crng)
        push!(A, aniso(ps))
    end
    _, E1 = set_totals(ps)
    # anisotropy relaxes toward isotropy (speed-dependent Coulomb rate)
    @test A[1] > 0.5
    @test A[end] < 0.1
    @test mean(A[2:16]) > mean(A[end-14:end])      # decreasing
    @test abs(E1 - E0) / E0 < 1e-10                # energy conserved throughout relaxation
end

@testset "TA-003 Coulomb edge cases, validation & determinism" begin
    T = Float64
    rng = MersenneTwister(1)
    ps = ParticleSet{1,T}(1000)
    load_maxwellian!(ps, rng, (0.1, 0.0, 0.0), (1.0, 1.0, 1.0))
    fill!(ps.weight, 1.0)
    snap = (copy(ps.v[1]), copy(ps.v[2]), copy(ps.v[3]))
    collide_coulomb!(ps, 0.0, 0.1)                 # gcoeff = 0 → untouched
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    collide_coulomb!(ps, 5.0, 0.0)                 # dt = 0 → untouched
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    @test_throws ArgumentError collide_coulomb!(ps, -1.0, 0.1)
    @test_throws ArgumentError collide_coulomb!(ps, 1.0, -0.1)
    @test_throws ArgumentError collide_coulomb!(ps, Inf, 0.1)
    @test_throws ArgumentError collide_coulomb!(ps, 1.0, Inf)
    @test_throws ArgumentError collide_coulomb!(ps, 1.0, 0.1; u_floor = 0.0)
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    # single particle: no pair → no-op, no error
    ps1 = ParticleSet{1,T}(1)
    ps1.v[1][1] = 0.7
    collide_coulomb!(ps1, 10.0, 1.0)
    @test ps1.v[1][1] == 0.7
    # deterministic for a fixed rng
    pa = ParticleSet{1,T}(2000)
    load_maxwellian!(pa, MersenneTwister(5), (0.0, 0.0, 0.0), (1.3, 0.6, 0.6))
    fill!(pa.weight, 1.0)
    pb = ParticleSet{1,T}(2000)
    load_maxwellian!(pb, MersenneTwister(5), (0.0, 0.0, 0.0), (1.3, 0.6, 0.6))
    fill!(pb.weight, 1.0)
    for _ = 1:5
        collide_coulomb!(pa, 2.0, 0.05; rng = MersenneTwister(42))
    end
    for _ = 1:5
        collide_coulomb!(pb, 2.0, 0.05; rng = MersenneTwister(42))
    end
    @test pa.v[1] == pb.v[1] && pa.v[2] == pb.v[2] && pa.v[3] == pb.v[3]
end

@testset "TA-004 Coulomb relaxation rate independent of macro-particle weights" begin
    # Discriminating rate oracle: the SAME physical anisotropic Maxwellian loaded twice —
    # uniform weights vs velocity-independent alternating weights {0.1, 1.9} — must relax
    # its weighted anisotropy at the same physical rate. Pre-fix (weight-weighted-COM
    # kinematics, no variance correction) the endpoints differed by ≈ 0.18 (log-decay
    # ratio ≈ 1.7-1.8); post-fix, measured across 4 seeds: |A60u - A60m| = 0.014-0.034
    # (finite-N estimator noise, dt-independent), log-decay ratio = 1.04-1.11.
    T = Float64
    N = 20_000
    waniso(ps) = (Tx = comp_temp(ps, 1);
    Tp = (comp_temp(ps, 2) + comp_temp(ps, 3)) / 2;
    (Tx - Tp) / (Tx + Tp))
    for seed in (555, 999)
        psu = ParticleSet{1,T}(N)
        load_maxwellian!(psu, MersenneTwister(31), (0.0, 0.0, 0.0), (2.0, 0.5, 0.5))
        fill!(psu.weight, 1.0)
        psm = ParticleSet{1,T}(N)
        load_maxwellian!(psm, MersenneTwister(31), (0.0, 0.0, 0.0), (2.0, 0.5, 0.5))
        for p = 1:N
            psm.weight[p] = isodd(p) ? 0.1 : 1.9  # velocity-independent mixed weights
        end
        A0u, A0m = waniso(psu), waniso(psm)
        ru = MersenneTwister(seed)
        rm = MersenneTwister(seed)
        for _ = 1:60
            collide_coulomb!(psu, 10.0, 0.05; rng = ru)
            collide_coulomb!(psm, 10.0, 0.05; rng = rm)
        end
        Au, Am = waniso(psu), waniso(psm)
        @test A0u > 0.5 && A0m > 0.5
        @test Au < 0.45 && Am < 0.45              # both actually relaxed (from ≈ 0.88)
        @test abs(Au - Am) < 0.06                 # rates agree (pre-fix diff ≈ 0.18)
        r = log(A0u / Au) / log(A0m / Am)
        @test 0.7 < r < 1.3                       # decay-rate ratio (pre-fix ≈ 1.75)
    end
end

@testset "TA-005 Coulomb isotropization monotone in collisionality" begin
    # Large-angle-fallback oracle: residual anisotropy after 3 identical steps must be
    # non-increasing in gcoeff. Pre-fix (Gaussian tan(Θ/2) with ⟨δ²⟩ ≫ 1, no fallback)
    # it was NON-monotonic and stalled: A3 ≈ 0.33, 1.10, 2.62, 2.97, 3.02 from A0 = 2.98
    # over gcoeff = 1e-2..1e8. Post-fix measured A3 ≈ 0.28-0.30 then 0.17-0.19 flat
    # (isotropic-scatter saturation), seed spread ≈ ±0.02. The slow population
    # (vth = 0.02/0.01, u_floor = 1e-3 default) makes ⟨δ²⟩ > 1 the dominant regime.
    T = Float64
    N = 20_000
    ar(ps) = comp_temp(ps, 1) / ((comp_temp(ps, 2) + comp_temp(ps, 3)) / 2) - 1
    ps0 = ParticleSet{1,T}(N)
    load_maxwellian!(ps0, MersenneTwister(7), (0.0, 0.0, 0.0), (0.02, 0.01, 0.01))
    fill!(ps0.weight, 1.0)
    A0 = ar(ps0)
    A3 = Float64[]
    for gc in (1e-2, 1.0, 1e2, 1e4, 1e8)
        ps = ParticleSet{1,T}(N)
        load_maxwellian!(ps, MersenneTwister(7), (0.0, 0.0, 0.0), (0.02, 0.01, 0.01))
        fill!(ps.weight, 1.0)
        crng = MersenneTwister(11)
        for _ = 1:3
            collide_coulomb!(ps, gc, 0.01; rng = crng)
        end
        push!(A3, ar(ps))
    end
    @test A0 > 2.5
    @test all(A3[k+1] <= A3[k] + 0.02 for k = 1:4) # monotone within MC slack
    @test A3[end] < 0.3                            # saturated fast isotropization
    @test A3[end] <= A3[1]                         # high gcoeff relaxes most (pre-fix: least)
end

@testset "TA-006 odd-N Takizuka-Abe triplet" begin
    T = Float64
    # (a) strong collisions: EVERY particle's velocity changes every call (pre-fix odd N
    # left exactly one particle per call bit-identical). Measured: 0 identical in 200
    # calls at N = 3 and N = 5.
    for N in (3, 5)
        rng = MersenneTwister(42)
        ps = ParticleSet{1,T}(N)
        for p = 1:N
            ps.v[1][p] = randn(rng)
            ps.v[2][p] = randn(rng)
            ps.v[3][p] = randn(rng)
        end
        fill!(ps.weight, 1.0)
        nident = 0
        for _ = 1:100
            snap = (copy(ps.v[1]), copy(ps.v[2]), copy(ps.v[3]))
            collide_coulomb!(ps, 1e4, 1.0; rng = rng)
            for p = 1:N
                if ps.v[1][p] == snap[1][p] && ps.v[2][p] == snap[2][p] && ps.v[3][p] == snap[3][p]
                    nident += 1
                end
            end
        end
        @test nident == 0
    end
    # (b) the halved triplet variances sum to the nominal per-particle small-angle rate:
    # with v1 = 0, v2 = x̂, v3 = ŷ (g12 = g31 = 1, g23 = √2) every permutation pairs all
    # three, so E|Δv_p|² = (gc·dt/2)(1/g_pq + 1/g_pr) + O(⟨δ²⟩²). Particle 1 gives
    # gc·dt — exactly the even-N pair rate at g = 1. Measured/predicted = 0.996-1.01
    # (MC SE ≈ 0.7-1% at M = 20_000): rtol 0.05 ≈ 5σ; triplet/pair measured 0.986.
    gcdt = 2e-4
    M = 20_000
    rng3 = MersenneTwister(2718)
    acc = zeros(3)
    ps3 = ParticleSet{1,T}(3)
    fill!(ps3.weight, 1.0)
    for _ = 1:M
        ps3.v[1] .= T[0.0, 1.0, 0.0]
        ps3.v[2] .= T[0.0, 0.0, 1.0]
        ps3.v[3] .= 0.0
        collide_coulomb!(ps3, gcdt, 1.0; rng = rng3)
        acc[1] += ps3.v[1][1]^2 + ps3.v[2][1]^2 + ps3.v[3][1]^2
        acc[2] += (ps3.v[1][2] - 1.0)^2 + ps3.v[2][2]^2 + ps3.v[3][2]^2
        acc[3] += ps3.v[1][3]^2 + (ps3.v[2][3] - 1.0)^2 + ps3.v[3][3]^2
    end
    acc ./= M
    @test isapprox(acc[1], gcdt; rtol = 0.05)
    @test isapprox(acc[2], gcdt * (1 + 1 / sqrt(2)) / 2; rtol = 0.05)
    @test isapprox(acc[3], gcdt * (1 + 1 / sqrt(2)) / 2; rtol = 0.05)
    rng2 = MersenneTwister(3141)
    acc2 = 0.0
    ps2 = ParticleSet{1,T}(2)
    fill!(ps2.weight, 1.0)
    for _ = 1:M
        ps2.v[1] .= T[0.0, 1.0]
        ps2.v[2] .= 0.0
        ps2.v[3] .= 0.0
        collide_coulomb!(ps2, gcdt, 1.0; rng = rng2)
        acc2 += ps2.v[1][1]^2 + ps2.v[2][1]^2 + ps2.v[3][1]^2
    end
    acc2 /= M
    @test isapprox(acc2, gcdt; rtol = 0.05)        # even-N pair baseline
    @test 0.93 < acc[1] / acc2 < 1.07              # odd-N rate matches even-N
end

# ---- Monte-Carlo neutral collisions (collide_neutral_mcc!) --------------------
# Oracle: the charged population relaxes toward the neutral reservoir — temperature
# → T_n, drift → u_n (full thermalization at equal mass).

@testset "MCC-001 neutral collisions thermalize toward the bath" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(3)
    ps = ParticleSet{1,T}(N)
    load_maxwellian!(ps, rng, (1.0, 0.0, 0.0), (sqrt(2.0), sqrt(2.0), sqrt(2.0)))  # hot T_p≈2, drift
    fill!(ps.weight, 1.0)
    Tof(ps) = (comp_temp(ps, 1) + comp_temp(ps, 2) + comp_temp(ps, 3)) / 3
    driftx(ps) = sum(ps.weight .* ps.v[1]) / sum(ps.weight)

    @test Tof(ps) > 1.5 && driftx(ps) > 0.8
    crng = MersenneTwister(9)
    for _ = 1:400
        collide_neutral_mcc!(ps, 0.1; nσ = 1.0, T_n = 0.2, m_n = 1.0, rng = crng)
    end
    # equal-mass elastic MCC → full thermalization: T → T_n = 0.2, drift → u_n = 0
    @test isapprox(Tof(ps), 0.2; atol = 0.03)
    @test abs(driftx(ps)) < 0.05
end

@testset "MCC-002 neutral collisions: edge cases, validation & determinism" begin
    T = Float64
    rng = MersenneTwister(1)
    ps = ParticleSet{1,T}(1000)
    load_maxwellian!(ps, rng, (0.3, 0.0, 0.0), (1.0, 1.0, 1.0))
    fill!(ps.weight, 1.0)
    snap = (copy(ps.v[1]), copy(ps.v[2]), copy(ps.v[3]))
    collide_neutral_mcc!(ps, 0.1; nσ = 0.0, T_n = 0.3)     # nσ = 0 → untouched
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    collide_neutral_mcc!(ps, 0.0; nσ = 1.0, T_n = 0.3)     # dt = 0 → untouched
    @test ps.v[1] == snap[1] && ps.v[2] == snap[2] && ps.v[3] == snap[3]
    @test_throws ArgumentError collide_neutral_mcc!(ps, 0.1; nσ = -1.0, T_n = 0.3)
    @test_throws ArgumentError collide_neutral_mcc!(ps, -0.1; nσ = 1.0, T_n = 0.3)
    @test_throws ArgumentError collide_neutral_mcc!(ps, 0.1; nσ = 1.0, T_n = -0.3)
    @test_throws ArgumentError collide_neutral_mcc!(ps, 0.1; nσ = 1.0, T_n = 0.3, m_n = 0.0)
    # deterministic for a fixed rng
    pa = ParticleSet{1,T}(2000)
    load_maxwellian!(pa, MersenneTwister(5), (0.0, 0.0, 0.0), (1.2, 1.2, 1.2))
    fill!(pa.weight, 1.0)
    pb = ParticleSet{1,T}(2000)
    load_maxwellian!(pb, MersenneTwister(5), (0.0, 0.0, 0.0), (1.2, 1.2, 1.2))
    fill!(pb.weight, 1.0)
    for _ = 1:4
        collide_neutral_mcc!(pa, 0.1; nσ = 1.0, T_n = 0.3, rng = MersenneTwister(42))
    end
    for _ = 1:4
        collide_neutral_mcc!(pb, 0.1; nσ = 1.0, T_n = 0.3, rng = MersenneTwister(42))
    end
    @test pa.v[1] == pb.v[1] && pa.v[2] == pb.v[2] && pa.v[3] == pb.v[3]
end

# ---- electron-impact ionization (ionize_mcc!) ---------------------------------
# Oracles: pair creation (Ne,Ni each +nb → net-neutral), and — with cold neutrals —
# an EXACT electron-population energy loss of nb·E_iz (each primary cooled by E_iz).

@testset "IZ-001 ionization: pair creation & exact energy cost" begin
    T = Float64
    Ne0 = 10_000
    el = ParticleSet{2,T}(Ne0; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(0; q = 1.0, m = 100.0)
    rng = MersenneTwister(2)
    for p = 1:Ne0
        el.x[1][p] = rand(rng)
        el.x[2][p] = rand(rng)
        el.v[1][p] = 2.0                            # KE = 2.0 each (> E_iz)
    end
    fill!(el.weight, 1.0)
    E_iz = 0.5
    eKE(ps) = 0.5 * ps.m * sum(@. ps.v[1]^2 + ps.v[2]^2 + ps.v[3]^2)
    E0 = eKE(el)
    nb = ionize_mcc!(
        el,
        ions,
        0.1;
        nσ_iz = 3.0,
        E_iz = E_iz,
        T_n = 0.0,
        m_n = 100.0,
        rng = MersenneTwister(7),
    )
    @test nb > 0
    @test nparticles(el) == Ne0 + nb               # a secondary electron per ionization
    @test nparticles(ions) == nb                   # an ion per ionization (net-neutral)
    # cold neutrals ⇒ secondaries born at rest ⇒ electrons lose EXACTLY nb·E_iz
    @test isapprox(E0 - eKE(el), nb * E_iz; rtol = 1e-12)
    # newborns get unique ids (no collision with the incident 1..Ne0 ids)
    @test length(unique(el.id)) == nparticles(el)
    @test length(unique(ions.id)) == nparticles(ions)
end

@testset "IZ-002 ionization: threshold, edge cases, validation & determinism" begin
    T = Float64
    # below threshold (KE = 0.125 < E_iz = 0.5) ⇒ no ionization, no growth
    el = ParticleSet{2,T}(5000; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(0; q = 1.0, m = 100.0)
    rng = MersenneTwister(1)
    for p = 1:5000
        el.x[1][p] = rand(rng)
        el.x[2][p] = rand(rng)
        el.v[1][p] = 0.5
    end
    fill!(el.weight, 1.0)
    @test ionize_mcc!(el, ions, 0.1; nσ_iz = 10.0, E_iz = 0.5, rng = MersenneTwister(1)) == 0
    @test nparticles(el) == 5000 && nparticles(ions) == 0
    # no-op: nσ_iz = 0, dt = 0
    el2 = ParticleSet{2,T}(1000; q = -1.0, m = 1.0)
    i2 = ParticleSet{2,T}(0; q = 1.0, m = 100.0)
    for p = 1:1000
        el2.v[1][p] = 3.0
    end
    fill!(el2.weight, 1.0)
    @test ionize_mcc!(el2, i2, 0.1; nσ_iz = 0.0, E_iz = 0.5) == 0
    @test ionize_mcc!(el2, i2, 0.0; nσ_iz = 5.0, E_iz = 0.5) == 0
    @test nparticles(el2) == 1000 && nparticles(i2) == 0
    # input validation
    @test_throws ArgumentError ionize_mcc!(el2, i2, 0.1; nσ_iz = -1.0, E_iz = 0.5)
    @test_throws ArgumentError ionize_mcc!(el2, i2, 0.1; nσ_iz = 1.0, E_iz = -0.5)
    @test_throws ArgumentError ionize_mcc!(el2, i2, 0.1; nσ_iz = 1.0, E_iz = 0.5, T_n = -1.0)
    @test_throws ArgumentError ionize_mcc!(el2, i2, 0.1; nσ_iz = 1.0, E_iz = 0.5, m_n = 0.0)
    # deterministic for a fixed rng (count + all created velocities)
    ea = ParticleSet{2,T}(3000; q = -1.0, m = 1.0)
    ia = ParticleSet{2,T}(0; q = 1.0, m = 50.0)
    eb = ParticleSet{2,T}(3000; q = -1.0, m = 1.0)
    ib = ParticleSet{2,T}(0; q = 1.0, m = 50.0)
    for p = 1:3000, s in (ea, eb)
        s.x[1][p] = 0.3
        s.x[2][p] = 0.4
        s.v[1][p] = 1.8
    end
    fill!(ea.weight, 1.0)
    fill!(eb.weight, 1.0)
    na = ionize_mcc!(ea, ia, 0.1; nσ_iz = 2.0, E_iz = 0.4, T_n = 0.1, rng = MersenneTwister(9))
    nbb = ionize_mcc!(eb, ib, 0.1; nσ_iz = 2.0, E_iz = 0.4, T_n = 0.1, rng = MersenneTwister(9))
    @test na == nbb && ea.v[1] == eb.v[1] && ia.v[1] == ib.v[1]
end

@testset "IZ-003 threaded id counter ⇒ globally-unique ids (self-heal + no reuse under removal)" begin
    T = Float64
    el = ParticleSet{2,T}(200; q = -1.0, m = 1.0)          # pre-existing ids 1..200
    ions = ParticleSet{2,T}(0; q = 1.0, m = 100.0)
    for p = 1:200
        el.v[1][p] = 2.0
    end
    fill!(el.weight, 1.0)
    ce = Ref(UInt64(1))   # deliberately below the live max: must self-heal above it
    ci = Ref(UInt64(1))
    n1 = ionize_mcc!(
        el,
        ions,
        0.1;
        nσ_iz = 5.0,
        E_iz = 0.5,
        e_nextid = ce,
        i_nextid = ci,
        rng = MersenneTwister(1),
    )
    @test n1 > 0
    @test length(unique(el.id)) == nparticles(el)          # newborns don't collide with 1..200
    @test ce[] > UInt64(200)                                # counter self-healed above the live max
    call1_ids = el.id[end-n1+1:end]
    # counter only increases ⇒ even if the current max-id particle is later removed, the next
    # batch's ids strictly exceed every previously-issued id, so a removed id is never reissued
    n2 = ionize_mcc!(
        el,
        ions,
        0.1;
        nσ_iz = 5.0,
        E_iz = 0.5,
        e_nextid = ce,
        i_nextid = ci,
        rng = MersenneTwister(2),
    )
    call2_ids = el.id[end-n2+1:end]
    @test all(id -> id > maximum(call1_ids), call2_ids)    # monotonic across calls
    @test length(unique(el.id)) == nparticles(el)          # still globally unique

    # a 0-valued counter must NOT underflow (e_nextid[]-1 would wrap to typemax) — it is
    # treated as "start above the live set"
    el0 = ParticleSet{2,T}(10; q = -1.0, m = 1.0)
    ions0 = ParticleSet{2,T}(0; q = 1.0, m = 100.0)
    for p = 1:10
        el0.v[1][p] = 2.0
    end
    fill!(el0.weight, 1.0)
    z = Ref(UInt64(0))
    ionize_mcc!(
        el0,
        ions0,
        0.1;
        nσ_iz = 5.0,
        E_iz = 0.5,
        e_nextid = z,
        i_nextid = Ref(UInt64(0)),
        rng = MersenneTwister(3),
    )
    @test length(unique(el0.id)) == nparticles(el0)        # no wrap/collision with live 1..10
    @test all(>(UInt64(10)), el0.id[11:end])               # newborns above the live max
    @test z[] > UInt64(10)                                  # counter self-healed, no wrap
end
