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
# Oracles: exact per-pair (hence whole-set) momentum + energy conservation for
# arbitrary weights (the weighted-centre-of-mass rotation), and physical relaxation
# of a temperature anisotropy toward isotropy.

@testset "TA-001 Coulomb momentum & energy conservation (weighted)" begin
    T = Float64
    N = 20_000
    rng = MersenneTwister(2024)
    ps = ParticleSet{1,T}(N)
    load_maxwellian!(ps, rng, (0.4, -0.3, 0.2), (1.0, 1.6, 0.7))
    for p = 1:N
        ps.weight[p] = 0.5 + rand(rng)            # heterogeneous weights
    end

    P0, E0 = set_totals(ps)
    crng = MersenneTwister(99)
    for _ = 1:20
        collide_coulomb!(ps, 0.8, 0.05; rng = crng)
    end
    P1, E1 = set_totals(ps)

    scale = sqrt(P0[1]^2 + P0[2]^2 + P0[3]^2) + 1.0
    @test abs(P1[1] - P0[1]) / scale < 1e-10       # weighted-CM update ⇒ exact momentum
    @test abs(P1[2] - P0[2]) / scale < 1e-10
    @test abs(P1[3] - P0[3]) / scale < 1e-10
    @test abs(E1 - E0) / E0 < 1e-10                # pure rotation |g'|=|g| ⇒ exact energy
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
    @test_throws ArgumentError collide_coulomb!(ps, 1.0, 0.1; u_floor = 0.0)
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
