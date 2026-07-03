# Phase-12b full electromagnetic 1D PIC. Oracles (all closed-form / exact):
#   (1) transverse EM-wave dispersion ω² = ω_pe² + c²k²  (≥2 k values, <3%);
#   (2) Esirkepov charge conservation  max|ρ^{n+1}−ρ^n + dt ∂xJx|/scale < 1e-10;
#   (3) total (field+kinetic) energy bounded over the run.

using HybridPlasmaPIC, FFTW, Test, Random

@testset "EMPIC1D parameter validation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    es = EMPIC1D(g, 4; n0 = 1.0, c = 5.0, mi = 1836.0)
    @test es.n0 == 1.0
    @test es.c == 5.0
    @test es.mi == 1836.0
    @test_throws ArgumentError EMPIC1D(g, 4; n0 = NaN)
    @test_throws ArgumentError EMPIC1D(g, 4; n0 = -1.0)
    @test_throws ArgumentError EMPIC1D(g, 4; c = NaN)
    @test_throws ArgumentError EMPIC1D(g, 4; c = 0.0)
    @test_throws ArgumentError EMPIC1D(g, 4; c = -1.0)
    @test_throws ArgumentError EMPIC1D(g, 4; mi = NaN)
    @test_throws ArgumentError EMPIC1D(g, 4; mi = 0.0)
end

@testset "EMPIC1D electron species validation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))

    badq = ParticleSet{1,T}(8; q = 0.0, m = 1.0)
    load_lattice_1d!(badq, 0.0, 2π)
    set_density_weight!(badq, 1.0, g)
    @test_throws ArgumentError init_empic!(EMPIC1D(g, 8; n0 = 1.0, c = 5.0, shape = CIC()), badq)

    badm = ParticleSet{1,T}(8; q = -1.0, m = 2.0)
    load_lattice_1d!(badm, 0.0, 2π)
    set_density_weight!(badm, 1.0, g)
    @test_throws ArgumentError init_empic!(EMPIC1D(g, 8; n0 = 1.0, c = 5.0, shape = CIC()), badm)

    e = ParticleSet{1,T}(8; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, 2π)
    set_density_weight!(e, 1.0, g)
    es = EMPIC1D(g, 8; n0 = 1.0, c = 5.0, shape = CIC())
    init_empic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 1)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    Ex0 = copy(es.Ex)
    Ey0 = copy(es.Ey)
    Bz0 = copy(es.Bz)
    rho_n0 = copy(es.rho_n)
    rho_np10 = copy(es.rho_np1)
    time0 = es.time[]
    step0 = es.step[]

    e.q = 0.0
    @test_throws ArgumentError step_empic!(es, e, 0.1)
    @test all(e.x[d] == x0[d] for d = 1:1)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test es.Ex == Ex0
    @test es.Ey == Ey0
    @test es.Bz == Bz0
    @test es.rho_n == rho_n0
    @test es.rho_np1 == rho_np10
    @test es.time[] == time0
    @test es.step[] == step0

    e.q = -1.0
    e.m = 2.0
    @test_throws ArgumentError step_empic!(es, e, 0.1)
    @test all(e.x[d] == x0[d] for d = 1:1)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test es.Ex == Ex0
    @test es.Ey == Ey0
    @test es.Bz == Bz0
    @test es.rho_n == rho_n0
    @test es.rho_np1 == rho_np10
    @test es.time[] == time0
    @test es.step[] == step0
end

@testset "step_empic! validates timestep before mutation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    e = ParticleSet{1,T}(8; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, 2π)
    set_density_weight!(e, 1.0, g)
    es = EMPIC1D(g, 8; n0 = 1.0, c = 5.0, shape = CIC())
    init_empic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 1)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    Ex0 = copy(es.Ex)
    Ey0 = copy(es.Ey)
    Bz0 = copy(es.Bz)
    rho_n0 = copy(es.rho_n)
    rho_np10 = copy(es.rho_np1)
    time0 = es.time[]
    step0 = es.step[]

    @test_throws ArgumentError step_empic!(es, e, NaN)
    @test_throws ArgumentError step_empic!(es, e, -0.1)
    @test all(e.x[d] == x0[d] for d = 1:1)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test es.Ex == Ex0
    @test es.Ey == Ey0
    @test es.Bz == Bz0
    @test es.rho_n == rho_n0
    @test es.rho_np1 == rho_np10
    @test es.time[] == time0
    @test es.step[] == step0
end

# dominant positive-frequency peak of a complex time series (parabolic interp)
function _em_peakfreq(series, dt)
    Nt = length(series)
    S = fft(series)
    mag = abs.(S)
    best = 2
    bm = -1.0
    for j = 2:Nt÷2
        mag[j+1] > bm && (bm = mag[j+1]; best = j)
    end
    a = mag[best]
    b = mag[best+1]
    c = mag[best+2]
    d = a - 2b + c
    δ = d != 0 ? 0.5 * (a - c) / d : 0.0
    return 2π * (best + δ) / (Nt * dt)
end

# Seed a small transverse plane wave of integer mode m on Ey/Bz consistent with
# ∂tBz = −∂xEy for a wave Ey = A cos(kx): the linearized dispersion gives Bz so
# both branches are excited cleanly. We simply seed Ey and let Bz grow self-
# consistently; a small amplitude keeps it linear.
@testset "Phase-12b EM PIC: transverse dispersion ω²=ω_pe²+c²k²" begin
    T = Float64
    n = 64
    L = 2π
    cc = 5.0
    n0 = 1.0
    ωpe = sqrt(n0)
    g = FourierGrid((n,), (L,))
    nppc = 200
    N = nppc * n

    measured = Tuple{Float64,Float64}[]   # (ω_measured, ω_theory)
    for m in (1, 2)
        e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
        load_lattice_1d!(e, 0.0, L)
        set_density_weight!(e, n0, g)
        # cold electrons (no thermal spread) so the dispersion is sharp
        es = EMPIC1D(g, N; n0 = n0, c = cc, shape = CIC())
        # seed Ey = A cos(k x), k = 2π m / L
        k = 2π * m / L
        A = 1e-3
        for i = 1:n
            x = (i - 1) * g.dx[1]
            es.Ey[i] = A * cos(k * x)
        end
        init_empic!(es, e)
        dt = 0.01
        Nt = 4000
        series = ComplexF64[]
        for _ = 1:Nt
            step_empic!(es, e, dt)
            push!(series, mode_amplitude(es.Ey, g, (m,)))
        end
        ω = _em_peakfreq(series, dt)
        ωth = sqrt(ωpe^2 + cc^2 * k^2)
        push!(measured, (ω, ωth))
        @test abs(ω - ωth) / ωth < 0.03
    end
    @info "EM dispersion (ω_measured vs ω_theory)" measured
end

@testset "Phase-12b EM PIC: Esirkepov charge conservation" begin
    T = Float64
    n = 48
    L = 2π
    g = FourierGrid((n,), (L,))
    nppc = 100
    N = nppc * n
    rng = MersenneTwister(7)
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_uniform!(e, rng, (0.0,), (L,))
    set_density_weight!(e, 1.0, g)
    # give the electrons a real spread of velocities so Jx is nontrivial
    for p = 1:N
        e.v[1][p] = 0.5 * randn(rng)
        e.v[2][p] = 0.3 * randn(rng)
        e.v[3][p] = 0.2 * randn(rng)
    end
    es = EMPIC1D(g, N; n0 = 1.0, c = 5.0, shape = CIC())
    # seed some transverse field too, to exercise the full step
    for i = 1:n
        es.Ey[i] = 0.01 * sin(2π * (i - 1) * g.dx[1] / L)
    end
    init_empic!(es, e)
    dt = 0.02
    maxres = 0.0
    for _ = 1:200
        step_empic!(es, e, dt)
        maxres = max(maxres, charge_conservation_residual(es, dt))
    end
    @info "charge-conservation residual (max over run)" maxres
    @test maxres < 1e-10
end

@testset "Phase-12b EM PIC: energy bounded" begin
    T = Float64
    n = 64
    L = 2π
    cc = 5.0
    g = FourierGrid((n,), (L,))
    nppc = 200
    N = nppc * n
    rng = MersenneTwister(3)
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, 1.0, g)
    for p = 1:N
        e.v[1][p] = 0.05 * randn(rng)
    end
    es = EMPIC1D(g, N; n0 = 1.0, c = cc, shape = CIC())
    for i = 1:n
        es.Ey[i] = 1e-3 * cos(2π * (i - 1) * g.dx[1] / L)
    end
    init_empic!(es, e)
    dt = 0.01
    W0 = em_field_energy(es) + kinetic_energy(e)
    Wmax = W0
    Wmin = W0
    for _ = 1:2000
        step_empic!(es, e, dt)
        W = em_field_energy(es) + kinetic_energy(e)
        Wmax = max(Wmax, W)
        Wmin = min(Wmin, W)
    end
    @info "energy bounds" W0 Wmin Wmax
    @test (Wmax - Wmin) / W0 < 0.10
end

# ---- helper: run the transverse dispersion for one (mode m) and return ω -----
# `n_sub` and `mobile` exercise the new full electron–ion PIC paths.
function _em_dispersion_omega(
    m;
    n = 64,
    L = 2π,
    cc = 5.0,
    n0 = 1.0,
    nppc = 200,
    dt = 0.01,
    Nt = 4000,
    mobile = false,
    mi = 1836.0,
    n_sub = 1,
)
    T = Float64
    g = FourierGrid((n,), (L,))
    N = nppc * n
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, n0, g)
    ions = nothing
    if mobile
        ions = ParticleSet{1,T}(N; q = 1.0, m = mi)
        load_lattice_1d!(ions, 0.0, L)
        set_density_weight!(ions, n0, g)
    end
    es = EMPIC1D(g, N; n0 = n0, c = cc, shape = CIC(), mobile = mobile, mi = mi, n_sub = n_sub)
    k = 2π * m / L
    A = 1e-3
    for i = 1:n
        x = (i - 1) * g.dx[1]
        es.Ey[i] = A * cos(k * x)
    end
    init_empic!(es, e, ions)
    series = ComplexF64[]
    for _ = 1:Nt
        step_empic!(es, e, ions, dt)
        push!(series, mode_amplitude(es.Ey, g, (m,)))
    end
    ω = _em_peakfreq(series, dt)
    ωth = sqrt(n0 + cc^2 * k^2)
    return ω, ωth
end

@testset "Full e–i PIC: mobile ions leave high-f EM dispersion unchanged" begin
    # ions (mi=1836) are far too heavy to respond at the EM wave frequency, so the
    # measured ω must still track ω²=ω_pe²+c²k² to within 3%.
    meas = Tuple{Float64,Float64}[]
    for m in (1, 2)
        ω, ωth = _em_dispersion_omega(m; mobile = true, mi = 1836.0)
        push!(meas, (ω, ωth))
        @test abs(ω - ωth) / ωth < 0.03
    end
    @info "mobile-ion EM dispersion (ω_measured vs ω_theory)" meas
end

@testset "Full e–i PIC: electron subcycling convergence (n_sub=1 vs 2)" begin
    agree = Float64[]
    for m in (1, 2)
        ω1, _ = _em_dispersion_omega(m; mobile = true, n_sub = 1)
        # n_sub=2 ⇒ ion step 2·dt, electron substep dt; halve Nt to keep wall-time
        ω2, _ = _em_dispersion_omega(m; mobile = true, n_sub = 2, dt = 0.02, Nt = 2000)
        rel = abs(ω1 - ω2) / ω1
        push!(agree, rel)
        @test rel < 0.05
    end
    @info "subcycle agreement |ω(n_sub=1)−ω(n_sub=2)|/ω(n_sub=1)" agree
end

@testset "Full e–i PIC: n_sub≥3 ion Jy independent of worki init (no stale read)" begin
    # For n_sub ≥ 3 the substeps before the ion push read es.worki before it is
    # written; the step must seed it from the ions' current positions so the
    # transverse-Jy deposit never depends on uninitialized/stale memory. Two runs
    # whose only difference is a garbage worki seed must produce identical fields.
    T = Float64
    n = 32
    L = 2π
    cc = 5.0
    n0 = 1.0
    g = FourierGrid((n,), (L,))
    nppc = 40
    N = nppc * n
    mmode = 1
    function run_poisoned(seed_val)
        e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
        load_lattice_1d!(e, 0.0, L)
        set_density_weight!(e, n0, g)
        ions = ParticleSet{1,T}(N; q = 1.0, m = 1836.0)
        load_lattice_1d!(ions, 0.0, L)
        set_density_weight!(ions, n0, g)
        ions.v[2] .= 0.1                        # transverse drift ⇒ nonzero ion Jy
        es = EMPIC1D(g, N; n0 = n0, c = cc, shape = CIC(), mobile = true, mi = 1836.0, n_sub = 3)
        k = 2π * mmode / L
        for i = 1:n
            es.Ey[i] = 1e-3 * cos(k * (i - 1) * g.dx[1])
        end
        init_empic!(es, e, ions)
        fill!(es.worki, seed_val)               # garbage: must NOT influence output
        step_empic!(es, e, ions, 0.01)
        return copy(es.Ey)
    end
    Ey_a = run_poisoned(0.0)
    Ey_b = run_poisoned(L / 2)
    Ey_c = run_poisoned(1.0e9)
    @test Ey_a == Ey_b
    @test Ey_a == Ey_c
end

@testset "Full e–i PIC: numerical-Cherenkov beam energy bounded (v≈0.9c)" begin
    T = Float64
    n = 64
    L = 2π
    cc = 5.0
    n0 = 1.0
    g = FourierGrid((n,), (L,))
    nppc = 100
    N = nppc * n
    # relativistic electron beam drifting at v≈0.9c; mobile ions form the neutral
    # return background. The numerical-Cherenkov instability would couple the beam
    # to grid EM modes and blow up the field energy — verify it stays bounded.
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, n0, g)
    vbeam = 0.9 * cc
    fill!(e.v[1], vbeam)
    ions = ParticleSet{1,T}(N; q = 1.0, m = 1836.0)
    load_lattice_1d!(ions, 0.0, L)
    set_density_weight!(ions, n0, g)
    es = EMPIC1D(g, N; n0 = n0, c = cc, shape = CIC(), relativistic = true, mobile = true)
    init_empic!(es, e, ions)
    dt = 0.005
    W0field = em_field_energy(es)
    # seed a tiny field floor so the ratio is meaningful (Poisson Ex ≈ 0 for a
    # neutral quiet start, so reference the run-max against the early-time field)
    Wmax = W0field
    Wearly = W0field
    finite = true
    for s = 1:3000
        step_empic!(es, e, ions, dt)
        W = em_field_energy(es)
        finite &= isfinite(W)
        Wmax = max(Wmax, W)
        s == 50 && (Wearly = max(W, eps()))
        # beam stays subluminal under the relativistic push
        @test maximum(abs, e.v[1]) < cc
    end
    @info "numerical-Cherenkov field energy" Wearly Wmax ratio = Wmax / Wearly
    @test finite
    # bounded: the field energy does not run away by orders of magnitude
    @test Wmax / Wearly < 100.0
end

@testset "Phase-12b EM PIC: relativistic push runs + bounded" begin
    T = Float64
    n = 64
    L = 2π
    cc = 5.0
    g = FourierGrid((n,), (L,))
    nppc = 100
    N = nppc * n
    rng = MersenneTwister(5)
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, 1.0, g)
    for p = 1:N
        e.v[1][p] = 0.3 * cc * randn(rng) * 0.1    # mildly relativistic
    end
    es = EMPIC1D(g, N; n0 = 1.0, c = cc, shape = CIC(), relativistic = true)
    for i = 1:n
        es.Ey[i] = 1e-3 * cos(2π * (i - 1) * g.dx[1] / L)
    end
    init_empic!(es, e)
    dt = 0.01
    W0 = em_field_energy(es) + kinetic_energy(e)
    Wmax = W0
    for _ = 1:1000
        step_empic!(es, e, dt)
        W = em_field_energy(es) + kinetic_energy(e)
        @test isfinite(W)
        Wmax = max(Wmax, W)
        # velocities stay subluminal
        @test maximum(abs, e.v[1]) < cc
    end
    @test Wmax / W0 < 2.0
end

@testset "EMPIC1D resizes electron gather buffers when the electron count grows" begin
    T = Float64
    g = FourierGrid((32,), (2π,))
    N = 64
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{1,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0,), g.L, (N,))
    load_lattice!(ions, (0.0,), g.L, (N,))
    es = EMPIC1D(g, N; mobile = true, mi = 25.0, c = 8.0)
    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.005)
    # grow the electron population (as ionize_mcc! secondaries would); the next step must not
    # DimensionMismatch on stale-sized electron gather/scratch buffers.
    xe = ParticleSet{1,T}(5; q = -1.0, m = 1.0)
    for k = 1:5
        xe.x[1][k] = 0.1k
    end
    HybridPlasmaPIC.append_particles!(e, xe)
    @test HybridPlasmaPIC.nparticles(e) == N + 5
    step_empic!(es, e, ions, 0.005)
    @test isfinite(em_field_energy(es) + kinetic_energy(e))
end

@testset "EMPIC1D-order: step_empic! is 2nd-order in dt (leapfrog v + Bz priming)" begin
    # The first step primes BOTH species' velocities (v^0→v^{-1/2}) and the seeded Bz
    # (Bz^0→Bz^{-1/2}); without it the integer-level fields are 1st-order. Cold-deterministic
    # longitudinal Langmuir self-convergence on Ex (integer level) — rate must sit near 2.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.02
    function ex(nsteps, Tf; lead_zero)
        n = 32
        N = 400 * n
        g = FourierGrid((n,), (T(L),))
        e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
        ions = ParticleSet{1,T}(N; q = 1.0, m = 100.0)
        load_lattice_1d!(e, 0.0, T(L))
        load_lattice_1d!(ions, 0.0, T(L))
        set_density_weight!(e, 1.0, g)
        set_density_weight!(ions, 1.0, g)
        for p = 1:N
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
        end
        es = EMPIC1D(g, N; mobile = true, mi = 100.0, c = 8.0)
        init_empic!(es, e, ions)
        lead_zero && step_empic!(es, e, ions, 0.0)         # dt=0 must not consume the priming
        for _ = 1:nsteps
            step_empic!(es, e, ions, Tf / nsteps)
        end
        abs(mode_amplitude(es.Ex, g, (1,)))
    end
    Tf = 1.0
    seq = (20, 40, 80, 160)
    for lz in (false, true)
        v = [ex(ns, Tf; lead_zero = lz) for ns in seq]
        r1 = log2(abs(v[1] - v[2]) / abs(v[2] - v[3]))
        r2 = log2(abs(v[2] - v[3]) / abs(v[3] - v[4]))
        @test r1 > 1.6        # ≈2.0 with priming; ≈1.0 (would fail) without it
        @test r2 > 1.6
    end
end

@testset "EMPIC1D relativistic priming is 2nd-order for a relativistic drifting beam" begin
    # The relativistic push works in momentum u=γv, so the prime must back up MOMENTUM
    # (u^{-1/2}=γ^0 v^0−h·a^0), not velocity. With a bulk drift the velocity-space prime leaves
    # an O(dt) momentum error → 1st-order; the momentum-space prime restores 2nd-order.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.02
    function ex(nsteps; vd)
        n = 32
        N = 400 * n
        g = FourierGrid((n,), (T(L),))
        e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
        ions = ParticleSet{1,T}(N; q = 1.0, m = 100.0)
        load_lattice_1d!(e, 0.0, T(L))
        load_lattice_1d!(ions, 0.0, T(L))
        set_density_weight!(e, 1.0, g)
        set_density_weight!(ions, 1.0, g)
        for p = 1:N
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
            e.v[1][p] = vd
        end
        es = EMPIC1D(g, N; mobile = true, mi = 100.0, c = 8.0, relativistic = true)
        init_empic!(es, e, ions)
        for _ = 1:nsteps
            step_empic!(es, e, ions, 1.0 / nsteps)
        end
        abs(mode_amplitude(es.Ex, g, (1,)))
    end
    v = [ex(ns; vd = 4.0) for ns in (20, 40, 80, 160)]   # drift 0.5c, γ≈1.15
    @test log2(abs(v[1] - v[2]) / abs(v[2] - v[3])) > 1.6
    @test log2(abs(v[2] - v[3]) / abs(v[3] - v[4])) > 1.6
end

@testset "EMPIC ion ParticleSet mass must match the constructor mi" begin
    T = Float64
    g = FourierGrid((16,), (2π,))
    es = EMPIC1D(g, 100; mobile = true, mi = 100.0, c = 5.0)
    e = ParticleSet{1,T}(100; q = -1.0, m = 1.0)
    bad = ParticleSet{1,T}(100; q = 1.0, m = 1.0)      # m ≠ mi=100 (was silently ignored)
    @test_throws ArgumentError init_empic!(es, e, bad)
end

@testset "EMPIC1D mobile subcycling (n_sub≥2) is 2nd-order and charge-conserving" begin
    # Splitting the ion drift across substeps (kick once on the first substep, drift dt_e each
    # substep) makes the mobile n_sub≥2 scheme 2nd-order (was 1st-order: the ion was drifted the
    # whole dt_ion at once, lumping its current into one substep) while keeping exact Esirkepov
    # charge conservation.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.05
    function ex(nsteps; n_sub)
        n = 32
        N = 400 * n
        g = FourierGrid((n,), (T(L),))
        e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
        ions = ParticleSet{1,T}(N; q = 1.0, m = 16.0)
        load_lattice_1d!(e, 0.0, T(L))
        load_lattice_1d!(ions, 0.0, T(L))
        set_density_weight!(e, 1.0, g)
        set_density_weight!(ions, 1.0, g)
        for p = 1:N
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
        end
        es = EMPIC1D(g, N; mobile = true, mi = 16.0, c = 8.0, n_sub = n_sub)
        init_empic!(es, e, ions)
        for _ = 1:nsteps
            step_empic!(es, e, ions, 1.0 / nsteps)
        end
        abs(mode_amplitude(es.Ex, g, (1,)))
    end
    for ns in (2, 3)      # n_sub=1 is covered elsewhere; ≥3 kicks first for all n_sub
        v = [ex(nsteps; n_sub = ns) for nsteps in (20, 40, 80, 160)]
        @test log2(abs(v[1] - v[2]) / abs(v[2] - v[3])) > 1.6
        @test log2(abs(v[2] - v[3]) / abs(v[3] - v[4])) > 1.6
    end
    # exact Esirkepov charge conservation per substep (dt_e) for a mobile n_sub=3 run
    n = 32
    N = 100 * n
    g = FourierGrid((n,), (T(L),))
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{1,T}(N; q = 1.0, m = 16.0)
    load_lattice_1d!(e, 0.0, T(L))
    load_lattice_1d!(ions, 0.0, T(L))
    set_density_weight!(e, 1.0, g)
    set_density_weight!(ions, 1.0, g)
    load_maxwellian!(e, MersenneTwister(2), (0.1, 0.0, 0.0), (0.3, 0.3, 0.3))
    es = EMPIC1D(g, N; mobile = true, mi = 16.0, c = 8.0, n_sub = 3)
    init_empic!(es, e, ions)
    dt = 0.02
    cmax = 0.0
    for _ = 1:20
        step_empic!(es, e, ions, dt)
        cmax = max(cmax, charge_conservation_residual(es, dt / 3))   # per-substep dt_e
    end
    @test cmax < 1e-10
end
