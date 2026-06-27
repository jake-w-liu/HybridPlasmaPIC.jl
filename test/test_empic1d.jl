# Phase-12b full electromagnetic 1D PIC. Oracles (all closed-form / exact):
#   (1) transverse EM-wave dispersion ω² = ω_pe² + c²k²  (≥2 k values, <3%);
#   (2) Esirkepov charge conservation  max|ρ^{n+1}−ρ^n + dt ∂xJx|/scale < 1e-10;
#   (3) total (field+kinetic) energy bounded over the run.

using HybridPlasmaPIC, FFTW, Test, Random

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
