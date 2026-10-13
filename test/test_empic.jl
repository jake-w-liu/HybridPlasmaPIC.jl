# Dimension-parametric electromagnetic PIC. These tests cover the 2D/3D paths
# that are intentionally separate from the legacy, more specialized EMPIC1D
# regression suite.

using HybridPlasmaPIC, FFTW, Test, Random

function _empic_peakfreq(series, dt)
    Nt = length(series)
    S = fft(series)
    mag = abs.(S)
    best = 2
    bm = -1.0
    for j = 2:(Nt÷2)
        if mag[j+1] > bm
            bm = mag[j+1]
            best = j
        end
    end
    a = mag[best]
    b = mag[best+1]
    c = mag[best+2]
    d = a - 2b + c
    δ = d != 0 ? 0.5 * (a - c) / d : 0.0
    return 2π * (best + δ) / (Nt * dt)
end

function _empic_gauss_residual(es, rho)
    g = es.g
    rhohat = complex.(rho)
    fft!(rhohat)
    for I in CartesianIndices(rhohat)
        k2 = zero(eltype(rho))
        for d = 1:length(g.n)
            k2 += g.kvec[d][I[d]]^2
        end
        k2 == 0 && (rhohat[I] = zero(eltype(rhohat)))
    end
    ifft!(rhohat)
    rhorepr = real.(rhohat)
    divE = similar(rho)
    divergence!(divE, es.E, g)
    scale = max(maximum(abs, rhorepr), 1.0)
    return maximum(abs, divE .- rhorepr) / scale
end

@testset "EMPIC constructor validation" begin
    g2 = FourierGrid((8, 6), (2π, 3π))
    es2 = EMPIC(g2, 5; n0 = 0.75, c = 4.0, mi = 100.0)
    @test es2.n0 == 0.75
    @test es2.c == 4.0
    @test es2.mi == 100.0
    @test size(es2.E[1]) == (8, 6)
    @test all(size(es2.B[c]) == (8, 6) for c = 1:3)
    @test all(length(es2.Ep[c]) == 5 for c = 1:3)
    @test all(length(es2.mide[d]) == 5 for d = 1:2)

    g3 = FourierGrid((6, 5, 4), (2π, 3π, 4π))
    es3 = EMPIC(g3, 7; n0 = 1.25)
    @test size(es3.rho_n) == (6, 5, 4)
    @test all(size(es3.J[c]) == (6, 5, 4) for c = 1:3)
    @test all(length(es3.mide[d]) == 7 for d = 1:3)

    @test_throws ArgumentError EMPIC(g2, 5; n0 = NaN)
    @test_throws ArgumentError EMPIC(g2, 5; n0 = -1.0)
    @test_throws ArgumentError EMPIC(g2, 5; c = NaN)
    @test_throws ArgumentError EMPIC(g2, 5; c = 0.0)
    @test_throws ArgumentError EMPIC(g2, 5; mi = 0.0)
    @test_throws ArgumentError EMPIC(g2, -1)

    g4 = FourierGrid((2, 2, 2, 2), (1.0, 1.0, 1.0, 1.0))
    @test_throws ArgumentError EMPIC(g4, 1)
end

@testset "EMPIC validates species and timestep before mutation" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 2π))
    e = ParticleSet{2,T}(48; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (8, 6))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, 48; n0 = 1.0, c = 5.0)
    init_empic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 2)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    E0 = ntuple(c -> copy(es.E[c]), 3)
    B0 = ntuple(c -> copy(es.B[c]), 3)
    rho_n0 = copy(es.rho_n)
    rho_np10 = copy(es.rho_np1)
    time0 = es.time[]
    step0 = es.step[]

    @test_throws ArgumentError step_empic!(es, e, NaN)
    @test_throws ArgumentError step_empic!(es, e, -0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.rho_n == rho_n0
    @test es.rho_np1 == rho_np10
    @test es.time[] == time0
    @test es.step[] == step0

    e.q = 0.0
    @test_throws ArgumentError step_empic!(es, e, 0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.time[] == time0
    @test es.step[] == step0
end

@testset "EMPIC rejects invalid relativistic velocities before mutation" begin
    g = FourierGrid((4, 4), (2π, 2π))
    e = ParticleSet{2,Float64}(16; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (4, 4))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, 16; c = 5.0, relativistic = true)
    init_empic!(es, e)
    e.v[1][end] = 5.5
    x0 = map(copy, e.x)
    v0 = map(copy, e.v)
    E0 = map(copy, es.E)
    B0 = map(copy, es.B)
    rho_n0, rho_np10 = copy(es.rho_n), copy(es.rho_np1)
    time0, step0 = es.time[], es.step[]

    @test_throws ArgumentError step_empic!(es, e, 0.01)

    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.rho_n == rho_n0 && es.rho_np1 == rho_np10
    @test es.time[] == time0 && es.step[] == step0
end

@testset "EMPIC init preserves transverse field and enforces Gauss law" begin
    T = Float64
    nx, ny = 16, 8
    L = 2π
    g = FourierGrid((nx, ny), (L, L))
    N = nx * ny
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (nx, ny))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, N; n0 = 1.0, c = 5.0)
    seed = zeros(T, nx, ny)
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        seed[i, j] = 1e-3 * cos(x)
        es.E[2][i, j] = seed[i, j]
    end

    init_empic!(es, e)
    @test maximum(abs, es.E[2] .- seed) < 1e-12
    @test _empic_gauss_residual(es, es.rho_n) < 1e-12
end

@testset "EMPIC 2D transverse dispersion ω²=ωpe²+c²k²" begin
    T = Float64
    nx, ny = 32, 4
    L = 2π
    cc = 5.0
    n0 = 1.0
    g = FourierGrid((nx, ny), (L, L))
    nppc = 80
    N = nppc * nx * ny
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (nx * nppc, ny))
    set_density_weight!(e, n0, g)
    es = EMPIC(g, N; n0 = n0, c = cc, shape = CIC())
    m = 1
    k = 2π * m / L
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        es.E[2][i, j] = 1e-3 * cos(k * x)
    end
    init_empic!(es, e)
    dt = 0.01
    series = ComplexF64[]
    for _ = 1:1000
        step_empic!(es, e, dt)
        push!(series, mode_amplitude(es.E[2], g, (m, 0)))
    end
    ω = _empic_peakfreq(series, dt)
    ωth = sqrt(n0 + cc^2 * k^2)
    @test abs(ω - ωth) / ωth < 0.03
    @test charge_conservation_residual(es, dt) < 1e-8
    @test _empic_gauss_residual(es, es.rho_np1) < 1e-8
end

@testset "EMPIC 2D/3D spectral current correction conserves charge" begin
    T = Float64

    rng2 = MersenneTwister(11)
    g2 = FourierGrid((16, 12), (2π, 2π))
    e2 = ParticleSet{2,T}(1200; q = -1.0, m = 1.0)
    load_uniform!(e2, rng2, (0.0, 0.0), g2.L)
    set_density_weight!(e2, 1.0, g2)
    for p in eachindex(e2.weight)
        e2.v[1][p] = 0.2 * randn(rng2)
        e2.v[2][p] = 0.2 * randn(rng2)
        e2.v[3][p] = 0.1 * randn(rng2)
    end
    es2 = EMPIC(g2, nparticles(e2); n0 = 1.0, c = 5.0)
    init_empic!(es2, e2)
    maxres2 = 0.0
    for _ = 1:20
        step_empic!(es2, e2, 0.005)
        maxres2 = max(maxres2, charge_conservation_residual(es2, 0.005))
    end
    @test maxres2 < 1e-10
    @test _empic_gauss_residual(es2, es2.rho_np1) < 1e-10

    rng3 = MersenneTwister(13)
    g3 = FourierGrid((8, 6, 4), (2π, 2π, 2π))
    e3 = ParticleSet{3,T}(384; q = -1.0, m = 1.0)
    load_uniform!(e3, rng3, (0.0, 0.0, 0.0), g3.L)
    set_density_weight!(e3, 1.0, g3)
    for p in eachindex(e3.weight)
        e3.v[1][p] = 0.15 * randn(rng3)
        e3.v[2][p] = 0.15 * randn(rng3)
        e3.v[3][p] = 0.15 * randn(rng3)
    end
    es3 = EMPIC(g3, nparticles(e3); n0 = 1.0, c = 4.0)
    init_empic!(es3, e3)
    maxres3 = 0.0
    for _ = 1:10
        step_empic!(es3, e3, 0.003)
        maxres3 = max(maxres3, charge_conservation_residual(es3, 0.003))
    end
    @test maxres3 < 1e-10
    @test _empic_gauss_residual(es3, es3.rho_np1) < 1e-10
end

@testset "EMPIC 2D mobile-ion path" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 2π))
    N = prod(g.n)
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0, 0.0), g.L, g.n)
    load_lattice!(ions, (0.0, 0.0), g.L, g.n)
    set_density_weight!(e, 1.0, g)
    set_density_weight!(ions, 1.0, g)
    es = EMPIC(g, N; mobile = true, mi = 25.0, c = 4.0)
    @test_throws ArgumentError init_empic!(es, e)

    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.01)
    @test charge_conservation_residual(es, 0.01) < 1e-12
    @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
    @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))

    badions = ParticleSet{2,T}(N; q = 0.0, m = 25.0)
    load_lattice!(badions, (0.0, 0.0), g.L, g.n)
    set_density_weight!(badions, 1.0, g)
    @test_throws ArgumentError step_empic!(es, e, badions, 0.01)
end

@testset "EMPIC resizes electron gather buffers when the electron count grows" begin
    T = Float64
    g = FourierGrid((16, 16), (2π, 2π))
    N = 64
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0, 0.0), g.L, (8, 8))
    load_lattice!(ions, (0.0, 0.0), g.L, (8, 8))
    es = EMPIC(g, N; mobile = true, mi = 25.0, c = 8.0)
    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.005)
    # grow the electron population (as ionize_mcc! secondaries would); the next step must not
    # DimensionMismatch on stale-sized electron gather buffers (they are ion-symmetric now).
    xe = ParticleSet{2,T}(5; q = -1.0, m = 1.0)
    for d = 1:2, k = 1:5
        xe.x[d][k] = 0.1k
    end
    HybridPlasmaPIC.append_particles!(e, xe)
    @test HybridPlasmaPIC.nparticles(e) == N + 5
    step_empic!(es, e, ions, 0.005)
    @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))
end

@testset "EMPIC-order: step_empic! is 2nd-order in dt (leapfrog v + B priming)" begin
    # The multi-dim EM PIC primes both species' velocities and the seeded B on the first step
    # (same fix as EMPIC1D). Cold 2D longitudinal Langmuir self-convergence on Ex — rate ≈ 2.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.02
    function ex(nsteps, Tf)
        n = 16
        pd = 80
        N = pd * pd
        g = FourierGrid((n, n), (T(L), T(L)))
        e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
        ions = ParticleSet{2,T}(N; q = 1.0, m = 100.0)
        load_lattice!(e, (0.0, 0.0), g.L, (pd, pd))
        load_lattice!(ions, (0.0, 0.0), g.L, (pd, pd))
        set_density_weight!(e, 1.0, g)
        set_density_weight!(ions, 1.0, g)
        for p = 1:N
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
        end
        es = EMPIC(g, N; mobile = true, mi = 100.0, c = 8.0)
        init_empic!(es, e, ions)
        for _ = 1:nsteps
            step_empic!(es, e, ions, Tf / nsteps)
        end
        abs(mode_amplitude(es.E[1], g, (1, 0)))
    end
    Tf = 1.0
    seq = (40, 80, 160, 320)
    v = [ex(ns, Tf) for ns in seq]
    r1 = log2(abs(v[1] - v[2]) / abs(v[2] - v[3]))
    r2 = log2(abs(v[2] - v[3]) / abs(v[3] - v[4]))
    @test r1 > 1.6
    @test r2 > 1.6
end

@testset "EMPIC ↔ EMPIC1D spectral equivalence incl. Nyquist current exclusion" begin
    # The continuity correction zeroes the pure-Nyquist current modes (as
    # _esirkepov_Jx! / _deposit_Jy_at_midpoints! do in EMPIC1D), so the two 1D
    # solvers must agree to roundoff on EVERY Fourier mode and Ex must never
    # accumulate a grid-Nyquist sawtooth. Pre-fix: the raw |Ĵx_Nyq| survived
    # only in EMPIC (1.5e-4 after one step here), Ex's Nyquist mode
    # random-walked to a few % of max|Ex|, and the solvers drifted apart by
    # ~2.5% within 200 steps.
    T = Float64
    n = 32
    L = 2π
    N = 1280
    dt = 0.01
    g1 = FourierGrid((n,), (L,))
    g2 = FourierGrid((n,), (L,))
    rng = MersenneTwister(1234)
    e1 = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_uniform!(e1, rng, (0.0,), (L,))
    set_density_weight!(e1, 1.0, g1)
    for p = 1:N
        e1.v[1][p] = 0.1 * randn(rng)
        e1.v[2][p] = 0.1 * randn(rng)     # exercise the transverse Jy channel
    end
    e2 = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    e2.x[1] .= e1.x[1]
    for c = 1:3
        e2.v[c] .= e1.v[c]
    end
    e2.weight .= e1.weight
    es1 = EMPIC1D(g1, N; n0 = 1.0, c = 5.0, shape = CIC())
    es2 = EMPIC(g2, N; n0 = 1.0, c = 5.0, shape = CIC())
    init_empic!(es1, e1)
    init_empic!(es2, e2)

    nyq(f) = abs(fft(complex.(f))[length(f)÷2+1]) / length(f)

    # one step: current spectra agree on ALL modes; Nyquist current is zero
    step_empic!(es1, e1, dt)
    step_empic!(es2, e2, dt)
    @test maximum(abs.(fft(complex.(es1.Jx)) .- fft(complex.(es2.J[1])))) / n < 1e-14
    @test maximum(abs.(fft(complex.(es1.Jy)) .- fft(complex.(es2.J[2])))) / n < 1e-14
    @test nyq(es2.J[1]) < 1e-15
    @test nyq(es2.J[2]) < 1e-15

    # 300 steps: solvers stay equivalent to roundoff and the Nyquist mode of
    # Ex stays at machine zero (measured 9e-15 / 5.5e-16 with the fix; the
    # pre-fix values 2.5e-2 / ≳1e-3 fail these bounds by many orders)
    nyqmax = 0.0
    for _ = 2:300
        step_empic!(es1, e1, dt)
        step_empic!(es2, e2, dt)
        nyqmax = max(nyqmax, nyq(es2.E[1]) / maximum(abs.(es2.E[1])))
    end
    @test maximum(abs.(es1.Ex .- es2.E[1])) / maximum(abs.(es1.Ex)) < 1e-12
    @test maximum(abs.(es1.Ey .- es2.E[2])) / maximum(abs.(es1.Ey)) < 1e-12
    @test nyqmax < 1e-13
    @test charge_conservation_residual(es2, dt) < 1e-10
end

@testset "EM PIC changing timestep preserves electron and ion half levels" begin
    # Constant E gives an exact momentum oracle, including relativistic motion.
    # Zero-weight tracers leave the imposed field unchanged through Maxwell.
    for (constructor, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
        relativistic in (false, true),
        mobile in (false, true)

        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0, D))
        e = ParticleSet{D,Float64}(1; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(1; q = 1.0, m = 16.0)
        e.weight .= 0.0
        ions.weight .= 0.0
        for d = 1:D
            e.x[d] .= 1.0
            ions.x[d] .= 2.0
        end
        es = constructor(g, 1; n0 = 0.0, c = 2.0, mobile, mi = 16.0, relativistic)
        species = mobile ? ions : nothing
        init_empic!(es, e, species)
        function set_field!()
            if es isa EMPIC1D
                es.Ex .= -0.7
                es.Ey .= 0.2
            else
                es.E[1] .= -0.7
                es.E[2] .= 0.2
            end
        end
        set_field!()
        elapsed = 0.0
        for (dt, n_sub) in ((0.1, 2), (0.2, 3), (0.05, 1))
            es isa EMPIC1D && (es.n_sub = n_sub)
            previous_dt = es.last_dt[]
            step_empic!(es, e, species, 0.0)
            @test es.last_dt[] == previous_dt
            step_empic!(es, e, species, dt)
            elapsed += dt
            electron_dt = es isa EMPIC1D ? dt / n_sub : dt
            for (s, acceleration, h) in
                ((e, (0.7, -0.2, 0.0), electron_dt), (ions, (-0.7 / 16, 0.2 / 16, 0.0), dt))
                s === ions && !mobile && continue
                u = ntuple(c -> acceleration[c] * (elapsed - h / 2), 3)
                gamma = relativistic ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0
                @test [s.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 2e-15
                if !relativistic
                    for d = 1:D
                        x0 = s === e ? 1.0 : 2.0
                        @test s.x[d][1] ≈ x0 + acceleration[d] * elapsed^2 / 2 atol = 2e-15
                    end
                end
            end
        end
        before = (deepcopy(e.x), deepcopy(e.v), es.time[], es.step[], es.last_dt[])
        for dt in (NaN, Inf, -0.1)
            @test_throws ArgumentError step_empic!(es, e, species, dt)
            @test (e.x, e.v, es.time[], es.step[], es.last_dt[]) == before
        end
        # Reinitialization takes freshly loaded physical velocities; a previous
        # run's unequal timesteps must not leak into the next run's prime.
        for c = 1:3
            e.v[c] .= 0.0
            ions.v[c] .= 0.0
        end
        init_empic!(es, e, species)
        @test es.last_dt[] == 0.0
        es isa EMPIC1D && @test es.last_sub_dt[] == 0.0
        set_field!()
        step_empic!(es, e, species, 0.08)
        u = (0.7 * 0.04, -0.2 * 0.04, 0.0)
        gamma = relativistic ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0
        @test [e.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 2e-15
    end
end

@testset "EM PIC changing timestep preserves Maxwell time staggering" begin
    function variable_dt_vacuum_error(constructor, refinement)
        g = FourierGrid((16,), (2π,))
        e = ParticleSet{1,Float64}(0; q = -1.0, m = 1.0)
        es = constructor(g, 0; n0 = 0.0, c = 2.0)
        init_empic!(es, e)
        E = es isa EMPIC1D ? es.Ey : es.E[2]
        B = es isa EMPIC1D ? es.Bz : es.B[3]
        x = range(0.0; step = g.dx[1], length = g.n[1])
        E .= sin.(x)
        for (count, dt) in ((10refinement, 0.04 / refinement), (20refinement, 0.02 / refinement))
            for _ = 1:count
                step_empic!(es, e, dt)
            end
        end
        # Analytic vacuum mode: Ey=sin(x)cos(ct), Bz=-cos(x)sin(ct)/c.
        # E is at t=0.8, while the stored B trails it by the last dt/2.
        Eerror = maximum(abs, E .- cos(1.6) .* sin.(x))
        Berror = maximum(abs, B .+ sin(2 * (0.8 - 0.01 / refinement)) / 2 .* cos.(x))
        return (Eerror, Berror)
    end
    for constructor in (EMPIC1D, EMPIC)
        errors = [variable_dt_vacuum_error(constructor, r) for r in (1, 2, 4)]
        for component = 1:2, level = 1:2
            @test log2(errors[level][component] / errors[level+1][component]) > 1.8
        end
    end
end

@testset "EM PIC changing timestep gyro-orbit convergence" begin
    function variable_dt_gyro_error(constructor, relativistic, refinement)
        g = FourierGrid((8,), (10.0,))
        e = ParticleSet{1,Float64}(1; q = -1.0, m = 1.0)
        e.weight .= 0.0
        e.x[1] .= 1.0
        e.v[1] .= 1.0
        es = constructor(g, 1; n0 = 0.0, c = 2.0, relativistic)
        init_empic!(es, e)
        if es isa EMPIC1D
            es.Bz .= 1.0
        else
            es.B[3] .= 1.0
        end
        for (count, dt) in ((10refinement, 0.04 / refinement), (20refinement, 0.02 / refinement))
            for _ = 1:count
                step_empic!(es, e, dt)
            end
        end
        omega = relativistic ? sqrt(1 - 1 / es.c^2) : 1.0
        return abs(e.x[1][1] - (1 + sin(omega * 0.8) / omega))
    end
    for constructor in (EMPIC1D, EMPIC), relativistic in (false, true)
        errors = [variable_dt_gyro_error(constructor, relativistic, r) for r in (1, 2, 4)]
        @test all(log2(errors[j] / errors[j+1]) > 1.8 for j = 1:2)
    end
end

@testset "EM PIC changing timestep retains charge and Gauss constraints" begin
    for (constructor, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3))
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> 2π, D))
        e = ParticleSet{D,Float64}(8^D; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(8^D; q = 1.0, m = 16.0)
        for s in (e, ions)
            load_lattice!(s, ntuple(_ -> 0.0, D), g.L, ntuple(_ -> 8, D))
            set_density_weight!(s, 1.0, g)
        end
        @. e.x[1] -= 0.02 * sin(e.x[1])
        e.v[1] .= 0.05
        ions.v[1] .= 0.02
        es = constructor(g, nparticles(e); mobile = true, mi = 16.0, c = 2.0)
        init_empic!(es, e, ions)
        for (dt, n_sub) in ((0.01, 1), (0.005, 3), (0.02, 2), (0.01, 1))
            es isa EMPIC1D && (es.n_sub = n_sub)
            step_empic!(es, e, ions, dt)
            sub_dt = es isa EMPIC1D ? dt / n_sub : dt
            delta_hat = fft(es.rho_np1 .- es.rho_n)
            dc = first(CartesianIndices(delta_hat))
            for I in CartesianIndices(delta_hat)
                if I != dc && all(iszero(g.kvec[d][I[d]]) for d = 1:D)
                    delta_hat[I] = 0
                end
            end
            delta = real.(ifft(delta_hat))
            divJ = similar(delta)
            if es isa EMPIC1D
                deriv!(divJ, es.Jx, g, 1)
            else
                divergence!(divJ, es.J, g)
            end
            # Density is normalized to one. Keep the DC conservation check,
            # whose deposition roundoff can dominate a tiny density increment;
            # bound absolute charge error rather than divide it by that increment.
            @test maximum(abs, delta .+ sub_dt .* divJ) < 32eps(Float64)
            @test charge_conservation_residual(es, sub_dt) * maximum(abs, delta) < 32eps(Float64)
            if es isa EMPIC1D
                rhohat = fft(es.rho_np1)
                rhohat[1] = rhohat[g.n[1]÷2+1] = 0
                rho = real.(ifft(rhohat))
                divE = similar(rho)
                deriv!(divE, es.Ex, g, 1)
                @test maximum(abs, divE .- rho) < 1e-12
            else
                @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
            end
        end
    end
end

function _empic_subcycling_fixture(
    D;
    n_sub = 1,
    shape = CIC(),
    mobile = false,
    relativistic = false,
    mi = 16.0,
)
    g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> 2π, D))
    e = ParticleSet{D,Float64}(8^D; q = -1.0, m = 1.0)
    ions = ParticleSet{D,Float64}(8^D; q = 1.0, m = mi)
    for s in (e, ions)
        load_lattice!(s, ntuple(_ -> 0.0, D), g.L, g.n)
        set_density_weight!(s, 1.0, g)
    end
    @. e.x[1] -= 0.02 * sin(e.x[1])
    @. e.v[1] = 0.1 + 0.02 * cos(e.x[1])
    @. e.v[2] = 0.01 * sin(e.x[1])
    ions.v[1] .= 0.03
    es = EMPIC(g, nparticles(e); n_sub, shape, mobile, relativistic, mi, c = 3.0)
    for I in CartesianIndices(es.E[1])
        x = (I[1] - 1) * g.dx[1]
        es.E[2][I] = 0.01 * sin(x)
        es.E[3][I] = 0.02 * cos(x)
        es.B[2][I] = 0.01 * cos(x)
        es.B[3][I] = 0.02 * sin(x)
    end
    init_empic!(es, e, mobile ? ions : nothing)
    return es, e, ions
end

@testset "EMPIC subcycling validates and resets timing state" begin
    for D = 1:3
        es, e, ions = _empic_subcycling_fixture(D; n_sub = 3, mobile = true)
        @test es.n_sub == 3
        @test_throws ArgumentError EMPIC(es.g, 1; n_sub = 0)
        @test_throws ArgumentError EMPIC(es.g, 1; n_sub = -1)
        step_empic!(es, e, ions, 0.03)
        before = (
            deepcopy(e.x),
            deepcopy(e.v),
            deepcopy(ions.x),
            deepcopy(ions.v),
            deepcopy(es.E),
            deepcopy(es.B),
            es.time[],
            es.step[],
            es.last_dt[],
            es.last_sub_dt[],
        )
        step_empic!(es, e, ions, 0.0)
        @test (
            e.x,
            e.v,
            ions.x,
            ions.v,
            es.E,
            es.B,
            es.time[],
            es.step[],
            es.last_dt[],
            es.last_sub_dt[],
        ) == before
        for n_sub in (0, -1)
            es.n_sub = n_sub
            @test_throws ArgumentError step_empic!(es, e, ions, 0.03)
            @test (
                e.x,
                e.v,
                ions.x,
                ions.v,
                es.E,
                es.B,
                es.time[],
                es.step[],
                es.last_dt[],
                es.last_sub_dt[],
            ) == before
        end
        es.n_sub = 3
        @test_throws ArgumentError step_empic!(es, e, ions, nextfloat(0.0))
        @test (
            e.x,
            e.v,
            ions.x,
            ions.v,
            es.E,
            es.B,
            es.time[],
            es.step[],
            es.last_dt[],
            es.last_sub_dt[],
        ) == before
        init_empic!(es, e, ions)
        @test (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]) == (0.0, 0, 0.0, 0.0)
    end
end

@testset "EMPIC immobile subcycling equals explicit electron steps" begin
    for D = 1:3, shape in (CIC(), TSC()), relativistic in (false, true), ns in (1, 2, 4)
        es, e, _ = _empic_subcycling_fixture(D; n_sub = ns, shape, relativistic)
        fine, efine, _ = _empic_subcycling_fixture(D; shape, relativistic)
        for (dt, n_sub) in ((0.01, ns), (0.02, 3), (0.015, ns))
            es.n_sub = n_sub
            step_empic!(es, e, dt)
            for _ = 1:n_sub
                step_empic!(fine, efine, dt / n_sub)
            end
            @test e.x == efine.x
            @test e.v == efine.v
            @test es.E == fine.E
            @test es.B == fine.B
            @test es.J == fine.J
            @test es.time[] == fine.time[]
            @test es.last_sub_dt[] == fine.last_sub_dt[]
        end
        @test es.step[] == 3
    end
end

@testset "EMPIC subcycling agrees with polarized EMPIC1D" begin
    for shape in (CIC(), TSC()), relativistic in (false, true), mobile in (false, true), ns = 1:4
        g = FourierGrid((16,), (2π,))
        e = ParticleSet{1,Float64}(256; q = -1.0, m = 1.0)
        ions = ParticleSet{1,Float64}(256; q = 1.0, m = 16.0)
        for s in (e, ions)
            load_lattice_1d!(s, 0.0, g.L[1])
            set_density_weight!(s, 1.0, g)
        end
        @. e.x[1] -= 0.02 * sin(e.x[1])
        @. e.v[2] = 0.02 * cos(e.x[1])
        e.v[1] .= relativistic ? 1.0 : 0.1
        es = EMPIC(g, nparticles(e); n_sub = ns, shape, relativistic, mobile, mi = 16.0, c = 3.0)
        ref = EMPIC1D(
            FourierGrid(g.n, g.L),
            nparticles(e);
            n_sub = ns,
            shape,
            relativistic,
            mobile,
            mi = 16.0,
            c = 3.0,
        )
        eref, iref = deepcopy(e), deepcopy(ions)
        for j = 1:g.n[1]
            x = (j - 1) * g.dx[1]
            es.E[2][j] = ref.Ey[j] = 0.01 * sin(x)
            es.B[3][j] = ref.Bz[j] = 0.02 * cos(x)
        end
        init_empic!(es, e, mobile ? ions : nothing)
        init_empic!(ref, eref, mobile ? iref : nothing)
        for (dt, n_sub) in ((0.01, ns), (0.02, 3), (0.015, 2))
            es.n_sub = ref.n_sub = n_sub
            step_empic!(es, e, mobile ? ions : nothing, dt)
            step_empic!(ref, eref, mobile ? iref : nothing, dt)
            @test all(isapprox.(e.x, eref.x; atol = 2e-13))
            @test all(isapprox.(e.v, eref.v; atol = 2e-13))
            @test all(isapprox.(ions.x, iref.x; atol = 2e-13))
            @test all(isapprox.(ions.v, iref.v; atol = 2e-13))
            @test es.E[1] ≈ ref.Ex atol = 2e-13
            @test es.E[2] ≈ ref.Ey atol = 2e-13
            @test es.B[3] ≈ ref.Bz atol = 2e-13
            @test es.time[] == ref.time[]
        end
    end
end

@testset "EMPIC multidimensional mobile subcycling conserves charge after population growth" begin
    for D in (2, 3), shape in (CIC(), TSC())
        es, e, ions = _empic_subcycling_fixture(D; n_sub = 3, shape, mobile = true)
        for (dt, ns) in ((0.01, 3), (0.02, 2), (0.005, 4))
            es.n_sub = ns
            step_empic!(es, e, ions, dt)
            @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
            delta_scale = maximum(abs, es.rho_np1 .- es.rho_n)
            @test charge_conservation_residual(es, dt / ns) * delta_scale < 128eps(Float64)
            @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))
            # An equal-position, equal-weight pair leaves the total charge
            # unchanged at insertion and exercises both lazy workspace resizes.
            for s in (e, ions)
                extra = ParticleSet{D,Float64}(2; q = s.q, m = s.m)
                for d = 1:D
                    extra.x[d] .= (0.1d, 0.2d)
                end
                extra.weight .= 0.001
                HybridPlasmaPIC.append_particles!(s, extra)
            end
        end
        step_empic!(es, e, ions, 0.005)
        @test all(length(v) == nparticles(e) for v in (es.Ep..., es.Bp..., es.mide..., es.worke))
        @test all(
            length(v) == nparticles(ions) for v in (es.Epi..., es.Bpi..., es.midi..., es.worki)
        )
        @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
    end
end

@testset "EMPIC multidimensional mobile subcycling has second-order temporal convergence" begin
    function subcycling_amplitude(D, mi, transverse, refinement)
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> 2π, D))
        e = ParticleSet{D,Float64}(16^D; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(16^D; q = 1.0, m = mi)
        for s in (e, ions)
            load_lattice!(s, ntuple(_ -> 0.0, D), g.L, ntuple(_ -> 16, D))
            set_density_weight!(s, 1.0, g)
        end
        es = EMPIC(g, nparticles(e); mobile = true, mi, c = 2.0, n_sub = 3)
        if transverse
            # The polarization (1,-1,0) is orthogonal to k=(1,...,1).
            for I in CartesianIndices(es.E[1])
                phase = sum((I[d] - 1) * g.dx[d] for d = 1:D)
                es.E[1][I] = 0.001 * sin(phase)
                es.E[2][I] = -es.E[1][I]
            end
        else
            for p = 1:nparticles(e)
                phase = sum(e.x[d][p] for d = 1:D)
                e.x[1][p] -= 0.01 * sin(phase)
            end
        end
        init_empic!(es, e, ions)
        for (count, dt, ns) in
            ((10refinement, 0.04 / refinement, 2), (20refinement, 0.02 / refinement, 3))
            es.n_sub = ns
            for _ = 1:count
                step_empic!(es, e, ions, dt)
            end
        end
        return abs(mode_amplitude(es.E[1], g, ntuple(_ -> 1, D)))
    end
    for D in (2, 3), mi in (16.0, 100.0), transverse in (false, true)
        values = [subcycling_amplitude(D, mi, transverse, r) for r in (1, 2, 4)]
        @test log2(abs(values[1] - values[2]) / abs(values[2] - values[3])) > 1.7
    end
end

@testset "EMPIC multidimensional subcycling matches constant-field particle oracles" begin
    for D in (2, 3), relativistic in (false, true)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0, D))
        e = ParticleSet{D,Float64}(1; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(1; q = 1.0, m = 16.0)
        for (s, x0) in ((e, 1.0), (ions, 2.0))
            s.weight .= 0.0
            for d = 1:D
                s.x[d] .= x0
            end
        end
        es = EMPIC(g, 1; n0 = 0.0, mobile = true, mi = 16.0, c = 2.0, n_sub = 2, relativistic)
        init_empic!(es, e, ions)
        force = (0.7, -0.2, 0.1)
        for c = 1:3
            es.E[c] .= -force[c]
        end
        elapsed = 0.0
        for (dt, ns) in ((0.1, 2), (0.2, 3), (0.05, 4))
            es.n_sub = ns
            step_empic!(es, e, ions, dt)
            elapsed += dt
            for (s, x0, half_dt) in ((e, 1.0, dt / (2ns)), (ions, 2.0, dt / 2))
                acceleration = ntuple(c -> -s.q / s.m * force[c], 3)
                u = ntuple(c -> acceleration[c] * (elapsed - half_dt), 3)
                gamma = relativistic ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0
                @test [s.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 5e-14
                if !relativistic
                    for d = 1:D
                        @test s.x[d][1] ≈ x0 + acceleration[d] * elapsed^2 / 2 atol = 5e-14
                    end
                end
            end
            @test es.last_sub_dt[] == dt / ns
            @test es.last_dt[] == dt
        end
    end
end

@testset "Float32 EM subcycling charge and phase constraints" begin
    for D = 1:3, mobile in (false, true), rel in (false, true), shape in (CIC(), TSC())
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> Float32(2π), D))
        e = ParticleSet{D,Float32}(8^D; q = -1.0f0, m = 1.0f0)
        ions = ParticleSet{D,Float32}(8^D; q = 1.0f0, m = 16.0f0)
        for s in (e, ions)
            load_lattice!(s, ntuple(_ -> 0.0f0, D), g.L, g.n)
            set_density_weight!(s, 1.0f0, g)
        end
        @. e.x[1] -= 0.02f0 * sin(e.x[1])
        e.v[1] .= 0.1f0
        @. e.v[2] = 0.02f0 * cos(e.x[1])
        es = EMPIC(
            g,
            nparticles(e);
            mobile,
            mi = 16.0f0,
            c = 2.0f0,
            relativistic = rel,
            shape,
            n_sub = 2,
        )
        init_empic!(es, e, mobile ? ions : nothing)
        for (dt, ns) in ((0.01, 2), (0.02, 3), (0.015, 1))
            es.n_sub = ns
            step_empic!(es, e, mobile ? ions : nothing, dt)
            @test es.last_dt[] === Float32(dt)
            @test es.last_sub_dt[] === Float32(dt) / ns
            @test _empic_gauss_residual(es, es.rho_np1) < 128eps(Float32)
            ds = maximum(abs, es.rho_np1 .- es.rho_n)
            @test charge_conservation_residual(es, Float32(dt) / ns) * ds < 128eps(Float32)
            @test all(all(isfinite, v) for v in (e.x..., e.v..., es.E..., es.B...))
            if rel
                @test maximum(sum(e.v[c] .^ 2 for c = 1:3)) < es.c^2
            end
        end
    end
end
@testset "Float32 EM subcycling constant-field momentum" begin
    for D = 1:3, mobile in (false, true), rel in (false, true)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0f0, D))
        e = ParticleSet{D,Float32}(1; q = -1.0f0, m = 1.0f0)
        ions = ParticleSet{D,Float32}(1; q = 1.0f0, m = 16.0f0)
        for s in (e, ions)
            s.weight .= 0.0f0
            for d = 1:D
                s.x[d] .= 1.0f0
            end
        end
        es = EMPIC(g, 1; n0 = 0.0f0, mobile, mi = 16.0f0, c = 2.0f0, relativistic = rel, n_sub = 2)
        init_empic!(es, e, mobile ? ions : nothing)
        force = (0.7f0, -0.2f0, 0.1f0)
        for c = 1:3
            es.E[c] .= -force[c]
        end
        elapsed = 0.0f0
        for (dt, ns) in ((0.1f0, 2), (0.2f0, 3), (0.05f0, 4))
            es.n_sub = ns
            step_empic!(es, e, mobile ? ions : nothing, dt)
            elapsed += dt
            for (s, h) in ((e, dt / (2ns)), (ions, dt / 2))
                s === ions && !mobile && continue
                u = ntuple(c -> -s.q / s.m * force[c] * (elapsed - h), 3)
                gamma = rel ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0f0
                @test [s.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 32eps(Float32)
            end
        end
    end
end
