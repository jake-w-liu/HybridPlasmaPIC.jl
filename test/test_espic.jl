# Phase-12 full-PIC (electrostatic): cold Langmuir oscillation ω = ω_pe, the
# two-stream instability growth rate γ = ω_pe/(2√2), and a falsifiable null
# (kv0 > ω_pe ⇒ no growth). Oracles are closed-form.

using HybridPlasmaPIC, FFTW, Test, Random

@testset "Electrostatic1D parameter validation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    es = Electrostatic1D(g, 4; n0 = 1.0)
    @test es.n0 == 1.0
    @test_throws ArgumentError Electrostatic1D(g, 4; n0 = NaN)
    @test_throws ArgumentError Electrostatic1D(g, 4; n0 = -1.0)
end

# dominant positive-frequency peak (parabolic interpolation)
function _peakfreq(series, dt)
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

# exponential growth rate of |E1| over its clean log-linear window
function _growth(amp, tt)
    amax = maximum(amp)
    lo = findfirst(a -> a > 20 * amp[1], amp)
    hi = findfirst(a -> a > 0.3 * amax, amp)
    lo = lo === nothing ? 1 : lo
    hi = hi === nothing ? length(amp) : hi
    hi <= lo && return 0.0
    return (log(amp[hi]) - log(amp[lo])) / (tt[hi] - tt[lo])
end

@testset "Phase-12 electrostatic PIC: Langmuir ω = ω_pe" begin
    T = Float64
    n = 64
    L = 2π
    k = 2π / L
    g = FourierGrid((n,), (L,))
    nppc = 400
    N = nppc * n
    e = ParticleSet{1,T}(N)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, 1.0, g)
    rng = MersenneTwister(1)
    A = 0.05
    for p = 1:N
        e.v[1][p] = 0.03 * randn(rng)
        e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)   # δn = A cos(kx)
    end
    es = Electrostatic1D(g, N; n0 = 1.0)
    init_espic!(es, e)
    dt = 0.05
    series = ComplexF64[]
    W = Float64[]
    for _ = 1:1200
        step_espic!(es, e, dt)
        push!(series, mode_amplitude(es.E, g, (1,)))
        push!(W, field_energy(es) + kinetic_energy(e))
    end
    ω = _peakfreq(series, dt)
    @test abs(ω - 1.0) < 0.04                       # ω_pe = 1
    @test (maximum(W) - minimum(W)) / W[1] < 0.05   # total energy bounded
end

@testset "Phase-12 electrostatic PIC: two-stream growth rate" begin
    T = Float64
    v0 = 1.0
    k = 0.6124
    L = 2π / k    # fastest-growing mode
    n = 64
    g = FourierGrid((n,), (L,))
    nppc = 800
    N = nppc * n
    e = ParticleSet{1,T}(N)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, 1.0, g)
    for p = 1:N
        e.v[1][p] = isodd(p) ? v0 : -v0
        e.x[1][p] = mod(e.x[1][p] - (0.001 / k) * sin(k * e.x[1][p]), L)
    end
    es = Electrostatic1D(g, N; n0 = 1.0)
    init_espic!(es, e)
    dt = 0.05
    amp = Float64[]
    tt = Float64[]
    for st = 1:800
        step_espic!(es, e, dt)
        push!(amp, abs(mode_amplitude(es.E, g, (1,))))
        push!(tt, st * dt)
    end
    γ = _growth(amp, tt)
    @test abs(γ - 0.35355) / 0.35355 < 0.08          # γ = ω_pe/(2√2) ≈ 0.3536
end

@testset "Phase-12 electrostatic PIC: stable null (kv0 > ω_pe)" begin
    T = Float64
    v0 = 1.0
    k = 1.2
    L = 2π / k        # kv0 = 1.2 > ω_pe ⇒ no instability
    n = 64
    g = FourierGrid((n,), (L,))
    nppc = 800
    N = nppc * n
    e = ParticleSet{1,T}(N)
    load_lattice_1d!(e, 0.0, L)
    set_density_weight!(e, 1.0, g)
    for p = 1:N
        e.v[1][p] = isodd(p) ? v0 : -v0
        e.x[1][p] = mod(e.x[1][p] - (0.001 / k) * sin(k * e.x[1][p]), L)
    end
    es = Electrostatic1D(g, N; n0 = 1.0)
    init_espic!(es, e)
    dt = 0.05
    a0 = abs(mode_amplitude(es.E, g, (1,)))
    amax = a0
    for _ = 1:800
        step_espic!(es, e, dt)
        amax = max(amax, abs(mode_amplitude(es.E, g, (1,))))
    end
    @test amax < 50 * (a0 + 1e-6)                     # bounded — no exponential growth
end
