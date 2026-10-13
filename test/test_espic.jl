# Phase-12 full-PIC (electrostatic): cold Langmuir oscillation ω = ω_pe, the
# two-stream instability growth rate γ = ω_pe/(2√2), and a falsifiable null
# (kv0 > ω_pe ⇒ no growth). Oracles are closed-form.

using HybridPlasmaPIC, FFTW, Test, Random

# This oracle expands real and imaginary components independently of the
# production complex projection and checks exact rounding cells without a
# Rational-to-float conversion in the expected-value path.
module SpectralProjectionOracle
using Test
exact(x::AbstractFloat) = Rational{BigInt}(x)
exact(z::Complex) = complex(exact(real(z)), exact(imag(z)))
function projection(a, k, source, divisor)
    kr, ar = map(exact, k), map(exact, a)
    norm = sum(x * x for x in kr)
    iszero(norm) && return ar
    return ntuple(length(k)) do d
        x =
            real(ar[d]) +
            kr[d] * (
                imag(exact(source)) / exact(divisor) -
                sum(kr[j] * real(ar[j]) for j in eachindex(kr))
            ) / norm
        y =
            imag(ar[d]) -
            kr[d] * (
                real(exact(source)) / exact(divisor) +
                sum(kr[j] * imag(ar[j]) for j in eachindex(kr))
            ) / norm
        complex(x, y)
    end
end
function in_rounding_cell(actual::T, expected) where {T}
    if isinf(actual)
        largest = exact(floatmax(T))
        threshold = largest + (largest - exact(prevfloat(floatmax(T)))) / 2
        return actual > 0 ? expected >= threshold : expected <= -threshold
    end
    isfinite(actual) || return false
    here = exact(actual)
    previous, following = prevfloat(actual), nextfloat(actual)
    low = isfinite(previous) ? (exact(previous) + here) / 2 : here - (exact(following) - here) / 2
    high = isfinite(following) ? (exact(following) + here) / 2 : here + (here - exact(previous)) / 2
    even = T === Float32 ? iseven(reinterpret(UInt32, actual)) : iseven(reinterpret(UInt64, actual))
    return low < expected < high || even && (expected == low || expected == high)
end
function check(actual, expected)
    for (value, reference) in zip(actual, expected)
        for (v, r) in ((real(value), real(reference)), (imag(value), imag(reference)))
            @test in_rounding_cell(v, r)
            if iszero(v) && !iszero(r)
                @test signbit(v) == (r < 0)
            end
        end
    end
end
end

function _spectral_projection_test_allocations(a, k, source, divisor)
    HybridPlasmaPIC._spectral_project_mode(a, k, source, divisor)
    HybridPlasmaPIC._spectral_poisson_mode(k, source)
    return (
        @allocated(HybridPlasmaPIC._spectral_project_mode(a, k, source, divisor)),
        @allocated(HybridPlasmaPIC._spectral_poisson_mode(k, source))
    )
end

@testset "Longitudinal projection complete-expression rounding" begin
    H = HybridPlasmaPIC
    oracle = SpectralProjectionOracle
    rng = MersenneTwister(531789)
    for T in (Float32, Float64), D = 1:3
        emin, emax = T === Float32 ? (-149, 125) : (-1074, 1021)
        for iteration = 1:60
            sample() = ldexp(
                T(rand(rng, (-1.25, -1.0, -0.75, 0.0, 0.75, 1.0, 1.25))),
                rand(rng, emin:emax),
            )
            # One square must lose range, independently of the production
            # selector. Other axes retain their own unrelated exponent scales.
            extreme_exponent =
                T === Float32 ? (isodd(iteration) ? -120 : 100) : (isodd(iteration) ? -1000 : 800)
            k = ntuple(
                d ->
                    d == 1 ? ldexp(T(isodd(iteration) ? -1.25 : 1.25), extreme_exponent) : sample(),
                D,
            )
            source = complex(sample(), sample())
            a = ntuple(_ -> complex(sample(), sample()), D)
            divisor = sample()
            iszero(divisor) && (divisor = one(T))
            actual = H._spectral_project_mode(a, k, source, divisor)
            oracle.check(actual, oracle.projection(a, k, source, divisor))
            zeros_mode = ntuple(_ -> zero(Complex{T}), D)
            oracle.check(
                H._spectral_poisson_mode(k, source),
                oracle.projection(zeros_mode, k, source, one(T)),
            )
        end
        # Ordinary helpers remain allocation-free; step allocations are
        # checked separately because the surrounding PIC path owns workspaces.
        k = ntuple(d -> T(d), D)
        a = ntuple(d -> complex(T(d) / 3, T(d) / 7), D)
        @test _spectral_projection_test_allocations(a, k, complex(T(0.1), T(0.2)), T(0.3)) == (0, 0)
        null = ntuple(_ -> zero(T), D)
        @test H._spectral_project_mode(a, null, complex(T(3)), one(T)) === a
        @test H._spectral_poisson_mode(null, complex(T(3))) == ntuple(_ -> zero(Complex{T}), D)
    end
    for T in (Float32, Float64)
        tiny = nextfloat(zero(T))
        for sign in (-one(T), one(T)),
            source in (complex(sign * tiny, zero(T)), complex(zero(T), sign * tiny))

            k = (T(2),)
            SpectralProjectionOracle.check(
                H._spectral_poisson_mode(k, source),
                SpectralProjectionOracle.projection((zero(Complex{T}),), k, source, one(T)),
            )
        end
        # Binade midpoint ties arise after cancellation even though each
        # supplied field component is exactly representable.
        for a in ((complex(T(2)), complex(-eps(T))), (complex(prevfloat(T(4))), complex(eps(T))))
            k = (tiny, tiny)
            SpectralProjectionOracle.check(
                H._spectral_project_mode(a, k, zero(Complex{T}), one(T)),
                SpectralProjectionOracle.projection(a, k, zero(Complex{T}), one(T)),
            )
        end
    end
end

@testset "Electrostatic Poisson preserves nonzero extreme wave numbers" begin
    for T in (Float32, Float64), D = 1:3, axis = 1:D, large_domain in (false, true)
        scale =
            T === Float32 ? (large_domain ? T(1e25) : T(1e-25)) :
            (large_domain ? T(1e200) : T(1e-200))
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(d -> d == axis ? T(2pi) * scale : T(4), D))
        electrons = ParticleSet{D,T}(4^(D - 1); q = -one(T), m = one(T))
        for (p, I) in enumerate(CartesianIndices(ntuple(d -> d == axis ? 1 : 4, D)))
            for d = 1:D
                electrons.x[d][p] = T(I[d] - 1) * g.dx[d]
            end
            electrons.weight[p] = prod(g.dx)
        end
        es = ElectrostaticPIC(g, length(electrons.weight); n0 = T(0.25))
        init_espic!(es, electrons)
        # Plane rho=(-3/4,1/4,1/4,1/4) has only the represented +/-1
        # coefficients -1 and the discarded Nyquist coefficient -1.
        amplitude = setprecision(BigFloat, 256) do
            T(-inv(2BigFloat(g.kvec[axis][2])))
        end
        for I in CartesianIndices(es.E[axis])
            @test es.E[axis][I] ≈ amplitude * (zero(T), one(T), zero(T), -one(T))[I[axis]] rtol =
                16eps(T) atol = zero(T)
        end
        @test all(d -> d == axis || all(iszero, es.E[d]), 1:3)
    end
    for T in (Float32, Float64), overflow in (false, true)
        kscale = T === Float32 ? (overflow ? T(1e15) : T(1e-15)) : (overflow ? T(1e150) : T(1e-150))
        amplitude =
            T === Float32 ? (overflow ? T(1e30) : T(1e-35)) : (overflow ? T(1e200) : T(1e-200))
        g = FourierGrid((4,), (T(2pi) / kscale,))
        es = ElectrostaticPIC(g, 0; n0 = zero(T))
        es.ne .= amplitude .* T[1, 0, 0, 0]
        poisson_E!(es)
        expected = setprecision(BigFloat, 256) do
            T(-BigFloat(amplitude) / (2BigFloat(g.kvec[1][2])))
        end
        @test es.E[1][2] ≈ expected rtol = 16eps(T)
    end
end


@testset "Electrostatic1D parameter validation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    es = Electrostatic1D(g, 4; n0 = 1.0)
    @test es.n0 == 1.0
    @test_throws ArgumentError Electrostatic1D(g, 4; n0 = NaN)
    @test_throws ArgumentError Electrostatic1D(g, 4; n0 = -1.0)
    @test_throws ArgumentError Electrostatic1D(g, -1; n0 = 1.0)
end

@testset "ElectrostaticPIC constructor validation" begin
    g2 = FourierGrid((8, 6), (2π, 3π))
    es2 = ElectrostaticPIC(g2, 5; n0 = 0.75)
    @test es2.n0 == 0.75
    @test size(es2.ne) == (8, 6)
    @test all(size(Ec) == (8, 6) for Ec in es2.E)
    @test all(length(Ep) == 5 for Ep in es2.Ep)

    g3 = FourierGrid((6, 5, 4), (2π, 3π, 4π))
    es3 = ElectrostaticPIC(g3, 7; n0 = 1.25)
    @test size(es3.ne) == (6, 5, 4)
    @test all(size(Ec) == (6, 5, 4) for Ec in es3.E)
    @test all(length(Ep) == 7 for Ep in es3.Ep)

    @test_throws ArgumentError ElectrostaticPIC(g2, 5; n0 = NaN)
    @test_throws ArgumentError ElectrostaticPIC(g2, 5; n0 = -1.0)
    @test_throws ArgumentError ElectrostaticPIC(g2, -1; n0 = 1.0)

    g4 = FourierGrid((2, 2, 2, 2), (1.0, 1.0, 1.0, 1.0))
    @test_throws ArgumentError ElectrostaticPIC(g4, 1; n0 = 1.0)
end

@testset "Electrostatic PIC resizes particle buffers after particle count changes" begin
    T = Float64

    g1 = FourierGrid((8,), (2π,))
    e1 = ParticleSet{1,T}(2; q = -1.0, m = 1.0)
    load_lattice_1d!(e1, 0.0, 2π)
    set_density_weight!(e1, 1.0, g1)
    es1 = Electrostatic1D(g1, 1; n0 = 1.0)
    init_espic!(es1, e1)
    @test length(es1.Ep) == nparticles(e1)
    @test step_espic!(es1, e1, 0.01) === es1

    g2 = FourierGrid((8, 8), (2π, 2π))
    e2 = ParticleSet{2,T}(4; q = -1.0, m = 1.0)
    load_lattice!(e2, (0.0, 0.0), g2.L, (2, 2))
    set_density_weight!(e2, 1.0, g2)
    es2 = ElectrostaticPIC(g2, 1; n0 = 1.0)
    init_espic!(es2, e2)
    @test all(length(Ep) == nparticles(e2) for Ep in es2.Ep)
    @test step_espic!(es2, e2, 0.01) === es2
end

@testset "ElectrostaticPIC 2D spectral Poisson oracle" begin
    nx, ny = 32, 24
    Lx, Ly = 2π, 2π
    g = FourierGrid((nx, ny), (Lx, Ly))
    es = ElectrostaticPIC(g, 0; n0 = 1.25)
    A = 0.2
    mx, my = 2, 3
    kx = 2π * mx / Lx
    ky = 2π * my / Ly
    k2 = kx^2 + ky^2
    ex = zeros(Float64, nx, ny)
    ey = zeros(Float64, nx, ny)

    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        y = (j - 1) * g.dx[2]
        rho = A * cos(kx * x) * cos(ky * y)
        es.ne[i, j] = es.n0 - rho
        ex[i, j] = A * kx / k2 * sin(kx * x) * cos(ky * y)
        ey[i, j] = A * ky / k2 * cos(kx * x) * sin(ky * y)
    end

    poisson_E!(es)
    @test maximum(abs, es.E[1] .- ex) < 1e-12
    @test maximum(abs, es.E[2] .- ey) < 1e-12
    @test maximum(abs, es.E[3]) < 1e-12
    @test field_energy(es) ≈ 0.5 * (sum(abs2, es.E[1]) + sum(abs2, es.E[2])) * prod(g.dx)
end

@testset "ElectrostaticPIC 3D spectral Poisson oracle" begin
    nx, ny, nz = 16, 12, 10
    Lx, Ly, Lz = 2π, 2π, 2π
    g = FourierGrid((nx, ny, nz), (Lx, Ly, Lz))
    es = ElectrostaticPIC(g, 0; n0 = 0.5)
    A = 0.15
    mx, my, mz = 1, 2, 3
    kx = 2π * mx / Lx
    ky = 2π * my / Ly
    kz = 2π * mz / Lz
    k2 = kx^2 + ky^2 + kz^2
    ex = zeros(Float64, nx, ny, nz)
    ey = zeros(Float64, nx, ny, nz)
    ez = zeros(Float64, nx, ny, nz)

    for k = 1:nz, j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        y = (j - 1) * g.dx[2]
        z = (k - 1) * g.dx[3]
        cx = cos(kx * x)
        cy = cos(ky * y)
        cz = cos(kz * z)
        sx = sin(kx * x)
        sy = sin(ky * y)
        sz = sin(kz * z)
        rho = A * cx * cy * cz
        es.ne[i, j, k] = es.n0 - rho
        ex[i, j, k] = A * kx / k2 * sx * cy * cz
        ey[i, j, k] = A * ky / k2 * cx * sy * cz
        ez[i, j, k] = A * kz / k2 * cx * cy * sz
    end

    poisson_E!(es)
    @test maximum(abs, es.E[1] .- ex) < 1e-12
    @test maximum(abs, es.E[2] .- ey) < 1e-12
    @test maximum(abs, es.E[3] .- ez) < 1e-12
end

@testset "ElectrostaticPIC 2D electron species validation" begin
    T = Float64
    g = FourierGrid((8, 8), (2π, 2π))

    badq = ParticleSet{2,T}(4; q = 0.0, m = 1.0)
    load_lattice!(badq, (0.0, 0.0), (2π, 2π), (2, 2))
    set_density_weight!(badq, 1.0, g)
    es = ElectrostaticPIC(g, 4; n0 = 1.0)
    @test_throws ArgumentError init_espic!(es, badq)
    for c = 1:3
        fill!(es.E[c], c)
    end
    fill!(es.ne, 0.25)
    x0 = ntuple(d -> copy(badq.x[d]), 2)
    v0 = ntuple(c -> copy(badq.v[c]), 3)
    E0 = ntuple(c -> copy(es.E[c]), 3)
    ne0 = copy(es.ne)
    @test_throws ArgumentError step_espic!(es, badq, 0.1)
    @test all(badq.x[d] == x0[d] for d = 1:2)
    @test all(badq.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test es.ne == ne0

    badm = ParticleSet{2,T}(4; q = -1.0, m = 2.0)
    load_lattice!(badm, (0.0, 0.0), (2π, 2π), (2, 2))
    set_density_weight!(badm, 1.0, g)
    @test_throws ArgumentError init_espic!(ElectrostaticPIC(g, 4; n0 = 1.0), badm)
end

@testset "ElectrostaticPIC 2D timestep validation before mutation" begin
    T = Float64
    g = FourierGrid((8, 8), (2π, 2π))
    e = ParticleSet{2,T}(64; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), (2π, 2π), (8, 8))
    set_density_weight!(e, 1.0, g)
    es = ElectrostaticPIC(g, 64; n0 = 1.0)
    init_espic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 2)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    E0 = ntuple(c -> copy(es.E[c]), 3)
    ne0 = copy(es.ne)

    @test_throws ArgumentError step_espic!(es, e, NaN)
    @test_throws ArgumentError step_espic!(es, e, -0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test es.ne == ne0
end

@testset "ElectrostaticPIC 2D uniform equilibrium is force-free" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 3π))
    counts = (8, 6)
    N = prod(counts)
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, counts)
    set_density_weight!(e, 1.0, g)
    es = ElectrostaticPIC(g, N; n0 = 1.0)
    init_espic!(es, e)

    @test maximum(abs, es.ne .- 1.0) < 1e-12
    @test all(maximum(abs, es.E[c]) < 1e-12 for c = 1:3)
    x0 = ntuple(d -> copy(e.x[d]), 2)
    for _ = 1:5
        step_espic!(es, e, 0.1)
    end
    @test maximum(abs, es.ne .- 1.0) < 1e-12
    @test all(maximum(abs, es.E[c]) < 1e-12 for c = 1:3)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test field_energy(es) < 1e-24
end

@testset "Electrostatic1D electron species validation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))

    badq = ParticleSet{1,T}(8; q = 0.0, m = 1.0)
    load_lattice_1d!(badq, 0.0, 2π)
    set_density_weight!(badq, 1.0, g)
    es = Electrostatic1D(g, 8; n0 = 1.0)
    @test_throws ArgumentError init_espic!(es, badq)
    fill!(es.E, 1.0)
    x0 = ntuple(d -> copy(badq.x[d]), 1)
    v0 = ntuple(c -> copy(badq.v[c]), 3)
    E0 = copy(es.E)
    ne0 = copy(es.ne)
    @test_throws ArgumentError step_espic!(es, badq, 0.1)
    @test all(badq.x[d] == x0[d] for d = 1:1)
    @test all(badq.v[c] == v0[c] for c = 1:3)
    @test es.E == E0
    @test es.ne == ne0

    badm = ParticleSet{1,T}(8; q = -1.0, m = 2.0)
    load_lattice_1d!(badm, 0.0, 2π)
    set_density_weight!(badm, 1.0, g)
    @test_throws ArgumentError init_espic!(Electrostatic1D(g, 8; n0 = 1.0), badm)
end

@testset "step_espic! validates timestep before mutation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    e = ParticleSet{1,T}(8; q = -1.0, m = 1.0)
    load_lattice_1d!(e, 0.0, 2π)
    set_density_weight!(e, 1.0, g)
    es = Electrostatic1D(g, 8; n0 = 1.0)
    init_espic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 1)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    E0 = copy(es.E)
    ne0 = copy(es.ne)

    @test_throws ArgumentError step_espic!(es, e, NaN)
    @test_throws ArgumentError step_espic!(es, e, -0.1)
    @test all(e.x[d] == x0[d] for d = 1:1)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test es.E == E0
    @test es.ne == ne0
end

# dominant positive-frequency peak (parabolic interpolation)
function _peakfreq(series, dt)
    Nt = length(series)
    S = fft(series)
    mag = abs.(S)
    best = 2
    bm = -1.0
    for j = 2:(Nt÷2)
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
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
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
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
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
    e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
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

@testset "ESPIC-order: step_espic! is 2nd-order in dt (leapfrog priming)" begin
    # The loaded velocity is physical v^0; step_espic! primes it once to v^{-1/2} so the
    # Boris+leapfrog scheme is 2nd-order (was 1st-order unprimed). Cold-Langmuir self-
    # convergence on the fundamental mode amplitude — deterministic (lattice load), so the
    # differences are pure temporal error; the order must sit near 2.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.01
    function amp(dt, Tf; lead_zero)
        n = 32
        N = 200 * n
        g = FourierGrid((n,), (T(L),))
        e = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
        load_lattice_1d!(e, 0.0, T(L))
        set_density_weight!(e, 1.0, g)
        for p = 1:N
            e.v[1][p] = 0.0
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
        end
        es = Electrostatic1D(g, N; n0 = 1.0)
        init_espic!(es, e)
        lead_zero && step_espic!(es, e, 0.0)          # dt=0 must not consume the priming
        for _ = 1:round(Int, Tf / dt)
            step_espic!(es, e, dt)
        end
        abs(mode_amplitude(es.E, g, (1,)))
    end
    Tf = 2.0
    dts = (0.02, 0.01, 0.005, 0.0025)
    for lz in (false, true)
        v = [amp(dt, Tf; lead_zero = lz) for dt in dts]
        r1 = log2(abs(v[1] - v[2]) / abs(v[2] - v[3]))
        r2 = log2(abs(v[2] - v[3]) / abs(v[3] - v[4]))
        @test r1 > 1.6        # ≈2.0 with priming; ≈1.0 (would fail) without it
        @test r2 > 1.6        # lead_zero=true also ≈2.0 (dt=0 short-circuit preserves priming)
    end
end

@testset "Electrostatic PIC preserves velocity staggering across timestep changes" begin
    # A zero-weight tracer sees a prescribed constant field, giving an exact
    # acceleration oracle independent of deposition and the Poisson solver.
    for (constructor, D) in
        ((Electrostatic1D, 1), (ElectrostaticPIC, 1), (ElectrostaticPIC, 2), (ElectrostaticPIC, 3))
        for dts in ((0.1, 0.2, 0.05), (0.1, 0.1, 0.1))
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0, D))
            e = ParticleSet{D,Float64}(1; q = -1.0, m = 1.0)
            e.weight .= 0.0
            for d = 1:D
                e.x[d] .= 1.0
            end
            es = constructor(g, 1; n0 = 0.0)
            init_espic!(es, e)
            acceleration = (0.7, -0.3, 0.2)
            function set_field!()
                if es isa Electrostatic1D
                    es.E .= -acceleration[1]
                else
                    for d = 1:D
                        es.E[d] .= -acceleration[d]
                    end
                end
            end
            elapsed = 0.0
            for dt in dts
                previous_dt = es.last_dt[]
                step_espic!(es, e, 0.0)
                @test es.last_dt[] == previous_dt
                set_field!()
                step_espic!(es, e, dt)
                elapsed += dt
                for d = 1:D
                    @test e.x[d][1] ≈ 1 + acceleration[d] * elapsed^2 / 2 atol = 2e-15
                    @test e.v[d][1] ≈ acceleration[d] * (elapsed - dt / 2) atol = 2e-15
                end
            end
            before = (deepcopy(e.x), deepcopy(e.v), es.last_dt[])
            for dt in (NaN, Inf, -0.1)
                @test_throws ArgumentError step_espic!(es, e, dt)
                @test (e.x, e.v, es.last_dt[]) == before
            end
            # Reload physical velocities and reinitialize: old timestep history
            # must not influence the next run's first half-step prime.
            for c = 1:3
                e.v[c] .= 0.0
            end
            init_espic!(es, e)
            @test es.last_dt[] == 0.0
            set_field!()
            step_espic!(es, e, 0.08)
            for d = 1:D
                @test e.x[d][1] ≈ before[1][d][1] + acceleration[d] * 0.08^2 / 2 atol = 2e-15
            end
        end
    end
end

@testset "Electrostatic PIC changing-timestep Langmuir convergence" begin
    function variable_dt_amplitude(constructor, refinement)
        g = FourierGrid((32,), (2π,))
        e = ParticleSet{1,Float64}(2048; q = -1.0, m = 1.0)
        load_lattice_1d!(e, 0.0, g.L[1])
        set_density_weight!(e, 1.0, g)
        @. e.x[1] -= 0.02 * sin(e.x[1])
        es = constructor(g, nparticles(e))
        init_espic!(es, e)
        # Change step at the same physical time under each refinement.
        for (count, dt) in ((10refinement, 0.04 / refinement), (20refinement, 0.02 / refinement))
            for _ = 1:count
                step_espic!(es, e, dt)
            end
        end
        field = es isa Electrostatic1D ? es.E : es.E[1]
        return abs(mode_amplitude(field, g, (1,)))
    end
    for constructor in (Electrostatic1D, ElectrostaticPIC)
        values = [variable_dt_amplitude(constructor, r) for r in (1, 2, 4)]
        @test log2(abs(values[1] - values[2]) / abs(values[2] - values[3])) > 1.7
    end
end
