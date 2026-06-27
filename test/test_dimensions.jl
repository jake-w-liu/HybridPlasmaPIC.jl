# DIM-001 — dimensional reduction. The D-parametric operators, field model, and
# integrator must reduce a y-invariant 2D configuration exactly to the
# corresponding 1D problem. First exercises the 2D code paths.

using HybridPlasmaPIC, Test, LinearAlgebra, Random

@testset "DIM operators: 1D ≡ y-uniform 2D (exact)" begin
    T = Float64
    nx = 16
    ny = 8
    L = 2π
    g1 = FourierGrid((nx,), (L,))
    g2 = FourierGrid((nx, ny), (L, L))
    rep(a) = repeat(a, 1, ny)                       # F(x,y) = f(x)

    f1 = randn(MersenneTwister(1), nx)
    f2 = rep(f1)

    # gradient: ∂x matches 1D, ∂y = 0
    gr1 = (zeros(T, nx),)
    gradient!(gr1, f1, g1)
    gr2 = (zeros(T, nx, ny), zeros(T, nx, ny))
    gradient!(gr2, f2, g2)
    @test maximum(abs, gr2[1] .- rep(gr1[1])) < 1e-12
    @test maximum(abs, gr2[2]) < 1e-12

    # curl of a y-uniform 3-vector field reduces to the 1D curl
    A1 = ntuple(_ -> randn(MersenneTwister(2), nx), 3)
    A2 = ntuple(c -> rep(A1[c]), 3)
    B1 = ntuple(_ -> zeros(T, nx), 3)
    curl!(B1, A1, g1)
    B2 = ntuple(_ -> zeros(T, nx, ny), 3)
    curl!(B2, A2, g2)
    for c = 1:3
        @test maximum(abs, B2[c] .- rep(B1[c])) < 1e-12
    end

    # Ohm's law reduces identically
    model = HybridModel(IsothermalElectrons(0.5))
    h1 = HybridFields{1,T}((nx,))
    h2 = HybridFields{2,T}((nx, ny))
    n1 = 1 .+ 0.1 .* randn(MersenneTwister(3), nx)
    h1.n .= n1
    h2.n .= rep(n1)
    for c = 1:3
        u = randn(MersenneTwister(10 + c), nx)
        b = randn(MersenneTwister(20 + c), nx)
        h1.ui[c] .= u
        h2.ui[c] .= rep(u)
        h1.B[c] .= b
        h2.B[c] .= rep(b)
    end
    ohms_law!(h1, model, g1)
    ohms_law!(h2, model, g2)
    for c = 1:3
        @test maximum(abs, h2.E[c] .- rep(h1.E[c])) < 1e-12
    end
end

@testset "DIM-001 2D integrator: y-invariant ion-acoustic ≡ 1D" begin
    T = Float64
    nx = 32
    ny = 4
    L = 2π
    Te = 1.0
    mx = 1
    k = 2π * mx / L
    g = FourierGrid((nx, ny), (L, L))
    counts = (nx * 8, ny * 8)               # cell-centered lattice: 64 particles/cell, y-uniform
    N = prod(counts)
    ps = ParticleSet{2,T}(N)
    load_lattice!(ps, (0.0, 0.0), (L, L), counts)
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(8), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    for p = 1:N
        ps.v[1][p] += 0.005 * sin(k * ps.x[1][p])   # x-only perturbation, no y dependence
    end
    st = HybridStepper(g, HybridModel(IsothermalElectrons(Te)), CIC(), N)  # B0=0
    init!(st, ps)
    dt = 0.02
    series = Float64[]
    ymode_max = 0.0
    for _ = 1:700
        step!(st, ps, dt)
        push!(series, real(mode_amplitude(st.fields.n, g, (mx, 0))))
        ymode_max = max(ymode_max, abs(mode_amplitude(st.fields.n, g, (0, 1))))  # spurious y-structure
    end
    @test all(isfinite, series)
    # x-mode oscillates at the 1D ion-acoustic frequency
    zc = Int[]
    for i = 2:length(series)
        series[i-1] < 0 && series[i] >= 0 && push!(zc, i)
    end
    @test length(zc) >= 2
    per = (zc[end] - zc[1]) / (length(zc) - 1) * dt
    @test abs(2π / per - k * sqrt(Te)) / (k * sqrt(Te)) < 0.04
    # y-structure stays at particle-noise level, far below the x-signal amplitude
    @test ymode_max < 0.1 * maximum(abs, series)
end

@testset "DIM-002 operators: 2D ≡ z-uniform 3D (exact)" begin
    T = Float64
    nx = 12
    ny = 10
    nz = 8
    L = 2π
    g2 = FourierGrid((nx, ny), (L, L))
    g3 = FourierGrid((nx, ny, nz), (L, L, L))
    rep(a) = repeat(reshape(a, nx, ny, 1), 1, 1, nz)     # F(x,y,z) = f(x,y)

    f2 = randn(MersenneTwister(1), nx, ny)
    f3 = rep(f2)
    gr2 = ntuple(_ -> zeros(T, nx, ny), 2)
    gradient!(gr2, f2, g2)
    gr3 = ntuple(_ -> zeros(T, nx, ny, nz), 3)
    gradient!(gr3, f3, g3)
    @test maximum(abs, gr3[1] .- rep(gr2[1])) < 1e-12
    @test maximum(abs, gr3[2] .- rep(gr2[2])) < 1e-12
    @test maximum(abs, gr3[3]) < 1e-12                   # ∂z = 0

    A2 = ntuple(_ -> randn(MersenneTwister(2), nx, ny), 3)
    A3 = ntuple(c -> rep(A2[c]), 3)
    B2 = ntuple(_ -> zeros(T, nx, ny), 3)
    curl!(B2, A2, g2)
    B3 = ntuple(_ -> zeros(T, nx, ny, nz), 3)
    curl!(B3, A3, g3)
    for c = 1:3
        @test maximum(abs, B3[c] .- rep(B2[c])) < 1e-12
    end

    model = HybridModel(IsothermalElectrons(0.5))
    h2 = HybridFields{2,T}((nx, ny))
    h3 = HybridFields{3,T}((nx, ny, nz))
    n2 = 1 .+ 0.1 .* randn(MersenneTwister(3), nx, ny)
    h2.n .= n2
    h3.n .= rep(n2)
    for c = 1:3
        u = randn(MersenneTwister(10 + c), nx, ny)
        b = randn(MersenneTwister(20 + c), nx, ny)
        h2.ui[c] .= u
        h3.ui[c] .= rep(u)
        h2.B[c] .= b
        h3.B[c] .= rep(b)
    end
    ohms_law!(h2, model, g2)
    ohms_law!(h3, model, g3)
    for c = 1:3
        @test maximum(abs, h3.E[c] .- rep(h2.E[c])) < 1e-12
    end
end
