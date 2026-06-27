# Phase-11 2D perpendicular shock. Robust checks (noise-averaged): the y-averaged
# profile reduces to the 1D shock (frozen-in flux freezing, supercritical
# compression), and a y-perturbed run develops more transverse structure than the
# y-uniform run. Per-column x_s(y) is provided (shock_surface) but, at this test
# resolution, is noise-limited so it is not asserted tightly.

using HybridPlasmaPIC, Test, Random, Statistics

function setup2d(nx, ny, Lx, Ly, U0; nppc = 32, vthi = 0.35, seed = 1, ripple = 0.0)
    T = Float64
    sh = PerpShock2D(nx, ny, T(Lx), T(Ly); Te = 0.125, γe = 5 / 3, η = 0.02, τ = U0, B0 = 1.0)
    Np = nppc * nx * ny
    ps = ParticleSet{2,T}(Np)
    rng = MersenneTwister(seed)
    for p = 1:Np
        ps.x[1][p] = Lx * rand(rng)
        ps.x[2][p] = Ly * rand(rng)
        ux = -U0 * (1 + ripple * cos(2π * ps.x[2][p] / Ly))   # y-dependent inflow seeds rippling
        ps.v[1][p] = ux + vthi * randn(rng)
        ps.v[2][p] = vthi * randn(rng)
        ps.v[3][p] = vthi * randn(rng)
    end
    ps.weight .= shock2d_density_weight(1.0, Lx, Ly, Np)
    init_shock2d!(sh, ps)
    return sh, ps
end

# coherent m_y Fourier amplitude of the downstream-averaged Bz(y) — isolates a
# seeded ripple (cos 2πy/Ly = m_y=1) from incoherent particle noise
function ripple_mode(sh, xlo, xhi, m)
    rows = findall((sh.x .> xlo) .& (sh.x .< xhi))
    acc = 0.0 + 0.0im
    for j = 1:sh.ny
        acc += mean(@view sh.Bz[rows, j]) * cis(-2π * m * (j - 1) / sh.ny)
    end
    return abs(acc) / sh.ny
end

@testset "step_shock2d! validates timestep and subcycles before mutation" begin
    T = Float64
    sh = PerpShock2D(8, 4, T(1.0), T(1.0))
    ps = ParticleSet{2,T}(2)
    ps.x[1] .= (T(0.25), T(0.75))
    ps.x[2] .= (T(0.10), T(0.80))
    ps.v[1] .= (T(-0.1), T(-0.2))
    ps.v[2] .= (T(0.05), T(0.06))
    ps.v[3] .= (T(0.0), T(0.0))
    ps.weight .= shock2d_density_weight(T(1.0), T(1.0), T(1.0), nparticles(ps))
    init_shock2d!(sh, ps)

    x0 = ntuple(d -> copy(ps.x[d]), 2)
    v0 = ntuple(c -> copy(ps.v[c]), 3)
    B0 = copy(sh.Bz)
    n0 = copy(sh.n)

    @test_throws ArgumentError step_shock2d!(sh, ps, 0.1; NB = 0)
    @test_throws ArgumentError step_shock2d!(sh, ps, NaN; NB = 1)
    @test_throws ArgumentError step_shock2d!(sh, ps, -0.1; NB = 1)
    @test all(ps.x[d] == x0[d] for d = 1:2)
    @test all(ps.v[c] == v0[c] for c = 1:3)
    @test sh.Bz == B0
    @test sh.n == n0
end

@testset "PerpShock2D transverse FFT workspace" begin
    sh = PerpShock2D(9, 8, 1.0, 2π)
    @test sh.ywork isa FourierDerivYWorkspace
    compute_E2d!(sh)
    @test (@allocated compute_E2d!(sh)) == 0
end

@testset "Phase-11 2D perpendicular shock" begin
    T = Float64
    nx = 160
    ny = 8
    Lx = 120.0
    Ly = 12.0
    U0 = 3.0

    # y-uniform load → reduces to the 1D shock on the y-averaged profile
    sh, ps = setup2d(nx, ny, Lx, Ly, U0; nppc = 32, seed = 1, ripple = 0.0)
    for _ = 1:700
        step_shock2d!(sh, ps, 0.02; NB = 2)
    end
    @test all(isfinite, sh.Bz)
    nprof = vec(mean(sh.n; dims = 2))            # y-averaged density profile
    bzprof = vec(mean(sh.Bz; dims = 2))
    # downstream plateau between the wall buffer and the shock ramp (well inside x<20)
    dmask = (sh.x .> 2.0) .& (sh.x .< 14.0)
    n2 = mean(nprof[dmask])
    bz2 = mean(bzprof[dmask])
    @test 2.0 < n2 < 4.0                          # supercritical compression
    @test abs((bz2 / 1.0) / n2 - 1) < 0.05        # frozen-in Bz2/Bz1 = n2/n1
    m1_uniform = ripple_mode(sh, 2.0, 14.0, 1)

    # y-perturbed load → a coherent m_y=1 ripple in the downstream compression
    sh2, ps2 = setup2d(nx, ny, Lx, Ly, U0; nppc = 32, seed = 1, ripple = 0.2)
    for _ = 1:700
        step_shock2d!(sh2, ps2, 0.02; NB = 2)
    end
    @test all(isfinite, sh2.Bz)
    m1_ripple = ripple_mode(sh2, 2.0, 14.0, 1)
    @test m1_ripple > 2 * m1_uniform             # coherent ripple ≫ uniform-run m=1 noise

    # shock_surface returns a finite front for every transverse column
    xs, mxs, σs = shock_surface(sh)
    @test all(isfinite, xs) && 0 < mxs < Lx
end
