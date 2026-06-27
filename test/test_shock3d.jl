# test_shock3d.jl — 3D perpendicular hybrid shock (Phase 11, 3-D).
# Tolerances are the MEASURED behavior (see comments), not guessed.

using HybridPlasmaPIC, Test
using HybridPlasmaPIC: _curl3d!, _fourier_d!, _sbp_dx3d!, _b_rhs3d!

@testset "step_shock3d! validates timestep and subcycles before mutation" begin
    T = Float64
    sh = PerpShock3D(6, 4, 4, T(1.0), T(1.0), T(1.0))
    ps = ParticleSet{3,T}(2)
    ps.x[1] .= (T(0.25), T(0.75))
    ps.x[2] .= (T(0.10), T(0.80))
    ps.x[3] .= (T(0.20), T(0.60))
    ps.v[1] .= (T(-0.1), T(-0.2))
    ps.v[2] .= (T(0.05), T(0.06))
    ps.v[3] .= (T(0.01), T(0.02))
    ps.weight .= shock3d_density_weight(T(1.0), T(1.0), T(1.0), T(1.0), nparticles(ps))
    init_shock3d!(sh, ps)

    x0 = ntuple(d -> copy(ps.x[d]), 3)
    v0 = ntuple(c -> copy(ps.v[c]), 3)
    B0 = ntuple(c -> copy(sh.B[c]), 3)
    n0 = copy(sh.n)

    @test_throws ArgumentError step_shock3d!(sh, ps, 0.1; NB = 0)
    @test_throws ArgumentError step_shock3d!(sh, ps, NaN; NB = 1)
    @test_throws ArgumentError step_shock3d!(sh, ps, -0.1; NB = 1)
    @test all(ps.x[d] == x0[d] for d = 1:3)
    @test all(ps.v[c] == v0[c] for c = 1:3)
    @test all(sh.B[c] == B0[c] for c = 1:3)
    @test sh.n == n0
end

@testset "Shock3D operators & div-B" begin
    nx, ny, nz = 24, 16, 12
    Lx, Ly, Lz = 3.0, 2.0, 2.5
    sh = PerpShock3D(nx, ny, nz, Lx, Ly, Lz)
    x, y, z = sh.x, sh.y, sh.z

    # Fourier derivatives are spectrally exact on a band-limited mode (~1e-14).
    f = [sin(2π * y[j] / Ly) * cos(2 * 2π * z[k] / Lz) for i = 1:nx, j = 1:ny, k = 1:nz]
    dyf = similar(f)
    _fourier_d!(dyf, f, Ly, 2)
    dy_exact =
        [(2π / Ly) * cos(2π * y[j] / Ly) * cos(2 * 2π * z[k] / Lz) for i = 1:nx, j = 1:ny, k = 1:nz]
    @test maximum(abs, dyf .- dy_exact) < 1e-10
    dzf = similar(f)
    _fourier_d!(dzf, f, Lz, 3)
    dz_exact = [
        sin(2π * y[j] / Ly) * (-2 * 2π / Lz) * sin(2 * 2π * z[k] / Lz) for i = 1:nx, j = 1:ny,
        k = 1:nz
    ]
    @test maximum(abs, dzf .- dz_exact) < 1e-10

    # SBP first derivative is EXACT for a linear field (central + one-sided).
    g = [2.0 * x[i] + 1.0 for i = 1:nx, j = 1:ny, k = 1:nz]
    dxg = similar(g)
    _sbp_dx3d!(dxg, g, sh.sbp)
    @test maximum(abs, dxg .- 2.0) < 1e-10

    # div(curl E) = 0 to machine precision: the SBP-x and Fourier-y,z operators
    # act on independent indices ⇒ commute exactly ⇒ the induction preserves ∇·B.
    E1 = [
        (0.3sin(2π * y[j] / Ly) + 0.2cos(2 * 2π * z[k] / Lz)) * (1 + 0.4cos(π * x[i] / Lx)) for
        i = 1:nx, j = 1:ny, k = 1:nz
    ]
    E2 = [
        (0.5cos(2 * 2π * y[j] / Ly)) * (1 + 0.3sin(2π * x[i] / Lx)) for i = 1:nx, j = 1:ny, k = 1:nz
    ]
    E3 = [
        (0.4sin(2π * z[k] / Lz) * cos(2π * y[j] / Ly)) * (1 + 0.2cos(2π * x[i] / Lx)) for
        i = 1:nx, j = 1:ny, k = 1:nz
    ]
    Ev = (E1, E2, E3)
    Jv = ntuple(_ -> zeros(nx, ny, nz), 3)
    _curl3d!(Jv, Ev, sh)
    divJ = zeros(nx, ny, nz)
    tmp = zeros(nx, ny, nz)
    _sbp_dx3d!(divJ, Jv[1], sh.sbp)
    _fourier_d!(tmp, Jv[2], Ly, 2)
    divJ .+= tmp
    _fourier_d!(tmp, Jv[3], Lz, 3)
    divJ .+= tmp
    @test maximum(abs, divJ) < 1e-9
end

@testset "Shock3D reduces to 1-D when transversally uniform" begin
    # A y,z-uniform state must yield a y,z-uniform ∂t B (machine zero variation).
    n = 32
    sh = PerpShock3D(n, 8, 8, 60.0, 8.0, 8.0)
    for i = 1:n, j = 1:8, k = 1:8
        sh.B[3][i, j, k] = 1.0 + 0.5exp(-((sh.x[i] - 20) / 5)^2)
        sh.u[1][i, j, k] = -2.0
    end
    K = ntuple(_ -> zeros(n, 8, 8), 3)
    _b_rhs3d!(K, sh.B, sh)
    for c = 1:3
        @test maximum(abs, diff(K[c]; dims = 2)) < 1e-12
        @test maximum(abs, diff(K[c]; dims = 3)) < 1e-12
    end
end

@testset "Shock3D physics: compression, flux-freezing, div-B, stability" begin
    @test_throws ArgumentError run_perp_shock3d(;
        MA = 3.0,
        nx = 2,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = 1,
        dt = 0.03,
        seed = 1,
    )
    @test_throws ArgumentError run_perp_shock3d(;
        MA = 3.0,
        nx = 8,
        ny = 0,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = 1,
        dt = 0.03,
        seed = 1,
    )
    @test_throws ArgumentError run_perp_shock3d(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 0,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = 1,
        dt = 0.03,
        seed = 1,
    )
    @test_throws ArgumentError run_perp_shock3d(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 0,
        nsteps = 1,
        dt = 0.03,
        seed = 1,
    )
    @test_throws ArgumentError run_perp_shock3d(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = -1,
        nsteps = 1,
        dt = 0.03,
        seed = 1,
    )
    @test_throws ArgumentError run_perp_shock3d(;
        MA = 3.0,
        nx = 8,
        ny = 4,
        nz = 4,
        Lx = 12.0,
        Ly = 4.0,
        Lz = 4.0,
        nppc = 1,
        nsteps = -1,
        dt = 0.03,
        seed = 1,
    )

    # Modest grid for a fast but physically converged perpendicular shock at MA=3.
    r = run_perp_shock3d(;
        MA = 3.0,
        nx = 40,
        ny = 8,
        nz = 8,
        Lx = 70.0,
        Ly = 10.0,
        Lz = 10.0,
        nppc = 8,
        nsteps = 500,
        dt = 0.03,
        seed = 1,
    )
    @test all(isfinite, r.sh.B[3])                 # stable
    @test isfinite(r.n2) && r.n2 > 1.8             # clear compression (measured ≈2.9)
    @test r.n2 < 1.1 * r.X_rh                       # kinetic ≤ fluid RH (measured 2.9 < 3.19)
    @test 0.9 < r.frozen_ratio < 1.1               # flux freezing (measured 0.989)
    @test isfinite(r.sigma_xs) && r.sigma_xs >= 0  # transverse ripple measured

    # ∇·B is machine-zero in the interior (the SAT artifact is confined to the
    # inflow boundary cells; the induction cannot propagate it inward).
    db = magnetic_divergence3d(r.sh)
    interior = @view db[1:r.sh.nx-2, :, :]
    @test maximum(abs, interior) < 1e-8
end
