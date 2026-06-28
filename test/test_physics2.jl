# Core physics extensions: spectral Laplacian, hyperresistivity term in Ohm's
# law, multi-species charge-weighted moments, and the adiabatic electron-pressure
# equation — each vs an analytic oracle (band-limited ⇒ spectrally exact).

using HybridPlasmaPIC, Test, Random, Statistics

@testset "spectral Laplacian ∇²(sin kx) = −k² sin kx" begin
    T = Float64
    for (nc, modes) in (((32,), (3,)), ((16, 16), (2, 3)), ((8, 8, 8), (1, 2, 1)))
        D = length(nc)
        L = ntuple(_ -> T(2π), D)
        g = FourierGrid(nc, L)
        f = Array{T,D}(undef, nc)
        ex = similar(f)
        ksum = sum(abs2, modes)
        for I in CartesianIndices(f)
            t = Tuple(I)
            val = one(T)
            for d = 1:D
                val *= sin(modes[d] * (t[d] - 1) * g.dx[d])
            end
            f[I] = val
            ex[I] = -ksum * val
        end
        out = similar(f)
        laplacian!(out, f, g)
        @test maximum(abs, out .- ex) < 1e-10
    end
end

@testset "hyperresistivity term −ηH ∇²J (1D analytic)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    f = HybridFields{1,T}((n,))
    fill!(f.n, 1.0)
    x = [(i - 1) * g.dx[1] for i = 1:n]
    k = 3.0
    f.B[3] .= cos.(k .* x)                      # Bz=cos(kx) ⇒ Jy=k sin(kx), ∇²Jy=−k³ sin(kx)
    ηH = 0.05
    ohms_law!(f, HybridModel(IsothermalElectrons(0.0); ηH = ηH), g)
    @test maximum(abs, f.E[1] .- (@. k * sin(k * x) * cos(k * x))) < 1e-9   # Hall (unchanged)
    @test maximum(abs, f.E[2] .- (@. ηH * k^3 * sin(k * x))) < 1e-8         # +ηH k³ sin(kx)
end

@testset "multi-species charge-weighted moments" begin
    T = Float64
    n = 32
    L = 2π
    g = FourierGrid((n,), (L,))
    nA = 1.0
    nB = 0.5
    uA = 0.4
    uB = -0.8
    NA = 64 * n
    NB = 64 * n
    psA = ParticleSet{1,T}(NA)
    load_lattice_1d!(psA, 0.0, L)
    set_density_weight!(psA, nA, g)
    load_quiet_velocities!(psA, MersenneTwister(1), (uA, 0.0, 0.0), (0.0, 0.0, 0.0))
    psB = ParticleSet{1,T}(NB)
    psB.q = 2.0
    load_lattice_1d!(psB, 0.0, L)
    set_density_weight!(psB, nB, g)
    load_quiet_velocities!(psB, MersenneTwister(2), (uB, 0.0, 0.0), (0.0, 0.0, 0.0))
    species = [psA, psB]
    f = HybridFields{1,T}((n,))
    ntmp = zeros(T, n)
    mtmp = ntuple(_ -> zeros(T, n), 3)
    works = [Vector{T}(undef, nparticles(psA)), Vector{T}(undef, nparticles(psB))]
    compute_moments_multi!(f, species, g, CIC(), 1e-6; ntmp, mtmp, works)
    charge_density = psA.q * nA + psB.q * nB
    charge_weighted_u = (psA.q * nA * uA + psB.q * nB * uB) / charge_density
    @test isapprox(mean(f.n), charge_density; rtol = 0.02)                       # Σ q_s n_s
    @test isapprox(mean(f.ui[1]), charge_weighted_u; rtol = 0.02)
    @test (@allocated compute_moments_multi!(f, species, g, CIC(), 1e-6; ntmp, mtmp, works)) <= 128
    @test_throws DimensionMismatch compute_moments_multi!(
        f,
        species,
        g,
        CIC(),
        1e-6;
        ntmp = zeros(T, n - 1),
        mtmp,
        works,
    )
    @test_throws DimensionMismatch compute_moments_multi!(
        f,
        species,
        g,
        CIC(),
        1e-6;
        ntmp,
        mtmp,
        works = works[1:1],
    )
    bad_works = [Vector{T}(undef, nparticles(psA) - 1), Vector{T}(undef, nparticles(psB))]
    @test_throws DimensionMismatch compute_moments_multi!(
        f,
        species,
        g,
        CIC(),
        1e-6;
        ntmp,
        mtmp,
        works = bad_works,
    )
end

@testset "compute_moments_multi! validates density floor before mutation" begin
    T = Float64
    g = FourierGrid((8,), (2π,))
    species = [ParticleSet{1,T}(0)]
    f = HybridFields{1,T}((8,))
    ntmp = zeros(T, 8)
    mtmp = ntuple(_ -> zeros(T, 8), 3)
    works = [Vector{T}(undef, 0)]
    fill!(f.n, 1.0)
    fill!(f.ui[1], 2.0)
    fill!(f.ui[2], 3.0)
    fill!(f.ui[3], 4.0)

    @test_throws ArgumentError compute_moments_multi!(
        f,
        species,
        g,
        NGP(),
        0.0;
        ntmp,
        mtmp,
        works,
    )
    @test all(==(1.0), f.n)
    @test all(==(2.0), f.ui[1])
    @test all(==(3.0), f.ui[2])
    @test all(==(4.0), f.ui[3])

    @test_throws ArgumentError compute_moments_multi!(
        f,
        species,
        g,
        NGP(),
        NaN;
        ntmp,
        mtmp,
        works,
    )
    @test all(==(1.0), f.n)
    @test all(==(2.0), f.ui[1])
    @test all(==(3.0), f.ui[2])
    @test all(==(4.0), f.ui[3])
end

@testset "electron pressure equation (analytic single step)" begin
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    k = 2π / L
    γe = 5 / 3
    dt = 1e-3
    # (a) uniform p_e, compressive flow ⇒ −γe p_e ∇·u_e
    p0 = 0.7
    A = 0.1
    pe = fill(p0, n)
    ue = (A .* sin.(k .* x), zeros(T, n), zeros(T, n))
    pe0 = copy(pe)
    advance_electron_pressure!(pe, ue, dt, γe, g)
    @test maximum(abs, (pe .- pe0) .- (@. -dt * γe * p0 * A * k * cos(k * x))) < 1e-12
    # (b) advection by a uniform flow ⇒ −u_e·∇p_e
    ε = 0.05
    u0 = 0.3
    pe = @. p0 * (1 + ε * cos(k * x))
    ue = (fill(u0, n), zeros(T, n), zeros(T, n))
    pe0 = copy(pe)
    advance_electron_pressure!(pe, ue, dt, γe, g)
    @test maximum(abs, (pe .- pe0) .- (@. dt * u0 * p0 * ε * k * sin(k * x))) < 1e-12
end

@testset "electron pressure: adiabatic invariant pe/n^γe conserved (energy closure)" begin
    # Integration-level check: evolving p_e (adiabatic equation) and n (continuity
    # ∂t n = −∇·(n u)) under the SAME compressive flow must materially conserve the
    # entropy s = p_e/n^γe — the adiabatic energy closure. A wrong compression
    # coefficient breaks it (measured drift: 1.4e-5 with γe=5/3, but 3.3e-2 / 2.2e-2
    # with γe=2/3 / 1.0).
    T = Float64
    n = 64
    L = 2π
    g = FourierGrid((n,), (L,))
    x = [(i - 1) * g.dx[1] for i = 1:n]
    k = 2π / L
    γe = 5 / 3
    A = 0.4
    ux = A .* sin.(k .* x)
    ue = (ux, zeros(T, n), zeros(T, n))
    dens = fill(1.0, n)
    pe = fill(0.5, n)
    s0 = pe ./ dens .^ γe                       # uniform initial entropy
    dt = 2e-3
    nsteps = 40
    divnu = similar(dens)
    for _ = 1:nsteps
        divergence!(divnu, (dens .* ux, zeros(T, n), zeros(T, n)), g)  # ∇·(n u)
        advance_electron_pressure!(pe, ue, dt, γe, g)
        @. dens += -dt * divnu
    end
    drift = maximum(abs, (pe ./ dens .^ γe) .- s0) / maximum(abs, s0)
    @test drift < 1e-3
end
