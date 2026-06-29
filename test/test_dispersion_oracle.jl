# HYB-005 (fast/slow MHD limits) + HYB-006 (independent warm dispersion oracle).
#
# HYB-006 builds a standalone linearized warm Hall-MHD eigenvalue solver
# (test/oracles/hybrid_dispersion_oracle.jl) that shares NO code path with the
# production spatial operator, and verifies it reproduces every known analytic
# branch (parallel whistler/ion-cyclotron, perpendicular fast magnetosonic,
# oblique fast/slow/Alfvén), plus polarization (eigenvectors) and group velocity.
# HYB-005 verifies the oblique fast/slow phase speeds c_{f,s} both in the oracle
# (closed form, exact) and in an actual PIC run (perpendicular fast magnetosonic,
# the cleanly-excitable fluid-regime mode), compared in the fluid regime where the
# checklist's fast/slow relation applies.

using HybridPlasmaPIC, FFTW, Test, Random, Statistics, LinearAlgebra

include("oracles/hybrid_dispersion_oracle.jl")
using .HybridDispersionOracle

# dominant frequency of a real series via early-window neg→pos zero crossings
function _freq_zerocross(series, dt)
    zc = Int[]
    for i = 2:length(series)
        if series[i-1] < 0 && series[i] >= 0
            push!(zc, i)
        end
    end
    length(zc) >= 2 || return NaN
    per = (zc[end] - zc[1]) / (length(zc) - 1) * dt
    return 2π / per
end

@testset "HYB-006 oracle reproduces the analytic branches" begin
    # Parallel cold (B0 ∥ k = x̂): whistler ωW and ion-cyclotron ωIC, exact.
    for K in (0.5, 1.0, 2.0)
        f = dispersion_frequencies((K, 0.0, 0.0), (1.0, 0.0, 0.0); cs = 0.0)
        ωW = (sqrt(K^4 + 4K^2) + K^2) / 2
        ωIC = (sqrt(K^4 + 4K^2) - K^2) / 2
        @test any(z -> isapprox(z, ωW; rtol = 1e-6), f)
        @test any(z -> isapprox(z, ωIC; rtol = 1e-6), f)
    end
    # Perpendicular (B0 = ẑ ⊥ k = x̂): fast magnetosonic ω = k√(vA²+cs²), exact.
    for (K, cs) in ((1.0, 0.0), (1.0, 0.5), (2.0, 0.7))
        f = dispersion_frequencies((K, 0.0, 0.0), (0.0, 0.0, 1.0); cs = cs)
        @test any(z -> isapprox(z, K * sqrt(1 + cs^2); rtol = 1e-6), f)
    end
    # Oblique MHD limit (small K): the three eigenfrequencies/K → c_fast, |cosθ| (Alfvén), c_slow.
    let K = 0.02, cs = 0.6
        for θ in (π / 6, π / 4, π / 3)
            B0 = (cos(θ), 0.0, sin(θ))
            f = dispersion_frequencies((K, 0.0, 0.0), B0; cs = cs)
            cf, csl = fast_slow_speeds(1.0, cs, θ)
            speeds = sort(f ./ K)
            @test isapprox(speeds[1], csl; rtol = 2e-3)             # slow
            @test isapprox(speeds[2], abs(cos(θ)); rtol = 2e-3)     # intermediate (Alfvén)
            @test isapprox(speeds[3], cf; rtol = 2e-3)              # fast
        end
    end
end

@testset "HYB-003 long-wavelength Alfvén limit (kd_i ≪ 1)" begin
    # As K = k d_i → 0 the parallel whistler and ion-cyclotron branches both
    # degenerate to the Alfvén wave ω = ±k v_A, with the Walén phase relation
    # |δu_⊥/δB_⊥| = 1/√(ρ0) = 1 (normalized). Verified on the independent oracle.
    vA = 1.0
    for K in (0.05, 0.01)
        f = dispersion_frequencies((K, 0.0, 0.0), (1.0, 0.0, 0.0); cs = 0.0)
        @test length(f) == 2
        for ω in f
            @test isapprox(ω / K, vA; rtol = 3 * K)     # both branches → k v_A, error O(K)
        end
        # branch splitting shrinks with K (convergence to the degenerate Alfvén mode)
        @test (f[2] - f[1]) / K < 3 * K
        # Walén relation for each branch
        vals, vecs = dispersion_eigen((K, 0.0, 0.0), (1.0, 0.0, 0.0); cs = 0.0)
        for ω in f
            i = argmin(abs.(vals .- ω))
            @test isapprox(abs(vecs[2, i] / vecs[5, i]), 1.0; rtol = 3 * K)
        end
    end
end

@testset "HYB-006 polarization + group velocity" begin
    # Parallel whistler eigenvector: transverse δB is circularly polarized, δB_z = ±i δB_y.
    K = 1.0
    vals, vecs = dispersion_eigen((K, 0.0, 0.0), (1.0, 0.0, 0.0); cs = 0.0)
    ωW = (sqrt(K^4 + 4K^2) + K^2) / 2
    iw = argmin(abs.(vals .- ωW))
    by, bz = vecs[5, iw], vecs[6, iw]                 # δB_y, δB_z components
    @test abs(by) > 1e-6
    @test isapprox(abs(bz / by), 1.0; rtol = 1e-4)    # equal transverse magnitude
    @test isapprox(abs(imag(bz / by)), 1.0; atol = 1e-4)  # ±90° phase ⇒ circular
    @test abs(real(bz / by)) < 1e-4
    # Group velocity: the perpendicular fast magnetosonic mode is non-dispersive in
    # the MHD limit, so dω/dk = phase speed = √(vA²+cs²), exactly.
    for cs in (0.0, 0.5, 1.0)
        vg = group_velocity(1, (1.0, 0.0, 0.0), 0.02, (0.0, 0.0, 1.0); cs = cs)
        @test isapprox(vg, sqrt(1 + cs^2); rtol = 1e-3)
    end
    # And a DISPERSIVE branch is handled correctly: the parallel ion-cyclotron group
    # velocity at K=0.02 matches its analytic derivative d/dK[(√(K⁴+4K²)−K²)/2] ≈ 0.980.
    ωIC(k) = (sqrt(k^4 + 4k^2) - k^2) / 2
    vg_ic = group_velocity(1, (1.0, 0.0, 0.0), 0.02, (1.0, 0.0, 0.0); cs = 0.0)
    @test isapprox(vg_ic, (ωIC(0.02 + 1e-6) - ωIC(0.02 - 1e-6)) / 2e-6; rtol = 1e-3)
end

@testset "HYB-005 fast/slow closed form matches the eigensolver" begin
    # The fast_slow_speeds closed form equals the eigensolver's fast & slow roots
    # at oblique angles in the MHD (small-K) limit, for several β (via cs).
    for cs in (0.3, 0.6, 1.0), θ in (π / 8, π / 4, 3π / 8)
        cf, csl = fast_slow_speeds(1.0, cs, θ)
        f = dispersion_frequencies((0.01, 0.0, 0.0), (cos(θ), 0.0, sin(θ)); cs = cs)
        sp = sort(f ./ 0.01)
        @test isapprox(sp[end], cf; rtol = 3e-3)      # fast
        @test isapprox(sp[1], csl; rtol = 3e-3)       # slow
    end
end

@testset "HYB-005 PIC fast magnetosonic matches the fluid limit" begin
    # Perpendicular fast magnetosonic in the actual hybrid PIC: B0 = ẑ, k = x̂,
    # warm isothermal electrons (cs² = Te, vA = 1) ⇒ ω = k√(1+Te). This is the
    # cleanly-excitable fluid-regime magnetosonic mode (no kinetic damping at
    # perpendicular propagation), per the checklist's "compare in a regime where
    # the fluid approximation is valid".
    T = Float64
    for (Te, tol) in ((0.25, 0.03), (0.5, 0.03))
        n, L, m, nppc, dt, nsteps = 64, 2π, 1, 1000, 0.02, 1500
        g = FourierGrid((n,), (T(L),))
        N = nppc * n
        ps = ParticleSet{1,T}(N)
        load_lattice_1d!(ps, 0.0, T(L))
        set_density_weight!(ps, 1.0, g)
        load_quiet_velocities!(ps, MersenneTwister(7), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        st = HybridStepper(g, HybridModel(IsothermalElectrons(T(Te))), CIC(), N)
        fill!(st.fields.B[3], 1.0)                       # B0 = ẑ
        k = 2π * m / L
        x = [(i - 1) * g.dx[1] for i = 1:n]
        st.fields.B[3] .+= 0.01 .* cos.(k .* x)          # compressional seed
        init!(st, ps)
        s = Float64[]
        for _ = 1:nsteps
            step!(st, ps, dt; NB = 2)
            push!(s, real(mode_amplitude(st.fields.B[3], g, (m,))))
        end
        ω = _freq_zerocross(s, dt)
        ωth = k * sqrt(1 + Te)                           # = oracle fast branch, perpendicular
        @test isapprox(ω, ωth; rtol = tol)
        # cross-check: the oracle's perpendicular fast branch equals ωth
        fo = dispersion_frequencies((k, 0.0, 0.0), (0.0, 0.0, 1.0); cs = sqrt(Te))
        @test any(z -> isapprox(z, ωth; rtol = 1e-6), fo)
    end
end
