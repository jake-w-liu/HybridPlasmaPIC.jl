# test_completeness.jl — closure / consistency diagnostics:
#   * particle_work! matches q (v·E) dt·Σw to roundoff on a uniform field
#   * mixed_divcurl_residual is small and DECREASES with x-refinement
#   * ion-acoustic frequency is independent of the startup transient amplitude

using HybridPlasmaPIC, Test, Random, Statistics

using FFTW

# dominant positive frequency of a real series via its DFT peak with parabolic
# interpolation (robust to the mild nonlinear steepening that breaks a
# zero-crossing estimate at larger seed amplitude).
function _freq_fft(series, dt)
    Nt = length(series)
    S = abs.(fft(series .- sum(series) / Nt))   # remove DC
    best = 0
    mb = -1.0
    for j = 1:Nt-1
        jj = j <= Nt ÷ 2 ? j : j - Nt
        ω = 2π * jj / (Nt * dt)
        if ω > 0.05 && S[j+1] > mb
            mb = S[j+1]
            best = j
        end
    end
    best == 0 && return NaN
    a = S[mod(best - 1, Nt)+1]
    b = S[best+1]
    c = S[mod(best + 1, Nt)+1]
    d = a - 2b + c
    δ = d != 0 ? 0.5 * (a - c) / d : 0.0       # parabolic peak refinement
    jr = best + δ
    jj = jr <= Nt / 2 ? jr : jr - Nt
    return 2π * jj / (Nt * dt)
end

@testset "particle_work! matches q v·E dt on a uniform field" begin
    T = Float64
    n = 16
    L = 2π
    g = FourierGrid((n,), (L,))
    N = 500
    ps = ParticleSet{1,T}(N; q = 1.7, m = 1.0)
    rng = MersenneTwister(42)
    load_uniform!(ps, rng, (0.0,), (L,))
    # nonuniform weights and velocities — the closure must not assume w=1 or v=const
    for p = 1:N
        ps.weight[p] = 0.5 + rand(rng)
        ps.v[1][p] = randn(rng)
        ps.v[2][p] = randn(rng)
        ps.v[3][p] = randn(rng)
    end
    # spatially UNIFORM E so gather(E)(x_p) = E0 for every particle (any shape, any x)
    E0 = (0.3, -0.45, 0.8)
    Egrid = ntuple(c -> fill(T(E0[c]), n), 3)
    dt = 0.05
    for shape in (NGP(), CIC(), TSC())
        work = zeros(T, N)
        particle_work!(work, ps, Egrid, g, shape, dt)
        # exact per-particle reference q (v·E) dt
        ref =
            [ps.q * (ps.v[1][p] * E0[1] + ps.v[2][p] * E0[2] + ps.v[3][p] * E0[3]) * dt for p = 1:N]
        @test maximum(abs.(work .- ref)) < 1e-12
        # accumulation: a second call doubles the work (adds, not overwrites)
        particle_work!(work, ps, Egrid, g, shape, dt)
        @test maximum(abs.(work .- 2 .* ref)) < 1e-12
    end
    poisoned = fill(T(7), N)
    snapshot = copy(poisoned)
    @test_throws ArgumentError particle_work!(poisoned, ps, Egrid, g, CIC(), NaN)
    @test poisoned == snapshot
    @test_throws ArgumentError particle_work!(poisoned, ps, Egrid, g, CIC(), Inf)
    @test poisoned == snapshot
    # population total weighted by particle weight: Σ_p w_p q (v_p·E) dt
    shape = CIC()
    wwork = zeros(T, N)
    # weight the per-particle work by w_p via scaling the velocity-independent
    # increment: build a weighted reference and compare against a weighted sum
    particle_work!(wwork, ps, Egrid, g, shape, dt)
    tot = sum(ps.weight[p] * wwork[p] for p = 1:N)
    ref_tot = sum(
        ps.weight[p] * ps.q * (ps.v[1][p] * E0[1] + ps.v[2][p] * E0[2] + ps.v[3][p] * E0[3]) * dt for p = 1:N
    )
    @test abs(tot - ref_tot) < 1e-10 * abs(ref_tot)
end

@testset "mixed_divcurl_residual small, not machine-zero, converges in nx" begin
    T = Float64
    Ly = 2π
    ny = 16
    ns = [33, 65, 129, 257]
    rs = Float64[]
    for nx in ns
        s = SBP1D(nx, 1.0)            # Lx = 1.0
        r = mixed_divcurl_residual(s, nx, ny, Ly)
        push!(rs, r)
        @test isfinite(r)
        @test r < 1e-1               # bounded by SBP truncation
        @test r > 1e-8               # NOT machine zero — genuine SBP-x truncation
    end
    # strictly decreasing under x-refinement
    for i = 2:length(rs)
        @test rs[i] < rs[i-1]
    end
    rate = log(rs[1] / rs[end]) / log((ns[end] - 1) / (ns[1] - 1))
    @info "mixed_divcurl_residual" residuals = rs convergence_rate = rate
    @test 1.3 < rate < 1.7           # diagonal-norm SBP-(2,1) H-norm rate ≈ 3/2
end

@testset "ion-acoustic frequency independent of startup amplitude" begin
    Te = 1.0
    m = 1
    L = 2π
    k = 2π * m / L
    dt = 0.02

    function ia_freq(amp; seed = 1)
        T = Float64
        n = 64
        nppc = 600
        N = nppc * n
        g = FourierGrid((n,), (T(L),))
        ps = ParticleSet{1,T}(N)
        load_lattice_1d!(ps, 0.0, T(L))
        set_density_weight!(ps, 1.0, g)
        load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        st = HybridStepper(g, HybridModel(IsothermalElectrons(T(Te))), CIC(), N)
        # B0 = 0 ⇒ electrostatic ion-acoustic branch
        for p = 1:nparticles(ps)
            ps.v[1][p] += amp * sin(k * ps.x[1][p])     # startup transient: scaled seed
        end
        init!(st, ps)
        series = Float64[]
        for _ = 1:700
            step!(st, ps, dt)
            push!(series, real(mode_amplitude(st.fields.n, g, (m,))))
        end
        return _freq_fft(series, dt)
    end

    ω_small = ia_freq(0.005; seed = 1)
    ω_large = ia_freq(0.02; seed = 7)   # 4× the seed amplitude, different seed
    ω_theory = k * sqrt(Te)
    @info "ion-acoustic startup independence" ω_small ω_large ω_theory
    @test isfinite(ω_small) && isfinite(ω_large)
    # crucially: frequency is INDEPENDENT of the startup transient (within 4%)
    @test abs(ω_small - ω_large) / ω_large < 0.04
    # and both track the analytic ω = k·c_s to within the FFT bin resolution
    # (bin width 2π/(700·0.02) ≈ 0.45 ⇒ ~one-bin, ~10%, offset at this window)
    @test abs(ω_small - ω_theory) / ω_theory < 0.10
    @test abs(ω_large - ω_theory) / ω_theory < 0.10
end
