# test_camcl.jl — CAM-CL hybrid integrator (Matthews 1994) verification.
#
# Oracle = the SAME analytic dispersion the existing integrator passes:
#   1-D ion-acoustic wave (B0=0, isothermal Te, cold quiet-start protons,
#   seed δv ∝ sin(kx)) ⇒ ω = k·√Te. Require |ω − k√Te| < 4%.
# Plus an integrator-comparison: run the SAME wave with the EXISTING step! and
# with CAM-CL; require the two measured frequencies agree within 5%.
# Plus an optional parallel whistler/ion-cyclotron check.

using HybridPlasmaPIC, FFTW, Test, Random, Statistics, LinearAlgebra

# dominant frequency of a real series via early-window neg→pos zero crossings
function camcl_freq_zerocross(series, dt)
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

# split a complex mode-amplitude series into its dominant +/− freq peaks
function camcl_branch_freqs(s, dt)
    Nt = length(s)
    S = fft(s)
    mag = abs.(S)
    refine(j) = begin
        a = mag[mod(j - 1, Nt)+1]
        b = mag[j+1]
        c = mag[mod(j + 1, Nt)+1]
        d = a - 2b + c
        δ = d != 0 ? 0.5 * (a - c) / d : 0.0
        jr = j + δ
        jj = jr <= Nt / 2 ? jr : jr - Nt
        2π * jj / (Nt * dt)
    end
    bp = 0
    mp = -1.0
    bn = 0
    mn = -1.0
    for j = 1:Nt-1
        jj = j <= Nt ÷ 2 ? j : j - Nt
        ω = 2π * jj / (Nt * dt)
        if ω > 0.05 && mag[j+1] > mp
            mp = mag[j+1]
            bp = j
        end
        if ω < -0.05 && mag[j+1] > mn
            mn = mag[j+1]
            bn = j
        end
    end
    f1 = abs(refine(bp))
    f2 = abs(refine(bn))
    return max(f1, f2), min(f1, f2)
end

# Build a fresh ion-acoustic setup (B0=0, isothermal Te, cold quiet protons,
# seed δv ∝ sin(kx)) and return the seeded ParticleSet + grid.
function ia_setup(; n = 64, L = 2π, seed = 1, nppc = 600, Te = 1.0, mode = 1, amp = 0.005)
    T = Float64
    g = FourierGrid((n,), (T(L),))
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    k = 2π * mode / L
    for p = 1:nparticles(ps)
        ps.v[1][p] += amp * sin(k * ps.x[1][p])
    end
    return g, ps
end

@testset "CAM-CL validates timestep and subcycles before mutation" begin
    T = Float64
    g = FourierGrid((4,), (2π,))
    ps = ParticleSet{1,T}(4)
    load_lattice_1d!(ps, 0.0, 2π)
    set_density_weight!(ps, 1.0, g)
    st = CAMCLStepper(g, HybridModel(IsothermalElectrons(0.0)), NGP(), nparticles(ps))
    fill!(st.fields.B[3], 1.0)
    init_camcl!(st, ps)

    x0 = ntuple(d -> copy(ps.x[d]), 1)
    v0 = ntuple(c -> copy(ps.v[c]), 3)
    B0 = ntuple(c -> copy(st.fields.B[c]), 3)
    time0 = st.time[]
    step0 = st.step[]

    @test_throws ArgumentError step_camcl!(st, ps, 0.1; NB = 1)
    @test_throws ArgumentError step_camcl!(st, ps, 0.1; NB = 0)
    @test_throws ArgumentError step_camcl!(st, ps, NaN; NB = 2)
    @test_throws ArgumentError step_camcl!(st, ps, -0.1; NB = 2)
    @test st.time[] == time0
    @test st.step[] == step0
    @test all(ps.x[d] == x0[d] for d = 1:1)
    @test all(ps.v[c] == v0[c] for c = 1:3)
    @test all(st.fields.B[c] == B0[c] for c = 1:3)
end

@testset "CAM-CL ion-acoustic ω = k·√Te (analytic oracle)" begin
    Te = 1.0
    m = 1
    L = 2π
    k = 2π * m / L
    g, ps = ia_setup(; n = 64, L = L, seed = 1, nppc = 600, Te = Te, mode = m, amp = 0.005)
    st = CAMCLStepper(g, HybridModel(IsothermalElectrons(Te)), CIC(), nparticles(ps))
    # B0 = 0 already (zeros); electrostatic limit
    init_camcl!(st, ps)
    dt = 0.02
    series = Float64[]
    for _ = 1:700
        step_camcl!(st, ps, dt)
        push!(series, real(mode_amplitude(st.fields.n, g, (m,))))
    end
    @test all(isfinite, series)
    ω = camcl_freq_zerocross(series, dt)
    rel = abs(ω - k * sqrt(Te)) / (k * sqrt(Te))
    @info "CAM-CL ion-acoustic" ω_measured = ω ω_analytic = k * sqrt(Te) rel_err = rel
    @test rel < 0.04
end

@testset "CAM-CL vs existing step! frequency agreement" begin
    Te = 1.0
    m = 1
    L = 2π
    k = 2π * m / L
    dt = 0.02
    nsteps = 700

    # existing integrator
    g1, ps1 = ia_setup(; n = 64, L = L, seed = 1, nppc = 600, Te = Te, mode = m, amp = 0.005)
    st1 = HybridStepper(g1, HybridModel(IsothermalElectrons(Te)), CIC(), nparticles(ps1))
    init!(st1, ps1)
    s1 = Float64[]
    for _ = 1:nsteps
        step!(st1, ps1, dt)
        push!(s1, real(mode_amplitude(st1.fields.n, g1, (m,))))
    end
    ω_old = camcl_freq_zerocross(s1, dt)

    # CAM-CL
    g2, ps2 = ia_setup(; n = 64, L = L, seed = 1, nppc = 600, Te = Te, mode = m, amp = 0.005)
    st2 = CAMCLStepper(g2, HybridModel(IsothermalElectrons(Te)), CIC(), nparticles(ps2))
    init_camcl!(st2, ps2)
    s2 = Float64[]
    for _ = 1:nsteps
        step_camcl!(st2, ps2, dt)
        push!(s2, real(mode_amplitude(st2.fields.n, g2, (m,))))
    end
    ω_new = camcl_freq_zerocross(s2, dt)

    rel = abs(ω_new - ω_old) / ω_old
    @info "integrator comparison" ω_existing = ω_old ω_camcl = ω_new rel_diff = rel
    @test rel < 0.05
end

@testset "CAM-CL parallel whistler + ion-cyclotron (optional)" begin
    m = 1
    L = 2π
    dt = 0.02
    ns = 2500
    k = 2π * m / L
    K = k
    ωW = (sqrt(K^4 + 4K^2) + K^2) / 2
    ωIC = (sqrt(K^4 + 4K^2) - K^2) / 2
    T = Float64
    n = 32
    g = FourierGrid((n,), (T(L),))
    N = 400 * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    st = CAMCLStepper(g, HybridModel(IsothermalElectrons(0.0)), CIC(), N)
    fill!(st.fields.B[1], 1.0)                  # B0 ∥ propagation
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st.fields.B[2] .= 0.005 .* cos.(k .* x)     # transverse seed: equal R + L
    init_camcl!(st, ps)
    s = ComplexF64[]
    for _ = 1:ns
        step_camcl!(st, ps, dt; NB = 24)        # cyclic-leapfrog whistler CFL (tighter than RK4)
        push!(
            s,
            mode_amplitude(st.fields.B[2], g, (m,)) + im * mode_amplitude(st.fields.B[3], g, (m,)),
        )
    end
    @test all(isfinite, real.(s))
    hi, lo = camcl_branch_freqs(s, dt)
    @info "CAM-CL parallel branches" whistler = hi ωW ion_cyclotron = lo ωIC
    @test abs(hi - ωW) / ωW < 0.06
    @test abs(lo - ωIC) / ωIC < 0.07
end
