# Linear-wave dispersion: HYB-001 stationarity, HYB-002 ion-acoustic, and
# HYB-003/004 parallel Alfvén/whistler/ion-cyclotron branches, each measured
# from an actual PIC run and compared to the analytic dispersion relation.
# Stochastic test: fixed seeds + tolerances set well above the measured error
# (FFT bin-resolution limited to ~1-2%; tightens with longer integration).

using HybridPlasmaPIC, FFTW, Test, Random, Statistics, LinearAlgebra

# 1-D cold/warm proton setup with n0=1 and B0
function setup1d(n, L, seed; nppc = 400, vth = (0.0, 0.0, 0.0), Te = 0.0, B0 = (0.0, 0.0, 0.0))
    T = Float64
    g = FourierGrid((n,), (T(L),))
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, T(L))
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), vth)
    st = HybridStepper(g, HybridModel(IsothermalElectrons(T(Te))), CIC(), N)
    for c = 1:3
        fill!(st.fields.B[c], T(B0[c]))
    end
    return g, ps, st
end

# dominant frequency of a real series via early-window neg→pos zero crossings
function freq_zerocross(series, dt)
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

# split a complex mode-amplitude series into its dominant +freq and −freq peaks
# (the two parallel branches counter-rotate), with parabolic interpolation
function branch_freqs(s, dt)
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

@testset "HYB-001 uniform equilibrium stays stationary" begin
    g, ps, st =
        setup1d(32, 2π, 5; nppc = 400, vth = (0.1, 0.1, 0.1), Te = 0.5, B0 = (0.0, 0.0, 1.0))
    init!(st, ps)
    E0 = magnetic_energy(st.fields.B, g)
    for _ = 1:300
        step!(st, ps, 0.05; NB = 4)
    end
    @test all(isfinite, st.fields.B[3])
    @test abs(mean(st.fields.B[3]) - 1.0) < 1e-10          # mean field exactly conserved (curl has no k=0)
    @test abs(magnetic_energy(st.fields.B, g) - E0) / E0 < 0.1
    @test maximum(abs, st.fields.E[1]) < 0.3              # no coherent E, only noise
end

@testset "HYB-002 ion-acoustic ω = k·c_s" begin
    Te = 1.0
    m = 1
    L = 2π
    k = 2π * m / L
    g, ps, st = setup1d(64, L, 1; nppc = 600, Te = Te)  # B0=0 ⇒ electrostatic
    for p = 1:nparticles(ps)
        ps.v[1][p] += 0.005 * sin(k * ps.x[1][p])
    end
    init!(st, ps)
    dt = 0.02
    series = Float64[]
    for _ = 1:700                                        # early linear window (cold IA steepens later)
        step!(st, ps, dt)
        push!(series, real(mode_amplitude(st.fields.n, g, (m,))))
    end
    ω = freq_zerocross(series, dt)
    @test abs(ω - k * sqrt(Te)) / (k * sqrt(Te)) < 0.03
end

@testset "HYB-003/004 parallel whistler + ion-cyclotron" begin
    for (m, dt, ns) in ((1, 0.02, 2500), (2, 0.015, 2500))
        L = 2π
        k = 2π * m / L
        K = k
        ωW = (sqrt(K^4 + 4K^2) + K^2) / 2
        ωIC = (sqrt(K^4 + 4K^2) - K^2) / 2
        g, ps, st = setup1d(32, L, 2; nppc = 400, B0 = (1.0, 0.0, 0.0))  # B0 ∥ propagation
        x = [(i - 1) * g.dx[1] for i = 1:g.n[1]]
        st.fields.B[2] .= 0.005 .* cos.(k .* x)          # transverse seed: equal R + L
        init!(st, ps)
        s = ComplexF64[]
        for _ = 1:ns
            step!(st, ps, dt; NB = 8)                    # subcycle: whistler CFL
            push!(
                s,
                mode_amplitude(st.fields.B[2], g, (m,)) +
                im * mode_amplitude(st.fields.B[3], g, (m,)),
            )
        end
        @test all(isfinite, real.(s))
        hi, lo = branch_freqs(s, dt)
        @test abs(hi - ωW) / ωW < 0.04
        @test abs(lo - ωIC) / ωIC < 0.05
    end
end

@testset "HYB-007 adiabatic energy convergence" begin
    clo = PolytropicElectrons(0.5, 1.0, 5 / 3)        # adiabatic, η=0
    function max_drift(dt; nsteps = 1200, NB = 4)
        T = Float64
        n = 32
        L = 2π
        k = 2π / L
        g = FourierGrid((n,), (L,))
        N = 300 * n
        ps = ParticleSet{1,T}(N)
        load_lattice_1d!(ps, 0.0, L)
        set_density_weight!(ps, 1.0, g)
        load_quiet_velocities!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
        st = HybridStepper(g, HybridModel(clo), CIC(), N)
        fill!(st.fields.B[1], 1.0)
        x = [(i - 1) * g.dx[1] for i = 1:n]
        st.fields.B[2] .= 0.05 .* cos.(k .* x)
        init!(st, ps)
        Etot() =
            kinetic_energy(ps) +
            magnetic_energy(st.fields.B, g) +
            electron_internal_energy(st.fields.n, clo, g)
        E0 = Etot()
        d = 0.0
        for _ = 1:nsteps
            step!(st, ps, dt; NB)
            d = max(d, abs(Etot() - E0) / E0)
        end
        return d
    end
    d_coarse = max_drift(0.04)
    d_fine = max_drift(0.01)
    @test d_coarse < 0.15                  # bounded (no secular blow-up)
    @test d_fine < d_coarse                # energy drift decreases with Δt
end

@testset "HYB-008 subcycling consistency" begin
    function b_after(NB; nsteps = 400)
        T = Float64
        n = 32
        L = 2π
        k = 2π / L
        g = FourierGrid((n,), (L,))
        N = 256 * n
        ps = ParticleSet{1,T}(N)
        load_lattice_1d!(ps, 0.0, L)
        set_density_weight!(ps, 1.0, g)
        load_quiet_velocities!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        st = HybridStepper(g, HybridModel(IsothermalElectrons(0.0)), CIC(), N)
        fill!(st.fields.B[1], 1.0)
        x = [(i - 1) * g.dx[1] for i = 1:n]
        st.fields.B[2] .= 0.01 .* cos.(k .* x)
        init!(st, ps)
        for _ = 1:nsteps
            step!(st, ps, 0.02; NB)
        end
        return copy(st.fields.B[2])
    end
    b4 = b_after(4)
    b8 = b_after(8)
    @test norm(b4 .- b8) / norm(b8) < 0.05   # field is converged in the subcycle count
end
