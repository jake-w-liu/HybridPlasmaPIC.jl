# Ion-acoustic wave: a minimal end-to-end HybridPlasmaPIC run that measures the wave
# frequency from a PIC simulation and compares it to ω = k·c_s.
#
#   julia --project=. examples/ion_acoustic.jl

using HybridPlasmaPIC, Random, Printf

const T = Float64
n = 64;
L = 2π;
g = FourierGrid((n,), (L,))
Te = 1.0;
m = 1;
k = 2π * m / L;                 # mode 1
cs = sqrt(Te)                                   # isothermal sound speed (γ_e = 1)

# cold, quiet-start protons with n0 = 1
N = n * 400
ps = ParticleSet{1,T}(N)
load_lattice_1d!(ps, 0.0, L)
set_density_weight!(ps, 1.0, g)
load_quiet_velocities!(ps, MersenneTwister(1), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))

# seed a small longitudinal velocity perturbation
for p = 1:N
    ps.v[1][p] += 0.005 * sin(k * ps.x[1][p])
end

# B0 = 0 ⇒ electrostatic; isothermal electrons
st = HybridStepper(g, HybridModel(IsothermalElectrons(Te)), CIC(), N)
init!(st, ps)

dt = 0.02
series = Float64[]
for _ = 1:700
    step!(st, ps, dt)
    push!(series, real(mode_amplitude(st.fields.n, g, (m,))))
end

# frequency from the first few zero crossings (early linear window)
zc = Int[]
for i = 2:length(series)
    series[i-1] < 0 && series[i] >= 0 && push!(zc, i)
end
period = (zc[end] - zc[1]) / (length(zc) - 1) * dt
ω = 2π / period

@printf("ion-acoustic wave  k = %.3f\n", k)
@printf("  ω measured = %.4f\n", ω)
@printf("  ω = k·c_s   = %.4f   (error %.2f%%)\n", k * cs, 100 * abs(ω - k * cs) / (k * cs))
