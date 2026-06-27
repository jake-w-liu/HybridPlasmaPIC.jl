# benchmarks.jl — performance regression harness for the hot kernels. Reports
# throughput in particles-advanced/s and cell-updates/s (the §21.4 / Phase-8
# performance metrics, measured on CPU) and a steady-state allocation check.
#
#   julia --project=. benchmark/benchmarks.jl
#
# Prints a table; exits 0. A CI job can diff the numbers against a stored
# baseline to catch regressions (see .github/workflows).

using HybridPlasmaPIC, Printf, Random

best(f, reps) = minimum(begin
    f()                       # warm up once before timing each rep
    (@elapsed f())
end for _ = 1:reps)

println("HybridPlasmaPIC kernel benchmarks (CPU, minimum of repeated runs)")
println("="^64)

# spectral derivative
let
    g = FourierGrid((256,), (2π,))
    f = randn(256)
    out = similar(f)
    deriv!(out, f, g, 1)
    t = best(() -> deriv!(out, f, g, 1), 50)
    a = @allocated deriv!(out, f, g, 1)
    @printf("deriv! 1D n=256        : %8.2f µs   alloc=%d B\n", 1e6t, a)
end

# CIC deposition throughput
let
    n = 64
    g = FourierGrid((n,), (2π,))
    N = 200 * n
    ps = ParticleSet{1,Float64}(N)
    load_uniform!(ps, MersenneTwister(1), (0.0,), (2π,))
    set_density_weight!(ps, 1.0, g)
    out = zeros(n)
    deposit_scalar!(out, ps, ps.weight, g, CIC())
    t = best(() -> deposit_scalar!(out, ps, ps.weight, g, CIC()), 20)
    @printf("CIC deposit  N=%d  : %8.2f µs   %.1f Mparticles/s\n", N, 1e6t, N / t / 1e6)
end

# full hybrid step throughput (particles advanced per second)
let
    n = 32
    L = 2π
    g = FourierGrid((n,), (L,))
    N = 100 * n
    ps = ParticleSet{1,Float64}(N)
    load_lattice_1d!(ps, 0.0, L)
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(2), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    fill!(st.fields.B[3], 1.0)
    init!(st, ps)
    step!(st, ps, 0.02; NB = 2)
    t = best(() -> step!(st, ps, 0.02; NB = 2), 10)
    @printf("hybrid step  N=%d   : %8.2f ms   %.2f Mparticle-steps/s\n", N, 1e3t, N / t / 1e6)
end
println("="^64)
println("done")
