# Steady-state cost of scalar, one-population, and two-ion hybrid evolution.
#   JULIA_NUM_THREADS=1 julia --project=. benchmark/multispecies_hybrid.jl
# Reports measured costs, not portable pass/fail latency thresholds. The
# two-ion case has 50% more particles and a second set of moment deposits.
using HybridPlasmaPIC, FFTW, Statistics, Printf

Threads.nthreads() == 1 || error("Run with JULIA_NUM_THREADS=1")
FFTW.set_num_threads(1)

benchmark_init!(st::HybridStepper, ps) = init!(st, ps)
benchmark_init!(st::CAMCLStepper, ps) = init_camcl!(st, ps)
benchmark_step!(st::HybridStepper, ps) = step!(st, ps, 1e-5; NB = 2)
benchmark_step!(st::CAMCLStepper, ps) = step_camcl!(st, ps, 1e-5; NB = 2)

function benchmark_population(g::FourierGrid{D,T}, per_cell, q, m, density, velocity) where {D,T}
    ps = ParticleSet{D,T}(prod(g.n) * per_cell; q = T(q), m = T(m))
    p = 0
    for I in CartesianIndices(g.n), j = 1:per_cell
        p += 1
        for d = 1:D
            fraction = d == 1 ? (T(j) - T(0.5)) / per_cell : T(0.5)
            ps.x[d][p] = (I[d] - 1 + fraction) * g.dx[d]
        end
        for c = 1:3
            ps.v[c][p] = T(velocity[c])
        end
    end
    set_density_weight!(ps, density, g)
    return ps
end

function measure_multispecies(St, D, n, mode)
    g = FourierGrid(ntuple(_ -> n, D), ntuple(_ -> 2π, D))
    proton = benchmark_population(g, 64, 1, 1, 0.8, (0.01, 0.12, -0.04))
    populations = if mode === :scalar
        proton
    elseif mode === :collection
        [proton]
    else
        [proton, benchmark_population(g, 32, 2, 4, 0.1, (-0.01, -0.08, 0.09))]
    end
    st = St(g, HybridModel(IsothermalElectrons(0.1)), CIC(), 0)
    fill!(st.fields.B[1], 1.0)
    benchmark_init!(st, populations)
    for _ = 1:3
        benchmark_step!(st, populations)
    end
    bytes = @allocated benchmark_step!(st, populations)
    samples = [@elapsed benchmark_step!(st, populations) for _ = 1:11]
    all(all(isfinite, component) for component in st.fields.B) || error("Nonfinite benchmark field")
    count = mode === :scalar ? nparticles(populations) : sum(nparticles, populations)
    return (; count, bytes, median_seconds = median(samples), minimum_seconds = minimum(samples))
end

println((
    julia = VERSION,
    julia_threads = Threads.nthreads(),
    fftw_threads = FFTW.get_num_threads(),
))
for St in (HybridStepper, CAMCLStepper), (D, n) in ((1, 32), (2, 16), (3, 8))
    for mode in (:scalar, :collection, :two_ions)
        result = measure_multispecies(St, D, n, mode)
        @printf(
            "%s %dD %-10s N=%6d median=%8.3f ms min=%8.3f ms alloc=%d B\n",
            nameof(St),
            D,
            mode,
            result.count,
            1e3result.median_seconds,
            1e3result.minimum_seconds,
            result.bytes
        )
    end
end
