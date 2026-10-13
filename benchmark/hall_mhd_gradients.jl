# Reproduce the Hall velocity-gradient tradeoff with the current full RHS.
#   JULIA_NUM_THREADS=1 julia --project=. benchmark/hall_mhd_gradients.jl
# Ratios > 1 favor the shared-gradient implementation. Process CPU time isolates
# this process's work; wall time remains sensitive to other running workloads.
# There are no portable latency thresholds. Compare the printed allocation cost
# and both timing distributions on the deployment machine.
using HybridPlasmaPIC, FFTW, Statistics

Threads.nthreads() == 1 || error("Run this benchmark with JULIA_NUM_THREADS=1")
Sys.isunix() || error("This CPU timing harness requires a POSIX clock() implementation")
FFTW.set_num_threads(1)
println((
    julia = VERSION,
    julia_threads = Threads.nthreads(),
    fftw_threads = FFTW.get_num_threads(),
))

# Derive the previous per-axis implementation from the current RHS, replacing
# only velocity-gradient evaluation. Assert the source shape so a future RHS
# refactor cannot silently turn this into a comparison of identical methods.
source = read(joinpath(@__DIR__, "..", "src", "models", "hall_mhd.jl"), String)
first_marker = findfirst("function _hall_mhd_rhs!(", source)
last_marker = findfirst("\n\"\"\"\n    hall_mhd_rhs!", source)
@assert first_marker !== nothing && last_marker !== nothing
baseline_source = source[first(first_marker):first(last_marker)-1]
shared = join(
    (
        "        if D == 1",
        "            deriv!(gradu[1], u[c], g, 1)",
        "        else",
        "            gradient!(gradu, u[c], g)",
        "        end",
    ),
    "\n",
)
per_axis =
    join(("        for d = 1:D", "            deriv!(gradu[d], u[c], g, d)", "        end"), "\n")
@assert occursin(shared, baseline_source)
baseline_source = replace(
    baseline_source,
    "function _hall_mhd_rhs!(" => "function _hall_mhd_rhs_per_axis_benchmark!(",
    shared => per_axis,
)
include_string(HybridPlasmaPIC, baseline_source, @__FILE__)

function benchmark_state(D, T, n)
    g = FourierGrid(ntuple(_ -> n, D), ntuple(_ -> T(2π), D))
    st = HallMHDState(g, HallMHDModel(IsothermalElectrons(0.2); Ti = 0.1, η = 0.01))
    for I in CartesianIndices(st.fields.n)
        phase = sum((I[d] - 1) * g.dx[d] for d = 1:D)
        st.fields.n[I] = 1 + 0.1cos(phase)
        for c = 1:3
            st.fields.ui[c][I] = 0.05sin(c * phase)
            st.fields.B[c][I] = 0.1cos(phase) + 0.2c
        end
    end
    return st
end

function call_rhs!(fn, st)
    return fn(st.rhs_n, st.rhs_u, st.rhs_B, st, st.fields.n, st.fields.ui, st.fields.B)
end

function measured_batch(fn, st, repetitions)
    cpu_start = ccall(:clock, Clong, ())
    cpu_start >= 0 || error("Process CPU clock is unavailable")
    wall_start = time_ns()
    for _ = 1:repetitions
        call_rhs!(fn, st)
    end
    wall = (time_ns() - wall_start) / 1e9
    ticks = ccall(:clock, Clong, ()) - cpu_start
    ticks > 0 || error("Process CPU clock did not advance")
    return (; wall, ticks)
end

function rhs_allocations(fn, st)
    call_rhs!(fn, st)
    return @allocated call_rhs!(fn, st)
end

baseline = HybridPlasmaPIC._hall_mhd_rhs_per_axis_benchmark!
current = HybridPlasmaPIC._hall_mhd_rhs!
for T in (Float32, Float64), D = 1:3
    st = benchmark_state(D, T, 8)
    call_rhs!(baseline, st)
    expected = deepcopy((st.rhs_n, st.rhs_u, st.rhs_B))
    call_rhs!(current, st)
    @assert expected == (st.rhs_n, st.rhs_u, st.rhs_B)
    println((
        T,
        D,
        baseline_bytes = rhs_allocations(baseline, st),
        current_bytes = rhs_allocations(current, st),
    ))
end

for (D, n) in ((2, 128), (3, 64))
    st = benchmark_state(D, Float64, n)
    call_rhs!(baseline, st)
    expected = deepcopy((st.rhs_n, st.rhs_u, st.rhs_B))
    call_rhs!(current, st)
    @assert expected == (st.rhs_n, st.rhs_u, st.rhs_B)
    for _ = 1:3
        measured_batch(baseline, st, 1)
        measured_batch(current, st, 1)
    end
    baseline_wall, current_wall, wall_ratios, cpu_ratios = (Float64[] for _ = 1:4)
    for pair = 1:12
        if isodd(pair)
            old = measured_batch(baseline, st, 3)
            new = measured_batch(current, st, 3)
        else
            new = measured_batch(current, st, 3)
            old = measured_batch(baseline, st, 3)
        end
        push!(baseline_wall, old.wall / 3)
        push!(current_wall, new.wall / 3)
        push!(wall_ratios, old.wall / new.wall)
        push!(cpu_ratios, old.ticks / new.ticks)
        println((D, n, pair, wall_ratio = last(wall_ratios), cpu_ratio = last(cpu_ratios)))
        flush(stdout)
    end
    println((
        D,
        n,
        baseline_wall = median(baseline_wall),
        current_wall = median(current_wall),
        wall_ratio = median(wall_ratios),
        wall_ratio_range = extrema(wall_ratios),
        cpu_ratio = median(cpu_ratios),
        cpu_ratio_range = extrema(cpu_ratios),
    ))
    flush(stdout)
end
