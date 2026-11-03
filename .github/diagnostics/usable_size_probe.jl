# Linux probe for the intermittent test_linear_cn_range.jl allocation equality.
#
# Julia counts malloc-backed array buffers and counted mallocs (GMP limbs) by
# memory_block_usable_size, which is malloc_usable_size on Linux
# (src/gc-common.c, src/gc-stock.c in v1.12.7 and v1.13.3). This probe checks,
# for identical repeated calls:
#   1. whether @allocated varies,
#   2. whether bytes minus the usable size of the reachable malloc-backed
#      buffers is constant (the variation is exactly allocator slack),
#   3. whether the profiled requested allocations (type, size) and the
#      allocation-event count are identical across calls.
using Profile
using Random
using HybridPlasmaPIC

const SAMPLES = 64
const B0 = Ref{Int64}(0)
const B1 = Ref{Int64}(0)

@noinline function counted(f::F) where {F}
    Base.gc_bytes(B0)
    x = f()
    Base.gc_bytes(B1)
    return B1[] - B0[], x
end

# The same query Julia's memory_block_usable_size makes on each platform.
usable(a::Array) = Int(@static if Sys.isapple()
    ccall(:malloc_size, Csize_t, (Ptr{Cvoid},), pointer(a))
else
    ccall(:malloc_usable_size, Csize_t, (Ptr{Cvoid},), pointer(a))
end)

# Arrays above the 2032-byte pool limit own a separately malloc'd buffer.
function malloc_arrays(x, depth = 2, out = Array[])
    if x isa Array && isbitstype(eltype(x)) && sizeof(x) > 2032
        push!(out, x)
    elseif depth > 0 && !(x isa Module) && isstructtype(typeof(x))
        for i = 1:nfields(x)
            isdefined(x, i) && malloc_arrays(getfield(x, i), depth - 1, out)
        end
    end
    return out
end

function distribution(xs)
    u = sort!(unique(xs))
    return join(("$(v)×$(count(==(v), xs))" for v in u), " ")
end

function report(name, f; keep = Any[])
    bytes, slack, residual = Int[], Int[], Int[]
    f()
    for i = 1:SAMPLES
        b, x = counted(f)
        arrays = malloc_arrays(x)
        s = sum((usable(a) - sizeof(a) for a in arrays); init = 0)
        push!(bytes, b)
        push!(slack, s)
        push!(residual, b - s)
        isodd(i) && push!(keep, x)
        i % 16 == 0 && GC.gc(false)
    end
    println("[usable] ", name)
    println("    bytes:              ", distribution(bytes))
    println("    usable slack:       ", distribution(slack))
    println("    bytes - slack:      ", distribution(residual))
    return bytes, residual
end

function inventory(f)
    f()
    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate = 1 f()
    allocs = Profile.Allocs.fetch().allocs
    return sort!([(string(a.type), a.size) for a in allocs]), sum(a.size for a in allocs)
end

function inventories(name, f)
    seen = Set{Vector{Tuple{String,Int}}}()
    sums = Int[]
    events = Int[]
    for _ = 1:8
        inv, total = inventory(f)
        push!(seen, inv)
        push!(sums, total)
        push!(events, @allocations f())
    end
    println("[usable] ", name, " distinct profiled inventories=", length(seen),
        " requested sums: ", distribution(sums), " allocation events: ", distribution(events))
end

println("[usable] Julia ", VERSION, " ", Sys.MACHINE, " threads=", Threads.nthreads(:default),
    ",", Threads.nthreads(:interactive), " check-bounds=", Base.JLOptions().check_bounds)

keep = Any[]
report("Vector{ComplexF64}(undef, 128)", () -> Vector{ComplexF64}(undef, 128); keep)
report("Vector{Float64}(undef, 1002)", () -> Vector{Float64}(undef, 1002); keep)
report("MersenneTwister(1)", () -> MersenneTwister(1); keep)

whistler(B0, n0, L) = () -> run_whistler(; n = 128, B0, n0, L, nsteps = 1, dt = 1e-6)
base = whistler(1.0, 1.0, 8.0)
variant = whistler(1, 1, 8)
report("run_whistler n=128 baseline Float64", base; keep)
report("run_whistler n=128 variant Int", variant; keep)
inventories("run_whistler n=128 baseline Float64", base)
inventories("run_whistler n=128 variant Int", variant)

# The permanent test's single-shot comparison, repeated.
let mismatches = 0, comparisons = 0
    for round = 1:20, n in (8, 128), L in (2π, 8, 8.0)
        b = () -> run_whistler(; n, B0 = 1.0, n0 = 1.0, L = Float64(L), nsteps = 1, dt = 1e-6)
        b()
        ceiling = @allocated b()
        for (B0, n0) in ((1, 1), (1, 1.0), (1.0, 1), (UInt128(1), Int8(1)))
            g = () -> run_whistler(; n, B0, n0, L, nsteps = 1, dt = 1e-6)
            g()
            bytes = @allocated g()
            comparisons += 1
            if bytes != ceiling
                mismatches += 1
                mismatches <= 10 && println("[usable] single-shot mismatch round=$round n=$n L=$(repr(L)) B0=$(repr(B0)) n0=$(repr(n0)) bytes=$bytes ceiling=$ceiling diff=$(bytes - ceiling)")
            end
        end
    end
    println("[usable] single-shot comparisons=", comparisons, " mismatches=", mismatches)
end
length(keep) > 0 || error("unreachable")
