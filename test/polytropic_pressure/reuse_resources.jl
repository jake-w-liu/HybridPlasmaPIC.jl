module ReuseResources
using HybridPlasmaPIC, Test, Random, Statistics
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function old_apply!(pe::Array{S}, n::Array{T}, c) where {S,T}
    @inbounds for i in eachindex(pe, n)
        pe[i] = H._polytropic_pressure(S, max(n[i], zero(T)), c)
    end
    return pe
end
new_apply!(pe, n, c) = electron_pressure!(pe, n, c)
sample(f, pe, n, c) = @timed f(pe, n, c)
allocation_signature(sample) = (
    sample.gcstats.malloc,
    sample.gcstats.realloc,
    sample.gcstats.poolalloc,
    sample.gcstats.bigalloc,
)
println((VERSION,))
flush(stdout)
@testset "Repeated exceptional pressure allocates one scalar evaluation per density run" begin
    for T in (Float32, Float64), sign in (-1, 1), kind in (:uniform, :blocks)
        k = T === Float32 ? 100 : 600
        base = ldexp(one(T), sign * k)
        palette = kind === :uniform ? T[base] : base .* T[1, 1.25, 1.5, 1.75]
        c = PolytropicElectrons(inv(base), one(T), T(2))
        expected = [HybridPlasmaPIC._from_exact(T, Q(c.pe0) * Q(x)^2) for x in palette]
        # The reference evaluates each density run once. Match every allocation
        # counter to that work; bytes include allocator-provided block capacity
        # and can differ even for identical BigInt requests (Julia gc-stock.c).
        scalar_n = copy(palette)
        scalar_pe = similar(scalar_n)
        old_apply!(scalar_pe, scalar_n, c)
        reference = sample(old_apply!, scalar_pe, scalar_n, c)
        @test isequal(scalar_pe, expected)
        expected_allocations = allocation_signature(reference)
        allocation_counts = typeof(expected_allocations)[]
        for count in (32, 256, 4194304)
            block = count ÷ length(palette)
            n = [palette[cld(i, block)] for i = 1:count]
            pe = similar(n)
            new_apply!(pe, n, c)
            samples = [sample(new_apply!, pe, n, c) for _ = 1:3]
            @test all(x -> x.value === pe, samples)
            @test all(i -> isequal(pe[i], expected[cld(i, block)]), eachindex(pe))
            @test all(i -> isequal(n[i], palette[cld(i, block)]), eachindex(n))
            @test all(x -> allocation_signature(x) == expected_allocations, samples)
            push!(allocation_counts, allocation_signature(samples[1]))
            if count <= 256
                previous = similar(n)
                old_apply!(previous, n, c)
                old = sample(old_apply!, previous, n, c)
                @test isequal(pe, previous)
                @test samples[1].bytes <= old.bytes
                println((T, sign, kind, count, old_bytes = old.bytes, old_seconds = old.time))
            end
            println((
                T,
                sign,
                kind,
                count,
                allocations = [allocation_signature(x) for x in samples],
                reference_allocations = expected_allocations,
                bytes = [x.bytes for x in samples],
                seconds = [x.time for x in samples],
            ))
            flush(stdout)
        end
        @test all(==(expected_allocations), allocation_counts)
    end
end
@testset "Varied ordinary pressure keeps scalar results and zero allocations" begin
    for T in (Float32, Float64), count in (4096, 4194304), gamma in (one(T), T(5 / 3))
        n = [T(0.5) + T(1.5) * T((i - 1) % 1024) / T(1024) for i = 1:count]
        c = PolytropicElectrons(T(0.7), T(1.25), gamma)
        old, new = similar(n), similar(n)
        old_apply!(old, n, c)
        new_apply!(new, n, c)
        @test isequal(old, new)
        old_times = Float64[]
        new_times = Float64[]
        for repetition = 1:5
            GC.gc()
            if isodd(repetition)
                a = sample(old_apply!, old, n, c)
                b = sample(new_apply!, new, n, c)
            else
                b = sample(new_apply!, new, n, c)
                a = sample(old_apply!, old, n, c)
            end
            @test a.bytes == 0
            @test b.bytes == 0
            @test isequal(old, new)
            push!(old_times, a.time)
            push!(new_times, b.time)
        end
        println((
            T,
            count,
            gamma,
            old_times,
            new_times,
            median_ratio = median(new_times) / median(old_times),
        ))
        flush(stdout)
    end
end
println("POLYTROPIC_ARRAY_REUSE_RESOURCES_COMPLETE")

end
