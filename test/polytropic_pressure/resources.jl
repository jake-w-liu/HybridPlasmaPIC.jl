module ResourceChecks
using HybridPlasmaPIC, Test, Random, Statistics
import ..Oracle
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function public_pressure(::Type{S}, density::T, closure) where {S,T}
    source = T[density]
    destination = Vector{S}(undef, 1)
    electron_pressure!(destination, source, closure)
    return only(destination)
end
public_candidate(density::T, closure) where {T} = public_pressure(T, density, closure)
apply!(pressure, density, closure) = electron_pressure!(pressure, density, closure)
function sample(f, pe, n, c)
    return @timed f(pe, n, c)
end
println((VERSION, threads = Threads.nthreads()))
flush(stdout)
@testset "Ordinary pressure whole-grid resources and sampled independent powers" begin
    for T in (Float32, Float64), count in (4096, 4194304), gamma in (one(T), T(5 / 3))
        n = [T(0.5) + T(1.5) * T((i - 1) % 1024) / T(1024) for i = 1:count]
        before = copy(n)
        closure = PolytropicElectrons(T(0.7), T(1.25), gamma)
        old, new = similar(n), similar(n)
        sample(H._legacy_polytropic_pressure!, old, n, closure)
        sample(apply!, new, n, closure)
        for i in (1, 256, 512, 1024)
            lo, hi, e = H._polytropic_power_bounds(Q(n[i]) / Q(closure.n0), Q(gamma), 160)
            lo *= Q(closure.pe0) * Q(2)^Int(e)
            hi *= Q(closure.pe0) * Q(2)^Int(e)
            reference = HybridPlasmaPIC._from_exact(T, (lo + hi) / 2)
            radius =
                4max(
                    abs(Q(nextfloat(reference)) - Q(reference)),
                    abs(Q(reference) - Q(prevfloat(reference))),
                )
            @test lo - radius <= Q(new[i]) <= hi + radius
        end
        old_times, new_times = Float64[], Float64[]
        for repetition = 1:3
            GC.gc()
            if isodd(repetition)
                a = sample(H._legacy_polytropic_pressure!, old, n, closure)
                b = sample(apply!, new, n, closure)
            else
                b = sample(apply!, new, n, closure)
                a = sample(H._legacy_polytropic_pressure!, old, n, closure)
            end
            @test a.bytes == 0
            @test b.bytes == 0
            push!(old_times, a.time)
            push!(new_times, b.time)
        end
        @test n == before
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
@testset "Complete exceptional pressure resources at actual array sizes" begin
    for T in (Float32, Float64), count in (32, 256), sign in (-1, 1)
        k = T === Float32 ? 100 : 600
        density = ldexp(one(T), sign * k)
        n = fill(density, count)
        pe = similar(n)
        closure = PolytropicElectrons(inv(density), one(T), T(2))
        apply!(pe, n, closure)
        measurements = [sample(apply!, pe, n, closure) for _ = 1:3]
        @test all(==(density), pe)
        @test all(==(density), n)
        @test all(x -> x.value === pe, measurements)
        println((
            T,
            count,
            sign,
            bytes = [x.bytes for x in measurements],
            seconds = [x.time for x in measurements],
        ))
        flush(stdout)
    end
end
println("POLYTROPIC_ORDINARY_RESOURCES_COMPLETE")

end
