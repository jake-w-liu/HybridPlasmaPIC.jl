module IsothermalLogRangeChecks
using HybridPlasmaPIC, Test, Random
include("isothermal_energy/oracles.jl")
using .Oracle: Q, oracle_log
@testset "Range-reduced logarithm intervals enclose independent exact bounds" begin
    rng = MersenneTwister(0x3189)
    for T in (Float32, Float64)
        U = T === Float32 ? UInt32 : UInt64
        values = T[
            nextfloat(zero(T)),
            prevfloat(floatmin(T)),
            floatmin(T),
            nextfloat(floatmin(T)),
            prevfloat(T(0.5)),
            T(0.5),
            nextfloat(T(0.5)),
            prevfloat(one(T)),
            one(T),
            nextfloat(one(T)),
            prevfloat(T(2)),
            T(2),
            nextfloat(T(2)),
            floatmax(T),
        ]
        append!(values, [reinterpret(T, rand(rng, U) & (typemax(U) >> 1)) for i = 1:64])
        filter!(x -> isfinite(x) && x > 0, values)
        for x in values
            valid, lo, hi = HybridPlasmaPIC._isothermal_logarithm_bounds(x)
            @test valid
            a, b = oracle_log(Q(x); bits = 192)
            @test Q(lo) <= a <= b <= Q(hi)
            @test isfinite(lo) && isfinite(hi)
            if x < T(0.5) || x > T(2)
                @test abs(hi - lo) < max(abs(lo), abs(hi)) / 2
            end
        end
        for x in (T(0), -one(T), T(NaN), T(Inf), -T(Inf))
            @test !first(HybridPlasmaPIC._isothermal_logarithm_bounds(x))
        end
    end
end
println("ISOTHERMAL_REDUCED_LOG_BOUNDS_COMPLETE")

function measured(n, te, dx)
    g = FourierGrid(size(n), ntuple(i -> eltype(n)(size(n, i)) * dx[i], ndims(n)))
    closure = IsothermalElectrons(te)
    electron_internal_energy(n, closure, g)
    return @timed electron_internal_energy(n, closure, g)
end
println((VERSION,))
flush(stdout)
@testset "Extreme-density normal integrals retain zero whole-call allocation" begin
    for T in (Float32, Float64), count in (256, 4194304)
        for (density, te) in (
            (T === Float32 ? T(1e20) : T(1e200), T === Float32 ? T(1e-20) : T(1e-200)),
            (floatmin(T), inv(floatmin(T))),
        )
            n = fill(density, count)
            dx = (one(T),)
            @test HybridPlasmaPIC._isothermal_general_probe(n, te, dx).certified
            result = measured(n, te, dx)
            lo, hi = oracle_log(Q(density); bits = 256)
            factor = Q(count) * Q(te) * Q(density)
            lo *= factor
            hi *= factor
            reference = HybridPlasmaPIC._from_exact(T, (lo + hi) / 2)
            if T === Float32
                @test isequal(result.value, reference)
            else
                radius =
                    2max(
                        abs(Q(nextfloat(reference)) - Q(reference)),
                        abs(Q(reference) - Q(prevfloat(reference))),
                    )
                @test lo - radius <= Q(result.value) <= hi + radius
            end
            @test result.bytes == 0
            @test first(n) == last(n) == density
            println((
                T,
                count,
                density,
                te,
                value = result.value,
                reference,
                bytes = result.bytes,
                seconds = result.time,
            ))
            flush(stdout)
        end
    end
end
println("ISOTHERMAL_EXTREME_RESOURCES_COMPLETE")

end
