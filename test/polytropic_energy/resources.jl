module IntegratedEnergyResourceChecks
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
include(joinpath(pkgdir(H), "test/isothermal_energy/oracles.jl"))
function reference(n, c, dx, rounding)
    setprecision(BigFloat, 512) do
        setrounding(BigFloat, rounding) do
            total = BigFloat(0)
            for density in n
                density > 0 || continue
                x = BigFloat(density)
                total += c.γ == 1 ? x * log(x) : (x / BigFloat(c.n0))^BigFloat(c.γ)
            end
            total *= BigFloat(c.pe0)
            total /= c.γ == 1 ? BigFloat(c.n0) : BigFloat(c.γ) - 1
            for spacing in dx
                total *= BigFloat(spacing)
            end
            return Q(total)
        end
    end
end
function accurate(value::T, lo, hi) where {T}
    target = H._from_exact(T, (lo + hi) / 2)
    radius = 4max(abs(Q(nextfloat(target)) - Q(target)), abs(Q(target) - Q(prevfloat(target))))
    return isfinite(value) && lo - radius <= Q(value) <= hi + radius
end
measure(n::N, c::C, g::G) where {N,C,G} = @timed electron_internal_energy(n, c, g)
println((VERSION, root = pkgdir(H)))
flush(stdout)
@testset "Public integrated energy correctness and resource measurements" begin
    for T in (Float32, Float64),
        kind in (:ordinary, :nearunit, :repeated, :distinct),
        count in (32, 256, 4096)

        g = FourierGrid((count,), (one(T),))
        c = PolytropicElectrons(T(0.7), T(3), kind === :ordinary ? T(1.5) : one(T))
        n =
            kind === :ordinary ? fill(T(1.25), count) :
            kind === :nearunit ?
            [one(T) + (isodd(i) ? -one(T) : one(T)) * ldexp(one(T), -10) for i = 1:count] :
            kind === :repeated ? [isodd(i) ? T(0.5) : T(1.5) for i = 1:count] :
            [one(T) + T(i) / T(2count) for i = 1:count]
        before = copy(n)
        lo = reference(n, c, g.dx, RoundDown)
        hi = reference(n, c, g.dx, RoundUp)
        @test lo <= hi
        value = electron_internal_energy(n, c, g)
        @test kind === :ordinary ? accurate(value, lo, hi) : Oracle.rounded(value, lo, hi)
        measure(n, c, g)
        sample = measure(n, c, g)
        @test isequal(sample.value, value)
        @test n == before
        @test sample.bytes == 0
        println((
            T,
            kind,
            count,
            value,
            bytes = sample.bytes,
            allocations = Base.gc_alloc_count(sample.gcstats),
            seconds = sample.time,
        ))
        flush(stdout)
    end
end
println("PUBLIC_POLYTROPIC_ENERGY_RESOURCES_COMPLETE")

end
