module OrdinaryChecks
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
function accurate(value::T, lower::Q, upper::Q) where {T}
    reference = HybridPlasmaPIC._from_exact(T, (lower + upper) / 2)
    isfinite(reference) && !iszero(reference) || return isequal(value, reference)
    isfinite(value) || return false
    radius =
        4max(
            abs(Q(nextfloat(reference)) - Q(reference)),
            abs(Q(reference) - Q(prevfloat(reference))),
        )
    return lower - radius <= Q(value) <= upper + radius
end
println((VERSION, threads = Threads.nthreads()))
flush(stdout)
@testset "First-order correction has a bounded nonlinear remainder" begin
    epsilon = Q(eps(Float64))
    u = epsilon / 2
    computed_limit = Q(sqrt(eps(Float64)) / 4)
    # Two divisions and one product, with a conservative allowance for a final
    # underflow, bound the exact gamma*relative-residual represented by the guard.
    t = (computed_limit + Q(nextfloat(0.0))) / (1 - u)^3
    log_remainder = t * epsilon / (2(1 - epsilon))
    y = t + log_remainder
    full_remainder = log_remainder + y^2 / (2(1 - y))
    @test y < 1
    @test full_remainder < epsilon / 16
    for gamma in (Q(1, 2), Q(1), Q(2), Q(3), Q(17)),
        e in (-epsilon, -epsilon / 3, epsilon / 3, epsilon)

        t = gamma * e
        @test abs(t) <= computed_limit
        if denominator(gamma) == 1
            exact = (1 + e)^Int(gamma) - 1 - t
            @test abs(exact) <= full_remainder
        else
            # Squaring checks the independent rational interval for sqrt(1+e).
            center = 1 + t
            @test (center - full_remainder)^2 <= 1 + e <= (center + full_remainder)^2
        end
    end
end
@testset "Ordinary pressures retain complete rational factors" begin
    for T in (Float32, Float64),
        C in (Float32, Float64),
        a = 1:6,
        b = 1:5,
        q in (1, 2),
        p in (1, 2, 3, 7)

        n = T(a^q)
        closure = PolytropicElectrons(C(0.7), C(b^q), C(p / q))
        expected = Q(closure.pe0) * Q(a, b)^p
        result = public_candidate(n, closure)
        @test H._polytropic_pressure_probe(T, n, closure) !== nothing
        @test accurate(result, expected, expected)
        @test result isa T
    end
end
@testset "Independent integer-root intervals cover nonrational powers" begin
    for T in (Float32, Float64), C in (Float32, Float64), m = 2:7, q in (2, 4, 8), p in (1, 3)
        scale = BigInt(1) << 192
        target = BigInt(m) * scale^q
        root = target
        for _ = 1:trailing_zeros(q)
            root = isqrt(root)
        end
        @test root^q <= target < (root + 1)^q
        c = PolytropicElectrons(C(0.7), one(C), C(p / q))
        lo, hi = Q(c.pe0) * (root // scale)^p, Q(c.pe0) * ((root + 1) // scale)^p
        @test accurate(public_candidate(T(m), c), lo, hi)
    end
end
@testset "Complete pointwise range and density flooring" begin
    for T in (Float32, Float64)
        k = T === Float32 ? 100 : 600
        large, small = ldexp(one(T), k), ldexp(one(T), -k)
        for (density, closure, expected) in (
            (large, PolytropicElectrons(small, one(T), T(2)), Q(large)),
            (small, PolytropicElectrons(large, one(T), T(2)), Q(small)),
            (large, PolytropicElectrons(one(T), small, T(0.5)), Q(large)),
        )
            @test Oracle.rounded(public_candidate(density, closure), expected, expected)
            @test public_candidate(zero(T), closure) == 0
            @test public_candidate(-one(T), closure) == 0
        end
        c = PolytropicElectrons(one(T), prevfloat(T(2)), ldexp(one(T), precision(T)))
        lo, hi, e = H._polytropic_power_bounds(Q(2) / Q(c.n0), Q(c.γ), 160)
        lo *= Q(2)^Int(e)
        hi *= Q(2)^Int(e)
        @test 2 < lo <= hi < 4
        @test Oracle.rounded(public_candidate(T(2), c), lo, hi)
        @test public_candidate(floatmax(T), PolytropicElectrons(zero(T), one(T), floatmax(T))) == 0
        for density in (T(NaN), T(Inf), -T(Inf))
            c = PolytropicElectrons(T(0.7), one(T), T(1.5))
            expected = c.pe0 * (max(density, zero(T)) / c.n0)^c.γ
            @test isequal(public_candidate(density, c), expected)
        end
    end
end
println("POLYTROPIC_ORDINARY_FOCUSED_COMPLETE")

end
