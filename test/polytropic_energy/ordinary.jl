module OrdinaryChecks
using Test, HybridPlasmaPIC, Random
include("../polytropic_pressure/binomial_oracle.jl")
const P = (
    probe = HybridPlasmaPIC._polytropic_energy_probe,
    exponent_bounds = HybridPlasmaPIC._polytropic_energy_exponent_bounds,
)
const B = BinomialOracle
const Q = Rational{BigInt}
function round_signed(::Type{T}, value::Q) where {T}
    return value < 0 ? -B.independent_round(T, -value) : B.independent_round(T, value)
end
function distance(a::T, b::T) where {T<:Union{Float32,Float64}}
    U = T === Float32 ? UInt32 : UInt64
    return abs(BigInt(reinterpret(U, abs(a))) - BigInt(reinterpret(U, abs(b))))
end
println((VERSION, root = pkgdir(HybridPlasmaPIC)))
flush(stdout)
@testset "Ordinary range certificate against exact rational powers" begin
    for T in (Float32, Float64),
        C in (Float32, Float64),
        D = 1:3,
        (power, root) in ((1, 2), (3, 2), (5, 4), (7, 4), (2, 1), (3, 1), (7, 1)),
        count in (3, 7, 16)

        values = [mod1(i, 3) for i = 1:count]
        n = reshape(T[a^root for a in values], ntuple(i -> i == 1 ? count : 1, D))
        c = PolytropicElectrons(C(0.7), C(2^root), C(power / root))
        dx = ntuple(_ -> T(0.5), D)
        exact = Q(c.pe0) * sum(a -> Q(a, 2)^power, values) * prod(Q, dx) / (Q(c.γ) - 1)
        lo, hi = P.exponent_bounds(maximum(n), count, c, dx)
        @test Q(2)^floor(Int, lo) <= abs(exact) <= Q(2)^ceil(Int, hi)
        value = P.probe(n, c, dx)
        @test value !== nothing
        if value !== nothing
            expected = round_signed(T, exact)
            @test signbit(value) == signbit(expected)
            @test distance(value, expected) <= 4
        end
    end
end
flush(stdout)
@testset "Stored-exponent binomial integral oracle" begin
    rng = MersenneTwister(613775)
    for T in (Float32, Float64),
        C in (Float32, Float64),
        D = 1:3,
        gamma in C[0.1, 0.7, prevfloat(one(C)), nextfloat(one(C)), 1.4, 5/3, 2.5, 3.7]

        n = reshape(T.(0.5 .+ rand(rng, 7)), ntuple(i -> i == 1 ? 7 : 1, D))
        c = PolytropicElectrons(C(0.7), one(C), gamma)
        dx = ntuple(_ -> T(0.5), D)
        lower = upper = zero(Q)
        for density in n
            lo, hi, _ = B.binomial_bounds(Q(density), Q(gamma), 168)
            lower += lo
            upper += hi
        end
        factor = Q(c.pe0) * prod(Q, dx) / (Q(gamma) - 1)
        lo, hi = factor > 0 ? (lower * factor, upper * factor) : (upper * factor, lower * factor)
        left, right = round_signed(T, lo), round_signed(T, hi)
        @test isequal(left, right)
        value = P.probe(n, c, dx)
        @test value !== nothing
        if value !== nothing
            @test signbit(value) == signbit(left)
            @test distance(value, left) <= 4
        end
    end
end
flush(stdout)
@testset "Range margins and complete-input fallback boundaries" begin
    for T in (Float32, Float64)
        c = PolytropicElectrons(one(T), one(T), T(2))
        k = T === Float32 ? 100 : 600
        for sign in (-1, 1)
            density = ldexp(one(T), sign * k)
            dx = (ldexp(one(T), -sign * k),)
            exact = Q(density)^2 * Q(only(dx))
            lo, hi = P.exponent_bounds(density, 1, c, dx)
            @test Q(2)^floor(Int, lo) <= exact <= Q(2)^ceil(Int, hi)
            # The native probe may reject even a finite complete result; it must
            # never accept a materialized zero or infinity as a range repair.
            value = P.probe(T[density], c, dx)
            @test value === nothing || isequal(value, round_signed(T, exact))
        end
        for density in (sqrt(floatmax(T)), nextfloat(zero(T)), floatmin(T))
            @test P.probe(T[density], c, (one(T),)) === nothing
        end
        for n in (T[Inf], T[NaN], T[-Inf], T[], T[-1, 0])
            @test P.probe(n, c, (one(T),)) === nothing
        end
        @test P.probe(T[1, 2], PolytropicElectrons(one(T), one(T), one(T)), (one(T),)) === nothing
        @test P.probe(T[1, 2], c, (zero(T),)) === nothing
        @test P.probe(T[1, 2], c, (T(Inf),)) === nothing
        for gamma in (prevfloat(1.0), nextfloat(1.0))
            value = P.probe(ones(T, 2), PolytropicElectrons(0.5, 1.0, gamma), (one(T),))
            @test value !== nothing && isequal(value, round_signed(T, Q(1) / (Q(gamma) - 1)))
        end
    end
    a = 1.0 + 2.0^-26
    b = 1.0 + 2.0^-27
    for k in (26, 27, 28)
        n = zeros(2, 2, 2)
        n[1] = 1
        n[2] = 0.5
        @test P.probe(n, PolytropicElectrons(a, 1.0, 2.0^(2k)), (2.0^k - 1, 2.0^k + 1, b)) ===
              nothing
    end
end
println("ORDINARY_POLYTROPIC_INTEGRAL_NUMERICAL_COMPLETE")
flush(stdout)


end
