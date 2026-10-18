module BinomialChecks

using Test, HybridPlasmaPIC, Random
const Q = Rational{BigInt}
const H = HybridPlasmaPIC

# NIST DLMF 4.6.7: https://dlmf.nist.gov/4.6#E7, for abs(base-1) < 1.
# After t_k, every ratio after the next term is bounded by
# abs(z)*max(1, gamma/(k+2)); a ratio below one bounds the remaining
# absolute tail geometrically. Rational arithmetic retains this enclosure.
# The requested 144-bit tail is checked by requiring both endpoints to round
# identically in each tested destination precision.
function binomial_bounds(base::Q, gamma::Q, bits::Int = 144)
    0 < base < 2 || throw(ArgumentError("oracle base must lie in (0,2)"))
    gamma >= 0 || throw(ArgumentError("oracle exponent must be nonnegative"))
    bits > 0 || throw(ArgumentError("oracle precision must be positive"))
    z = base - 1
    tolerance = BigInt(1) // (BigInt(1) << bits)
    partial = term = one(Q)
    k = 0
    while true
        next = term * (gamma - k) * z / (k + 1)
        iszero(next) && return partial, partial, k
        ratio = abs(z) * max(one(Q), gamma / (k + 2))
        if ratio < 1
            tail = abs(next) / (1 - ratio)
            tail <= tolerance && return max(zero(Q), partial - tail), partial + tail, k
        end
        partial += next
        term = next
        k = Base.checked_add(k, 1)
    end
end

# Positive float encodings are numerically ordered. Search the complete finite
# range; compare rational distances and encoding parity at the final midpoint.
function independent_round(::Type{T}, x::Q) where {T<:Union{Float32,Float64}}
    x >= 0 || throw(ArgumentError("positive oracle output required"))
    U = T === Float32 ? UInt32 : UInt64
    largest = reinterpret(U, floatmax(T))
    maximum = Q(floatmax(T))
    overflow = maximum + (maximum - Q(prevfloat(floatmax(T)))) / 2
    x >= overflow && return T(Inf)
    x >= maximum && return floatmax(T)
    lo, hi = zero(U), largest
    while hi - lo > 1
        mid = lo + ((hi - lo) >> 1)
        if Q(reinterpret(T, mid)) <= x
            lo = mid
        else
            hi = mid
        end
    end
    a, b = reinterpret(T, lo), reinterpret(T, hi)
    midpoint = (Q(a) + Q(b)) / 2
    return x < midpoint ? a : x > midpoint ? b : iseven(lo) ? a : b
end

@testset "Independent binomial and rounding oracle controls" begin
    for x in Q[1//2, 3//4, 1, 5//4, 3//2], p in (0, 1, 2, 3, 8, 97)
        lo, hi, k = binomial_bounds(x, Q(p))
        @test lo <= x^p <= hi
        @test hi - lo <= Q(1, BigInt(1) << 143)
    end
    for x in Q[1//2, 3//4, 1, 5//4, 3//2]
        lo, hi, k = binomial_bounds(x, Q(1, 2))
        @test 0 <= lo <= hi
        @test lo^2 <= x <= hi^2
    end
    for T in (Float32, Float64),
        value in T[0, nextfloat(zero(T)), floatmin(T), 0.5, 1, 1.5, prevfloat(floatmax(T))]

        adjacent = nextfloat(value)
        midpoint = (Q(value) + Q(adjacent)) / 2
        U = T === Float32 ? UInt32 : UInt64
        @test independent_round(T, Q(value)) === value
        @test independent_round(T, midpoint) === (iseven(reinterpret(U, value)) ? value : adjacent)
        @test independent_round(T, (Q(value) + midpoint) / 2) === value
        @test independent_round(T, (midpoint + Q(adjacent)) / 2) === adjacent
    end
    for T in (Float32, Float64)
        maximum = Q(floatmax(T))
        gap = maximum - Q(prevfloat(floatmax(T)))
        @test independent_round(T, maximum + gap / 4) === floatmax(T)
        @test independent_round(T, maximum + gap / 2) === T(Inf)
    end
end
flush(stdout)

function check_pressure(::Type{S}, density::T, c::PolytropicElectrons{C}) where {S,T,C}
    lower, upper, terms = binomial_bounds(Q(density) / Q(c.n0), Q(c.γ))
    lower *= Q(c.pe0)
    upper *= Q(c.pe0)
    left, right = independent_round(S, lower), independent_round(S, upper)
    @test isequal(left, right)
    destination = S[-1]
    @test electron_pressure!(destination, T[density], c) === destination
    result = only(destination)
    if iszero(left) || !isfinite(left)
        @test isequal(result, left)
    else
        radius = 4max(Q(nextfloat(left)) - Q(left), Q(left) - Q(prevfloat(left)))
        @test isfinite(result)
        @test lower - radius <= Q(result) <= upper + radius
    end
    return terms
end
@testset "Stored non-integer gamma against independent binomial intervals" begin
    maxterms = 0
    for S in (Float32, Float64),
        T in (Float32, Float64),
        C in (Float32, Float64),
        density in T[0.75, 1, 1.25, 1.5],
        gamma in C[0.1, prevfloat(one(C)), nextfloat(one(C)), 1.4, 5/3, 2.5, 3.7]

        c = PolytropicElectrons(C(0.7), C(1.25), gamma)
        maxterms = max(maxterms, check_pressure(S, density, c))
    end
    rng = MersenneTwister(682411)
    for T in (Float32, Float64), i = 1:32
        reference = T(0.75) + T(rand(rng)) / 2
        density = reference * (T(0.75) + T(rand(rng)) / 2)
        gamma = T(0.125) + T(8) * T(rand(rng))
        maxterms =
            max(maxterms, check_pressure(T, density, PolytropicElectrons(T(0.7), reference, gamma)))
    end
    for S in (Float32, Float64), T in (Float32, Float64), direction in (-1, 1)
        density = direction < 0 ? prevfloat(one(T)) : nextfloat(one(T))
        gamma = ldexp(one(T), precision(T))
        maxterms =
            max(maxterms, check_pressure(S, density, PolytropicElectrons(T(0.7), one(T), gamma)))
    end
    println((maximum_series_terms = maxterms,))
    flush(stdout)
end
println("INDEPENDENT_BINOMIAL_PRESSURE_COMPLETE")
flush(stdout)

end
