module Oracle
using HybridPlasmaPIC
const Q = Rational{BigInt}
# Independent alternating log(1+x) series: consecutive partial sums enclose
# the answer for 0 <= x <= 1/2. No atanh series or floating logarithm is used.
function alternating_log(x::Q, bits::Int)
    @assert 0 <= x <= Q(1, 2)
    iszero(x) && return zero(Q), zero(Q)
    target = BigInt(1) // (BigInt(1) << bits)
    partial = zero(Q)
    power = x
    k = 1
    while true
        partial += (isodd(k) ? power : -power) / k
        nextpower = power * x
        error = nextpower / (k + 1)
        if error <= target
            return isodd(k) ? (partial - error, partial) : (partial, partial + error)
        end
        power = nextpower
        k += 1
    end
end

function oracle_log(q::Q; bits = 768)
    @assert q > 0
    l3, u3 = alternating_log(Q(1, 2), bits)
    l4, u4 = alternating_log(Q(1, 3), bits)
    e = ndigits(numerator(q); base = 2) - ndigits(denominator(q); base = 2)
    m = q / Q(2)^e
    if m < 1
        m *= 2
        e -= 1
    end
    @assert 1 <= m < 2
    if m <= Q(3, 2)
        lo, hi = alternating_log(m - 1, bits)
    else
        lo, hi = alternating_log(m / Q(3, 2) - 1, bits)
        lo += l3
        hi += u3
    end
    return lo + e * (e >= 0 ? l3 + l4 : u3 + u4), hi + e * (e >= 0 ? u3 + u4 : l3 + l4)
end

function oracle_integral(n, temperature, dx)
    factor = Q(temperature) * prod(Q, dx)
    lower = upper = zero(Q)
    for value in n
        value > 0 || continue
        a = Q(value) * factor
        lo, hi = oracle_log(Q(value))
        lower += a * (a >= 0 ? lo : hi)
        upper += a * (a >= 0 ? hi : lo)
    end
    return lower, upper
end

function rounded(value::T, lo, hi) where {T}
    lo <= hi || return false
    isfinite(value) || isinf(value) || return false
    edge = Q(floatmax(T)) + (Q(floatmax(T)) - Q(prevfloat(floatmax(T)))) / 2
    if isinf(value)
        return signbit(value) ? hi <= -edge : lo >= edge
    end
    if signbit(value)
        value = abs(value)
        lo, hi = -hi, -lo
    end
    center = Q(value)
    lower = (Q(prevfloat(value)) + center) / 2
    upper = isfinite(nextfloat(value)) ? (center + Q(nextfloat(value))) / 2 : edge
    U = T === Float32 ? UInt32 : UInt64
    even = iseven(reinterpret(U, value))
    return (lower < lo || even && lower == lo) && (hi < upper || even && hi == upper)
end
function exact_density_moments(n::Vector{Float32})
    scale = Int128(1) << 24
    s1 = s2 = s3 = s4 = s5 = s6 = Int128(0)
    largest = Int128(0)
    @assert length(n) <= 1 << 22
    for value in n
        @assert 0.5f0 <= value <= 2.0f0
        d = Int128(ldexp(value, 24)) - scale
        # At most 2^22 terms of magnitude (2^15)^6 fit within Int128.
        @assert abs(d) <= Int128(1) << 15
        largest = max(largest, abs(d))
        d2 = d * d
        d3 = d2 * d
        d4 = d3 * d
        d5 = d4 * d
        d6 = d5 * d
        s1 += d
        s2 += d2
        s3 += d3
        s4 += d4
        s5 += d5
        s6 += d6
    end
    return (s1, s2, s3, s4, s5, s6), largest
end
function moment_eight_bounds(n::Vector{Float32}, temperature, dx)
    moments, maximum_delta = exact_density_moments(n)
    seventh = eighth = BigInt(0)
    scale = BigInt(1) << 24
    for first = 1:64:length(n)
        s7 = s8 = Int128(0)
        for i = first:min(first + 63, length(n))
            d = Int128(ldexp(n[i], 24)) - (Int128(1) << 24)
            @assert abs(d) <= Int128(1) << 15
            # Each eighth power is <=2^120; at most 64 fit below 2^127.
            d7 = d^7
            s7 += d7
            s8 += d7 * d
        end
        seventh += s7
        eighth += s8
    end
    exact_moments = (BigInt.(moments)..., seventh, eighth)
    polynomial = Q(exact_moments[1], scale)
    for k = 2:8
        polynomial += (iseven(k) ? 1 : -1) * Q(exact_moments[k], scale^k * k * (k - 1))
    end
    maximum = Q(maximum_delta, scale)
    tail = maximum * Q(eighth, scale^8) / (72 * (1 - maximum))
    factor = Q(temperature) * prod(Q, dx)
    @assert factor >= 0
    return factor * (polynomial - tail), factor * (polynomial + tail)
end

function complete_bounds(
    n::AbstractArray{T},
    temperature,
    dx,
    state,
    bits,
    degree::Val{K},
) where {T,K}
    @assert 40 <= bits <= state.cap
    lo, hi = HybridPlasmaPIC._isothermal_scaled_polynomial(n, state, bits, degree)
    magnitude = Q(2)^(state.q - state.p)
    tail = Q(state.active) * magnitude^(K - 1) / (K * (K + 1) * (1 - magnitude))
    scale = BigInt(1) << bits
    linear = Q(state.delta, state.scale)
    lower = linear + magnitude^2 * (Q(lo, scale) - tail)
    upper = linear + magnitude^2 * (Q(hi, scale) + tail)
    multiplier = Q(temperature) * prod(Q, dx)
    return multiplier >= 0 ? (multiplier * lower, multiplier * upper) :
           (multiplier * upper, multiplier * lower)
end

function interval(
    n::AbstractArray{T},
    temperature::C,
    dx;
    bits = precision(T) + 16,
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    bits >= 40 || throw(ArgumentError("near-unit interval requires at least 40 fraction bits"))
    isfinite(temperature) && all(isfinite, dx) || return nothing
    state = HybridPlasmaPIC._isothermal_statistics(n)
    state === nothing && return nothing
    (iszero(state.active) || iszero(temperature) || any(iszero, dx)) && return (zero(Q), zero(Q))
    bits <= state.cap || return nothing
    return complete_bounds(n, temperature, dx, state, bits, Val(T === Float32 ? 8 : 12))
end


function evaluate(n, temperature, dx)
    isfinite(temperature) && all(isfinite, dx) || return nothing
    state = HybridPlasmaPIC._isothermal_statistics(n)
    state === nothing && return nothing
    return HybridPlasmaPIC._isothermal_nearunit_energy(n, temperature, dx, state)
end
end
