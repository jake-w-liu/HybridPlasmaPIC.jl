module BinomialOracle
const Q = Rational{BigInt}

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

end
