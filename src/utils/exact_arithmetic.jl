# Exact conversions used by cold numerical fallbacks. Callers supply finite
# floating-point operands and finite Rational{BigInt} results.

_exact_rational(value) = Rational{BigInt}(value)
function _exact_rational(value::BigFloat)
    # Julia's Rational conversion scales at the current BigFloat precision.
    # Preserve stored bits if the caller has lowered its working precision.
    return setprecision(BigFloat, precision(value)) do
        Rational{BigInt}(value)
    end
end

_from_exact(::Type{T}, value::Rational{BigInt}) where {T<:AbstractFloat} = T(value)

function _from_exact(::Type{BigFloat}, value::Rational{BigInt})
    n, d = numerator(value), denominator(value)
    # Both operands are exact. MPFR division then rounds once at the caller's
    # working precision and rounding mode, without changing either setting.
    nfloat = BigFloat(n; precision = max(2, ndigits(n; base = 2)))
    dfloat = BigFloat(d; precision = max(2, ndigits(d; base = 2)))
    return nfloat / dfloat
end

function _from_exact(::Type{T}, value::Rational{BigInt}) where {T<:Union{Float16,Float32,Float64}}
    signed_n, d = numerator(value), denominator(value)
    iszero(signed_n) && return zero(T)
    negative = signed_n < 0
    n = abs(signed_n)
    p = precision(T)
    maximum_exponent = exponent(floatmax(T))
    minimum_quantum = exponent(floatmin(T)) - p + 1

    # Bit lengths bracket the binade exponent within one. Resolve it with an
    # exact comparison; early extremes avoid unnecessary large integer shifts.
    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    if e > maximum_exponent + 1
        return negative ? -T(Inf) : T(Inf)
    elseif e < minimum_quantum - 1
        return negative ? -zero(T) : zero(T)
    end
    if e >= 0 ? n < d << e : n << -e < d
        e -= 1
    end
    e > maximum_exponent && return negative ? -T(Inf) : T(Inf)

    # Subnormals use one fixed quantum. Quotient/remainder comparison rounds
    # to nearest, ties to even, without a BigFloat or narrower-float intermediate.
    quantum = max(e - p + 1, minimum_quantum)
    scaled_n, scaled_d = quantum >= 0 ? (n, d << quantum) : (n << -quantum, d)
    quotient, remainder = divrem(scaled_n, scaled_d)
    twice_remainder = remainder << 1
    if twice_remainder > scaled_d || (twice_remainder == scaled_d && isodd(quotient))
        quotient += 1
    end
    # The rounded integer is at most 2^p and therefore exactly representable.
    # Rescaling is exact unless the rounded result is the overflow endpoint.
    magnitude = ldexp(T(quotient), quantum)
    return negative ? -magnitude : magnitude
end
