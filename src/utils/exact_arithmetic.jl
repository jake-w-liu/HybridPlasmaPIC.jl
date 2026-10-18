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

# Fixed-width dyadic conversion used by bounded exact diagnostic intervals.
@inline _exact_int128_magnitude(x::Int128) = x < 0 ? -unsigned(x) : unsigned(x)
function _from_exact(::Type{T}, value::Int128, exponent::Int) where {T<:Union{Float32,Float64}}
    iszero(value) && return zero(T)
    integer = _exact_int128_magnitude(value)
    count = 128 - leading_zeros(integer)
    actual_exponent = Int128(exponent) + count - 1
    minimum_quantum = Base.exponent(floatmin(T)) - precision(T) + 1
    actual_exponent > Base.exponent(floatmax(T)) && return copysign(T(Inf), value)
    actual_exponent < minimum_quantum - 1 && return copysign(zero(T), value)
    quantum = max(exponent + count - precision(T), minimum_quantum)
    shift = quantum - exponent
    rounded = if shift <= 0
        integer << -shift
    elseif shift >= 128
        # Int128 magnitude is at most 2^127; the only half-way case rounds to zero.
        UInt128(0)
    else
        whole = integer >> shift
        rest = integer & ((UInt128(1) << shift) - 1)
        half = UInt128(1) << (shift - 1)
        whole + (rest > half || (rest == half && isodd(whole)))
    end
    result = ldexp(T(rounded), quantum)
    return value < 0 ? -result : result
end
