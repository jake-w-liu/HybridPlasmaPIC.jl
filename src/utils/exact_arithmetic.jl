# Exact conversions used by cold numerical fallbacks. Callers supply finite
# floating-point operands and finite Rational{BigInt} results.

_exact_rational(value) = Rational{BigInt}(value)
function _exact_rational(value::BigFloat)
    isfinite(value) || return Rational{BigInt}(value)
    # Read stored bits directly. Floating scaling can overflow under a compact
    # MPFR exponent range, even when both the input and rational are finite.
    # Integer scaling also preserves inputs wider than the caller's precision.
    significand, power, sign = Base.decompose(value)
    n = BigInt(sign) * significand
    return power >= 0 ? (n << power) // BigInt(1) : n // (BigInt(1) << -power)
end

_from_exact(::Type{T}, value::Rational{BigInt}) where {T<:AbstractFloat} = T(value)

function _from_exact(::Type{BigFloat}, value::Rational{BigInt})
    n, d = numerator(value), denominator(value)
    iszero(n) && return zero(BigFloat)
    nbits, dbits = ndigits(n; base = 2), ndigits(d; base = 2)
    minimum = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
    maximum = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
    if iszero(d) || (minimum <= nbits <= maximum && minimum <= dbits <= maximum)
        # Exact integer operands followed by one division retain the ordinary
        # result and the caller's precision and rounding mode.
        nfloat = BigFloat(n; precision = max(2, nbits))
        dfloat = BigFloat(d; precision = max(2, dbits))
        return nfloat / dfloat
    end

    negative = n < 0
    magnitude = abs(n)
    binade = nbits - dbits
    if binade >= 0 ? magnitude < d << binade : magnitude << -binade < d
        binade -= 1
    end
    # Retain p significant bits plus guard/sticky bits without first converting
    # either integer. A nonzero remainder distinguishes a tie from its upper
    # side; folding it into the low bit preserves all directed rounding choices.
    quantum = binade - precision(BigFloat) - 1
    scaled_n, scaled_d = quantum >= 0 ? (magnitude, d << quantum) : (magnitude << -quantum, d)
    whole, remainder = divrem(scaled_n, scaled_d)
    iszero(remainder) || (whole |= 1)
    negative && (whole = -whole)
    result = zero(BigFloat)
    # Significand and scale enter MPFR together, so compact exponent limits
    # apply to the complete result, including overflow and underflow endpoints.
    ccall(
        (:mpfr_set_z_2exp, Base.MPFR.libmpfr),
        Cint,
        (Ref{BigFloat}, Ref{BigInt}, Clong, Base.MPFR.MPFRRoundingMode),
        result,
        whole,
        quantum,
        Base.MPFR.rounding_raw(BigFloat),
    )
    return result
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

# Round a bounded integer ratio and binary scale without a floating intermediate.
function _from_exact(
    ::Type{T},
    signed_n::Int128,
    d::UInt64,
    exponent::Int,
) where {T<:Union{Float32,Float64}}
    iszero(d) && throw(ArgumentError("exact scaled conversion denominator must be positive"))
    iszero(signed_n) && return zero(T)
    negative = signed_n < 0
    n = negative ? -unsigned(signed_n) : unsigned(signed_n)
    denominator = UInt128(d)
    nw = 128 - leading_zeros(n)
    dw = 64 - leading_zeros(d)
    binade = nw - dw
    if binade >= 0 ? n < (denominator << binade) : (n << -binade) < denominator
        binade -= 1
    end
    actual = Int128(exponent) + binade
    minimum_quantum = Base.exponent(floatmin(T)) - precision(T) + 1
    actual > Base.exponent(floatmax(T)) && return negative ? -T(Inf) : T(Inf)
    actual < minimum_quantum - 1 && return negative ? -zero(T) : zero(T)
    quantum = max(Int(actual) - precision(T) + 1, minimum_quantum)
    shift = Int128(exponent) - quantum
    if shift >= 0
        # The output exponent guards imply nw+shift <= dw+precision(T) <= 117.
        numerator = n << Int(shift)
        scaled_denominator = denominator
    else
        denominator_bits = Int128(dw) - shift
        denominator_bits > 128 && return negative ? -zero(T) : zero(T)
        numerator = n
        scaled_denominator = denominator << Int(-shift)
    end
    whole, rest = divrem(numerator, scaled_denominator)
    complement = scaled_denominator - rest
    rounded = whole + (rest > complement || (rest == complement && isodd(whole)))
    magnitude = ldexp(T(rounded), quantum)
    return negative ? -magnitude : magnitude
end

# Round the square root directly from a nonnegative exact rational. Integer
# midpoint comparisons implement nearest, ties to even, independently of the
# caller's BigFloat rounding mode and without a transcendental intermediate.
function _sqrt_from_exact(::Type{T}, value::Rational{BigInt}) where {T<:AbstractFloat}
    n, d = numerator(value), denominator(value)
    n < 0 && throw(DomainError(value, "exact square root requires a nonnegative value"))
    iszero(n) && return zero(T)
    p = precision(T)
    maximum_exponent = exponent(floatmax(T))
    minimum_quantum = exponent(floatmin(T)) - p + 1
    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    if e >= 0 ? n < d << e : n << -e < d
        e -= 1
    end
    root_exponent = fld(e, 2)
    root_exponent > maximum_exponent && return T(Inf)

    # MPFR has no subnormals. Decide its half-minimum endpoint before rounding
    # a significand; a rounded intermediate can otherwise erase the tie side.
    minimum_positive = nextfloat(zero(T))
    minimum_exponent = exponent(minimum_positive)
    root_exponent < minimum_exponent - 1 && return zero(T)
    if root_exponent < minimum_exponent
        square_exponent = 2 * (minimum_exponent - 1)
        above_half = square_exponent >= 0 ? n > d << square_exponent : n << -square_exponent > d
        return above_half ? minimum_positive : zero(T)
    end

    quantum = max(root_exponent - p + 1, minimum_quantum)
    square_quantum = 2 * quantum
    scaled_n, scaled_d = square_quantum >= 0 ? (n, d << square_quantum) : (n << -square_quantum, d)
    whole = isqrt(div(scaled_n, scaled_d))
    midpoint_square = scaled_d * (2 * whole + 1)^2
    four_n = 4 * scaled_n
    if four_n > midpoint_square || (four_n == midpoint_square && isodd(whole))
        whole += 1
    end
    root_exponent == maximum_exponent && ndigits(whole; base = 2) > p && return T(Inf)
    if T === BigFloat
        # Construct significand and scale together: under a compact MPFR
        # exponent range the integer alone need not be representable.
        result = BigFloat(0; precision = p)
        ccall(
            (:mpfr_set_z_2exp, Base.MPFR.libmpfr),
            Cint,
            (Ref{BigFloat}, Ref{BigInt}, Clong, Base.MPFR.MPFRRoundingMode),
            result,
            whole,
            quantum,
            Base.MPFR.MPFRRoundNearest,
        )
        return result
    end
    return ldexp(T(whole), quantum)
end
