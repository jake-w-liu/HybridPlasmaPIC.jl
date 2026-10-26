# Exact conversions used by cold numerical fallbacks. Callers supply finite
# floating-point operands and finite Rational{BigInt} results.

_exact_rational(value) = Rational{BigInt}(value)

# Compare positive integer dyadics without shifting by their absolute scale.
# Equal effective bit lengths imply the alignment is bounded by operand widths.
function _compare_scaled(a::BigInt, ea::Int128, b::BigInt, eb::Int128)
    iszero(a) && return iszero(b) ? 0 : -1
    iszero(b) && return 1
    ba, bb = Int128(ndigits(a; base = 2)) + ea, Int128(ndigits(b; base = 2)) + eb
    ba != bb && return cmp(ba, bb)
    return ea >= eb ? cmp(a << Int(ea - eb), b) : cmp(a, b << Int(eb - ea))
end

# Compare one positive dyadic against a sum of two. Only overlapping binades
# need subtraction; a remote small term is compared to the exact remainder.
function _compare_scaled_sum(a::BigInt, ea::Int128, b::BigInt, eb::Int128, c::BigInt, ec::Int128)
    iszero(b) && return _compare_scaled(a, ea, c, ec)
    iszero(c) && return _compare_scaled(a, ea, b, eb)
    iszero(a) && return -1
    if _compare_scaled(b, eb, c, ec) < 0
        b, eb, c, ec = c, ec, b, eb
    end
    ba, bb = Int128(ndigits(a; base = 2)) + ea, Int128(ndigits(b; base = 2)) + eb
    ba < bb && return -1
    ba > bb + 1 && return 1
    scale = min(ea, eb)
    remainder = (a << Int(ea - scale)) - (b << Int(eb - scale))
    remainder <= 0 && return -1
    return _compare_scaled(remainder, scale, c, ec)
end

# Round a positive exact value whose binade and comparisons with q*2^e are
# supplied by its owner. The binary search uses p-bit integers, not absolute
# input exponents. This also handles MPFR's abrupt underflow (no subnormals).
function _round_exact_comparison(
    ::Type{T},
    binade::Int128,
    compare,
    negative::Bool,
) where {T<:AbstractFloat}
    mode = T === BigFloat ? rounding(BigFloat) : RoundNearest
    nearest = mode === RoundNearest
    away = mode === RoundFromZero || (negative ? mode === RoundDown : mode === RoundUp)
    p = precision(T)
    minimum = nextfloat(zero(T))
    minimum_exponent = exponent(minimum)
    maximum_exponent = exponent(floatmax(T))
    if binade < minimum_exponent
        above_half = nearest && compare(BigInt(1), Int128(minimum_exponent) - 1) > 0
        magnitude = away || above_half ? minimum : zero(T)
        return negative ? -magnitude : magnitude
    elseif binade > maximum_exponent
        magnitude = nearest || away ? T(Inf) : floatmax(T)
        return negative ? -magnitude : magnitude
    end
    quantum = binade - p + 1
    if T !== BigFloat
        quantum = max(quantum, Int128(minimum_exponent))
    end
    lower, upper = BigInt(0), BigInt(1) << p
    while upper - lower > 1
        middle = (lower + upper) >> 1
        if compare(middle, quantum) >= 0
            lower = middle
        else
            upper = middle
        end
    end
    if nearest
        side = compare(2lower + 1, quantum - 1)
        (side > 0 || (side == 0 && isodd(lower))) && (lower += 1)
    elseif away && compare(lower, quantum) > 0
        lower += 1
    end
    negative && (lower = -lower)
    if T === BigFloat
        result = zero(BigFloat)
        ccall(
            (:mpfr_set_z_2exp, Base.MPFR.libmpfr),
            Cint,
            (Ref{BigFloat}, Ref{BigInt}, Clong, Base.MPFR.MPFRRoundingMode),
            result,
            lower,
            Clong(quantum),
            Base.MPFR.rounding_raw(BigFloat),
        )
        return result
    end
    return ldexp(T(lower), Int(quantum))
end

function _from_exact(::Type{T}, value::Rational{BigInt}, power::Int128) where {T<:AbstractFloat}
    n, d = numerator(value), denominator(value)
    iszero(n) && return zero(T)
    negative = n < 0
    n = abs(n)
    binade = Int128(ndigits(n; base = 2) - ndigits(d; base = 2)) + power
    _compare_scaled(n, power, d, binade) < 0 && (binade -= 1)
    compare = (q, e) -> _compare_scaled(n, power, d * q, e)
    return _round_exact_comparison(T, binade, compare, negative)
end
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

# Existing collision callers require nearest, ties to even, independently of
# the caller's BigFloat mode. Scaled callers select rounding explicitly.
_sqrt_from_exact(::Type{T}, value::Rational{BigInt}) where {T<:AbstractFloat} =
    _sqrt_from_exact(T, value, Int128(0), RoundNearest, false)

# Round sign*sqrt(value*2^power) without expanding the absolute binary scale.
# Integer work depends on the rational's bit widths and destination precision.
function _sqrt_from_exact(
    ::Type{T},
    value::Rational{BigInt},
    power::Int128,
    mode::RoundingMode,
    negative::Bool,
) where {T<:AbstractFloat}
    n, d = numerator(value), denominator(value)
    n < 0 && throw(DomainError(value, "exact square root requires a nonnegative value"))
    iszero(d) && throw(ArgumentError("exact square root requires a finite value"))
    nearest = mode === RoundNearest
    away = mode === RoundFromZero || (negative ? mode === RoundDown : mode === RoundUp)
    (nearest || away || mode === RoundToZero || mode === RoundUp || mode === RoundDown) ||
        throw(ArgumentError("unsupported exact square-root rounding mode"))
    iszero(n) && return negative ? -zero(T) : zero(T)
    p = precision(T)
    maximum_exponent = exponent(floatmax(T))
    minimum_quantum = exponent(floatmin(T)) - p + 1

    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    if e >= 0 ? n < d << e : n << -e < d
        e -= 1
    end
    # Int128 holds sums of the bounded number of native MPFR input exponents,
    # including cancellation near their Int64 limits.
    root_exponent = fld(Int128(e) + power, 2)
    if root_exponent > maximum_exponent
        magnitude = nearest || away ? T(Inf) : floatmax(T)
        return negative ? -magnitude : magnitude
    end

    minimum_positive = nextfloat(zero(T))
    minimum_exponent = exponent(minimum_positive)
    if root_exponent < minimum_exponent
        above_half = false
        if nearest && root_exponent == minimum_exponent - 1
            # Resolve the half-minimum tie before significand rounding. MPFR
            # has no subnormals; rounding first could erase the side of a tie.
            shift = power - 2 * Int128(minimum_exponent - 1)
            above_half = shift >= 0 ? n << Int(shift) > d : n > d << Int(-shift)
        end
        magnitude = away || above_half ? minimum_positive : zero(T)
        return negative ? -magnitude : magnitude
    end

    quantum = max(root_exponent - p + 1, Int128(minimum_quantum))
    shift = power - 2 * quantum
    # The selected root binade cancels power. These shifts are bounded by the
    # rational bit widths plus twice the output precision, even for huge power.
    scaled_n, scaled_d = shift >= 0 ? (n << Int(shift), d) : (n, d << Int(-shift))
    whole = isqrt(div(scaled_n, scaled_d))
    if nearest
        midpoint_square = scaled_d * (2 * whole + 1)^2
        four_n = 4 * scaled_n
        if four_n > midpoint_square || (four_n == midpoint_square && isodd(whole))
            whole += 1
        end
    elseif away && whole^2 * scaled_d != scaled_n
        whole += 1
    end
    if root_exponent == maximum_exponent && ndigits(whole; base = 2) > p
        return negative ? -T(Inf) : T(Inf)
    end
    if T === BigFloat
        # The integer alone may lie outside the caller's exponent range. This
        # signed significand and quantum already lie exactly on the output grid.
        result = BigFloat(0; precision = p)
        ccall(
            (:mpfr_set_z_2exp, Base.MPFR.libmpfr),
            Cint,
            (Ref{BigFloat}, Ref{BigInt}, Clong, Base.MPFR.MPFRRoundingMode),
            result,
            negative ? -whole : whole,
            Clong(quantum),
            Base.MPFR.MPFRRoundNearest,
        )
        return result
    end
    magnitude = ldexp(T(whole), Int(quantum))
    return negative ? -magnitude : magnitude
end

# Exact totals use bounded integer accumulation when the represented exponents
# permit it; the wider path retains every finite stored operand.
function _exact_float_sum(values::AbstractArray{T}) where {T<:Union{Float16,Float32,Float64}}
    Q = Rational{BigInt}
    isempty(values) && return Q(0)
    reference = first(values)
    isfinite(reference) || throw(ArgumentError("exact dyadic sum requires finite values"))
    minimum_exponent = iszero(reference) ? typemax(Int) : last(frexp(reference))
    maximum_exponent = iszero(reference) ? typemin(Int) : minimum_exponent
    all_equal = true
    for value in values
        if value != reference
            all_equal = false
            isfinite(value) || throw(ArgumentError("exact dyadic sum requires finite values"))
            if !iszero(value)
                _, e = frexp(value)
                minimum_exponent = min(minimum_exponent, e)
                maximum_exponent = max(maximum_exponent, e)
            end
        end
    end
    if all_equal
        return length(values) == 1 ? Q(reference) : Q(reference) * length(values)
    end
    minimum_exponent == typemax(Int) && return Q(0)
    p = precision(T)
    # n < 2^ndigits(n), so this bounds every absolute partial sum strictly
    # below 2^127 in the selected quantum, including arbitrary sign order.
    if p + maximum_exponent - minimum_exponent + ndigits(length(values); base = 2) > 127
        return _exact_repeated_float_sum(values)
    end
    total = zero(Int128)
    for value in values
        if !iszero(value)
            fraction, e = frexp(value)
            significand = Int128(ldexp(fraction, p))
            total += significand << (e - minimum_exponent)
        end
    end
    quantum = minimum_exponent - p
    integer = BigInt(total)
    return quantum >= 0 ? Q(integer << quantum) : integer // (BigInt(1) << -quantum)
end

function _exact_repeated_float_sum(values)
    Q = Rational{BigInt}
    total = Q(0)
    previous = zero(eltype(values))
    count = 0
    for value in values
        if count > 0 && value != previous
            total += count == 1 ? Q(previous) : Q(previous) * count
            count = 0
        end
        previous = value
        count += 1
    end
    count > 0 && (total += count == 1 ? Q(previous) : Q(previous) * count)
    return total
end
