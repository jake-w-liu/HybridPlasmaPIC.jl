# Shared complete-input arithmetic for polytropic pressure.
# Ordinary powers retain native evaluation with a bounded ratio correction;
# unresolved range/correction cases use complete scaled intervals and tie checks.

function _polytropic_positive_exp_bounds(x::Rational{BigInt}, bits::Int)
    zero(Rational{BigInt}) <= x <= one(Rational{BigInt}) ||
        throw(ArgumentError("exponent must lie in [0, 1]"))
    bits > 0 || throw(ArgumentError("bound precision must be positive"))
    iszero(x) && return one(Rational{BigInt}), one(Rational{BigInt})
    tolerance = BigInt(1) // (BigInt(1) << bits)
    partial = term = one(Rational{BigInt})
    k = 0
    while true
        following = term * x / (k + 1)
        # Every later term ratio is <= x/(k+2), so this bounds the full tail.
        tail = following / (1 - x / (k + 2))
        tail <= tolerance && return partial, partial + tail
        partial += following
        term = following
        k = Base.checked_add(k, 1)
    end
end

function _polytropic_exp_bounds(x::Rational{BigInt}, bits::Int)
    -one(Rational{BigInt}) <= x <= one(Rational{BigInt}) ||
        throw(ArgumentError("reduced exponent must lie in [-1, 1]"))
    if x < 0
        lower, upper = _polytropic_positive_exp_bounds(-x, bits)
        return inv(upper), inv(lower)
    end
    return _polytropic_positive_exp_bounds(x, bits)
end

function _polytropic_scaled_exp_bounds(lower::Rational{BigInt}, upper::Rational{BigInt}, bits::Int)
    lower <= upper || throw(ArgumentError("exponent interval must be ordered"))
    upper - lower <= Rational{BigInt}(1, 8) || throw(ArgumentError("exponent interval is too wide"))
    bits > 0 || throw(ArgumentError("bound precision must be positive"))
    midpoint = (lower + upper) / 2
    magnitude_bits =
        iszero(midpoint) ? 0 :
        max(
            0,
            ndigits(abs(numerator(midpoint)); base = 2) - ndigits(denominator(midpoint); base = 2) +
            1,
        )
    logbits = Base.checked_add(Base.checked_add(bits, magnitude_bits), 8)
    loglo, loghi = _isothermal_atanh_log_bounds(Rational{BigInt}(1, 3), logbits)
    exponent = floor(BigInt, 2midpoint / (loglo + loghi))
    reduced_lower = lower - exponent * (exponent >= 0 ? loghi : loglo)
    reduced_upper = upper - exponent * (exponent >= 0 ? loglo : loghi)
    # An unresolved reduction asks the caller for tighter logarithm bounds.
    -one(Rational{BigInt}) <= reduced_lower <= reduced_upper <= one(Rational{BigInt}) ||
        return nothing
    lo, _ = _polytropic_exp_bounds(reduced_lower, bits)
    _, hi = _polytropic_exp_bounds(reduced_upper, bits)
    return lo, hi, exponent
end

"""
Return the exact nonnegative q-th root of `base`, or `nothing` when it is
irrational, for a positive reduced binary-rational exponent p/q. The returned
root has not been raised to p, so even very large stored exponents stay bounded.
"""
function _polytropic_rational_root(base::Rational{BigInt}, gamma::Rational{BigInt})
    base >= 0 || throw(ArgumentError("power base must be nonnegative"))
    gamma > 0 || throw(ArgumentError("power exponent must be positive"))
    q = denominator(gamma)
    ispow2(q) || throw(ArgumentError("power exponent denominator must be a power of two"))
    iszero(base) && return zero(Rational{BigInt})
    base == 1 && return one(Rational{BigInt})
    a, b = numerator(base), denominator(base)
    remaining = ndigits(q; base = 2) - 1
    while remaining > 0
        ar, br = isqrt(a), isqrt(b)
        ar * ar == a && br * br == b || return nothing
        a, b = ar, br
        remaining -= 1
    end
    return a // b
end

function _polytropic_integer_trim(lower::BigInt, upper::BigInt, exponent::BigInt, bits::Int)
    shift = max(0, ndigits(upper; base = 2) - bits)
    iszero(shift) && return lower, upper, exponent
    lo = lower >> shift
    hi = upper >> shift
    hi += (hi << shift) != upper
    return lo, hi, exponent + shift
end

function _polytropic_integer_product(a, b, bits::Int)
    return _polytropic_integer_trim(a[1] * b[1], a[2] * b[2], a[3] + b[3], bits)
end

function _polytropic_integer_initial(base::Rational{BigInt}, bits::Int)
    n, d = numerator(base), denominator(base)
    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    if e >= 0 ? n < d << e : n << -e < d
        e -= 1
    end
    scale = e - bits + 1
    a, b = scale >= 0 ? (n, d << scale) : (n << -scale, d)
    lower, remainder = divrem(a, b)
    return lower, lower + !iszero(remainder), BigInt(scale)
end

"""Enclose a nonnegative rational raised to a nonnegative integer, as lo/hi * 2^e."""
function _polytropic_integer_bounds(base::Rational{BigInt}, exponent::BigInt, bits::Int)
    base >= 0 || throw(ArgumentError("power base must be nonnegative"))
    exponent >= 0 || throw(ArgumentError("integer power must be nonnegative"))
    bits > 0 || throw(ArgumentError("power precision must be positive"))
    iszero(exponent) && return one(Rational{BigInt}), one(Rational{BigInt}), BigInt(0)
    iszero(base) && return zero(Rational{BigInt}), zero(Rational{BigInt}), BigInt(0)
    base == 1 && return one(Rational{BigInt}), one(Rational{BigInt}), BigInt(0)
    # Initial relative uncertainty can be amplified by the integer exponent.
    # Its bit length supplies that precision; extra bits cover subsequent
    # outward products. Correct enclosure does not depend on a width estimate.
    working_bits = Base.checked_add(bits, Base.checked_add(ndigits(exponent; base = 2), 4))
    factor = _polytropic_integer_initial(base, working_bits)
    result = (BigInt(1), BigInt(1), BigInt(0))
    remaining = exponent
    while !iszero(remaining)
        isodd(remaining) && (result = _polytropic_integer_product(result, factor, working_bits))
        remaining >>= 1
        iszero(remaining) || (factor = _polytropic_integer_product(factor, factor, working_bits))
    end
    lower, upper, e = result
    shift = ndigits(upper; base = 2) - 1
    denominator = BigInt(1) << shift
    return lower // denominator, upper // denominator, e + shift
end

function _polytropic_magnitude_bits(x::Rational{BigInt})
    iszero(x) && return 0
    return max(0, ndigits(abs(numerator(x)); base = 2) - ndigits(denominator(x); base = 2) + 1)
end

function _polytropic_logarithm_bounds(x::Rational{BigInt}, bits::Int)
    x > 0 || throw(ArgumentError("logarithm input must be positive"))
    logtwo = _isothermal_atanh_log_bounds(Rational{BigInt}(1, 3), bits)
    nlo, nhi = _isothermal_integer_log_bounds(numerator(x), bits, logtwo)
    dlo, dhi = _isothermal_integer_log_bounds(denominator(x), bits, logtwo)
    return nlo - dhi, nhi - dlo
end

function _polytropic_power_bounds(base::Rational{BigInt}, gamma::Rational{BigInt}, bits::Int)
    base >= 0 || throw(ArgumentError("power base must be nonnegative"))
    gamma > 0 || throw(ArgumentError("power exponent must be positive"))
    bits > 0 || throw(ArgumentError("bound precision must be positive"))
    iszero(base) && return zero(Rational{BigInt}), zero(Rational{BigInt}), BigInt(0)
    base == 1 && return one(Rational{BigInt}), one(Rational{BigInt}), BigInt(0)
    if ispow2(denominator(gamma))
        root = _polytropic_rational_root(base, gamma)
        root === nothing || return _polytropic_integer_bounds(root, numerator(gamma), bits)
    end
    size = max(ndigits(numerator(base); base = 2), ndigits(denominator(base); base = 2))
    size_bits = ndigits(size; base = 2)
    logbits =
        Base.checked_add(Base.checked_add(bits, _polytropic_magnitude_bits(gamma)), size_bits + 8)
    while true
        lower, upper = _polytropic_logarithm_bounds(base, logbits)
        lower *= gamma
        upper *= gamma
        if upper - lower <= Rational{BigInt}(1, 8)
            interval = _polytropic_scaled_exp_bounds(lower, upper, bits)
            interval === nothing || return interval
        end
        logbits = Base.checked_mul(logbits, 2)
    end
end

function _polytropic_quantize(value::Rational{BigInt}, shift::BigInt, left_limit::Int)
    value >= 0 || throw(ArgumentError("scaled interval endpoints must be nonnegative"))
    shift <= left_limit || throw(ArgumentError("left shift exceeds the precision bound"))
    iszero(value) && return BigInt(0), BigInt(0)
    n, d = numerator(value), denominator(value)
    # The positive value is strictly between zero and one; avoid a huge shift.
    shift <= -BigInt(ndigits(n; base = 2)) && return BigInt(0), BigInt(1)
    count = Int(shift)
    a, b = count >= 0 ? (n << count, d) : (n, d << -count)
    whole, remainder = divrem(a, b)
    return whole, whole + !iszero(remainder)
end

function _polytropic_scaled_sum_bounds(intervals, bits::Int)
    bits > 0 || throw(ArgumentError("sum precision must be positive"))
    largest = nothing
    for (lower, upper, exponent) in intervals
        0 <= lower <= upper ||
            throw(ArgumentError("scaled intervals must be nonnegative and ordered"))
        upper > 0 || continue
        largest = largest === nothing ? exponent : max(largest, exponent)
    end
    largest === nothing && return BigInt(0), BigInt(0), BigInt(0)
    scale = BigInt(largest) - bits
    lower = upper = BigInt(0)
    for (lo, hi, exponent) in intervals
        iszero(hi) && continue
        shift = BigInt(exponent) - scale
        a, _ = _polytropic_quantize(lo, shift, bits)
        _, b = _polytropic_quantize(hi, shift, bits)
        lower += a
        upper += b
    end
    return lower, upper, scale
end

function _polytropic_round_scaled(
    ::Type{T},
    coefficient::BigInt,
    exponent::BigInt,
) where {T<:Union{Float32,Float64}}
    iszero(coefficient) && return zero(T)
    negative = coefficient < 0
    n = abs(coefficient)
    actual_exponent = exponent + ndigits(n; base = 2) - 1
    minimum_quantum = Base.exponent(floatmin(T)) - precision(T) + 1
    actual_exponent > Base.exponent(floatmax(T)) && return negative ? -T(Inf) : T(Inf)
    actual_exponent < minimum_quantum - 1 && return negative ? -zero(T) : zero(T)
    # In the remaining range, the shift costs at most the existing coefficient's
    # bit length plus the target's finite exponent range.
    e = Int(exponent)
    value = e >= 0 ? Rational{BigInt}(coefficient << e) : coefficient // (BigInt(1) << -e)
    return _from_exact(T, value)
end

function _polytropic_normalized_factor(value::Rational{BigInt})
    iszero(value) && return false, zero(Rational{BigInt}), BigInt(0)
    negative = value < 0
    n, d = abs(numerator(value)), denominator(value)
    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    if e >= 0 ? n < d << e : n << -e < d
        e -= 1
    end
    mantissa = e >= 0 ? n // (d << e) : (n << -e) // d
    return negative, mantissa, BigInt(e)
end

function _polytropic_integer_power_equal(base::BigInt, power::BigInt, target::BigInt)
    base >= 1 && power >= 1 && target >= 1 ||
        throw(ArgumentError("integer power equality requires positive integers"))
    base == 1 && return target == 1
    target_bits = ndigits(target; base = 2)
    base_bits = ndigits(base; base = 2)
    power * (base_bits - 1) >= target_bits && return false
    # This branch proves power < target_bits and bounds the constructed
    # power to fewer than 2*target_bits bits, independently of the input power.
    return base^Int(power) == target
end

function _polytropic_power_equal(
    base::Rational{BigInt},
    gamma::Rational{BigInt},
    target::Rational{BigInt},
)
    base >= 0 && gamma > 0 && ispow2(denominator(gamma)) || throw(
        ArgumentError("power equality requires a nonnegative base and positive binary exponent"),
    )
    target < 0 && return false
    iszero(base) && return iszero(target)
    iszero(target) && return false
    root = _polytropic_rational_root(base, gamma)
    root === nothing && return false
    power = numerator(gamma)
    return _polytropic_integer_power_equal(numerator(root), power, numerator(target)) &&
           _polytropic_integer_power_equal(denominator(root), power, denominator(target))
end

function _polytropic_midpoint_after(value::T) where {T<:Union{Float32,Float64}}
    adjacent = nextfloat(value)
    return isfinite(adjacent) ? (Rational{BigInt}(value) + Rational{BigInt}(adjacent)) / 2 :
           Rational{BigInt}(value) +
           (Rational{BigInt}(value) - Rational{BigInt}(prevfloat(value))) / 2
end

function _polytropic_complete_power(
    ::Type{T},
    base::Rational{BigInt},
    gamma::Rational{BigInt},
    coefficient::Rational{BigInt},
) where {T<:Union{Float32,Float64}}
    base >= 0 && coefficient >= 0 && gamma > 0 && ispow2(denominator(gamma)) || throw(
        ArgumentError(
            "complete pressure requires nonnegative factors and a positive binary exponent",
        ),
    )
    (iszero(base) || iszero(coefficient)) && return zero(T)
    _, multiplier, exponent = _polytropic_normalized_factor(coefficient)
    bits = precision(T) + 16
    while true
        lower, upper, e = _polytropic_power_bounds(base, gamma, bits)
        lo, hi, scale = _polytropic_scaled_sum_bounds(
            ((lower * multiplier, upper * multiplier, e + exponent),),
            bits,
        )
        left = _polytropic_round_scaled(T, lo, scale)
        right = _polytropic_round_scaled(T, hi, scale)
        isequal(left, right) && return left
        if isfinite(left) && nextfloat(left) == right
            boundary = _polytropic_midpoint_after(left)
            _polytropic_power_equal(base, gamma, boundary / coefficient) &&
                return _from_exact(T, boundary)
        end
        bits = Base.checked_mul(bits, 2)
    end
end

@inline _polytropic_normal_positive(x::Float64) = isfinite(x) && x >= 2floatmin(Float64)

function _polytropic_pressure_probe_power(
    ::Type{S},
    density::T,
    closure::PolytropicElectrons{C},
    ::Type{P},
) where {
    S<:Union{Float32,Float64},
    T<:Union{Float32,Float64},
    C<:Union{Float32,Float64},
    P<:Union{Float32,Float64},
}
    value, reference = Float64(density), Float64(closure.n0)
    gamma, coefficient = Float64(closure.γ), Float64(closure.pe0)
    ratio = Float64(P(density) / P(closure.n0))
    isfinite(ratio) && ratio >= 2Float64(floatmin(P)) || return nothing
    a, ae = frexp(value)
    b, be = frexp(reference)
    normalized_ratio = ldexp(ratio, be - ae)
    residual = fma(-normalized_ratio, b, a) / b
    relative = residual / normalized_ratio
    correction = gamma * relative
    # For |gamma*e| <= sqrt(eps)/4 and |e| <= eps, the logarithm and
    # exponential remainder bounds make the omitted nonlinear term < eps/16.
    # An unresolved probe advances to wider arithmetic or complete evaluation.
    abs(correction) <= Float64(sqrt(eps(P)) / 4) || return nothing
    power = Float64(P(ratio)^P(gamma))
    isfinite(power) && power >= 2Float64(floatmin(P)) || return nothing
    adjusted = fma(power, correction, power)
    _polytropic_normal_positive(adjusted) || return nothing
    product = coefficient * adjusted
    _polytropic_normal_positive(product) || return nothing
    result = S(product)
    2floatmin(S) <= result <= floatmax(S) / 2 || return nothing
    return result
end
function _polytropic_pressure_probe(
    ::Type{S},
    density::T,
    closure::PolytropicElectrons{C},
) where {S<:Union{Float32,Float64},T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    if S === Float32 && T === Float32 && C === Float32
        result = _polytropic_pressure_probe_power(S, density, closure, Float32)
        result === nothing || return result
    end
    return _polytropic_pressure_probe_power(S, density, closure, Float64)
end

function _polytropic_pressure(
    ::Type{S},
    density::T,
    closure::PolytropicElectrons{C},
) where {S<:Union{Float32,Float64},T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    if !(
        isfinite(density) &&
        isfinite(closure.pe0) &&
        isfinite(closure.n0) &&
        isfinite(closure.γ) &&
        closure.pe0 >= 0 &&
        closure.n0 > 0 &&
        closure.γ > 0
    ) || density <= 0
        # Raw shock callers retain signed powers and their existing failures.
        # The closure caller applies its density floor before entering here.
        return S(closure.pe0 * (density / closure.n0)^closure.γ)
    end
    iszero(closure.pe0) && return S(closure.pe0)
    value = _polytropic_pressure_probe(S, density, closure)
    value === nothing || return value
    return _polytropic_complete_power(
        S,
        Rational{BigInt}(density) / Rational{BigInt}(closure.n0),
        Rational{BigInt}(closure.γ),
        Rational{BigInt}(closure.pe0),
    )
end


# Only known CPU storage reaches the complete scalar evaluator. Built-in
# wrappers delegate ownership to their parent; unknown storage keeps its backend.
_polytropic_host_storage(::Any) = false
_polytropic_host_storage(::Array) = true
_polytropic_host_storage(
    a::Union{SubArray,Base.ReshapedArray,Base.ReinterpretArray,PermutedDimsArray,Transpose,Adjoint},
) = _polytropic_host_storage(parent(a))

function _polytropic_pressure_array!(
    pe::Array{S},
    n::Array{T},
    c::PolytropicElectrons{C},
    ::Val{Floor},
) where {S<:Union{Float32,Float64},T<:Union{Float32,Float64},C<:Union{Float32,Float64},Floor}
    axes(pe) == axes(n) ||
        throw(DimensionMismatch("repeated-density traversal requires equal axes"))
    pe === n ||
        !Base.mightalias(pe, n) ||
        throw(
            ArgumentError(
                "repeated-density traversal requires independent arrays or exact self-alias",
            ),
        )
    previous_density = zero(T)
    previous_pressure = zero(S)
    have_previous = false
    @inbounds for i in eachindex(pe, n)
        density = Floor ? max(n[i], zero(T)) : n[i]
        if have_previous && isequal(density, previous_density)
            pe[i] = previous_pressure
        else
            value = _polytropic_pressure(S, density, c)
            pe[i] = value
            previous_density = density
            previous_pressure = value
            have_previous = isfinite(density)
        end
    end
    return pe
end

_polytropic_array_compatible(pe, n) =
    pe isa Array && n isa Array && axes(pe) == axes(n) && (pe === n || !Base.mightalias(pe, n))

function electron_pressure!(
    pe::AbstractArray{S},
    n::AbstractArray{T},
    closure::PolytropicElectrons{C},
) where {S<:Union{Float32,Float64},T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    # Union element storage belongs to the generic broadcast owner; it cannot
    # supply the single concrete output/input precision used by this path.
    (S === Float32 || S === Float64) && (T === Float32 || T === Float64) ||
        return _legacy_polytropic_pressure!(pe, n, closure)
    if _polytropic_host_storage(pe) && _polytropic_host_storage(n)
        if pe isa Array && n isa Array && _polytropic_array_compatible(pe, n)
            return _polytropic_pressure_array!(pe, n, closure, Val(true))
        end
        broadcast!(pe, n) do density
            _polytropic_pressure(S, max(density, zero(T)), closure)
        end
        return pe
    end
    return _legacy_polytropic_pressure!(pe, n, closure)
end

function _polytropic_unreferenced_pressure(
    ::Type{T},
    density::T,
    coefficient::T,
    gamma::T,
) where {T}
    if T === Float32 || T === Float64
        closure = PolytropicElectrons{T}(coefficient, one(T), gamma)
        return _polytropic_pressure(T, density, closure)
    end
    # Retain stored-input arithmetic for other precisions without dividing by one.
    return coefficient * density^gamma
end

function _polytropic_unfloored_pressure!(
    pe::Array{T},
    n::Array{T},
    coefficient::T,
    gamma::T,
) where {T}
    if T === Float32 || T === Float64
        closure = PolytropicElectrons{T}(coefficient, one(T), gamma)
        if _polytropic_array_compatible(pe, n)
            return _polytropic_pressure_array!(pe, n, closure, Val(false))
        end
        # Preserve the singleton destination type in scalar broadcast storage.
        broadcast!(_polytropic_pressure, pe, Ref{Type{T}}(T), n, Ref(closure))
    else
        @. pe = coefficient * n^gamma
    end
    return pe
end
