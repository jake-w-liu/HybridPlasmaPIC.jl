# Complete and bounded reductions for the isothermal free-energy diagnostic.


function _isothermal_add_log!(
    terms::Dict{UInt64,Rational{BigInt}},
    base::BigInt,
    coefficient::Rational{BigInt},
)
    base >= 1 || throw(ArgumentError("logarithm base must be a positive integer"))
    (base == 1 || iszero(coefficient)) && return terms
    # A Float32/64 density's odd significand has at most 53 bits.
    key = UInt64(base)
    value = get(terms, key, zero(Rational{BigInt})) + coefficient
    if iszero(value)
        delete!(terms, key)
    else
        terms[key] = value
    end
    return terms
end

function _isothermal_atanh_log_bounds(z::Rational{BigInt}, bits::Int)
    zero(Rational{BigInt}) <= z <= Rational{BigInt}(1, 3) ||
        throw(ArgumentError("reduced logarithm argument is outside [0, 1/3]"))
    bits > 0 || throw(ArgumentError("bound precision must be positive"))
    iszero(z) && return zero(Rational{BigInt}), zero(Rational{BigInt})
    tolerance = BigInt(1) // (BigInt(1) << bits)
    square = z * z
    power = z
    partial = zero(Rational{BigInt})
    denominator = BigInt(1)
    while true
        partial += power / denominator
        power *= square
        denominator += 2
        tail = 2power / (denominator * (1 - square))
        tail <= tolerance && return 2partial, 2partial + tail
    end
end

function _isothermal_integer_log_bounds(base::BigInt, bits::Int, logtwo)
    base >= 1 || throw(ArgumentError("logarithm base must be a positive integer"))
    exponent = ndigits(base; base = 2) - 1
    scale = BigInt(1) << exponent
    lo, hi = _isothermal_atanh_log_bounds((base - scale) // (base + scale), bits)
    return lo + exponent * logtwo[1], hi + exponent * logtwo[2]
end

function _isothermal_finish_log_energy(
    ::Type{T},
    terms::Dict{UInt64,Rational{BigInt}},
) where {T<:Union{Float32,Float64}}
    isempty(terms) && return zero(T)
    # All odd-base coefficients share the sign of Te * volume. Clearing their
    # rational denominators gives an odd integer product on one side of any
    # putative zero identity and a power of two on the other. A nonempty odd
    # part therefore cannot cancel; a pure log(2) part is zero only if empty.
    odd_sign = 0
    for (base, coefficient) in terms
        base > 1 && !iszero(coefficient) ||
            throw(ArgumentError("logarithm terms must have nonunit bases and nonzero coefficients"))
        base == 2 && continue
        isodd(base) || throw(ArgumentError("density logarithm bases must be odd or two"))
        current_sign = coefficient > 0 ? 1 : -1
        odd_sign == 0 ||
            odd_sign == current_sign ||
            throw(ArgumentError("odd density logarithms must share a coefficient sign"))
        odd_sign = current_sign
    end
    bits = precision(T) + 16
    while true
        logtwo = _isothermal_atanh_log_bounds(Rational{BigInt}(1, 3), bits)
        lower = upper = zero(Rational{BigInt})
        for (base, coefficient) in terms
            lo, hi = _isothermal_integer_log_bounds(BigInt(base), bits, logtwo)
            lower += coefficient * (coefficient > 0 ? lo : hi)
            upper += coefficient * (coefficient > 0 ? hi : lo)
        end
        lo, hi = _from_exact(T, lower), _from_exact(T, upper)
        isequal(lo, hi) && return lo
        bits = Base.checked_mul(bits, 2)
    end
end

function _isothermal_exact_energy(
    n::AbstractArray{T},
    temperature::C,
    dx,
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    all(isfinite, n) && isfinite(temperature) && all(isfinite, dx) ||
        throw(ArgumentError("the complete isothermal evaluator requires finite inputs"))
    multiplier = Rational{BigInt}(temperature) * prod(Rational{BigInt}, dx)
    iszero(multiplier) && return zero(T)
    terms = Dict{UInt64,Rational{BigInt}}()
    for value in n
        value > zero(T) || continue
        density = Rational{BigInt}(value)
        coefficient = density * multiplier
        odd = numerator(density)
        exponent = -(ndigits(denominator(density); base = 2) - 1)
        while iseven(odd)
            odd >>= 1
            exponent += 1
        end
        _isothermal_add_log!(terms, BigInt(2), coefficient * exponent)
        _isothermal_add_log!(terms, odd, coefficient)
    end
    return _isothermal_finish_log_energy(T, terms)
end




_isothermal_ceil_log2(x::Integer) = 128 - leading_zeros(Int128(x) - 1)

function _isothermal_statistics(n::AbstractArray{T}) where {T<:Union{Float32,Float64}}
    p = precision(T)
    scale = Int128(1) << p
    scale_float = ldexp(one(T), p)
    delta = maximum_delta = Int128(0)
    active = 0
    for value in n
        isfinite(value) || return nothing
        value > zero(T) || continue
        abs(value - one(T)) <= ldexp(one(T), -8) || return nothing
        d = Int128(Int64(value * scale_float)) - scale
        delta += d
        maximum_delta = max(maximum_delta, abs(d))
        active += !iszero(d)
    end
    q = iszero(maximum_delta) ? 0 : _isothermal_ceil_log2(maximum_delta)
    cap = iszero(active) ? 124 : min(124 - q, 124 - _isothermal_ceil_log2(active))
    return (; p, scale, delta, maximum_delta, active, q, cap)
end

function _isothermal_coefficient_bounds(k::Int, exponent::Int, bits::Int)
    power = bits + exponent * (k - 2)
    if power < 0
        lo, hi = Int128(0), Int128(1)
    else
        numerator = Int128(1) << power
        denominator = k * (k - 1)
        lo, remainder = divrem(numerator, denominator)
        hi = lo + !iszero(remainder)
    end
    return iseven(k) ? (lo, hi) : (-hi, -lo)
end

@inline function _isothermal_multiply_y(
    lo::I,
    hi::I,
    d::I,
    q::Int,
    mask::I,
) where {I<:Union{Int64,Int128}}
    a, b = d >= 0 ? (lo * d, hi * d) : (hi * d, lo * d)
    lower = a >> q
    upper = (b >> q) + !iszero(b & mask)
    return lower, upper
end

function _isothermal_scaled_polynomial_storage(
    n::AbstractArray{T},
    state,
    bits::Int,
    ::Val{K},
    ::Type{I},
) where {T,K,I<:Union{Int64,Int128}}
    exponent = state.q - state.p
    coefficients = ntuple(Val(K - 1)) do i
        lo, hi = _isothermal_coefficient_bounds(i + 1, exponent, bits)
        return I(lo), I(hi)
    end
    mask = (I(1) << state.q) - 1
    scale = I(state.scale)
    scale_float = T(state.scale)
    lower = upper = Int128(0)
    for value in n
        value > zero(T) || continue
        d = I(Int64(value * scale_float)) - scale
        iszero(d) && continue
        lo, hi = coefficients[K-1]
        for k = (K-1):-1:2
            lo, hi = _isothermal_multiply_y(lo, hi, d, state.q, mask)
            c = coefficients[k-1]
            lo += c[1]
            hi += c[2]
        end
        lo, hi = _isothermal_multiply_y(lo, hi, d, state.q, mask)
        lo, hi = _isothermal_multiply_y(lo, hi, d, state.q, mask)
        lower += lo
        upper += hi
    end
    return lower, upper
end

function _isothermal_scaled_polynomial(n, state, bits, degree)
    # The coefficient absolute sum is <1; three sign/carry bits cover each
    # product and outward-rounding increment. Whole-grid totals remain Int128.
    if bits + state.q <= 60
        return _isothermal_scaled_polynomial_storage(n, state, bits, degree, Int64)
    end
    return _isothermal_scaled_polynomial_storage(n, state, bits, degree, Int128)
end

@inline _isothermal_bitwidth(x::Int128) = 128 - leading_zeros(_exact_int128_magnitude(x))

# Return floor/ceil(x / 2^shift). Callers bound left shifts before reaching here.
@inline function _isothermal_outward_shift(x::Int128, shift::Int)
    shift <= 0 && return (x << -shift, x << -shift)
    shift >= 128 && return (Int128(x < 0 ? -1 : 0), Int128(x > 0 ? 1 : 0))
    mask = (UInt128(1) << shift) - 1
    lo = x >> shift
    return lo, lo + !iszero(unsigned(x) & mask)
end

function _isothermal_add_dyadic(a, b)
    iszero(a[1]) && iszero(a[2]) && return b
    iszero(b[1]) && iszero(b[2]) && return a
    abits = max(_isothermal_bitwidth(a[1]), _isothermal_bitwidth(a[2]))
    bbits = max(_isothermal_bitwidth(b[1]), _isothermal_bitwidth(b[2]))
    # Both aligned magnitudes fit within 122 bits before outward rounding;
    # adding the rounded endpoints therefore remains below the signed limit.
    exponent = max(min(a[3], b[3]), max(a[3] + abits, b[3] + bbits) - 122)
    alo, _ = _isothermal_outward_shift(a[1], exponent - a[3])
    _, ahi = _isothermal_outward_shift(a[2], exponent - a[3])
    blo, _ = _isothermal_outward_shift(b[1], exponent - b[3])
    _, bhi = _isothermal_outward_shift(b[2], exponent - b[3])
    return alo + blo, ahi + bhi, exponent
end

function _isothermal_multiply_dyadic(a, factor::T) where {T<:Union{Float32,Float64}}
    iszero(factor) && return (Int128(0), Int128(0), 0)
    fraction, exponent = frexp(factor)
    integer = Int128(ldexp(fraction, precision(T)))
    exponent -= precision(T)
    zeros = trailing_zeros(integer)
    integer >>= zeros
    exponent += zeros
    # An outward ceil can add one bit. The product then needs at most 123 bits.
    shift = max(
        0,
        max(_isothermal_bitwidth(a[1]), _isothermal_bitwidth(a[2])) +
        _isothermal_bitwidth(integer) - 122,
    )
    lo, _ = _isothermal_outward_shift(a[1], shift)
    _, hi = _isothermal_outward_shift(a[2], shift)
    lower, upper = integer > 0 ? (lo * integer, hi * integer) : (hi * integer, lo * integer)
    return lower, upper, a[3] + shift + exponent
end

function _isothermal_tail_units(active::Int, r::Int, bits::Int, ::Val{K}) where {K}
    # ceil(active * 2^(bits-r*(K-2)) / (K*(K+1)*(2^r-1))).
    iszero(active) && return Int128(0)
    numerator = Int128(active)
    denominator = Int128(K * (K + 1)) * ((Int128(1) << r) - 1)
    power = bits - r * (K - 2)
    if power >= 0
        numerator <<= power
    else
        # A denominator this wide exceeds active (an Int), so the positive ceil is one.
        _isothermal_bitwidth(denominator) - power > 126 && return Int128(1)
        denominator <<= -power
    end
    whole, rest = divrem(numerator, denominator)
    return whole + !iszero(rest)
end

function _isothermal_complete_dyadic(
    n::AbstractArray{T},
    temperature,
    dx,
    state,
    bits,
    degree,
) where {T}
    @assert 40 <= bits <= state.cap
    lower, upper = _isothermal_scaled_polynomial(n, state, bits, degree)
    r = state.p - state.q
    tail = _isothermal_tail_units(state.active, r, bits, degree)
    entropy = (lower - tail, upper + tail, -2r - bits)
    bounds = _isothermal_add_dyadic((state.delta, state.delta, -state.p), entropy)
    bounds = _isothermal_multiply_dyadic(bounds, temperature)
    for spacing in dx
        bounds = _isothermal_multiply_dyadic(bounds, spacing)
    end
    return bounds
end

function _isothermal_nearunit_energy(
    n::AbstractArray{T},
    temperature::C,
    dx,
    state,
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    (iszero(state.active) || iszero(temperature) || any(iszero, dx)) && return zero(T)
    bits = precision(T) + 16
    bits <= state.cap || return nothing
    degree = Val(T === Float32 ? 8 : 12)
    while true
        lower, upper, exponent =
            _isothermal_complete_dyadic(n, temperature, dx, state, bits, degree)
        lo, hi = _from_exact(T, lower, exponent), _from_exact(T, upper, exponent)
        isequal(lo, hi) && return lo
        bits == state.cap && return nothing
        bits = min(bits + 8, state.cap)
    end
end



function _isothermal_direct_logarithm_bounds(x::T) where {T<:Union{Float32,Float64}}
    invalid = (false, T(NaN), T(NaN))
    isfinite(x) && x > zero(T) || return invalid
    x == one(T) && return true, zero(T), zero(T)
    numerator_lo = prevfloat(x - one(T))
    numerator_hi = nextfloat(x - one(T))
    denominator_lo = prevfloat(x + one(T))
    denominator_hi = nextfloat(x + one(T))
    denominator_lo > zero(T) &&
        isfinite(numerator_lo) &&
        isfinite(numerator_hi) &&
        isfinite(denominator_hi) || return invalid
    zlo = prevfloat(numerator_lo / (numerator_lo >= 0 ? denominator_hi : denominator_lo))
    zhi = nextfloat(numerator_hi / (numerator_hi >= 0 ? denominator_lo : denominator_hi))
    -one(T) < zlo <= zhi < one(T) || return invalid
    if zlo >= zero(T)
        lo = prevfloat(T(2) * zlo)
    else
        denominator = prevfloat(one(T) - nextfloat(zlo * zlo))
        denominator > zero(T) || return invalid
        lo = prevfloat(prevfloat(T(2) * zlo) / denominator)
    end
    if zhi <= zero(T)
        hi = nextfloat(T(2) * zhi)
    else
        denominator = prevfloat(one(T) - nextfloat(zhi * zhi))
        denominator > zero(T) || return invalid
        hi = nextfloat(nextfloat(T(2) * zhi) / denominator)
    end
    return isfinite(lo) && isfinite(hi), lo, hi
end

function _isothermal_logarithm_bounds(x::T) where {T<:Union{Float32,Float64}}
    isfinite(x) && x > zero(T) || return false, T(NaN), T(NaN)
    T(0.5) <= x <= T(2) && return _isothermal_direct_logarithm_bounds(x)
    fraction, exponent = frexp(x)
    mantissa = T(2) * fraction
    exponent -= 1
    good, lower, upper = _isothermal_direct_logarithm_bounds(mantissa)
    logtwo_good, logtwo_lower, logtwo_upper = _isothermal_direct_logarithm_bounds(T(2))
    good && logtwo_good || return false, T(NaN), T(NaN)
    # frexp is exact; the integer exponent fits exactly in Float32/64.
    scaled, elo, ehi = _isothermal_scale_bounds(logtwo_lower, logtwo_upper, T(exponent))
    scaled || return false, T(NaN), T(NaN)
    lo, hi = prevfloat(lower + elo), nextfloat(upper + ehi)
    return isfinite(lo) && isfinite(hi), lo, hi
end

function _isothermal_scale_bounds(lo::T, hi::T, factor) where {T}
    converted = T(factor)
    if converted == factor
        flo = fhi = converted
    else
        flo, fhi = prevfloat(converted), nextfloat(converted)
    end
    isfinite(flo) && isfinite(fhi) || return false, T(NaN), T(NaN)
    a, b, c, d = lo * flo, lo * fhi, hi * flo, hi * fhi
    lower = prevfloat(min(a, b, c, d))
    upper = nextfloat(max(a, b, c, d))
    return isfinite(lower) && isfinite(upper), lower, upper
end

function _isothermal_legacy_energy(n::AbstractArray{T}, temperature, dx) where {T}
    total = zero(T)
    for density in n
        density > zero(T) || continue
        pressure = T(temperature * density)
        total += pressure * log(density)
    end
    return total * prod(dx)
end

function _isothermal_general_probe(
    n::AbstractArray{T},
    temperature::C,
    dx::NTuple{D,T},
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64},D}
    total = correction = lower = upper = 0.0
    finite = isfinite(temperature) && all(isfinite, dx)
    operations = valid = true
    active = false
    te = Float64(temperature)
    for density in n
        finite &= isfinite(density)
        density > zero(T) || continue
        density == one(T) && continue
        iszero(temperature) && continue
        active = true
        value = Float64(density)
        pressure = te * value
        logarithm = log(value)
        term = pressure * logarithm
        total, correction = _amr_compensated_add(total, correction, term)
        operations &=
            _kinetic_normal_product(te, value, pressure) &&
            _kinetic_normal_product(pressure, logarithm, term) &&
            isfinite(total) &&
            isfinite(correction)
        if valid
            good, lo, hi = _isothermal_logarithm_bounds(value)
            if good && isfinite(pressure)
                plo, phi = prevfloat(pressure), nextfloat(pressure)
                a, b, c, d = plo * lo, plo * hi, phi * lo, phi * hi
                lower = prevfloat(lower + prevfloat(min(a, b, c, d)))
                upper = nextfloat(upper + nextfloat(max(a, b, c, d)))
                valid &= isfinite(lower) && isfinite(upper)
            else
                valid = false
            end
        end
    end
    if !finite
        return (;
            value = _isothermal_legacy_energy(n, temperature, dx),
            certified = false,
            finite,
            lower = NaN,
            upper = NaN,
        )
    end
    if !active
        return (; value = zero(T), certified = true, finite, lower = 0.0, upper = 0.0)
    end
    volume = 1.0
    for spacing in dx
        updated = volume * Float64(spacing)
        operations &= _kinetic_normal_product(volume, Float64(spacing), updated)
        volume = updated
        if valid
            valid, lower, upper = _isothermal_scale_bounds(lower, upper, spacing)
        end
    end
    reduced = total + correction
    product = reduced * volume
    result = T(product)
    operations &= isfinite(reduced) && _kinetic_normal_product(reduced, volume, product)
    certified =
        valid &&
        operations &&
        isfinite(result) &&
        lower <= Float64(result) <= upper &&
        -Float64(floatmax(T)) <= lower <= upper <= Float64(floatmax(T)) &&
        (lower >= Float64(floatmin(T)) || upper <= -Float64(floatmin(T)))
    return (; value = result, certified, finite, lower, upper)
end

function _isothermal_general_energy(n, temperature, dx)
    state = _isothermal_general_probe(n, temperature, dx)
    return state.certified || !state.finite ? state.value :
           _isothermal_exact_energy(n, temperature, dx)
end

function _isothermal_energy(
    n::AbstractArray{T},
    temperature::C,
    dx::NTuple{D,T},
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64},D}
    if isfinite(temperature) && all(isfinite, dx)
        state = _isothermal_statistics(n)
        if state !== nothing
            value = _isothermal_nearunit_energy(n, temperature, dx, state)
            return value === nothing ? _isothermal_exact_energy(n, temperature, dx) : value
        end
    end
    return _isothermal_general_energy(n, temperature, dx)
end
