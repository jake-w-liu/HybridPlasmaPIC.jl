# Complete algebraic polytropic energy, retaining the stored coefficient and exponent.

function _polytropic_dyadic_equal(coefficient::BigInt, exponent::BigInt, target::Rational{BigInt})
    coefficient >= 0 || throw(ArgumentError("dyadic comparison requires a nonnegative coefficient"))
    target < 0 && return false
    iszero(coefficient) && return iszero(target)
    iszero(target) && return false
    ispow2(denominator(target)) || return false
    a = trailing_zeros(coefficient)
    b = trailing_zeros(numerator(target))
    return (coefficient >> a) == (numerator(target) >> b) &&
           exponent + a == BigInt(b) - (ndigits(denominator(target); base = 2) - 1)
end

function _polytropic_round_sum(::Type{T}, state) where {T<:Union{Float32,Float64}}
    lo, hi, exponent, lower_strict, upper_strict = state
    left = _polytropic_round_scaled(T, lo, exponent)
    right = _polytropic_round_scaled(T, hi, exponent)
    isequal(left, right) && return left
    if isfinite(left) && nextfloat(left) == right
        boundary = _polytropic_midpoint_after(left)
        lower_strict && _polytropic_dyadic_equal(lo, exponent, boundary) && return right
        upper_strict && _polytropic_dyadic_equal(hi, exponent, boundary) && return left
    end
    return nothing
end

function _polytropic_sum_isprime(n::BigInt)
    n < 2 && return false
    n == 2 && return true
    iseven(n) && return false
    limit = isqrt(n)
    divisor = BigInt(3)
    while divisor <= limit
        iszero(rem(n, divisor)) && return false
        divisor += 2
    end
    return true
end
function _polytropic_sum_nextprime(n::BigInt)
    n < 2 && return BigInt(2)
    candidate = n == 2 ? BigInt(3) : n + 2
    while !_polytropic_sum_isprime(candidate)
        candidate += 2
    end
    return candidate
end

function _polytropic_sum_equal(bases, power::BigInt, target::Rational{BigInt})
    power > 0 || throw(ArgumentError("sum equality requires a positive integer exponent"))
    count = 0
    denominator_lcm = BigInt(1)
    largest = firstbase = zero(Rational{BigInt})
    uniform = true
    for base in bases
        base >= 0 || throw(ArgumentError("sum equality requires nonnegative bases"))
        iszero(base) && continue
        count = Base.checked_add(count, 1)
        count == 1 && (firstbase = base)
        uniform &= base == firstbase
        largest = max(largest, base)
        denominator_lcm = lcm(denominator_lcm, denominator(base))
    end
    if iszero(count) || target <= 0
        return (
            equal = iszero(count) && iszero(target),
            primes = 0,
            certified_bits = BigInt(0),
            bound_bits = BigInt(0),
            max_modulus_bits = 0,
        )
    end
    if uniform
        result = _polytropic_power_equal(firstbase, Rational{BigInt}(power), target / count)
        return (
            equal = result,
            primes = 0,
            certified_bits = BigInt(0),
            bound_bits = BigInt(0),
            max_modulus_bits = 0,
        )
    end
    c, d = numerator(target), denominator(target)
    maximum_numerator = numerator(largest) * div(denominator_lcm, denominator(largest))
    magnitude = max(denominator_lcm, maximum_numerator)
    bound_bits = BigInt(ndigits(d * count + c; base = 2)) + power * ndigits(magnitude; base = 2)
    modulus = BigInt(1)
    certified_bits = BigInt(0)
    primes = 0
    while true
        modulus = _polytropic_sum_nextprime(modulus)
        primes = Base.checked_add(primes, 1)
        sum = BigInt(0)
        for base in bases
            iszero(base) && continue
            integer = numerator(base) * div(denominator_lcm, denominator(base))
            sum = mod(sum + powermod(integer, power, modulus), modulus)
        end
        residue = mod(
            mod(d, modulus) * sum - mod(c, modulus) * powermod(denominator_lcm, power, modulus),
            modulus,
        )
        bits = ndigits(modulus; base = 2)
        if !iszero(residue)
            return (
                equal = false,
                primes = primes,
                certified_bits = certified_bits,
                bound_bits = bound_bits,
                max_modulus_bits = bits,
            )
        end
        certified_bits += bits - 1
        if certified_bits >= bound_bits
            return (
                equal = true,
                primes = primes,
                certified_bits = certified_bits,
                bound_bits = bound_bits,
                max_modulus_bits = bits,
            )
        end
    end
end

function _polytropic_energy_term(
    value,
    reference::Rational{BigInt},
    gamma::Rational{BigInt},
    multiplier::Rational{BigInt},
    exponent::BigInt,
    bits::Int,
)
    lo, hi, e = _polytropic_power_bounds(Rational{BigInt}(value) / reference, gamma, bits)
    return lo * multiplier, hi * multiplier, e + exponent
end
function _polytropic_energy_rational_roots(n, reference::Rational{BigInt}, gamma::Rational{BigInt})
    for value in n
        value > 0 || continue
        _polytropic_rational_root(Rational{BigInt}(value) / reference, gamma) === nothing &&
            return false
    end
    return true
end
function _polytropic_energy_root(value, reference::Rational{BigInt}, gamma::Rational{BigInt})
    return _polytropic_rational_root(Rational{BigInt}(value) / reference, gamma)::Rational{BigInt}
end

function _polytropic_exact_energy(
    n::AbstractArray{T},
    c::PolytropicElectrons{C},
    dx,
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64}}
    all(isfinite, n) && all(isfinite, dx) && all(isfinite, (c.pe0, c.n0, c.γ)) ||
        throw(ArgumentError("complete integral component requires finite inputs"))
    c.pe0 >= 0 && c.n0 > 0 && c.γ > 0 ||
        throw(ArgumentError("complete integral component requires valid physical coefficients"))
    c.γ != 1 || throw(ArgumentError("gamma one requires the separate logarithmic integral"))
    count = 0
    first = zero(T)
    uniform = true
    for value in n
        value > 0 || continue
        count = Base.checked_add(count, 1)
        count == 1 && (first = value)
        uniform &= value == first
    end
    if iszero(count) || iszero(c.pe0) || any(iszero, dx)
        negative = (c.γ < 1) ⊻ (!isempty(n) && signbit(c.pe0))
        for spacing in dx
            negative ⊻= signbit(spacing)
        end
        return negative ? -zero(T) : zero(T)
    end
    gamma, reference = Rational{BigInt}(c.γ), Rational{BigInt}(c.n0)
    factor = Rational{BigInt}(c.pe0) * prod(Rational{BigInt}, dx) / (gamma - 1)
    negative = factor < 0
    magnitude = abs(factor)
    if uniform
        value = _polytropic_complete_power(
            T,
            Rational{BigInt}(first) / reference,
            gamma,
            magnitude * count,
        )
        return negative ? -value : value
    end
    _, multiplier, exponent = _polytropic_normalized_factor(magnitude)
    bits = precision(T) + 16
    all_rational = nothing
    rejected_boundary = nothing
    while true
        intervals = (
            _polytropic_energy_term(value, reference, gamma, multiplier, exponent, bits) for
            value in n if value > 0
        )
        state = _polytropic_scaled_sum_state(intervals, bits)
        value = _polytropic_round_sum(T, state)
        value === nothing || return negative ? -value : value
        left = _polytropic_round_scaled(T, state[1], state[3])
        right = _polytropic_round_scaled(T, state[2], state[3])
        if isfinite(left) && nextfloat(left) == right
            boundary = _polytropic_midpoint_after(left)
            if rejected_boundary === nothing || boundary != rejected_boundary
                all_rational === nothing &&
                    (all_rational = _polytropic_energy_rational_roots(n, reference, gamma))
                if all_rational
                    roots = (
                        _polytropic_energy_root(value, reference, gamma) for
                        value in n if value > 0
                    )
                    decision = _polytropic_sum_equal(roots, numerator(gamma), boundary / magnitude)
                    if decision.equal
                        rounded = _from_exact(T, boundary)
                        return negative ? -rounded : rounded
                    end
                end
                rejected_boundary = boundary
            end
        end
        bits = Base.checked_mul(bits, 2)
    end
end

function _polytropic_energy_exponent_bounds(
    maximum_density::T,
    count::Int,
    c::PolytropicElectrons{C},
    dx::NTuple{D,T},
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64},D}
    count > 0 && isfinite(maximum_density) && maximum_density > 0 || return nothing
    all(isfinite, (c.pe0, c.n0, c.γ)) && c.pe0 > 0 && c.n0 > 0 && c.γ > 0 && c.γ != 1 ||
        return nothing
    a, ae = frexp(Float64(maximum_density))
    b, be = frexp(Float64(c.n0))
    k = ae - be - (a < b)
    gap = abs(Float64(c.γ) - 1.0)
    dl, du = prevfloat(gap), nextfloat(gap)
    dl > 0 && isfinite(du) || return nothing
    _, ep = frexp(Float64(c.pe0))
    _, edl = frexp(dl)
    _, edu = frexp(du)
    lower, upper = ep - 1 - edu, ep - edl + 1
    for spacing in dx
        isfinite(spacing) && !iszero(spacing) || return nothing
        _, e = frexp(abs(Float64(spacing)))
        lower = Base.checked_add(lower, e - 1)
        upper = Base.checked_add(upper, e)
    end
    # ceil(log2(count)) from the exact positive integer, including count == 1.
    count_bits = 8sizeof(Int) - leading_zeros(count - 1)
    lo = prevfloat(prevfloat(Float64(c.γ) * k) + lower)
    hi = nextfloat(nextfloat(Float64(c.γ) * (k + 1)) + Base.checked_add(upper, count_bits))
    return lo, hi
end

function _polytropic_energy_probe(
    n::AbstractArray{T},
    c::PolytropicElectrons{C},
    dx::NTuple{D,T},
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64},D}
    all(isfinite, (c.pe0, c.n0, c.γ)) && c.pe0 > 0 && c.n0 > 0 && c.γ > 0 && c.γ != 1 ||
        return nothing
    all(x -> isfinite(x) && !iszero(x), dx) || return nothing
    total = correction = 0.0
    count = 0
    largest = previous = zero(T)
    term = 0.0
    for density in n
        isfinite(density) || return nothing
        density > 0 || continue
        if density != previous
            value = _polytropic_pressure_probe_power(Float64, density, c, Float64)
            value === nothing && return nothing
            term = value
            previous = density
        end
        total, correction = _amr_compensated_add(total, correction, term)
        isfinite(total) && isfinite(correction) || return nothing
        largest = max(largest, density)
        count = Base.checked_add(count, 1)
    end
    count > 0 || return nothing
    bounds = _polytropic_energy_exponent_bounds(largest, count, c, dx)
    bounds === nothing && return nothing
    bounds[1] >= exponent(floatmin(T)) + 2 && bounds[2] <= exponent(floatmax(T)) - 2 ||
        return nothing
    value = (total + correction) / (Float64(c.γ) - 1.0)
    isfinite(value) && abs(value) >= 2floatmin(Float64) || return nothing
    for spacing in dx
        value *= Float64(spacing)
        isfinite(value) && abs(value) >= 2floatmin(Float64) || return nothing
    end
    result = T(value)
    isfinite(result) && abs(result) >= 2floatmin(T) || return nothing
    return result
end

function _polytropic_log_energy(
    n::AbstractArray{T},
    c::PolytropicElectrons{C},
    dx::NTuple{D,T},
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64},D}
    all(isfinite, n) && all(isfinite, dx) && all(isfinite, (c.pe0, c.n0, c.γ)) ||
        throw(ArgumentError("gamma-one component requires finite inputs"))
    c.γ == 1 && c.n0 > 0 && c.pe0 >= 0 ||
        throw(ArgumentError("gamma-one component requires valid stored coefficients"))
    state = _isothermal_statistics(n)
    if state !== nothing
        value = _isothermal_nearunit_energy(n, c.pe0, dx, state, c.n0)
        value === nothing || return value
    else
        probe = _isothermal_general_probe(n, c.pe0, dx, c.n0)
        if probe.certified
            lower, upper = T(probe.lower), T(probe.upper)
            isequal(lower, upper) && return lower
        end
    end
    value = _isothermal_fixed_probe(n, c.pe0, dx, c.n0)
    value === nothing || return value
    return _isothermal_exact_energy(n, Rational{BigInt}(c.pe0) / Rational{BigInt}(c.n0), dx)
end

function _polytropic_energy(
    n::AbstractArray{T},
    c::PolytropicElectrons{C},
    dx::NTuple{D,T},
) where {T<:Union{Float32,Float64},C<:Union{Float32,Float64},D}
    c.γ == 1 && return _polytropic_log_energy(n, c, dx)
    value = _polytropic_energy_probe(n, c, dx)
    return value === nothing ? _polytropic_exact_energy(n, c, dx) : value
end
