# Residual arithmetic for the one canonical relativistic Boris expression.
# This private payload is used only by that expression, in nearest/gradual
# IEEE arithmetic. Its DD tier reuses the existing priming residual owners.
struct _RelBorisBound{T<:Union{_RangeGuardedFloat,_DoubleFloat64}} <: Number
    center::T
    radius::Float64
end

@inline _rel_boris_ball(a::_RelBorisBound{T}) where {T<:_RangeGuardedFloat} =
    ((Float64(a.center), 0.0), a.radius)
@inline _rel_boris_ball(a::_RelBorisBound{_DoubleFloat64}) = (a.center, a.radius)
@inline _rel_boris_dd(x) = _RelBorisBound(_rel_prime_word(x), 0.0)
@inline _rel_boris_dd_ball(a::_RelPrimeBall) = _RelBorisBound(a[1], a[2])
@inline _rel_boris_valid(a) = _rel_prime_ball_valid(_rel_boris_ball(a))
@inline _rel_boris_zero(a) = iszero(a.radius) && all(iszero, _rel_boris_ball(a)[1])
@inline _rel_boris_upper(a) =
    _eft_add_upper(_rel_prime_abs_upper(_rel_boris_ball(a)[1]), a.radius)
@inline _rel_boris_lower(a) =
    max(0.0, prevfloat(_rel_prime_abs_lower(_rel_boris_ball(a)[1]) - a.radius))

@inline _rel_boris_constant(::Type{T}, x) where {T<:_RangeGuardedFloat} =
    _RelBorisBound(T(x), 0.0)
@inline _rel_boris_constant(::Type{_DoubleFloat64}, x) = _rel_boris_dd(x)
Base.one(::Type{_RelBorisBound{T}}) where {T} = _rel_boris_constant(T, 1)

@inline function Base.:+(a::_RelBorisBound{T}, b::_RelBorisBound{T}) where {T<:_RangeGuardedFloat}
    z = a.center + b.center
    residual = _eft_abs_sum_bound((Float64(a.center), Float64(b.center), -Float64(z)))
    return _RelBorisBound(z, _boris_sum_error(a.radius, b.radius, residual))
end
@inline Base.:-(a::_RelBorisBound) = _RelBorisBound(-a.center, a.radius)
@inline Base.:-(a::_RelBorisBound, b::_RelBorisBound) = a + (-b)

@inline function Base.:*(a::_RelBorisBound{T}, b::_RelBorisBound{T}) where {T<:_RangeGuardedFloat}
    z = a.center * b.center
    (_rel_boris_zero(a) || _rel_boris_zero(b)) && return _RelBorisBound(z, 0.0)
    x, y = Float64(a.center), Float64(b.center)
    product = _eft_two_product(x, y)
    residual = _eft_add_upper(
        _eft_abs_sum_bound((product..., -Float64(z))), nextfloat(0.0),
    )
    return _RelBorisBound(z, _boris_product_error(x, y, a.radius, b.radius, residual))
end

@inline function Base.:/(a::_RelBorisBound{T}, b::_RelBorisBound{T}) where {T<:_RangeGuardedFloat}
    z = a.center / b.center
    lower = _rel_boris_lower(b)
    lower > 0 || return _RelBorisBound(z, Inf)
    _rel_boris_zero(a) && return _RelBorisBound(z, 0.0)
    product = _eft_two_product(Float64(b.center), Float64(z))
    residual = _eft_add_upper(
        _eft_abs_sum_bound((Float64(a.center), -product[1], -product[2])), nextfloat(0.0),
    )
    numerator = _eft_add_upper(
        residual, _eft_add_upper(a.radius, _eft_mul_upper(abs(Float64(z)), b.radius)),
    )
    return _RelBorisBound(z, _eft_div_upper(numerator, lower))
end

@inline function Base.sqrt(a::_RelBorisBound{T}) where {T<:_RangeGuardedFloat}
    z = sqrt(a.center)
    _rel_boris_zero(a) && return _RelBorisBound(z, 0.0)
    z > 0 || return _RelBorisBound(z, Inf)
    product = _eft_two_product(Float64(z), Float64(z))
    residual = _eft_add_upper(
        _eft_abs_sum_bound((Float64(a.center), -product[1], -product[2])), nextfloat(0.0),
    )
    # sqrt(A)-z = (A-z²)/(sqrt(A)+z), whose denominator is at least z.
    return _RelBorisBound(z, _eft_div_upper(_eft_add_upper(residual, a.radius), Float64(z)))
end

@inline Base.:+(a::_RelBorisBound{_DoubleFloat64}, b::_RelBorisBound{_DoubleFloat64}) =
    _rel_boris_dd_ball(_rel_prime_ball_add(_rel_boris_ball(a), _rel_boris_ball(b)))
@inline Base.:-(a::_RelBorisBound{_DoubleFloat64}) =
    _rel_boris_dd_ball(_rel_prime_ball_negate(_rel_boris_ball(a)))
@inline function Base.:*(a::_RelBorisBound{_DoubleFloat64}, b::_RelBorisBound{_DoubleFloat64})
    (_rel_boris_zero(a) || _rel_boris_zero(b)) && return _rel_boris_dd(0)
    return _rel_boris_dd_ball(_rel_prime_ball_mul(_rel_boris_ball(a), _rel_boris_ball(b)))
end
@inline function Base.:/(a::_RelBorisBound{_DoubleFloat64}, b::_RelBorisBound{_DoubleFloat64})
    _rel_boris_lower(b) > 0 && _rel_boris_zero(a) && return _rel_boris_dd(0)
    return _rel_boris_dd_ball(_rel_prime_ball_div(_rel_boris_ball(a), _rel_boris_ball(b)))
end
@inline Base.sqrt(a::_RelBorisBound{_DoubleFloat64}) =
    _rel_boris_zero(a) ? _rel_boris_dd(0) :
    _rel_boris_dd_ball(_rel_prime_ball_sqrt(_rel_boris_ball(a)))

for operation in (:+, :-, :*, :/)
    @eval begin
        @inline Base.$operation(a::_RelBorisBound{T}, b::Integer) where {T} =
            $operation(a, _rel_boris_constant(T, b))
        @inline Base.$operation(a::Integer, b::_RelBorisBound{T}) where {T} =
            $operation(_rel_boris_constant(T, a), b)
    end
end

@inline function _rel_boris_factor(factor::_BorisSpeciesFactors, make)
    return make(factor.q) / make(factor.m)
end
@inline _rel_boris_factor(factor::Number, make) = make(factor)

@inline function _rel_boris_distance(a, x)
    isfinite(x) || return Inf
    delta = _rel_prime_ball_sub(_rel_boris_ball(a), _rel_prime_ball(x))
    return _eft_add_upper(_rel_prime_abs_upper(delta[1]), delta[2])
end
@inline function _rel_boris_distance(a::_RelBorisBound{T}, x::T) where {T<:_RangeGuardedFloat}
    isequal(a.center, x) && return a.radius
    isfinite(x) || return Inf
    return _eft_add_upper(_eft_abs_sum_bound((Float64(a.center), -Float64(x))), a.radius)
end

@inline function _rel_boris_norm(a)
    squared = a[1] * a[1] + a[2] * a[2] + a[3] * a[3]
    return sqrt(squared)
end

# Fixed aggregate magnitude/error pairs bound the same finite-c expression.
# Denominator lower bounds belong to the particular positive Boris stages;
# no cancellation or generic interval policy is inferred from a magnitude pair.
@inline function _rel_boris_batch_divide(a, b, native_lower, exact_lower, u, D)
    native_lower > 0 && exact_lower > 0 || return (Inf, Inf)
    quotient = _boris_upper_div(a[1], native_lower)
    rounding = iszero(a[1]) ? 0.0 : _boris_upper_add(_boris_upper_mul(u, quotient), D)
    propagated = _boris_upper_div(
        _boris_upper_add(a[2], _boris_upper_mul(quotient, b[2])), exact_lower,
    )
    return _boris_upper_add(quotient, rounding), _boris_upper_add(propagated, rounding)
end

@inline function _rel_boris_batch_sqrt(a, native_lower, exact_lower, u, D)
    native_lower > 0 && exact_lower > 0 || return (Inf, Inf)
    upper = nextfloat(sqrt(a[1]))
    rounding = _boris_upper_add(_boris_upper_mul(u, upper), D)
    denominator = max(0.0, prevfloat(prevfloat(sqrt(native_lower)) + prevfloat(sqrt(exact_lower))))
    return _boris_upper_add(upper, rounding),
        _boris_upper_add(_boris_upper_div(a[2], denominator), rounding)
end

@inline function _boris_rel_batch_accuracy(h::T, factor, dt::T, c::T, vmax, emax, bmax, smin) where {T}
    _boris_accuracy_mode() || return false
    eligible, coefficient, _ = _boris_coefficient(factor, dt, h)
    eligible || return false
    # The reduction omits native-zero scales. Such a row has v=0 and
    # RN(hE_i)=0, hence |H E_i|<=(1+kappa)D/2. Exact magnetic rotation is
    # orthogonal and final gamma>=1, so norm(out)<=sqrt(3)(1+kappa)D<4D.
    # It therefore satisfies the quantum-floor budget independently of smin.
    one_plus = _boris_upper_add(1.0, coefficient)
    _boris_upper_mul(3.0, _boris_upper_mul(one_plus, one_plus)) < 16.0 || return false
    isinf(smin) && return true
    u, D = ldexp(1.0, -precision(T)), Float64(nextfloat(zero(T)))
    velocity, electric, magnetic = (Float64(vmax), 0.0), (Float64(emax), 0.0), (Float64(bmax), 0.0)
    C = (Float64(c), 0.0)
    c2 = c * c
    C2 = (Float64(c2), _boris_local_rounding(Float64(c2), u, D))
    exact_c2_lower = _boris_lower_product(Float64(c), Float64(c))
    beta = _rel_boris_batch_divide(velocity, C, C[1], C[1], u, D)
    beta_square = _boris_batch_product(beta, beta, u, D)
    beta_sum = _boris_batch_sum(_boris_batch_sum(beta_square, beta_square, u, D), beta_square, u, D)
    # The range admission already requires vmax<=c/4. Both exact and native
    # gaps are positive; subtraction rounds a value no larger than one.
    subtraction_rounding = _boris_upper_add(u, D)
    gap_radius = _boris_upper_add(beta_sum[2], subtraction_rounding)
    native_gap_lower = _boris_lower_difference(
        _boris_lower_difference(1.0, beta_sum[1]), subtraction_rounding,
    )
    exact_gap_lower = _boris_lower_difference(native_gap_lower, gap_radius)
    root = _rel_boris_batch_sqrt((1.0, gap_radius), native_gap_lower, exact_gap_lower, u, D)
    native_root_lower = _boris_lower_difference(
        _boris_lower_product(max(0.0, prevfloat(sqrt(native_gap_lower))), 1.0 - u), D,
    )
    exact_root_lower = max(0.0, prevfloat(sqrt(exact_gap_lower)))
    gamma0 = _rel_boris_batch_divide((1.0, 0.0), root, native_root_lower, exact_root_lower, u, D)
    initial = _boris_batch_product(gamma0, velocity, u, D)
    H = (abs(Float64(h)), _boris_upper_mul(coefficient, abs(Float64(h))))
    kick = _boris_batch_product(H, electric, u, D)
    minus = _boris_batch_sum(initial, kick, u, D)
    square = _boris_batch_product(minus, minus, u, D)
    squares = _boris_batch_sum(_boris_batch_sum(square, square, u, D), square, u, D)
    ratio = _rel_boris_batch_divide(squares, C2, C2[1], exact_c2_lower, u, D)
    gamma_minus = _rel_boris_batch_sqrt(_boris_batch_sum((1.0, 0.0), ratio, u, D), 1.0, 1.0, u, D)
    t = _rel_boris_batch_divide(_boris_batch_product(H, magnetic, u, D), gamma_minus, 1.0, 1.0, u, D)
    t_square = _boris_batch_product(t, t, u, D)
    t2 = _boris_batch_sum(_boris_batch_sum(t_square, t_square, u, D), t_square, u, D)
    denominator = _boris_batch_sum((1.0, 0.0), t2, u, D)
    f = _rel_boris_batch_divide((2.0, 0.0), denominator, 1.0, 1.0, u, D)
    s = _boris_batch_product(f, t, u, D)
    cross_product = _boris_batch_product(minus, t, u, D)
    prime = _boris_batch_sum(minus, _boris_batch_sum(cross_product, cross_product, u, D), u, D)
    rotated_product = _boris_batch_product(prime, s, u, D)
    rotated = _boris_batch_sum(minus, _boris_batch_sum(rotated_product, rotated_product, u, D), u, D)
    final = _boris_batch_sum(rotated, kick, u, D)
    final_square = _boris_batch_product(final, final, u, D)
    final_squares = _boris_batch_sum(_boris_batch_sum(final_square, final_square, u, D), final_square, u, D)
    final_ratio = _rel_boris_batch_divide(final_squares, C2, C2[1], exact_c2_lower, u, D)
    gamma_final = _rel_boris_batch_sqrt(_boris_batch_sum((1.0, 0.0), final_ratio, u, D), 1.0, 1.0, u, D)
    output = _rel_boris_batch_divide(final, gamma_final, 1.0, 1.0, u, D)
    # smin bounds max(max|v|,max|RN(hE)|) for each nonzero row.
    # Exact gamma0>=1, so subtracting electric error gives a lower bound on
    # max(norm(gamma0*v),norm(hE)); divide by an upper exact gamma_final.
    lower_term = _boris_lower_difference(Float64(smin), kick[2])
    gamma_upper = _boris_upper_add(gamma_final[1], gamma_final[2])
    lower = isfinite(gamma_upper) && gamma_upper > 0 && lower_term > 0 ?
        max(0.0, prevfloat(lower_term / gamma_upper)) : 0.0
    return _euler_accuracy_finish(output, lower, T)
end

# Preserve the historical range predicate's truthful contract. Only this
# joint owner authorizes native bypass in actual EMPIC kick loops.
@inline function _boris_rel_batch_admission(h::T, factor, dt::T, c::T, vmax, emax, bmax, smin) where {T}
    _boris_rel_batch_certificate(h, factor, dt, c, vmax, emax, bmax, smin) || return false
    T <: _RangeGuardedFloat || return true
    return _boris_rel_batch_accuracy(h, factor, dt, c, vmax, emax, bmax, smin)
end

@inline function _rel_boris_scale_lower(values, initial, parts, ::Val{full}) where {full}
    electric = (parts[2], parts[3], parts[4])
    gamma = _rel_boris_upper(parts[11])
    isfinite(gamma) && gamma > 0 || return 0.0
    if full
        output_lower = _rel_boris_lower(_rel_boris_norm(values))
        initial_lower = _rel_boris_lower(_rel_boris_norm(initial))
        electric_lower = _rel_boris_lower(_rel_boris_norm(electric))
    else
        output_lower = maximum(_rel_boris_lower, values)
        initial_lower = maximum(_rel_boris_lower, initial)
        electric_lower = maximum(_rel_boris_lower, electric)
    end
    lower = max(initial_lower, electric_lower)
    divided = iszero(lower) ? 0.0 : max(0.0, prevfloat(lower / gamma))
    return max(output_lower, divided)
end

@inline function _rel_boris_accuracy_ok(out::NTuple{3,T}, values, lower) where {T}
    all(_rel_boris_valid, values) || return false
    isfinite(lower) && lower >= 0 || return false
    e = iszero(lower) ? exponent(floatmin(T)) : max(exponent(lower), exponent(floatmin(T)))
    tolerance = ldexp(1.0, e - precision(T) + 3)
    return all(1:3) do i
        error = _rel_boris_distance(values[i], out[i])
        isfinite(error) && error <= tolerance &&
            _boris_upper_add(abs(Float64(out[i])), error) <= Float64(floatmax(T))
    end
end

# Radius around the actual historical initial gamma, from the original v,c.
@inline function _boris_rel_accuracy_native(out, ginit::T, v, E, B, factor, dt::T, c::T) where {T}
    _boris_accuracy_mode() || return false
    all(isfinite, (v..., E..., B..., _boris_factor_operands(factor)..., dt, c)) || return false
    # Retain the existing gamma-ball input proof. Binary scaling outside
    # this window needs a separate scaled-input loss bound, not an assumption.
    all(_rel_prime_dd_operand, (v..., c)) || return false
    gamma = _rel_prime_gamma_ball(v, c)
    _rel_prime_ball_valid(gamma) || return false
    G = _RelBorisBound(ginit, _rel_boris_distance(_rel_boris_dd_ball(gamma), ginit))
    make(x) = _RelBorisBound(x, 0.0)
    initial = map(x -> G * make(x), v)
    fields = map(make, (E..., B...))
    coefficient = _rel_boris_factor(factor, make)
    C = make(c)
    values, parts = _boris_rel_momentum(initial..., fields..., coefficient, make(dt), C * C)
    # The adapter must preserve every native operation, including literal
    # promotion and the original three-term addition order.
    isequal(map(x -> x.center, values), out) || return false
    return _rel_boris_accuracy_ok(out, values, _rel_boris_scale_lower(values, initial, parts, Val(false)))
end

# This is a bounded, residual-certified repair; the inherited extreme replay
# remains separate. Each value has an enclosure, not an assumed DD error bound.
@noinline function _boris_rel_dd_attempt(out::NTuple{3,T}, v, E, B, factor, dt::T, c::T) where {T}
    _boris_accuracy_mode() || return out, false
    iszero(c) && return out, false
    all(_rel_prime_dd_operand, (v..., E..., B..., _boris_factor_operands(factor)..., dt, c)) ||
        return out, false
    gamma = _rel_prime_gamma_ball(v, c)
    _rel_prime_ball_valid(gamma) || return out, false
    G = _rel_boris_dd_ball(gamma)
    initial = map(x -> G * _rel_boris_dd(x), v)
    fields = map(_rel_boris_dd, (E..., B...))
    coefficient = _rel_boris_factor(factor, _rel_boris_dd)
    C = _rel_boris_dd(c)
    values, parts = _boris_rel_momentum(initial..., fields..., coefficient, _rel_boris_dd(dt), C * C)
    all(_rel_boris_valid, (initial..., values..., parts...)) || return out, false
    lower = _rel_boris_scale_lower(values, initial, parts, Val(true))
    _supported_relativistic_output(out, c) && _rel_boris_accuracy_ok(out, values, lower) &&
        return out, true
    candidate = map(x -> T(x.center[1] + x.center[2]), values)
    return candidate,
        _supported_relativistic_output(candidate, c) && _rel_boris_accuracy_ok(candidate, values, lower)
end
