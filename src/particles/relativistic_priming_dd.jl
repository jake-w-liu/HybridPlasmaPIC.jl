# Accurately rounded publication tier of _rel_prime_repair. Native priming
# can meet a normwise bound yet lose the next stage's Lorentz gap. Every fast
# publication is certified by a residual interval; unresolved cells replay.
@inline _rel_prime_dd_operand(x) = iszero(x) || 0x1p-100 <= abs(x) <= 0x1p100
@inline _rel_prime_word(x) = (Float64(x), 0.0)
const _RelPrimeBall = Tuple{_DoubleFloat64,Float64}
@inline _rel_prime_ball(x) = (_rel_prime_word(x), 0.0)
@inline _rel_prime_unresolved() = ((0.0, 0.0), Inf)
@inline _rel_prime_ball_valid(a::_RelPrimeBall) = all(isfinite, a[1]) && isfinite(a[2]) && a[2] >= 0
@inline _rel_prime_abs_upper(a::_DoubleFloat64) = _eft_add_upper(abs(a[1]), abs(a[2]))
@inline _rel_prime_abs_lower(a::_DoubleFloat64) = prevfloat(abs(a[1]) - abs(a[2]))

# The four products of the exact sums of two word pairs. Each FMA residual
# can lose at most one minimum-subnormal quantum; callers include 4*D.
@inline function _rel_prime_product_words(a::_DoubleFloat64, b::_DoubleFloat64)
    p11 = _eft_two_product(a[1], b[1])
    p12 = _eft_two_product(a[1], b[2])
    p21 = _eft_two_product(a[2], b[1])
    p22 = _eft_two_product(a[2], b[2])
    return (p11..., p12..., p21..., p22...)
end
const _REL_PRIME_PRODUCT_LOSS = 4nextfloat(0.0)

@inline function _rel_prime_ball_add(a::_RelPrimeBall, b::_RelPrimeBall)
    _rel_prime_ball_valid(a) && _rel_prime_ball_valid(b) || return _rel_prime_unresolved()
    z = _dd_add(a[1], b[1])
    residual = _eft_abs_sum_bound((a[1]..., b[1]..., -z[1], -z[2]))
    radius = _eft_add_upper(_eft_add_upper(a[2], b[2]), residual)
    return z, radius
end
@inline _rel_prime_ball_negate(a::_RelPrimeBall) = ((-a[1][1], -a[1][2]), a[2])
@inline _rel_prime_ball_sub(a::_RelPrimeBall, b::_RelPrimeBall) =
    _rel_prime_ball_add(a, _rel_prime_ball_negate(b))

@inline function _rel_prime_ball_mul(a::_RelPrimeBall, b::_RelPrimeBall)
    _rel_prime_ball_valid(a) && _rel_prime_ball_valid(b) || return _rel_prime_unresolved()
    z = _dd_mul(a[1], b[1])
    products = _rel_prime_product_words(a[1], b[1])
    residual =
        _eft_add_upper(_eft_abs_sum_bound((products..., -z[1], -z[2])), _REL_PRIME_PRODUCT_LOSS)
    uncertainty = _eft_add_upper(
        _eft_mul_upper(_rel_prime_abs_upper(a[1]), b[2]),
        _eft_add_upper(
            _eft_mul_upper(_rel_prime_abs_upper(b[1]), a[2]),
            _eft_mul_upper(a[2], b[2]),
        ),
    )
    return z, _eft_add_upper(residual, uncertainty)
end

@inline function _rel_prime_ball_div(a::_RelPrimeBall, b::_RelPrimeBall)
    _rel_prime_ball_valid(a) && _rel_prime_ball_valid(b) || return _rel_prime_unresolved()
    lower = prevfloat(_rel_prime_abs_lower(b[1]) - b[2])
    lower > 0 || return _rel_prime_unresolved()
    z = _dd_div(a[1], b[1])
    products = _rel_prime_product_words(b[1], z)
    residual =
        _eft_add_upper(_eft_abs_sum_bound((a[1]..., map(-, products)...)), _REL_PRIME_PRODUCT_LOSS)
    # For exact inputs A,B in the two balls, A/B-z = (A-B*z)/B.
    numerator = _eft_add_upper(
        residual,
        _eft_add_upper(a[2], _eft_mul_upper(b[2], _rel_prime_abs_upper(z))),
    )
    return z, _eft_div_upper(numerator, lower)
end

@inline function _rel_prime_ball_sqrt(a::_RelPrimeBall)
    _rel_prime_ball_valid(a) && a[1][1] > 0 || return _rel_prime_unresolved()
    prevfloat(_rel_prime_abs_lower(a[1]) - a[2]) > 0 || return _rel_prime_unresolved()
    z = _dd_sqrt(a[1])
    lower = _rel_prime_abs_lower(z)
    z[1] > 0 && lower > 0 || return _rel_prime_unresolved()
    products = _rel_prime_product_words(z, z)
    residual =
        _eft_add_upper(_eft_abs_sum_bound((a[1]..., map(-, products)...)), _REL_PRIME_PRODUCT_LOSS)
    # sqrt(A)-z = (A-z²)/(sqrt(A)+z), and sqrt(A)+z >= z > 0.
    return z, _eft_div_upper(_eft_add_upper(residual, a[2]), lower)
end

@inline function _rel_prime_gamma_ball(v, c)
    scaled_c, gap, squares = _lorentz_gap_words(v, c)
    gap[1] >= 0x1p-900 || return _rel_prime_unresolved()
    residual =
        _eft_add_upper(_eft_abs_sum_bound((squares..., -gap[1], -gap[2])), _REL_PRIME_PRODUCT_LOSS)
    return _rel_prime_ball_div(_rel_prime_ball(scaled_c), _rel_prime_ball_sqrt((gap, residual)))
end

@inline function _rel_prime_structural_zero(v, e, b, qm, h, i)
    j, k = mod1(i + 1, 3), mod1(i + 2, 3)
    return iszero(v[i]) && (
        iszero(h) ||
        iszero(qm) ||
        (iszero(e[i]) && (iszero(v[j]) || iszero(b[k])) && (iszero(v[k]) || iszero(b[j])))
    )
end

# A strict interior-cell comparison is safe with rounded endpoints: a
# representable point strictly below RN(boundary) is below the true boundary.
# The subtraction residual is included; this does not assume DD error bounds.
@inline function _rel_prime_round_cell(::Type{T}, value::_DoubleFloat64, radius) where {T}
    rounded = T(value[1] + value[2])
    isfinite(rounded) || return rounded, false
    below, above = prevfloat(rounded), nextfloat(rounded)
    isfinite(below) && isfinite(above) || return rounded, false
    center = Float64(rounded)
    lower = (center - Float64(below)) / 2
    upper = (Float64(above) - center) / 2
    lower > 0 && upper > 0 || return rounded, false
    delta = _dd_sub(value, (center, 0.0))
    residual = _eft_abs_sum_bound((value..., -center, -delta[1], -delta[2]))
    uncertainty = _eft_add_upper(abs(delta[2]), _eft_add_upper(radius, residual))
    return rounded, -lower + uncertainty < delta[1] < upper - uncertainty
end

@noinline function _rel_prime_dd_attempt(out::NTuple{3,T}, v, e, b, qm::T, h::T, c::T) where {T}
    all(_rel_prime_dd_operand, (v..., e..., b..., qm, h, c)) || return out, false
    gamma0 = _rel_prime_gamma_ball(v, c)
    _rel_prime_ball_valid(gamma0) || return out, false
    vmax, emax, bmax = _abs3(Float64.(v)...), _abs3(Float64.(e)...), _abs3(Float64.(b)...)
    hq = abs(Float64(h) * Float64(qm))
    scale = abs(gamma0[1][1]) * vmax + hq * (emax + 2vmax * bmax)
    isfinite(scale) && scale <= 0x1p300 || return out, false
    V, E, B = map(_rel_prime_ball, v), map(_rel_prime_ball, e), map(_rel_prime_ball, b)
    Qm, H, C = _rel_prime_ball(qm), _rel_prime_ball(h), _rel_prime_ball(c)
    momentum = ntuple(3) do i
        j, k = mod1(i + 1, 3), mod1(i + 2, 3)
        cross = _rel_prime_ball_sub(
            _rel_prime_ball_mul(V[j], B[k]),
            _rel_prime_ball_mul(V[k], B[j]),
        )
        acceleration = _rel_prime_ball_mul(Qm, _rel_prime_ball_add(E[i], cross))
        _rel_prime_ball_sub(_rel_prime_ball_mul(gamma0, V[i]), _rel_prime_ball_mul(H, acceleration))
    end
    magnitude = _rel_prime_ball_add(
        _rel_prime_ball_add(
            _rel_prime_ball_mul(momentum[1], momentum[1]),
            _rel_prime_ball_mul(momentum[2], momentum[2]),
        ),
        _rel_prime_ball_mul(momentum[3], momentum[3]),
    )
    gamma = _rel_prime_ball_sqrt(
        _rel_prime_ball_add(
            _rel_prime_ball_div(magnitude, _rel_prime_ball_mul(C, C)),
            _rel_prime_ball(1.0),
        ),
    )
    _rel_prime_ball_valid(gamma) && gamma[1][1] >= 1 || return out, false
    attempts = ntuple(3) do i
        if _rel_prime_structural_zero(v, e, b, qm, h, i)
            (zero(T), true)
        else
            component = _rel_prime_ball_div(momentum[i], gamma)
            _rel_prime_round_cell(T, component[1], component[2])
        end
    end
    all(last, attempts) || return out, false
    return map(first, attempts), true
end
