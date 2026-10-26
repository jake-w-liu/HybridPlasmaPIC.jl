# Pure error-free and Float64 double-word arithmetic shared by particle
# kinematics and radiation. Callers own scaling, admission and physical policy.

@inline function _eft_two_sum(a::W, b::W) where {W<:AbstractFloat}
    s = a + b
    t = s - a
    return s, (a - (s - t)) + (b - t)
end

# Requires |a| ≥ |b| or a == 0.
@inline function _eft_fast_two_sum(a::W, b::W) where {W<:AbstractFloat}
    s = a + b
    return s, b - (s - a)
end

@inline function _eft_two_product(a::W, b::W) where {W<:AbstractFloat}
    p = a * b
    return p, fma(a, b, -p)
end


const _DoubleFloat64 = NTuple{2,Float64}

@inline function _dd_add(a::_DoubleFloat64, b::Float64)
    s, e = _eft_two_sum(a[1], b)
    return _eft_fast_two_sum(s, e + a[2])
end

# Accurate double-word addition for normalized operands.
@inline function _dd_add(a::_DoubleFloat64, b::_DoubleFloat64)
    s, e = _eft_two_sum(a[1], b[1])
    t, f = _eft_two_sum(a[2], b[2])
    s, e = _eft_fast_two_sum(s, e + t)
    return _eft_fast_two_sum(s, e + f)
end

@inline _dd_sub(a::_DoubleFloat64, b::_DoubleFloat64) = _dd_add(a, (-b[1], -b[2]))

@inline function _dd_mul(a::_DoubleFloat64, b::_DoubleFloat64)
    p, e = _eft_two_product(a[1], b[1])
    return _eft_fast_two_sum(p, fma(a[1], b[2], fma(a[2], b[1], e)))
end

@inline function _dd_mul(a::_DoubleFloat64, b::Float64)
    p, e = _eft_two_product(a[1], b)
    return _eft_fast_two_sum(p, fma(a[2], b, e))
end

@inline function _dd_div(a::_DoubleFloat64, b::_DoubleFloat64)
    q1 = a[1] / b[1]
    rem = _dd_sub(a, _dd_mul(b, q1))
    q2 = rem[1] / b[1]
    rem = _dd_sub(rem, _dd_mul(b, q2))
    q3 = rem[1] / b[1]
    return _dd_add(_eft_fast_two_sum(q1, q2), q3)
end

@inline function _dd_sqrt(a::_DoubleFloat64)   # a[1] > 0
    x = sqrt(a[1])
    p, e = _eft_two_product(x, x)
    return _eft_fast_two_sum(x, (((a[1] - p) - e) + a[2]) / (2x))
end

# Exact sum of Float64 values as a nonoverlapping expansion (Shewchuk's grow-expansion),
# then rounded to a double-double by accumulating from the smallest component.
@inline _eft_grow(::Tuple{}, b::Float64) = (b,)
@inline function _eft_grow(e::Tuple, b::Float64)
    q, h = _eft_two_sum(b, first(e))
    return (h, _eft_grow(Base.tail(e), q)...)
end
@inline _eft_expansion(t::Tuple{Float64}) = t
@inline _eft_expansion(t::Tuple) = _eft_grow(_eft_expansion(Base.front(t)), last(t))
@inline function _eft_sum_dd(t::Tuple)
    acc = (0.0, 0.0)
    for x in _eft_expansion(t)
        acc = _dd_add(acc, x)
    end
    return acc
end

# Outward-rounded nonnegative bounds. One nextfloat also covers a rounded
# subnormal result; exact zero factors need no artificial uncertainty.
@inline _eft_add_upper(a::Float64, b::Float64) = iszero(a) && iszero(b) ? 0.0 : nextfloat(a + b)
@inline _eft_mul_upper(a::Float64, b::Float64) = iszero(a) || iszero(b) ? 0.0 : nextfloat(a * b)
@inline _eft_div_upper(a::Float64, b::Float64) = iszero(a) ? 0.0 : nextfloat(a / b)

# TwoSum is exact for finite IEEE additions without overflow, including
# subnormals. A nonfinite intermediate makes the bound unresolved (Inf).
# Expansion first cancels leading terms, so a small residual does not inherit
# the scale of the operands that produced it. No DD error theorem is assumed.
@inline function _eft_abs_sum_bound(terms::Tuple)
    all(isfinite, terms) || return Inf
    expansion = _eft_expansion(terms)
    all(isfinite, expansion) || return Inf
    bound = 0.0
    for x in expansion
        bound = _eft_add_upper(bound, abs(x))
    end
    return bound
end
