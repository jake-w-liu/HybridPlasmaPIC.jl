# Sufficient CPU accuracy bounds for the canonical Boris and Euler expressions.
# All arithmetic below is positive Float64 outward bounding, not a general
# interval API. The fixed Boris nodes must stay DAG-equal to boris.jl.
struct _BorisSpeciesFactors{T<:_RangeGuardedFloat}
    q::T
    m::T
    qm::T
end
@inline _boris_native_factor(qm::Number) = qm
@inline _boris_native_factor(f::_BorisSpeciesFactors) = f.qm
@inline _boris_zero_factor(qm::Number) = iszero(qm)
@inline _boris_zero_factor(f::_BorisSpeciesFactors) = iszero(f.q)
_boris_exact_factor(qm::Number) = _exact_rational(qm)
_boris_exact_factor(f::_BorisSpeciesFactors) = _exact_rational(f.q) / _exact_rational(f.m)
@inline _boris_accuracy_mode() = rounding(Float64) === RoundNearest && !get_zero_subnormals()

# Euler range certificates must distinguish an original zero charge from a
# quotient rounded to zero. Normal quotients keep the historical native path.
@inline _euler_factor_in_range(qm::Number, ::Type) = true
@inline function _euler_factor_in_range(f::_BorisSpeciesFactors{T}, ::Type{T}) where {T}
    return _boris_accuracy_mode() && (iszero(f.q) || (isfinite(f.qm) && abs(f.qm) >= floatmin(T)))
end

@inline _boris_bound_operands(a, b) = isfinite(a) && isfinite(b) && a >= 0 && b >= 0
@inline function _boris_upper_add(a::Float64, b::Float64)
    _boris_bound_operands(a, b) || return Inf
    iszero(a) && return b
    iszero(b) && return a
    return nextfloat(a + b)
end
@inline function _boris_upper_mul(a::Float64, b::Float64)
    _boris_bound_operands(a, b) || return Inf
    (iszero(a) || iszero(b)) && return 0.0
    return nextfloat(a * b)
end
@inline function _boris_upper_div(a::Float64, b::Float64)
    (_boris_bound_operands(a, b) && b > 0) || return Inf
    iszero(a) && return 0.0
    return nextfloat(a / b)
end
@inline function _boris_lower_difference(a::Float64, b::Float64)
    _boris_bound_operands(a, b) || return 0.0
    iszero(b) && return a
    a <= b && return 0.0
    return max(0.0, prevfloat(a - b))
end
@inline function _boris_lower_product(a::Float64, b::Float64)
    _boris_bound_operands(a, b) || return 0.0
    (iszero(a) || iszero(b)) && return 0.0
    return max(0.0, prevfloat(a * b))
end
@inline function _boris_relative_coefficient(factor, ::Type{T}) where {T}
    u = ldexp(1.0, -precision(T))
    return u
end
@inline function _boris_relative_coefficient(factor::_BorisSpeciesFactors, ::Type{T}) where {T}
    u = ldexp(1.0, -precision(T))
    return _boris_upper_add(_boris_upper_mul(2.0, u), _boris_upper_mul(u, u))
end
@inline function _boris_coefficient(factor, dt::T, h::T) where {T}
    qm = _boris_native_factor(factor)
    truezero = _boris_zero_factor(factor) || iszero(dt)
    normal_quotient = !(factor isa _BorisSpeciesFactors) || (isfinite(qm) && abs(qm) >= floatmin(T))
    eligible =
        isfinite(qm) && isfinite(h) && ((abs(h) >= floatmin(T) && normal_quotient) || truezero)
    return eligible, truezero ? 0.0 : _boris_relative_coefficient(factor, T), truezero
end
@inline _boris_local_rounding(center, u, D) = _boris_upper_add(_boris_upper_mul(u, abs(center)), D)
@inline _boris_sum_error(ea, eb, rho) = _boris_upper_add(_boris_upper_add(ea, eb), rho)
@inline function _boris_product_error(a, b, ea, eb, rho)
    return _boris_upper_add(
        _boris_upper_add(
            _boris_upper_add(_boris_upper_mul(abs(a), eb), _boris_upper_mul(abs(b), ea)),
            _boris_upper_mul(ea, eb),
        ),
        rho,
    )
end
@inline function _boris_reciprocal_error(d, ed, rho)
    lower = max(1.0, _boris_lower_difference(d, ed))
    raw_product = d * lower
    denominator = isfinite(raw_product) ? _boris_lower_product(d, lower) : 0.0
    return denominator > 0 ?
           _boris_upper_add(_boris_upper_div(_boris_upper_mul(2.0, ed), denominator), rho) : Inf
end

# Every radius corresponds to the same numbered actual native operation.
@inline function _boris_scalar_radii(nodes::NTuple{48,T}, v, E, B, coefficient, truezero) where {T}
    u, D = ldexp(1.0, -precision(T)), Float64(nextfloat(zero(T)))
    vx, vy, vz = Float64.(v)
    Ex, Ey, Ez = Float64.(E)
    Bx, By, Bz = Float64.(B)
    z1 = truezero
    r1 = _boris_upper_mul(coefficient, abs(Float64(nodes[1])))
    z2 = truezero
    r2 = _boris_upper_mul(coefficient, abs(Float64(nodes[2])))
    z3 = (z2 || iszero(Ex))
    r3 =
        (z2 || iszero(Ex)) ? 0.0 :
        _boris_product_error(
            Float64(nodes[2]),
            Ex,
            r2,
            0.0,
            ((z2 || iszero(Ex)) ? 0.0 : _boris_local_rounding(Float64(nodes[3]), u, D)),
        )
    z4 = (z2 || iszero(Ey))
    r4 =
        (z2 || iszero(Ey)) ? 0.0 :
        _boris_product_error(
            Float64(nodes[2]),
            Ey,
            r2,
            0.0,
            ((z2 || iszero(Ey)) ? 0.0 : _boris_local_rounding(Float64(nodes[4]), u, D)),
        )
    z5 = (z2 || iszero(Ez))
    r5 =
        (z2 || iszero(Ez)) ? 0.0 :
        _boris_product_error(
            Float64(nodes[2]),
            Ez,
            r2,
            0.0,
            ((z2 || iszero(Ez)) ? 0.0 : _boris_local_rounding(Float64(nodes[5]), u, D)),
        )
    z6 = (iszero(vx) && z3)
    r6 = _boris_sum_error(
        0.0,
        r3,
        ((iszero(vx) || z3) ? 0.0 : _boris_local_rounding(Float64(nodes[6]), u, D)),
    )
    z7 = (iszero(vy) && z4)
    r7 = _boris_sum_error(
        0.0,
        r4,
        ((iszero(vy) || z4) ? 0.0 : _boris_local_rounding(Float64(nodes[7]), u, D)),
    )
    z8 = (iszero(vz) && z5)
    r8 = _boris_sum_error(
        0.0,
        r5,
        ((iszero(vz) || z5) ? 0.0 : _boris_local_rounding(Float64(nodes[8]), u, D)),
    )
    z9 = (z2 || iszero(Bx))
    r9 =
        (z2 || iszero(Bx)) ? 0.0 :
        _boris_product_error(
            Float64(nodes[2]),
            Bx,
            r2,
            0.0,
            ((z2 || iszero(Bx)) ? 0.0 : _boris_local_rounding(Float64(nodes[9]), u, D)),
        )
    z10 = (z2 || iszero(By))
    r10 =
        (z2 || iszero(By)) ? 0.0 :
        _boris_product_error(
            Float64(nodes[2]),
            By,
            r2,
            0.0,
            ((z2 || iszero(By)) ? 0.0 : _boris_local_rounding(Float64(nodes[10]), u, D)),
        )
    z11 = (z2 || iszero(Bz))
    r11 =
        (z2 || iszero(Bz)) ? 0.0 :
        _boris_product_error(
            Float64(nodes[2]),
            Bz,
            r2,
            0.0,
            ((z2 || iszero(Bz)) ? 0.0 : _boris_local_rounding(Float64(nodes[11]), u, D)),
        )
    z12 = (z9 || z9)
    r12 =
        (z9 || z9) ? 0.0 :
        _boris_product_error(
            Float64(nodes[9]),
            Float64(nodes[9]),
            r9,
            r9,
            ((z9 || z9) ? 0.0 : _boris_local_rounding(Float64(nodes[12]), u, D)),
        )
    z13 = (z10 || z10)
    r13 =
        (z10 || z10) ? 0.0 :
        _boris_product_error(
            Float64(nodes[10]),
            Float64(nodes[10]),
            r10,
            r10,
            ((z10 || z10) ? 0.0 : _boris_local_rounding(Float64(nodes[13]), u, D)),
        )
    z14 = (z11 || z11)
    r14 =
        (z11 || z11) ? 0.0 :
        _boris_product_error(
            Float64(nodes[11]),
            Float64(nodes[11]),
            r11,
            r11,
            ((z11 || z11) ? 0.0 : _boris_local_rounding(Float64(nodes[14]), u, D)),
        )
    z15 = (z12 && z13)
    r15 = _boris_sum_error(
        r12,
        r13,
        ((z12 || z13) ? 0.0 : _boris_local_rounding(Float64(nodes[15]), u, D)),
    )
    z16 = (z15 && z14)
    r16 = _boris_sum_error(
        r15,
        r14,
        ((z15 || z14) ? 0.0 : _boris_local_rounding(Float64(nodes[16]), u, D)),
    )
    z17 = (false && z16)
    r17 = _boris_sum_error(
        0.0,
        r16,
        ((false || z16) ? 0.0 : _boris_local_rounding(Float64(nodes[17]), u, D)),
    )
    z18 = false
    r18 = _boris_reciprocal_error(
        Float64(nodes[17]),
        r17,
        _boris_local_rounding(Float64(nodes[18]), u, D),
    )
    z19 = (z18 || z9)
    r19 =
        (z18 || z9) ? 0.0 :
        _boris_product_error(
            Float64(nodes[18]),
            Float64(nodes[9]),
            r18,
            r9,
            ((z18 || z9) ? 0.0 : _boris_local_rounding(Float64(nodes[19]), u, D)),
        )
    z20 = (z18 || z10)
    r20 =
        (z18 || z10) ? 0.0 :
        _boris_product_error(
            Float64(nodes[18]),
            Float64(nodes[10]),
            r18,
            r10,
            ((z18 || z10) ? 0.0 : _boris_local_rounding(Float64(nodes[20]), u, D)),
        )
    z21 = (z18 || z11)
    r21 =
        (z18 || z11) ? 0.0 :
        _boris_product_error(
            Float64(nodes[18]),
            Float64(nodes[11]),
            r18,
            r11,
            ((z18 || z11) ? 0.0 : _boris_local_rounding(Float64(nodes[21]), u, D)),
        )
    z22 = (z7 || z11)
    r22 =
        (z7 || z11) ? 0.0 :
        _boris_product_error(
            Float64(nodes[7]),
            Float64(nodes[11]),
            r7,
            r11,
            ((z7 || z11) ? 0.0 : _boris_local_rounding(Float64(nodes[22]), u, D)),
        )
    z23 = (z8 || z10)
    r23 =
        (z8 || z10) ? 0.0 :
        _boris_product_error(
            Float64(nodes[8]),
            Float64(nodes[10]),
            r8,
            r10,
            ((z8 || z10) ? 0.0 : _boris_local_rounding(Float64(nodes[23]), u, D)),
        )
    z24 = (z22 && z23)
    r24 = _boris_sum_error(
        r22,
        r23,
        ((z22 || z23) ? 0.0 : _boris_local_rounding(Float64(nodes[24]), u, D)),
    )
    z25 = (z6 && z24)
    r25 = _boris_sum_error(
        r6,
        r24,
        ((z6 || z24) ? 0.0 : _boris_local_rounding(Float64(nodes[25]), u, D)),
    )
    z26 = (z8 || z9)
    r26 =
        (z8 || z9) ? 0.0 :
        _boris_product_error(
            Float64(nodes[8]),
            Float64(nodes[9]),
            r8,
            r9,
            ((z8 || z9) ? 0.0 : _boris_local_rounding(Float64(nodes[26]), u, D)),
        )
    z27 = (z6 || z11)
    r27 =
        (z6 || z11) ? 0.0 :
        _boris_product_error(
            Float64(nodes[6]),
            Float64(nodes[11]),
            r6,
            r11,
            ((z6 || z11) ? 0.0 : _boris_local_rounding(Float64(nodes[27]), u, D)),
        )
    z28 = (z26 && z27)
    r28 = _boris_sum_error(
        r26,
        r27,
        ((z26 || z27) ? 0.0 : _boris_local_rounding(Float64(nodes[28]), u, D)),
    )
    z29 = (z7 && z28)
    r29 = _boris_sum_error(
        r7,
        r28,
        ((z7 || z28) ? 0.0 : _boris_local_rounding(Float64(nodes[29]), u, D)),
    )
    z30 = (z6 || z10)
    r30 =
        (z6 || z10) ? 0.0 :
        _boris_product_error(
            Float64(nodes[6]),
            Float64(nodes[10]),
            r6,
            r10,
            ((z6 || z10) ? 0.0 : _boris_local_rounding(Float64(nodes[30]), u, D)),
        )
    z31 = (z7 || z9)
    r31 =
        (z7 || z9) ? 0.0 :
        _boris_product_error(
            Float64(nodes[7]),
            Float64(nodes[9]),
            r7,
            r9,
            ((z7 || z9) ? 0.0 : _boris_local_rounding(Float64(nodes[31]), u, D)),
        )
    z32 = (z30 && z31)
    r32 = _boris_sum_error(
        r30,
        r31,
        ((z30 || z31) ? 0.0 : _boris_local_rounding(Float64(nodes[32]), u, D)),
    )
    z33 = (z8 && z32)
    r33 = _boris_sum_error(
        r8,
        r32,
        ((z8 || z32) ? 0.0 : _boris_local_rounding(Float64(nodes[33]), u, D)),
    )
    z34 = (z29 || z21)
    r34 =
        (z29 || z21) ? 0.0 :
        _boris_product_error(
            Float64(nodes[29]),
            Float64(nodes[21]),
            r29,
            r21,
            ((z29 || z21) ? 0.0 : _boris_local_rounding(Float64(nodes[34]), u, D)),
        )
    z35 = (z33 || z20)
    r35 =
        (z33 || z20) ? 0.0 :
        _boris_product_error(
            Float64(nodes[33]),
            Float64(nodes[20]),
            r33,
            r20,
            ((z33 || z20) ? 0.0 : _boris_local_rounding(Float64(nodes[35]), u, D)),
        )
    z36 = (z34 && z35)
    r36 = _boris_sum_error(
        r34,
        r35,
        ((z34 || z35) ? 0.0 : _boris_local_rounding(Float64(nodes[36]), u, D)),
    )
    z37 = (z6 && z36)
    r37 = _boris_sum_error(
        r6,
        r36,
        ((z6 || z36) ? 0.0 : _boris_local_rounding(Float64(nodes[37]), u, D)),
    )
    z38 = (z33 || z19)
    r38 =
        (z33 || z19) ? 0.0 :
        _boris_product_error(
            Float64(nodes[33]),
            Float64(nodes[19]),
            r33,
            r19,
            ((z33 || z19) ? 0.0 : _boris_local_rounding(Float64(nodes[38]), u, D)),
        )
    z39 = (z25 || z21)
    r39 =
        (z25 || z21) ? 0.0 :
        _boris_product_error(
            Float64(nodes[25]),
            Float64(nodes[21]),
            r25,
            r21,
            ((z25 || z21) ? 0.0 : _boris_local_rounding(Float64(nodes[39]), u, D)),
        )
    z40 = (z38 && z39)
    r40 = _boris_sum_error(
        r38,
        r39,
        ((z38 || z39) ? 0.0 : _boris_local_rounding(Float64(nodes[40]), u, D)),
    )
    z41 = (z7 && z40)
    r41 = _boris_sum_error(
        r7,
        r40,
        ((z7 || z40) ? 0.0 : _boris_local_rounding(Float64(nodes[41]), u, D)),
    )
    z42 = (z25 || z20)
    r42 =
        (z25 || z20) ? 0.0 :
        _boris_product_error(
            Float64(nodes[25]),
            Float64(nodes[20]),
            r25,
            r20,
            ((z25 || z20) ? 0.0 : _boris_local_rounding(Float64(nodes[42]), u, D)),
        )
    z43 = (z29 || z19)
    r43 =
        (z29 || z19) ? 0.0 :
        _boris_product_error(
            Float64(nodes[29]),
            Float64(nodes[19]),
            r29,
            r19,
            ((z29 || z19) ? 0.0 : _boris_local_rounding(Float64(nodes[43]), u, D)),
        )
    z44 = (z42 && z43)
    r44 = _boris_sum_error(
        r42,
        r43,
        ((z42 || z43) ? 0.0 : _boris_local_rounding(Float64(nodes[44]), u, D)),
    )
    z45 = (z8 && z44)
    r45 = _boris_sum_error(
        r8,
        r44,
        ((z8 || z44) ? 0.0 : _boris_local_rounding(Float64(nodes[45]), u, D)),
    )
    z46 = (z37 && z3)
    r46 = _boris_sum_error(
        r37,
        r3,
        ((z37 || z3) ? 0.0 : _boris_local_rounding(Float64(nodes[46]), u, D)),
    )
    z47 = (z41 && z4)
    r47 = _boris_sum_error(
        r41,
        r4,
        ((z41 || z4) ? 0.0 : _boris_local_rounding(Float64(nodes[47]), u, D)),
    )
    z48 = (z45 && z5)
    r48 = _boris_sum_error(
        r45,
        r5,
        ((z45 || z5) ? 0.0 : _boris_local_rounding(Float64(nodes[48]), u, D)),
    )
    return (
        r1,
        r2,
        r3,
        r4,
        r5,
        r6,
        r7,
        r8,
        r9,
        r10,
        r11,
        r12,
        r13,
        r14,
        r15,
        r16,
        r17,
        r18,
        r19,
        r20,
        r21,
        r22,
        r23,
        r24,
        r25,
        r26,
        r27,
        r28,
        r29,
        r30,
        r31,
        r32,
        r33,
        r34,
        r35,
        r36,
        r37,
        r38,
        r39,
        r40,
        r41,
        r42,
        r43,
        r44,
        r45,
        r46,
        r47,
        r48,
    )
end

@inline function _boris_accuracy_scalar(
    out,
    parts,
    nodes::NTuple{48,T},
    v,
    E,
    B,
    factor,
    dt::T,
) where {T}
    _boris_accuracy_mode() || return false
    all(isfinite, nodes) || return false
    eligible, coefficient, truezero = _boris_coefficient(factor, dt, parts[1])
    eligible || return false
    radii = _boris_scalar_radii(nodes, v, E, B, coefficient, truezero)
    all(x -> isfinite(x) && x >= 0, radii) || return false
    lower = max(
        maximum(abs, Float64.(v)),
        _boris_lower_difference(abs(Float64(nodes[3])), radii[3]),
        _boris_lower_difference(abs(Float64(nodes[4])), radii[4]),
        _boris_lower_difference(abs(Float64(nodes[5])), radii[5]),
    )
    e = iszero(lower) ? exponent(floatmin(T)) : max(exponent(lower), exponent(floatmin(T)))
    tolerance = ldexp(1.0, e - precision(T) + 3)
    margin =
        _boris_upper_add(abs(Float64(out[1])), radii[46]) <= Float64(floatmax(T)) &&
        _boris_upper_add(abs(Float64(out[2])), radii[47]) <= Float64(floatmax(T)) &&
        _boris_upper_add(abs(Float64(out[3])), radii[48]) <= Float64(floatmax(T))
    return margin && max(radii[46], radii[47], radii[48]) <= tolerance
end

@inline function _boris_batch_sum(a::NTuple{2,Float64}, b::NTuple{2,Float64}, u, D)
    P = _boris_upper_add(a[1], b[1])
    R = (iszero(a[1]) || iszero(b[1])) ? 0.0 : _boris_upper_add(_boris_upper_mul(u, P), D)
    return _boris_upper_add(P, R), _boris_sum_error(a[2], b[2], R)
end
@inline function _boris_batch_product(a::NTuple{2,Float64}, b::NTuple{2,Float64}, u, D)
    P = _boris_upper_mul(a[1], b[1])
    R = (iszero(a[1]) || iszero(b[1])) ? 0.0 : _boris_upper_add(_boris_upper_mul(u, P), D)
    return _boris_upper_add(P, R), _boris_product_error(a[1], b[1], a[2], b[2], R)
end

# Fixed Euler magnitude/error pairs for u ± H(E + u×B). Maxima bound each
# component; subtraction has the same absolute bound as addition. The native
# coefficient is hq, and coefficient bounds its relative error against H.
@inline function _euler_accuracy_bounds(hq::T, umax, emax, bmax, coefficient) where {T}
    u, D = ldexp(1.0, -precision(T)), Float64(nextfloat(zero(T)))
    velocity = (Float64(umax), 0.0)
    electric = (Float64(emax), 0.0)
    magnetic = (Float64(bmax), 0.0)
    radius = _boris_upper_mul(coefficient, abs(Float64(hq)))
    product = _boris_batch_product(velocity, magnetic, u, D)
    cross = _boris_batch_sum(product, product, u, D)
    force = _boris_batch_sum(electric, cross, u, D)
    kick = _boris_batch_product((abs(Float64(hq)), radius), force, u, D)
    output = _boris_batch_sum(velocity, kick, u, D)
    return output, kick, radius
end

@inline function _euler_accuracy_finish(bound, lower, ::Type{T}) where {T}
    all(x -> isfinite(x) && x >= 0, bound) || return false
    isfinite(lower) && lower >= 0 || return false
    e = iszero(lower) ? exponent(floatmin(T)) : max(exponent(lower), exponent(floatmin(T)))
    tolerance = ldexp(1.0, e - precision(T) + 3)
    return _boris_upper_add(bound[1], bound[2]) <= Float64(floatmax(T)) &&
           bound[2] <= tolerance
end

@inline function _euler_accuracy_scalar(out, v, E, B, h::T, factor, hq::T) where {T}
    T <: _RangeGuardedFloat || return true
    _boris_accuracy_mode() || return false
    eligible, coefficient, _ = _boris_coefficient(factor, h, hq)
    eligible || return false
    V, A, C = map(x -> maximum(abs, Float64.(x)), (v, E, B))
    bound, _, radius = _euler_accuracy_bounds(hq, V, A, C, coefficient)
    all(isfinite, bound) || return false
    # Each term is a lower bound on the retained exact Euclidean term scale
    # max(norm(v), |H|norm(E), |H|norm(v)norm(B)); no cancelled force is used.
    lower_h = _boris_lower_difference(abs(Float64(hq)), radius)
    lower = max(
        V,
        _boris_lower_product(lower_h, A),
        _boris_lower_product(_boris_lower_product(lower_h, V), C),
    )
    return _euler_accuracy_finish((maximum(abs, Float64.(out)), bound[2]), lower, T)
end

@inline function _euler_batch_accuracy(h::T, factor, hq::T, umax, emax, bmax, smin) where {T}
    T <: _RangeGuardedFloat || return true
    _boris_accuracy_mode() || return false
    eligible, coefficient, _ = _boris_coefficient(factor, h, hq)
    eligible || return false
    bound, kick, _ = _euler_accuracy_bounds(hq, umax, emax, bmax, coefficient)
    # The Euler reduction chooses max|v| for a moving row. A resting row has
    # v=0 exactly and chooses max|RN(hq*E)|, so subtracting the kick error gives
    # a lower bound on |H|norm(E). ES chooses the max of these two quantities
    # per component. Inf means no nonzero native scale: use the quantum floor.
    lower = _boris_lower_difference(Float64(smin), kick[2])
    return _euler_accuracy_finish(bound, lower, T)
end

@inline function _boris_batch_reciprocal(d::NTuple{2,Float64}, u, D)
    return 2.0,
    _boris_upper_add(_boris_upper_mul(2.0, d[2]), _boris_upper_add(_boris_upper_mul(2.0, u), D))
end

# Fixed pairs bound native magnitudes and exact-stage errors, respectively.
@inline function _boris_batch_bounds(h::T, factor, dt::T, vmax, emax, bmax, coefficient) where {T}
    u, D = ldexp(1.0, -precision(T)), Float64(nextfloat(zero(T)))
    v = (Float64(vmax), 0.0)
    E = (Float64(emax), 0.0)
    B = (Float64(bmax), 0.0)
    product = _boris_native_factor(factor) * dt
    c1 = (abs(Float64(product)), _boris_upper_mul(coefficient, abs(Float64(product))))
    c2 = (abs(Float64(h)), _boris_upper_mul(coefficient, abs(Float64(h))))
    c3 = _boris_batch_product(c2, E, u, D)
    c4 = _boris_batch_product(c2, E, u, D)
    c5 = _boris_batch_product(c2, E, u, D)
    c6 = _boris_batch_sum(v, c3, u, D)
    c7 = _boris_batch_sum(v, c4, u, D)
    c8 = _boris_batch_sum(v, c5, u, D)
    c9 = _boris_batch_product(c2, B, u, D)
    c10 = _boris_batch_product(c2, B, u, D)
    c11 = _boris_batch_product(c2, B, u, D)
    c12 = _boris_batch_product(c9, c9, u, D)
    c13 = _boris_batch_product(c10, c10, u, D)
    c14 = _boris_batch_product(c11, c11, u, D)
    c15 = _boris_batch_sum(c12, c13, u, D)
    c16 = _boris_batch_sum(c15, c14, u, D)
    c17 = _boris_batch_sum((1.0, 0.0), c16, u, D)
    c18 = _boris_batch_reciprocal(c17, u, D)
    c19 = _boris_batch_product(c18, c9, u, D)
    c20 = _boris_batch_product(c18, c10, u, D)
    c21 = _boris_batch_product(c18, c11, u, D)
    c22 = _boris_batch_product(c7, c11, u, D)
    c23 = _boris_batch_product(c8, c10, u, D)
    c24 = _boris_batch_sum(c22, c23, u, D)
    c25 = _boris_batch_sum(c6, c24, u, D)
    c26 = _boris_batch_product(c8, c9, u, D)
    c27 = _boris_batch_product(c6, c11, u, D)
    c28 = _boris_batch_sum(c26, c27, u, D)
    c29 = _boris_batch_sum(c7, c28, u, D)
    c30 = _boris_batch_product(c6, c10, u, D)
    c31 = _boris_batch_product(c7, c9, u, D)
    c32 = _boris_batch_sum(c30, c31, u, D)
    c33 = _boris_batch_sum(c8, c32, u, D)
    c34 = _boris_batch_product(c29, c21, u, D)
    c35 = _boris_batch_product(c33, c20, u, D)
    c36 = _boris_batch_sum(c34, c35, u, D)
    c37 = _boris_batch_sum(c6, c36, u, D)
    c38 = _boris_batch_product(c33, c19, u, D)
    c39 = _boris_batch_product(c25, c21, u, D)
    c40 = _boris_batch_sum(c38, c39, u, D)
    c41 = _boris_batch_sum(c7, c40, u, D)
    c42 = _boris_batch_product(c25, c20, u, D)
    c43 = _boris_batch_product(c29, c19, u, D)
    c44 = _boris_batch_sum(c42, c43, u, D)
    c45 = _boris_batch_sum(c8, c44, u, D)
    c46 = _boris_batch_sum(c37, c3, u, D)
    c47 = _boris_batch_sum(c41, c4, u, D)
    c48 = _boris_batch_sum(c45, c5, u, D)
    return (
        c1,
        c2,
        c3,
        c4,
        c5,
        c6,
        c7,
        c8,
        c9,
        c10,
        c11,
        c12,
        c13,
        c14,
        c15,
        c16,
        c17,
        c18,
        c19,
        c20,
        c21,
        c22,
        c23,
        c24,
        c25,
        c26,
        c27,
        c28,
        c29,
        c30,
        c31,
        c32,
        c33,
        c34,
        c35,
        c36,
        c37,
        c38,
        c39,
        c40,
        c41,
        c42,
        c43,
        c44,
        c45,
        c46,
        c47,
        c48,
    )
end

@inline function _boris_batch_accuracy(
    h::T,
    factor,
    dt::T,
    vmax::T,
    emax::T,
    bmax::T,
    smin::T,
) where {T}
    T <: _RangeGuardedFloat || return true
    _boris_accuracy_mode() || return false
    eligible, coefficient, _ = _boris_coefficient(factor, dt, h)
    eligible || return false
    bounds = _boris_batch_bounds(h, factor, dt, vmax, emax, bmax, coefficient)
    maximum_value = Float64(floatmax(T))
    all(
        c -> isfinite(c[1]) && isfinite(c[2]) && c[1] >= 0 && c[2] >= 0 && c[1] <= maximum_value,
        bounds,
    ) || return false
    Ek = max(bounds[3][2], bounds[4][2], bounds[5][2])
    # Inf means every native scale is zero (or there is no particle). Each
    # such particle has v=0 and RN(hE)=0, so exact |aE_i|<=(1+kappa)D/2.
    # Exact Cayley orthogonality gives |out|<4D when 3(1+kappa)^2<16.
    one_plus = _boris_upper_add(1.0, coefficient)
    _boris_upper_mul(3.0, _boris_upper_mul(one_plus, one_plus)) < 16.0 || return false
    if !isinf(smin)
        isfinite(smin) || return false
        lower = _boris_lower_difference(Float64(smin), Ek)
        e = iszero(lower) ? exponent(floatmin(T)) : max(exponent(lower), exponent(floatmin(T)))
        tolerance = ldexp(1.0, e - precision(T) + 3)
        max(bounds[46][2], bounds[47][2], bounds[48][2]) <= tolerance || return false
    end
    return all(
        c -> _boris_upper_add(c[1], c[2]) <= maximum_value,
        (bounds[46], bounds[47], bounds[48]),
    )
end

# The old predicate remains a range predicate. This owner alone authorizes
# production native batch bypass after independently checking full accuracy.
@inline function _boris_batch_admission(h, factor, dt, vmax, emax, bmax, smin)
    range, vbound =
        _boris_batch_certificate(h, _boris_native_factor(factor), dt, vmax, emax, bmax, smin)
    return range && _boris_batch_accuracy(h, factor, dt, vmax, emax, bmax, smin), vbound
end
@inline function _boris_batch_admission(ps::ParticleSet, h, factor, dt, vmax, emax, bmax, smin)
    admitted, vbound = _boris_batch_admission(h, factor, dt, vmax, emax, bmax, smin)
    # Repeated-index storage can change later velocity samples within one component.
    stable = admitted && all(_particle_batch_injective, ps.v) && !_particle_components_mightalias(ps)
    return stable, vbound
end
