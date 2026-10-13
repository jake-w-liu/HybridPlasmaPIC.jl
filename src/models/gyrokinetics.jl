# gyrokinetics.jl -- guiding-centre and gyro-averaged particle dynamics.

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross3(a, b) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline _scale3(s, a) = (s * a[1], s * a[2], s * a[3])
@inline _add3(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])

function _require_finite_vector3(name::AbstractString, a)
    length(a) == 3 || throw(DimensionMismatch("$name must have exactly three components"))
    @inbounds for i = 1:3
        a[i] isa Real || throw(ArgumentError("$name[$i] must be real"))
        isfinite(a[i]) || throw(ArgumentError("$name[$i] must be finite"))
    end
    return nothing
end

function _convert_finite_vector3(name::AbstractString, a, ::Type{T}) where {T<:AbstractFloat}
    _require_finite_vector3(name, a)
    return ntuple(3) do i
        value = T(a[i])
        isfinite(value) || throw(ArgumentError("$name[$i] must be finite and representable as $T"))
        return value
    end
end

function _require_finite_scalar(name::AbstractString, x)
    x isa Real || throw(ArgumentError("$name must be real"))
    isfinite(x) || throw(ArgumentError("$name must be finite"))
    return x
end

function _require_positive_scalar(name::AbstractString, x)
    _require_finite_scalar(name, x)
    x > zero(x) || throw(ArgumentError("$name must be positive"))
    return x
end

function _require_nonnegative_scalar(name::AbstractString, x)
    _require_finite_scalar(name, x)
    x >= zero(x) || throw(ArgumentError("$name must be nonnegative"))
    return x
end

@inline function _gyro_scaled_vector(a)
    values = promote(float(a[1]), float(a[2]), float(a[3]))
    scale = max(abs(values[1]), abs(values[2]), abs(values[3]))
    if iszero(scale)
        z = zero(scale)
        return scale, (z, z, z)
    end
    return scale, (values[1] / scale, values[2] / scale, values[3] / scale)
end

@noinline function _gyro_exact_ratio(
    numerators::Tuple,
    denominators::Tuple,
    ::Type{T},
) where {T<:AbstractFloat}
    exact = _exact_rational(one(T))
    for x in numerators
        exact *= _exact_rational(x)
    end
    for x in denominators
        exact /= _exact_rational(x)
    end
    return _from_exact(T, exact)
end

@inline function _gyro_ratio_needs_exact(result::T) where {T<:Union{Float16,Float32,Float64}}
    return iszero(result) ||
           !isfinite(result) ||
           issubnormal(result) ||
           abs(result) >= floatmax(T) / T(2)
end

@inline _gyro_ratio_needs_exact(result) = iszero(result) || !isfinite(result)

# Evaluate a product ratio through binary mantissas so finite final values are
# retained even when a direct numerator or denominator would overflow or underflow.
@inline function _gyro_scaled_ratio(numerators::Tuple, denominators::Tuple)
    promoted = promote_type(map(typeof, numerators)..., map(typeof, denominators)...)
    T = float(promoted)
    converted_numerators = map(T, numerators)
    converted_denominators = map(T, denominators)
    mantissa = one(T)
    exponent = 0
    for x in converted_numerators
        iszero(x) && return zero(T)
        m, e = frexp(x)
        mantissa *= m
        exponent += e
    end
    for x in converted_denominators
        m, e = frexp(x)
        mantissa /= m
        exponent -= e
    end
    result = ldexp(mantissa, exponent)
    if _gyro_ratio_needs_exact(result) &&
       all(isfinite, converted_numerators) &&
       all(isfinite, converted_denominators)
        return _gyro_exact_ratio(converted_numerators, converted_denominators, T)
    end
    return result
end

@inline function _gyro_compensated_add(total, correction, term)
    updated = total + term
    correction += abs(total) >= abs(term) ? (total - updated) + term : (term - updated) + total
    return updated, correction
end

@inline function _gyro_compensated_sum3(a, b, c)
    values = promote(a, b, c)
    total = zero(values[1])
    correction = zero(values[1])
    total, correction = _gyro_compensated_add(total, correction, values[1])
    total, correction = _gyro_compensated_add(total, correction, values[2])
    total, correction = _gyro_compensated_add(total, correction, values[3])
    return total + correction
end

# Combine three finite drift terms without an order-dependent overflow. Large
# opposite-signed terms are combined first; otherwise a common scale is used.
@inline function _gyro_scaled_sum3(a, b, c)
    x, y, z = promote(float(a), float(b), float(c))
    if abs(x) < abs(y)
        x, y = y, x
    end
    if abs(y) < abs(z)
        y, z = z, y
    end
    if abs(x) < abs(y)
        x, y = y, x
    end

    if !iszero(y) && signbit(x) != signbit(y)
        result = _gyro_compensated_sum3(x, y, z)
        isfinite(result) && return result
    elseif !iszero(z) && signbit(x) != signbit(z)
        result = _gyro_compensated_sum3(x, z, y)
        isfinite(result) && return result
    end

    scale = max(abs(x), abs(y), abs(z))
    iszero(scale) && return zero(scale)
    scaled_sum = _gyro_compensated_sum3(x / scale, y / scale, z / scale)
    return _gyro_scaled_ratio((scale, scaled_sum), (one(scale),))
end

@inline function _gyro_two_sum(a, b)
    total = a + b
    recovered_b = total - a
    error = (a - (total - recovered_b)) + (b - recovered_b)
    return total, error
end

# An FMA product residual fits p bits but can still underflow. Each p-bit
# operand is an integer multiple of 2^(exponent(x) - p + 1); require the
# product's smallest possible bit to remain representable. This also guards
# a zero quotient residual against silently losing a nonzero remainder.
@inline function _gyro_product_residual_is_exact(
    a::T,
    b::T,
    residual::T,
) where {T<:Union{Float16,Float32,Float64}}
    isfinite(residual) || return false
    (iszero(a) || iszero(b)) && return true
    return exponent(abs(a)) + exponent(abs(b)) - 2(precision(T) - 1) >= exponent(nextfloat(zero(T)))
end

@inline _gyro_product_residual_is_exact(a, b, residual) = false

@inline _gyro_product_is_exact(a, b, residual) =
    iszero(residual) && _gyro_product_residual_is_exact(a, b, residual)

@inline _gyro_sum_is_exact(a, b) = iszero(_gyro_two_sum(a, b)[2])

# Keep a product ratio as a compensated mantissa and a separate binary
# exponent. The certificate is deliberately stronger than compensated accuracy:
# true means the high-plus-low expansion is exact, without discarded residuals.
@inline function _gyro_ratio_parts(numerators::NTuple{N,T}, denominators::NTuple{D,T}) where {N,D,T}
    any(iszero, numerators) && return (zero(T), zero(T), 0, true)
    hi, lo, exponent = one(T), zero(T), 0
    exact = true
    # A hardware float's frexp mantissa is exact. Start from the first
    # factor directly; generic floats retain their working-precision multiply.
    if N > 0 && T <: Union{Float16,Float32,Float64}
        hi, exponent = frexp(first(numerators))
        factors = Base.tail(numerators)
    else
        factors = numerators
    end
    for x in factors
        mantissa, power = frexp(x)
        product = hi * mantissa
        residual = fma(hi, mantissa, -product)
        tail = lo * mantissa
        exact =
            exact &&
            _gyro_product_residual_is_exact(hi, mantissa, residual) &&
            _gyro_product_is_exact(lo, mantissa, fma(lo, mantissa, -tail)) &&
            _gyro_sum_is_exact(residual, tail)
        error = residual + tail
        hi, lo = _gyro_two_sum(product, error)
        exact = exact && isfinite(hi) && isfinite(lo)
        exponent += power
    end
    for x in denominators
        mantissa, power = frexp(x)
        quotient = hi / mantissa
        residual = fma(-quotient, mantissa, hi)
        # A nonzero low component would require proving its division too;
        # conservatively leave that quotient uncertified.
        exact = exact && iszero(lo) && _gyro_product_is_exact(quotient, mantissa, residual)
        error = (residual + lo) / mantissa
        hi, lo = _gyro_two_sum(quotient, error)
        exponent -= power
    end
    return hi, lo, exponent, exact
end

@inline function _gyro_sum_parts(terms, factor_operations = 7)
    T = typeof(terms[1][1])
    common_exponent = typemin(Int)
    all_exact = true
    for (hi, lo, exponent, exact) in terms
        all_exact &= exact
        if !iszero(hi) || !iszero(lo)
            common_exponent = max(common_exponent, exponent)
        end
    end
    common_exponent == typemin(Int) && return (zero(T), zero(T), 0, all_exact), !all_exact
    hi, lo, magnitude = zero(T), zero(T), zero(T)
    lost_component = false
    for (term_hi, term_lo, exponent, _) in terms
        iszero(term_hi) && iszero(term_lo) && continue
        shift = exponent - common_exponent
        scaled_hi = ldexp(term_hi, shift)
        scaled_lo = ldexp(term_lo, shift)
        lost_component |=
            (!iszero(term_hi) && iszero(scaled_hi)) || (!iszero(term_lo) && iszero(scaled_lo))
        sum_hi, error_hi = _gyro_two_sum(hi, scaled_hi)
        sum_lo, error_lo = _gyro_two_sum(lo, scaled_lo)
        # Both aligned components and each compression of the error terms
        # must be exact; TwoSum itself retains the residual of each addition.
        all_exact =
            all_exact &&
            ldexp(scaled_hi, -shift) == term_hi &&
            ldexp(scaled_lo, -shift) == term_lo &&
            _gyro_sum_is_exact(error_hi, sum_lo)
        hi, lo = _gyro_two_sum(sum_hi, error_hi + sum_lo)
        all_exact = all_exact && _gyro_sum_is_exact(lo, error_lo)
        hi, lo = _gyro_two_sum(hi, lo + error_lo)
        all_exact = all_exact && isfinite(hi) && isfinite(lo)
        magnitude += abs(scaled_hi) + abs(scaled_lo)
    end
    combined = hi + lo
    # The caller bounds the factor operations in each term. Allow a conservative
    # two rounding contributions per operation plus summation, and use exact
    # arithmetic when cancellation could expose that second-order roundoff.
    cancellation_limit = T((2factor_operations + 2) * length(terms)) * eps(T) * magnitude
    needs_exact = lost_component || abs(combined) <= cancellation_limit
    needs_exact &= !(all_exact && iszero(combined))
    return (hi, lo, common_exponent, all_exact), needs_exact
end

@inline function _gyro_sum_scaled_parts(terms, factor_operations = 7)
    (hi, lo, exponent, exact), needs_exact = _gyro_sum_parts(terms, factor_operations)
    combined = hi + lo
    result = ldexp(combined, exponent)
    needs_exact |= _gyro_ratio_needs_exact(result) && !(exact && iszero(combined))
    return result, needs_exact
end

@inline function _gyro_multiply_parts(left, right)
    left_hi, left_lo, left_exponent, left_exact = left
    right_hi, right_lo, right_exponent, right_exact = right
    if (left_exact && iszero(left_hi) && iszero(left_lo)) ||
       (right_exact && iszero(right_hi) && iszero(right_lo))
        return zero(left_hi), zero(left_hi), 0, true
    end
    product = left_hi * right_hi
    residual = fma(left_hi, right_hi, -product)
    high_low = left_hi * right_lo
    low_high = left_lo * right_hi
    low_low = left_lo * right_lo
    correction = residual + high_low
    next_correction = correction + low_high
    exact =
        left_exact &&
        right_exact &&
        _gyro_product_residual_is_exact(left_hi, right_hi, residual) &&
        _gyro_product_is_exact(left_hi, right_lo, fma(left_hi, right_lo, -high_low)) &&
        _gyro_product_is_exact(left_lo, right_hi, fma(left_lo, right_hi, -low_high)) &&
        _gyro_product_is_exact(left_lo, right_lo, fma(left_lo, right_lo, -low_low)) &&
        _gyro_sum_is_exact(residual, high_low) &&
        _gyro_sum_is_exact(correction, low_high) &&
        _gyro_sum_is_exact(next_correction, low_low)
    error = next_correction + low_low
    hi, lo = _gyro_two_sum(product, error)
    exact = exact && isfinite(hi) && isfinite(lo)
    return hi, lo, left_exponent + right_exponent, exact
end

# Keep both products of a cross component scaled until its caller chooses the
# materialization boundary (a drift velocity or a complete position update).
@inline function _gyro_cross_parts(a, b, coefficient, i)
    j, k = mod1(i + 1, 3), mod1(i + 2, 3)
    left = _gyro_ratio_parts((a[j], b[k]), ())
    right = _gyro_ratio_parts((a[k], b[j]), ())
    T = typeof(left[1])
    # Each normalized two-factor product has an exact high/low expansion in
    # these formats: its residual quantum 2^(-2precision(T)) exceeds minsub.
    # Equal expansions therefore prove a zero cross component before an
    # inexact coefficient can obscure this common parallel-vector cancellation.
    if T <: Union{Float16,Float32,Float64} &&
       left[1] == right[1] &&
       left[2] == right[2] &&
       left[3] == right[3]
        z = (zero(T), zero(T), 0, true)
        return z, z
    end
    negative_right = (-right[1], -right[2], right[3], right[4])
    return _gyro_multiply_parts(left, coefficient),
    _gyro_multiply_parts(negative_right, coefficient)
end

@noinline function _gyro_exact_position_component(
    X,
    vpar,
    μ,
    q,
    m,
    dt,
    E,
    B,
    gradB,
    κ,
    Bscale,
    Bscaled_norm,
    i,
)
    exact(x) = _exact_rational(x)
    j, k = mod1(i + 1, 3), mod1(i + 2, 3)
    cross(a, b) = exact(a[j]) * exact(b[k]) - exact(a[k]) * exact(b[j])
    magnitude = exact(Bscale) * exact(Bscaled_norm)
    perpendicular =
        (
            cross(E, B) +
            exact(μ) / exact(q) * cross(B, gradB) +
            exact(m) * exact(vpar)^2 / exact(q) * cross(B, κ)
        ) / magnitude^2
    velocity = exact(vpar) * exact(B[i]) / magnitude + perpendicular
    return _from_exact(typeof(vpar), exact(X[i]) + exact(dt) * velocity)
end

@inline function _gyro_position(X, vpar, μ, q, m, dt, E, B, gradB, κ, Bscale, Bscaled_norm)
    iszero(dt) && return X
    streaming = _gyro_ratio_parts((dt, vpar), (Bscale, Bscaled_norm))
    displacement = _gyro_ratio_parts((dt,), (q, Bscale, Bscale, Bscaled_norm, Bscaled_norm))
    electric = _gyro_ratio_parts((q,), ())
    gradient = _gyro_ratio_parts((μ,), ())
    curvature = _gyro_ratio_parts((m, vpar, vpar), ())
    return ntuple(3) do i
        initial = _gyro_ratio_parts((X[i],), ())
        parallel = _gyro_multiply_parts(_gyro_ratio_parts((B[i],), ()), streaming)
        exb = _gyro_cross_parts(E, B, electric, i)
        grad = _gyro_cross_parts(B, gradB, gradient, i)
        curv = _gyro_cross_parts(B, κ, curvature, i)
        # Retain the common displacement coefficient outside the numerator:
        # an exact drift balance can then be certified before inexact scaling.
        perpendicular, uncertain_perpendicular = _gyro_sum_parts((exb..., grad..., curv...))
        perpendicular_step = _gyro_multiply_parts(perpendicular, displacement)
        # A composed term has at most six numerator operations, six shared
        # ratio factors, and the final multiplication; allow fourteen in total.
        result, needs_exact = _gyro_sum_scaled_parts((initial, parallel, perpendicular_step), 14)
        return needs_exact || uncertain_perpendicular ?
               _gyro_exact_position_component(
            X,
            vpar,
            μ,
            q,
            m,
            dt,
            E,
            B,
            gradB,
            κ,
            Bscale,
            Bscaled_norm,
            i,
        ) : result
    end
end

@noinline function _gyro_exact_parallel_velocity(
    vpar::T,
    dt,
    q,
    E,
    B,
    Bscale,
    Bscaled_norm,
    μ,
    gradpar_B,
    m,
) where {T}
    exact(x) = _exact_rational(x)
    field_dot = sum(exact(E[i]) * exact(B[i]) for i = 1:3)
    Epar = field_dot / (exact(Bscale) * exact(Bscaled_norm))
    return _from_exact(
        T,
        exact(vpar) + exact(dt) / exact(m) * (exact(q) * Epar - exact(μ) * exact(gradpar_B)),
    )
end

@inline function _gyro_parallel_velocity(vpar, dt, q, E, B, Bscale, Bscaled_norm, μ, gradpar_B, m)
    iszero(dt) && return vpar
    electric_scale = _gyro_ratio_parts((dt, q), (m, Bscale, Bscaled_norm))
    electric =
        ntuple(i -> _gyro_multiply_parts(_gyro_ratio_parts((E[i], B[i]), ()), electric_scale), 3)
    mirror = _gyro_ratio_parts((-dt, μ, gradpar_B), (m,))
    initial = _gyro_ratio_parts((vpar,), ())
    result, needs_exact = _gyro_sum_scaled_parts((electric..., mirror, initial))
    if needs_exact
        return _gyro_exact_parallel_velocity(
            vpar,
            dt,
            q,
            E,
            B,
            Bscale,
            Bscaled_norm,
            μ,
            gradpar_B,
            m,
        )
    end
    return result
end

@inline function _gyro_nonzero_field(B)
    _require_finite_vector3("B", B)
    scale, scaled = _gyro_scaled_vector(B)
    scale > zero(scale) ||
        throw(ArgumentError("magnetic field magnitude must be finite and positive"))
    scaled_norm = hypot(scaled[1], scaled[2], scaled[3])
    direction = ntuple(i -> scaled[i] / scaled_norm, 3)
    return scale, scaled_norm, direction
end

function _normalize3(a)
    return _gyro_nonzero_field(a)[3]
end

"""
    exb_drift(E, B) -> NTuple{3}

The `E×B` drift `v_E = (E × B)/B^2`. `E` and `B` are 3-component vectors.
"""
function exb_drift(E, B)
    _require_finite_vector3("E", E)
    Bscale, Bscaled_norm, b = _gyro_nonzero_field(B)
    Escale, Escaled = _gyro_scaled_vector(E)
    direction = _cross3(Escaled, b)
    drift = ntuple(i -> _gyro_scaled_ratio((Escale, direction[i]), (Bscale, Bscaled_norm)), 3)
    _require_finite_vector3("E×B drift", drift)
    return drift
end

"""
    gradb_drift(vperp, q, m, B, gradB) -> NTuple{3}

Grad-B drift `v_∇B = (m v⊥^2)/(2 q B^3) (B × ∇B)`, where `gradB = ∇|B|`.
"""
function gradb_drift(vperp, q, m, B, gradB)
    _require_nonnegative_scalar("vperp", vperp)
    _require_finite_scalar("q", q)
    q != zero(q) || throw(ArgumentError("charge q must be nonzero"))
    _require_positive_scalar("m", m)
    _require_finite_vector3("gradB", gradB)
    Bscale, Bscaled_norm, b = _gyro_nonzero_field(B)
    grad_scale, grad_scaled = _gyro_scaled_vector(gradB)
    direction = _cross3(b, grad_scaled)
    drift = ntuple(
        i -> _gyro_scaled_ratio(
            (m, vperp, vperp, grad_scale, direction[i]),
            (2, q, Bscale, Bscale, Bscaled_norm, Bscaled_norm),
        ),
        3,
    )
    _require_finite_vector3("grad-B drift", drift)
    return drift
end

"""
    curvature_drift(vpar, q, m, B, κ) -> NTuple{3}

Curvature drift `v_κ = (m v∥^2)/(q B^2) (B × κ)`.
"""
function curvature_drift(vpar, q, m, B, κ)
    _require_finite_scalar("vpar", vpar)
    _require_finite_scalar("q", q)
    q != zero(q) || throw(ArgumentError("charge q must be nonzero"))
    _require_positive_scalar("m", m)
    _require_finite_vector3("κ", κ)
    Bscale, Bscaled_norm, b = _gyro_nonzero_field(B)
    κscale, κscaled = _gyro_scaled_vector(κ)
    direction = _cross3(b, κscaled)
    drift = ntuple(
        i -> _gyro_scaled_ratio((m, vpar, vpar, κscale, direction[i]), (q, Bscale, Bscaled_norm)),
        3,
    )
    _require_finite_vector3("curvature drift", drift)
    return drift
end

"""
    drift_velocity(; vpar, vperp, q, m, E, B, gradB, κ) -> NTuple{3}

Total perpendicular guiding-centre drift `v_E + v_∇B + v_κ`.
"""
function drift_velocity(; vpar, vperp, q, m, E, B, gradB, κ)
    vE = exb_drift(E, B)
    vgradB = gradb_drift(vperp, q, m, B, gradB)
    vcurvature = curvature_drift(vpar, q, m, B, κ)
    drift = ntuple(i -> _gyro_scaled_sum3(vE[i], vgradB[i], vcurvature[i]), 3)
    _require_finite_vector3("total drift", drift)
    return drift
end

"""
    GuidingCentre(X, vpar, μ, q, m)

Guiding-centre state with position `X`, parallel velocity `vpar`, magnetic
moment `μ = m v⊥^2/(2B)`, charge `q`, and mass `m`.
"""
mutable struct GuidingCentre{T<:AbstractFloat}
    X::NTuple{3,T}
    vpar::T
    μ::T
    q::T
    m::T

    function GuidingCentre{T}(X::NTuple{3,T}, vpar::T, μ::T, q::T, m::T) where {T<:AbstractFloat}
        _require_finite_vector3("X", X)
        _require_finite_scalar("vpar", vpar)
        _require_nonnegative_scalar("μ", μ)
        _require_finite_scalar("q", q)
        q != zero(T) || throw(ArgumentError("charge q must be nonzero"))
        _require_positive_scalar("m", m)
        return new{T}(X, vpar, μ, q, m)
    end
end

function GuidingCentre(X::NTuple{3,S}, vpar, μ, q, m) where {S<:Real}
    T = S <: AbstractFloat ? S : typeof(float(zero(S)))
    XT = ntuple(i -> _require_finite_real("X[$i]", X[i], T), 3)
    return GuidingCentre{T}(XT, T(vpar), T(μ), T(q), T(m))
end

function _validate_guiding_centre(gc::GuidingCentre)
    _require_finite_vector3("X", gc.X)
    _require_finite_scalar("vpar", gc.vpar)
    _require_nonnegative_scalar("μ", gc.μ)
    _require_finite_scalar("q", gc.q)
    gc.q != zero(gc.q) || throw(ArgumentError("charge q must be nonzero"))
    _require_positive_scalar("m", gc.m)
    return nothing
end

"""
    push_guiding_centre!(gc; dt, E, B, gradB, κ, gradpar_B) -> gc

Advance the guiding centre one local-field step. The position advances by
perpendicular drifts plus parallel streaming `v∥ b`; `v∥` advances by
`m dv∥/dt = q E∥ - μ ∂∥B`. `μ` is conserved.
"""
function push_guiding_centre!(gc::GuidingCentre{T}; dt, E, B, gradB, κ, gradpar_B) where {T}
    _validate_guiding_centre(gc)
    dtT = _require_finite_real("dt", dt, T)
    (iszero(dt) || !iszero(dtT)) ||
        throw(ArgumentError("dt must be finite and representable as $T"))
    ET = _convert_finite_vector3("E", E, T)
    BT = _convert_finite_vector3("B", B, T)
    gradBT = _convert_finite_vector3("gradB", gradB, T)
    κT = _convert_finite_vector3("κ", κ, T)
    gradparT = _require_finite_real("gradpar_B", gradpar_B, T)
    Bscale, Bscaled_norm, _ = _gyro_nonzero_field(BT)
    # Use μ directly and retain X plus all displacement terms until the final
    # rounding; individual drift velocities need not fit in the state type.
    new_X = _gyro_position(
        gc.X,
        gc.vpar,
        gc.μ,
        gc.q,
        gc.m,
        dtT,
        ET,
        BT,
        gradBT,
        κT,
        Bscale,
        Bscaled_norm,
    )
    new_vpar = _gyro_parallel_velocity(
        gc.vpar,
        dtT,
        gc.q,
        ET,
        BT,
        Bscale,
        Bscaled_norm,
        gc.μ,
        gradparT,
        gc.m,
    )
    _require_finite_vector3("updated X", new_X)
    _require_finite_scalar("updated vpar", new_vpar)

    # Commit only after all direct, derived, and updated values are valid.
    gc.X = new_X
    gc.vpar = new_vpar
    return gc
end

function _perp_basis(b::NTuple{3,T}) where {T}
    ref = abs(b[1]) < T(0.9) ? (one(T), zero(T), zero(T)) : (zero(T), one(T), zero(T))
    e1 = _normalize3(_sub3(ref, _scale3(_dot3(ref, b), b)))
    e2 = _cross3(b, e1)
    return e1, e2
end

"""
    gyroaverage(f, X, ρ, B; n=16) -> value

Gyro-average of scalar callable `f(x)` over a ring of radius `ρ` centred at `X`
in the plane perpendicular to `B`, sampled at `n` points.
"""
function gyroaverage(f, X::NTuple{3,T}, ρ, B; n::Integer = 16) where {T<:AbstractFloat}
    nring = _require_positive_intlike("n", n)
    nring >= 3 || throw(ArgumentError("need n >= 3 ring points"))
    XT = _convert_finite_vector3("X", X, T)
    BT = _convert_finite_vector3("B", B, T)
    ρT = _require_finite_nonnegative_real("ρ", ρ, T)
    b = _normalize3(BT)
    e1, e2 = _perp_basis(b)
    scale = zero(T)
    total = zero(T)
    correction = zero(T)
    @inbounds for k = 1:nring
        φ = 2 * T(π) * (k - 1) / nring
        pt = _add3(XT, _add3(_scale3(ρT * cos(φ), e1), _scale3(ρT * sin(φ), e2)))
        raw_sample = f(pt)
        raw_sample isa Real || throw(ArgumentError("f must return real scalar values"))
        isfinite(raw_sample) || throw(ArgumentError("f must return finite values"))
        sample = float(raw_sample) + zero(T)
        isfinite(sample) ||
            throw(ArgumentError("f return values must be representable as a finite float"))
        magnitude = abs(sample)
        if magnitude > scale
            if iszero(scale)
                total = zero(sample)
                correction = zero(sample)
            else
                factor = scale / magnitude
                scaled_total = total * factor
                scaled_correction = correction * factor
                (
                    !iszero(total) && iszero(scaled_total) ||
                    !iszero(correction) && iszero(scaled_correction)
                ) && throw(
                    ArgumentError("sample dynamic range is too wide for lossless normalisation"),
                )
                total = scaled_total
                correction = scaled_correction
            end
            scale = magnitude
        end
        if !iszero(scale)
            normalised = sample / scale
            !iszero(sample) &&
                iszero(normalised) &&
                throw(ArgumentError("sample dynamic range is too wide for lossless normalisation"))
            total, correction = _gyro_compensated_add(total, correction, normalised)
        else
            total += sample
        end
    end
    iszero(scale) && return total
    return _gyro_scaled_ratio((scale, total + correction), (nring,))
end
