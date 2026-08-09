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
    exact = Rational{BigInt}(one(T))
    for x in numerators
        exact *= Rational{BigInt}(x)
    end
    for x in denominators
        exact /= Rational{BigInt}(x)
    end
    return T(exact)
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

@inline function _gyro_compensated_sum4(a, b, c, d)
    values = promote(a, b, c, d)
    total = zero(values[1])
    correction = zero(values[1])
    total, correction = _gyro_compensated_add(total, correction, values[1])
    total, correction = _gyro_compensated_add(total, correction, values[2])
    total, correction = _gyro_compensated_add(total, correction, values[3])
    total, correction = _gyro_compensated_add(total, correction, values[4])
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

@inline function _gyro_product_parts(a, b)
    ma, ea = frexp(a)
    mb, eb = frexp(b)
    product = ma * mb
    error = fma(ma, mb, -product)
    return product, error, ea + eb
end

@inline function _gyro_scaled_product(a, b)
    (iszero(a) || iszero(b)) && return zero(a * b)
    product, error, exponent = _gyro_product_parts(a, b)
    return ldexp(product + error, exponent)
end

# Evaluate `a*b - c*d` from aligned exact product expansions, avoiding `Inf-Inf`
# and crossed-scale underflow when the mathematical result is representable.
@inline function _gyro_scaled_product_difference(a, b, c, d)
    aT, bT, cT, dT = promote(float(a), float(b), float(c), float(d))
    left_zero = iszero(aT) || iszero(bT)
    right_zero = iszero(cT) || iszero(dT)
    left_zero && right_zero && return zero(aT)
    left_zero && return -_gyro_scaled_product(cT, dT)
    right_zero && return _gyro_scaled_product(aT, bT)

    left, left_error, left_exponent = _gyro_product_parts(aT, bT)
    right, right_error, right_exponent = _gyro_product_parts(cT, dT)
    common_exponent = max(left_exponent, right_exponent)
    left_shift = left_exponent - common_exponent
    right_shift = right_exponent - common_exponent
    difference = _gyro_compensated_sum4(
        ldexp(left, left_shift),
        -ldexp(right, right_shift),
        ldexp(left_error, left_shift),
        -ldexp(right_error, right_shift),
    )
    return ldexp(difference, common_exponent)
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
    Bscale, Bscaled_norm, b = _gyro_nonzero_field(BT)
    vperp2 = _gyro_scaled_ratio((T(2), gc.μ, Bscale, Bscaled_norm), (gc.m,))
    isfinite(vperp2) && vperp2 >= zero(T) ||
        throw(ArgumentError("local perpendicular speed is not finite and representable as $T"))
    vperp = sqrt(vperp2)
    vd = drift_velocity(;
        vpar = gc.vpar,
        vperp,
        q = gc.q,
        m = gc.m,
        E = ET,
        B = BT,
        gradB = gradBT,
        κ = κT,
    )
    velocity = ntuple(i -> _gyro_scaled_sum3(vd[i], gc.vpar * b[i], zero(T)), 3)
    step = ntuple(i -> _gyro_scaled_ratio((dtT, velocity[i]), (one(T),)), 3)
    new_X = _add3(gc.X, step)
    Escale, Escaled = _gyro_scaled_vector(ET)
    Epar = _gyro_scaled_ratio((Escale, _dot3(Escaled, b)), (one(T),))
    force = _gyro_scaled_product_difference(gc.q, Epar, gc.μ, gradparT)
    delta_vpar = _gyro_scaled_ratio((dtT, force), (gc.m,))
    new_vpar = gc.vpar + delta_vpar
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
