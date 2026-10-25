# Canonical sampled-angle and complete Coulomb velocity arithmetic.
# Extracted from the local multispecies candidate so both public owners share
# one scattering law; pairing, rejection and publication stay with each caller.

@inline _coulomb_product_needs_complete(a::T, b::T) where {T} =
    !iszero(a) && !iszero(b) && abs(a * b) < floatmin(T)

@inline function _coulomb_increment_needs_complete(value::T, a::T, b::T, c::T) where {T}
    scale = max(abs(a), abs(b), abs(c))
    # Re-evaluate cancellation before it loses a leading bit. Waiting until
    # half the significand is lost permits large relative component errors.
    return !iszero(scale) && (abs(value) < floatmin(T) || abs(value) <= scale / 2)
end

# All inputs are finite positive rate factors from validated public inputs.
# Their binary fractions stay near unity while integer exponents retain a rate
# that would overflow or underflow if multiplied before taking its square root.
function _coulomb_scaled_variance(rate_factors, wmax::T, wbar, speed::T) where {T}
    numerator, power = frexp(wmax)
    for factor in rate_factors
        fraction, exponent = frexp(factor)
        numerator *= fraction
        power += exponent
    end
    average, average_power = frexp(wbar)
    velocity, velocity_power = frexp(speed)
    power -= average_power + 3velocity_power
    ratio = T(numerator / (average * velocity * velocity * velocity))
    return ratio, power
end

function _coulomb_scaled_tangent(ratio::T, power::Int, draw::T) where {T}
    isodd(power) && (ratio *= 2)
    fraction, draw_power = frexp(draw)
    return sqrt(ratio) * fraction, fld(power, 2) + draw_power
end

# Resolve the model's sharp variance-one boundary from the original finite
# factors. Only their bounded binary fractions become rationals; large input
# exponents cancel as integers instead of creating huge powers of two.
@noinline function _coulomb_large_angle_exact(rate_factors, wmax, wbar, speed)
    fraction, exponent = frexp(wmax)
    numerator = _exact_rational(fraction)
    power = BigInt(exponent)
    for factor in rate_factors
        fraction, exponent = frexp(factor)
        numerator *= _exact_rational(fraction)
        power += exponent
    end
    average, average_power = frexp(wbar)
    velocity, velocity_power = frexp(speed)
    denominator = _exact_rational(average) * _exact_rational(velocity)^3
    power -= BigInt(average_power) + 3BigInt(velocity_power)
    # Each positive fraction lies in [1/2,1). Thus the ratio lies strictly
    # between 2^(-N-1) and 16, bounding any required final binary shift.
    power > length(rate_factors) + 1 && return true
    power < -4 && return false
    p = Int(power)
    return p >= 0 ? numerator * (BigInt(1) << p) > denominator :
           numerator > denominator * (BigInt(1) << -p)
end

@inline function _coulomb_sample_increment(
    gx::T,
    gy::T,
    gz::T,
    gmag::T,
    gcdt::T,
    wmax::T,
    wbar,
    uf::T,
    twoπ::T,
    rng,
    rate_factors,
) where {T}
    # Takizuka-Abe ⟨δ²⟩, weight-corrected: × max(wᵢ,wⱼ)/w̄ (Higginson 2020 eq. 31).
    # The ordinary-range path avoids transcendental work in this hot loop.
    # Binary scaling retains the original factors when a product leaves the
    # normal range. A narrow exact comparison resolves the model boundary.
    speed = max(gmag, uf)
    weight_ratio = T(wmax / wbar)
    numerator = gcdt * weight_ratio
    speed2 = speed * speed
    denominator = speed2 * speed
    use_fast =
        isfinite(weight_ratio) &&
        weight_ratio >= floatmin(T) &&
        isfinite(numerator) &&
        numerator >= floatmin(T) &&
        isfinite(denominator) &&
        denominator >= floatmin(T) &&
        isfinite(gcdt) &&
        gcdt >= floatmin(T)
    ratio, rate_power = gcdt, 0
    if use_fast
        var = numerator / denominator
    else
        ratio, rate_power = _coulomb_scaled_variance(rate_factors, wmax, wbar, speed)
        var = ldexp(ratio, rate_power)
    end
    # At most eight rounded products/divisions contribute to the variance.
    # Sixteen eps covers their forward error for the supported hardware floats
    # and ordinary BigFloat precision. Coarser BigFloat always takes the exact
    # comparison, where a relative forward-error bound would be too wide.
    large_angle =
        (eps(T) >= one(T) / 16 || abs(var - one(T)) <= 16eps(T)) ?
        _coulomb_large_angle_exact(rate_factors, wmax, wbar, speed) : var > one(T)
    # Isotropic sampling does not use the tangent; reuse a scalar for that field.
    tangent, tangent_power, complete = gx, 0, false
    if large_angle                           # large-angle regime → isotropic scatter
        cosθ = 2 * rand(rng, T) - one(T)
        cosθ_minus_one = cosθ - one(T)
        sinθ = sqrt(max(zero(T), one(T) - cosθ * cosθ))
        identity = cosθ == one(T) && iszero(sinθ)
    else
        draw = randn(rng, T)
        δ = sqrt(var) * draw
        isfinite(δ) ||
            throw(ArgumentError("collide_coulomb!: sampled scattering angle must be finite"))
        identity = iszero(draw)
        tangent = δ
        if !iszero(draw) && (!use_fast || var < floatmin(T) || abs(δ) < floatmin(T))
            # A rounded variance or tangent can vanish while the final
            # physical kick remains representable. Retain its binary scale.
            if use_fast
                ratio, rate_power = _coulomb_scaled_variance(rate_factors, wmax, wbar, speed)
            end
            tangent, tangent_power = _coulomb_scaled_tangent(ratio, rate_power, draw)
            complete = true
        elseif !iszero(δ)
            # Both tangent and reciprocal-tangent forms must retain their
            # squares; the lost term may survive after multiplication by g.
            small = abs(δ) <= one(T) ? abs(δ) : inv(abs(δ))
            complete = abs(δ) > one(T) || small * small < floatmin(T)
        end
        if abs(δ) <= one(T)
            δ2 = δ * δ
            denom = one(T) + δ2
            cosθ = (one(T) - δ2) / denom
            cosθ_minus_one = -2δ2 / denom
            sinθ = 2δ / denom
        else
            invδ = one(T) / δ
            invδ2 = invδ * invδ
            denom = one(T) + invδ2
            cosθ = (invδ2 - one(T)) / denom
            cosθ_minus_one = -2 / denom
            sinθ = 2invδ / denom
        end
    end
    φ = twoπ * rand(rng, T)
    (isfinite(cosθ) && isfinite(sinθ) && isfinite(φ)) ||
        throw(ArgumentError("collide_coulomb!: sampled scattering direction must be finite"))
    sinφ, cosφ = sincos(φ)

    # Build the physical increment directly; reconstructing a normalized
    # rotated vector would add angle-independent noise before subtraction.
    invgmag = inv(gmag)
    if isfinite(invgmag)
        ugx = gx * invgmag
        ugy = gy * invgmag
        ugz = gz * invgmag
    else
        ugx = gx / gmag
        ugy = gy / gmag
        ugz = gz / gmag
    end
    dugx, dugy, dugz = ugx, ugy, ugz
    ugperp = hypot(ugx, ugy)
    if ugperp > 0                           # Takizuka & Abe 1977 rotation
        invugperp = inv(ugperp)
        if isfinite(invugperp)
            xperp = ugx * invugperp
            yperp = ugy * invugperp
        else
            xperp = ugx / ugperp
            yperp = ugy / ugperp
        end
        polar_term = ugz * sinθ * cosφ
        azimuthal_term = sinθ * sinφ
        lx, ly, lz = ugx * cosθ_minus_one, ugy * cosθ_minus_one, ugz * cosθ_minus_one
        px, py = xperp * polar_term, yperp * polar_term
        ax, ay = yperp * azimuthal_term, xperp * azimuthal_term
        pz = ugperp * sinθ * cosφ
        dugx, dugy, dugz = lx + px - ax, ly + py + ay, lz - pz
        # Every factor can be normal while a unit-scale product vanishes.
        # The later multiplication by physical speed can recover that term,
        # so preserve the sampled expression before information is lost.
        complete =
            complete ||
            _coulomb_product_needs_complete(ugx, cosθ_minus_one) ||
            _coulomb_product_needs_complete(ugy, cosθ_minus_one) ||
            _coulomb_product_needs_complete(ugz, cosθ_minus_one) ||
            _coulomb_product_needs_complete(ugz, sinθ) ||
            _coulomb_product_needs_complete(ugz * sinθ, cosφ) ||
            _coulomb_product_needs_complete(sinθ, sinφ) ||
            _coulomb_product_needs_complete(xperp, polar_term) ||
            _coulomb_product_needs_complete(yperp, polar_term) ||
            _coulomb_product_needs_complete(yperp, azimuthal_term) ||
            _coulomb_product_needs_complete(xperp, azimuthal_term) ||
            _coulomb_product_needs_complete(ugperp, sinθ) ||
            _coulomb_product_needs_complete(ugperp * sinθ, cosφ) ||
            _coulomb_increment_needs_complete(dugx, lx, px, ax) ||
            _coulomb_increment_needs_complete(dugy, ly, py, ay) ||
            _coulomb_increment_needs_complete(dugz, lz, pz, zero(T))
    else                                    # g along ±z
        dugx, dugy, dugz = sinθ * cosφ, sinθ * sinφ, ugz * cosθ_minus_one
        complete =
            complete ||
            _coulomb_product_needs_complete(sinθ, cosφ) ||
            _coulomb_product_needs_complete(sinθ, sinφ) ||
            _coulomb_product_needs_complete(ugz, cosθ_minus_one)
    end
    # Normalizing a very oblique vector can discard a physical component.
    complete |=
        (!iszero(gx) && abs(ugx) < floatmin(T)) ||
        (!iszero(gy) && abs(ugy) < floatmin(T)) ||
        (!iszero(gz) && abs(ugz) < floatmin(T))
    increment = (gmag * dugx, gmag * dugy, gmag * dugz)
    return (;
        increment,
        tangent,
        tangent_power,
        isotropic = large_angle,
        cosine = cosθ,
        sine = sinθ,
        cosφ,
        sinφ,
        complete,
        identity,
    )
end

@noinline function _coulomb_mass_update_exact(va::T, vb::T, ma::T, mb::T, increment::T) where {T}
    a, b = _exact_rational(ma), _exact_rational(mb)
    # The caller routes nonfinite increments to the complete sampled-angle
    # owner; this path completes a finite mass-weighted addition.
    change = _exact_rational(increment)
    return _from_exact(T, _exact_rational(va) + b * change / (a + b)),
    _from_exact(T, _exact_rational(vb) - a * change / (a + b))
end

@inline function _coulomb_mass_update(va::T, vb::T, ma::T, mb::T, fa::T, fb::T, change::T) where {T}
    # A null relative change is the identity even when the two original
    # velocities have very different magnitudes. Reconstructing a midpoint
    # first would discard the smaller velocity through cancellation.
    iszero(change) && return va, vb
    kick_a, kick_b = fb * change, fa * change
    newa, newb = va + kick_a, vb - kick_b
    if isfinite(newa) &&
       isfinite(newb) &&
       abs(kick_a) >= floatmin(T) &&
       abs(kick_b) >= floatmin(T) &&
       abs(newa) >= floatmin(T) &&
       abs(newb) >= floatmin(T) &&
       fa >= floatmin(T) &&
       fb >= floatmin(T)
        return newa, newb
    end
    return _coulomb_mass_update_exact(va, vb, ma, mb, change)
end

# One physical expression for widened and exact complete-angle arithmetic.
@inline function _coulomb_complete_expression(
    ga,
    gb,
    a,
    b,
    magnitude,
    radius,
    phase_cos,
    phase_sin,
    cosine_change,
    sine,
)
    g = map(-, ga, gb)
    rotation = if iszero(radius)
        (magnitude * phase_cos, magnitude * phase_sin, zero(magnitude))
    else
        (
            (g[1] * g[3] * phase_cos - g[2] * magnitude * phase_sin) / radius,
            (g[2] * g[3] * phase_cos + g[1] * magnitude * phase_sin) / radius,
            -radius * phase_cos,
        )
    end
    change = ntuple(c -> g[c] * cosine_change + rotation[c] * sine, 3)
    newa = ntuple(c -> ga[c] + b * change[c] / (a + b), 3)
    newb = ntuple(c -> gb[c] - a * change[c] / (a + b), 3)
    return newa, newb, change
end

@inline function _coulomb_widened_update(va, vb, ma, mb, speed, sample)
    if sample.isotropic
        cosine_change = Float64(sample.cosine) - 1
        sine = Float64(sample.sine)
    else
        tangent = ldexp(Float64(sample.tangent), sample.tangent_power)
        iszero(sample.tangent) && return va, vb
        squared = tangent * tangent
        (isfinite(squared) && squared >= floatmin(Float64)) || return nothing
        cosine_change = -2squared / (1 + squared)
        sine = 2tangent / (1 + squared)
    end
    ga, gb = map(Float64, va), map(Float64, vb)
    radius = Float64(hypot(va[1] - vb[1], va[2] - vb[2]))
    newa, newb, change = _coulomb_complete_expression(
        ga,
        gb,
        Float64(ma),
        Float64(mb),
        Float64(speed),
        radius,
        Float64(sample.cosφ),
        Float64(sample.sinφ),
        cosine_change,
        sine,
    )
    # Cancellation can occur inside the rotation before the final sum. Bound
    # acceptance by the pair scale, including coordinates whose rounded change
    # is zero. Reserve four extra bits beyond half of Float64's significand
    # for accumulated arithmetic error before rounding back to Float32.
    scale = max(Float64(speed), maximum(abs, ga), maximum(abs, gb))
    threshold = 16sqrt(eps(Float64)) * scale
    # A widened sum can round across Float32's overflow midpoint even when
    # its exact result rounds to a finite Float32. Complete that boundary
    # exactly instead of converting an already rounded midpoint to infinity.
    limit = Float64(floatmax(Float32))
    all(x -> isfinite(x) && threshold <= abs(x) <= limit, (newa..., newb...)) || return nothing
    return ntuple(c -> iszero(change[c]) ? va[c] : Float32(newa[c]), 3),
    ntuple(c -> iszero(change[c]) ? vb[c] : Float32(newb[c]), 3)
end

# Bounded compensated operands reuse the same complete physical expression.
@inline function _coulomb_compensated_update(
    va::NTuple{3,Float64},
    vb::NTuple{3,Float64},
    ma,
    mb,
    speed,
    sample,
)
    DP = Base.TwicePrecision{Float64}
    margin =
        fld(min(exponent(floatmax(Float64)), -exponent(floatmin(Float64))) - 4precision(Float64), 8)
    low, high = ldexp(1.0, -margin), ldexp(1.0, margin)
    inrange(x) = iszero(x) || (low <= abs(x) <= high)
    radius = hypot(va[1] - vb[1], va[2] - vb[2])
    all(inrange, (va..., vb..., ma, mb, speed, radius, sample.cosφ, sample.sinφ)) || return nothing
    if sample.isotropic
        all(inrange, (sample.cosine, sample.sine)) || return nothing
        cosine_change = DP(sample.cosine) - 1
        sine = DP(sample.sine)
    else
        iszero(sample.tangent) && return va, vb
        sample.tangent_power == 0 && inrange(sample.tangent) || return nothing
        tangent = DP(sample.tangent)
        squared = tangent * tangent
        cosine_change = -2squared / (1 + squared)
        sine = 2tangent / (1 + squared)
    end
    newa, newb, change = _coulomb_complete_expression(
        map(DP, va),
        map(DP, vb),
        DP(ma),
        DP(mb),
        DP(speed),
        DP(radius),
        DP(sample.cosφ),
        DP(sample.sinφ),
        cosine_change,
        sine,
    )
    all(x -> isfinite(x.hi) && isfinite(x.lo), (newa..., newb..., change...)) || return nothing
    any(x -> iszero(x.hi) && iszero(x.lo), change) && return nothing
    a, b = map(Float64, newa), map(Float64, newb)
    scale = max(speed, maximum(abs, va), maximum(abs, vb))
    all(x -> abs(x) > 256eps(Float64) * scale, (a..., b...)) || return nothing
    return a, b
end

# Retain the sampled angle through the final physical velocity, instead of
# trying to recover information from a rounded cosine or increment.
@noinline function _coulomb_complete_update(
    va::NTuple{3,T},
    vb::NTuple{3,T},
    ma::T,
    mb::T,
    speed::T,
    sample,
) where {T}
    if T === Float32
        widened = _coulomb_widened_update(va, vb, ma, mb, speed, sample)
        widened === nothing || return widened
    elseif T === Float64
        compensated = _coulomb_compensated_update(va, vb, ma, mb, speed, sample)
        compensated === nothing || return compensated
    end
    a, b = _exact_rational(ma), _exact_rational(mb)
    ga = map(_exact_rational, va)
    gb = map(_exact_rational, vb)
    if sample.isotropic
        cosine_change = _exact_rational(sample.cosine) - 1
        sine = _exact_rational(sample.sine)
    else
        tangent = _exact_rational(sample.tangent)
        power = sample.tangent_power
        tangent = power >= 0 ? tangent * (BigInt(1) << power) : tangent / (BigInt(1) << -power)
        squared = tangent * tangent
        cosine_change = -2squared / (1 + squared)
        sine = 2tangent / (1 + squared)
    end
    phase_cos, phase_sin = _exact_rational(sample.cosφ), _exact_rational(sample.sinφ)
    magnitude = _exact_rational(speed)
    perpendicular = hypot(va[1] - vb[1], va[2] - vb[2])
    result_a, result_b, change = _coulomb_complete_expression(
        ga,
        gb,
        a,
        b,
        magnitude,
        _exact_rational(perpendicular),
        phase_cos,
        phase_sin,
        cosine_change,
        sine,
    )
    # Preserve exact identities, including the sign of an unchanged zero.
    newa = ntuple(c -> iszero(change[c]) ? va[c] : _from_exact(T, result_a[c]), 3)
    newb = ntuple(c -> iszero(change[c]) ? vb[c] : _from_exact(T, result_b[c]), 3)
    return newa, newb
end

@inline function _coulomb_complete_needed(va::T, vb::T, change::T, fa::T, fb::T) where {T}
    isfinite(change) || return true
    iszero(change) && return false
    kick_a, kick_b = fb * change, fa * change
    newa, newb = va + kick_a, vb - kick_b
    # Preserve component accuracy when the final addition loses a leading bit,
    # including the small near-backscatter remainder.
    threshold = one(T) / 2
    return fa < floatmin(T) ||
           fb < floatmin(T) ||
           abs(kick_a) < floatmin(T) ||
           abs(kick_b) < floatmin(T) ||
           !isfinite(newa) ||
           !isfinite(newb) ||
           (signbit(va) != signbit(kick_a) && abs(newa) <= threshold * max(abs(va), abs(kick_a))) ||
           (signbit(vb) == signbit(kick_b) && abs(newb) <= threshold * max(abs(vb), abs(kick_b)))
end
