# EnergyEquation.jl — electron pressure evolution (from hybrid.jl)

"""
    advance_electron_pressure!(pe, ue, dt, γe, g)

One forward-Euler step of the adiabatic electron-pressure equation
`∂p_e/∂t = −u_e·∇p_e − γ_e p_e ∇·u_e` on a periodic grid (the energy-equation
closure; source terms −∇·q_e, resistive heating are optional and omitted here).
This standalone building block has no stage-positivity or rollback guarantee.
Use [`EvolvingElectrons`](@ref) with `HybridStepper` for coupled B/pressure RK4
subcycles and validated persistent pressure state.
"""
function advance_electron_pressure!(
    pe::Array{T,D},
    ue::NTuple{3,<:Array{T,D}},
    dt,
    γe,
    g::FourierGrid{D,T},
) where {D,T}
    dtT = _validated_nonnegative_dt(T, dt; name = "advance_electron_pressure!")
    γ = _require_finite_positive_real("γe", γe, T)
    size(pe) == g.n ||
        throw(DimensionMismatch("pe size $(size(pe)) does not match grid size $(g.n)"))
    for c = 1:3
        size(ue[c]) == g.n ||
            throw(DimensionMismatch("ue[$c] size $(size(ue[c])) does not match grid size $(g.n)"))
    end
    iszero(dtT) && return pe

    gradpe = ntuple(_ -> similar(pe), D)
    gradient!(gradpe, pe, g)
    divu = similar(pe)
    divergence!(divu, ue, g)
    @inbounds for I in eachindex(pe)
        pe[I] += dtT * _electron_pressure_rate(pe, ue, gradpe, divu, γ, I, Val(D))
    end
    return pe
end


# The primitive equation has one arithmetic owner for the standalone Euler
# update and the coupled RK stages. Derivatives are caller-owned scratch.
@inline function _electron_pressure_rate(pe, ue, gradpe, divu, γ, I, ::Val{D}) where {D}
    adv = zero(eltype(pe))
    @inbounds for d = 1:D
        adv += ue[d][I] * gradpe[d][I]
    end
    @inbounds return -adv - γ * pe[I] * divu[I]
end

function _electron_pressure_rhs!(out, pe, ue, gradpe, divu, γ, ::Val{D}) where {D}
    @inbounds for I in eachindex(out, pe, divu)
        out[I] = _electron_pressure_rate(pe, ue, gradpe, divu, γ, I, Val(D))
    end
    return out
end

# Capture possible represented pressure sources before derivative scratch is
# reused for heating. A contributing factor pair is only a routing condition;
# the complete cold equation decides exact cancellation and physical range.
function _transport_pressure_rhs_with_activity!(
    out,
    pe,
    ue,
    gradpe,
    divu,
    γ,
    J,
    η,
    active::Bool,
    ::Val{D},
) where {D}
    if !active && !iszero(η)
        for component in J
            for value in component
                if !iszero(value)
                    active = true
                    break
                end
            end
            active && break
        end
    end
    @inbounds for I in eachindex(out, pe, divu)
        out[I] = _electron_pressure_rate(pe, ue, gradpe, divu, γ, I, Val(D))
        if !active
            active = !iszero(pe[I]) && !iszero(divu[I])
            for d = 1:D
                active |= !iszero(ue[d][I]) && !iszero(gradpe[d][I])
            end
        end
    end
    return active
end

# This scan requests complete signed replay; it never certifies a nonzero
# source or rejects a derivative channel independently of the pressure equation.
function _scalar_pressure_rhs_needs_completion(
    out::Array{T},
    pe,
    ue,
    gradpe,
    divu,
    J,
    γ::T,
    η::T,
    active::Bool,
    ::Val{D},
) where {T,D}
    magnitude = zero(T)
    γm1 = _pressure_heat_gamma(γ)
    coefficient = γm1 * η
    @inbounds for I in eachindex(out, pe, divu)
        isfinite(out[I]) && isfinite(divu[I]) || return true
        magnitude = max(magnitude, abs(out[I]))
        weighted = γ * pe[I]
        _ohm_product_needs_completion(γ, pe[I], weighted) && return true
        _ohm_product_needs_completion(weighted, divu[I], weighted * divu[I]) && return true
        for d = 1:D
            isfinite(gradpe[d][I]) || return true
            _ohm_product_needs_completion(ue[d][I], gradpe[d][I], ue[d][I] * gradpe[d][I]) &&
                return true
        end
        if !iszero(η)
            square = zero(T)
            nonzero_current = false
            for c = 1:3
                current = J[c][I]
                nonzero_current |= !iszero(current)
                term = current * current
                _ohm_product_needs_completion(current, current, term) && return true
                square += term
            end
            nonzero_current && _ohm_product_needs_completion(γm1, η, coefficient) && return true
            isfinite(square) || return true
            _ohm_product_needs_completion(coefficient, square, coefficient * square) && return true
        end
    end
    return active && magnitude < 32floatmin(T)
end

function _complete_scalar_pressure_if_needed!(
    out::Array{T,D},
    pe,
    ue,
    gradpe,
    divu,
    J,
    γ::T,
    η::T,
    active::Bool,
    g::FourierGrid{D,T};
    derivative_risk::Bool = false,
) where {D,T}
    derivative_risk ||
        _scalar_pressure_rhs_needs_completion(out, pe, ue, gradpe, divu, J, γ, η, active, Val(D)) ||
        return false
    # Only original represented inputs enter exact extraction. Rounded
    # derivative/rate scratch may be nonfinite and is recomputed by its owner.
    _validate_pressure_finite(pe, "complete electron pressure input")
    _validate_pressure_finite((γ, η), "complete electron pressure coefficients")
    for d = 1:D
        _validate_pressure_finite(ue[d], "complete electron pressure velocity")
    end
    if !iszero(η)
        for component in J
            _validate_pressure_finite(component, "complete electron pressure current")
        end
    end
    _complete_transport_pressure_rhs!(
        out,
        pe,
        ue,
        J,
        γ,
        η,
        zero(T),
        g;
        context = "complete electron pressure rate",
    )
    return true
end

# (γ−1)η|J|² from the exact power-of-two factors of the coefficient: the
# largest current component is scaled to [1, 2), so only the final ldexp can
# underflow (within one subnormal ulp) or overflow (rejected as non-finite).
@inline function _scaled_ohmic_heating(mantissa::T, exponent::Int, jx::T, jy::T, jz::T) where {T}
    s = max(abs(jx), abs(jy), abs(jz))
    iszero(s) && return zero(T)
    e = Base.exponent(s)
    x, y, z = ldexp(jx, -e), ldexp(jy, -e), ldexp(jz, -e)
    return ldexp(mantissa * (x * x + y * y + z * z), exponent + 2e)
end

# Add the Ohmic pressure source (γ−1)η|J|² of the same RK stage current `J`.
# The ordinary product is used when the coefficient and squared current are
# normal and finite: then underflowed component squares or a subnormal result
# cost at most a few ulps of the operation scale. Other cells take the scaled
# evaluation of the same expression; a non-finite result is left for the
# caller's finite-rate rejection.
@inline function _pressure_heating_factors(γ::T, η::T) where {T}
    γm1 = _pressure_heat_gamma(γ)
    coefficient = γm1 * η
    fγ, eγ = frexp(γm1)
    fη, eη = frexp(η)
    mantissa = fγ * fη
    exponent = eγ + eη
    ordinary = isfinite(coefficient) && coefficient >= floatmin(T)
    return coefficient, mantissa, exponent, ordinary
end

@inline _pressure_heat_gamma(γ::T) where {T} = γ - one(T)

@inline function _squared_pressure_heating(jx, jy, jz, coefficient, mantissa, exponent, ordinary)
    square = jx * jx + jy * jy + jz * jz
    return ordinary && isfinite(square) && square >= floatmin(typeof(jx)) ? coefficient * square :
           _scaled_ohmic_heating(mantissa, exponent, jx, jy, jz)
end

function _add_ohmic_pressure_heating!(out, J, γ::T, η::T) where {T}
    factors = _pressure_heating_factors(γ, η)
    Jx, Jy, Jz = J
    @inbounds for I in eachindex(out, Jx, Jy, Jz)
        jx, jy, jz = Jx[I], Jy[I], Jz[I]
        out[I] += _squared_pressure_heating(jx, jy, jz, factors...)
    end
    return out
end

# Each real channel contributes a positive square through the same range-safe
# numerical owner as Ohmic heating. Both scratch owners have been consumed
# by advection; gradient! shares one forward transform across spatial axes.
function _add_hyperresistive_pressure_heating!(
    out,
    J,
    gradJ,
    scratch,
    γ::T,
    ηH::T,
    g::FourierGrid{D,T},
) where {D,T}
    _add_represented_hyperresistive_heating!(out, J, gradJ, scratch, γ, ηH, g) || throw(
        ArgumentError("hyperresistive pressure channels require complete scaled rate evaluation"),
    )
    return out
end

# The coupled owner may replay the complete equation when a physical derivative
# channel cannot be stored. Never feed Inf to the scaled square's exponent owner.
function _add_represented_hyperresistive_heating!(
    out,
    J,
    gradJ,
    scratch,
    γ::T,
    ηH::T,
    g::FourierGrid{D,T},
) where {D,T}
    iszero(ηH) && return true
    factors = _pressure_heating_factors(γ, ηH)
    for c = 1:3
        gradient!(gradJ, J[c], g)
        for d = 1:D
            all(isfinite, gradJ[d]) || return false
            _add_pressure_heating_channel!(out, gradJ[d], factors)
            if iseven(g.n[d])
                nyquist_channel!(scratch, J[c], g, d)
                all(isfinite, scratch) || return false
                _add_pressure_heating_channel!(out, scratch, factors)
            end
        end
    end
    return true
end

function _add_pressure_heating_channel!(out::Array{T}, channel, factors) where {T}
    @inbounds for I in eachindex(out, channel)
        out[I] += _squared_pressure_heating(channel[I], zero(T), zero(T), factors...)
    end
    return out
end

# Classification follows the complete signed sum, never an individual channel.
# Magnetic components use the returned flags jointly before norm rejection.
function _publish_complete_rate!(out::Array{T}, exact, context) where {T}
    nonzero, represented = false, false
    for I in eachindex(out, exact)
        value = _from_exact(T, exact[I])
        _validate_pressure_finite((value,), context)
        nonzero |= !iszero(exact[I])
        represented |= !iszero(value)
        out[I] = value
    end
    return nonzero, represented
end

# The physical temperature is a derived stage operand, not a carried reservoir.
# Original division rounds once; reciprocal multiplication would define a
# different represented source. Scratch may change on rejection.
function _conduction_temperature!(temperature::Array{T,D}, pe, n, g::FourierGrid{D,T}) where {D,T}
    _boris_accuracy_mode() ||
        throw(ArgumentError("active electron conduction requires nearest rounding and gradual underflow"))
    _require_grid_array(:temperature, temperature, g)
    _require_grid_array(:pe, pe, g)
    _require_grid_array(:n, n, g)
    (_particle_arrays_mightalias(temperature, pe) || _particle_arrays_mightalias(temperature, n)) &&
        throw(ArgumentError("electron temperature scratch must be independent of pressure and density"))
    for I in eachindex(temperature, pe, n)
        pressure, density = pe[I], n[I]
        isfinite(pressure) && pressure >= zero(T) && isfinite(density) && density > zero(T) ||
            throw(ArgumentError("electron temperature requires finite nonnegative pressure and positive density"))
        exact = _exact_rational(pressure) / _exact_rational(density)
        value = _from_exact(T, exact)
        isfinite(value) || throw(ArgumentError("electron temperature does not round finite"))
        (!iszero(exact) && iszero(value)) &&
            throw(ArgumentError("nonzero electron temperature rounds to zero"))
        temperature[I] = value
    end
    return temperature
end

@noinline function _complete_transport_pressure_rhs!(
    out::Array{T,D},
    pe,
    ue,
    J,
    γ::T,
    η::T,
    ηH::T,
    g::FourierGrid{D,T};
    context = "complete transport pressure rate",
    κ::T = zero(T),
    temperature = nothing,
) where {D,T}
    exact = zeros(Rational{BigInt}, size(out))
    gamma = _exact_rational(γ)
    for d = 1:D
        data = SpectralOperators.deriv_scaled_data!(g, pe, d)
        factor = _exact_binary_scale(data.input_scale, data.exponent)
        for I in eachindex(exact, pe, ue[d], data.values)
            exact[I] -= _exact_rational(ue[d][I]) * _exact_rational(real(data.values[I])) * factor
        end
        data = SpectralOperators.deriv_scaled_data!(g, ue[d], d)
        factor = gamma * _exact_binary_scale(data.input_scale, data.exponent)
        for I in eachindex(exact, pe, data.values)
            exact[I] -= _exact_rational(pe[I]) * _exact_rational(real(data.values[I])) * factor
        end
    end
    gamma_minus_one = _exact_rational(_pressure_heat_gamma(γ))
    if !iszero(η)
        coefficient = gamma_minus_one * _exact_rational(η)
        for c = 1:3, I in eachindex(exact, J[c])
            current = _exact_rational(J[c][I])
            exact[I] += coefficient * current * current
        end
    end
    if !iszero(ηH)
        coefficient = gamma_minus_one * _exact_rational(ηH)
        for c = 1:3, d = 1:D
            data = SpectralOperators.deriv_scaled_data!(g, J[c], d)
            _add_exact_pressure_channel!(exact, data, coefficient)
            if iseven(g.n[d])
                data = SpectralOperators.nyquist_scaled_data!(g, J[c], d)
                _add_exact_pressure_channel!(exact, data, coefficient)
            end
        end
    end
    if !iszero(κ)
        temperature isa Array{T,D} ||
            throw(ArgumentError("active conduction requires grid-precision temperature storage"))
        data = SpectralOperators.laplacian_scaled_data!(g, temperature)
        coefficient = gamma_minus_one * _exact_rational(κ) *
            _exact_binary_scale(data.input_scale, data.exponent)
        # Consume borrowed Fourier data before any other operator uses g.
        # No physical Laplacian or partial conductive rate is rounded here.
        for I in eachindex(exact, data.values)
            exact[I] += coefficient * _exact_rational(real(data.values[I]))
        end
    end
    nonzero, represented = _publish_complete_rate!(out, exact, context)
    nonzero && !represented && throw(ArgumentError("$context rounds entirely to zero in $T"))
    return nothing
end

function _add_exact_pressure_channel!(exact, data, coefficient)
    factor = _exact_binary_scale(data.input_scale, data.exponent)
    scaled_coefficient = coefficient * factor * factor
    for I in eachindex(exact, data.values)
        value = _exact_rational(real(data.values[I]))
        exact[I] += scaled_coefficient * value * value
    end
    return nothing
end

function _validate_evolving_pressure(pe, g; context = "electron pressure", subcycle = 0, stage = 0)
    _require_grid_array(:pe, pe, g)
    for value in pe
        isfinite(value) && value >= zero(value) || throw(
            ArgumentError(
                "$context must be finite and nonnegative (subcycle $subcycle, stage $stage)",
            ),
        )
    end
    return nothing
end

function _validate_pressure_finite(values, context)
    for value in values
        isfinite(value) || throw(ArgumentError("$context must be finite"))
    end
    return nothing
end

# Stored pressure is a physical input to Ohm's law. Validate its independence
# once at the public field boundary, also shared by stepping and persistence.
function _validate_evolving_pressure_storage(f::HybridFields{D}, g) where {D}
    for (name, a) in ((:n, f.n), (:pe, f.pe), (:ninv, f.ninv))
        _require_grid_array(name, a, g)
    end
    for (name, arrays) in
        ((:ui, f.ui), (:B, f.B), (:E, f.E), (:J, f.J), (:lapJ, f.lapJ), (:gradp, f.gradp))
        _require_grid_tuple(name, arrays, g)
    end
    scalar_shape = ntuple(_ -> 0, D)
    for a in (f.pforce..., f.cgl_dp)
        size(a) == scalar_shape ||
            throw(DimensionMismatch("evolving scalar pressure requires empty CGL scratch"))
    end
    for a in
        (f.n, f.ui..., f.B..., f.E..., f.J..., f.ninv, f.gradp..., f.lapJ..., f.pforce..., f.cgl_dp)
        _particle_arrays_mightalias(f.pe, a) &&
            throw(ArgumentError("evolving pressure must have independent storage"))
    end
    return nothing
end
