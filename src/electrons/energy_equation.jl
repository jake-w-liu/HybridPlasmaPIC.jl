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
function _add_ohmic_pressure_heating!(out, J, γ::T, η::T) where {T}
    γm1 = γ - one(T)
    coefficient = γm1 * η
    fγ, eγ = frexp(γm1)
    fη, eη = frexp(η)
    mantissa = fγ * fη
    exponent = eγ + eη
    ordinary = isfinite(coefficient) && coefficient >= floatmin(T)
    Jx, Jy, Jz = J
    @inbounds for I in eachindex(out, Jx, Jy, Jz)
        jx, jy, jz = Jx[I], Jy[I], Jz[I]
        square = jx * jx + jy * jy + jz * jz
        heating =
            ordinary && isfinite(square) && square >= floatmin(T) ? coefficient * square :
            _scaled_ohmic_heating(mantissa, exponent, jx, jy, jz)
        out[I] += heating
    end
    return out
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
