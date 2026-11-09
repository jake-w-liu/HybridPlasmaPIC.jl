# Derived scratch only: neither a pressure reservoir nor restart carry.
struct _FieldAlignedConductionWorkspace{D}
    norm2::Array{Rational{BigInt},D}
    gradient::NTuple{D,Array{Rational{BigInt},D}}
end

_field_aligned_storage(::FourierGrid, ::ElectronClosure) = nothing
function _field_aligned_storage(g::FourierGrid{D}, ::FieldAlignedConductingElectrons) where {D}
    return _FieldAlignedConductionWorkspace{D}(
        zeros(Rational{BigInt}, g.n),
        ntuple(_ -> zeros(Rational{BigInt}, g.n), D),
    )
end

_field_aligned_arrays(::Nothing) = ()
_field_aligned_arrays(w::_FieldAlignedConductionWorkspace) = (w.norm2, w.gradient...)

_validate_field_aligned_direction(::ElectronClosure, B, ::Type) = nothing
function _validate_field_aligned_direction(
    closure::FieldAlignedConductingElectrons,
    B,
    ::Type{T},
    norm2 = nothing,
) where {T}
    parallel, perpendicular = _field_aligned_conductivities(closure, T)
    parallel == perpendicular && return nothing
    for I in eachindex(B[1], B[2], B[3])
        all(c -> isfinite(B[c][I]), 1:3) ||
            throw(ArgumentError("field-aligned conduction requires finite magnetic direction"))
        magnitude = sum(c -> _exact_rational(B[c][I])^2, 1:3)
        iszero(magnitude) &&
            throw(ArgumentError("unequal conductivities require nonzero magnetic direction"))
        norm2 === nothing || (norm2[I] = magnitude)
    end
    return nothing
end

# Project a computed real channel to exact zero mean on each axis line. This
# is an error-nonincreasing projection for the ideal derivative's range.
function _field_aligned_project!(out, values, d, g::FourierGrid{D}) where {D}
    for origin in CartesianIndices(ntuple(c -> c == d ? 1 : g.n[c], D))
        total = zero(Rational{BigInt})
        for j = 1:g.n[d]
            I = CartesianIndex(ntuple(c -> c == d ? j : origin[c], D))
            total += real(values[I])
        end
        mean = total / g.n[d]
        for j = 1:g.n[d]
            I = CartesianIndex(ntuple(c -> c == d ? j : origin[c], D))
            out[I] = real(values[I]) - mean
        end
    end
    return out
end

# One real odd-derivative or positive Nyquist channel. Both the original
# quotient and the generally nondyadic tensor flux use this same bounded
# preparation; only dyadic terms enter the banded FFT owner.
function _field_aligned_channel!(out, values, d, nyquist, g::FourierGrid{D,T}) where {D,T}
    if nyquist && isodd(g.n[d])
        fill!(out, zero(Rational{BigInt}))
        return out
    end
    Q = Rational{BigInt}
    input = Array{Complex{Q},D}(undef, g.n)
    active = false
    for I in CartesianIndices(values)
        anchor = CartesianIndex(ntuple(c -> c == d ? 1 : I[c], D))
        residual = values[I] - values[anchor]
        active |= !iszero(residual)
        input[I] = complex(_conduction_dyadic(T, residual))
    end
    if !active
        fill!(out, zero(Q))
        return out
    end
    spectrum = _em_dyadic_transform(input, g)
    nyquist_index = g.n[d] ÷ 2 + 1
    for I in CartesianIndices(spectrum)
        factor = if nyquist
            complex(I[d] == nyquist_index ? abs(_exact_rational(g.kfull[d][I[d]])) : zero(Q))
        else
            complex(zero(Q), _exact_rational(g.kvec[d][I[d]]))
        end
        spectrum[I] *= factor
    end
    values_out = _em_dyadic_transform(spectrum, g, true)
    return _field_aligned_project!(out, values_out, d, g)
end

# Exact arithmetic here retains the original stored B components and the
# resolved grid coefficients. The norm is never formed in native arithmetic.
function _field_aligned_diagonal(B, norm2, I, d, perpendicular, difference)
    return perpendicular + difference * _exact_rational(B[d][I])^2 / norm2[I]
end

function _field_aligned_conduction!(
    exact,
    quotient,
    coefficient,
    closure::FieldAlignedConductingElectrons,
    B,
    w::_FieldAlignedConductionWorkspace{D},
    g::FourierGrid{D,T};
    budget = false,
) where {D,T}
    _require_conduction_quotient(quotient, g)
    parallelT, perpendicularT = _field_aligned_conductivities(closure, T)
    Q = Rational{BigInt}
    parallel, perpendicular = _exact_rational(parallelT), _exact_rational(perpendicularT)
    difference = parallel - perpendicular
    iszero(parallel) && iszero(perpendicular) && return zero(Q)
    SpectralOperators._laplacian_gain(g.kfull)
    _validate_field_aligned_direction(closure, B, T, w.norm2)
    if iszero(difference) && !budget
        _add_conduction_pressure!(exact, quotient, coefficient * perpendicular, g)
        return zero(Q)
    end
    for d = 1:D
        _field_aligned_channel!(w.gradient[d], quotient, d, false, g)
    end
    dissipation = zero(Q)
    if budget
        for I in eachindex(quotient)
            square = sum(d -> w.gradient[d][I]^2, 1:D)
            aligned =
                iszero(difference) ? zero(Q) :
                sum(d -> _exact_rational(B[d][I]) * w.gradient[d][I], 1:D)
            dissipation +=
                perpendicular * square +
                (iszero(difference) ? zero(Q) : difference * aligned^2 / w.norm2[I])
        end
    end
    flux, channel = similar(quotient), similar(quotient)
    for d = 1:D
        if !iszero(difference)
            for I in eachindex(flux)
                aligned = sum(j -> _exact_rational(B[j][I]) * w.gradient[j][I], 1:D)
                flux[I] =
                    perpendicular * w.gradient[d][I] +
                    difference * _exact_rational(B[d][I]) * aligned / w.norm2[I]
            end
            _field_aligned_channel!(channel, flux, d, false, g)
            for I in eachindex(exact, channel)
                exact[I] += coefficient * channel[I]
            end
        end
        if iseven(g.n[d])
            _field_aligned_channel!(channel, quotient, d, true, g)
            for I in eachindex(flux, channel)
                diagonal =
                    iszero(difference) ? perpendicular :
                    _field_aligned_diagonal(B, w.norm2, I, d, perpendicular, difference)
                budget && (dissipation += diagonal * channel[I]^2)
                flux[I] = diagonal * channel[I]
            end
            if !iszero(difference)
                _field_aligned_channel!(channel, flux, d, true, g)
                for I in eachindex(exact, channel)
                    exact[I] -= coefficient * channel[I]
                end
            end
        end
    end
    # Equal coefficients retain policy-2 arithmetic even in a diagnostic call.
    iszero(difference) && _add_conduction_pressure!(exact, quotient, coefficient * perpendicular, g)
    return dissipation
end

function _add_transport_conduction!(
    exact,
    quotient,
    coefficient,
    κ,
    g,
    closure::FieldAlignedConductingElectrons,
    B,
    work,
)
    _field_aligned_conduction!(exact, quotient, coefficient, closure, B, work, g)
    return nothing
end
