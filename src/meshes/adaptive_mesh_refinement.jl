# adaptive_mesh_refinement.jl -- block-structured AMR primitives.
#
# These utilities operate on cell-centred finite-volume levels. They are kept in
# the mesh subsystem because they define mesh resolution, cell locations, and
# conservative inter-level transfer operators.

"""
    AMRGrid(u, dx; x0=0.0)

One refinement level: cell-centred data `u` (length `N`) on `[x0, x0+N*dx]`,
with cells of width `dx` centred at `x0 + (i-1/2)dx`.
"""
struct AMRGrid{T<:AbstractFloat}
    u::Vector{T}
    dx::T
    x0::T

    function AMRGrid{T}(u::Vector{T}, dx::T, x0::T) where {T<:AbstractFloat}
        isempty(u) && throw(ArgumentError("AMRGrid requires at least one cell"))
        isfinite(dx) && dx > zero(T) ||
            throw(ArgumentError("AMRGrid dx must be finite and positive"))
        isfinite(x0) || throw(ArgumentError("AMRGrid x0 must be finite"))
        return new{T}(u, dx, x0)
    end
end

function AMRGrid(u::AbstractVector{S}, dx::Real; x0::Real = 0.0) where {S<:Real}
    T = S <: AbstractFloat ? S : typeof(float(zero(S)))
    values = u isa Vector{T} ? u : T.(u)
    dxT = _require_finite_positive_real("AMRGrid dx", dx, T)
    x0T = _require_finite_real("AMRGrid x0", x0, T)
    return AMRGrid{T}(values, dxT, x0T)
end

ncells(g::AMRGrid) = length(g.u)
effective_resolution(g::AMRGrid) = g.dx

function cell_center(g::AMRGrid{T}, i::Integer) where {T}
    W = promote_type(T, typeof(i))
    if W <: Union{Float16,Float32,Float64}
        F = W === Float16 ? Float64 : W
        bound = Int64(1) << (precision(F) - 1)
        if 1 - bound <= i <= bound
            # Within these bounds the half-integer coefficient is exact.
            factor = F(i) - F(0.5)
            # Float16 fma can double-round through Float32. Its exact half-cell
            # coordinate lies on the 2^-25 lattice, so Float64 represents every
            # result near Float16's finite range exactly, including its midpoints.
            return W(fma(factor, F(g.dx), F(g.x0)))
        end
    end
    return _amr_exact_cell_center(g, i, W)
end

@noinline function _amr_exact_cell_center(g::AMRGrid, i::Integer, ::Type{W}) where {W}
    exact = _amr_exact_rational(g.x0) + (2BigInt(i) - 1) * _amr_exact_rational(g.dx) / 2
    return W(exact)
end

@inline function _amr_compensated_add(total::T, correction::T, term::T) where {T}
    updated = total + term
    correction += abs(total) >= abs(term) ? (total - updated) + term : (term - updated) + total
    return updated, correction
end

@inline function _amr_two_sum(a::T, b::T) where {T}
    total = a + b
    recovered_b = total - a
    error = (a - (total - recovered_b)) + (b - recovered_b)
    return total, error
end

# The two components are exact and their addition is safe at the working scale.
# A result above the smallest normal (or overflow) preserves this single
# rounding under binary rescaling. Subnormals and the smallest normal share a
# rounding boundary where rescaling can round again and need certification.
@inline function _amr_round_exact_pair(high::T, low::T, power::Int) where {T}
    result = ldexp(high + low, power)
    iszero(low) && return result, true
    ordinary =
        T <: Union{Float16,Float32,Float64} ? abs(result) > floatmin(T) :
        T === BigFloat && !iszero(result)
    return result, ordinary
end

# Callers supply exact, safely scaled components under round-to-nearest.
@inline function _amr_round_scaled_sum3(a::T, b::T, c::T, power::Int) where {T}
    first, e1 = _amr_two_sum(a, b)
    high, e2 = _amr_two_sum(first, c)
    low, e3 = _amr_two_sum(e1, e2)
    if iszero(e3)
        result, certified = _amr_round_exact_pair(high, low, power)
        certified && return result, true
    end
    value, e4 = _amr_two_sum(high, low)
    error_hi, error_lo = _amr_two_sum(e4, e3)
    if iszero(error_lo)
        result, certified = _amr_round_exact_pair(value, error_hi, power)
        certified && return result, true
    end

    result = ldexp(value, power)
    # The exact normalized result is value + error_hi + error_lo. Reversible
    # rescaling and an error strictly inside both rounding midpoints certify it.
    if isfinite(result) && ldexp(result, -power) == value
        radius = min(value - prevfloat(value), nextfloat(value) - value) / T(2)
        abs(error_hi) + abs(error_lo) < radius && return result, true
    end
    return result, false
end

@inline function _amr_compensated_sum4(a::T, b::T, c::T, d::T) where {T}
    total = zero(T)
    correction = zero(T)
    total, correction = _amr_compensated_add(total, correction, a)
    total, correction = _amr_compensated_add(total, correction, b)
    total, correction = _amr_compensated_add(total, correction, c)
    total, correction = _amr_compensated_add(total, correction, d)
    return total + correction
end

@inline function _amr_product_parts(a::T, b::T) where {T<:AbstractFloat}
    ma, ea = frexp(a)
    mb, eb = frexp(b)
    product = ma * mb
    error = fma(ma, mb, -product)
    return product, error, ea + eb
end

@inline function _amr_scaled_product(a::T, b::T) where {T<:AbstractFloat}
    (iszero(a) || iszero(b)) && return zero(T)
    product, error, exponent = _amr_product_parts(a, b)
    return ldexp(product + error, exponent)
end

# Evaluate `a*b - c*d` after aligning two exact product expansions. This avoids
# both overflow in a representable difference and underflow of either product.
@inline function _amr_scaled_product_difference(a::T, b::T, c::T, d::T) where {T<:AbstractFloat}
    left_zero = iszero(a) || iszero(b)
    right_zero = iszero(c) || iszero(d)
    left_zero && right_zero && return zero(T)
    left_zero && return -_amr_scaled_product(c, d)
    right_zero && return _amr_scaled_product(a, b)

    left, left_error, left_exponent = _amr_product_parts(a, b)
    right, right_error, right_exponent = _amr_product_parts(c, d)
    common_exponent = max(left_exponent, right_exponent)
    left_shift = left_exponent - common_exponent
    right_shift = right_exponent - common_exponent
    difference = _amr_compensated_sum4(
        ldexp(left, left_shift),
        -ldexp(right, right_shift),
        ldexp(left_error, left_shift),
        -ldexp(right_error, right_shift),
    )
    return ldexp(difference, common_exponent)
end

# Widen integer inputs before forming the second difference. Int128 and UInt128
# widen to BigInt, so the public refinement criterion cannot wrap at any width.
@inline function _amr_second_difference(u::AbstractVector{T}, i::Int) where {T<:Integer}
    W = widen(signed(T))
    return W(u[i+1]) - W(2) * W(u[i]) + W(u[i-1])
end

@inline function _amr_normalisation_loses_precision(
    value::T,
    normalised::T,
) where {T<:Union{Float16,Float32,Float64}}
    return !iszero(value) && (iszero(normalised) || issubnormal(normalised))
end

@inline _amr_normalisation_loses_precision(value, normalised) = !iszero(value) && iszero(normalised)

@inline function _amr_normalisation_loses_precision(value::BigFloat, normalised::BigFloat)
    return precision(value) > precision(BigFloat) || (!iszero(value) && iszero(normalised))
end

_amr_exact_rational(value) = Rational{BigInt}(value)
function _amr_exact_rational(value::BigFloat)
    # Julia's Rational conversion scales in the current BigFloat precision.
    # Preserve stored input bits even when the caller has lowered that precision.
    return setprecision(BigFloat, precision(value)) do
        Rational{BigInt}(value)
    end
end

# Rational{BigInt} exactly represents every finite binary floating-point input.
# This cold fallback handles scaling loss and uncertified final rounding.
@noinline function _amr_exact_second_difference(a::T, b::T, c::T) where {T<:AbstractFloat}
    ar = _amr_exact_rational(a)
    br = _amr_exact_rational(b)
    cr = _amr_exact_rational(c)
    exact = if signbit(a) != signbit(c)
        (ar + cr) - 2br
    elseif !iszero(b) && signbit(a) == signbit(b)
        (ar - br) + (cr - br)
    else
        (ar + cr) - 2br
    end
    return T(exact)
end

@inline function _amr_second_difference(u::AbstractVector{T}, i::Int) where {T<:AbstractFloat}
    a, b, c = u[i-1], u[i], u[i+1]
    scale = max(abs(a), abs(b), abs(c))
    iszero(scale) && return zero(T)
    isfinite(scale) || return c - T(2) * b + a
    # TwoSum and the midpoint certificate below require round-to-nearest.
    if T === BigFloat && rounding(BigFloat) != RoundNearest
        return _amr_exact_second_difference(a, b, c)
    end

    # A general divisor rounds each input independently and can turn an exactly
    # affine profile into nonzero curvature. A power of two preserves the inputs
    # unless normalization reaches the subnormal range, checked below.
    _, power = frexp(scale)
    an, bn, cn = ldexp(a, -power), ldexp(b, -power), ldexp(c, -power)
    if _amr_normalisation_loses_precision(a, an) ||
       _amr_normalisation_loses_precision(b, bn) ||
       _amr_normalisation_loses_precision(c, cn)
        return _amr_exact_second_difference(a, b, c)
    end
    result, certified = _amr_round_scaled_sum3(an, -T(2) * bn, cn, power)
    return certified ? result : _amr_exact_second_difference(a, b, c)
end

@inline _amr_second_difference(u::AbstractVector, i::Int) = u[i+1] - 2u[i] + u[i-1]

@noinline function _amr_exact_reconstruction(a::T, b::T, c::T, shift::Int) where {T}
    exact =
        _amr_exact_rational(a) +
        (_amr_exact_rational(b) - _amr_exact_rational(c)) / (BigInt(1) << -shift)
    return T(exact)
end

# Reconstruct both children from the same unrounded slope. Rounding the slope
# first can erase a representable child when it cancels against the cell mean.
@inline function _amr_reconstruction_pair(a::T, b::T, c::T, shift::Int) where {T}
    scale = max(abs(a), abs(b), abs(c))
    iszero(scale) && return a - zero(T), a + zero(T)
    if !isfinite(scale)
        adjustment = ldexp(c - b, shift)
        return a - adjustment, a + adjustment
    end
    if T === BigFloat && rounding(BigFloat) != RoundNearest
        return _amr_exact_reconstruction(a, b, c, shift), _amr_exact_reconstruction(a, c, b, shift)
    end

    fraction = ldexp(one(T), shift)
    inverse_fraction = inv(fraction)
    bn, cn = b * fraction, c * fraction
    # With |a,b,c| <= max/2 and fraction <= 1/4, all intermediate sums are
    # bounded by 3max/4. Reversible scaling proves no input bits were discarded.
    direct =
        T <: Union{Float16,Float32,Float64} &&
        scale <= floatmax(T) / T(2) &&
        bn * inverse_fraction == b &&
        cn * inverse_fraction == c
    if direct
        an, power = a, 0
    else
        _, power = frexp(scale)
        an, bn, cn = ldexp(a, -power), ldexp(b, shift - power), ldexp(c, shift - power)
        if _amr_normalisation_loses_precision(a, an) ||
           _amr_normalisation_loses_precision(b, bn) ||
           _amr_normalisation_loses_precision(c, cn)
            return _amr_exact_reconstruction(a, b, c, shift),
            _amr_exact_reconstruction(a, c, b, shift)
        end
    end
    left, left_certified = _amr_round_scaled_sum3(an, bn, -cn, power)
    right, right_certified = _amr_round_scaled_sum3(an, cn, -bn, power)
    return (
        left_certified ? left : _amr_exact_reconstruction(a, b, c, shift),
        right_certified ? right : _amr_exact_reconstruction(a, c, b, shift),
    )
end

# `(a+b)/2` overflows for two large equal children, whereas evaluating each half
# separately underflows for two minimum subnormals.
@inline function _amr_average(a::T, b::T) where {T<:AbstractFloat}
    if isfinite(a) && isfinite(b)
        if T === BigFloat &&
           (precision(a) > precision(BigFloat) || precision(b) > precision(BigFloat))
            # Normalization would round stored input bits before cancellation.
            # Form the complete dyadic average before rounding to working precision.
            return T((_amr_exact_rational(a) + _amr_exact_rational(b)) / 2)
        end
        half = inv(T(2))
        return _amr_scaled_product_difference(a, half, -b, half)
    end
    return (a + b) / T(2)
end

"""
    refine_flags(u, threshold) -> BitVector

Tag cells for refinement where the undivided second difference exceeds
`threshold`: `abs(u[i+1] - 2u[i] + u[i-1]) > threshold`. Boundary cells are
never tagged.
"""
function refine_flags(u::AbstractVector{T}, threshold::Real) where {T}
    isfinite(threshold) && threshold >= 0 ||
        throw(ArgumentError("threshold must be finite and >= 0"))
    n = length(u)
    flags = falses(n)
    @inbounds for i = 2:n-1
        if abs(_amr_second_difference(u, i)) > threshold
            flags[i] = true
        end
    end
    return flags
end

function _validate_refinement_pair(fine::AMRGrid{T}, coarse::AMRGrid{T}) where {T}
    n = ncells(coarse)
    ncells(fine) == 2n || throw(DimensionMismatch("fine must have 2*$(n) = $(2n) cells"))
    fine.x0 == coarse.x0 || throw(ArgumentError("fine and coarse grids must have the same x0"))
    expected_dx = coarse.dx / T(2)
    fine.dx == expected_dx || throw(
        ArgumentError("fine dx must equal coarse dx / 2 (expected $expected_dx, got $(fine.dx))"),
    )
    return n
end

"""
    prolong!(fine, coarse) -> fine

Piecewise-linear prolongation from a coarse grid to a 2x-refined grid. Each
coarse cell `i` maps to fine cells `2i-1, 2i`; endpoint slopes are one-sided.
"""
function prolong!(fine::AMRGrid{T}, coarse::AMRGrid{T}) where {T}
    n = _validate_refinement_pair(fine, coarse)
    u = coarse.u
    if n == 1
        fine.u[1] = u[1] - zero(T)
        fine.u[2] = u[1] + zero(T)
        return fine
    end
    @inbounds for i = 1:n
        left = i == 1 ? u[i] : u[i-1]
        right = i == n ? u[i] : u[i+1]
        shift = i == 1 || i == n ? -2 : -3
        fine.u[2i-1], fine.u[2i] = _amr_reconstruction_pair(u[i], left, right, shift)
    end
    return fine
end

"""
    restrict!(coarse, fine) -> coarse

Conservative restriction from a 2x-refined grid to a coarse grid:
`coarse[i] = (fine[2i-1] + fine[2i]) / 2`.
"""
function restrict!(coarse::AMRGrid{T}, fine::AMRGrid{T}) where {T}
    n = _validate_refinement_pair(fine, coarse)
    @inbounds for i = 1:n
        coarse.u[i] = _amr_average(fine.u[2i-1], fine.u[2i])
    end
    return coarse
end

"""
    refine(coarse) -> AMRGrid

Return a uniformly 2x-refined level obtained by prolonging `coarse`.
"""
function refine(coarse::AMRGrid{T}) where {T}
    n = ncells(coarse)
    n <= typemax(Int) ÷ 2 || throw(OverflowError("refined AMR cell count does not fit in Int"))
    fine = AMRGrid(Vector{T}(undef, 2n), coarse.dx / 2; x0 = coarse.x0)
    return prolong!(fine, coarse)
end
