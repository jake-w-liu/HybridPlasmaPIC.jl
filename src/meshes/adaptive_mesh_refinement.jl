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
cell_center(g::AMRGrid, i::Integer) = g.x0 + (i - oftype(g.dx, 0.5)) * g.dx

@inline function _amr_compensated_add(total::T, correction::T, term::T) where {T}
    updated = total + term
    correction +=
        abs(total) >= abs(term) ?
        (total - updated) + term :
        (term - updated) + total
    return updated, correction
end

@inline function _amr_compensated_sum3(a::T, b::T, c::T) where {T}
    total = zero(T)
    correction = zero(T)
    total, correction = _amr_compensated_add(total, correction, a)
    total, correction = _amr_compensated_add(total, correction, b)
    total, correction = _amr_compensated_add(total, correction, c)
    return total + correction
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
@inline function _amr_scaled_product_difference(
    a::T,
    b::T,
    c::T,
    d::T,
) where {T<:AbstractFloat}
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

@inline function _amr_normalisation_loses_precision(value::T, normalised::T) where {
    T<:Union{Float16,Float32,Float64},
}
    return !iszero(value) && (iszero(normalised) || issubnormal(normalised))
end

@inline _amr_normalisation_loses_precision(value, normalised) =
    !iszero(value) && iszero(normalised)

# Rational{BigInt} exactly represents every finite binary floating-point input.
# This cold fallback is used only when a shared floating scale erases information.
@noinline function _amr_exact_second_difference(a::T, b::T, c::T) where {T<:AbstractFloat}
    ar = Rational{BigInt}(a)
    br = Rational{BigInt}(b)
    cr = Rational{BigInt}(c)
    exact = if signbit(a) != signbit(c)
        (ar + cr) - 2br
    elseif !iszero(b) && signbit(a) == signbit(b)
        (ar - br) + (cr - br)
    else
        (ar + cr) - 2br
    end
    return T(exact)
end

@inline function _amr_second_difference(
    u::AbstractVector{T},
    i::Int,
) where {T<:AbstractFloat}
    a, b, c = u[i-1], u[i], u[i+1]
    scale = max(abs(a), abs(b), abs(c))
    iszero(scale) && return zero(T)
    isfinite(scale) || return c - T(2) * b + a

    an, bn, cn = a / scale, b / scale, c / scale
    if _amr_normalisation_loses_precision(a, an) ||
       _amr_normalisation_loses_precision(b, bn) ||
       _amr_normalisation_loses_precision(c, cn)
        return _amr_exact_second_difference(a, b, c)
    end
    scaled = _amr_compensated_sum3(an, -T(2) * bn, cn)
    return _amr_scaled_product(scale, scaled)
end

@inline _amr_second_difference(u::AbstractVector, i::Int) =
    u[i+1] - 2u[i] + u[i-1]

# Form `(a-b)*fraction` without overflowing the subtraction. The fractions used
# by prolongation are exact powers of two.
@inline function _amr_scaled_difference(a::T, b::T, fraction::T) where {T<:AbstractFloat}
    if isfinite(a) && isfinite(b)
        return _amr_scaled_product_difference(a, fraction, b, fraction)
    end
    return (a - b) * fraction
end

# `(a+b)/2` overflows for two large equal children, whereas evaluating each half
# separately underflows for two minimum subnormals.
@inline function _amr_average(a::T, b::T) where {T<:AbstractFloat}
    if isfinite(a) && isfinite(b)
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
    ncells(fine) == 2n ||
        throw(DimensionMismatch("fine must have 2*$(n) = $(2n) cells"))
    fine.x0 == coarse.x0 ||
        throw(ArgumentError("fine and coarse grids must have the same x0"))
    expected_dx = coarse.dx / T(2)
    fine.dx == expected_dx || throw(
        ArgumentError(
            "fine dx must equal coarse dx / 2 (expected $expected_dx, got $(fine.dx))",
        ),
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
    quarter = inv(T(4))
    eighth = inv(T(8))
    @inbounds for i = 1:n
        adjustment =
            n == 1 ? zero(T) :
            i == 1 ? _amr_scaled_difference(u[2], u[1], quarter) :
            i == n ? _amr_scaled_difference(u[n], u[n-1], quarter) :
            _amr_scaled_difference(u[i+1], u[i-1], eighth)
        fine.u[2i-1] = u[i] - adjustment
        fine.u[2i] = u[i] + adjustment
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
    n <= typemax(Int) ÷ 2 ||
        throw(OverflowError("refined AMR cell count does not fit in Int"))
    fine = AMRGrid(Vector{T}(undef, 2n), coarse.dx / 2; x0 = coarse.x0)
    return prolong!(fine, coarse)
end
