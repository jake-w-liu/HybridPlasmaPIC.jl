# deposit.jl — particle ⇄ grid coupling on the collocated periodic mesh.
#
# Deposition and gathering share ONE shape function (`_stencil1d`), so they are
# exact transposes of each other (the deposit/gather adjoint identity holds to
# roundoff). Grid nodes sit at (i−1)·dx, i=1..n, periodic; the FourierGrid
# supplies n and dx (no FFT needed here, just the mesh geometry).

"Tensor-product particle shapes. Touched nodes per axis: NGP 1, CIC 2, TSC 3."
abstract type ShapeFunction end
struct NGP <: ShapeFunction end
struct CIC <: ShapeFunction end
struct TSC <: ShapeFunction end

"""
    YeeTriangularCloud(half_width)

Optional finite physical-width triangular charge cloud for the closed 1D Yee
charge/current/voltage operators and [`yee_longitudinal_field`](@ref). The
positive finite Float32/Float64 half-width is retained without conversion to
grid precision. Hold it fixed under mesh refinement. This is a separate shape
from the default point/CIC model, with particle-bound endpoint tail support;
it is not a periodic deposit/gather or device shape. See the bounded-particle
contract for its distinct wall response and remaining coupled-model work.
"""
struct YeeTriangularCloud{T<:Union{Float32,Float64}}
    half_width::T
    function YeeTriangularCloud(a::T) where {T<:Union{Float32,Float64}}
        isfinite(a) && _exact_rational(a) > 0 ||
            throw(ArgumentError("Yee cloud half-width must be finite and positive"))
        return new{T}(a)
    end
end

# Exact primitive R and cumulative P of the triangular density. Every cloud
# operator uses this same cumulative A and derivative g, in physical geometry.
@inline function _yee_cloud_primitives(t, a)
    if t <= -a
        return zero(t), zero(t)
    elseif t <= 0
        r = (t + a)^2 / (2a^2)
        return r * (t + a) / 3, r
    elseif t < a
        r = (a - t)^2 / (2a^2)
        return t + r * (a - t) / 3, one(t) - r
    end
    return t, one(t)
end

@inline function _yee_cloud_edge(x, lower, h, a)
    left = _yee_cloud_primitives(x - lower, a)
    right = _yee_cloud_primitives(x - lower - h, a)
    return (left[1] - right[1]) / h, (left[2] - right[2]) / h
end

# One further primitive U'=R supplies exact nodal path averages. Existing
# A/g keep their original arithmetic; inputs here are original rationals.
@inline function _yee_cloud_path_primitive(t, a)
    positive(x) = max(zero(x), x)
    return (positive(t + a)^4 - 2 * positive(t)^4 + positive(t - a)^4) / (24a^2)
end
@inline function _yee_cloud_integrated_edge(x, lower, h, a)
    return (_yee_cloud_path_primitive(x - lower, a) -
        _yee_cloud_path_primitive(x - lower - h, a)) / h
end

# Outside this exact span both endpoint cumulatives equal the same saturated
# constant and their derivatives vanish. Clip before Int conversion: widths can
# exceed the domain, while positive a keeps the first retained edge below n.
@inline function _yee_cloud_edge_span(x, y, h, a, n)
    lower = max(zero(h), min(x, y) - a)
    upper = min(n * h, max(x, y) + a)
    return floor(Int, lower / h):min(n - 1, floor(Int, upper / h))
end

@inline width(::NGP) = 1
@inline width(::CIC) = 2
@inline width(::TSC) = 3

@inline _periodic_cell_position(x, L, dx) = mod(x, L) / dx

Base.@propagate_inbounds @inline function _particle_cell_position(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    d::Int,
    p::Int,
) where {D,T}
    x = ps.x[d][p]
    isfinite(x) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
    return _periodic_cell_position(x, g.L[d], g.dx[d])
end

function _validate_particle_positions(ps::ParticleSet{D}) where {D}
    @inbounds for p in eachindex(ps.weight), d = 1:D
        isfinite(ps.x[d][p]) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
    end
    return nothing
end

# 1-D stencil at fractional cell position s: returns (base, weights), where the
# touched nodes are base, base+1, … (0-based, to be wrapped mod n). Weights sum
# to 1 (partition of unity) for every s.
# Round before integer conversion so CPU and device stencils share the same
# upward-tie policy without rounding a just-below-half coordinate onto the tie.
@inline _ngp_node(s) = round(s, RoundNearestTiesUp)

@inline function _stencil1d(::NGP, s::T) where {T}
    c = Int(_ngp_node(s))
    return (c, (one(T),))
end
@inline function _stencil1d(::CIC, s::T) where {T}
    i0 = floor(Int, s)
    f = s - i0
    return (i0, (one(T) - f, f))
end
@inline function _stencil1d(::TSC, s::T) where {T}
    c = round(Int, s)
    δ = s - c
    half = T(0.5)
    return (c - 1, (half * (half - δ)^2, T(0.75) - δ^2, half * (half + δ)^2))
end

# Bounded CIC uses the closed physical interval [0,L]. Its exact coordinate is
# N*x/L, independently of the rounded stored dx used by the Maxwell field rows.
# These endpoint shapes do not select a particle force or trajectory.
@inline function _yee_cic_fraction(terms::Tuple, denominator::Float64)
    value = 0.0
    for term in _eft_expansion(terms)
        value += term
    end
    error = _eft_abs_sum_bound((terms..., -value))
    iszero(value) && iszero(error) && return 0.0, true
    threshold = max(0.0, prevfloat(0x1p-53 * abs(value)))
    value > error && (iszero(error) || error <= threshold) || return 0.0, false
    fraction = value / denominator
    return fraction, floatmin(Float64) <= fraction <= 1.0
end

# Validation admits IEEE coordinates even when a collection or vector load is
# abstractly typed. Express that scalar domain before conversion/scaling.
@inline function _yee_cic_stencil(x::Union{Float32,Float64}, g::YeeGrid{1}, ::Val{false})
    n, length = g.n[1], Float64(g.L[1])
    iszero(x) && return 0, (1.0, 0.0), true
    x == g.L[1] && return n - 1, (0.0, 1.0), true
    n <= 0x1p53 || return 0, (0.0, 0.0), false
    scale = exponent(length)
    normalized_length = ldexp(length, -scale)
    normalized_x = ldexp(Float64(x), -scale)
    ldexp(normalized_x, scale) == x || return 0, (0.0, 0.0), false
    product, product_error = _eft_two_product(normalized_x, Float64(n))
    # Two 53-bit factors have at most 106 bits. At this floor their last
    # possible product/residual bit is normal, so the FMA product is exact.
    product >= 0x1p-916 || return 0, (0.0, 0.0), false
    cell, _ = _stencil1d(CIC(), product / normalized_length)
    0 <= cell <= n || return 0, (0.0, 0.0), false
    cell == n && (cell -= 1)
    lower, lower_error = _eft_two_product(Float64(cell), normalized_length)
    upper, upper_error = _eft_two_product(Float64(cell + 1), normalized_length)
    right, right_ok =
        _yee_cic_fraction((product, product_error, -lower, -lower_error), normalized_length)
    left, left_ok =
        _yee_cic_fraction((upper, upper_error, -product, -product_error), normalized_length)
    return cell, (left, right), left_ok && right_ok
end

@inline function _yee_cic_stencil(x::Union{Float32,Float64}, g::YeeGrid{1}, ::Val{true})
    n = g.n[1]
    s = _exact_rational(n) * _exact_rational(x) / _exact_rational(g.L[1])
    if s == n
        return n - 1, (zero(s), one(s)), true
    end
    cell, weights = _stencil1d(CIC(), s)
    return cell, weights, true
end

@inline function _yee_cic_displacement(a, b, g::YeeGrid{1}, ::Val{false})
    scale = exponent(Float64(g.L[1]))
    difference = Float64(b) - Float64(a)
    normalized = ldexp(difference, -scale)
    ldexp(normalized, scale) == difference || return 0.0, false
    product = normalized * Float64(g.n[1])
    coefficient = product / ldexp(Float64(g.L[1]), -scale)
    return coefficient, isfinite(coefficient) && floatmin(Float64) <= coefficient <= 1.0
end

@inline function _yee_cic_displacement(a, b, g::YeeGrid{1}, ::Val{true})
    value =
        _exact_rational(g.n[1]) * (_exact_rational(b) - _exact_rational(a)) /
        _exact_rational(g.L[1])
    return value, true
end

# Nonnegative normalized cell overlaps for an ordered path a <= b. Interior
# crossed cells have weight one; the two returned weights own its endpoint
# cells. When both endpoints share a cell, only the first weight is used.
@inline function _yee_cic_path_stencil(a, b, g::YeeGrid{1}, exact)
    lower, left, lower_ok = _yee_cic_stencil(a, g, exact)
    upper, right, upper_ok = _yee_cic_stencil(b, g, exact)
    lower_ok && upper_ok || return lower, upper, (zero(left[1]), zero(right[2])), false
    if lower == upper
        displacement, admitted = _yee_cic_displacement(a, b, g, exact)
        return lower, upper, (displacement, zero(displacement)), admitted
    end
    return lower, upper, (left[1], right[2]), true
end
