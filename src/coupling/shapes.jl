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

@inline width(::NGP) = 1
@inline width(::CIC) = 2
@inline width(::TSC) = 3

Base.@propagate_inbounds @inline function _particle_cell_position(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    d::Int,
    p::Int,
) where {D,T}
    x = ps.x[d][p]
    isfinite(x) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
    return x / g.dx[d]
end

# 1-D stencil at fractional cell position s: returns (base, weights), where the
# touched nodes are base, base+1, … (0-based, to be wrapped mod n). Weights sum
# to 1 (partition of unity) for every s.
@inline function _stencil1d(::NGP, s::T) where {T}
    c = round(Int, s)
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
