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
struct AMRGrid{T}
    u::Vector{T}
    dx::T
    x0::T
end

AMRGrid(u::Vector{T}, dx::Real; x0::Real = 0.0) where {T} = AMRGrid{T}(u, T(dx), T(x0))

ncells(g::AMRGrid) = length(g.u)
effective_resolution(g::AMRGrid) = g.dx
cell_center(g::AMRGrid, i::Integer) = g.x0 + (i - oftype(g.dx, 0.5)) * g.dx

"""
    refine_flags(u, threshold) -> BitVector

Tag cells for refinement where the undivided second difference exceeds
`threshold`: `abs(u[i+1] - 2u[i] + u[i-1]) > threshold`. Boundary cells are
never tagged.
"""
function refine_flags(u::AbstractVector{T}, threshold::Real) where {T}
    threshold >= 0 || throw(ArgumentError("threshold must be >= 0"))
    n = length(u)
    flags = falses(n)
    thr = T(threshold)
    @inbounds for i = 2:n-1
        if abs(u[i+1] - 2u[i] + u[i-1]) > thr
            flags[i] = true
        end
    end
    return flags
end

"""
    prolong!(fine, coarse) -> fine

Piecewise-linear prolongation from a coarse grid to a 2x-refined grid. Each
coarse cell `i` maps to fine cells `2i-1, 2i`; endpoint slopes are one-sided.
"""
function prolong!(fine::AMRGrid{T}, coarse::AMRGrid{T}) where {T}
    n = ncells(coarse)
    ncells(fine) == 2n || throw(DimensionMismatch("fine must have 2*$(n) = $(2n) cells"))
    u = coarse.u
    @inbounds for i = 1:n
        s =
            n == 1 ? zero(T) :
            i == 1 ? (u[2] - u[1]) : i == n ? (u[n] - u[n-1]) : (u[i+1] - u[i-1]) / 2
        fine.u[2i-1] = u[i] - s / 4
        fine.u[2i] = u[i] + s / 4
    end
    return fine
end

"""
    restrict!(coarse, fine) -> coarse

Conservative restriction from a 2x-refined grid to a coarse grid:
`coarse[i] = (fine[2i-1] + fine[2i]) / 2`.
"""
function restrict!(coarse::AMRGrid{T}, fine::AMRGrid{T}) where {T}
    n = ncells(coarse)
    ncells(fine) == 2n || throw(DimensionMismatch("fine must have 2*$(n) = $(2n) cells"))
    @inbounds for i = 1:n
        coarse.u[i] = (fine.u[2i-1] + fine.u[2i]) / 2
    end
    return coarse
end

"""
    refine(coarse) -> AMRGrid

Return a uniformly 2x-refined level obtained by prolonging `coarse`.
"""
function refine(coarse::AMRGrid{T}) where {T}
    n = ncells(coarse)
    fine = AMRGrid(Vector{T}(undef, 2n), coarse.dx / 2; x0 = coarse.x0)
    return prolong!(fine, coarse)
end
