# cartesian_mesh.jl (§5.3 meshes/) — the concrete periodic Fourier and mixed
# SBP-x/Fourier-y,z operators live in the sibling SpectralOperators.jl package.
# Both are Cartesian. The bounded Maxwell solver owns the staggered mesh below.

"""
    YeeGrid(n, lengths; boundaries=ntuple(_ -> :pec, length(n)))

Uniform Cartesian Yee mesh in one to three spatial dimensions. `n` counts cells;
each axis has either two homogeneous perfectly conducting (`:pec`) walls or is
`:periodic`. Coordinates start at zero. Electric components live on edges and
magnetic components on faces, including all three vector components in 1D/2D.
Periodic nodes omit the duplicate endpoint; conducting axes retain both walls.
`lengths` must promote to Float32 or Float64.
"""
struct YeeGrid{D,T<:Union{Float32,Float64}}
    n::NTuple{D,Int}
    L::NTuple{D,T}
    dx::NTuple{D,T}
    periodic::NTuple{D,Bool}
    function YeeGrid(
        n::NTuple{D,Int},
        L::NTuple{D,T},
        periodic::NTuple{D,Bool},
    ) where {D,T<:Union{Float32,Float64}}
        _check_spatial_dimension(D)
        all(>=(2), n) || throw(ArgumentError("YeeGrid requires at least two cells per axis"))
        all(d -> periodic[d] || n[d] < typemax(Int), 1:D) ||
            throw(ArgumentError("conducting-axis node counts must fit Int"))
        all(x -> isfinite(x) && x > 0, L) ||
            throw(ArgumentError("grid lengths must be finite and positive"))
        dx = ntuple(d -> L[d] / n[d], D)
        all(x -> isfinite(x) && x > 0 && isfinite(inv(x)), dx) ||
            throw(ArgumentError("grid spacings and their reciprocals must be finite and positive"))
        new{D,T}(n, L, dx, periodic)
    end
end

function YeeGrid(
    n::NTuple{D,Integer},
    lengths::NTuple{D,Real};
    boundaries = ntuple(_ -> :pec, D),
) where {D}
    _check_spatial_dimension(D)
    length(boundaries) == D || throw(DimensionMismatch("one boundary pair is required per axis"))
    all(b -> b === :pec || b === :periodic, boundaries) ||
        throw(ArgumentError("each boundary pair must be :pec or :periodic"))
    L = float.(promote(lengths...))
    eltype(L) <: Union{Float32,Float64} ||
        throw(ArgumentError("YeeGrid supports Float32 and Float64"))
    return YeeGrid(Int.(n), L, ntuple(d -> boundaries[d] === :periodic, D))
end

@inline _yee_node(electric::Bool, component::Int, axis::Int) =
    electric ? component != axis : component == axis
@inline _yee_shape(g::YeeGrid{D}, electric::Bool, component::Int) where {D} =
    ntuple(d -> g.n[d] + Int(_yee_node(electric, component, d) && !g.periodic[d]), D)

function _yee_electric(field::Symbol)
    field === :electric && return true
    field === :magnetic && return false
    throw(ArgumentError("field must be :electric or :magnetic"))
end

"""
    yee_coordinates(g, field, component)

Return one coordinate vector per axis for `field = :electric` or `:magnetic`
and vector component 1, 2, or 3. Components have different array shapes.
"""
function yee_coordinates(g::YeeGrid{D,T}, field::Symbol, component::Integer) where {D,T}
    electric = _yee_electric(field)
    1 <= component <= 3 || throw(ArgumentError("component must be 1, 2, or 3"))
    c = Int(component)
    sizes = _yee_shape(g, electric, c)
    return ntuple(D) do d
        offset = _yee_node(electric, c, d) ? zero(T) : T(0.5)
        [(T(i - 1) + offset) * g.dx[d] for i = 1:sizes[d]]
    end
end
