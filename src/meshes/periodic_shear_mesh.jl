# One fixed nonsingular periodic shear, not a toroidal or arbitrary-map backend.
struct PeriodicShearGrid{T,G}
    logical::G
    amplitude::T
    tau::T
    _s::Vector{T}
end

function _shear_round(::Type{T}, value, context; preserve_nonzero = false) where {T}
    out = _from_exact(T, value)
    isfinite(out) && (!preserve_nonzero || iszero(value) || !iszero(out)) ||
        throw(ArgumentError("$context exceeds its materialized range"))
    return out
end

# Keep the original ordinary order; replay original factors on a range event.
@inline function _shear_linear(a::T, b::T, c::T) where {T}
    all(isfinite, (a, b, c)) ||
        throw(ArgumentError("shear transform requires finite original operands"))
    product = b * c
    value = a + product
    if _deposition_product_lost(b, c, product) ||
       _deposition_range_lost(value) ||
       (iszero(value) && !iszero(a) && !iszero(product))
        Q = _exact_rational
        return _shear_round(T, Q(a) + Q(b) * Q(c), "shear component")
    end
    return value
end

function _shear_fourier_grid(g::FourierGrid{3,T}) where {T}
    all(n -> n >= 3 && isodd(n), g.n) ||
        throw(ArgumentError("periodic shear requires odd grid sizes of at least three"))
    expected = SpectralOperators._fourier_metadata(g.n, g.L)
    all(
        name -> isequal(getproperty(g, name), getproperty(expected, name)),
        (:dx, :midx, :kfull, :kvec, :ik),
    ) || throw(ArgumentError("periodic shear Fourier metadata changed"))
    all(a -> a isa Array{Complex{T},3} && size(a) == g.n, (g.cbuf, g.tbuf, g.abuf, g.sbuf)) ||
        throw(ArgumentError("periodic shear requires matching ordinary Fourier buffers"))
    return nothing
end

function _shear_slope(mesh::PeriodicShearGrid{T}, y::T) where {T}
    L = mesh.logical.L[2]
    isfinite(y) || throw(ArgumentError("shear coordinate must be finite"))
    wrapped = _periodic_coordinate(y, zero(T), L, L)
    _, cosine = sincospi(T(2) * (wrapped / L))
    Q = _exact_rational
    return _shear_round(T, Q(mesh.amplitude) * Q(mesh.tau) * Q(cosine) / Q(L), "shear metric")
end

"""
    PeriodicShearGrid(logical; amplitude)

Fixed periodic map x=(xi1+A*sin(2pi*xi2/L2),xi2,xi3), with J=1 and a
nonorthogonal basis. Serial CPU Float32/64, odd sizes >=3. No toroidal axis/wall.
"""
function PeriodicShearGrid(g::FourierGrid{3,T}; amplitude) where {T<:Union{Float32,Float64}}
    _boris_accuracy_mode() ||
        throw(ArgumentError("periodic shear requires nearest rounding and gradual underflow"))
    _shear_fourier_grid(g)
    A = _require_finite_real("shear amplitude", amplitude, T)
    (iszero(amplitude) || !iszero(A)) ||
        throw(ArgumentError("nonzero shear amplitude must survive conversion"))
    mesh = PeriodicShearGrid{T,typeof(g)}(g, A, T(2pi), zeros(T, g.n[2]))
    for j in eachindex(mesh._s)
        mesh._s[j] = _shear_slope(mesh, T(j - 1) * g.dx[2])
    end
    return mesh
end

gridsize(mesh::PeriodicShearGrid) = mesh.logical.n
function _shear_geometry(mesh::PeriodicShearGrid{T}) where {T}
    _shear_fourier_grid(mesh.logical)
    isfinite(mesh.amplitude) && isequal(mesh.tau, T(2pi)) && length(mesh._s) == mesh.logical.n[2] ||
        throw(ArgumentError("fixed shear geometry changed"))
    for j in eachindex(mesh._s)
        isequal(mesh._s[j], _shear_slope(mesh, T(j - 1) * mesh.logical.dx[2])) ||
            throw(ArgumentError("fixed shear metric cache changed"))
    end
    return nothing
end

function _shear_map(mesh::PeriodicShearGrid{T}, coordinates, inverse::Bool) where {T}
    _boris_accuracy_mode() ||
        throw(ArgumentError("periodic shear requires nearest rounding and gradual underflow"))
    _shear_geometry(mesh)
    point = ntuple(d -> _require_finite_real("coordinate", coordinates[d], T), 3)
    wrapped = ntuple(
        d -> _periodic_coordinate(point[d], zero(T), mesh.logical.L[d], mesh.logical.L[d]),
        3,
    )
    sine, _ = sincospi(T(2) * (wrapped[2] / mesh.logical.L[2]))
    first = _shear_linear(wrapped[1], inverse ? -mesh.amplitude : mesh.amplitude, T(sine))
    return (
        _periodic_coordinate(first, zero(T), mesh.logical.L[1], mesh.logical.L[1]),
        wrapped[2],
        wrapped[3],
    )
end

"Map periodic logical coordinates to physical Cartesian coordinates."
to_cartesian(mesh::PeriodicShearGrid, coordinates::NTuple{3,Real}) =
    _shear_map(mesh, coordinates, false)
"Invert the fixed periodic shear map; return half-open logical representatives."
to_logical(mesh::PeriodicShearGrid, coordinates::NTuple{3,Real}) =
    _shear_map(mesh, coordinates, true)
