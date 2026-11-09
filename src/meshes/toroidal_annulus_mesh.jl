# The logical radial interval is an annulus, with normalized periodic turns.
# This does not reinterpret the cell-centred solid-axis ToroidalGrid.
struct ToroidalAnnulusGrid{T,O}
    op::O
    n::NTuple{3,Int}
    L::NTuple{3,T}
    dx::NTuple{3,T}
    R0::T
    r_inner::T
    r_outer::T
    tau::T
    r::Vector{T}
    R::Matrix{T}
    su::Vector{T}
    cu::Vector{T}
    sv::Vector{T}
    cv::Vector{T}
end

function _tah_round(::Type{T}, value, context; preserve_nonzero = false) where {T}
    result = _from_exact(T, value)
    isfinite(result) && (!preserve_nonzero || iszero(value) || !iszero(result)) ||
        throw(ArgumentError("$context is outside its materialized range"))
    return result
end

function _tah_trig(turn::T) where {T}
    isfinite(turn) || throw(ArgumentError("toroidal turn must be finite"))
    wrapped = _periodic_coordinate(turn, zero(T), one(T), one(T))
    return sincospi(T(2) * wrapped)
end

function _tah_basis(su, cu, sv, cv)
    Q = _exact_rational
    s, c, v, w = Q(su), Q(cu), Q(sv), Q(cv)
    return ((c * w, c * v, s), (-s * w, -s * v, c), (-v, w, zero(_CompositeQ)))
end
_tah_basis(mesh::ToroidalAnnulusGrid, I::CartesianIndex{3}) =
    _tah_basis(mesh.su[I[2]], mesh.cu[I[2]], mesh.sv[I[3]], mesh.cv[I[3]])

function _tah_factors(tau, r, R)
    Q = _exact_rational
    h = (one(_CompositeQ), Q(tau) * Q(r), Q(tau) * Q(R))
    return h, prod(h)
end
_tah_metric(mesh::ToroidalAnnulusGrid, I::CartesianIndex{3}) =
    _tah_factors(mesh.tau, mesh.r[I[1]], mesh.R[I[1], I[2]])
function _tah_volume(mesh::ToroidalAnnulusGrid, I::CartesianIndex{3})
    _, J = _tah_metric(mesh, I)
    return J * _exact_rational(mesh.op.sbp.H[I[1]]) *
        _exact_rational(mesh.dx[2]) * _exact_rational(mesh.dx[3])
end

"""
    ToroidalAnnulusGrid((Nr, Nu, Nv); R0, r_inner, r_outer, T=Float64)

Nonsingular circular toroidal annulus with an SBP radial interval and periodic
normalized turns. Positive physical volume has orientation sigma=-1 for this
coordinate order. This is not a solid-axis or material-wall mesh.
"""
function ToroidalAnnulusGrid(n::NTuple{3,Integer}; R0, r_inner, r_outer, T::Type = Float64)
    T in (Float32, Float64) || throw(ArgumentError("toroidal annulus requires Float32 or Float64"))
    _mixed_mode(T)
    all(k -> k >= 3, n) && all(isodd, Base.tail(n)) ||
        throw(ArgumentError("annular radial size must be >=3; periodic sizes must be odd and >=3"))
    major = _require_finite_positive_real("toroidal major radius", R0, T)
    inner = _require_finite_positive_real("annular inner radius", r_inner, T)
    outer = _require_finite_positive_real("annular outer radius", r_outer, T)
    inner < outer < major || throw(ArgumentError("require 0 < r_inner < r_outer < R0"))
    ell = _tah_round(T, _exact_rational(outer) - _exact_rational(inner), "annular width";
        preserve_nonzero = true)
    op = MixedShockOperators(n, (ell, one(T), one(T)))
    tau = T(2pi)
    r = T[_tah_round(T, _exact_rational(inner) + _exact_rational(x), "annular node")
        for x in op.coordinates[1]]
    all(i -> r[i] < r[i+1], 1:(length(r)-1)) && isequal(r[end], outer) ||
        throw(ArgumentError("annular radius nodes must remain distinct and retain the outer radius"))
    su, cu = Vector{T}(undef, op.n[2]), Vector{T}(undef, op.n[2])
    sv, cv = Vector{T}(undef, op.n[3]), Vector{T}(undef, op.n[3])
    for j in eachindex(su)
        su[j], cu[j] = _tah_trig(op.coordinates[2][j])
    end
    for k in eachindex(sv)
        sv[k], cv[k] = _tah_trig(op.coordinates[3][k])
    end
    R = Matrix{T}(undef, op.n[1], op.n[2])
    for j in axes(R, 2), i in axes(R, 1)
        R[i, j] = _tah_round(T, _exact_rational(major) +
            _exact_rational(r[i]) * _exact_rational(cu[j]), "annular major coordinate";
            preserve_nonzero = true)
        R[i, j] > zero(T) || throw(ArgumentError("annular major coordinate must be positive"))
    end
    mesh = ToroidalAnnulusGrid(op, op.n, op.L, op.dx, major, inner, outer, tau,
        r, R, su, cu, sv, cv)
    _tah_geometry(mesh)
    return mesh
end

gridsize(mesh::ToroidalAnnulusGrid) = mesh.n

function _tah_geometry(mesh::ToroidalAnnulusGrid{T}) where {T}
    _mixed_storage(mesh.op)
    mesh.n == mesh.op.n && isequal(mesh.L, mesh.op.L) && isequal(mesh.dx, mesh.op.dx) &&
        0 < mesh.r_inner < mesh.r_outer < mesh.R0 && isequal(mesh.tau, T(2pi)) &&
        isequal(mesh.L[1], _tah_round(T, _exact_rational(mesh.r_outer) -
            _exact_rational(mesh.r_inner), "annular width")) &&
        mesh.L[2] == one(T) && mesh.L[3] == one(T) ||
        throw(ArgumentError("annular immutable geometry changed"))
    length(mesh.r) == mesh.n[1] && size(mesh.R) == (mesh.n[1], mesh.n[2]) &&
        length(mesh.su) == length(mesh.cu) == mesh.n[2] &&
        length(mesh.sv) == length(mesh.cv) == mesh.n[3] ||
        throw(DimensionMismatch("annular metric cache sizes changed"))
    for i in eachindex(mesh.r)
        expected = _tah_round(T, _exact_rational(mesh.r_inner) +
            _exact_rational(mesh.op.coordinates[1][i]), "annular radius")
        isequal(mesh.r[i], expected) || throw(ArgumentError("annular radius cache changed"))
        i == 1 || mesh.r[i-1] < mesh.r[i] || throw(ArgumentError("annular nodes collapsed"))
    end
    isequal(mesh.r[end], mesh.r_outer) || throw(ArgumentError("annular outer node changed"))
    for j in eachindex(mesh.su)
        isequal((mesh.su[j], mesh.cu[j]), _tah_trig(mesh.op.coordinates[2][j])) ||
            throw(ArgumentError("annular poloidal cache changed"))
        for i in eachindex(mesh.r)
            expected = _tah_round(T, _exact_rational(mesh.R0) +
                _exact_rational(mesh.r[i]) * _exact_rational(mesh.cu[j]), "annular major coordinate")
            isequal(mesh.R[i, j], expected) && expected > zero(T) ||
                throw(ArgumentError("annular major-coordinate cache changed"))
        end
    end
    for k in eachindex(mesh.sv)
        isequal((mesh.sv[k], mesh.cv[k]), _tah_trig(mesh.op.coordinates[3][k])) ||
            throw(ArgumentError("annular toroidal cache changed"))
    end
    return nothing
end

function _tah_mesh_arrays(mesh::ToroidalAnnulusGrid)
    arrays = Any[mesh.r, mesh.R, mesh.su, mesh.cu, mesh.sv, mesh.cv,
        _mixed_operator_arrays(mesh.op)...]
    for work in mesh.op.transverse
        _visit_checkpoint_objects((work.plan, work.iplan)) do object
            object isa AbstractArray && !any(a -> a === object, arrays) && push!(arrays, object)
        end
    end
    return arrays
end

function _tah_point_geometry(mesh::ToroidalAnnulusGrid{T}, point) where {T}
    length(point) == 3 || throw(DimensionMismatch("annular logical point requires three coordinates"))
    xi, u, v = ntuple(d -> _require_finite_real("annular logical coordinate", point[d], T), 3)
    0 <= xi <= mesh.L[1] || throw(ArgumentError("annular logical radius is outside [0,ell]"))
    su, cu = _tah_trig(u)
    sv, cv = _tah_trig(v)
    r = _tah_round(T, _exact_rational(mesh.r_inner) + _exact_rational(xi), "annular particle radius")
    R = _tah_round(T, _exact_rational(mesh.R0) + _exact_rational(r) * _exact_rational(cu),
        "annular particle major coordinate"; preserve_nonzero = true)
    r > 0 && R > 0 || throw(ArgumentError("annular point is singular"))
    h, J = _tah_factors(mesh.tau, r, R)
    return (r = r, R = R, basis = _tah_basis(su, cu, sv, cv), h = h, J = J,
        trig = (su, cu, sv, cv))
end

function to_cartesian(mesh::ToroidalAnnulusGrid{T}, point::NTuple{3,Real}) where {T}
    _tah_geometry(mesh)
    geo = _tah_point_geometry(mesh, point)
    su, cu, sv, cv = map(_exact_rational, geo.trig)
    values = _toroidal_cartesian_point(_exact_rational(mesh.R0), _exact_rational(geo.r), su, cu, sv, cv)
    return map(x -> _tah_round(T, x, "annular Cartesian coordinate"), values)
end

function to_logical(mesh::ToroidalAnnulusGrid{T}, point::NTuple{3,Real}) where {T}
    _tah_geometry(mesh)
    x, y, z = ntuple(d -> _require_finite_real("annular Cartesian coordinate", point[d], T), 3)
    major = hypot(x, y)
    r = hypot(major - mesh.R0, z)
    isfinite(major) && isfinite(r) && major > 0 && mesh.r_inner <= r <= mesh.r_outer ||
        throw(ArgumentError("Cartesian point is outside the admitted toroidal annulus"))
    xi = _tah_round(T, _exact_rational(r) - _exact_rational(mesh.r_inner), "annular inverse radius")
    u = _periodic_coordinate(atan(z, major - mesh.R0) / mesh.tau, zero(T), one(T), one(T))
    v = _periodic_coordinate(atan(y, x) / mesh.tau, zero(T), one(T), one(T))
    return (xi, u, v)
end
