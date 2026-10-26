# Constrained entries lie on the two end planes of each conducting axis along
# which the component is nodal (tangential E, normal B); only those O(surface)
# planes are visited.
@inline _yee_constrained_axis(g::YeeGrid, electric::Bool, c::Int, d::Int) =
    !g.periodic[d] && _yee_node(electric, c, d)

# Concrete ranges also for callers whose array types are inferred abstractly.
@inline function _yee_wall_planes(a::AbstractArray{T,D}, g::YeeGrid{D}, d::Int) where {T,D}
    return ntuple(2) do side
        index = side == 1 ? 1 : g.n[d] + 1
        CartesianIndices(ntuple(e -> e == d ? (index:index) : UnitRange{Int}(1, size(a, e)), D))
    end
end

function _yee_walls_zero(a::AbstractArray{T,D}, g::YeeGrid{D}, electric::Bool, c::Int) where {T,D}
    for d = 1:D
        _yee_constrained_axis(g, electric, c, d) || continue
        for plane in _yee_wall_planes(a, g, d), I in plane
            iszero(a[I]) || return false
        end
    end
    return true
end

"""
    apply_yee_pec!(fields::YeeFields, g)

Set tangential E and normal B to zero at every conducting face, edge and corner.
This explicit initialization operation changes boundary values. The vacuum
stepper instead rejects inconsistent boundary input before changing any state.
The homogeneous normal B condition describes a cavity with no trapped flux.
"""
function apply_yee_pec!(fields::YeeFields, g::YeeGrid{D}) where {D}
    _yee_check_shape(fields.E, g, true)
    _yee_check_shape(fields.B, g, false)
    _yee_check_independent_arrays((fields.E..., fields.B...))
    for (components, electric) in ((fields.E, true), (fields.B, false)), c = 1:3
        a = components[c]
        for d = 1:D
            _yee_constrained_axis(g, electric, c, d) || continue
            for plane in _yee_wall_planes(a, g, d), I in plane
                a[I] = zero(eltype(a))
            end
        end
    end
    return fields
end

function _yee_validate_electric_pec(components, g::YeeGrid)
    for c = 1:3
        _yee_walls_zero(components[c], g, true, c) ||
            throw(ArgumentError("conducting walls require zero tangential electric data"))
    end
    return nothing
end

function _yee_validate_pec(fields::YeeFields, g::YeeGrid)
    for (components, electric) in ((fields.E, true), (fields.B, false)), c = 1:3
        _yee_walls_zero(components[c], g, electric, c) ||
            throw(ArgumentError("conducting walls require zero tangential E and normal B"))
    end
    return nothing
end
