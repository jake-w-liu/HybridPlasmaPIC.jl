"""
    YeeFields(g::YeeGrid)

Zero electric edge fields `E` and magnetic face fields `B`. Each is a tuple of
three ordinary CPU arrays. Use [`yee_coordinates`](@ref) to initialize them and
[`apply_yee_pec!`](@ref) to impose homogeneous conducting-wall data explicitly.
"""
struct YeeFields{D,T}
    E::NTuple{3,Array{T,D}}
    B::NTuple{3,Array{T,D}}
end

function YeeFields(g::YeeGrid{D,T}) where {D,T}
    return YeeFields{D,T}(
        ntuple(c -> zeros(T, _yee_shape(g, true, c)), 3),
        ntuple(c -> zeros(T, _yee_shape(g, false, c)), 3),
    )
end

# Branch-free reduction: x*0 is a signed zero for finite x and NaN otherwise,
# and NaN survives any (SIMD-reassociated) order of the additions.
function _yee_all_finite(a::AbstractArray{T}) where {T<:AbstractFloat}
    flag = zero(T)
    @inbounds @simd for i in eachindex(a)
        flag += a[i] * zero(T)
    end
    return isfinite(flag)
end

function _yee_validate_impulse(impulse, g::YeeGrid{1,T}) where {T}
    impulse isa Tuple && length(impulse) == 3 ||
        throw(ArgumentError("impulse must be a tuple of three current arrays"))
    all(a -> a isa Vector{T}, impulse) || throw(
        ArgumentError(
            "Yee current arrays must be ordinary CPU vectors matching the grid precision",
        ),
    )
    _yee_validate_electric_data(impulse, g)
    return nothing
end

function _yee_validate_electric_data(components, g::YeeGrid)
    _yee_check_shape(components, g, true)
    for a in components
        _yee_all_finite(a) || throw(ArgumentError("Yee electric data must be finite"))
    end
    _yee_validate_electric_pec(components, g)
    return nothing
end

function _yee_validate_fields(fields::YeeFields, g::YeeGrid)
    _yee_check_shape(fields.E, g, true)
    _yee_check_shape(fields.B, g, false)
    for a in (fields.E..., fields.B...)
        _yee_all_finite(a) || throw(ArgumentError("Yee fields must be finite"))
    end
    return nothing
end

function _yee_check_independent_arrays(arrays)
    for j in eachindex(arrays), i = 1:j-1
        Base.mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("Yee physical and staging arrays must not alias"))
    end
    return nothing
end
