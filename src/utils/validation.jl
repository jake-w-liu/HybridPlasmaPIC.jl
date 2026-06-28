# validation.jl - runtime guards that must remain active in optimized runs.

function _require_all_finite(name::AbstractString, a, context::AbstractString)
    all(isfinite, a) || error("$(name) went non-finite ($(context))")
    return nothing
end

function _require_finite_real(name::AbstractString, value::Real, ::Type{T}) where {T<:AbstractFloat}
    v = T(value)
    isfinite(v) || throw(ArgumentError("$(name) must be finite"))
    return v
end

function _require_finite_nonnegative_real(name::AbstractString, value::Real, ::Type{T}) where {T<:AbstractFloat}
    v = _require_finite_real(name, value, T)
    v >= zero(T) || throw(ArgumentError("$(name) must be finite and non-negative"))
    return v
end

function _require_finite_positive_real(name::AbstractString, value::Real, ::Type{T}) where {T<:AbstractFloat}
    v = _require_finite_real(name, value, T)
    v > zero(T) || throw(ArgumentError("$(name) must be finite and positive"))
    return v
end

function _require_positive_intlike(name::AbstractString, value)
    value isa Integer || value isa Real || throw(ArgumentError("$(name) must be a positive integer"))
    if value isa Real
        isfinite(value) || throw(ArgumentError("$(name) must be a finite positive integer"))
        isinteger(value) || throw(ArgumentError("$(name) must be an integer-valued positive count"))
    end
    value > 0 || throw(ArgumentError("$(name) must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$(name) must fit in Int"))
    return Int(value)
end

function _validated_nonnegative_dt(::Type{T}, dt::Real; name::AbstractString) where {T<:AbstractFloat}
    dtT = T(dt)
    isfinite(dtT) || throw(ArgumentError("$name requires finite dt (got $dt)"))
    dtT >= zero(T) || throw(ArgumentError("$name requires nonnegative dt (got $dt)"))
    return dtT
end
