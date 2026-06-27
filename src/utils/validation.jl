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
