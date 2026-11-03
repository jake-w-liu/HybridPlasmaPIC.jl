# validation.jl - runtime guards that must remain active in optimized runs.

# Structural zeros use stored IEEE bits even when the caller enables FTZ.
# The masked test runs in a @noinline integer helper: LLVM folds a float-domain
# masked reinterpret into `fcmp oeq x, 0.0`, which honors caller FTZ/DAZ modes
# and would flush a stored subnormal to zero.
@noinline _stored_bits_zero(bits::U, mask::U) where {U<:Unsigned} = iszero(bits & mask)
_particle_stored_zero(x) = iszero(x)
_particle_stored_zero(x::Float16) = _stored_bits_zero(reinterpret(UInt16, x), UInt16(0x7fff))
_particle_stored_zero(x::Float32) = _stored_bits_zero(reinterpret(UInt32, x), UInt32(0x7fffffff))
_particle_stored_zero(x::Float64) =
    _stored_bits_zero(reinterpret(UInt64, x), UInt64(0x7fffffffffffffff))

# Callers validate finiteness separately; these predicates classify stored sign/zero.
_stored_positive(x) = x > zero(x)
_stored_nonnegative(x) = x >= zero(x)
_stored_positive(x::Union{Float16,Float32,Float64}) =
    !isnan(x) & !signbit(x) & !_particle_stored_zero(x)
_stored_nonnegative(x::Union{Float16,Float32,Float64}) =
    !isnan(x) & (!signbit(x) | _particle_stored_zero(x))

function _require_all_finite(name::AbstractString, a, context::AbstractString)
    all(isfinite, a) || error("$(name) went non-finite ($(context))")
    return nothing
end

function _require_finite_real(
    name::AbstractString,
    value::Real,
    ::Type{T},
    index::Union{Nothing,Integer} = nothing,
) where {T<:AbstractFloat}
    v = T(value)
    if !isfinite(v)
        label = index === nothing ? name : "$name[$index]"
        throw(ArgumentError("$label must be finite"))
    end
    return v
end

function _require_finite_nonnegative_real(
    name::AbstractString,
    value::Real,
    ::Type{T},
) where {T<:AbstractFloat}
    v = _require_finite_real(name, value, T)
    v >= zero(T) || throw(ArgumentError("$(name) must be finite and non-negative"))
    return v
end

function _require_finite_positive_real(
    name::AbstractString,
    value::Real,
    ::Type{T},
) where {T<:AbstractFloat}
    v = _require_finite_real(name, value, T)
    v > zero(T) || throw(ArgumentError("$(name) must be finite and positive"))
    return v
end

function _require_positive_intlike(name::AbstractString, value)
    value isa Integer ||
        value isa Real ||
        throw(ArgumentError("$(name) must be a positive integer"))
    if value isa Real
        isfinite(value) || throw(ArgumentError("$(name) must be a finite positive integer"))
        isinteger(value) || throw(ArgumentError("$(name) must be an integer-valued positive count"))
    end
    value > 0 || throw(ArgumentError("$(name) must be a positive integer"))
    value <= typemax(Int) || throw(ArgumentError("$(name) must fit in Int"))
    return Int(value)
end

function _check_spatial_dimension(D)
    D isa Int || throw(ArgumentError("D must be an Int spatial dimension"))
    1 <= D <= 3 || throw(ArgumentError("D must be 1, 2, or 3 for 1D/2D/3D spatial models"))
    return nothing
end

function _validated_nonnegative_dt(
    ::Type{T},
    dt::Real;
    name::AbstractString,
) where {T<:AbstractFloat}
    dtT = T(dt)
    isfinite(dtT) || throw(ArgumentError("$name requires finite dt (got $dt)"))
    _stored_nonnegative(dtT) || throw(ArgumentError("$name requires nonnegative dt (got $dt)"))
    return dtT
end

# Particle pushers and steppers form dt/2 for leapfrog priming, midpoints and
# staggered kicks. A nonzero |dt| < 2*floatmin(T) is subnormal or has a
# subnormal half step, which silently erases or truncates those terms, so it is
# rejected before any state changes. dt = 0 keeps each caller's own semantics.
function _require_resolved_dt(name::AbstractString, dt::T) where {T<:AbstractFloat}
    _particle_stored_zero(dt) ||
        abs(dt) >= 2 * floatmin(T) ||
        throw(
            ArgumentError(
                "$name requires dt = 0 or |dt| >= 2*floatmin($T) (got $dt): subnormal timesteps and half steps are not resolved",
            ),
        )
    return dt
end

function _yl_source_scalar(name, value, ::Type{T}) where {T}
    if value isa Union{Float32,Float64}
        isfinite(value) || throw(ArgumentError("$name must be finite"))
        return value
    end
    return _require_finite_real(name, value, T)
end
