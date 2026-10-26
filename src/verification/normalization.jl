# normalization.jl
#
# SI <-> normalized unit conversion for the hybrid plasma solver (checklist §7).
#
# The solver works in Ω_ci-normalized units (proton q/m = 1, n0 = 1, μ0 = 1).
# `PlasmaUnits` carries the three independent SI reference scales (density n0,
# magnetic field B0, ion mass mi) plus the physical constants needed to build
# every derived scale. All derived normalization scales follow from these.

"""
    PlasmaUnits{T}

Reference SI scales for converting between SI and Ω_ci-normalized hybrid-plasma
units. Constructed from three independent references:

- `n0` — reference number density [m^-3]
- `B0` — reference magnetic-field strength [T]
- `mi` — ion mass [kg]

Stored constants (SI): elementary charge `e`, vacuum permeability `mu0`, and the
speed of light `c`. Derived scales are obtained through the accessors
[`alfven_speed`](@ref), [`gyrofrequency`](@ref) and [`inertial_length`](@ref).

Construct via the keyword constructor:

    PlasmaUnits(; n0, B0, mi)

All six stored references must be finite and positive. The full positional
constructor accepts six values of the same floating-point type;
`PlasmaUnits{T}(n0, B0, mi, e, mu0, c)` converts them to `T` before validation.
Passing existing `BigFloat` references to the typed constructor retains their
stored precision. Arithmetic results use the caller's working precision.
"""
struct PlasmaUnits{T<:AbstractFloat}
    n0::T   # reference number density [m^-3]
    B0::T   # reference magnetic field [T]
    mi::T   # ion mass [kg]
    e::T    # elementary charge [C]
    mu0::T  # vacuum permeability [H/m]
    c::T    # speed of light [m/s]

    function PlasmaUnits{T}(n0, B0, mi, e, mu0, c) where {T<:AbstractFloat}
        values = map(x -> convert(T, x), (n0, B0, mi, e, mu0, c))
        for (name, value) in zip((:n0, :B0, :mi, :e, :mu0, :c), values)
            isfinite(value) && value > 0 ||
                throw(ArgumentError("$name must be finite and positive"))
        end
        return new{T}(values...)
    end
end

PlasmaUnits(n0::T, B0::T, mi::T, e::T, mu0::T, c::T) where {T<:AbstractFloat} =
    PlasmaUnits{T}(n0, B0, mi, e, mu0, c)

"""
    PlasmaUnits(; n0, B0, mi)

Build a [`PlasmaUnits`](@ref) from reference density `n0` [m^-3], magnetic field
`B0` [T] and ion mass `mi` [kg]. The physical constants `e`, `mu0`, `c` are set
to their CODATA SI values. The three references are promoted to a common
floating-point type.
"""
function PlasmaUnits(; n0, B0, mi)
    e = 1.602176634e-19
    mu0 = 1.25663706212e-6
    c = 2.99792458e8
    T = float(promote_type(typeof(n0), typeof(B0), typeof(mi), typeof(e), typeof(mu0), typeof(c)))
    n0T = _require_finite_positive_real("n0", n0, T)
    B0T = _require_finite_positive_real("B0", B0, T)
    miT = _require_finite_positive_real("mi", mi, T)
    return PlasmaUnits{T}(n0T, B0T, miT, T(e), T(mu0), T(c))
end

"""
    alfven_speed(u::PlasmaUnits) -> T

Reference Alfvén speed `v_A = B0 / sqrt(mu0 * n0 * mi)` [m/s].
"""
alfven_speed(u::PlasmaUnits) = _scale(:velocity, u)

"""
    gyrofrequency(u::PlasmaUnits) -> T

Reference ion gyrofrequency `Ω_ci = e * B0 / mi` [rad/s].
"""
function gyrofrequency(u::PlasmaUnits{T}) where {T}
    value, ordinary = _gyrofrequency_scale(u)
    ordinary && return value
    squared, power = _normalization_squared_scale((0, 2, -2, 2, 0), u)
    return _normalization_materialize(T, squared, power, false)
end

"""
    inertial_length(u::PlasmaUnits) -> T

Reference ion inertial length `d_i = v_A / Ω_ci` [m].
"""
inertial_length(u::PlasmaUnits) = _scale(:length, u)

# Keep ordinary arithmetic and its rounding order. A zero, subnormal or
# nonfinite intermediate instead requires the complete physical expression.
@inline _ordinary_unit_scale(x::T) where {T<:AbstractFloat} = isfinite(x) && x >= floatmin(T)
# MPFR has no subnormals. Directed overflow/underflow can return finite
# endpoints. Interior exponents certify the same strict endpoint predicate
# without constructing two BigFloat bounds for every ordinary arithmetic result.
@inline function _normalization_interior_magnitude(x::BigFloat)
    isfinite(x) && !iszero(x) || return false
    e = exponent(x)
    minimum = ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ())
    maximum = ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ())
    minimum - 1 <= e < maximum || return false
    if e == minimum - 1
        bound = floatmin(BigFloat)
        (signbit(x) ? x < -bound : x > bound) || return false
    end
    if e == maximum - 1
        bound = floatmax(BigFloat)
        (signbit(x) ? x > -bound : x < bound) || return false
    end
    return true
end

@inline _ordinary_unit_scale(x::BigFloat) = x > 0 && _normalization_interior_magnitude(x)

@inline function _alfven_scale(u::PlasmaUnits)
    density = u.mu0 * u.n0
    product = density * u.mi
    root = sqrt(product)
    value = u.B0 / root
    ordinary =
        _ordinary_unit_scale(density) &&
        _ordinary_unit_scale(product) &&
        _ordinary_unit_scale(root) &&
        _ordinary_unit_scale(value)
    return value, ordinary
end

@inline function _gyrofrequency_scale(u::PlasmaUnits)
    product = u.e * u.B0
    value = product / u.mi
    return value, _ordinary_unit_scale(product) && _ordinary_unit_scale(value)
end

@inline function _inertial_scale(u::PlasmaUnits)
    velocity, safe_velocity = _alfven_scale(u)
    frequency, safe_frequency = _gyrofrequency_scale(u)
    value = velocity / frequency
    return value, safe_velocity && safe_frequency && _ordinary_unit_scale(value)
end

@noinline function _unknown_normalization_kind(kind)
    throw(
        ArgumentError(
            "unknown normalization kind :$(kind); valid kinds are " *
            ":length, :time, :velocity, :magnetic, :electric, :density, :current, :pressure",
        ),
    )
end

@inline function _ordinary_normalization_scale(kind::Symbol, u::PlasmaUnits)
    if kind === :length
        return _inertial_scale(u)
    elseif kind === :time
        frequency, ordinary = _gyrofrequency_scale(u)
        value = one(frequency) / frequency
        return value, ordinary && _ordinary_unit_scale(value)
    elseif kind === :velocity
        return _alfven_scale(u)
    elseif kind === :magnetic
        return u.B0, _ordinary_unit_scale(u.B0)
    elseif kind === :electric
        velocity, ordinary = _alfven_scale(u)
        value = velocity * u.B0
        return value, ordinary && _ordinary_unit_scale(value)
    elseif kind === :density
        return u.n0, _ordinary_unit_scale(u.n0)
    elseif kind === :current
        length, ordinary = _inertial_scale(u)
        product = u.mu0 * length
        value = u.B0 / product
        return value, ordinary && _ordinary_unit_scale(product) && _ordinary_unit_scale(value)
    elseif kind === :pressure
        product = u.B0^2
        value = product / u.mu0
        return value, _ordinary_unit_scale(product) && _ordinary_unit_scale(value)
    end
    return _unknown_normalization_kind(kind)
end

# Powers of (n0, B0, mi, e, mu0) in the SQUARE of the physical scale.
# Squaring expresses every scale as a rational monomial, including v_A and d_i.
function _normalization_powers(kind::Symbol)
    kind === :length && return (-1, 0, 1, -2, -1)
    kind === :time && return (0, -2, 2, -2, 0)
    kind === :velocity && return (-1, 2, -1, 0, -1)
    kind === :magnetic && return (0, 2, 0, 0, 0)
    kind === :electric && return (-1, 4, -1, 0, -1)
    kind === :density && return (2, 0, 0, 0, 0)
    kind === :current && return (1, 2, -1, 2, -1)
    kind === :pressure && return (0, 4, 0, 0, -2)
    return _unknown_normalization_kind(kind)
end

function _normalization_squared_scale(powers::NTuple{5,Int}, u::PlasmaUnits)
    squared = BigInt(1) // BigInt(1)
    power = Int128(0)
    for (value, count) in zip((u.n0, u.B0, u.mi, u.e, u.mu0), powers)
        iszero(count) && continue
        n, e, d = Base.decompose(value)
        squared *= (BigInt(n) // BigInt(d))^count
        power += Int128(e) * count
    end
    return squared, power
end

function _normalization_materialize(::Type{T}, squared, power, negative) where {T<:AbstractFloat}
    mode = T === BigFloat ? rounding(BigFloat) : RoundNearest
    return _sqrt_from_exact(T, squared, power, mode, negative)
end

function _scale(kind::Symbol, u::PlasmaUnits{T}) where {T}
    value, ordinary = _ordinary_normalization_scale(kind, u)
    ordinary && return value
    squared, power = _normalization_squared_scale(_normalization_powers(kind), u)
    return _normalization_materialize(T, squared, power, false)
end

function _normalization_exact_value(value::Real, ::Type{U}, squared, power, inverse) where {U}
    T = promote_type(typeof(value), U)
    # A physical scale is strictly positive even if it lies outside T's range.
    # Multiplying/dividing zero or infinity by that scale retains its sign.
    (!isfinite(value) || iszero(value)) && return T(value)
    n, e, d = Base.decompose(value)
    factor = (BigInt(n) // BigInt(d))^2
    combined = inverse ? factor / squared : factor * squared
    combined_power = 2 * Int128(e) + (inverse ? -power : power)
    return _normalization_materialize(T, combined, combined_power, signbit(value))
end

# Irrational values follow Julia's existing promotion to a floating operand;
# unlike stored integers, rationals and floats, they have no exact decomposition.
function _normalization_exact_value(
    value::AbstractIrrational,
    ::Type{T},
    squared,
    power,
    inverse,
) where {T}
    promoted = convert(promote_type(typeof(value), T), value)
    return _normalization_exact_value(promoted, T, squared, power, inverse)
end

function _normalization_exact_value(value::Complex, T, squared, power, inverse)
    return complex(
        _normalization_exact_value(real(value), T, squared, power, inverse),
        _normalization_exact_value(imag(value), T, squared, power, inverse),
    )
end

_normalization_exact_value(::Missing, T, squared, power, inverse) = missing

@inline _normalization_finite_nonzero(value::Real) = isfinite(value) && !iszero(value)
@inline function _normalization_finite_nonzero(value::T) where {T<:Union{Float16,Float32,Float64}}
    # Device floating comparisons can treat stored subnormals as zero. Inspect
    # IEEE magnitude bits so a lost finite input still reaches exact conversion.
    U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
    magnitude = reinterpret(U, value) & (typemax(U) >> 1)
    return zero(U) < magnitude < reinterpret(U, T(Inf))
end

@inline _needs_complete_normalization(value, result) = false
@inline function _needs_complete_normalization(value::Real, result::AbstractFloat)
    return _normalization_finite_nonzero(value) && !_ordinary_unit_scale(abs(result))
end
@inline function _needs_complete_normalization(value::Real, result::BigFloat)
    return _normalization_finite_nonzero(value) && !_normalization_interior_magnitude(result)
end
@inline function _needs_complete_normalization(value::Complex, result::Complex)
    return _needs_complete_normalization(real(value), real(result)) ||
           _needs_complete_normalization(imag(value), imag(result))
end

function _normalization_convert_value(value, kind, u::PlasmaUnits{T}, scale, inverse) where {T}
    result = inverse ? value / scale : value * scale
    _needs_complete_normalization(value, result) || return result
    squared, power = _normalization_squared_scale(_normalization_powers(kind), u)
    return _normalization_exact_value(value, T, squared, power, inverse)
end

function _convert_normalization(value, kind::Symbol, u::PlasmaUnits{T}, inverse) where {T}
    scale, ordinary = _ordinary_normalization_scale(kind, u)
    if ordinary
        return broadcast(x -> _normalization_convert_value(x, kind, u, scale, inverse), value)
    end
    squared, power = _normalization_squared_scale(_normalization_powers(kind), u)
    return broadcast(x -> _normalization_exact_value(x, T, squared, power, inverse), value)
end

# Device extensions delegate here. Ordinary arithmetic and its inspection run
# on the device; only exceptional expressions need arbitrary-precision host work.
_normalization_device_needs_complete(value, result) =
    mapreduce(_needs_complete_normalization, |, value, result; init = false)

function _convert_normalization_device(value, kind, u, inverse)
    # Julia's zero-dimensional broadcast returns a scalar. Stage its single
    # input explicitly rather than invoking device scalar getindex/setindex.
    ndims(value) == 0 && return _convert_normalization(Array(value), kind, u, inverse)
    scale, ordinary = _ordinary_normalization_scale(kind, u)
    if ordinary
        result = inverse ? value ./ scale : value .* scale
        needs_exact = _normalization_device_needs_complete(value, result)
        needs_exact || return result
    end
    host = _convert_normalization(Array(value), kind, u, inverse)
    result = similar(value, eltype(host), size(host))
    copyto!(result, host)
    return result
end

"""
    to_SI(value, kind::Symbol, u::PlasmaUnits)

Convert a normalized `value` of physical quantity `kind` to SI units by
multiplying by the reference scale.

`kind` ∈ (`:length`, `:time`, `:velocity`, `:magnetic`, `:electric`,
`:density`, `:current`, `:pressure`). `value` may be a scalar or an array.
Ordinary arithmetic retains its evaluation order. When an intermediate scale or
final operation loses range, conversion evaluates the complete physical
expression before rounding; a representable result can therefore be recovered
even when the reference scale itself is outside the output type's range.
`BigFloat` arithmetic uses the current working precision and rounding mode.
Zero signs and nonfinite values are preserved under the positive physical scale;
complex values are converted componentwise.

CUDA/Metal array extensions inspect ordinary results on the device, with a
synchronization to decide whether exact host arithmetic is needed. Exceptional
conversions stage through host memory and return the original backend and shape.
Zero-dimensional arrays follow Julia broadcast behavior and return a scalar.
"""
to_SI(value, kind::Symbol, u::PlasmaUnits) = _convert_normalization(value, kind, u, false)

"""
    to_normalized(value, kind::Symbol, u::PlasmaUnits)

Convert an SI `value` of physical quantity `kind` to Ω_ci-normalized units by
dividing by the reference scale. Uses the same complete-expression, precision
and device-staging rules as [`to_SI`](@ref). Floating-point rounding can prevent
an exact round trip.
"""
to_normalized(value, kind::Symbol, u::PlasmaUnits) = _convert_normalization(value, kind, u, true)

"""
    coulomb_normalization(u::PlasmaUnits) -> (; epsilon0, weight_to_number)

Explicit collision parameters for Ω_ci-normalized positions, velocities, time,
charge `q/e`, mass `m/mi`, and weights from [`set_density_weight!`](@ref) with
density in units `u.n0`. Pass the result as keywords to cell-local
[`collide_coulomb!`](@ref).
The collision velocity floor is specified separately in normalized velocity
units; its SI counterpart is `to_SI(u_floor, :velocity, u)`.

With physical reference scales `L0 = d_i` and `T0 = 1/Ω_ci`, the returned
permittivity is `epsilon0_SI*mi*L0^3/(e^2*T0^2)`, with
`epsilon0_SI = 1/(mu0*c^2)`, and `weight_to_number = n0*L0^3`. For D<3,
the normalized bin has unit extent in each omitted dimension: its physical
volume is `prod(grid.dx)*L0^3`. A stored weight thus counts particles in that
reference transverse area/length, not an unspecified physical cross section.

Both parameters are evaluated from the complete stored-reference expressions,
without materializing an intermediate reference volume. BigFloat outputs use
the caller's working precision and rounding mode. The results must be finite
and positive to configure the collision operator.

These parameters require the same reference plasma as the field and particle
normalization. Neither a default EM permittivity nor a Coulomb logarithm alone
determines the physical binary collision rate.
"""
function coulomb_normalization(u::PlasmaUnits{T}) where {T}
    # Squared physical scales, with powers of (n0, B0, mi, e, mu0).
    count_squared, count_power = _normalization_squared_scale((-1, 0, 3, -6, -3), u)
    epsilon_squared, epsilon_power = _normalization_squared_scale((-3, 4, 1, -6, -5), u)
    # The squared permittivity also contains c^-4. Keep its binary scale
    # separate, including when reference lengths/volumes are unrepresentable.
    n, e, d = Base.decompose(u.c)
    epsilon_squared *= (BigInt(n) // BigInt(d))^-4
    epsilon_power -= 4 * Int128(e)
    epsilon0 = _normalization_materialize(T, epsilon_squared, epsilon_power, false)
    weight_to_number = _normalization_materialize(T, count_squared, count_power, false)
    _require_finite_positive_real("normalized collision permittivity", epsilon0, T)
    _require_finite_positive_real("normalized collision particle count", weight_to_number, T)
    return (; epsilon0, weight_to_number)
end
