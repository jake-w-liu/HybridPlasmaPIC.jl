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
"""
struct PlasmaUnits{T<:AbstractFloat}
    n0::T   # reference number density [m^-3]
    B0::T   # reference magnetic field [T]
    mi::T   # ion mass [kg]
    e::T    # elementary charge [C]
    mu0::T  # vacuum permeability [H/m]
    c::T    # speed of light [m/s]
end

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
    return PlasmaUnits{T}(T(n0), T(B0), T(mi), T(e), T(mu0), T(c))
end

"""
    alfven_speed(u::PlasmaUnits) -> T

Reference Alfvén speed `v_A = B0 / sqrt(mu0 * n0 * mi)` [m/s].
"""
alfven_speed(u::PlasmaUnits) = u.B0 / sqrt(u.mu0 * u.n0 * u.mi)

"""
    gyrofrequency(u::PlasmaUnits) -> T

Reference ion gyrofrequency `Ω_ci = e * B0 / mi` [rad/s].
"""
gyrofrequency(u::PlasmaUnits) = u.e * u.B0 / u.mi

"""
    inertial_length(u::PlasmaUnits) -> T

Reference ion inertial length `d_i = v_A / Ω_ci` [m].
"""
inertial_length(u::PlasmaUnits) = alfven_speed(u) / gyrofrequency(u)

# Reference scale (in SI units) for each normalization `kind`.
# A normalized value of 1.0 corresponds to exactly this many SI units.
@inline function _scale(kind::Symbol, u::PlasmaUnits)
    vA = alfven_speed(u)
    di = inertial_length(u)
    if kind === :length
        return di                       # [m]            × d_i
    elseif kind === :time
        return one(vA) / gyrofrequency(u)   # [s]        × Ω_ci^-1
    elseif kind === :velocity
        return vA                       # [m/s]          × v_A
    elseif kind === :magnetic
        return u.B0                     # [T]            × B0
    elseif kind === :electric
        return vA * u.B0                # [V/m]          × v_A B0
    elseif kind === :density
        return u.n0                     # [m^-3]         × n0
    elseif kind === :current
        return u.B0 / (u.mu0 * di)      # [A/m^2]        × B0/(mu0 d_i)
    elseif kind === :pressure
        return u.B0^2 / u.mu0           # [Pa]           × B0^2/mu0
    else
        throw(
            ArgumentError(
                "unknown normalization kind :$(kind); valid kinds are " *
                ":length, :time, :velocity, :magnetic, :electric, :density, :current, :pressure",
            ),
        )
    end
end

"""
    to_SI(value, kind::Symbol, u::PlasmaUnits)

Convert a normalized `value` of physical quantity `kind` to SI units by
multiplying by the reference scale.

`kind` ∈ (`:length`, `:time`, `:velocity`, `:magnetic`, `:electric`,
`:density`, `:current`, `:pressure`). `value` may be a scalar or an array.
"""
to_SI(value, kind::Symbol, u::PlasmaUnits) = value .* _scale(kind, u)

"""
    to_normalized(value, kind::Symbol, u::PlasmaUnits)

Convert an SI `value` of physical quantity `kind` to Ω_ci-normalized units by
dividing by the reference scale. Inverse of [`to_SI`](@ref).
"""
to_normalized(value, kind::Symbol, u::PlasmaUnits) = value ./ _scale(kind, u)
