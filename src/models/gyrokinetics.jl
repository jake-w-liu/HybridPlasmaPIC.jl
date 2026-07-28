# gyrokinetics.jl -- guiding-centre and gyro-averaged particle dynamics.

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross3(a, b) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline _norm3(a) = sqrt(_dot3(a, a))
@inline _scale3(s, a) = (s * a[1], s * a[2], s * a[3])
@inline _add3(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])

function _require_finite_vector3(name::AbstractString, a)
    length(a) == 3 || throw(DimensionMismatch("$name must have exactly three components"))
    @inbounds for i = 1:3
        a[i] isa Real || throw(ArgumentError("$name[$i] must be real"))
        isfinite(a[i]) || throw(ArgumentError("$name[$i] must be finite"))
    end
    return nothing
end

function _convert_finite_vector3(name::AbstractString, a, ::Type{T}) where {T<:AbstractFloat}
    _require_finite_vector3(name, a)
    return ntuple(i -> _require_finite_real("$name[$i]", a[i], T), 3)
end

function _require_finite_scalar(name::AbstractString, x)
    x isa Real || throw(ArgumentError("$name must be real"))
    isfinite(x) || throw(ArgumentError("$name must be finite"))
    return x
end

function _require_positive_scalar(name::AbstractString, x)
    _require_finite_scalar(name, x)
    x > zero(x) || throw(ArgumentError("$name must be positive"))
    return x
end

function _require_nonnegative_scalar(name::AbstractString, x)
    _require_finite_scalar(name, x)
    x >= zero(x) || throw(ArgumentError("$name must be nonnegative"))
    return x
end

function _normalize3(a)
    _require_finite_vector3("vector", a)
    n = _norm3(a)
    isfinite(n) && n > zero(n) ||
        throw(ArgumentError("cannot normalize a zero or non-finite vector"))
    return _scale3(inv(n), a)
end

"""
    exb_drift(E, B) -> NTuple{3}

The `E×B` drift `v_E = (E × B)/B^2`. `E` and `B` are 3-component vectors.
"""
function exb_drift(E, B)
    _require_finite_vector3("E", E)
    _require_finite_vector3("B", B)
    B2 = _dot3(B, B)
    isfinite(B2) && B2 > zero(B2) ||
        throw(ArgumentError("magnetic field magnitude must be finite and positive"))
    drift = _scale3(inv(B2), _cross3(E, B))
    _require_finite_vector3("E×B drift", drift)
    return drift
end

"""
    gradb_drift(vperp, q, m, B, gradB) -> NTuple{3}

Grad-B drift `v_∇B = (m v⊥^2)/(2 q B^3) (B × ∇B)`, where `gradB = ∇|B|`.
"""
function gradb_drift(vperp, q, m, B, gradB)
    _require_nonnegative_scalar("vperp", vperp)
    _require_finite_scalar("q", q)
    q != zero(q) || throw(ArgumentError("charge q must be nonzero"))
    _require_positive_scalar("m", m)
    _require_finite_vector3("B", B)
    _require_finite_vector3("gradB", gradB)
    Bmag = _norm3(B)
    isfinite(Bmag) && Bmag > zero(Bmag) ||
        throw(ArgumentError("magnetic field magnitude must be finite and positive"))
    drift = _scale3(m * vperp^2 / (2 * q * Bmag^3), _cross3(B, gradB))
    _require_finite_vector3("grad-B drift", drift)
    return drift
end

"""
    curvature_drift(vpar, q, m, B, κ) -> NTuple{3}

Curvature drift `v_κ = (m v∥^2)/(q B^2) (B × κ)`.
"""
function curvature_drift(vpar, q, m, B, κ)
    _require_finite_scalar("vpar", vpar)
    _require_finite_scalar("q", q)
    q != zero(q) || throw(ArgumentError("charge q must be nonzero"))
    _require_positive_scalar("m", m)
    _require_finite_vector3("B", B)
    _require_finite_vector3("κ", κ)
    B2 = _dot3(B, B)
    isfinite(B2) && B2 > zero(B2) ||
        throw(ArgumentError("magnetic field magnitude must be finite and positive"))
    drift = _scale3(m * vpar^2 / (q * B2), _cross3(B, κ))
    _require_finite_vector3("curvature drift", drift)
    return drift
end

"""
    drift_velocity(; vpar, vperp, q, m, E, B, gradB, κ) -> NTuple{3}

Total perpendicular guiding-centre drift `v_E + v_∇B + v_κ`.
"""
drift_velocity(; vpar, vperp, q, m, E, B, gradB, κ) = _add3(
    _add3(exb_drift(E, B), gradb_drift(vperp, q, m, B, gradB)),
    curvature_drift(vpar, q, m, B, κ),
)

"""
    GuidingCentre(X, vpar, μ, q, m)

Guiding-centre state with position `X`, parallel velocity `vpar`, magnetic
moment `μ = m v⊥^2/(2B)`, charge `q`, and mass `m`.
"""
mutable struct GuidingCentre{T<:AbstractFloat}
    X::NTuple{3,T}
    vpar::T
    μ::T
    q::T
    m::T

    function GuidingCentre{T}(X::NTuple{3,T}, vpar::T, μ::T, q::T, m::T) where {T<:AbstractFloat}
        _require_finite_vector3("X", X)
        _require_finite_scalar("vpar", vpar)
        _require_nonnegative_scalar("μ", μ)
        _require_finite_scalar("q", q)
        q != zero(T) || throw(ArgumentError("charge q must be nonzero"))
        _require_positive_scalar("m", m)
        return new{T}(X, vpar, μ, q, m)
    end
end

function GuidingCentre(X::NTuple{3,S}, vpar, μ, q, m) where {S<:Real}
    T = S <: AbstractFloat ? S : typeof(float(zero(S)))
    XT = ntuple(i -> _require_finite_real("X[$i]", X[i], T), 3)
    return GuidingCentre{T}(XT, T(vpar), T(μ), T(q), T(m))
end

function _validate_guiding_centre(gc::GuidingCentre)
    _require_finite_vector3("X", gc.X)
    _require_finite_scalar("vpar", gc.vpar)
    _require_nonnegative_scalar("μ", gc.μ)
    _require_finite_scalar("q", gc.q)
    gc.q != zero(gc.q) || throw(ArgumentError("charge q must be nonzero"))
    _require_positive_scalar("m", gc.m)
    return nothing
end

"""
    push_guiding_centre!(gc; dt, E, B, gradB, κ, gradpar_B) -> gc

Advance the guiding centre one local-field step. The position advances by
perpendicular drifts plus parallel streaming `v∥ b`; `v∥` advances by
`m dv∥/dt = q E∥ - μ ∂∥B`. `μ` is conserved.
"""
function push_guiding_centre!(gc::GuidingCentre{T}; dt, E, B, gradB, κ, gradpar_B) where {T}
    _validate_guiding_centre(gc)
    dtT = _require_finite_real("dt", dt, T)
    ET = _convert_finite_vector3("E", E, T)
    BT = _convert_finite_vector3("B", B, T)
    gradBT = _convert_finite_vector3("gradB", gradB, T)
    κT = _convert_finite_vector3("κ", κ, T)
    gradparT = _require_finite_real("gradpar_B", gradpar_B, T)
    Bmag = _norm3(BT)
    isfinite(Bmag) && Bmag > zero(T) ||
        throw(ArgumentError("magnetic field magnitude must be finite and positive"))
    b = _scale3(inv(Bmag), BT)
    vperp = sqrt(2 * gc.μ * Bmag / gc.m)
    vd = drift_velocity(;
        vpar = gc.vpar,
        vperp,
        q = gc.q,
        m = gc.m,
        E = ET,
        B = BT,
        gradB = gradBT,
        κ = κT,
    )
    new_X = _add3(gc.X, _scale3(dtT, _add3(vd, _scale3(gc.vpar, b))))
    Epar = _dot3(ET, b)
    new_vpar = gc.vpar + dtT * (gc.q * Epar - gc.μ * gradparT) / gc.m
    _require_finite_vector3("updated X", new_X)
    _require_finite_scalar("updated vpar", new_vpar)

    # Commit only after all direct, derived, and updated values are valid.
    gc.X = new_X
    gc.vpar = new_vpar
    return gc
end

function _perp_basis(b::NTuple{3,T}) where {T}
    ref = abs(b[1]) < T(0.9) ? (one(T), zero(T), zero(T)) : (zero(T), one(T), zero(T))
    e1 = _normalize3(_sub3(ref, _scale3(_dot3(ref, b), b)))
    e2 = _cross3(b, e1)
    return e1, e2
end

"""
    gyroaverage(f, X, ρ, B; n=16) -> value

Gyro-average of scalar callable `f(x)` over a ring of radius `ρ` centred at `X`
in the plane perpendicular to `B`, sampled at `n` points.
"""
function gyroaverage(f, X::NTuple{3,T}, ρ, B; n::Integer = 16) where {T<:AbstractFloat}
    nring = _require_positive_intlike("n", n)
    nring >= 3 || throw(ArgumentError("need n >= 3 ring points"))
    XT = _convert_finite_vector3("X", X, T)
    BT = _convert_finite_vector3("B", B, T)
    ρT = _require_finite_nonnegative_real("ρ", ρ, T)
    b = _normalize3(BT)
    e1, e2 = _perp_basis(b)
    pt = _add3(XT, _scale3(ρT, e1))
    s = f(pt)
    @inbounds for k = 2:nring
        φ = 2 * T(π) * (k - 1) / nring
        pt = _add3(XT, _add3(_scale3(ρT * cos(φ), e1), _scale3(ρT * sin(φ), e2)))
        s += f(pt)
    end
    return s / nring
end
