# gyrokinetics.jl -- guiding-centre and gyro-averaged particle dynamics.

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross3(a, b) =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline _norm3(a) = sqrt(_dot3(a, a))
@inline _scale3(s, a) = (s * a[1], s * a[2], s * a[3])
@inline _add3(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])

function _normalize3(a)
    T = eltype(a)
    n = _norm3(a)
    n > 0 || throw(ArgumentError("cannot normalize a zero vector"))
    return _scale3(one(T) / n, a)
end

"""
    exb_drift(E, B) -> NTuple{3}

The `E×B` drift `v_E = (E × B)/B^2`. `E` and `B` are 3-component vectors.
"""
function exb_drift(E, B)
    B2 = _dot3(B, B)
    B2 > 0 || throw(ArgumentError("magnetic field magnitude must be positive"))
    return _scale3(one(eltype(B)) / B2, _cross3(E, B))
end

"""
    gradb_drift(vperp, q, m, B, gradB) -> NTuple{3}

Grad-B drift `v_∇B = (m v⊥^2)/(2 q B^3) (B × ∇B)`, where `gradB = ∇|B|`.
"""
function gradb_drift(vperp, q, m, B, gradB)
    Bmag = _norm3(B)
    Bmag > 0 || throw(ArgumentError("magnetic field magnitude must be positive"))
    q != 0 || throw(ArgumentError("charge q must be nonzero"))
    return _scale3(m * vperp^2 / (2 * q * Bmag^3), _cross3(B, gradB))
end

"""
    curvature_drift(vpar, q, m, B, κ) -> NTuple{3}

Curvature drift `v_κ = (m v∥^2)/(q B^2) (B × κ)`.
"""
function curvature_drift(vpar, q, m, B, κ)
    B2 = _dot3(B, B)
    B2 > 0 || throw(ArgumentError("magnetic field magnitude must be positive"))
    q != 0 || throw(ArgumentError("charge q must be nonzero"))
    return _scale3(m * vpar^2 / (q * B2), _cross3(B, κ))
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
mutable struct GuidingCentre{T}
    X::NTuple{3,T}
    vpar::T
    μ::T
    q::T
    m::T
end

GuidingCentre(X::NTuple{3,T}, vpar, μ, q, m) where {T} =
    GuidingCentre{T}(X, T(vpar), T(μ), T(q), T(m))

"""
    push_guiding_centre!(gc; dt, E, B, gradB, κ, gradpar_B) -> gc

Advance the guiding centre one local-field step. The position advances by
perpendicular drifts plus parallel streaming `v∥ b`; `v∥` advances by
`m dv∥/dt = q E∥ - μ ∂∥B`. `μ` is conserved.
"""
function push_guiding_centre!(gc::GuidingCentre{T}; dt, E, B, gradB, κ, gradpar_B) where {T}
    gc.m != 0 || throw(ArgumentError("mass m must be nonzero"))
    Bmag = _norm3(B)
    Bmag > 0 || throw(ArgumentError("magnetic field magnitude must be positive"))
    b = _scale3(one(T) / Bmag, B)
    vperp = sqrt(2 * gc.μ * Bmag / gc.m)
    vd = drift_velocity(; vpar = gc.vpar, vperp, q = gc.q, m = gc.m, E, B, gradB, κ)
    gc.X = _add3(gc.X, _scale3(T(dt), _add3(vd, _scale3(gc.vpar, b))))
    Epar = _dot3(E, b)
    gc.vpar += T(dt) * (gc.q * Epar - gc.μ * gradpar_B) / gc.m
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
function gyroaverage(f, X::NTuple{3,T}, ρ, B; n::Integer = 16) where {T}
    n >= 3 || throw(ArgumentError("need n >= 3 ring points"))
    b = _normalize3(B)
    e1, e2 = _perp_basis(b)
    ρT = T(ρ)
    s = zero(T)
    @inbounds for k = 1:n
        φ = 2 * T(π) * (k - 1) / n
        pt = _add3(X, _add3(_scale3(ρT * cos(φ), e1), _scale3(ρT * sin(φ), e2)))
        s += f(pt)
    end
    return s / n
end
