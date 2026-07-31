# curvilinear_mesh.jl -- metric-aware toroidal mesh operators.
#
# The spectral Cartesian field engine uses uniform-grid wavenumber multipliers.
# Toroidal geometry needs an explicit metric/Jacobian layer, so these operators
# live with the other mesh backends and provide finite-difference derivatives in
# orthogonal toroidal coordinates.
#
# Toroidal map (r,θ,φ) -> Cartesian:
#     x = (R0 + r cosθ) cosφ,  y = (R0 + r cosθ) sinφ,  z = r sinθ.
# It is orthogonal with scale factors h_r=1, h_θ=r, h_φ=R0+r cosθ and
# Jacobian J = r(R0+r cosθ).

"""
    ToroidalGrid(R0, a, Nr, Nθ, Nφ; T=Float64)

A toroidal `(r,θ,φ)` grid on the torus of major radius `R0` and minor radius `a`.
`r` uses `Nr` cell centres on `(0, a)` (offset `dr/2` to avoid the `r=0` axis
singularity); `θ, φ` are periodic with `Nθ, Nφ` nodes on `[0, 2π)`.
"""
struct ToroidalGrid{T<:AbstractFloat}
    R0::T
    r::Vector{T}
    θ::Vector{T}
    φ::Vector{T}
    dr::T
    dθ::T
    dφ::T
end

function ToroidalGrid(R0::Real, a::Real, Nr::Integer, Nθ::Integer, Nφ::Integer; T::Type = Float64)
    isconcretetype(T) && T <: AbstractFloat ||
        throw(ArgumentError("ToroidalGrid T must be a concrete AbstractFloat type"))
    counts = (Nr, Nθ, Nφ)
    all(n -> 3 <= n <= typemax(Int), counts) ||
        throw(ArgumentError("need Nr >= 3 and Nθ,Nφ >= 3 for the finite-difference stencils"))
    Nri, Nθi, Nφi = Int.(counts)
    aT = _require_finite_positive_real("minor radius a", a, T)
    R0T = _require_finite_positive_real("major radius R0", R0, T)
    R0T > aT || throw(ArgumentError("R0 must exceed a (no self-intersecting torus)"))
    dr = aT / Nri
    dθ = T(2) * T(π) / Nθi
    dφ = T(2) * T(π) / Nφi
    all(x -> isfinite(x) && x > zero(T), (dr, dθ, dφ)) || throw(
        ArgumentError("ToroidalGrid spacings must remain finite and positive at precision $T"),
    )
    r = T[dr * (i - T(0.5)) for i = 1:Nri]
    θ = T[dθ * (j - 1) for j = 1:Nθi]
    φ = T[dφ * (k - 1) for k = 1:Nφi]
    return ToroidalGrid{T}(R0T, r, θ, φ, dr, dθ, dφ)
end

gridsize(g::ToroidalGrid) = (length(g.r), length(g.θ), length(g.φ))

"""
    scale_factors(g, i, j) -> (h_r, h_θ, h_φ)

Metric scale factors at grid node `(i,j,·)` (φ-independent): `h_r=1`,
`h_θ=r_i`, `h_φ = R0 + r_i cosθ_j`.
"""
@inline function scale_factors(g::ToroidalGrid{T}, i::Integer, j::Integer) where {T}
    r = g.r[i]
    R = g.R0 + r * cos(g.θ[j])
    return (one(T), r, R)
end

"Jacobian `J = r*R` at node `(i,j)`."
@inline function jacobian(g::ToroidalGrid{T}, i::Integer, j::Integer) where {T}
    h = scale_factors(g, i, j)
    return h[1] * h[2] * h[3]
end

"Cartesian position of `(r,θ,φ)`."
function to_cartesian(g::ToroidalGrid{T}, r::Real, θ::Real, φ::Real) where {T}
    rT = _require_finite_real("r", r, T)
    θT = _require_finite_real("θ", θ, T)
    φT = _require_finite_real("φ", φ, T)
    R = g.R0 + rT * cos(θT)
    return (R * cos(φT), R * sin(φT), rT * sin(θT))
end

@inline _dperiodic(f, lo, hi, h) = (f[hi] - f[lo]) / (2h)

function _require_toroidal_field(g::ToroidalGrid, name::Symbol, A::AbstractArray{<:Any,3})
    expected_size = gridsize(g)
    size(A) == expected_size || throw(DimensionMismatch("$name must be $expected_size"))
    expected_axes = ntuple(d -> Base.OneTo(expected_size[d]), 3)
    axes(A) == expected_axes ||
        throw(ArgumentError("$name must use one-based axes $expected_axes, got $(axes(A))"))
    return nothing
end

"""
    metric_gradient(g, f) -> (gr, gθ, gφ)

Physical components of `∇f` for a scalar field `f[i,j,k]` on a toroidal grid:
`(∇f)_r = ∂_r f`, `(∇f)_θ = (1/r)∂_θ f`, `(∇f)_φ = (1/R)∂_φ f`.
Central differences are used in the interior, second-order one-sided stencils at
the `r` boundaries, and periodic central differences in `θ,φ`.
"""
function metric_gradient(g::ToroidalGrid{T}, f::AbstractArray{T,3}) where {T}
    Nr, Nθ, Nφ = gridsize(g)
    _require_toroidal_field(g, :f, f)
    gr = similar(f)
    gθ = similar(f)
    gφ = similar(f)
    @inbounds for k = 1:Nφ
        kp = k == Nφ ? 1 : k + 1
        km = k == 1 ? Nφ : k - 1
        for j = 1:Nθ
            jp = j == Nθ ? 1 : j + 1
            jm = j == 1 ? Nθ : j - 1
            for i = 1:Nr
                r = g.r[i]
                R = g.R0 + r * cos(g.θ[j])
                if i == 1
                    dfr = (-3f[1, j, k] + 4f[2, j, k] - f[3, j, k]) / (2 * g.dr)
                elseif i == Nr
                    dfr = (3f[Nr, j, k] - 4f[Nr-1, j, k] + f[Nr-2, j, k]) / (2 * g.dr)
                else
                    dfr = (f[i+1, j, k] - f[i-1, j, k]) / (2 * g.dr)
                end
                gr[i, j, k] = dfr
                gθ[i, j, k] = _dperiodic(@view(f[i, :, k]), jm, jp, g.dθ) / r
                gφ[i, j, k] = _dperiodic(@view(f[i, j, :]), km, kp, g.dφ) / R
            end
        end
    end
    return gr, gθ, gφ
end

@inline function _radial_metric_flux(g, Ar, i, j, k, cosθ)
    r = g.r[i]
    return r * (g.R0 + r * cosθ) * Ar[i, j, k]
end

@inline function _radial_flux_derivative(g, Ar, i, j, k, cosθ, Nr)
    dr = g.dr
    if i == 1
        # At the coordinate axis, Fr = r*R*Ar is exactly zero for every
        # bounded physical Ar.  Including that known value with the first
        # three cell-centred fluxes gives a cubic-exact derivative at r=dr/2.
        return (
            _radial_metric_flux(g, Ar, 1, j, k, cosθ) / 2 +
            2 * _radial_metric_flux(g, Ar, 2, j, k, cosθ) / 3 -
            _radial_metric_flux(g, Ar, 3, j, k, cosθ) / 10
        ) / dr
    elseif i == Nr
        # J is bounded away from zero at the outer edge, so the usual
        # second-order one-sided derivative retains second-order accuracy.
        return (
            3 * _radial_metric_flux(g, Ar, Nr, j, k, cosθ) -
            4 * _radial_metric_flux(g, Ar, Nr - 1, j, k, cosθ) +
            _radial_metric_flux(g, Ar, Nr - 2, j, k, cosθ)
        ) / (2dr)
    elseif Nr == 3
        # Cubic interpolation through Fr(0)=0 and all three cell centres.
        return (
            -3 * _radial_metric_flux(g, Ar, 1, j, k, cosθ) / 2 +
            2 * _radial_metric_flux(g, Ar, 2, j, k, cosθ) / 3 +
            3 * _radial_metric_flux(g, Ar, 3, j, k, cosθ) / 10
        ) / dr
    elseif i == 2
        # Cubic-exact forward-biased derivative.
        return (
            -2 * _radial_metric_flux(g, Ar, 1, j, k, cosθ) -
            3 * _radial_metric_flux(g, Ar, 2, j, k, cosθ) +
            6 * _radial_metric_flux(g, Ar, 3, j, k, cosθ) -
            _radial_metric_flux(g, Ar, 4, j, k, cosθ)
        ) / (6dr)
    elseif i == Nr - 1
        # Cubic-exact backward-biased derivative.
        return (
            _radial_metric_flux(g, Ar, Nr - 3, j, k, cosθ) -
            6 * _radial_metric_flux(g, Ar, Nr - 2, j, k, cosθ) +
            3 * _radial_metric_flux(g, Ar, Nr - 1, j, k, cosθ) +
            2 * _radial_metric_flux(g, Ar, Nr, j, k, cosθ)
        ) / (6dr)
    else
        # The fourth-order centred stencil keeps the flux-derivative error
        # below O(dr^3), which is needed near J=rR=O(dr) for the divergence
        # itself to remain at least second-order accurate.
        return (
            _radial_metric_flux(g, Ar, i - 2, j, k, cosθ) -
            8 * _radial_metric_flux(g, Ar, i - 1, j, k, cosθ) +
            8 * _radial_metric_flux(g, Ar, i + 1, j, k, cosθ) -
            _radial_metric_flux(g, Ar, i + 2, j, k, cosθ)
        ) / (12dr)
    end
end

"""
    metric_divergence(g, Ar, Aθ, Aφ) -> divA

Divergence of a vector field with physical components `(Ar,Aθ,Aφ)` on the
toroidal grid, using
`∇⋅A = (1/J)[∂_r(J A_r) + ∂_θ(R A_θ) + ∂_φ(r A_φ)]` with `J = rR`.
"""
function metric_divergence(
    g::ToroidalGrid{T},
    Ar::AbstractArray{T,3},
    Aθ::AbstractArray{T,3},
    Aφ::AbstractArray{T,3},
) where {T}
    Nr, Nθ, Nφ = gridsize(g)
    for (nm, A) in ((:Ar, Ar), (:Aθ, Aθ), (:Aφ, Aφ))
        _require_toroidal_field(g, nm, A)
    end
    out = similar(Ar)
    @inbounds for k = 1:Nφ
        kp = k == Nφ ? 1 : k + 1
        km = k == 1 ? Nφ : k - 1
        for j = 1:Nθ
            jp = j == Nθ ? 1 : j + 1
            jm = j == 1 ? Nθ : j - 1
            cosθ = cos(g.θ[j])
            cosθp = cos(g.θ[jp])
            cosθm = cos(g.θ[jm])
            for i = 1:Nr
                r = g.r[i]
                R = g.R0 + r * cosθ
                J = r * R
                dFr = _radial_flux_derivative(g, Ar, i, j, k, cosθ, Nr)
                Rjp = g.R0 + r * cosθp
                Rjm = g.R0 + r * cosθm
                dFθ = (Rjp * Aθ[i, jp, k] - Rjm * Aθ[i, jm, k]) / (2 * g.dθ)
                dFφ = (r * Aφ[i, j, kp] - r * Aφ[i, j, km]) / (2 * g.dφ)
                out[i, j, k] = (dFr + dFθ + dFφ) / J
            end
        end
    end
    return out
end
