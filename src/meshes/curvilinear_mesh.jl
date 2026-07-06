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
`r` uses `Nr` cell centres on `(0, a]` (offset `dr/2` to avoid the `r=0` axis
singularity); `θ, φ` are periodic with `Nθ, Nφ` nodes on `[0, 2π)`.
"""
struct ToroidalGrid{T}
    R0::T
    r::Vector{T}
    θ::Vector{T}
    φ::Vector{T}
    dr::T
    dθ::T
    dφ::T
end

function ToroidalGrid(R0::Real, a::Real, Nr::Integer, Nθ::Integer, Nφ::Integer; T::Type = Float64)
    (Nr >= 3 && Nθ >= 3 && Nφ >= 3) ||
        throw(ArgumentError("need Nr >= 3 and Nθ,Nφ >= 3 for the finite-difference stencils"))
    a > 0 || throw(ArgumentError("minor radius a must be > 0"))
    R0 > a || throw(ArgumentError("R0 must exceed a (no self-intersecting torus)"))
    dr = T(a) / Nr
    dθ = 2 * T(π) / Nθ
    dφ = 2 * T(π) / Nφ
    r = T[dr * (i - T(0.5)) for i = 1:Nr]
    θ = T[dθ * (j - 1) for j = 1:Nθ]
    φ = T[dφ * (k - 1) for k = 1:Nφ]
    return ToroidalGrid{T}(T(R0), r, θ, φ, dr, dθ, dφ)
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
    R = g.R0 + T(r) * cos(T(θ))
    return (R * cos(T(φ)), R * sin(T(φ)), T(r) * sin(T(θ)))
end

@inline _dperiodic(f, lo, hi, h) = (f[hi] - f[lo]) / (2h)

"""
    metric_gradient(g, f) -> (gr, gθ, gφ)

Physical components of `∇f` for a scalar field `f[i,j,k]` on a toroidal grid:
`(∇f)_r = ∂_r f`, `(∇f)_θ = (1/r)∂_θ f`, `(∇f)_φ = (1/R)∂_φ f`.
Central differences are used in the interior, second-order one-sided stencils at
the `r` boundaries, and periodic central differences in `θ,φ`.
"""
function metric_gradient(g::ToroidalGrid{T}, f::AbstractArray{T,3}) where {T}
    Nr, Nθ, Nφ = gridsize(g)
    size(f) == (Nr, Nθ, Nφ) || throw(DimensionMismatch("f must be $(gridsize(g))"))
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
        size(A) == (Nr, Nθ, Nφ) || throw(DimensionMismatch("$nm must be $(gridsize(g))"))
    end
    out = similar(Ar)
    Fr = similar(Ar)
    Fθ = similar(Aθ)
    Fφ = similar(Aφ)
    @inbounds for k = 1:Nφ, j = 1:Nθ, i = 1:Nr
        r = g.r[i]
        R = g.R0 + r * cos(g.θ[j])
        Fr[i, j, k] = r * R * Ar[i, j, k]
        Fθ[i, j, k] = R * Aθ[i, j, k]
        Fφ[i, j, k] = r * Aφ[i, j, k]
    end
    @inbounds for k = 1:Nφ
        kp = k == Nφ ? 1 : k + 1
        km = k == 1 ? Nφ : k - 1
        for j = 1:Nθ
            jp = j == Nθ ? 1 : j + 1
            jm = j == 1 ? Nθ : j - 1
            for i = 1:Nr
                r = g.r[i]
                R = g.R0 + r * cos(g.θ[j])
                J = r * R
                if i == 1
                    dFr = (-3Fr[1, j, k] + 4Fr[2, j, k] - Fr[3, j, k]) / (2 * g.dr)
                elseif i == Nr
                    dFr = (3Fr[Nr, j, k] - 4Fr[Nr-1, j, k] + Fr[Nr-2, j, k]) / (2 * g.dr)
                else
                    dFr = (Fr[i+1, j, k] - Fr[i-1, j, k]) / (2 * g.dr)
                end
                dFθ = (Fθ[i, jp, k] - Fθ[i, jm, k]) / (2 * g.dθ)
                dFφ = (Fφ[i, j, kp] - Fφ[i, j, km]) / (2 * g.dφ)
                out[i, j, k] = (dFr + dFθ + dFφ) / J
            end
        end
    end
    return out
end
