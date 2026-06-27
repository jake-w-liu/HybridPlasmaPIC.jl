# PressureTensor.jl — ∥/⊥ temperatures + pressure-strain (from diagnostics.jl)

@inline _grid_axes(g::FourierGrid{D}) where {D} = ntuple(d -> Base.OneTo(g.n[d]), Val(D))

function _require_axes(name::Symbol, a::AbstractArray, ref_axes)
    axes(a) == ref_axes ||
        throw(DimensionMismatch("$(name) axes $(axes(a)) do not match expected axes $(ref_axes)"))
    return nothing
end

function _require_axes(name::Symbol, c::Int, a::AbstractArray, ref_axes)
    axes(a) == ref_axes || throw(
        DimensionMismatch("$(name)[$c] axes $(axes(a)) do not match expected axes $(ref_axes)"),
    )
    return nothing
end

function _require_grid_array(name::Symbol, a::AbstractArray, g::FourierGrid{D}) where {D}
    ndims(a) == D ||
        throw(DimensionMismatch("$(name) has ndims=$(ndims(a)); expected $D for grid size $(g.n)"))
    size(a) == g.n ||
        throw(DimensionMismatch("$(name) size $(size(a)) does not match grid size $(g.n)"))
    _require_axes(name, a, _grid_axes(g))
    return nothing
end

function _require_grid_array(name::Symbol, c::Int, a::AbstractArray, g::FourierGrid{D}) where {D}
    ndims(a) == D || throw(
        DimensionMismatch("$(name)[$c] has ndims=$(ndims(a)); expected $D for grid size $(g.n)"),
    )
    size(a) == g.n ||
        throw(DimensionMismatch("$(name)[$c] size $(size(a)) does not match grid size $(g.n)"))
    _require_axes(name, c, a, _grid_axes(g))
    return nothing
end

function _require_grid_tuple(name::Symbol, arrays::Tuple, g::FourierGrid)
    for c = 1:length(arrays)
        _require_grid_array(name, c, arrays[c], g)
    end
    return nothing
end

function _require_same_axes(name::Symbol, arrays::Tuple)
    isempty(arrays) && return nothing
    ref_axes = axes(arrays[1])
    for c = 2:length(arrays)
        _require_axes(name, c, arrays[c], ref_axes)
    end
    return nothing
end

"""
    temperatures_par_perp(P, n, B; nfloor=1e-6)

Per-cell parallel and perpendicular temperatures from the pressure tensor `P`
(6 components xx,yy,zz,xy,xz,yz) and magnetic field `B=(Bx,By,Bz)`:
`T∥ = b̂·P·b̂ / n`, `T⊥ = (trP − b̂·P·b̂)/(2n)`.
"""
function temperatures_par_perp(
    P::NTuple{6,<:AbstractArray{T}},
    n::AbstractArray{T},
    B::NTuple{3,<:AbstractArray{T}};
    nfloor = 1e-6,
) where {T}
    _require_same_axes(:temperature_input, (n, P..., B...))
    Tpar = similar(n)
    Tperp = similar(n)
    nf = T(nfloor)
    isfinite(nf) && nf > zero(T) ||
        throw(ArgumentError("nfloor must be finite and positive"))
    Pxx, Pyy, Pzz, Pxy, Pxz, Pyz = P
    Bx, By, Bz = B
    @inbounds for I in eachindex(n)
        bx = Bx[I]
        by = By[I]
        bz = Bz[I]
        b2 = bx * bx + by * by + bz * bz
        ninv = one(T) / max(n[I], nf)
        tr = Pxx[I] + Pyy[I] + Pzz[I]
        if b2 > 0
            ppar =
                (
                    bx * bx * Pxx[I] +
                    by * by * Pyy[I] +
                    bz * bz * Pzz[I] +
                    2 * (bx * by * Pxy[I] + bx * bz * Pxz[I] + by * bz * Pyz[I])
                ) / b2
            Tpar[I] = ppar * ninv
            Tperp[I] = (tr - ppar) / 2 * ninv
        else
            Tpar[I] = tr / 3 * ninv
            Tperp[I] = tr / 3 * ninv
        end
    end
    return Tpar, Tperp
end

# ---------------------------------------------------------------- distributions


"""
    pressure_strain(P, u, g)

Pressure–strain interaction `−Σ_ij P_ij ∂_j u_i` per cell (the rate of internal-
energy exchange), using the spectral gradient of the bulk velocity `u`.
"""
function pressure_strain(
    P::NTuple{6,<:AbstractArray{T,D}},
    u::NTuple{3,<:AbstractArray{T,D}},
    g::FourierGrid{D,T},
) where {T,D}
    _require_grid_tuple(:P, P, g)
    _require_grid_tuple(:u, u, g)
    Pidx = ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))  # symmetric component map
    Pmat(i, j) = begin
        for (k, (a, b)) in enumerate(Pidx)
            ((a == i && b == j) || (a == j && b == i)) && return P[k]
        end
    end
    out = zeros(T, g.n)
    duij = similar(out)
    for j = 1:D, i = 1:3
        deriv!(duij, u[i], g, j)             # ∂_j u_i
        Pij = Pmat(i, j)
        @inbounds for I in eachindex(out)
            out[I] -= Pij[I] * duij[I]
        end
    end
    return out
end

# ---------------------------------------------------------------- shock front
