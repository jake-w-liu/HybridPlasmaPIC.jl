# HybridPIC.jl — HybridModel + hybrid moment computation (from hybrid.jl)

"""
    HybridModel(closure; η=0.0, ηH=0.0, nfloor=1e-6)

Hybrid model parameters: electron `closure`, resistivity `η` and hyperresistivity
`ηH` (both off by default; Ohm's law adds `+η J − ηH ∇²J`), and a density floor
used in the 1/n divisions of Ohm's law.
"""
struct HybridModel{C<:ElectronClosure}
    closure::C
    η::Float64
    ηH::Float64
    nfloor::Float64
end
function HybridModel(closure::ElectronClosure; η = 0.0, ηH = 0.0, nfloor = 1e-6)
    ηT = _require_finite_nonnegative_real("η", η, Float64)
    ηHT = _require_finite_nonnegative_real("ηH", ηH, Float64)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, Float64)
    return HybridModel(closure, ηT, ηHT, nfloorT)
end


"""
    compute_moments!(f, ps, g, shape, nfloor; work=...)

Deposit number density into `f.n` and ion bulk velocity into `f.ui`
(= momentum density / max(n, nfloor)).
"""
function compute_moments!(
    f::HybridFields{D,T},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor;
    work::AbstractVector{T} = Vector{T}(undef, nparticles(ps)),
) where {D,T}
    density!(f.n, ps, g, shape)
    momentum!(f.ui, ps, g, shape; work)          # holds (n u) on entry
    nf = T(nfloor)
    @inbounds for I in eachindex(f.n)
        inv = one(T) / max(f.n[I], nf)
        f.ui[1][I] *= inv
        f.ui[2][I] *= inv
        f.ui[3][I] *= inv
    end
    return f
end

# ---------------------------------------------------------------- Ohm's law


"""
    compute_moments_multi!(f, species, g, shape, nfloor)

Quasineutral moments from several ion species (each a `ParticleSet` with its own
charge `q`): `f.n` ← e·n_e = Σ_s q_s n_s, and the charge-weighted ion bulk
velocity `f.ui` ← (Σ_s q_s n_s u_s)/(Σ_s q_s n_s). Reduces to the single-species
result for one proton population (q=1).
"""
function compute_moments_multi!(
    f::HybridFields{D,T},
    species::AbstractVector{<:ParticleSet{D,T}},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor,
    ;
    ntmp::Array{T,D} = similar(f.n),
    mtmp::NTuple{3,<:Array{T,D}} = ntuple(_ -> similar(f.n), 3),
    works = nothing,
) where {D,T}
    size(ntmp) == g.n ||
        throw(DimensionMismatch("ntmp size $(size(ntmp)) does not match grid size $(g.n)"))
    for c = 1:3
        size(mtmp[c]) == g.n || throw(
            DimensionMismatch("mtmp[$c] size $(size(mtmp[c])) does not match grid size $(g.n)"),
        )
    end
    if works !== nothing && length(works) != length(species)
        throw(
            DimensionMismatch(
                "works length $(length(works)) must equal species count $(length(species))",
            ),
        )
    end
    fill!(f.n, zero(T))
    for c = 1:3
        fill!(f.ui[c], zero(T))
    end
    for (is, s) in enumerate(species)
        work = if works === nothing
            Vector{T}(undef, nparticles(s))
        else
            w = works[is]
            eltype(w) === T ||
                throw(ArgumentError("works[$is] eltype $(eltype(w)) must match $T"))
            length(w) == nparticles(s) || throw(
                DimensionMismatch(
                    "works[$is] length $(length(w)) must equal particle count $(nparticles(s))",
                ),
            )
            w
        end
        density!(ntmp, s, g, shape)
        momentum!(mtmp, s, g, shape; work)           # (n u)_s
        @. f.n += s.q * ntmp                          # Σ q_s n_s
        for c = 1:3
            @. f.ui[c] += s.q * mtmp[c]               # Σ q_s (n u)_s
        end
    end
    nf = T(nfloor)
    @inbounds for I in eachindex(f.n)
        inv = one(T) / max(f.n[I], nf)
        f.ui[1][I] *= inv
        f.ui[2][I] *= inv
        f.ui[3][I] *= inv
    end
    return f
end

# ---------------------------------------------------------------- electron energy
