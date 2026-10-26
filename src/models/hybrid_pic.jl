# HybridPIC.jl — HybridModel + hybrid moment computation (from hybrid.jl)

"""
    HybridModel(closure; η=0.0, ηH=0.0, nfloor=1e-6, de2=0.0)

Hybrid model parameters: electron `closure`, resistivity `η` and hyperresistivity
`ηH` (both off by default; Ohm's law adds `+η J − ηH ∇²J`), a density floor used in
the 1/n divisions of Ohm's law, and the **electron inertial length squared** `de2 = d_e²`
(finite electron mass). With `de2 > 0` the generalized Ohm's law keeps the leading
electron-inertia term `(d_e²/n) ∂_t J`, which — for the transverse field on a uniform
background — is the spectral multiplier `E ← E/(1 + d_e² k²)`: it caps the whistler
frequency at `ω_pe` scales and enables collisionless reconnection at the `d_e` scale.
`de2 = 0` recovers the massless (ideal) hybrid exactly.

[`EvolvingElectrons`](@ref) requires `η = ηH = de2 = 0`; its independent
pressure is initialized in `HybridStepper.fields.pe` before calling `init!`.
[`TransportElectrons`](@ref) accepts constant `η ≥ 0` with `ηH = de2 = 0`; a
nonzero `η` must stay finite and nonzero after conversion to grid precision.
"""
struct HybridModel{C<:ElectronClosure}
    closure::C
    η::Float64
    ηH::Float64
    nfloor::Float64
    de2::Float64        # electron inertial length squared d_e²; 0 ⇒ massless electrons
end
function HybridModel(closure::ElectronClosure; η = 0.0, ηH = 0.0, nfloor = 1e-6, de2 = 0.0)
    ηT = _require_finite_nonnegative_real("η", η, Float64)
    ηHT = _require_finite_nonnegative_real("ηH", ηH, Float64)
    nfloorT = _require_finite_positive_real("nfloor", nfloor, Float64)
    de2T = _require_finite_nonnegative_real("de2", de2, Float64)
    model = HybridModel(closure, ηT, ηHT, nfloorT, de2T)
    _validate_evolving_model(model, Float64)
    return model
end


# Also used at stepping and restore boundaries: positional model construction
# and deserialization must not bypass the evolving closure's physical contract.
function _validate_evolving_model(model::HybridModel, ::Type{T}) where {T}
    if _evolves_pressure(model.closure)
        γ = _require_finite_real("γ", model.closure.γ, T)
        γ > one(T) ||
            throw(ArgumentError("evolving electron γ must be greater than one at grid precision"))
        if _transports_electrons(model.closure)
            _validate_transport_resistivity(model, T)
        else
            (iszero(model.η) && iszero(model.ηH) && iszero(model.de2)) ||
                throw(ArgumentError("EvolvingElectrons requires η = ηH = de2 = 0"))
        end
        _require_finite_positive_real("nfloor", model.nfloor, T)
    end
    return nothing
end

# Constant Ohmic transport only. A nonzero coefficient must survive conversion:
# overflow or disappearance to zero would silently change the model.
function _validate_transport_resistivity(model::HybridModel, ::Type{T}) where {T}
    (iszero(model.ηH) && iszero(model.de2)) ||
        throw(ArgumentError("TransportElectrons supports constant η only; requires ηH = de2 = 0"))
    η = _require_finite_nonnegative_real("η", model.η, T)
    (iszero(model.η) || !iszero(η)) ||
        throw(ArgumentError("TransportElectrons η underflows to zero at grid precision"))
    return η
end

# Transport requires the actual deposited density, never its floor: a floor
# would hide a missing population and break the paired force/heat balance.
function _validate_transport_density(n, nfloor::T, context) where {T}
    for value in n
        isfinite(value) && value >= nfloor || throw(
            ArgumentError(
                "TransportElectrons $context density must be finite and at least nfloor = $nfloor",
            ),
        )
    end
    return nothing
end


"""
    compute_moments!(f, ps, g, shape, nfloor; work=...)

Deposit number density into `f.n` and ion bulk velocity into `f.ui`
(= momentum density / max(n, nfloor)). Outputs and optional one-based scratch
must not overlap each other or particle storage; all are checked before writing.
"""
function compute_moments!(
    f::HybridFields{D,T},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor;
    work::Union{Nothing,AbstractVector{T}} = nothing,
) where {D,T}
    nf = _require_finite_positive_real("nfloor", nfloor, T)
    _deposit_density_momentum_moments!(f.n, f.ui, ps, g, shape, work)
    _normalize_ion_moments!(f.n, f.ui, nf)
    return f
end

function _normalize_ion_moments!(n::Array{T}, u, nfloor::T) where {T}
    @inbounds for I in eachindex(n)
        inv = one(T) / max(n[I], nfloor)
        u[1][I] *= inv
        u[2][I] *= inv
        u[3][I] *= inv
    end
    return n
end

function _zero_ion_moments!(n::Array{T}, u) where {T}
    fill!(n, zero(T))
    foreach(a -> fill!(a, zero(T)), u)
    return n
end

# Accumulate raw charge moments; normalize only after every species contributes.
function _accumulate_charge_moments!(n, u, ps, g, shape, ntmp, mtmp, work)
    _deposit_density_momentum_moments!(ntmp, mtmp, ps, g, shape, work)
    @. n += ps.q * ntmp
    for c = 1:3
        @. u[c] += ps.q * mtmp[c]
    end
    return n
end

# ---------------------------------------------------------------- Ohm's law


"""
    compute_moments_multi!(f, species, g, shape, nfloor)

Quasineutral moments from several ion species (each a `ParticleSet` with its own
charge `q`): `f.n` ← e·n_e = Σ_s q_s n_s, and the charge-weighted ion bulk
velocity `f.ui` ← (Σ_s q_s n_s u_s)/(Σ_s q_s n_s). Reduces to the single-species
result for one proton population (q=1). Grid outputs and scratch must be
independent and must not overlap any population. Each optional `works` entry
is a one-based particle vector; entries may share scratch reused sequentially,
but must not overlap particle or grid storage. Every population is checked
before any output or scratch is changed.
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
    nf = _require_finite_positive_real("nfloor", nfloor, T)
    _validate_moment_storage((f.n, f.ui..., ntmp, mtmp...), species, g, works)
    _zero_ion_moments!(f.n, f.ui)
    for (is, s) in enumerate(species)
        work = works === nothing ? nothing : works[is]
        _accumulate_charge_moments!(f.n, f.ui, s, g, shape, ntmp, mtmp, work)
    end
    _normalize_ion_moments!(f.n, f.ui, nf)
    return f
end

# ---------------------------------------------------------------- electron energy
