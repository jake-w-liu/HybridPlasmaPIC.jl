# FieldState.jl — HybridFields container (from hybrid.jl)

"""
    HybridFields{D,T}(n::NTuple{D,Int})

Grid-resident hybrid fields and scratch: density `n`, ion bulk velocity `ui`,
magnetic `B`, electric `E`, current `J`, electron pressure `pe`, its gradient
`gradp`, reciprocal density `ninv`, and a density-floor activation counter.
"""
mutable struct HybridFields{D,T}
    n::Array{T,D}
    ui::NTuple{3,Array{T,D}}
    B::NTuple{3,Array{T,D}}
    E::NTuple{3,Array{T,D}}
    J::NTuple{3,Array{T,D}}
    pe::Array{T,D}
    gradp::NTuple{D,Array{T,D}}
    ninv::Array{T,D}
    lapJ::NTuple{3,Array{T,D}}        # ∇²J workspace (hyperresistivity)
    pforce::NTuple{3,Array{T,D}}      # ∇·P_e output (anisotropic/CGL only; 0-length for scalar)
    floor_count::Base.RefValue{Int}
end

function _hybrid_fields(::Val{D}, ::Type{T}, nc::NTuple{D,Int}, anisotropic::Bool) where {D,T}
    z() = zeros(T, nc)
    pf() = anisotropic ? zeros(T, nc) : zeros(T, ntuple(_ -> 0, D))
    HybridFields{D,T}(
        z(),
        ntuple(_ -> z(), 3),
        ntuple(_ -> z(), 3),
        ntuple(_ -> z(), 3),
        ntuple(_ -> z(), 3),
        z(),
        ntuple(_ -> z(), D),
        z(),
        ntuple(_ -> z(), 3),
        ntuple(_ -> pf(), 3),
        Ref(0),
    )
end

# `anisotropic=true` (a CGL closure) allocates the `pforce` ∇·P_e buffer full-size; scalar
# closures get 0-length pforce (type-stable, never indexed on the scalar Ohm path) so the
# common case carries no dead weight.
function HybridFields{D,T}(nc::NTuple{D,Int}; anisotropic::Bool = false) where {D,T}
    _check_spatial_dimension(D)
    return _hybrid_fields(Val(D), T, nc, anisotropic)
end

function HybridFields{D,T}(nc::Tuple{Vararg{Int}}; anisotropic::Bool = false) where {D,T}
    _check_spatial_dimension(D)
    length(nc) == D ||
        throw(DimensionMismatch("grid size tuple length $(length(nc)) must equal D=$D"))
    nct = ntuple(d -> nc[d], Val(D))
    return _hybrid_fields(Val(D), T, nct, anisotropic)
end

# ---------------------------------------------------------------- moments
