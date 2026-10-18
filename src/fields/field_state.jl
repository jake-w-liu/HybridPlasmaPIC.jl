# FieldState.jl — HybridFields container (from hybrid.jl)

"""
    HybridFields{D,T}(n::NTuple{D,Int})

Grid-resident hybrid fields and scratch: density `n`, ion bulk velocity `ui`,
magnetic `B`, electric `E`, current `J`, electron pressure `pe`, its gradient
`gradp`, reciprocal density `ninv`, and a density-floor activation counter.
For [`EvolvingElectrons`](@ref), `pe` is persistent integer-time state alongside
B; algebraic closures use it as derived pressure storage.
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
    cgl_dp::Array{T,D}                # p∥-p⊥ scratch (anisotropic/CGL only; 0-length for scalar)
    floor_count::Base.RefValue{Int}
end

function _hybrid_fields(::Val{D}, ::Type{T}, nc::NTuple{D,Int}, anisotropic::Bool) where {D,T}
    z() = zeros(T, nc)
    cgl() = anisotropic ? zeros(T, nc) : zeros(T, ntuple(_ -> 0, D))
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
        ntuple(_ -> cgl(), 3),
        cgl(),
        Ref(0),
    )
end

# `anisotropic=true` (a CGL closure) allocates the `pforce` ∇·P_e buffers and `cgl_dp`
# scratch full-size; scalar closures get 0-length CGL arrays (type-stable, never indexed on
# the scalar Ohm path) so the common case carries no dead weight.
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

function _copy_hybrid_fields!(destination::HybridFields, candidate::HybridFields)
    # Preserve the arrays exposed to callers while committing the physical and
    # derived fields together. All fallible calculations happen on the candidate.
    for (dst, src) in (
        (destination.n, candidate.n),
        (destination.pe, candidate.pe),
        (destination.ninv, candidate.ninv),
        (destination.cgl_dp, candidate.cgl_dp),
    )
        copyto!(dst, src)
    end
    for (dst, src) in (
        (destination.ui, candidate.ui),
        (destination.B, candidate.B),
        (destination.E, candidate.E),
        (destination.J, candidate.J),
        (destination.lapJ, candidate.lapJ),
        (destination.pforce, candidate.pforce),
    )
        for c in eachindex(dst)
            copyto!(dst[c], src[c])
        end
    end
    for d in eachindex(destination.gradp)
        copyto!(destination.gradp[d], candidate.gradp[d])
    end
    destination.floor_count[] = candidate.floor_count[]
    return nothing
end
