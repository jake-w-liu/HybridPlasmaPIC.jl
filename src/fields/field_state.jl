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
    pforce::NTuple{3,Array{T,D}}      # frozen ∇·P_e (anisotropic/CGL closure; unused/0 for scalar)
    floor_count::Base.RefValue{Int}
end

function HybridFields{D,T}(nc::NTuple{D,Int}) where {D,T}
    z() = zeros(T, nc)
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
        ntuple(_ -> z(), 3),
        Ref(0),
    )
end

# ---------------------------------------------------------------- moments
