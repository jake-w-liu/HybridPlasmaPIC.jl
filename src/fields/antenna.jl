# antenna.jl — external-antenna / RF field source (Phase-4 forcing).
#
# An antenna drives waves by imposing an external electric field E_ant(x,t) (e.g. a
# localized oscillating profile). Its effect on the magnetic field over dt is the Faraday
# contribution ∂B/∂t = −∇×E_ant, so this adds −dt·∇×E_ant to B. Being a curl, the
# injection is divergence-free — it preserves ∇·B = 0 exactly. (Uniform / DC applied
# fields and an RF drive on the particles are already covered by `push_uniform!` called
# with a time-varying E; this is the spatially-structured, wave-launching antenna.)

"""
    apply_antenna!(B, E_ant, dt, g; work=nothing) -> B

Inject one substep of an external antenna into the magnetic field: `B ← B − dt·∇×E_ant`,
the Faraday response to a prescribed external electric field `E_ant` (a 3-tuple of grid
arrays; supply e.g. `E_ant = amp·profile(x)·sin(ω·t)` afresh each step). Because the
update is a curl it is divergence-free, so `∇·B` is unchanged. Combined with a field
solver that propagates `B`, a localized oscillating `E_ant` launches waves at its
frequency. `dt ≥ 0`. `work` is an optional 3-tuple of scratch arrays (allocated if
omitted). Returns `B`.
"""
function apply_antenna!(
    B::NTuple{3,<:Array{T,D}},
    E_ant::NTuple{3,<:Array{T,D}},
    dt::Real,
    g::FourierGrid{D,T};
    work::Union{Nothing,NTuple{3,Array{T,D}}} = nothing,
) where {D,T}
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    dtT = _require_finite_nonnegative_real("dt", dt, T)
    dtT == 0 && return B
    for c = 1:3
        size(B[c]) == g.n ||
            throw(DimensionMismatch("B[$c] size $(size(B[c])) does not match grid size $(g.n)"))
    end
    curlE = work === nothing ? ntuple(_ -> similar(B[1]), 3) : work
    curl!(curlE, E_ant, g)                       # ∇×E_ant
    @inbounds for c = 1:3
        Bc = B[c]
        cc = curlE[c]
        for I in eachindex(Bc)
            Bc[I] -= dtT * cc[I]                 # B ← B − dt ∇×E_ant  (Faraday, div-free)
        end
    end
    return B
end
