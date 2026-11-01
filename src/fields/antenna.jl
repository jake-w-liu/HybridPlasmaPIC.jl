# antenna.jl — external-antenna / RF field source (Phase-4 forcing).
#
# An antenna drives waves by imposing an external electric field E_ant(x,t) (e.g. a
# localized oscillating profile). Its effect on the magnetic field over dt is the Faraday
# contribution ∂B/∂t = −∇×E_ant, so this adds −dt·∇×E_ant to B. Being a curl, the
# injection preserves discrete divergence to roundoff. (Uniform / DC applied
# fields and an RF drive on the particles are already covered by `push_uniform!` called
# with a time-varying E; this is the spatially-structured, wave-launching antenna.)

"""
    apply_antenna!(B, E_ant, dt, g; work=nothing) -> B

Inject one substep of an external antenna into the magnetic field: `B ← B − dt·∇×E_ant`,
the Faraday response to a prescribed external electric field `E_ant` (a 3-tuple of grid
arrays; supply e.g. `E_ant = amp·profile(x)·sin(ω·t)` afresh each step). Because the
update is a curl it preserves discrete divergence to roundoff. Combined with a field
solver that propagates `B`, a localized oscillating `E_ant` launches waves at its
frequency. `dt ≥ 0`. `work` is an optional 3-tuple of scratch arrays (allocated if
omitted). Its components must be independent of every `B` component, and `B`
components must be pairwise independent. Neither `B` nor `work` may overwrite
read-only Fourier metadata. These storage conflicts reject before mutation.
`B` may share storage with `E_ant`: the complete curl is formed before updating
`B`. Returns `B`.
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
    for c = 1:2, d = (c+1):3
        Base.mightalias(B[c], B[d]) &&
            throw(ArgumentError("antenna magnetic field components must not alias each other"))
    end
    for scratch in curlE, field in B
        Base.mightalias(scratch, field) &&
            throw(ArgumentError("antenna workspace must not alias the magnetic field"))
    end
    for output in (B..., curlE...),
        metadata in (g.midx, g.kfull, g.ik, g.kvec), axis in metadata

        Base.mightalias(output, axis) &&
            throw(ArgumentError("antenna outputs must not alias Fourier metadata"))
    end
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
