# antenna.jl — external-antenna / RF field source (Phase-4 forcing).
#
# An antenna drives waves by imposing an external electric field E_ant(x,t) (e.g. a
# localized oscillating profile). Its effect on the magnetic field over dt is the Faraday
# contribution ∂B/∂t = −∇×E_ant, so this adds −dt·∇×E_ant to B. Being a curl, the
# injection preserves discrete divergence to roundoff. (Uniform / DC applied
# fields and an RF drive on the particles are already covered by `push_uniform!` called
# with a time-varying E; this is the spatially-structured, wave-launching antenna.)

@inline _antenna_update_native(b, curl, dt) = b - dt * curl

# Public admission supplies finite operands in nearest/gradual arithmetic.
# The half-maximum margin certifies only the overflow domain, not exact rounding.
@inline function _antenna_update_value(b::T, curl::T, dt::T) where {T}
    value = _antenna_update_native(b, curl, dt)
    isfinite(value) && (abs(value) <= floatmax(T) / 2 || iszero(curl) || iszero(dt)) && return value
    return _antenna_update_repair(value, b, curl, dt)
end

@noinline function _antenna_update_repair(value::T, b::T, curl::T, dt::T) where {T}
    Q = _exact_rational
    rounded = _from_exact(T, _antenna_update_native(Q(b), Q(curl), Q(dt)))
    isfinite(rounded) || throw(ArgumentError("complete antenna magnetic field is not finite in $T"))
    return isfinite(value) ? value : rounded
end

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
read-only Fourier metadata. Storage checks use `Base.mightalias`; callers must
supply arrays whose shared storage it reports. Detected conflicts reject before mutation.
`B` may share storage with `E_ant`: the complete curl is formed before updating
`B`. A nonzero converted timestep requires nearest rounding, gradual underflow,
finite `B` and a finite materialized curl. The update retains ordinary native
bits and repairs affine range events from the stored `B`, `dt` and curl operands.
All finite candidates are staged in `work` before any `B` write; numerical
failure may change scratch but leaves `B` unchanged. An unrepresentable curl is
rejected even if a later timestep factor could yield a finite update. A true
converted signed-zero timestep remains a no-op. Returns `B`.
"""
function apply_antenna!(
    B::NTuple{3,<:Array{T,D}},
    E_ant::NTuple{3,<:Array{T,D}},
    dt::Real,
    g::FourierGrid{D,T};
    work::Union{Nothing,NTuple{3,Array{T,D}}} = nothing,
) where {D,T}
    dt >= 0 || throw(ArgumentError("dt must be ≥ 0"))
    dtT = _validated_nonnegative_dt(T, dt; name = "apply_antenna!")
    _particle_stored_zero(dtT) && return B
    rounding(T) === RoundNearest && !get_zero_subnormals() ||
        throw(ArgumentError("apply_antenna! requires nearest rounding and gradual underflow"))
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
    for output in (B..., curlE...), metadata in (g.midx, g.kfull, g.ik, g.kvec), axis in metadata

        Base.mightalias(output, axis) &&
            throw(ArgumentError("antenna outputs must not alias Fourier metadata"))
    end
    curl!(curlE, E_ant, g)                       # ∇×E_ant
    @inbounds for c = 1:3
        Bc = B[c]
        cc = curlE[c]
        for I in eachindex(Bc)
            b, curl = Bc[I], cc[I]
            isfinite(b) || throw(ArgumentError("antenna magnetic field B[$c][$I] must be finite"))
            isfinite(curl) || throw(
                ArgumentError(
                    "antenna requires a finite materialized curl in $T; " *
                    "completion through a later timestep factor is not supported",
                ),
            )
            cc[I] = _antenna_update_value(b, curl, dtT)
        end
    end
    for c = 1:3
        copyto!(B[c], curlE[c])
    end
    return B
end
