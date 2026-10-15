# threaded.jl — CPU-threaded particle→grid deposition (§21.4).
#
# `deposit_scalar_threaded!` reproduces `deposit_scalar!` exactly (same stencil,
# same periodic wrap, same partition of unity) but spreads the particle loop over
# `Base.Threads.nthreads()` threads. Each thread accumulates into its OWN private
# grid array (no shared writes ⇒ no data races, no atomics); the per-thread
# partials are reduced into `out` at the end. With a single thread this is the
# serial algorithm with one extra (no-op) reduction. The only numerical
# difference vs serial is the order of floating-point additions in the per-cell
# sums, which is why callers compare with a small rtol rather than bit-equality.

"""
    deposit_scalar_threaded!(out, ps, vals, g, shape)

Threaded equivalent of [`deposit_scalar!`](@ref): accumulate
`out_g = Σ_p vals[p]·S_g(x_p)` over all particles using `Base.Threads`.
`out` is zeroed first. The result equals the serial deposition up to
floating-point reduction order (partition of unity ⇒ `Σ_g out_g = Σ_p vals[p]`).
The output must not overlap `vals` or any particle storage; validation precedes
all writes, as in the serial routine.
"""
function deposit_scalar_threaded!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    vals::AbstractVector,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    _validate_deposition_inputs(out, ps, vals, g)
    _deposit_scalar_threaded_kernel!(out, ps, (vals,), g, shape, Val(false))
    return out
end

function _deposit_scalar_threaded_kernel!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    factors::Tuple,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    checked::Val{check},
) where {D,T,check}
    Np = nparticles(ps)
    nth = Threads.nthreads()
    if nth == 1 || Np == 0
        return _deposit_scalar_kernel!(out, ps, factors, g, shape, checked)
    end
    fill!(out, zero(T))
    partials = [Array{T}(undef, g.n) for _ = 1:nth]
    # Each worker owns one grid and one flag record; reduction follows all jobs.
    states = check ? Vector{NamedTuple{(:eligible, :lost),Tuple{Bool,Bool}}}(undef, nth) : nothing
    chunks = _partition_ranges(Np, nth)
    Threads.@threads for t = 1:nth
        state = _deposit_scalar_kernel!(partials[t], ps, factors, g, shape, checked, chunks[t])
        states === nothing || (states[t] = state)
    end
    eligible = true
    lost = false
    for t = 1:nth
        if states !== nothing
            eligible &= states[t].eligible
            lost |= states[t].lost
        end
        acc = partials[t]
        for I in eachindex(out)
            out[I] += acc[I]
            check && (lost |= _deposition_range_lost(out[I]))
        end
    end
    return (; eligible, lost)
end

# Split 1:Np into `nparts` near-equal contiguous index ranges (the first
# `Np % nparts` ranges are one longer). Empty ranges (when nparts > Np) are
# returned as `1:0` so the worker simply skips them.
function _partition_ranges(Np::Int, nparts::Int)
    ranges = Vector{UnitRange{Int}}(undef, nparts)
    base, rem = divrem(Np, nparts)
    start = 1
    @inbounds for t = 1:nparts
        len = base + (t <= rem ? 1 : 0)
        ranges[t] = start:(start+len-1)
        start += len
    end
    return ranges
end

"""
    density_threaded!(nout, ps, g, shape)

Threaded equivalent of [`density!`](@ref): number density
`n_g = (Σ_p w_p S_g)/ΔV`, computed with [`deposit_scalar_threaded!`](@ref).
"""
function density_threaded!(
    nout::Array{T,D},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    _validate_deposition_inputs(nout, ps, ps.weight, g)
    factors = (ps.weight,)
    state = _deposit_scalar_threaded_kernel!(nout, ps, factors, g, shape, Val(true))
    _finish_moment!(nout, nout, ps, factors, g, shape, one(T), state)
    return nout
end
