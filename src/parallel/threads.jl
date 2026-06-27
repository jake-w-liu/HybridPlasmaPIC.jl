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
"""
function deposit_scalar_threaded!(
    out::Array{T,D},
    ps::ParticleSet{D,T},
    vals::AbstractVector,
    g::FourierGrid{D,T},
    shape::ShapeFunction,
) where {D,T}
    length(vals) == nparticles(ps) || throw(DimensionMismatch("vals length ≠ particle count"))
    size(out) == g.n ||
        throw(DimensionMismatch("out size $(size(out)) does not match grid size $(g.n)"))
    fill!(out, zero(T))
    Np = length(vals)
    val_first = firstindex(vals)
    nth = Threads.nthreads()
    n = g.n
    stamp = CartesianIndices(ntuple(_ -> width(shape), D))

    # Fall back to the in-place serial path when threading buys nothing: a single
    # thread or a particle count too small to bother allocating partials for.
    if nth == 1 || Np == 0
        @inbounds for p = 1:Np
            st = ntuple(d -> _stencil1d(shape, ps.x[d][p] / g.dx[d]), D)
            val = vals[val_first+p-1]
            for c in stamp
                o = Tuple(c)
                w = val
                for d = 1:D
                    w *= st[d][2][o[d]]
                end
                idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
                out[idx...] += w
            end
        end
        return out
    end

    # One private accumulator per thread, zero-initialised.
    partials = [zeros(T, n) for _ = 1:nth]

    # Static contiguous partition of the particle index range across threads, so
    # every particle is deposited exactly once. `@threads :static` would also
    # work but we partition explicitly so the mapping is independent of the
    # scheduler and each thread writes only its own `partials[tid]`.
    chunks = _partition_ranges(Np, nth)
    Threads.@threads for t = 1:nth
        rng = chunks[t]
        isempty(rng) && continue
        acc = partials[t]   # private to this thread: no other thread touches it
        @inbounds for p in rng
            st = ntuple(d -> _stencil1d(shape, ps.x[d][p] / g.dx[d]), D)
            val = vals[val_first+p-1]
            for c in stamp
                o = Tuple(c)
                w = val
                for d = 1:D
                    w *= st[d][2][o[d]]
                end
                idx = ntuple(d -> mod(st[d][1] + o[d] - 1, n[d]) + 1, D)
                acc[idx...] += w
            end
        end
    end

    # Reduce the per-thread partials into `out`.
    @inbounds for t = 1:nth
        acc = partials[t]
        for I in eachindex(out)
            out[I] += acc[I]
        end
    end
    return out
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
    deposit_scalar_threaded!(nout, ps, ps.weight, g, shape)
    nout ./= prod(g.dx)
    return nout
end
