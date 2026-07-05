# Periodic.jl — periodic particle boundary (from particles.jl)

"""
    apply_periodic!(ps, lo::NTuple{D}, hi::NTuple{D})

Wrap positions into `[lo, hi)` on every axis. Particle count is unchanged.
"""
function apply_periodic!(ps::ParticleSet{D,T}, lo::NTuple{D}, hi::NTuple{D}) where {D,T}
    loT, hiT = _validated_open_interval(lo, hi, T)
    @inbounds for d = 1:D
        xd = ps.x[d]
        for p in eachindex(xd)
            isfinite(xd[p]) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
        end
    end
    @inbounds for d = 1:D
        l = loT[d]
        h = hiT[d]
        L = h - l
        xd = ps.x[d]
        for p in eachindex(xd)
            xd[p] = l + mod(xd[p] - l, L)
            # float rounding can land exactly on hi (mod can return L, and l + m
            # can round up to hi); x==hi is the same physical point as lo — fold
            # it back so the half-open box [lo,hi) contract holds.
            xd[p] >= h && (xd[p] = l)
        end
    end
    return ps
end
