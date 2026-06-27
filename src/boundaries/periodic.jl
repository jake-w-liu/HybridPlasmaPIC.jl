# Periodic.jl — periodic particle boundary (from particles.jl)

"""
    apply_periodic!(ps, lo::NTuple{D}, hi::NTuple{D})

Wrap positions into `[lo, hi)` on every axis. Particle count is unchanged.
"""
function apply_periodic!(ps::ParticleSet{D,T}, lo::NTuple{D}, hi::NTuple{D}) where {D,T}
    loT, hiT = _validated_open_interval(lo, hi, T)
    @inbounds for d = 1:D
        l = loT[d]
        L = hiT[d] - loT[d]
        xd = ps.x[d]
        for p in eachindex(xd)
            xd[p] = l + mod(xd[p] - l, L)
        end
    end
    return ps
end
