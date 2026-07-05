# OpenFieldBoundary.jl — absorbing particle boundary (from particles.jl)

"""
    apply_absorbing!(ps, lo::NTuple{D}, hi::NTuple{D}) -> nremoved

Remove particles that left the box `[lo, hi)` on any axis (compacting all
arrays). Returns the number removed — the only sanctioned particle loss.
"""
function apply_absorbing!(ps::ParticleSet{D,T}, lo::NTuple{D}, hi::NTuple{D}) where {D,T}
    N = nparticles(ps)
    loT, hiT = _validated_open_interval(lo, hi, T)
    @inbounds for d = 1:D
        xd = ps.x[d]
        for p in eachindex(xd)
            isfinite(xd[p]) || throw(ArgumentError("particle position x[$d][$p] must be finite"))
        end
    end
    write = 0
    @inbounds for p = 1:N
        inside = true
        for d = 1:D
            xp = ps.x[d][p]
            inside &= loT[d] <= xp < hiT[d]
        end
        if inside
            write += 1
            if write != p
                for d = 1:D
                    ps.x[d][write] = ps.x[d][p]
                end
                for c = 1:3
                    ps.v[c][write] = ps.v[c][p]
                end
                ps.weight[write] = ps.weight[p]
                ps.id[write] = ps.id[p]
                ps.tag[write] = ps.tag[p]
            end
        end
    end
    if write < N
        for d = 1:D
            resize!(ps.x[d], write)
        end
        for c = 1:3
            resize!(ps.v[c], write)
        end
        resize!(ps.weight, write)
        resize!(ps.id, write)
        resize!(ps.tag, write)
    end
    return N - write
end
