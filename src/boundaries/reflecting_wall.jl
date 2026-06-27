# ReflectingWall.jl — reflecting particle boundary (from particles.jl)

"""
    apply_reflecting!(ps, lo::NTuple{D}, hi::NTuple{D})

Specularly reflect particles that crossed a wall: position is mirrored back into
the box and the normal velocity component is flipped. Count is unchanged.
Assumes at most one reflection per step (true when |v·dt| < box length).
"""
function apply_reflecting!(ps::ParticleSet{D,T}, lo::NTuple{D}, hi::NTuple{D}) where {D,T}
    @inbounds for d = 1:D
        l = T(lo[d])
        h = T(hi[d])
        xd = ps.x[d]
        vd = ps.v[d]
        for p in eachindex(xd)
            if xd[p] < l
                xd[p] = 2l - xd[p]
                vd[p] = -vd[p]
            elseif xd[p] >= h
                xd[p] = 2h - xd[p]
                # exact x==h maps back to h (outside the half-open box [l,h));
                # nudge strictly inside so cell_index and deposit_scalar! agree.
                xd[p] >= h && (xd[p] = prevfloat(h))
                vd[p] = -vd[p]
            end
        end
    end
    return ps
end
