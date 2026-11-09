# Periodic.jl — periodic particle boundary (from particles.jl)

# The caller supplies a validated interval and its loop-invariant width.
# Floating-point mod or the translated result can round to hi; all periodic
# consumers use the same half-open interval convention in that case.
@inline function _periodic_coordinate(x, lo, hi, width)
    wrapped = lo + mod(x - lo, width)
    return wrapped >= hi ? lo : wrapped
end

# Canonical original-Q periodic endpoint materialization. Transport bounds and
# history capacity belong to the caller, not to this scalar wrap/lift map.
function _periodic_represented_lift(
    ideal::Rational{BigInt},
    period::T,
) where {T<:Union{Float32,Float64}}
    isfinite(period) && _stored_positive(period) ||
        throw(ArgumentError("periodic represented lift requires a finite positive period"))
    L = _exact_rational(period)
    wrapped = _periodic_coordinate(ideal, zero(L), L, L)
    winding = BigInt((ideal - wrapped) / L)
    rounded = _from_exact(T, wrapped)
    rounded == period && (winding += 1)
    represented = _periodic_coordinate(rounded, zero(T), period, period)
    lifted = _exact_rational(represented) + winding * L
    return (wrapped = represented, lifted = lifted, winding = winding, defect = lifted - ideal)
end

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
            xd[p] = _periodic_coordinate(xd[p], l, h, L)
        end
    end
    return ps
end
