module KineticRangeOracle
using HybridPlasmaPIC
Q = Rational{BigInt}
function exact(species)
    result = Q(0)
    for p in species, i in eachindex(p.weight)
        result += Q(p.m) * Q(p.weight[i]) * sum(Q(v[i])^2 for v in p.v) / 2
    end
    return result
end
function ordinary(p::ParticleSet{D,T}) where {D,T}
    total = zero(T)
    for i in eachindex(p.weight)
        total += p.weight[i] * (p.v[1][i]^2 + p.v[2][i]^2 + p.v[3][i]^2)
    end
    return T(0.5) * p.m * total
end
function rounded(x::T, value) where {T}
    isfinite(x) || isinf(x) || return false
    largest = Q(floatmax(T))
    edge = largest + (largest - Q(prevfloat(floatmax(T)))) / 2
    if isinf(x)
        return signbit(x) ? value <= -edge : value >= edge
    end
    value = signbit(x) ? -value : value
    x = abs(x)
    if iszero(x)
        return 0 <= value <= Q(nextfloat(zero(T))) / 2
    end
    low = (Q(prevfloat(x)) + Q(x)) / 2
    high = isfinite(nextfloat(x)) ? (Q(x) + Q(nextfloat(x))) / 2 : edge
    U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
    even = iseven(reinterpret(U, x))
    return low < value < high || even && (low == value || high == value)
end
end
function kinetic_range_alloc(ps)
    kinetic_energy(ps)
    return @allocated kinetic_energy(ps)
end
