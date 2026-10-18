using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function original(J, de2::T, dx) where {T}
    s = zero(T)
    for I in eachindex(J[1], J[2], J[3])
        s += J[1][I]^2 + J[2][I]^2 + J[3][I]^2
    end
    return T(0.5) * de2 * s * prod(dx)
end
exact(J, de2, dx) = Q(de2) * sum(sum(Q(x)^2 for x in a) for a in J) * prod(Q, dx) / 2
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
    lo = (Q(prevfloat(x)) + Q(x)) / 2
    hi = isfinite(nextfloat(x)) ? (Q(x) + Q(nextfloat(x))) / 2 : edge
    even = T === Float32 ? iseven(reinterpret(UInt32, x)) : iseven(reinterpret(UInt64, x))
    return lo < value < hi || even && (lo == value || hi == value)
end
cost(J, d, dx) = @allocated H._electron_inertia_energy(J, d, dx)
