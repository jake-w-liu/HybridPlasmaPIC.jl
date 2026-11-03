module CGLEnergyRangeOracle
using HybridPlasmaPIC
const Q = Rational{BigInt}
function interval(n, B, closure, dx; bits = 2048)
    volume = prod(Q, dx)
    lower = upper = Q(0)
    scale = BigInt(1) << bits
    for I in eachindex(n)
        density = Q(n[I])
        iszero(density) && continue
        r = sum(Q(a[I])^2 for a in B)
        if iszero(r)
            iszero(closure.Cpar) || throw(DomainError(r, "singular parallel pressure"))
            continue
        end
        rational = Q(closure.Cpar) * density^3 * volume / (2r)
        lower += rational
        upper += rational
        a = Q(closure.Cperp) * density * volume
        iszero(a) && continue
        sn, sd = isqrt(numerator(r)), isqrt(denominator(r))
        if sn^2 == numerator(r) && sd^2 == denominator(r)
            exact = a * (sn // sd)
            lower += exact
            upper += exact
        else
            k = isqrt(div(numerator(r) << (2bits), denominator(r)))
            lo, hi = k // scale, (k + 1) // scale
            lower += a * (a > 0 ? lo : hi)
            upper += a * (a > 0 ? hi : lo)
        end
    end
    return lower, upper
end
function rounded(value::T, lo, hi) where {T}
    lo <= hi || return false
    isfinite(value) || isinf(value) || return false
    edge = Q(floatmax(T)) + (Q(floatmax(T)) - Q(prevfloat(floatmax(T)))) / 2
    if isinf(value)
        return signbit(value) ? hi <= -edge : lo >= edge
    end
    if signbit(value)
        value = abs(value)
        lo, hi = -hi, -lo
    end
    center = Q(value)
    lower = (Q(prevfloat(value)) + center) / 2
    upper = isfinite(nextfloat(value)) ? (center + Q(nextfloat(value))) / 2 : edge
    U = T === Float32 ? UInt32 : UInt64
    even = iseven(reinterpret(U, value))
    return (lower < lo || even && lower == lo) && (hi < upper || even && hi == upper)
end
function ordinary(n::Array{T}, B, closure, dx) where {T}
    s = zero(T)
    for I in eachindex(n, B[1], B[2], B[3])
        bm = sqrt(B[1][I]^2 + B[2][I]^2 + B[3][I]^2)
        s += T(cgl_pperp(closure, n[I], bm)) + T(0.5) * T(cgl_ppar(closure, n[I], bm))
    end
    return s * prod(dx)
end
function allocated(n, B, c, g)
    electron_internal_energy(n, B, c, g)
    return @allocated electron_internal_energy(n, B, c, g)
end
end
