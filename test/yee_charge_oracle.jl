module YeeChargeOracle

# Independent physical hat functions and clipped interval overlaps. Expected
# values do not invoke a package stencil, replay or rounding converter.
const Q = Rational{BigInt}
rat(x) = Q(x)

function charge(species, positions, n, length)
    values = fill(zero(Q), n + 1)
    magnitudes = fill(zero(Q), n + 1)
    for (ps, x) in zip(species, positions), p in eachindex(ps.weight), j = 0:n
        shape = max(zero(Q), one(Q) - abs(n * rat(x[p]) / rat(length) - j))
        contribution = rat(ps.q) * rat(ps.weight[p]) * shape
        values[j+1] += contribution
        magnitudes[j+1] += abs(contribution)
    end
    return values, magnitudes
end

function impulse(species, old, new, n, length)
    values = fill(zero(Q), n)
    magnitudes = fill(zero(Q), n)
    for (ps, x0, x1) in zip(species, old, new), p in eachindex(ps.weight), e = 0:(n-1)
        a, b = minmax(rat(x0[p]), rat(x1[p]))
        overlap = max(zero(Q), min(b, rat(length) * (e + 1) / n) - max(a, rat(length) * e / n))
        direction = x1[p] < x0[p] ? -1 : 1
        contribution = rat(ps.q) * rat(ps.weight[p]) * direction * overlap * n / rat(length)
        values[e+1] += contribution
        magnitudes[e+1] += abs(contribution)
    end
    return values, magnitudes
end

incidence(x) = [rat(x[1]); [rat(x[i]) - rat(x[i-1]) for i = 2:length(x)]; -rat(x[end])]

rounding_error(::Type{T}, exact) where {T} =
    rat(eps(T)) / 2 * abs(exact) + rat(nextfloat(zero(T))) / 2

function error_bound(::Type{T}, exact, magnitude, particles) where {T}
    u = rat(eps(Float64)) / 2
    gamma = (particles + 10) * u / (1 - (particles + 10) * u)
    native = gamma * magnitude
    return native + rounding_error(T, abs(exact) + native)
end

# Check a nearest-even rounding cell independently, including minsubnormal and
# finite/overflow neighbors. A nonzero reference also owns the sign of zero.
function rounded(value::T, exact::Q) where {T}
    isfinite(value) || return false
    iszero(exact) && return iszero(value)
    signbit(value) == (exact < 0) || return false
    x, reference = abs(value), abs(exact)
    lower = (rat(prevfloat(x)) + rat(x)) / 2
    upper =
        x == floatmax(T) ? rat(x) + (rat(x) - rat(prevfloat(x))) / 2 :
        (rat(x) + rat(nextfloat(x))) / 2
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, x))
    return (reference > lower || (even && reference == lower)) &&
           (reference < upper || (even && reference == upper))
end

end
