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

# Independent piecewise-polynomial quadrature (Boole, exact through degree5).
# Cloud references integrate density against physical hats/edges rather than
# evaluating a cumulative or using any package shape/replay/converter.
function polynomial_integral(f, knots)
    points = sort!(unique!(collect(knots)))
    total = zero(Q)
    for (left, right) in zip(points[1:end-1], points[2:end])
        h = (right - left) / 4
        total += (right - left) / 90 *
            (7 * f(left) + 32 * f(left + h) + 12 * f(left + 2h) + 32 * f(left + 3h) + 7 * f(right))
    end
    return total
end

density(r, x, a) = max(zero(Q), a - abs(r - x)) / a^2

function cloud_mass(r, x, a)
    hi = clamp(r, x - a, x + a)
    knots = [x-a, hi]
    x - a < x < hi && push!(knots, x)
    return polynomial_integral(t -> density(t, x, a), knots)
end

function cloud_nodes(x, n, length, a)
    x, length, a = rat(x), rat(length), rat(a)
    knots = [x-a, x, x+a]
    append!(knots, filter(r -> x - a < r < x + a, [length * j / n for j = 0:n]))
    return [polynomial_integral(knots) do r
        s = n * clamp(r, zero(Q), length) / length
        density(r, x, a) * max(zero(Q), 1 - abs(s - j))
    end for j = 0:n]
end

function cloud_edge(x, e, n, length, a)
    x, length, a = rat(x), rat(length), rat(a)
    lo, hi = length * e / n, length * (e + 1) / n
    knots = [lo, hi]
    append!(knots, filter(r -> lo < r < hi, [x-a, x, x+a]))
    return n / length * polynomial_integral(r -> density(r, x, a), knots)
end

function cloud_path(old, new, e, n, length, a)
    old, new, length, a = rat(old), rat(new), rat(length), rat(a)
    lo, hi = minmax(old, new)
    knots = [lo, hi]
    for wall in (length * e / n, length * (e + 1) / n), offset in (-a, zero(Q), a)
        point = wall + offset
        lo < point < hi && push!(knots, point)
    end
    value = polynomial_integral(x -> cloud_edge(x, e, n, length, a), knots)
    return new < old ? -value : value
end

function cloud_charge(species, positions, n, length, a)
    values = fill(zero(Q), n + 1)
    for (ps, x) in zip(species, positions), p in eachindex(ps.weight)
        values .+= rat(ps.q) * rat(ps.weight[p]) .* cloud_nodes(x[p], n, length, a)
    end
    return values
end

function cloud_impulse(species, old, new, n, length, a)
    return [sum((rat(ps.q) * rat(ps.weight[p]) *
        cloud_path(x0[p], x1[p], e, n, length, a)
        for (ps, x0, x1) in zip(species, old, new) for p in eachindex(ps.weight)); init = zero(Q))
        for e = 0:(n-1)]
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

function voltage(field, old, new, n, length, dx)
    a, b = minmax(rat(old), rat(new))
    total = zero(Q)
    for e = 0:(n-1)
        overlap = max(zero(Q), min(b, rat(length) * (e + 1) / n) - max(a, rat(length) * e / n))
        total += rat(field[e+1]) * overlap * n / rat(length)
    end
    return rat(dx) * (new < old ? -total : total)
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
