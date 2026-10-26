module SourcedYeeOracle

# Independent matrix topology and rational expressions; no package-private calls.
rat(x) = Rational{BigInt}(x)

function curl(v, n, periodic, electric, dx)
    out = ntuple(3) do c
        nodes = electric ? c == 1 : c != 1
        zeros(Rational{BigInt}, n + Int(nodes && !periodic))
    end
    for (c, source, sign) in ((2, 3, -1), (3, 2, 1)), i in eachindex(out[c])
        if electric
            j = periodic && i == n ? 1 : i + 1
            out[c][i] = sign * (rat(v[source][j]) - rat(v[source][i])) / rat(dx)
        elseif periodic
            j = i == 1 ? n : i - 1
            out[c][i] = sign * (rat(v[source][i]) - rat(v[source][j])) / rat(dx)
        elseif 1 < i <= n
            out[c][i] = sign * (rat(v[source][i]) - rat(v[source][i-1])) / rat(dx)
        end
    end
    return out
end

function step(E, B, impulse, n, periodic, dx, h, speed)
    T = eltype(E[1])
    ce = curl(E, n, periodic, true, dx)
    half = ntuple(c -> T.(rat.(B[c]) .- rat(h) / 2 .* ce[c]), 3)
    cb = curl(half, n, periodic, false, dx)
    en = ntuple(c -> T.(rat.(E[c]) .+ rat(h) * rat(speed)^2 .* cb[c] .- rat.(impulse[c])), 3)
    cen = curl(en, n, periodic, true, dx)
    bn = ntuple(c -> T.(rat.(half[c]) .- rat(h) / 2 .* cen[c]), 3)
    return en, bn, half
end

function charge(E, periodic)
    x = rat.(E[1])
    return periodic ? [x[i] - x[mod1(i - 1, length(x))] for i in eachindex(x)] :
           [x[1]; diff(x); -x[end]]
end

function work(old, new, impulse, dx)
    return rat(dx) / 2 * sum(
        rat(impulse[c][i]) * (rat(old[c][i]) + rat(new[c][i])) for c = 1:3 for
        i in eachindex(old[c])
    )
end

# Normwise scales: the same sums with every contribution taken in magnitude.
function work_scale(old, new, impulse, dx)
    return rat(dx) / 2 * sum(
        abs(rat(impulse[c][i])) * (abs(rat(old[c][i])) + abs(rat(new[c][i]))) for c = 1:3 for
        i in eachindex(old[c])
    )
end

function energy_scale(E, B, n, periodic, dx, h, speed)
    ce = curl(E, n, periodic, true, dx)
    e2 = sum(rat(x)^2 for a in E for x in a)
    b2 = sum(rat(x)^2 for a in B for x in a)
    c2 = sum(x^2 for a in ce for x in a)
    return rat(dx) / 2 * (e2 + rat(speed)^2 * b2 + rat(h)^2 * rat(speed)^2 * c2 / 4)
end

function energy(E, B, n, periodic, dx, h, speed)
    ce = curl(E, n, periodic, true, dx)
    e2 = sum(rat(x)^2 for a in E for x in a)
    b2 = sum(rat(x)^2 for a in B for x in a)
    c2 = sum(x^2 for a in ce for x in a)
    return rat(dx) / 2 * (e2 + rat(speed)^2 * b2 - rat(h)^2 * rat(speed)^2 * c2 / 4)
end

end
