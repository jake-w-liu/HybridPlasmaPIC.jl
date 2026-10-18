using HybridPlasmaPIC, Test, Random

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

@testset "CGL complete norm, pressure and volume ranges" begin
    O = CGLEnergyRangeOracle
    for T in (Float32, Float64), D = 1:3
        k = T === Float32 ? 100 : 600
        large, small = ldexp(one(T), k), ldexp(one(T), -k)
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> one(T), D))
        for amplitude in (large, small)
            B = (fill(amplitude, g.n), zeros(T, g.n), zeros(T, g.n))
            for (n, c) in (
                (ones(T, g.n), CGLElectrons(inv(amplitude), zero(T), one(T), one(T))),
                (fill(amplitude, g.n), CGLElectrons(zero(T), one(T), one(T), one(T))),
            )
                exact = O.interval(n, B, c, g.dx)
                before = map(copy, (n, B...))
                @test !isfinite(O.ordinary(n, B, c, g.dx)) || iszero(O.ordinary(n, B, c, g.dx))
                value = electron_internal_energy(n, B, c, g)
                @test isfinite(value) && value > 0
                @test O.rounded(value, exact...)
                @test all(isequal.(before, (n, B...)))
            end
        end
        # Two nonzero components require an irrational norm; volume and
        # invariant coefficients compensate their square-range excursion.
        for sign in (-one(T), one(T)), inverse in (false, true)
            amplitude = inverse ? small : large
            n = fill(sign, g.n)
            B = (fill(amplitude, g.n), fill(amplitude, g.n), zeros(T, g.n))
            c = CGLElectrons(inv(amplitude), zero(T), one(T), one(T))
            @test O.rounded(electron_internal_energy(n, B, c, g), O.interval(n, B, c, g.dx)...)
        end
    end
    for T in (Float32, Float64)
        k = T === Float32 ? 80 : 600
        large, small = ldexp(one(T), k), ldexp(one(T), -k)
        for lengths in ((large, large, small), (small, small, large), (small, large, small))
            g = FourierGrid((2, 2, 2), lengths)
            n = ones(T, g.n)
            B = (ones(T, g.n), zeros(T, g.n), zeros(T, g.n))
            c = CGLElectrons(T(1), T(1), T(1), T(1))
            @test CGLEnergyRangeOracle.rounded(
                electron_internal_energy(n, B, c, g),
                CGLEnergyRangeOracle.interval(n, B, c, g.dx)...,
            )
        end
    end
end

@testset "CGL equivalent radicals, midpoint tails and output overflow" begin
    O = CGLEnergyRangeOracle
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(4),))
        n = T[2, -1, 0, 0]
        k = T === Float32 ? 80 : 550
        a = ldexp(one(T), k)
        B = (T[a, 2a, 0, 0], T[a, 2a, 0, 0], zeros(T, 4))
        c = CGLElectrons(inv(a), zero(T), one(T), one(T))
        @test electron_internal_energy(n, B, c, g) === zero(T)
        c2 = CGLElectrons(inv(a), one(T), one(T), one(T))
        @test O.rounded(electron_internal_energy(n, B, c2, g), O.interval(n, B, c2, g.dx)...)
        # A rational norm makes exact rounding-cell boundaries independent of
        # the square-root interval. Very small fields exercise norm underflow.
        q = nextfloat(zero(T))
        B = (fill(q, 4), zeros(T, 4), zeros(T, 4))
        c = CGLElectrons(one(T), zero(T), one(T), one(T))
        for density in (T(0.5), T(1.5), T(2.5))
            fill!(n, zero(T))
            n[1] = density
            @test O.rounded(electron_internal_energy(n, B, c, g), O.interval(n, B, c, g.dx)...)
        end
        fill!(n, one(T))
        B = (fill(floatmax(T), 4), zeros(T, 4), zeros(T, 4))
        @test electron_internal_energy(n, B, c, g) === T(Inf)
        fill!(n, -one(T))
        @test electron_internal_energy(n, B, c, g) === T(-Inf)
    end
end

@testset "CGL ordinary output, allocation, singularity and nonfinite behavior" begin
    O = CGLEnergyRangeOracle
    rng = MersenneTwister(192841)
    for T in (Float32, Float64), C in (Float32, Float64), D = 1:3, mixed in (false, true)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        n = T(0.5) .+ rand(rng, T, g.n)
        mixed && (n[1:2:end] .*= -one(T))
        B = ntuple(_ -> T(0.5) .+ rand(rng, T, g.n), 3)
        c = CGLElectrons(C(0.7), C(1.1), one(C), one(C))
        expected = O.ordinary(n, B, c, g.dx)
        before = map(copy, (n, B...))
        @test isequal(electron_internal_energy(n, B, c, g), expected)
        @test O.allocated(n, B, c, g) == 0
        @test all(isequal.(before, (n, B...)))
    end
    for T in (Float32, Float64)
        g = FourierGrid((2,), (one(T),))
        n = ones(T, 2)
        B = ntuple(_ -> zeros(T, 2), 3)
        for ppar in (zero(T), one(T)), density in (zero(T), one(T), -one(T))
            c = CGLElectrons(one(T), ppar, one(T), one(T))
            fill!(n, density)
            @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
        end
        c = CGLElectrons(one(T), one(T), one(T), one(T))
        n .= (one(T), -one(T))
        @test isnan(electron_internal_energy(n, B, c, g))
        for value in (T(Inf), T(-Inf), T(NaN))
            fill!(n, one(T))
            fill!(B[1], one(T))
            n[1] = value
            @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
            fill!(n, one(T))
            B[1][1] = value
            @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
        end
        @test_throws ArgumentError electron_internal_energy(n, B, c, g; pe = ones(T, 2))
        @test isnan(electron_internal_energy(n, c, g))
    end
end

@testset "CGL mixed-precision complete factors and independent error enclosures" begin
    O = CGLEnergyRangeOracle
    Q = O.Q
    rng = MersenneTwister(391887)
    checked = 0
    accepted = 0
    g = FourierGrid((2,), (1.0,))
    n = fill(ldexp(1.0, 200), 2)
    B = (ones(2), zeros(2), zeros(2))
    c = CGLElectrons(ldexp(1.0f0, -100), 0.0f0, 1.0f0, 1.0f0)
    @test O.rounded(electron_internal_energy(n, B, c, g), O.interval(n, B, c, g.dx)...)
    for T in (Float32, Float64), C in (Float32, Float64), count in (2, 4, 16), trial = 1:16
        g = FourierGrid((count,), (one(T),))
        n = zeros(T, count)
        B = ntuple(_ -> zeros(T, count), 3)
        for i = 1:count
            n[i] = ldexp(T(0.75) + rand(rng, T) / T(2), rand(rng, (-20, -5, 0, 5, 20)))
            isodd(trial) && rand(rng, Bool) && (n[i] = -n[i])
            for a in B
                a[i] = ldexp(T(0.75) + rand(rng, T) / T(2), rand(rng, (-30, -5, 0, 5, 30)))
            end
        end
        c = CGLElectrons{C}(
            ldexp(one(C), rand(rng, (-30, 0, 30))),
            ldexp(one(C), rand(rng, (-30, 0, 30))),
        )
        value, absolute, normal, active, largest, nonnegative =
            HybridPlasmaPIC._cgl_energy_probe(n, B, c, g.dx, Val(true))
        lo, hi = O.interval(n, B, c, g.dx; bits = 512)
        certified =
            normal && HybridPlasmaPIC._cgl_energy_normal_certified(value, absolute, count, C)
        if normal && nonnegative && !iszero(value)
            volume = prod(Q, g.dx)
            @test Q(largest) * volume / (BigInt(1) << 19) <= lo
            @test hi <= Q(largest) * volume * (2count) * (BigInt(1) << 19)
            if HybridPlasmaPIC._cgl_energy_positive_certified(largest, g.dx, count)
                @test Q(floatmin(T)) <= lo <= hi < Q(floatmax(T))
                @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
            end
        end
        if normal
            checked += 1
            p = min(precision(T), precision(C))
            delta = Q(count + 22) / (BigInt(1) << p)
            error = delta / (1 - 2delta) * Q(absolute)
            @test Q(value) - error <= lo <= hi <= Q(value) + error
        end
        if certified
            accepted += 1
            @test isequal(electron_internal_energy(n, B, c, g), O.ordinary(n, B, c, g.dx))
            @test Q(floatmin(T)) <= min(abs(lo), abs(hi))
            @test max(abs(lo), abs(hi)) < Q(floatmax(T))
        else
            @test O.rounded(electron_internal_energy(n, B, c, g), lo, hi)
        end
    end
    @test checked > 0
    @test accepted > 0
end
