using HybridPlasmaPIC, Test, Random

module EMFieldEnergyOracle
using HybridPlasmaPIC: EMPIC1D, EMPIC
exact(x) = Rational{BigInt}(x)
components(es::EMPIC1D) = ((es.Ex, es.Ey), (es.Bz,))
components(es::EMPIC) = (es.E, es.B)
function integral(es)
    electric, magnetic = components(es)
    volume = prod(exact.(es.g.dx))
    value = zero(Rational{BigInt})
    for a in electric, x in a
        value += exact(x)^2 * volume / 2
    end
    for a in magnetic, x in a
        value += (exact(es.c) * exact(x))^2 * volume / 2
    end
    return value
end
function rounded(actual::T, expected) where {T}
    if isinf(actual)
        largest = exact(floatmax(T))
        return actual > 0 && expected >= largest + (largest - exact(prevfloat(floatmax(T)))) / 2
    end
    isfinite(actual) && actual >= 0 || return false
    center = exact(actual)
    previous, following = prevfloat(actual), nextfloat(actual)
    low = (exact(previous) + center) / 2
    high =
        isfinite(following) ? (exact(following) + center) / 2 :
        center + (center - exact(previous)) / 2
    even = T === Float32 ? iseven(reinterpret(UInt32, actual)) : iseven(reinterpret(UInt64, actual))
    return low < expected < high || even && (expected == low || expected == high)
end
end

function _em_energy_old(es::EMPIC1D{T}) where {T}
    c2, total = es.c * es.c, zero(T)
    for i in eachindex(es.Ex)
        total += es.Ex[i]^2 + es.Ey[i]^2 + c2 * es.Bz[i]^2
    end
    return T(0.5) * total * es.g.dx[1]
end
function _em_energy_old(es::EMPIC{D,T}) where {D,T}
    c2, total = es.c * es.c, zero(T)
    for c = 1:3, I in eachindex(es.E[c])
        total += es.E[c][I]^2 + c2 * es.B[c][I]^2
    end
    return T(0.5) * total * prod(es.g.dx)
end
function _em_energy_allocation(es)
    em_field_energy(es)
    return @allocated em_field_energy(es)
end

@testset "EM field energy retains complete finite integrals" begin
    oracle = EMFieldEnergyOracle
    for T in (Float32, Float64), D = 1:3, legacy in (false, true)
        legacy && D != 1 && continue
        exponent = T === Float32 ? 30 : 200
        for large in (false, true)
            c = T(10)^(large ? exponent : -exponent)
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
            es = legacy ? EMPIC1D(g, 0; c) : EMPIC(g, 0; c)
            E, B = oracle.components(es)
            @test em_field_energy(es) === zero(T)
            for a in B
                fill!(a, -inv(c))
                @test oracle.rounded(em_field_energy(es), oracle.integral(es))
                fill!(a, zero(T))
            end
            for a in E
                fill!(a, one(T))
                @test oracle.rounded(em_field_energy(es), oracle.integral(es))
                fill!(a, zero(T))
            end
        end
        for large in (false, true)
            exponent = T === Float32 ? 25 : 180
            amplitude = T(10)^(large ? -exponent : exponent)
            length_exponent = D == 3 ? exponent ÷ 2 : exponent
            length_scale = T(10)^(large ? length_exponent : -length_exponent)
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> length_scale, D))
            es = legacy ? EMPIC1D(g, 0) : EMPIC(g, 0)
            E, B = oracle.components(es)
            for a in (E..., B...)
                fill!(a, -amplitude)
                expected = oracle.integral(es)
                @test expected > 0
                @test oracle.rounded(em_field_energy(es), expected)
                @test isfinite(em_field_energy(es))
                fill!(a, zero(T))
            end
        end
    end
end

@testset "EM field energy rounding, volume associations and nonfinite inputs" begin
    oracle = EMFieldEnergyOracle
    for T in (Float32, Float64)
        # A tiny square beyond a midpoint still decides final rounding. A
        # single global normalization that erases this tail chooses the wrong
        # neighbor, even though the complete energy is ordinary-sized.
        for legacy in (false, true)
            g = FourierGrid((4,), (T(4),))
            es = legacy ? EMPIC1D(g, 0) : EMPIC(g, 0)
            E, _ = oracle.components(es)
            E[1][1] = one(T)
            if T === Float32
                E[1][2] = ldexp(one(T), -12)
            else
                E[1][2:3] .= ldexp(one(T), -27)
            end
            E[1][4] = nextfloat(zero(T))
            @test oracle.rounded(em_field_energy(es), oracle.integral(es))
            @test em_field_energy(es) === nextfloat(T(0.5))
        end
        # Every cell square fits, but their unscaled sum does not.
        g = FourierGrid((4, 4), (one(T), one(T)))
        es = EMPIC(g, 0)
        fill!(es.E[1], sqrt(floatmax(T) / 8))
        @test isfinite(em_field_energy(es))
        @test oracle.rounded(em_field_energy(es), oracle.integral(es))
        # Exact half of the smallest subnormal rounds to even zero, while the
        # next count crosses that tie. Here the field's square itself underflows.
        exponent = T === Float32 ? -75 : -538
        amplitude = ldexp(one(T), exponent)
        g = FourierGrid((4,), (T(8),))
        es = EMPIC(g, 0)
        for count = 1:4
            es.E[1][1:count] .= amplitude
            @test oracle.rounded(em_field_energy(es), oracle.integral(es))
        end
        # Small spacings can lose bits before a later large spacing restores
        # them; exercise all placements independently of prod's association.
        small = T === Float32 ? T(1e-30) : T(1e-200)
        for lengths in
            ((small, small, inv(small)), (small, inv(small), small), (inv(small), small, small))
            g = FourierGrid((4, 4, 4), lengths)
            es = EMPIC(g, 0)
            fill!(es.E[2], one(T))
            @test oracle.rounded(em_field_energy(es), oracle.integral(es))
            @test em_field_energy(es) > 0
        end
        for legacy in (false, true)
            g = FourierGrid((4,), (T(4),))
            es = legacy ? EMPIC1D(g, 0) : EMPIC(g, 0)
            E, B = oracle.components(es)
            for a in (E..., B...)
                a[1] = floatmax(T)
                @test em_field_energy(es) == T(Inf)
                a[1] = T(Inf)
                @test isequal(em_field_energy(es), _em_energy_old(es))
                a[1] = T(NaN)
                @test isnan(em_field_energy(es))
                a[1] = zero(T)
            end
        end
    end
end

@testset "EM field energy ordinary state and allocation compatibility" begin
    rng = MersenneTwister(120415)
    for T in (Float32, Float64), D = 1:3, legacy in (false, true)
        legacy && D != 1 && continue
        g = FourierGrid(ntuple(_ -> 6, D), ntuple(d -> T(d + 1), D))
        es = legacy ? EMPIC1D(g, 0; c = T(1.7)) : EMPIC(g, 0; c = T(1.7))
        E, B = EMFieldEnergyOracle.components(es)
        for a in (E..., B...)
            randn!(rng, a)
        end
        before = map(copy, (E..., B...))
        @test em_field_energy(es) === _em_energy_old(es)
        @test _em_energy_allocation(es) == 0
        @test all(isequal.(before, (E..., B...)))
    end
end
