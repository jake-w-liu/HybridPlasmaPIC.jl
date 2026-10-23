using HybridPlasmaPIC, Test, Random

module HybridFieldEnergyRangeOracle
using HybridPlasmaPIC: magnetic_energy
exact(x) = Rational{BigInt}(x)
integral(B, dx) = sum(sum(x -> exact(x)^2, a) for a in B) * prod(exact, dx) / 2
function rounded(actual::T, expected) where {T}
    largest = exact(floatmax(T))
    if isinf(actual)
        return actual > 0 && expected >= largest + (largest - exact(prevfloat(floatmax(T)))) / 2
    end
    isfinite(actual) && actual >= 0 || return false
    center = exact(actual)
    previous, following = prevfloat(actual), nextfloat(actual)
    lower = (exact(previous) + center) / 2
    upper =
        isfinite(following) ? (exact(following) + center) / 2 :
        center + (center - exact(previous)) / 2
    even = T === Float32 ? iseven(reinterpret(UInt32, actual)) : iseven(reinterpret(UInt64, actual))
    return lower < expected < upper || even && (expected == lower || expected == upper)
end
function original(B::NTuple{3,<:Array{T}}, dx) where {T}
    total = zero(T)
    for component in B, value in component
        total += value^2
    end
    return T(0.5) * total * prod(dx)
end
function allocations(B, g)
    magnetic_energy(B, g)
    return @allocated magnetic_energy(B, g)
end
end

@testset "Hybrid magnetic complete integrals and spacing associations" begin
    oracle = HybridFieldEnergyRangeOracle
    for T in (Float32, Float64), D = 1:3
        power = T === Float32 ? 25 : 180
        for inverse in (false, true)
            amplitude = T(10)^(inverse ? -power : power)
            scale = T(10)^(inverse ? power ÷ D : -(power ÷ D))
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> scale, D))
            B = ntuple(_ -> zeros(T, g.n), 3)
            for component in B
                fill!(component, -amplitude)
                expected = oracle.integral(B, g.dx)
                before = map(copy, B)
                @test expected > 0
                @test oracle.rounded(magnetic_energy(B, g), expected)
                @test all(isequal.(B, before))
                fill!(component, zero(T))
            end
        end
    end
    for T in (Float32, Float64)
        power = T === Float32 ? 80 : 600
        small, large = ldexp(one(T), -power), ldexp(one(T), power)
        for lengths in (
            (small, small, large),
            (small, large, small),
            (large, small, small),
            (large, large, small),
            (large, small, large),
            (small, large, large),
        )
            g = FourierGrid((4, 4, 4), lengths)
            B = (ones(T, g.n), zeros(T, g.n), zeros(T, g.n))
            value = magnetic_energy(B, g)
            @test isfinite(value) && value > 0
            @test oracle.rounded(value, oracle.integral(B, g.dx))
        end
        g = FourierGrid((2, 2, 2), (large, large, large))
        B = ntuple(_ -> zeros(T, g.n), 3)
        @test magnetic_energy(B, g) === zero(T)
        @test oracle.rounded(magnetic_energy(B, g), oracle.integral(B, g.dx))
    end
end

@testset "Hybrid magnetic subnormal ties and decisive small tails" begin
    oracle = HybridFieldEnergyRangeOracle
    for T in (Float32, Float64)
        amplitude = ldexp(one(T), T === Float32 ? -75 : -538)
        g = FourierGrid((4,), (T(8),))
        B = ntuple(_ -> zeros(T, g.n), 3)
        for count = 1:4
            B[1][1:count] .= amplitude
            @test oracle.rounded(magnetic_energy(B, g), oracle.integral(B, g.dx))
        end
        g = FourierGrid((4,), (T(4),))
        B = ntuple(_ -> zeros(T, g.n), 3)
        B[1][1] = one(T)
        if T === Float32
            B[1][2] = ldexp(one(T), -12)
        else
            B[1][2:3] .= ldexp(one(T), -27)
        end
        B[1][4] = nextfloat(zero(T))
        @test oracle.rounded(magnetic_energy(B, g), oracle.integral(B, g.dx))
        @test magnetic_energy(B, g) === nextfloat(T(0.5))
    end
end

@testset "Field energy finite-looking output overflow boundaries" begin
    oracle = HybridFieldEnergyRangeOracle
    for T in (Float32, Float64)
        k = (exponent(floatmax(T)) - 1) ÷ 2
        tail = fld(2k - precision(T) - 2, 2)
        g = FourierGrid((32,), (T(256),))
        B = ntuple(_ -> zeros(T, g.n), 3)
        legacy = EMPIC1D(g, 0; c = one(T))
        general = EMPIC(g, 0; c = one(T))
        for count in (0, 5, 6, 11, 12, 31)
            fill!(B[1], zero(T))
            B[1][1] = prevfloat(ldexp(one(T), k))
            B[1][2:count+1] .= ldexp(one(T), tail)
            expected = oracle.integral(B, g.dx)
            @test oracle.rounded(magnetic_energy(B, g), expected)
            legacy.Ex .= B[1]
            general.E[1] .= B[1]
            @test oracle.rounded(em_field_energy(legacy), expected)
            @test oracle.rounded(em_field_energy(general), expected)
        end
        for es in (legacy, general)
            if es isa EMPIC1D
                es.Bz .= es.Ex
                fill!(es.Ex, zero(T))
            else
                es.B[1] .= es.E[1]
                fill!(es.E[1], zero(T))
            end
            @test oracle.rounded(em_field_energy(es), oracle.integral(B, g.dx))
        end
    end
end

@testset "Hybrid magnetic ordinary output, allocation and nonfinite compatibility" begin
    oracle = HybridFieldEnergyRangeOracle
    rng = MersenneTwister(932167)
    for T in (Float32, Float64), D = 1:3, count in (4, 8, 16)
        g = FourierGrid(ntuple(_ -> count, D), ntuple(d -> T(d + 1), D))
        B = ntuple(_ -> randn(rng, T, g.n), 3)
        before = map(copy, B)
        if T === Float32
            @test oracle.rounded(magnetic_energy(B, g), oracle.integral(B, g.dx))
        else
            @test magnetic_energy(B, g) === oracle.original(B, g.dx)
        end
        @test oracle.allocations(B, g) == 0
        @test all(isequal.(B, before))
    end
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        B = ntuple(_ -> zeros(T, g.n), 3)
        for component in B, value in (T(Inf), T(-Inf), T(NaN), floatmax(T))
            component[1] = value
            observed = magnetic_energy(B, g)
            @test isequal(observed, oracle.original(B, g.dx))
            @test isequal(component[1], value)
            component[1] = zero(T)
        end
        wrong = (zeros(T, ntuple(_ -> 3, D)), B[2], B[3])
        @test_throws DimensionMismatch magnetic_energy(wrong, g)
    end
end

@testset "Field overflow bound includes magnetic light-speed factors" begin
    oracle = HybridFieldEnergyRangeOracle
    for T in (Float32, Float64), legacy in (false, true), c in (T(0.5), T(1.7), T(2))
        k = (exponent(floatmax(T)) - 1) ÷ 2
        tail = fld(2k - precision(T) - 2, 2)
        g = FourierGrid((32,), (T(256),))
        es = legacy ? EMPIC1D(g, 0; c) : EMPIC(g, 0; c)
        B = legacy ? (es.Bz,) : es.B
        B[1] .= ldexp(one(T), tail) / c
        B[1][1] = prevfloat(ldexp(one(T), k)) / c
        expected = oracle.integral(B, g.dx) * oracle.exact(c)^2
        before = map(copy, B)
        @test oracle.rounded(em_field_energy(es), expected)
        @test all(isequal.(B, before))
    end
    # A true zero spacing is a private-owner compatibility control. Public
    # FourierGrid construction requires positive spacings.
    for T in (Float32, Float64)
        E = (ones(T, 1),)
        value = HybridPlasmaPIC._em_field_energy_finish(
            E,
            (),
            one(T),
            (zero(T),),
            one(T),
            one(T),
            T(Inf),
        )
        @test value === zero(T)
    end
end
