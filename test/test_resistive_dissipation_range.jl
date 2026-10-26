module ResistiveDissipationRangeTests
using HybridPlasmaPIC, Test, Random

const Q = Rational{BigInt}

function original(J, eta::T, g) where {T}
    total = zero(T)
    for I in eachindex(J...)
        total += J[1][I] * J[1][I] + J[2][I] * J[2][I] + J[3][I] * J[3][I]
    end
    return eta * total * prod(g.dx)
end

exact(J, eta, g) = Q(eta) * sum(sum(Q(x)^2 for x in a) for a in J) * prod(Q, g.dx)

# Check the nearest-even rounding cell with exact rational endpoints, including
# the half-minimum-subnormal and overflow ties. No production helper is used.
function rounded(value::T, reference) where {T}
    isfinite(value) || isinf(value) || return false
    largest = Q(floatmax(T))
    overflow = largest + (largest - Q(prevfloat(floatmax(T)))) / 2
    isinf(value) && return value > 0 && reference >= overflow
    value < 0 && return false
    iszero(value) && return 0 <= reference <= Q(nextfloat(zero(T))) / 2
    lower = (Q(prevfloat(value)) + Q(value)) / 2
    upper = isfinite(nextfloat(value)) ? (Q(value) + Q(nextfloat(value))) / 2 : overflow
    even = T === Float32 ? iseven(reinterpret(UInt32, value)) : iseven(reinterpret(UInt64, value))
    return lower < reference < upper || even && (lower == reference || upper == reference)
end

cost(J, eta, g) = @allocated resistive_dissipation(J, eta, g)

@testset "Ohmic integral retains finite current/coefficient/volume range" begin
    for T in (Float32, Float64), D = 1:3
        e = T === Float32 ? 80 : 600
        small, large = ldexp(one(T), -e), ldexp(one(T), e)
        lengths = (
            ntuple(_ -> one(T), D),
            ntuple(_ -> small, D),
            ntuple(_ -> large, D),
            ntuple(d -> isodd(d) ? small : large, D),
            ntuple(d -> isodd(d) ? large : small, D),
        )
        for L in lengths
            g = FourierGrid(ntuple(_ -> 1, D), L)
            for amplitude in (zero(T), nextfloat(zero(T)), small, one(T), large, floatmax(T)),
                eta in (zero(T), nextfloat(zero(T)), floatmin(T), small, one(T), large, floatmax(T))

                J = (fill(amplitude, g.n), fill(-amplitude, g.n), zeros(T, g.n))
                snapshot = map(copy, J)
                result = resistive_dissipation(J, eta, g)
                # All nonzero arithmetic here consists of powers of two except
                # floatmax. Its square loses range and requires exact completion.
                @test rounded(result, exact(J, eta, g))
                @test result isa T
                @test all(isequal.(J, snapshot))
            end
        end
        # D=3 spacing products can overflow/underflow and then return to a
        # representable volume; all six orderings exercise prod association.
        if D == 3
            for L in (
                (small, small, large),
                (small, large, small),
                (large, small, small),
                (large, large, small),
                (large, small, large),
                (small, large, large),
            )
                g = FourierGrid((1, 1, 1), L)
                J = (ones(T, g.n), zeros(T, g.n), zeros(T, g.n))
                @test rounded(resistive_dissipation(J, one(T), g), exact(J, one(T), g))
            end
        end
    end
end

@testset "Ohmic subnormal and overflow rounding neighbors" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(4),))
        amplitude = ldexp(one(T), T === Float32 ? -75 : -538)
        J = ntuple(_ -> zeros(T, 4), 3)
        for count = 1:4, eta in (prevfloat(one(T)), one(T), nextfloat(one(T)))
            fill!(J[1], zero(T))
            J[1][1:count] .= amplitude
            @test rounded(resistive_dissipation(J, eta, g), exact(J, eta, g))
        end
        k = (exponent(floatmax(T)) - 1) ÷ 2
        tail = fld(2k - precision(T) - 2, 2)
        g = FourierGrid((32,), (T(256),))
        J = ntuple(_ -> zeros(T, 32), 3)
        for count in (0, 5, 6, 11, 12, 31), eta in (prevfloat(one(T)), one(T), nextfloat(one(T)))
            fill!(J[1], zero(T))
            J[1][1] = prevfloat(ldexp(one(T), k))
            J[1][2:count+1] .= ldexp(one(T), tail)
            @test rounded(resistive_dissipation(J, eta, g), exact(J, eta, g))
        end
    end
end

@testset "Ordinary Ohmic integral preserves bytes allocation and ownership" begin
    rng = MersenneTwister(72913)
    for T in (Float32, Float64), D = 1:3, count in (1, 4, 16), eta in (zero(T), T(0.07), one(T))
        g = FourierGrid(ntuple(_ -> count, D), ntuple(d -> T(count * d) / 2, D))
        for J in (
            ntuple(_ -> randn(rng, T, g.n), 3),
            ntuple(_ -> view(randn(rng, T, ntuple(_ -> 2count, D)), ntuple(_ -> 1:count, D)...), 3),
        )
            snapshot = map(copy, J)
            @test isequal(resistive_dissipation(J, eta, g), original(J, eta, g))
            cost(J, eta, g)
            @test cost(J, eta, g) == 0
            @test all(isequal.(J, snapshot))
        end
    end
end

@testset "Ohmic validation and nonfinite current propagation remain explicit" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(4),))
        for value in (T(NaN), T(Inf), T(-Inf)), eta in (zero(T), one(T))
            J = (T[value, 0, 0, 0], zeros(T, 4), zeros(T, 4))
            snapshot = map(copy, J)
            @test isequal(resistive_dissipation(J, eta, g), original(J, eta, g))
            @test all(isequal.(J, snapshot))
        end
        J = ntuple(_ -> ones(T, 4), 3)
        snapshot = map(copy, J)
        for eta in (T(-1), T(NaN), T(Inf))
            @test_throws ArgumentError resistive_dissipation(J, eta, g)
        end
        @test_throws DimensionMismatch resistive_dissipation((J[1], zeros(T, 3), J[3]), one(T), g)
        @test all(isequal.(J, snapshot))
        # Heating and its diagnostic must agree when the unscaled square
        # overflows or underflows but their common physical source is normal.
        e = T === Float32 ? 80 : 600
        for sign in (-1, 1)
            amplitude, eta = ldexp(one(T), sign * e), ldexp(one(T), -sign * e)
            J = (fill(amplitude, g.n), zeros(T, 4), zeros(T, 4))
            pressure_rate = zeros(T, 4)
            HybridPlasmaPIC._add_ohmic_pressure_heating!(pressure_rate, J, T(2), eta)
            @test isfinite(sum(pressure_rate)) && sum(pressure_rate) > 0
            @test resistive_dissipation(J, eta, g) == sum(pressure_rate)
        end
    end
end
end
