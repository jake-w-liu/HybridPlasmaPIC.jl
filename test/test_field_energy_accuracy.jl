module FieldEnergyAccuracyTests
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# The oracle converts each stored operand separately. It does not share the
# production bins, wide reduction, normalization, or exact-output conversion.
function squares(fields)
    result = zero(Q)
    for component in fields, value in component
        result += Q(value)^2
    end
    return result
end
integral(E, B, coefficient, dx) = (squares(E) + Q(coefficient) * squares(B)) * prod(Q, dx) / 2
components(es::EMPIC1D) = ((es.Ex, es.Ey), (es.Bz,))
components(es::EMPIC) = (es.E, es.B)
function rounded(value::Float32, expected::Q)
    isnan(value) && return false
    edge = Q(floatmax(Float32)) + (Q(floatmax(Float32)) - Q(prevfloat(floatmax(Float32)))) / 2
    isinf(value) && return value > 0 && expected >= edge
    signbit(value) && return false
    iszero(value) && return 0 <= expected <= Q(nextfloat(0.0f0)) / 2
    lower = (Q(prevfloat(value)) + Q(value)) / 2
    upper = isfinite(nextfloat(value)) ? (Q(value) + Q(nextfloat(value))) / 2 : edge
    return lower < expected < upper ||
           iseven(reinterpret(UInt32, value)) && (expected == lower || expected == upper)
end
magnetic_cost(B, g) = @allocated magnetic_energy(B, g)
field_cost(es) = @allocated em_field_energy(es)
inertia_cost(J, d, dx) = @allocated H._electron_inertia_energy(J, d, dx)

@testset "Exact Float32 square bins cover every stored exponent" begin
    rng = MersenneTwister(800217)
    for scale = 0:254,
        fraction in (UInt32(0), UInt32(1), UInt32(0x007fffff), rand(rng, UInt32) & 0x007fffff)

        bits = (UInt32(scale) << 23) | fraction
        value = reinterpret(Float32, bits)
        fields = (Float32[value, -value, -0.0],)
        before = copy(fields[1])
        @test H._field_square_sum_exact(fields) == squares(fields)
        @test isequal(fields[1], before)
        total, count = H._field_square_sum_wide(fields)
        bounds = H._field_roundoff_bounds(total, count)
        @test bounds !== nothing
        @test Q(bounds[1]) <= squares(fields) <= Q(bounds[2])
    end
    @test H._field_square_sum_exact(()) == 0
    @test H._field_square_sum_exact((Float32[],)) == 0
    largest = (big(2)^24 - 1)^2
    @test 3big(typemax(Int)) * largest < big(2)^113
    capacity = big(typemax(UInt128)) ÷ largest
    @test capacity * largest <= typemax(UInt128)
    @test (capacity + 1) * largest > typemax(UInt128)
    @test H._field_roundoff_bounds(1.0, Int128(1) << 51) === nothing
    for D = 1:3, trial = 1:32
        fields = ntuple(
            _ -> reshape(
                [reinterpret(Float32, rand(rng, UInt32) % UInt32(0x7f800000)) for _ = 1:4^D],
                ntuple(_ -> 4, D),
            ),
            3,
        )
        before = map(copy, fields)
        @test H._field_square_sum_exact(fields) == squares(fields)
        @test H._field_square_sum_exact(map(a -> view(a, ntuple(_ -> Colon(), D)...), fields)) ==
              squares(fields)
        total, count = H._field_square_sum_wide(fields)
        bounds = H._field_roundoff_bounds(total, count)
        @test Q(bounds[1]) <= squares(fields) <= Q(bounds[2])
        @test all(isequal.(fields, before))
    end
end

@testset "Complete Float32 field owners round independent exact integrals" begin
    rng = MersenneTwister(820317)
    for D = 1:3, trial = 1:32
        lengths = ntuple(_ -> ldexp(Float32(rand(rng, 4:8)) / 4, rand(rng, -90:90)), D)
        g = FourierGrid(ntuple(_ -> 4, D), lengths)
        c = ldexp(Float32(rand(rng, 4:8)) / 4, rand(rng, -100:100))
        d = ldexp(Float32(rand(rng, 4:8)) / 4, rand(rng, -120:120))
        B = ntuple(
            _ -> [
                ldexp(Float32(rand(rng, -8:8)) / 8, rand(rng, -120:120)) for
                _ in CartesianIndices(g.n)
            ],
            3,
        )
        before = map(copy, B)
        @test rounded(magnetic_energy(B, g), integral((), B, 1, g.dx))
        @test rounded(H._electron_inertia_energy(B, d, g.dx), integral((), B, d, g.dx))
        @test all(isequal.(B, before))
        for legacy in (false, true)
            legacy && D != 1 && continue
            es = legacy ? EMPIC1D(g, 0; c) : EMPIC(g, 0; c)
            E, M = components(es)
            for (index, a) in enumerate((E..., M...))
                a .= B[mod1(index, 3)]
            end
            saved = map(copy, (E..., M...))
            @test rounded(em_field_energy(es), integral(E, M, Q(c)^2, g.dx))
            @test all(isequal.((E..., M...), saved))
        end
    end
end

@testset "Float32 midpoint tails, independent coefficients and volume factors" begin
    for tail in (0.0f0, nextfloat(0.0f0)), field_sign in (-1, 1)
        g = FourierGrid((4,), (4.0f0,))
        B = (field_sign .* Float32[1, 2^-12, tail, 0], zeros(Float32, 4), zeros(Float32, 4))
        expected = integral((), B, 1, g.dx)
        @test rounded(magnetic_energy(B, g), expected)
        @test magnetic_energy(B, g) === (iszero(tail) ? 0.5f0 : nextfloat(0.5f0))
        # A Float64-only reduction erases the tiny positive square beyond the tie.
        iszero(tail) ||
            @test !rounded(Float32(sum(sum(abs2, Float64.(a)) for a in B) / 2), expected)
        for d in (prevfloat(1.0f0), 1.0f0, nextfloat(1.0f0), 1.7f0)
            @test rounded(H._electron_inertia_energy(B, d, g.dx), integral((), B, d, g.dx))
        end
        for c in (prevfloat(1.0f0), 1.0f0, nextfloat(1.0f0), 1.7f0), legacy in (false, true)
            es = legacy ? EMPIC1D(g, 0; c) : EMPIC(g, 0; c)
            E, M = components(es)
            E[1][4] = 2.0f0^-13
            M[1] .= B[1]
            @test rounded(em_field_energy(es), integral(E, M, Q(c)^2, g.dx))
        end
    end
    for dx in (
        (1.7f0,),
        (1.7f0, 1.3f0),
        (1.7f0, 1.3f0, 0.9f0),
        (2.0f0^-100, 2.0f0^-100, 2.0f0^100),
        (2.0f0^100, 2.0f0^100, 2.0f0^-100),
    )
        J = (Float32[1, 1.25, 2^-12, nextfloat(0.0f0)], zeros(Float32, 4), zeros(Float32, 4))
        for coefficient in (nextfloat(0.0f0), 1.7f0, floatmax(Float32))
            @test rounded(
                H._electron_inertia_energy(J, coefficient, dx),
                integral((), J, coefficient, dx),
            )
        end
    end
end

@testset "Ordinary Float32 component allocations and signed zero" begin
    rng = MersenneTwister(830417)
    for D = 1:3
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(i -> Float32(i + 1), D))
        B = ntuple(_ -> randn(rng, Float32, g.n), 3)
        magnetic_cost(B, g)
        inertia_cost(B, 1.7f0, g.dx)
        @test magnetic_cost(B, g) == 0
        @test inertia_cost(B, 1.7f0, g.dx) == 0
        for legacy in (false, true)
            legacy && D != 1 && continue
            es = legacy ? EMPIC1D(g, 0; c = 1.7f0) : EMPIC(g, 0; c = 1.7f0)
            E, M = components(es)
            foreach(a -> randn!(rng, a), (E..., M...))
            field_cost(es)
            @test field_cost(es) == 0
            for a in (E..., M...)
                fill!(a, -0.0f0)
            end
            @test em_field_energy(es) === 0.0f0
        end
        foreach(a -> fill!(a, -0.0f0), B)
        @test magnetic_energy(B, g) === 0.0f0
        @test H._electron_inertia_energy(B, 1.7f0, g.dx) === 0.0f0
        @test H._electron_inertia_energy(B, -0.0f0, g.dx) === -0.0f0
    end
    # The public budget converts coefficients to the field type. Its private
    # component owner previously rejected mismatched stored types; the new
    # Float32 certificate must not silently accept Float64 field products.
    for S in (Float32, Float64)
        C = S === Float32 ? Float64 : Float32
        J = (S[1, 1.1], zeros(S, 2), zeros(S, 2))
        @test_throws MethodError H._electron_inertia_energy(J, C(0.7), (S(1.3),))
    end
end

@testset "Large Float32 uniform alternating and exact midpoint fields" begin
    for side in (32, 80, 128, 160)
        g = FourierGrid((side, side, side), (Float32(side), Float32(side), Float32(side)))
        B = ntuple(_ -> zeros(Float32, g.n), 3)
        count = length(B[1])
        for alternating in (false, true)
            fill!(B[3], 1.25f0)
            alternating && (B[3][1:2:end] .= 1.0f0)
            expected =
                alternating ? (count ÷ 2) * (Q(1) + Q(1.25f0)^2) / 2 : count * Q(1.25f0)^2 / 2
            @test rounded(magnetic_energy(B, g), expected)
            magnetic_cost(B, g)
            @test magnetic_cost(B, g) == 0
            @test all(iszero, B[1]) && all(iszero, B[2])
            @test all(i -> B[3][i] === (alternating && isodd(i) ? 1.0f0 : 1.25f0), eachindex(B[3]))
        end
        fill!(B[3], 1.0f0)
        base = Float32(count ÷ 2)
        e = exponent(base)
        # Add exactly half a binary32 ulp using one or two exact power-of-two
        # squares. This exercises bounded fallback storage on every large grid.
        tails = iseven(e - 23) ? 1 : 2
        tail = ldexp(1.0f0, fld(e - 23, 2))
        B[1][1:tails] .= tail
        expected = Q(count) / 2 + tails * Q(tail)^2 / 2
        @test expected == (Q(base) + Q(nextfloat(base))) / 2
        @test rounded(magnetic_energy(B, g), expected)
        magnetic_cost(B, g)
        @test magnetic_cost(B, g) <= 256 * 1024
        B[2][1] = nextfloat(0.0f0)
        @test magnetic_energy(B, g) === nextfloat(base)
        @test rounded(magnetic_energy(B, g), expected + Q(B[2][1])^2 / 2)
    end
end

@testset "Shared Float32 exact completion has bounded owner allocations" begin
    for legacy in (false, true)
        g =
            legacy ? FourierGrid((4096,), (4096.0f0,)) :
            FourierGrid((64, 64, 64), (64.0f0, 64.0f0, 64.0f0))
        es = legacy ? EMPIC1D(g, 0) : EMPIC(g, 0)
        E, B = components(es)
        fill!(E[1], 1.0f0)
        base = Float32(length(E[1]) ÷ 2)
        e = exponent(base)
        tails = iseven(e - 23) ? 1 : 2
        tail = ldexp(1.0f0, fld(e - 23, 2))
        E[2][1:tails] .= tail
        expected = Q(length(E[1])) / 2 + tails * Q(tail)^2 / 2
        before = map(copy, (E..., B...))
        @test expected == (Q(base) + Q(nextfloat(base))) / 2
        @test rounded(em_field_energy(es), expected)
        field_cost(es)
        @test field_cost(es) <= 256 * 1024
        @test all(isequal.((E..., B...), before))
        J = (E[1], E[2], first(B))
        @test rounded(H._electron_inertia_energy(J, 1.0f0, g.dx), expected)
        inertia_cost(J, 1.0f0, g.dx)
        @test inertia_cost(J, 1.0f0, g.dx) <= 256 * 1024
        @test all(isequal.((E..., B...), before))
    end
end

@testset "Concurrent Float32 ordinary and exact component evaluations" begin
    for tail in (0.0f0, nextfloat(0.0f0), 0.125f0)
        g = FourierGrid((4,), (4.0f0,))
        B = (Float32[1, 2^-12, tail, 0], zeros(Float32, 4), zeros(Float32, 4))
        expected = integral((), B, 1, g.dx)
        before = map(copy, B)
        results = Vector{Float32}(undef, 64)
        Threads.@threads for i in eachindex(results)
            results[i] = magnetic_energy(B, g)
        end
        @test all(v -> rounded(v, expected), results)
        @test all(isequal.(B, before))
    end
end
end
