module NGPBoundaryTests
using HybridPlasmaPIC, Test, Random, LinearAlgebra
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
exact_node(s) = floor(Int, Q(s) + Q(1, 2))

@testset "NGP nearest-node rounding boundaries" begin
    rng = MersenneTwister(918211)
    for T in (Float32, Float64)
        points = T[zero(T), -zero(T), nextfloat(zero(T)), -nextfloat(zero(T))]
        for center = -32:32
            tie = T(center) + T(0.5)
            append!(points, (prevfloat(tie), tie, nextfloat(tie)))
        end
        # Cover transitions where fractional nodes cease to be representable.
        for exponent = 0:precision(T)
            x = ldexp(one(T), exponent)
            append!(points, (prevfloat(x), x, nextfloat(x), -prevfloat(x), -x, -nextfloat(x)))
        end
        for _ = 1:1024
            push!(points, ldexp(T(2rand(rng) - 1), rand(rng, -32:32)))
        end
        for s in points
            node, weights = H._stencil1d(NGP(), s)
            @test node == exact_node(s)
            @test weights === (one(T),)
        end
        for invalid in (T(Inf), T(-Inf), T(NaN), floatmax(T))
            @test_throws InexactError H._stencil1d(NGP(), invalid)
        end
        # A tie at an even cell distinguishes upward ties from nearest-even.
        @test exact_node(T(0.5)) != round(Int, T(0.5))
        @test exact_node(prevfloat(T(0.5))) != floor(Int, prevfloat(T(0.5)) + T(0.5))
    end
end

@testset "NGP public coupling shares nearest-node policy" begin
    for T in (Float32, Float64), D = 1:3
        cells = ntuple(d -> 4d, D)
        g = FourierGrid(cells, T.(cells)) # Unit cell volume makes moment oracles exact.
        ps = ParticleSet{D,T}(2; q = T(-2), m = T(2))
        ps.weight .= one(T)
        field = reshape(T.(1:prod(cells)), cells)
        vector_field = ntuple(c -> T(c) .* field, 3)
        out = zeros(T, 2)
        vector_out = ntuple(_ -> zeros(T, 2), 3)
        deposited = zeros(T, cells)
        expected = zeros(T, cells)
        moment = ntuple(_ -> zeros(T, cells), 3)
        pressure = ntuple(_ -> zeros(T, cells), 6)
        values = T[1, 2]
        for center = 0:3, neighbor in (-1, 0, 1)
            tie = T(center) + T(0.5)
            s = neighbor < 0 ? prevfloat(tie) : neighbor > 0 ? nextfloat(tie) : tie
            for x in ps.x
                fill!(x, s)
            end
            for c = 1:3
                fill!(ps.v[c], T(c))
            end
            before = (map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(field), copy(values))
            index = CartesianIndex(ntuple(d -> mod(exact_node(s), cells[d]) + 1, D))
            fill!(expected, zero(T))
            expected[index] = T(3)
            deposit_scalar!(deposited, ps, values, g, NGP())
            @test deposited == expected
            @test sum(deposited) == sum(values)
            deposit_scalar_threaded!(deposited, ps, values, g, NGP())
            @test deposited == expected
            gather_scalar!(out, field, ps, g, NGP())
            @test out == fill(field[index], 2)
            @test dot(vec(field), vec(deposited)) == dot(values, out)
            gather_vector!(vector_out, vector_field, ps, g, NGP())
            for c = 1:3
                @test vector_out[c] == fill(T(c) * field[index], 2)
            end
            density!(deposited, ps, g, NGP())
            expected[index] = T(2)
            @test deposited == expected
            momentum!(moment, ps, g, NGP())
            for c = 1:3
                @test moment[c] == T(c) .* expected
            end
            current!(moment, ps, g, NGP())
            for c = 1:3
                @test moment[c] == ps.q * T(c) .* expected
            end
            @test isequal(before, (ps.x, ps.v, ps.weight, field, values))
            for c = 1:3
                ps.v[c][2] = -T(c)
            end
            pressure_tensor!(pressure, ps, g, NGP())
            for (component, (i, j)) in enumerate(((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)))
                @test pressure[component] == T(2i * j) .* expected
            end
        end
    end
end

@testset "NGP half-cell translation and validation" begin
    for T in (Float32, Float64), D = 1:3
        n = ntuple(_ -> 4, D)
        g = FourierGrid(n, ntuple(_ -> T(4), D))
        ps = ParticleSet{D,T}(1)
        for x in ps.x
            x[1] = T(0.5)
        end
        out = zeros(T, n)
        values = T[1]
        deposit_scalar!(out, ps, values, g, NGP())
        original = copy(out)
        for axis = 1:D
            ps.x[axis][1] += one(T)
            deposit_scalar!(out, ps, values, g, NGP())
            shift = ntuple(d -> d <= axis ? 1 : 0, D)
            @test out == circshift(original, shift)
        end
        ps.x[1][1] = T(NaN)
        before = copy(out)
        @test_throws ArgumentError deposit_scalar!(out, ps, values, g, NGP())
        @test isequal(out, before)
        gathered = T[7]
        @test_throws ArgumentError gather_scalar!(gathered, out, ps, g, NGP())
        @test gathered == T[7]
    end
end
end
