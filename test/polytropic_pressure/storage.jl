module AlgebraicStorageChecks
using Test, HybridPlasmaPIC
const H = HybridPlasmaPIC
const U = Union{Float32,Float64}
@testset "Union storage retains real broadcast behavior" begin
    for S in (Float32, Float64, U), T in (Float32, Float64, U), C in (Float32, Float64, U)
        (S === U || T === U) || continue
        c = PolytropicElectrons{C}(0.7, 1.25, 5 / 3)
        for view_storage in (false, true), count in (0, 1, 4)
            n = T[(0.75f0, 1.5)[mod1(i, 2)] for i = 1:count]
            a = S[0.0f0 for i = 1:count]
            b = copy(a)
            before = copy(n)
            source = view_storage ? view(n, :) : n
            dest = view_storage ? view(a, :) : a
            reference = view_storage ? view(b, :) : b
            H._legacy_polytropic_pressure!(reference, source, c)
            @test electron_pressure!(dest, source, c) === dest
            @test isequal(dest, reference)
            @test isequal(n, before)
        end
    end
    for T in (Union{}, U)
        n = T[]
        pe = Float64[]
        c = PolytropicElectrons(0.7, 1.0, 5 / 3)
        @test electron_pressure!(pe, n, c) === pe
        @test isempty(pe)
    end
end

@testset "Raw pressure broadcast preserves Array alias and signed contracts" begin
    for T in (Float32, Float64), gamma in (T(2), T(3))
        coefficient = T(0.5)
        original = T[-2, -1, -0.0, 0.5, 1.5, 2]
        source = reshape(original, 2, 3)
        destination = reshape(original, 2, 3)
        @test source isa Array && destination isa Array
        @test source !== destination && Base.mightalias(source, destination)
        expected = coefficient .* copy(source) .^ gamma
        @test H._polytropic_unfloored_pressure!(destination, source, coefficient, gamma) ===
              destination
        @test isequal(destination, expected)
        wide = zeros(T, 2, 3)
        column = T[-1, 0.5]
        @test H._polytropic_unfloored_pressure!(wide, column, coefficient, gamma) === wide
        @test isequal(wide, repeat(coefficient .* column .^ gamma, 1, 3))
        n = T[-2, -1, -0.0, 0.5, 1.5, 2]
        before = copy(n)
        @test H._polytropic_unfloored_pressure!(n, n, coefficient, gamma) === n
        @test isequal(n, coefficient .* before .^ gamma)
    end
    for T in (Float32, Float64)
        n = T[-1]
        out = zeros(T, 1, 2)
        @test_throws DomainError H._polytropic_unfloored_pressure!(out, n, T(0.5), T(1.5))
        @test_throws MethodError H._polytropic_unfloored_pressure!(
            view(out, :),
            view(n, :),
            T(0.5),
            T(1.5),
        )
        n = reshape(T[-Inf, Inf, NaN, -0.0], 1, 4)
        out = zeros(T, 2, 4)
        @test H._polytropic_unfloored_pressure!(out, n, T(0.5), T(2)) === out
        @test isequal(out, repeat(T(0.5) .* n .^ T(2), 2, 1))
        empty = T[]
        @test H._polytropic_unfloored_pressure!(empty, empty, T(0.5), T(2)) === empty
    end
end
end
