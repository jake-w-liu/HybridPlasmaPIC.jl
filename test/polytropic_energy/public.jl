module PublicChecks
using Test, HybridPlasmaPIC
using ..GammaOneChecks: complete_oracle
using ..GammaOneChecks.Oracle: rounded
const H = HybridPlasmaPIC
function legacy(n, c, g)
    gamma = eltype(n)(H.closure_gamma(c))
    pe = similar(n)
    electron_pressure!(pe, n, c)
    if gamma == one(eltype(n))
        total = zero(eltype(n))
        for i in eachindex(pe, n)
            n[i] > 0 && (total += pe[i] * log(n[i]))
        end
        return total * prod(g.dx)
    end
    return sum(pe) / (gamma - one(eltype(n))) * prod(g.dx)
end
function outcome(f)
    try
        return (:value, f())
    catch error
        return (:error, typeof(error), sprint(showerror, error))
    end
end
struct CustomClosure{T} <: ElectronClosure
    coefficient::T
    gamma::T
end
H.closure_gamma(c::CustomClosure) = c.gamma
function H.electron_pressure!(pe, n, c::CustomClosure)
    pe .= c.coefficient .* max.(n, zero(eltype(n)))
    return pe
end
@testset "Public dimensions, complete range and stored gamma" begin
    for T in (Float32, Float64), D = 1:3
        k = T === Float32 ? 100 : 600
        for sign in (-1, 1)
            spacing = ldexp(one(T), -sign * k)
            extents = ntuple(i -> i == 1 ? 2spacing : T(2), D)
            g = FourierGrid(ntuple(_ -> 2, D), extents)
            n = fill(ldexp(one(T), sign * k), g.n)
            before = copy(n)
            c = PolytropicElectrons(one(T), one(T), T(2))
            @test electron_internal_energy(n, c, g) === ldexp(T(2^D), sign * k)
            @test isequal(n, before)
        end
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> one(T), D))
        n = ones(T, g.n)
        for C in (Float32, Float64), gamma in (prevfloat(one(C)), nextfloat(one(C)))
            c = PolytropicElectrons(C(0.5), one(C), gamma)
            @test electron_internal_energy(n, c, g) === T(0.5 / (Float64(gamma) - 1))
        end
        for C in (Float32, Float64)
            n = fill(one(T) + eps(T), g.n)
            before = copy(n)
            c = PolytropicElectrons(C(0.7), C(3), one(C))
            lo, hi = complete_oracle(n, c, g.dx)
            @test rounded(electron_internal_energy(n, c, g), lo, hi)
            B = ntuple(_ -> zeros(T, g.n), 3)
            @test isequal(electron_internal_energy(n, B, c, g), electron_internal_energy(n, c, g))
            @test isequal(n, before)
            @test all(b -> all(iszero, b), B)
        end
    end
end
@testset "Public tiny-term midpoint direction" begin
    a = 1.0 + 2.0^-26
    b = 1.0 + 2.0^-27
    lower = 1.0 + 3 * 2.0^-27
    for k in (26, 27, 28)
        dx = (2.0^k - 1, 2.0^k + 1, b)
        g = FourierGrid((2, 2, 2), 2 .* dx)
        n = zeros(2, 2, 2)
        n[1] = 1
        n[2] = 0.5
        c = PolytropicElectrons(a, 1.0, 2.0^(2k))
        @test electron_internal_energy(n, c, g) === nextfloat(lower)
    end
end
@testset "Public validation, nonfinite propagation and generic closures" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(4),))
        valid = PolytropicElectrons(T(0.7), one(T), T(2))
        n = ones(T, g.n)
        @test_throws ArgumentError electron_internal_energy(n, valid, g; pe = copy(n))
        bad = ones(T, 3)
        expected = outcome(() -> H._require_grid_array(:n, bad, g))
        @test first(expected) === :error
        @test isequal(
            outcome(() -> electron_internal_energy(bad, valid, g; pe = copy(bad))),
            expected,
        )
        for gamma in T[0.5, 1, 2], value in T[NaN, Inf, -Inf]
            n = T[1, value, 0.5, -1]
            c = PolytropicElectrons(T(0.7), one(T), gamma)
            before = copy(n)
            @test isequal(
                outcome(() -> electron_internal_energy(n, c, g)),
                outcome(() -> legacy(n, c, g)),
            )
            @test isequal(n, before)
        end
        n = T[0.5, 1, 2, -1]
        for c in (
            PolytropicElectrons{T}(T(-1), one(T), T(2)),
            PolytropicElectrons{T}(one(T), zero(T), T(2)),
            PolytropicElectrons{T}(one(T), one(T), zero(T)),
            PolytropicElectrons{T}(T(Inf), one(T), one(T)),
            PolytropicElectrons{T}(one(T), one(T), T(NaN)),
        )
            @test isequal(
                outcome(() -> electron_internal_energy(n, c, g)),
                outcome(() -> legacy(n, c, g)),
            )
        end
        for gamma in T[1, 2]
            c = CustomClosure(T(0.7), gamma)
            @test isequal(electron_internal_energy(n, c, g), legacy(n, c, g))
        end
    end
end
end
