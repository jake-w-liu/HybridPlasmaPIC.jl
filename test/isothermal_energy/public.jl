module PublicChecks
using HybridPlasmaPIC, Test
import ..Oracle
using ..Oracle: Q, oracle_integral
const H = HybridPlasmaPIC
function original(n::Array{T}, closure, g) where {T}
    pe = similar(n)
    electron_pressure!(pe, n, closure)
    gamma = T(H.closure_gamma(closure))
    if gamma == one(T)
        total = zero(T)
        for i in eachindex(n, pe)
            n[i] > zero(T) && (total += pe[i] * log(n[i]))
        end
        return total * prod(g.dx)
    end
    return sum(pe) / (gamma - one(T)) * prod(g.dx)
end
function public_measure(n, closure, g)
    electron_internal_energy(n, closure, g)
    return @timed electron_internal_energy(n, closure, g)
end
@testset "Public isothermal owners, dimensions and energy budgets" begin
    for T in (Float32, Float64), C in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 2, D), ntuple(_ -> one(T), D))
        B = ntuple(_ -> zeros(T, g.n), 3)
        for density in (T(0.5), T(1.5), one(T) + ldexp(one(T), -10)), te in (C(0), C(0.7))
            n = fill(density, g.n)
            closure = IsothermalElectrons(te)
            lo, hi = oracle_integral(n, te, g.dx)
            before = deepcopy((n, B))
            value = electron_internal_energy(n, closure, g)
            ref = H._from_exact(T, (lo + hi) / 2)
            if density == one(T) + ldexp(one(T), -10) || T === Float32 || iszero(te)
                @test Oracle.rounded(value, lo, hi)
            else
                radius = 2max(abs(Q(nextfloat(ref)) - Q(ref)), abs(Q(ref) - Q(prevfloat(ref))))
                @test lo - radius <= Q(value) <= hi + radius
            end
            @test isequal(electron_internal_energy(n, B, closure, g), value)
            ions = ParticleSet{D,T}(0)
            for particles in (ions, [ions], view([ions], :))
                budget = energy_budget(particles, B, n, closure, g)
                @test isequal(budget.electron_internal, value)
                @test isequal(budget.total, value)
            end
            @test isequal((n, B), before)
            @test public_measure(n, closure, g).bytes == 0
        end
        n = ones(T, g.n)
        closure = IsothermalElectrons(C(0.7))
        @test_throws DimensionMismatch electron_internal_energy(
            zeros(T, ntuple(_ -> 3, D)),
            closure,
            g,
        )
        @test_throws ArgumentError electron_internal_energy(n, closure, g; pe = ones(T, g.n))
        wrongB = ntuple(_ -> zeros(T, ntuple(_ -> 3, D)), 3)
        @test_throws DimensionMismatch electron_internal_energy(n, wrongB, closure, g)
        @test isnan(electron_internal_energy(n, CGLElectrons(T(1), T(1), T(1), T(1)), g))
        for bad in (T(NaN), T(Inf), -T(Inf))
            altered = copy(n)
            altered[1] = bad
            @test isequal(
                electron_internal_energy(altered, closure, g),
                original(altered, closure, g),
            )
        end
    end
end
@testset "Unsupported coefficients and gamma-rounded polytropes keep their owner" begin
    for T in (Float32, Float64), C in (Float16, BigFloat)
        g = FourierGrid((4,), (T(4),))
        n = T[0.5, 0.75, 1.5, 2]
        closure = IsothermalElectrons(C(0.7))
        @test isequal(electron_internal_energy(n, closure, g), original(n, closure, g))
    end
    g = FourierGrid((4,), (4.0f0,))
    n = Float32[0.5, 0.75, 1.5, 2]
    closure = PolytropicElectrons(0.7, 1.0, nextfloat(1.0))
    @test Float32(H.closure_gamma(closure)) == 1.0f0
    @test isequal(electron_internal_energy(n, closure, g), original(n, closure, g))
end
@testset "Public large weak-wave independent eighth-moment accuracy" begin
    for T in (Float32, Float64), count in (256, 4096, 4194304)
        stored32 =
            [1.0f0 + 0.001f0 * sin(Float32(2pi) * Float32(i) / Float32(count)) for i = 1:count]
        n = T.(stored32)
        g = FourierGrid((count,), (T(count),))
        closure = IsothermalElectrons(T(0.7))
        before = copy(n)
        lo, hi = Oracle.moment_eight_bounds(stored32, closure.Te, g.dx)
        result = public_measure(n, closure, g)
        @test Oracle.rounded(result.value, lo, hi)
        @test result.bytes == 0
        @test isequal(n, before)
        println((
            T,
            count,
            case = :weak_wave,
            value = result.value,
            bytes = result.bytes,
            seconds = result.time,
        ))
        flush(stdout)
    end
end
@testset "Public uniform states avoid large-grid reduction error" begin
    for T in (Float32, Float64), count in (256, 4096, 4194304), density in (T(0.5), T(0.75), T(1.5))
        n = fill(density, count)
        g = FourierGrid((count,), (T(count),))
        closure = IsothermalElectrons(T(0.7))
        lo, hi = Oracle.oracle_log(Q(density); bits = 256)
        factor = Q(count) * Q(closure.Te) * Q(density) * prod(Q, g.dx)
        lo *= factor
        hi *= factor
        result = public_measure(n, closure, g)
        if T === Float32
            @test Oracle.rounded(result.value, lo, hi)
        else
            ref = H._from_exact(T, (lo + hi) / 2)
            radius = 2max(abs(Q(nextfloat(ref)) - Q(ref)), abs(Q(ref) - Q(prevfloat(ref))))
            @test lo - radius <= Q(result.value) <= hi + radius
        end
        @test result.bytes == 0
        @test first(n) == last(n) == density
        println((
            T,
            count,
            density,
            value = result.value,
            bytes = result.bytes,
            seconds = result.time,
        ))
        flush(stdout)
    end
end
end
