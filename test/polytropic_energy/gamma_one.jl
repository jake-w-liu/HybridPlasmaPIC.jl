module GammaOneChecks
using HybridPlasmaPIC, Test

include("../isothermal_energy/oracles.jl")
const P = (
    energy = HybridPlasmaPIC._polytropic_log_energy,
    _isothermal_nearunit_energy = HybridPlasmaPIC._isothermal_nearunit_energy,
    _isothermal_statistics = HybridPlasmaPIC._isothermal_statistics,
)
const Q = Rational{BigInt}

function complete_oracle(n, closure, dx)
    factor = Q(closure.pe0) * prod(Q, dx) / Q(closure.n0)
    # For exact binary powers, combine the log(2) coefficient algebraically.
    # Independent per-cell enclosures cannot certify exact zero at finite width.
    if all(value <= 0 || first(frexp(value)) == eltype(n)(0.5) for value in n)
        coefficient =
            sum((Q(value) * exponent(value) for value in n if value > 0); init = zero(Q)) * factor
        iszero(coefficient) && return zero(Q), zero(Q)
        lo, hi = Oracle.oracle_log(Q(2); bits = 768)
        return coefficient > 0 ? (coefficient * lo, coefficient * hi) :
               (coefficient * hi, coefficient * lo)
    end
    lower = upper = zero(Q)
    for value in n
        value > 0 || continue
        lo, hi = Oracle.oracle_log(Q(value); bits = 768)
        lower += Q(value) * lo
        upper += Q(value) * hi
    end
    return factor >= 0 ? (factor * lower, factor * upper) : (factor * upper, factor * lower)
end

println((VERSION,))
flush(stdout)
@testset "Gamma-one coefficient remains an exact complete ratio" begin
    for T in (Float32, Float64), C in (Float32, Float64), D = 1:3
        for (coefficient, reference) in ((C(0.7), C(3)), (C(0.75), C(1.5)), (zero(C), C(3)))
            closure = PolytropicElectrons(coefficient, reference, one(C))
            dx = ntuple(_ -> T(0.5), D)
            for n in (T[0.5, 0.75, 1.5], T[2, 0.5, 0.5, 0.5, 0.5], T[-1, 0, 1], T[])
                before = copy(n)
                lo, hi = complete_oracle(n, closure, dx)
                value = P.energy(n, closure, dx)
                @test Oracle.rounded(value, lo, hi)
                @test isequal(P.energy(view(n, :), closure, dx), value)
                @test isequal(n, before)
            end
        end
    end
end
@testset "Coefficient range loss precedes normal complete integrals" begin
    for T in (Float32, Float64)
        shift = T === Float32 ? 120 : 1000
        for (coefficient, reference, spacing) in (
            (floatmax(T), floatmin(T), ldexp(one(T), -shift)),
            (nextfloat(zero(T)), floatmax(T), ldexp(one(T), shift)),
        )
            closure = PolytropicElectrons(coefficient, reference, one(T))
            @test !isfinite(coefficient / reference) || iszero(coefficient / reference)
            for n in (T[0.5], T[2], T[0.75, 1.5])
                dx = (spacing, spacing)
                lo, hi = complete_oracle(n, closure, dx)
                value = P.energy(n, closure, dx)
                @test isfinite(value) && !iszero(value)
                @test Oracle.rounded(value, lo, hi)
            end
        end
        @test_throws ArgumentError P.energy(T[Inf], PolytropicElectrons(T(1), T(1), T(1)), (T(1),))
        @test_throws ArgumentError P.energy(
            T[1],
            PolytropicElectrons(T(1), T(1), nextfloat(T(1))),
            (T(1),),
        )
        @test_throws ArgumentError P.energy(T[1], PolytropicElectrons{T}(T(1), T(0), T(1)), (T(1),))
    end
end
println("POLYTROPIC_GAMMA_ONE_COMPLETE_COEFFICIENT_CHECKS_COMPLETE")

@testset "Near-unit reference division precedes rounding" begin
    for T in (Float32, Float64),
        C in (Float32, Float64),
        D = 1:3,
        reference in C[1, 3, 0.7],
        sign in (-1, 1)

        e = eps(T)
        n = T[1+e, 1-e, 1+2e, 1-2e, 1+sign*e]
        c = PolytropicElectrons(C(0.7), reference, one(C))
        dx = ntuple(_ -> T(0.5), D)
        lo, hi = complete_oracle(n, c, dx)
        @test Oracle.rounded(P.energy(n, c, dx), lo, hi)
        state = P._isothermal_statistics(n)
        value = P._isothermal_nearunit_energy(n, c.pe0, dx, state, c.n0)
        @test value !== nothing
        @test value === nothing || Oracle.rounded(value, lo, hi)
    end
end

@testset "Repeated densities retain complete logarithmic coefficients" begin
    for T in (Float32, Float64),
        count in (0, 1, 3, 17, 64),
        temperature in (Q(0), Q(7, 30), Q(-3, 7)),
        profile in (:repeated, :alternating, :distinct)

        n =
            profile === :repeated ? fill(T(0.75), count) :
            profile === :alternating ? T[(0.5, 0.75, 1.5)[mod1(i, 3)] for i = 1:count] :
            T[0.75 + i / 256 for i = 1:count]
        dx = (T(0.5),)
        before = copy(n)
        lo = hi = zero(Q)
        for density in n
            a, b = Oracle.oracle_log(Q(density); bits = 768)
            lo += Q(density) * a
            hi += Q(density) * b
        end
        factor = temperature * prod(Q, dx)
        lo, hi = factor >= 0 ? (factor * lo, factor * hi) : (factor * hi, factor * lo)
        value = HybridPlasmaPIC._isothermal_exact_energy(n, temperature, dx)
        @test Oracle.rounded(value, lo, hi)
        @test isequal(n, before)
    end
end

end
