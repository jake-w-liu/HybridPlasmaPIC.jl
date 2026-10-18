module GeneralChecks
using HybridPlasmaPIC, Test, Random
import ..Oracle
using ..Oracle: Q, oracle_log, oracle_integral
const H = HybridPlasmaPIC
const O = Oracle
function accurate(value, lo, hi, ::Type{T}) where {T}
    isfinite(value) || return isequal(value, HybridPlasmaPIC._from_exact(T, (lo + hi) / 2))
    ref = HybridPlasmaPIC._from_exact(T, (lo + hi) / 2)
    iszero(ref) && return iszero(value)
    # Two output ulps are a finite ordinary-path test contract, not a claim of
    # universal correct rounding. Cold-path cases below still require exact bits.
    radius = 2max(abs(Q(nextfloat(ref)) - Q(ref)), abs(Q(ref) - Q(prevfloat(ref))))
    return lo - radius <= Q(value) <= hi + radius
end
@testset "Independent finite general isothermal accuracy" begin
    for T in (Float32, Float64), C in (Float32, Float64), D = 1:3
        dx = ntuple(_ -> T(0.5), D)
        for n in (
            T[0.5],
            T[0.75],
            T[1.5],
            T[3],
            T[0.5, 0.75, 0, -1],
            T[1.5, 3, 8],
            T[0.75, 1.5],
            T[2, 0.5, 0.5, 0.5, 0.5],
        )
            for te in (C(-0.7), zero(C), C(0.7))
                before = copy(n)
                lo, hi = oracle_integral(n, te, dx)
                state = H._isothermal_general_probe(n, te, dx)
                result = H._isothermal_energy(n, te, dx)
                @test accurate(result, lo, hi, T)
                @test isequal(n, before)
                @test isequal(H._isothermal_energy(view(n, :), te, dx), result)
                if state.certified
                    @test Q(state.lower) <= lo <= hi <= Q(state.upper)
                else
                    @test isequal(result, H._isothermal_exact_energy(n, te, dx))
                end
            end
        end
    end
end
@testset "Nonfinite legacy compatibility" begin
    for T in (Float32, Float64), C in (Float32, Float64)
        for bad in (T(NaN), T(Inf), -T(Inf))
            n = T[1.5, bad, 0.75]
            grid = FourierGrid((3,), (T(3),))
            closure = IsothermalElectrons(C(0.7))
            @test isequal(
                H._isothermal_energy(n, closure.Te, grid.dx),
                electron_internal_energy(n, closure, grid),
            )
        end
        for te in (C(NaN), C(Inf), -C(Inf)), n in (T[1], T[0, -1], T[0.75, 1.5])
            @test isequal(
                H._isothermal_energy(n, te, (T(1),)),
                H._isothermal_legacy_energy(n, te, (T(1),)),
            )
        end
    end
end
@testset "Range controls retain independent complete integrals" begin
    for T in (Float32, Float64), C in (Float32, Float64)
        a = ldexp(one(T), T === Float32 ? 90 : 700)
        b = inv(a)
        cases = Any[
            (T[1.5, 0.75], C(1), dx) for
            dx in ((a, a, b), (a, b, a), (b, a, a), (b, b, a), (b, a, b), (a, b, b))
        ]
        append!(
            cases,
            [
                (T[2], floatmax(C), (T(0.25),)),
                (T[0.5], nextfloat(zero(C)), (floatmax(T),)),
                (T[1], floatmax(C), (floatmax(T), floatmax(T))),
                (T[], C(1), (floatmax(T), floatmax(T))),
                (T[2], zero(C), (floatmax(T), floatmax(T))),
            ],
        )
        for multiplier = 1:4, density in T[0.75, 0.875, 1.125, 1.25]
            push!(
                cases,
                (
                    fill(density, 2),
                    C(T(multiplier) * floatmin(T)),
                    (ldexp(one(T), T === Float32 ? 79 : 799),),
                ),
            )
        end
        for (n, te, dx) in cases
            lo, hi = oracle_integral(n, te, dx)
            state = H._isothermal_general_probe(n, te, dx)
            value = H._isothermal_energy(n, te, dx)
            @test accurate(value, lo, hi, T)
            if !state.certified
                @test isequal(value, H._isothermal_exact_energy(n, te, dx))
            else
                @test Q(state.lower) <= lo <= hi <= Q(state.upper)
            end
        end
    end
end
println("GENERAL_ISOTHERMAL_FOCUSED_COMPLETE")

end
