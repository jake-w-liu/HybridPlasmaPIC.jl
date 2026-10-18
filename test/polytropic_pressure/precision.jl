module PrecisionChecks
using HybridPlasmaPIC, Test, Random, Statistics
import ..Oracle
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
function public_pressure(::Type{S}, density::T, closure) where {S,T}
    source = T[density]
    destination = Vector{S}(undef, 1)
    electron_pressure!(destination, source, closure)
    return only(destination)
end
public_candidate(density::T, closure) where {T} = public_pressure(T, density, closure)
function root_bounds(base::Q, p)
    scale = BigInt(1) << 192
    integer = isqrt(fld(numerator(base) * scale^2, denominator(base)))
    @test integer^2 * denominator(base) <=
          numerator(base) * scale^2 <
          (integer + 1)^2 * denominator(base)
    return (integer // scale)^p, ((integer + 1) // scale)^p
end
println((VERSION,))
flush(stdout)
@testset "Destination precision and complete stored inputs" begin
    for S in (Float32, Float64),
        T in (Float32, Float64),
        C in (Float32, Float64),
        n in T[0.5, 0.75, 1, 1.5, 2],
        gamma in C[0.5, 1, 1.5, 2, 3],
        reference in C[1, 1.25]

        c = PolytropicElectrons(C(0.7), reference, gamma)
        base = Q(n) / Q(c.n0)
        if isinteger(gamma)
            lo = hi = base^Int(gamma)
        else
            lo, hi = root_bounds(base, Int(2gamma))
        end
        lo *= Q(c.pe0)
        hi *= Q(c.pe0)
        rounded = HybridPlasmaPIC._from_exact(S, (lo + hi) / 2)
        radius =
            4max(abs(Q(nextfloat(rounded)) - Q(rounded)), abs(Q(rounded) - Q(prevfloat(rounded))))
        value = public_pressure(S, n, c)
        @test value isa S
        @test lo - radius <= Q(value) <= hi + radius
    end
    for S in (Float32, Float64)
        for (n, c, exact) in (
            (2.0f0^100, PolytropicElectrons(0.5f0, 1.0f0, 2.0f0), Q(2)^199),
            (2.0f0^-100, PolytropicElectrons(2.0f0^100, 1.0f0, 2.0f0), Q(2)^-100),
            (1.0f0, PolytropicElectrons(2.0^600, 1.0, 2.0), Q(2)^600),
            (1.0f0, PolytropicElectrons(nextfloat(0.0), 1.0, 2.0), Q(nextfloat(0.0))),
        )
            @test isequal(public_pressure(S, n, c), HybridPlasmaPIC._from_exact(S, exact))
        end
    end
end
@testset "Raw signed powers and closure flooring retain their separate contracts" begin
    for S in (Float32, Float64), T in (Float32, Float64), C in (Float32, Float64)
        for density in T[-2, -0.5], gamma in C[1, 2, 3]
            c = PolytropicElectrons(C(0.7), C(1.25), gamma)
            @test isequal(H._polytropic_pressure(S, density, c), S(c.pe0 * (density / c.n0)^c.γ))
            @test iszero(public_candidate(density, c))
        end
        c = PolytropicElectrons(C(0.7), C(1.25), C(1.5))
        @test_throws DomainError H._polytropic_pressure(S, -one(T), c)
        @test iszero(public_candidate(-one(T), c))
        for density in T[0, -0.0, 1.5], coefficient in C[0, -0.0, 0.7], gamma in C[1, 1.5, 2, 3]
            # Positive ordinary values are covered by the independent bounds above.
            (iszero(density) || iszero(coefficient)) || continue
            c = PolytropicElectrons(coefficient, C(1.25), gamma)
            @test isequal(H._polytropic_pressure(S, density, c), S(c.pe0 * (density / c.n0)^c.γ))
            if iszero(density)
                @test isequal(
                    public_candidate(density, c),
                    T(c.pe0 * (max(density, zero(T)) / c.n0)^c.γ),
                )
            end
        end
    end
end
println("POLYTROPIC_SCALAR_OWNERSHIP_COMPLETE")

end
