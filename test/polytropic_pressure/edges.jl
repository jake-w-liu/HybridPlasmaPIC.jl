module EdgeChecks
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
function captured(f)
    try
        return (:value, f())
    catch error
        # This test compares actual legacy failures with candidate failures.
        return (:error, typeof(error), sprint(showerror, error))
    end
end
println((VERSION,))
flush(stdout)
@testset "Nonfinite and directly constructed invalid closures retain legacy outcomes" begin
    for T in (Float32, Float64), C in (Float32, Float64)
        closures = [PolytropicElectrons{C}(x, one(C), C(1.5)) for x in (C(-1), C(Inf), C(NaN))]
        append!(
            closures,
            [PolytropicElectrons{C}(C(0.7), x, C(1.5)) for x in (C(-1), C(0), C(Inf), C(NaN))],
        )
        append!(
            closures,
            [PolytropicElectrons{C}(C(0.7), one(C), x) for x in (C(-1), C(0), C(Inf), C(NaN))],
        )
        for c in closures, density in (T(NaN), T(Inf), -T(Inf), zero(T), -one(T), T(0.5), T(2))
            n = T[density]
            pe = similar(n)
            old = captured(() -> first(H._legacy_polytropic_pressure!(pe, n, c)))
            new = captured(() -> public_candidate(density, c))
            @test isequal(old, new)
            @test isequal(n[1], density)
        end
    end
end
@testset "Normalized FMA residual retains the exact stored density ratio" begin
    rng = MersenneTwister(129705)
    checked = normal_outputs = 0
    for T in (Float32, Float64)
        U = T === Float32 ? UInt32 : UInt64
        signmask = typemax(U) >> 1
        values = T[
            nextfloat(zero(T)),
            floatmin(T),
            floatmax(T),
            prevfloat(one(T)),
            one(T),
            nextfloat(one(T)),
        ]
        for repetition = 1:1024
            value = reinterpret(T, rand(rng, U) & signmask)
            isfinite(value) && value > 0 && push!(values, value)
        end
        for (density, reference) in zip(values, reverse(values))
            value, base = Float64(density), Float64(reference)
            ratio = value / base
            H._polytropic_normal_positive(ratio) || continue
            a, ae = frexp(value)
            b, be = frexp(base)
            r = ldexp(ratio, be - ae)
            residual = fma(-r, b, a) / b
            relative = residual / r
            exact = Q(value) / Q(base) / Q(ratio) - 1
            @test abs(exact) <= Q(eps(Float64))
            @test abs(Q(relative) - exact) <= 4Q(eps(Float64)) * abs(exact) + Q(nextfloat(0.0))
            @test iszero(exact) ? iszero(relative) : sign(relative) == sign(exact)
            closure = PolytropicElectrons(T(0.7), reference, T(2))
            result = H._polytropic_pressure_probe(T, density, closure)
            if result !== nothing
                expected = Q(closure.pe0) * (Q(density) / Q(reference))^2
                rounded = HybridPlasmaPIC._from_exact(T, expected)
                @test isfinite(rounded) && !iszero(rounded)
                radius =
                    4max(
                        abs(Q(nextfloat(rounded)) - Q(rounded)),
                        abs(Q(rounded) - Q(prevfloat(rounded))),
                    )
                @test abs(Q(result) - expected) <= radius
                normal_outputs += 1
            end
            checked += 1
        end
    end
    @test checked > 1000
    @test normal_outputs > 100
    println((checked, normal_outputs))
    flush(stdout)
end
println("POLYTROPIC_ORDINARY_EDGE_CHECKS_COMPLETE")

end
