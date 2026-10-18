module Float32Checks
using HybridPlasmaPIC, Test, Random, Statistics
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
println((VERSION,))
flush(stdout)
@testset "Power-precision correction bounds include Float64 residual arithmetic" begin
    u = Q(eps(Float64)) / 2
    for T in (Float32, Float64)
        epsilon = Q(eps(T))
        limit = Q(Float64(sqrt(eps(T)) / 4))
        t = (limit + Q(nextfloat(0.0))) / (1 - u)^3
        logarithm_remainder = t * epsilon / (2(1 - epsilon))
        y = t + logarithm_remainder
        total = logarithm_remainder + y^2 / (2(1 - y))
        @test y < 1
        @test total < epsilon / 16
        for exponent in (1, 2, 3, 17), residual in (-epsilon, -epsilon / 3, epsilon / 3, epsilon)
            @test abs(exponent * residual) <= limit
            exact = (1 + residual)^exponent - 1 - exponent * residual
            @test abs(exact) <= total
        end
    end
end
@testset "Native Float32 residuals and pressures satisfy independent rational bounds" begin
    rng = MersenneTwister(971263)
    values = Float32[
        nextfloat(0.0f0),
        floatmin(Float32),
        floatmax(Float32),
        prevfloat(1.0f0),
        1.0f0,
        nextfloat(1.0f0),
    ]
    for _ = 1:4096
        x = reinterpret(Float32, rand(rng, UInt32) & 0x7fffffff)
        isfinite(x) && x > 0 && push!(values, x)
    end
    residuals = outputs = 0
    for (density, reference) in zip(values, reverse(values))
        ratio = Float64(density / reference)
        isfinite(ratio) && ratio >= 2Float64(floatmin(Float32)) || continue
        a, ae = frexp(Float64(density))
        b, be = frexp(Float64(reference))
        normalized = ldexp(ratio, be - ae)
        residual = fma(-normalized, b, a) / b / normalized
        exact = Q(density) / Q(reference) / Q(ratio) - 1
        @test abs(exact) <= Q(eps(Float32))
        @test abs(Q(residual) - exact) <= 4Q(eps(Float64)) * abs(exact) + Q(nextfloat(0.0))
        @test iszero(exact) ? iszero(residual) : sign(residual) == sign(exact)
        residuals += 1
        c = PolytropicElectrons(0.7f0, reference, 2.0f0)
        result = H._polytropic_pressure_probe_power(Float32, density, c, Float32)
        if result !== nothing
            truth = Q(c.pe0) * (Q(density) / Q(reference))^2
            rounded = HybridPlasmaPIC._from_exact(Float32, truth)
            @test isfinite(rounded) && rounded > 0
            radius = 4max(Q(nextfloat(rounded)) - Q(rounded), Q(rounded) - Q(prevfloat(rounded)))
            @test abs(Q(result) - truth) <= radius
            outputs += 1
        end
    end
    @test residuals > 2000
    @test outputs > 1000
    println((residuals, outputs))
    flush(stdout)
end
@testset "Float32 range and correction rejection preserve the wider probe" begin
    for (n, c) in (
        (2.0f0^100, PolytropicElectrons(2.0f0^-100, 1.0f0, 2.0f0)),
        (2.0f0^-100, PolytropicElectrons(2.0f0^100, 1.0f0, 2.0f0)),
        (2.0f0, PolytropicElectrons(1.0f0, prevfloat(2.0f0), 2.0f0^16)),
    )
        @test H._polytropic_pressure_probe_power(Float32, n, c, Float32) === nothing
        wide = H._polytropic_pressure_probe_power(Float32, n, c, Float64)
        @test wide !== nothing
        @test isequal(H._polytropic_pressure_probe(Float32, n, c), wide)
    end
    c = PolytropicElectrons(0.7f0, 1.25f0, 1.5f0)
    @test isequal(
        H._polytropic_pressure_probe(Float64, 1.5f0, c),
        H._polytropic_pressure_probe_power(Float64, 1.5f0, c, Float64),
    )
end
println("POLYTROPIC_FLOAT32_PROBE_COMPLETE")

end
