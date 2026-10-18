module ShockFluxBoundsTests
using HybridPlasmaPIC, Test, Random


const H = HybridPlasmaPIC
const Q = Rational{BigInt}
println((VERSION,))
flush(stdout)
@testset "Float32 seven-operation and ULP error bound" begin
    rng = MersenneTwister(661732)
    u = Q(1, BigInt(1) << 24)
    @test Q(nextfloat((7 * 2.0^-24) / (1 - 14 * 2.0^-24))) >= 7u / (1 - 14u)
    for k = 1:4000
        x, y, z, w, c = ntuple(_ -> ldexp(Float32(randn(rng)), rand(rng, -25:25)), 5)
        c = abs(c)
        x2, y2, z2 = x * x, y * y, z * z
        speed = x2 + y2 + z2
        half = 0.5f0 * speed
        directed = half * x
        flux = directed * w
        actual = c * flux
        certified =
            H._kinetic_normal_product(x, x, x2) &&
            H._kinetic_normal_product(y, y, y2) &&
            H._kinetic_normal_product(z, z, z2) &&
            isfinite(speed) &&
            H._kinetic_normal_product(0.5f0, speed, half) &&
            H._kinetic_normal_product(half, x, directed) &&
            H._kinetic_normal_product(directed, w, flux) &&
            H._kinetic_normal_product(c, flux, actual)
        if certified
            exact = Q(c) * Q(w) * Q(x) * (Q(x)^2 + Q(y)^2 + Q(z)^2) / 2
            @test Q(H._shock_flux_float32_product_error(actual)) >= abs(Q(actual) - exact)
        end
        left, right = ntuple(_ -> reinterpret(Float32, rand(rng, UInt32)), 2)
        result = left + right
        if isfinite(left) && isfinite(right) && isfinite(result)
            @test Q(H._shock_flux_float32_sum_error(result)) >= abs(Q(result) - Q(left) - Q(right))
        end
    end
    for result in (
        0.0f0,
        -0.0f0,
        nextfloat(0.0f0),
        -nextfloat(0.0f0),
        floatmin(Float32),
        floatmax(Float32),
        -floatmax(Float32),
    )
        @test isfinite(H._shock_flux_float32_sum_error(result))
        @test H._shock_flux_float32_sum_error(result) > 0
    end
    @test H._shock_flux_float32_sum_error(Inf32) == Inf
    @test H._shock_flux_float32_sum_error(-Inf32) == Inf
end
println("FAST_BOUND_CHECKS_COMPLETE")
flush(stdout)

end
