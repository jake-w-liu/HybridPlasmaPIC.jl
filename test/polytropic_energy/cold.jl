module ColdChecks
using Test, HybridPlasmaPIC

const P = (complete = HybridPlasmaPIC._polytropic_exact_energy,)
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
println((VERSION, root = pkgdir(H)))
flush(stdout)
@testset "Complete cold integrals against independent rational powers" begin
    for T in (Float32, Float64),
        C in (Float32, Float64),
        D = 1:3,
        (p, q) in ((1, 2), (3, 2), (5, 4), (7, 4), (2, 1), (3, 1), (7, 1))

        n = reshape(T[a^q for a in (1, 2, 3)], ntuple(i -> i == 1 ? 3 : 1, D))
        c = PolytropicElectrons(C(0.7), C(2^q), C(p / q))
        dx = ntuple(_ -> T(0.5), D)
        exact = Q(c.pe0) * sum(a -> Q(a, 2)^p, (1, 2, 3)) * prod(Q, dx) / (Q(c.γ) - 1)
        before = copy(n)
        @test isequal(P.complete(n, c, dx), H._from_exact(T, exact))
        @test n == before
    end
end
flush(stdout)
@testset "Complete range and stored gamma remain intact" begin
    for T in (Float32, Float64), sign in (-1, 1)
        k = T === Float32 ? 100 : 600
        density = ldexp(one(T), sign * k)
        spacing = ldexp(one(T), -sign * k)
        c = PolytropicElectrons(one(T), one(T), T(2))
        for n in (fill(density, 2), T[density, density/2])
            exact = sum(x -> Q(x)^2, n) * Q(spacing)
            @test isequal(P.complete(n, c, (spacing,)), H._from_exact(T, exact))
        end
    end
    for gamma in (prevfloat(1.0), nextfloat(1.0)), D = 1:3
        n = ones(Float32, ntuple(_ -> 2, D))
        dx = ntuple(_ -> 1.0f0, D)
        c = PolytropicElectrons(0.5, 1.0, gamma)
        exact = Q(c.pe0) * length(n) / (Q(gamma) - 1)
        @test isequal(P.complete(n, c, dx), H._from_exact(Float32, exact))
    end
end
flush(stdout)
@testset "Non-dyadic positive and negative sums at exact midpoints" begin
    for T in (Float32, Float64)
        a = T === Float32 ? one(T) + T(2)^-12 : one(T) + T(2)^-26
        b = T === Float32 ? one(T) + T(2)^-12 : one(T) + T(2)^-27
        boundary = Q(a) * Q(b)
        expected = H._from_exact(T, boundary)
        @test H._polytropic_midpoint_after(expected) == boundary
        @test P.complete(T[3, 4], PolytropicElectrons(a, T(5), T(2)), (b,)) === expected
        @test P.complete(T[9, 4], PolytropicElectrons(a / 2, T(25), T(0.5)), (b,)) === -expected
    end
end
flush(stdout)
@testset "Tiny positive terms retain midpoint direction without enormous shifts" begin
    a = 1.0 + 2.0^-26
    b = 1.0 + 2.0^-27
    lower = 1.0 + 3 * 2.0^-27
    expected = nextfloat(lower)
    for k in (26, 27, 28)
        gamma = 2.0^(2k)
        dx = (2.0^k - 1, 2.0^k + 1, b)
        n = zeros(2, 2, 2)
        n[1] = 1
        n[2] = 0.5
        c = PolytropicElectrons(a, 1.0, gamma)
        @test Q(c.pe0) * prod(Q, dx) / (Q(gamma) - 1) == (Q(lower) + Q(expected)) / 2
        @test P.complete(n, c, dx) === expected
        trial = @timed P.complete(n, c, dx)
        @test trial.value === expected
        println((k, bytes = trial.bytes, seconds = trial.time))
        flush(stdout)
    end
end
@testset "Finite zero signs and component rejection boundaries" begin
    for T in (Float32, Float64),
        count in (1, 2, 8),
        pe0 in (zero(T), -zero(T)),
        gamma in T[0.5, 2],
        density in T[0, 1, 2]

        c = PolytropicElectrons(pe0, one(T), gamma)
        n = fill(density, count)
        expected = signbit(pe0) ⊻ (gamma < 1) ? -zero(T) : zero(T)
        @test P.complete(n, c, (one(T),)) === expected
    end
    @test P.complete(Float64[], PolytropicElectrons(-0.0, 1.0, 0.5), (1.0,)) === -0.0
    @test_throws ArgumentError P.complete([1.0], PolytropicElectrons(1.0, 1.0, 1.0), (1.0,))
    @test_throws ArgumentError P.complete([Inf], PolytropicElectrons(1.0, 1.0, 2.0), (1.0,))
    @test_throws ArgumentError P.complete(
        [1.0],
        PolytropicElectrons{Float64}(1.0, -1.0, 2.0),
        (1.0,),
    )
end
println("COMPLETE_POLYTROPIC_COLD_ROUNDING_COMPLETE")
flush(stdout)

end
