module InertiaEnergyPropertyTests
include("support.jl")
@testset "Mixed extreme inertia factors use complete rational rounding cells" begin
    rng = MersenneTwister(981503)
    for T in (Float32, Float64), D = 1:3, iteration = 1:100
        e = T === Float32 ? 80 : 600
        ce = T === Float32 ? 100 : 800
        J = ntuple(_ -> [ldexp(T(rand(rng, -8:8)) / 8, rand(rng, -e:e)) for _ = 1:16], 3)
        J[1][1] = ldexp(one(T), e)
        coefficient = ldexp(one(T), rand(rng, -ce:0))
        dx = ntuple(_ -> ldexp(T(rand(rng, 1:8)) / 8, rand(rng, -e:e)), D)
        before = map(copy, J)
        expected = exact(J, coefficient, dx)
        @test isinf(J[1][1]^2)
        @test rounded(H._electron_inertia_energy(J, coefficient, dx), expected)
        @test all(isequal.(J, before))
    end
end
end
