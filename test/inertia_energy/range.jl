module InertiaEnergyRangeTests
include("support.jl")
@testset "Complete inertia range and retained finite spectral current" begin
    for T in (Float32, Float64), sign in (-1, 1)
        e = T === Float32 ? 80 : 600
        ce = T === Float32 ? 100 : 800
        amplitude = ldexp(one(T), sign * e)
        de2 = ldexp(one(T), -sign * ce)
        g = FourierGrid((4,), (T(2pi),))
        B = (zeros(T, 4), zeros(T, 4), amplitude .* T[1, 0, -1, 0])
        J = ntuple(_ -> zeros(T, 4), 3)
        H.curl!(J, B, g)
        before = map(copy, J)
        value = H._electron_inertia_energy(J, de2, g.dx)
        expected = exact(J, de2, g.dx)
        @test Q(floatmin(T)) < expected < Q(floatmax(T))
        @test rounded(value, expected)
        @test all(isequal.(J, before))
    end
    for T in (Float32, Float64), D = 1:3
        e = T === Float32 ? 80 : 600
        small, large = ldexp(one(T), -e), ldexp(one(T), e)
        for dx in (ntuple(_ -> small, D), ntuple(_ -> large, D)),
            coefficient in (small, one(T), large)

            for value in (zero(T), small, one(T), large)
                J = (fill(value, 4), zeros(T, 4), zeros(T, 4))
                result = H._electron_inertia_energy(J, coefficient, dx)
                @test rounded(result, exact(J, coefficient, dx))
            end
        end
    end
    for T in (Float32, Float64)
        e = T === Float32 ? 80 : 600
        small, large = ldexp(one(T), -e), ldexp(one(T), e)
        for dx in (
            (small, small, large),
            (small, large, small),
            (large, small, small),
            (large, large, small),
            (large, small, large),
            (small, large, large),
        )
            J = (ones(T, 4), zeros(T, 4), zeros(T, 4))
            @test rounded(H._electron_inertia_energy(J, one(T), dx), exact(J, one(T), dx))
        end
    end
end
@testset "Inertia subnormal and finite overflow rounding neighbors" begin
    for T in (Float32, Float64)
        amplitude = ldexp(one(T), T === Float32 ? -75 : -538)
        J = ntuple(_ -> zeros(T, 4), 3)
        for count = 1:4, coefficient in (prevfloat(one(T)), one(T), nextfloat(one(T)))
            fill!(J[1], zero(T))
            J[1][1:count] .= amplitude
            @test rounded(
                H._electron_inertia_energy(J, coefficient, (T(2),)),
                exact(J, coefficient, (T(2),)),
            )
        end
        k = (exponent(floatmax(T)) - 1) ÷ 2
        tail = fld(2k - precision(T) - 2, 2)
        J = ntuple(_ -> zeros(T, 32), 3)
        for count in (0, 5, 6, 11, 12, 31),
            coefficient in (prevfloat(one(T)), one(T), nextfloat(one(T)))

            fill!(J[1], zero(T))
            J[1][1] = prevfloat(ldexp(one(T), k))
            J[1][2:count+1] .= ldexp(one(T), tail)
            @test rounded(
                H._electron_inertia_energy(J, coefficient, (T(8),)),
                exact(J, coefficient, (T(8),)),
            )
        end
    end
end
@testset "Ordinary inertia order allocation and nonfinite propagation" begin
    rng = MersenneTwister(902013)
    for T in (Float32, Float64),
        D = 1:3,
        count in (4, 16, 64),
        d in (T(0.1), one(T), T(2)^-20, T(2)^20)

        J = ntuple(_ -> randn(rng, T, ntuple(_ -> count, D)), 3)
        dx = ntuple(i -> T(i) / 2, D)
        before = map(copy, J)
        if T === Float32
            @test rounded(H._electron_inertia_energy(J, d, dx), exact(J, d, dx))
        else
            @test isequal(H._electron_inertia_energy(J, d, dx), original(J, d, dx))
        end
        cost(J, d, dx)
        @test cost(J, d, dx) == 0
        @test all(isequal.(J, before))
    end
    for T in (Float32, Float64), v in (T(Inf), T(-Inf), T(NaN))
        J = (T[v, 0], zeros(T, 2), zeros(T, 2))
        before = map(copy, J)
        @test isequal(
            H._electron_inertia_energy(J, T(0.25), (one(T),)),
            original(J, T(0.25), (one(T),)),
        )
        @test all(isequal.(J, before))
    end
end
println("INERTIA_COMPONENT_INITIAL_COMPLETE")
flush(stdout)

end
