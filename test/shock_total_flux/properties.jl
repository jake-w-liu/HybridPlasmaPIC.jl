module ShockTotalPropertiesTests
include("support.jl")
println((VERSION,))
flush(stdout)
@testset "Complete total neighboring cancellation and rounding boundaries" begin
    for T in (Float32, Float64),
        sign in (-1, 1),
        u in (T(2), T(4), T(16)),
        transverse in (zero(T), one(T), T(2))

        sh = PerpShock(8, T(7))
        sh.ux .= sign * u
        sh.uy .= transverse
        n = ldexp(one(T), exponent(floatmax(T)) - 1)
        sh.n .= n
        sh.Ey .= -n
        sh.Bz .= (u * u + transverse * transverse) * (sign * u) / 2
        for pressure in (nextfloat(zero(T)), T(0.125), T(1.5)), gamma in (T(1.5), T(2), T(5 / 3))
            sh.pe .= pressure
            sh.γe = gamma
            expected = reference(sh, nothing)
            actual = H.boundary_energy_flux(sh)
            @test all(isequal(actual.total[i], rounded(T, expected[i])) for i = 1:2)
        end
    end
    for T in (Float32, Float64), sign in (-1, 1)
        sh = PerpShock(8, T(7); γe = T(2))
        sh.ux .= T(sign)
        sh.Ey .= sign * floatmax(T)
        sh.Bz .= one(T)
        midpoint_pressure = ldexp(one(T), exponent(floatmax(T)) - precision(T) - 1)
        for pressure in
            (prevfloat(midpoint_pressure), midpoint_pressure, nextfloat(midpoint_pressure))
            sh.pe .= pressure
            exact = reference(sh, nothing)
            actual = H.boundary_energy_flux(sh)
            @test all(isequal(actual.total[i], rounded(T, exact[i])) for i = 1:2)
        end
        q = nextfloat(zero(T))
        sh.n .= 64q
        sh.ux .= T(sign * 0.25)
        sh.pe .= q
        sh.Ey .= sign * q
        for bz in (prevfloat(T(0.5)), T(0.5), nextfloat(T(0.5)))
            sh.Bz .= bz
            exact = reference(sh, nothing)
            actual = H.boundary_energy_flux(sh)
            @test all(isequal(actual.total[i], rounded(T, exact[i])) for i = 1:2)
        end
    end
end

@testset "Independent complete-field error bound" begin
    rng = MersenneTwister(291032)
    for T in (Float32, Float64), iteration = 1:500
        sh = PerpShock(8, T(7))
        span = fld(min(exponent(floatmax(T)), -exponent(floatmin(T))) - 4, 4)
        scale = T === Float32 ? 20 : 100
        for field in (sh.n, sh.ux, sh.uy, sh.pe, sh.Ey, sh.Bz)
            field .= ldexp(T(rand(rng, -16:16)) / 4, rand(rng, -scale:scale))
        end
        sh.γe = rand(rng, T[0.25, 0.5, 0.75, 1.25, 1.5, 5/3, 2, ldexp(one(T), span)])
        ordinary = Legacy.boundary_energy_flux(sh)
        exact = reference(sh, nothing)
        factors = (sh.γe, sh.γe - one(T), sh.pe[1], sh.ux[1], sh.Ey[1], sh.Bz[1], sh.n[1], sh.uy[1])
        @test H._shock_flux_factor_range_certified(T, factors)
        absolute = abs(ordinary.magnetic[1]) + abs(ordinary.kinetic[1]) + abs(ordinary.enthalpy[1])
        error = H._kinetic_energy_error_bound(absolute, Int128(4))
        @test Q(error) >= abs(Q(ordinary.total[1]) - exact[1])
    end
end

@testset "Ion state bound encloses the complete boundary expression" begin
    rng = MersenneTwister(718291)
    for T in (Float32, Float64), profile in (:ordinary, :extreme, :normalization), iteration = 1:80
        sh = PerpShock(8, T(7))
        ps = ParticleSet{1,T}(32)
        ps.x[1] .= T.(rand(rng, -8:64, 32)) ./ 8
        span = profile === :ordinary ? 3 : T === Float32 ? 70 : 400
        for a in (ps.weight, ps.v...)
            for i in eachindex(a)
                a[i] = ldexp(T(rand(rng, -16:16)) / 4, rand(rng, -span:span))
            end
        end
        if profile === :normalization
            he = T === Float32 ? 100 : 800
            sh.s.H .= ldexp(one(T), rand(rng, (-he, he)))
        end
        before = deepcopy((sh.s.H, ps.x, ps.v, ps.weight))
        actual, error = H._shock_ion_flux_state(sh, ps)
        @test isfinite(error) && error >= 0
        exact = reference(sh, ps)
        for (j, i) in ((1, 8), (2, 1))
            if isfinite(actual[i])
                normalization_error = nextfloat(0.5 * Float64(eps(actual[i])))
                radius = nextfloat(nextfloat(error / Float64(sh.s.H[i])) + normalization_error)
                @test isinf(radius) || Q(radius) >= abs(Q(actual[i]) - exact[j])
            else
                @test isequal(actual[i], rounded(T, exact[j]))
            end
        end
        @test isequal((sh.s.H, ps.x, ps.v, ps.weight), before)
    end
end
println("SHOCK_TOTAL_PROPERTIES_COMPLETE")
flush(stdout)

end
