module ShockTotalRangeTests
include("support.jl")
println((VERSION,))
flush(stdout)
@testset "Complete public shock totals with unchanged components" begin
    for T in (Float32, Float64), sign in (-1, 1)
        sh = PerpShock(8, T(7); γe = T(2))
        n = ldexp(one(T), exponent(floatmax(T)) - 1)
        sh.n .= n
        sh.ux .= T(2sign)
        sh.pe .= T(0.25)
        sh.Bz .= T(4)
        sh.Ey .= -sign * n
        ps = ParticleSet{1,T}(2)
        ps.x[1] .= T[0, 7]
        ps.v[1] .= T(2sign)
        ps.weight .= n / 2
        for p in (nothing, ps)
            before = deepcopy((sh.n, sh.ux, sh.pe, sh.Bz, sh.Ey, ps.x, ps.v, ps.weight))
            expected = reference(sh, p)
            baseline = Legacy.boundary_energy_flux(sh; ps = p)
            actual = H.boundary_energy_flux(sh; ps = p)
            @test all(isequal(actual.total[i], rounded(T, expected[i])) for i = 1:2)
            @test isequal(
                (actual.magnetic, actual.kinetic, actual.enthalpy),
                (baseline.magnetic, baseline.kinetic, baseline.enthalpy),
            )
            @test isequal((sh.n, sh.ux, sh.pe, sh.Bz, sh.Ey, ps.x, ps.v, ps.weight), before)
        end
        q = nextfloat(zero(T))
        sh.n .= 64q
        sh.ux .= T(sign * 0.25)
        sh.pe .= q
        sh.Ey .= sign * q
        sh.Bz .= T(0.5)
        ps.weight .= 32q
        ps.v[1] .= T(sign * 0.25)
        for p in (nothing, ps)
            actual = H.boundary_energy_flux(sh; ps = p)
            expected = reference(sh, p)
            @test all(isequal(actual.total[i], rounded(T, expected[i])) for i = 1:2)
        end
    end
end
@testset "Ordinary public totals and nodal allocation remain compatible" begin
    rng = MersenneTwister(89231)
    for T in (Float32, Float64), nodes in (8, 64, 512), count in (0, 1, 64, 8192)
        sh = PerpShock(nodes, T(nodes - 1); γe = T(1.5))
        for field in (sh.n, sh.pe, sh.Bz)
            field .= T.(1 .+ rand(rng, nodes))
        end
        for field in (sh.ux, sh.uy, sh.Ey)
            field .= T.(randn(rng, nodes))
        end
        ps = ParticleSet{1,T}(count)
        ps.x[1] .= T.(rand(rng, count) .* (nodes - 1))
        ps.weight .= T.(rand(rng, count))
        for v in ps.v
            v .= T.(randn(rng, count))
        end
        for p in (nothing, ps)
            old(sh, p)
            new(sh, p)
            a, b = old(sh, p), new(sh, p)
            @test isequal(a.value, b.value)
            @test a.bytes == b.bytes
            println((
                T,
                nodes,
                count,
                particles = p !== nothing,
                old_bytes = a.bytes,
                new_bytes = b.bytes,
            ))
            flush(stdout)
        end
        node_old(sh, ps)
        node_new(sh, ps)
        @test node_old(sh, ps) == node_new(sh, ps)
    end
    for T in (Float32, Float64), count in (0, 1, 64)
        sh = PerpShock(8, T(7))
        ps = ParticleSet{1,T}(count)
        ps.x[1] .= T(3.5)
        for v in ps.v
            v .= one(T)
        end
        for p in (nothing, ps)
            old(sh, p)
            new(sh, p)
            a, b = old(sh, p), new(sh, p)
            @test isequal(a.value, b.value)
            @test a.bytes == b.bytes
        end
    end
end
println("SHOCK_TOTAL_INITIAL_COMPLETE")
flush(stdout)

end
