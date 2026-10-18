module ShockTotalThreadedTests
include("support.jl")
println((version = VERSION, threads = Threads.nthreads()))
flush(stdout)
@testset "Parent-checked concurrent complete boundary totals" begin
    for T in (Float32, Float64),
        kind in (:ordinary, :overflow, :subnormal),
        particles in (false, true)

        sh = PerpShock(8, T(7); γe = T(2))
        ps = ParticleSet{1,T}(2)
        ps.x[1] .= T[0, 7]
        if kind === :overflow
            n = ldexp(one(T), exponent(floatmax(T)) - 1)
            sh.n .= n
            sh.ux .= T(2)
            sh.pe .= T(0.25)
            sh.Bz .= T(4)
            sh.Ey .= -n
            ps.v[1] .= T(2)
            ps.weight .= n / 2
        elseif kind === :subnormal
            q = nextfloat(zero(T))
            sh.n .= 64q
            sh.ux .= T(0.25)
            sh.pe .= q
            sh.Ey .= q
            sh.Bz .= T(0.5)
            ps.weight .= 32q
            ps.v[1] .= T(0.25)
        else
            sh.n .= T(1.5)
            sh.ux .= T(-1.25)
            sh.uy .= T(0.5)
            sh.pe .= T(0.75)
            sh.Ey .= T(0.25)
            ps.v[1] .= T(-1.25)
            ps.weight .= T(0.75)
        end
        p = particles ? ps : nothing
        before = deepcopy((
            ps.x,
            ps.v,
            ps.weight,
            sh.n,
            sh.ux,
            sh.uy,
            sh.pe,
            sh.Ey,
            sh.Bz,
            sh.s.H,
            sh.γe,
        ))
        baseline = Legacy.boundary_energy_flux(sh; ps = p)
        expected =
            kind === :ordinary ? baseline :
            merge(baseline, (total = map(q -> rounded(T, q), reference(sh, p)),))
        output = Vector{typeof(expected)}(undef, 64)
        Threads.@threads for i in eachindex(output)
            output[i] = H.boundary_energy_flux(sh; ps = p)
        end
        for value in output
            @test isequal(value, expected)
        end
        @test isequal(
            (ps.x, ps.v, ps.weight, sh.n, sh.ux, sh.uy, sh.pe, sh.Ey, sh.Bz, sh.s.H, sh.γe),
            before,
        )
    end
end
println("SHOCK_TOTAL_THREADED_COMPLETE")
flush(stdout)

end
