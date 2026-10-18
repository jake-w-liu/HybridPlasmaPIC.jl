module ShockExactCacheTests
include("support.jl")
import JET
@testset "Partial exact nodes retain ownership and uncached boundary behavior" begin
    for T in (Float32, Float64)
        sh = PerpShock(8, T(7); γe = T(2))
        ps = ParticleSet{1,T}(3)
        ps.x[1] .= T[0, 0, 7]
        ps.v[1] .= T[1, -1, 2]
        sh.ux .= T(0.25)
        sh.pe .= 3nextfloat(zero(T))
        sh.Ey[8] = T(-8)
        state = H._shock_ion_flux_state(sh, ps)
        @test state.exact !== nothing
        @test haskey(state.exact, 1)
        @test state.exact[1] == 0
        @test !haskey(state.exact, 8)
        @test state.flux[8] == 8
        actual = H.boundary_energy_flux(sh; ps = ps)
        expected = reference(sh, ps)
        @test all(isequal(actual.total[i], rounded(T, expected[i])) for i = 1:2)
        @test actual.total == (2nextfloat(zero(T)), 2nextfloat(zero(T)))
        another = H._shock_ion_flux_state(sh, ps)
        @test another.exact !== state.exact
        @test another.flux !== state.flux
        state.exact[1] = Q(19)
        @test another.exact[1] == 0
        @test isequal(H.boundary_energy_flux(sh; ps = ps), actual)
        for p in (ParticleSet{1,T}(0), ps)
            report = JET.report_call(H._shock_ion_flux_state, Tuple{typeof(sh),typeof(p)})
            if !isempty(JET.get_reports(report))
                show(stdout, MIME"text/plain"(), report)
                println()
            end
            @test isempty(JET.get_reports(report))
        end
        # A concrete immutable result contains the optional exact nodes; it does
        # not introduce a shared cache or an ordinary output workspace.
        empty_state = H._shock_ion_flux_state(sh, ParticleSet{1,T}(0))
        @test empty_state.exact === nothing
        @test empty_state.error == 0.0
        @test all(iszero, empty_state.flux)
        output = Vector{typeof(another)}(undef, 32)
        Threads.@threads for i in eachindex(output)
            output[i] = H._shock_ion_flux_state(sh, ps)
        end
        for value in output
            @test isequal(value.flux, another.flux)
            @test value.exact == another.exact
            @test value.exact !== another.exact
        end
        @test all(output[i].exact !== output[j].exact for i in eachindex(output) for j = 1:i-1)
    end
end
println("SHOCK_EXACT_CACHE_COMPLETE")
flush(stdout)

end
