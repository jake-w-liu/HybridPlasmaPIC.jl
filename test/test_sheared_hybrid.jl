using Test
using HybridPlasmaPIC

@testset "Coupled periodic shear essentials" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64)
        @testset "$T" begin
            g = FourierGrid((5, 5, 5), (T(5), T(5), T(5)))
            mesh = PeriodicShearGrid(g; amplitude = T(1 / 8))
            point = (T(1 / 4), T(1 / 2), T(3 / 4))
            returned = to_logical(mesh, to_cartesian(mesh, point))
            @test all(d -> isapprox(returned[d], point[d]; atol = T(16) * eps(T), rtol = zero(T)), 1:3)
            @test_throws ArgumentError PeriodicShearGrid(FourierGrid((4, 5, 5), g.L); amplitude = zero(T))
            p = ParticleSet{3,T}(prod(g.n); q = one(T), m = one(T))
            load_sheared_uniform!(p, mesh; number_density = one(T), velocity = (T(1 / 8), zero(T), T(1 / 8)))
            species = [p, ParticleSet{3,T}(0; q = one(T), m = one(T))]
            st = ShearedPeriodicHybridPIC(mesh, species; Te = T(1 / 8), B0 = (one(T), zero(T), one(T)))
            init!(st, species)
            @test all(isone, st.fields.n)
            @test maximum(a -> maximum(abs, a), st.fields.E) <= T(128) * eps(T)
            @test st.last_report.divergence <= st.last_report.tolerance
            @test_throws ArgumentError H._sh_report(ShearedHybridFieldReport{T}(T(2), one(T), one(T), false, true))
            @test H._sh_report(ShearedHybridFieldReport{T}(T(Inf), T(Inf), T(Inf), false, false)) === nothing
            rate = H._sheared_lorentz_rate(T(2), T(4), (one(T), T(2), T(3)), (zero(T), zero(T), zero(T)), (zero(T), zero(T), one(T)))
            @test rate == (T(1 / 2), one(T), T(3 / 2))
            before = H._checkpoint_state(st, species)
            @test step!(st, species, -zero(T)) === st
            @test isequal(H._checkpoint_state(st, species), before)
            @test_throws ArgumentError step!(st, species, -one(T))
            @test isequal(H._checkpoint_state(st, species), before)
            old_x = map(copy, p.x)
            old_v = map(copy, p.v)
            identity = p.x
            h = T(1 / 8)
            step!(st, species, h)
            @test p.x === identity
            @test st.time == h && st.step == 1 && st.last_dt == h
            for d in 1:3
                @test isapprox(p.x[d], mod.(old_x[d] .+ h .* old_v[d], g.L[d]); atol = T(128) * eps(T), rtol = zero(T))
                @test isapprox(p.v[d], old_v[d]; atol = T(128) * eps(T), rtol = zero(T))
            end
            diagnostics = sheared_hybrid_diagnostics(st)
            @test diagnostics.range_resolved
            @test diagnostics.charge == [big(125)//1, big(0)//1]
            @test diagnostics.quadrature_volume_defect == 0
            mktempdir() do directory
                checkpoint = joinpath(directory, "sheared.bin")
                archive = joinpath(directory, "sheared.run")
                save_checkpoint(checkpoint, st, species)
                saved = H._checkpoint_state(st, species)
                step!(st, species, h)
                load_checkpoint!(st, checkpoint)
                @test isequal(H._checkpoint_state(st, species), saved)
                archive_run(archive, st; rng_seed = 0)
                step!(st, species, h)
                load_archive!(st, archive)
                @test isequal(H._checkpoint_state(st, species), saved)
                @test p.x === identity
            end
        end
    end
end
