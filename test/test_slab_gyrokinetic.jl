using Test
using HybridPlasmaPIC

@testset "Full-f slab gyrokinetic essentials" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64)
        @testset "$T" begin
            g = FourierGrid((4, 3, 3), (T(4), T(3), T(3)))
            p = GyrocentreSet{T}(prod(g.n))
            for (i, I) in enumerate(CartesianIndices(g.n)), d in 1:3
                p.X[d][i] = T(I[d] - 1)
            end
            # Base.mod may round a negative near-zero coordinate to L. The
            # public boundary must still use the canonical half-open interval.
            p.X[1][1] = -nextfloat(zero(T))
            fill!(p.vpar, T(1 / 8))
            species = [p, GyrocentreSet{T}(0)]
            st = SlabGyrokineticPIC(g, species; B0 = one(T), Te = one(T), nring = 4)
            init!(st, species)
            @test iszero(p.X[1][1])
            @test all(x -> 0 <= x < g.L[1], p.X[1])
            @test maximum(abs, st.rho_gc .- one(T)) <= T(16) * eps(T)
            @test maximum(abs, st.phi) <= T(64) * eps(T)

            # A zonal even-grid checkerboard has nonzero full-k polarization
            # despite the zero centered first derivative at Nyquist.
            for I in CartesianIndices(g.n)
                st.work.rho[I] = one(T) + (isodd(I[1]) ? T(1 / 8) : -T(1 / 8))
            end
            receipt = H._sg_field!(st)
            amplitude = T(1 / 8) / (st.chi * (T(2π) / T(4) * T(2))^2)
            expected = [isodd(I[1]) ? amplitude : -amplitude for I in CartesianIndices(g.n)]
            @test isapprox(st.work.phi, expected; rtol = T(64) * eps(T), atol = T(64) * eps(T) * amplitude)
            @test receipt.residual <= receipt.tolerance
            sample = H._sg_sample(fill(T(3), g.n), (T(0.2), T(0.3), T(0.4)), g)
            @test all(iszero, sample[2:4])

            before = H._checkpoint_state(st, species)
            @test step!(st, species, zero(T)) === st
            @test isequal(H._checkpoint_state(st, species), before)
            @test_throws ArgumentError step!(st, species, -one(T))
            @test isequal(H._checkpoint_state(st, species), before)
            old_z = copy(p.X[3])
            identity = p.X
            h = T(1 / 8)
            step!(st, species, h)
            @test p.X === identity
            @test st.time == h && st.step == 1 && st.last_dt == h
            @test isapprox(p.X[3], mod.(old_z .+ h * T(1 / 8), g.L[3]); atol = T(64) * eps(T), rtol = zero(T))
            @test isapprox(p.vpar, fill(T(1 / 8), length(p.vpar)); atol = T(64) * eps(T), rtol = zero(T))
            @test isapprox(st.rho_gc, fill(one(T), g.n); atol = T(64) * eps(T), rtol = zero(T))
            diagnostics = slab_gyrokinetic_diagnostics(st)
            @test diagnostics.range_resolved
            @test diagnostics.charge[1] == prod(g.n)

            mktempdir() do directory
                checkpoint = joinpath(directory, "slab.bin")
                archive = joinpath(directory, "slab.run")
                save_checkpoint(checkpoint, st, species)
                saved = H._checkpoint_state(st, species)
                step!(st, species, h)
                load_checkpoint!(st, species, checkpoint)
                @test isequal(H._checkpoint_state(st, species), saved)
                archive_run(archive, st, species; rng_seed = 0)
                step!(st, species, h)
                load_archive!(st, archive)
                @test isequal(H._checkpoint_state(st, species), saved)
                @test p.X === identity
            end
        end
    end
end
