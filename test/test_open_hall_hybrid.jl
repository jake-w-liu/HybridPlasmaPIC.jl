using Test
using Random
using HybridPlasmaPIC

@testset "Synchronous open Hall essentials" begin
    for T in (Float32, Float64)
        @testset "$T" begin
            op = MixedShockOperators((3, 3), (T(2), T(3)))
            p = ParticleSet{2,T}(prod(op.n))
            for (i, I) in enumerate(CartesianIndices(op.n))
                for d in 1:2
                    p.x[d][i] = op.coordinates[d][I[d]]
                end
                p.weight[i] = op.sbp.H[I[1]] * op.dx[2]
            end
            # Unit inward speed, density 8, face area 3 and dt=1/8 give
            # exactly three births. No old particle moves or exits initially.
            left = OpenHallReservoir(T(8), one(T), zero(T))
            right = OpenHallReservoir(zero(T), zero(T), zero(T))
            species = [p]
            rng = T === Float32 ? Xoshiro(41) : MersenneTwister(41)
            st = OpenHallHybrid(op, HybridModel(IsothermalElectrons(zero(T))),
                species, [(left, right)]; rng, max_levels = 0, max_births = 3,
                max_particles = 32, max_draws = 1024)
            init!(st, species)
            snapshot() = (x = map(copy, p.x), v = map(copy, p.v),
                weight = copy(p.weight), id = copy(p.id), tag = copy(p.tag),
                fields = deepcopy(HybridPlasmaPIC._hybrid_field_state(st.fields)),
                time = st.time[], step = st.step[], last_dt = st.last_dt[],
                alpha = copy(st.alpha), next_index = copy(st.next_index),
                elapsed = st.source_elapsed, future = rand(deepcopy(st.rng), UInt64, 8))
            initial = snapshot()
            @test step!(st, species, zero(T)) === st
            @test isequal(snapshot(), initial)
            @test_throws ArgumentError step!(st, species, -one(T))
            @test isequal(snapshot(), initial)
            h = T(1 / 8)
            step!(st, species, h)
            @test st.time[] == h && st.step[] == 1 && st.last_dt[] == h
            @test nparticles(p) == prod(op.n) + 3
            @test st.births[1].count == st.last_births[1].count == 3
            @test all(e -> iszero(e.count), st.exits)
            @test st.last_report.births == 3 && st.last_report.exits == 0
            @test 0 < st.last_report.uniform_requests <= st.controls.max_draws
            @test st.next_index[1] == prod(op.n) + 4
            @test all(a -> all(isfinite, a), st.fields.B)
            diagnostics = open_hall_diagnostics(st)
            @test diagnostics.particle_totals.count == nparticles(p)
            @test diagnostics.source.source_elapsed == h
            mktempdir() do directory
                checkpoint = joinpath(directory, "hall.bin")
                archive = joinpath(directory, "hall.run")
                save_checkpoint(checkpoint, st, species)
                saved = snapshot()
                step!(st, species, h)
                @test st.births[1].count == 6 && st.last_births[1].count == 3
                load_checkpoint!(st, species, checkpoint)
                @test isequal(snapshot(), saved)
                archive_run(archive, st, species; rng_seed = 41)
                step!(st, species, h)
                load_archive!(st, species, archive)
                @test isequal(snapshot(), saved)
            end
        end
    end
end
