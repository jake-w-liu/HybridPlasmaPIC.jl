using HybridPlasmaPIC, Test, Serialization

function evolving_io_fixture(::Type{T}, ::Val{D}, collection) where {T,D}
    g = FourierGrid(ntuple(d -> d == 1 ? 8 : 4, D), ntuple(_ -> T(2π), D))
    counts = ntuple(d -> 2g.n[d], D)
    ps = ParticleSet{D,T}(prod(counts))
    load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
    set_density_weight!(ps, one(T), g)
    ps.v[1] .= T(0.02) .* sin.(ps.x[1])
    ps.v[2] .= T(0.01)
    particles = collection ? [ps, ParticleSet{D,T}(0; q = T(2), m = T(4))] : ps
    st = HybridStepper(g, HybridModel(EvolvingElectrons(T(5 / 3))), CIC(), 0)
    fill!(st.fields.B[1], T(0.5))
    for I in CartesianIndices(g.n)
        st.fields.pe[I] = T(0.2) * (one(T) + T(0.1) * sin(T(I[1] - 1) * g.dx[1]))
    end
    init!(st, particles)
    return st, particles
end

function evolving_io_public_state(st, particles)
    fields = map(fieldnames(typeof(st.fields))) do name
        value = getfield(st.fields, name)
        value isa Ref ? value[] : value
    end
    populations = particles isa ParticleSet ? (particles,) : particles
    return deepcopy((
        st.model,
        st.shape,
        fields,
        [(ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m) for ps in populations],
        st.time[],
        st.step[],
        st.last_dt[],
    ))
end

@testset "evolving pressure raw and archive restart" begin
    for T in (Float32, Float64), D = 1:3, collection in (false, true)
        source, particles = evolving_io_fixture(T, Val(D), collection)
        initial = copy(source.fields.pe)
        step!(source, particles, T(0.003); NB = 2)
        @test source.fields.pe != initial
        mktempdir() do dir
            checkpoint = joinpath(dir, "state.ser")
            archive = joinpath(dir, "archive.ser")
            save_checkpoint(checkpoint, source, particles)
            archive_run(archive, source, particles; rng_seed = 17)
            raw = deserialize(checkpoint)
            archived = load_archive(archive).state
            @test haskey(raw, :pe)
            @test isequal(raw, archived)
            for saved in (raw, archived)
                target, restored = evolving_io_fixture(T, Val(D), collection)
                expected, full = deepcopy((source, particles))
                pressure_identity = target.fields.pe
                fill!(target.fields.pe, T(NaN)) # Valid saved state can replace invalid values.
                HybridPlasmaPIC._restore_checkpoint_state!(target, restored, saved)
                @test target.fields.pe === pressure_identity
                @test isequal(HybridPlasmaPIC._checkpoint_state(target, restored), raw)
                for (dt, NB) in ((T(0.002), 1), (T(0.004), 3))
                    step!(target, restored, dt; NB)
                    step!(expected, full, dt; NB)
                    @test isequal(
                        evolving_io_public_state(target, restored),
                        evolving_io_public_state(expected, full),
                    )
                end
            end
        end
    end
end

@testset "evolving pressure checkpoint rejects malformed state before mutation" begin
    for T in (Float32, Float64), collection in (false, true)
        source, particles = evolving_io_fixture(T, Val(1), collection)
        step!(source, particles, T(0.001))
        state = deepcopy(HybridPlasmaPIC._checkpoint_state(source, particles))
        # Serialization can recreate isbits fields without invoking constructors.
        # Reinterpret produces that malformed gamma without relaxing validation.
        bad_gamma = reinterpret(EvolvingElectrons{T}, T[one(T)])[1]
        raw_model(closure; η = 0.0, ηH = 0.0, nfloor = 1e-6, de2 = 0.0) =
            HybridModel(closure, η, ηH, nfloor, de2)
        malformed = (
            (; (k => v for (k, v) in pairs(state) if k != :pe)...),
            merge(state, (pe = zeros(T, 7),)),
            merge(state, (pe = zeros(T === Float32 ? Float64 : Float32, 8),)),
            merge(state, (pe = view(copy(state.pe), :),)),
            merge(state, (pe = fill(-one(T), 8),)),
            merge(state, (pe = fill(T(NaN), 8),)),
            merge(state, (pe = fill(T(Inf), 8),)),
            merge(state, (model = raw_model(bad_gamma),)),
            merge(state, (model = raw_model(state.model.closure; η = 0.1),)),
            merge(state, (model = raw_model(state.model.closure; ηH = 0.1),)),
            merge(state, (model = raw_model(state.model.closure; de2 = 0.1),)),
            merge(state, (model = raw_model(state.model.closure; nfloor = NaN),)),
        )
        target, restored = evolving_io_fixture(T, Val(1), collection)
        before = evolving_io_public_state(target, restored)
        mktempdir() do dir
            path = joinpath(dir, "bad.ser")
            for saved in malformed
                serialize(path, saved)
                @test_throws Union{ArgumentError,DimensionMismatch} load_checkpoint!(
                    target,
                    restored,
                    path,
                )
                @test isequal(evolving_io_public_state(target, restored), before)
            end
            good = joinpath(dir, "good.ser")
            save_checkpoint(good, target, restored)
            old_bytes = read(good)
            fill!(target.fields.pe, T(NaN))
            @test_throws ArgumentError save_checkpoint(good, target, restored)
            @test read(good) == old_bytes
            @test_throws ArgumentError archive_run(
                joinpath(dir, "bad-archive.ser"),
                target,
                restored;
                rng_seed = 17,
            )
            @test !isfile(joinpath(dir, "bad-archive.ser"))
        end
        target.fields.pe = target.fields.B[1]
        aliased = evolving_io_public_state(target, restored)
        @test_throws ArgumentError HybridPlasmaPIC._restore_checkpoint_state!(
            target,
            restored,
            state,
        )
        @test isequal(evolving_io_public_state(target, restored), aliased)
    end
end

@testset "algebraic checkpoint payload remains pressure-free" begin
    g = FourierGrid((8,), (2π,))
    ps = ParticleSet{1,Float64}(8)
    load_lattice_1d!(ps, 0.0, 2π)
    set_density_weight!(ps, 1.0, g)
    for closure in (IsothermalElectrons(0.2), PolytropicElectrons(0.2, 1.0, 5 / 3))
        st = HybridStepper(g, HybridModel(closure), CIC(), 8)
        init!(st, ps)
        saved = deepcopy(HybridPlasmaPIC._checkpoint_state(st, ps))
        @test !haskey(saved, :pe)
        target = HybridStepper(g, st.model, st.shape, 0)
        restored = ParticleSet{1,Float64}(0)
        HybridPlasmaPIC._restore_checkpoint_state!(target, restored, saved)
        @test isequal(HybridPlasmaPIC._checkpoint_state(target, restored), saved)
    end
end
