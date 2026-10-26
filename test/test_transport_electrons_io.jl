using HybridPlasmaPIC, Test, Serialization

function transport_io_fixture(::Type{T}, ::Val{D}, collection; η = 0.03) where {T,D}
    g = FourierGrid(ntuple(d -> d == 1 ? 8 : 4, D), ntuple(_ -> T(2π), D))
    counts = ntuple(d -> 2g.n[d], D)
    ps = ParticleSet{D,T}(prod(counts))
    load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
    set_density_weight!(ps, one(T), g)
    ps.v[1] .= T(0.02) .* sin.(ps.x[1])
    ps.v[2] .= T(0.01)
    particles = collection ? [ps, ParticleSet{D,T}(0; q = T(2), m = T(4))] : ps
    st = HybridStepper(g, HybridModel(TransportElectrons(T(5 / 3)); η), CIC(), 0)
    fill!(st.fields.B[1], T(0.5))
    for I in CartesianIndices(g.n)
        x = T(I[1] - 1) * g.dx[1]
        st.fields.B[3][I] = T(0.1) * sin(x)
        st.fields.pe[I] = T(0.2) * (one(T) + T(0.1) * sin(x))
    end
    init!(st, particles)
    return st, particles
end

function transport_io_public_state(st, particles)
    fields = map(fieldnames(typeof(st.fields))) do name
        value = getfield(st.fields, name)
        value isa Ref ? value[] : value
    end
    populations = particles isa ParticleSet ? (particles,) : particles
    return deepcopy((
        st.model,
        st.shape,
        fields,
        st.pressure_work.Eforce,
        [(ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m) for ps in populations],
        st.time[],
        st.step[],
        st.last_dt[],
    ))
end

@testset "transport raw and archive restart" begin
    for T in (Float32, Float64), D = 1:3, collection in (false, true), evolved in (false, true)
        source, particles = transport_io_fixture(T, Val(D), collection)
        # Unstepped state checks the first-step priming phase after restore.
        evolved && step!(source, particles, T(0.003); NB = 2)
        mktempdir() do dir
            checkpoint = joinpath(dir, "state.ser")
            archive = joinpath(dir, "archive.ser")
            save_checkpoint(checkpoint, source, particles)
            archive_run(archive, source, particles; rng_seed = 17)
            raw = deserialize(checkpoint)
            archived = load_archive(archive).state
            @test raw.transport_policy_version === 1
            @test isequal(raw.Eforce, source.pressure_work.Eforce)
            @test isequal(raw.E, source.fields.E)
            @test isequal(raw, archived)
            for saved in (:path, raw, archived)
                # The restored model replaces the destination coefficient.
                target, restored = transport_io_fixture(T, Val(D), collection; η = 0.07)
                expected, full = deepcopy((source, particles))
                force_identity = target.pressure_work.Eforce
                foreach(a -> fill!(a, T(NaN)), force_identity)
                if saved === :path
                    @test load_checkpoint!(target, restored, checkpoint) === target
                else
                    HybridPlasmaPIC._restore_checkpoint_state!(target, restored, saved)
                end
                @test target.pressure_work.Eforce === force_identity
                @test target.model.η == 0.03
                @test isequal(HybridPlasmaPIC._checkpoint_state(target, restored), raw)
                for (dt, NB) in ((T(0.002), 1), (T(0.004), 3))
                    step!(target, restored, dt; NB)
                    step!(expected, full, dt; NB)
                    @test isequal(
                        transport_io_public_state(target, restored),
                        transport_io_public_state(expected, full),
                    )
                end
            end
            # The force is restored verbatim, not rebuilt from particle moments.
            perturbed = merge(raw, (Eforce = map(copy, raw.Eforce),))
            perturbed.Eforce[1][1] = nextfloat(perturbed.Eforce[1][1])
            target, restored = transport_io_fixture(T, Val(D), collection)
            HybridPlasmaPIC._restore_checkpoint_state!(target, restored, perturbed)
            @test isequal(target.pressure_work.Eforce, perturbed.Eforce)
            @test !isequal(target.pressure_work.Eforce, source.pressure_work.Eforce)
        end
    end
end

@testset "transport checkpoint rejects malformed state before mutation" begin
    for T in (Float32, Float64), collection in (false, true)
        source, particles = transport_io_fixture(T, Val(1), collection)
        step!(source, particles, T(0.001))
        state = deepcopy(HybridPlasmaPIC._checkpoint_state(source, particles))
        # Serialization can recreate isbits fields without invoking constructors.
        bad_gamma = reinterpret(TransportElectrons{T}, T[one(T)])[1]
        raw_model(closure; η = 0.03, ηH = 0.0, nfloor = 1e-6, de2 = 0.0) =
            HybridModel(closure, η, ηH, nfloor, de2)
        other = T === Float32 ? Float64 : Float32
        drop(key) = (; (k => v for (k, v) in pairs(state) if k != key)...)
        force = state.Eforce
        malformed = Any[
            drop(:Eforce),
            drop(:transport_policy_version),
            drop(:pe),
            merge(state, (transport_policy_version = 2,)),
            merge(state, (transport_policy_version = 1.0,)),
            merge(state, (transport_policy_version = Int32(1),)),
            merge(state, (Eforce = force[1:2],)),
            merge(state, (Eforce = (zeros(T, 7), force[2], force[3]),)),
            merge(state, (Eforce = map(a -> other.(a), force),)),
            merge(state, (Eforce = map(a -> view(copy(a), :), force),)),
            merge(state, (Eforce = (fill(T(NaN), 8), force[2], force[3]),)),
            merge(state, (Eforce = (force[1], fill(T(Inf), 8), force[3]),)),
            merge(state, (Eforce = (force[1], force[1], force[3]),)),
            merge(state, (Eforce = (state.E[1], force[2], force[3]),)),
            merge(state, (model = raw_model(state.model.closure; ηH = 0.1),)),
            merge(state, (model = raw_model(state.model.closure; de2 = 0.1),)),
            merge(state, (model = raw_model(state.model.closure; η = -0.1),)),
            merge(state, (model = raw_model(state.model.closure; η = NaN),)),
            merge(state, (model = raw_model(bad_gamma),)),
        ]
        collection || push!(malformed, merge(state, (q = T(2),)))
        target, restored = transport_io_fixture(T, Val(1), collection)
        before = transport_io_public_state(target, restored)
        mktempdir() do dir
            path = joinpath(dir, "bad.ser")
            for saved in malformed
                serialize(path, saved)
                @test_throws Union{ArgumentError,DimensionMismatch} load_checkpoint!(
                    target,
                    restored,
                    path,
                )
                @test isequal(transport_io_public_state(target, restored), before)
            end
            good = joinpath(dir, "good.ser")
            save_checkpoint(good, target, restored)
            old_bytes = read(good)
            saved_force = target.pressure_work.Eforce[3][2]
            target.pressure_work.Eforce[3][2] = T(NaN)
            @test_throws ArgumentError save_checkpoint(good, target, restored)
            @test read(good) == old_bytes
            @test_throws ArgumentError archive_run(
                joinpath(dir, "bad-archive.ser"),
                target,
                restored;
                rng_seed = 17,
            )
            @test !isfile(joinpath(dir, "bad-archive.ser"))
            target.pressure_work.Eforce[3][2] = saved_force
        end
        # A destination whose carried force aliases its fields is rejected.
        force_storage = target.pressure_work.Eforce
        target.pressure_work.Eforce = (target.fields.E[1], force_storage[2], force_storage[3])
        aliased = transport_io_public_state(target, restored)
        @test_throws ArgumentError HybridPlasmaPIC._restore_checkpoint_state!(
            target,
            restored,
            state,
        )
        @test isequal(transport_io_public_state(target, restored), aliased)
        @test_throws ArgumentError save_checkpoint(tempname(), target, restored)
    end
end

@testset "non-transport checkpoint payloads carry no force" begin
    g = FourierGrid((8,), (2π,))
    ps = ParticleSet{1,Float64}(8)
    load_lattice_1d!(ps, 0.0, 2π)
    set_density_weight!(ps, 1.0, g)
    for closure in (IsothermalElectrons(0.2), EvolvingElectrons(5 / 3))
        st = HybridStepper(g, HybridModel(closure), CIC(), 8)
        closure isa EvolvingElectrons && fill!(st.fields.pe, 0.3)
        init!(st, ps)
        saved = HybridPlasmaPIC._checkpoint_state(st, ps)
        @test !haskey(saved, :Eforce)
        @test !haskey(saved, :transport_policy_version)
        @test haskey(saved, :pe) == (closure isa EvolvingElectrons)
    end
end
