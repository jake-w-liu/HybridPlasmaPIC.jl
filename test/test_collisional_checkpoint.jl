using Test, Random, Serialization
using HybridPlasmaPIC
isdefined(@__MODULE__, :_composition_fixture) || include("collisional_fixtures.jl")

@testset "Collisional checkpoint validation and failure ownership" begin
    H = HybridPlasmaPIC
    @testset "transport force pressure and advanced RNG restart together $R etaH=$ηH kind=$kind" for R in (
            Xoshiro,
            MersenneTwister,
        ),
        ηH in (0.0, 0.01), kind in (:transport, :conduction, :wide_conduction)

        # Matching particle/grid lengths allow a real force/particle alias control.
        conducting = kind !== :transport
        closure = if kind === :wide_conduction
            setprecision(BigFloat, 256) do
                ConductingElectrons(BigFloat(5) / 3; κ = BigFloat(0.07))
            end
        else
            conducting ? ConductingElectrons(5 / 3; κ = 0.07) : TransportElectrons()
        end
        st, ps = _composition_transport_fixture(Float64, R; count = 4, ηH, closure)
        init!(st, ps)
        randn(st.rng, 3)
        step!(st, ps, st.dt)
        saved = deepcopy(H._checkpoint_state(st, ps))
        @test saved.schema === 1
        @test saved.inner.transport_policy_version === 1
        @test saved.inner.model.ηH == ηH
        conducting && @test saved.inner.conduction_policy_version === 2
        @test isequal(saved.inner.Eforce, st.inner.pressure_work.Eforce)
        dest, dps = _composition_transport_fixture(Float64, R; count = 4, ηH, closure)
        init!(dest, dps)
        entry = _composition_snapshot(dest, dps)
        drop(key) = (; (k => v for (k, v) in pairs(saved.inner) if k !== key)...)
        force = saved.inner.Eforce
        corruptions = Any[
            drop(:Eforce),
            drop(:transport_policy_version),
            merge(saved.inner, (transport_policy_version = 2,)),
            merge(
                saved.inner,
                (
                    model = HybridModel(
                        saved.inner.model.closure;
                        η = saved.inner.model.η,
                        ηH = ηH + 0.01,
                    ),
                ),
            ),
            merge(saved.inner, (Eforce = (fill(NaN, 4), force[2], force[3]),)),
            merge(saved.inner, (Eforce = (saved.fields.ui[1], force[2], force[3]),)),
            merge(saved.inner, (Eforce = (saved.inner.species[1].x[1], force[2], force[3]),)),
        ]
        if conducting
            changed_κ = kind === :wide_conduction ? nextfloat(closure.κ) : 0.08
            if kind === :wide_conduction
                @test changed_κ != closure.κ
                @test Float64(changed_κ) == Float64(closure.κ)
            end
            changed_closure = typeof(closure)(closure.γ; κ = changed_κ)
            changed_model = HybridModel(changed_closure;
                η = saved.inner.model.η, ηH)
            append!(corruptions, (
                drop(:conduction_policy_version),
                merge(saved.inner, (conduction_policy_version = 1,)),
                merge(saved.inner, (model = changed_model,)),
            ))
        end
        for inner in corruptions
            @test_throws ArgumentError H._restore_checkpoint_state!(
                dest,
                dps,
                merge(saved, (; inner)),
            )
            @test isequal(entry, _composition_snapshot(dest, dps))
        end
        corrupt = merge(saved, (fields = merge(saved.fields, (pe = 2 .* saved.inner.pe,)),))
        @test_throws ArgumentError H._restore_checkpoint_state!(dest, dps, corrupt)
        @test isequal(entry, _composition_snapshot(dest, dps))
        if R === MersenneTwister
            corrupt = deepcopy(saved)
            corrupt.rng.seed = corrupt.inner.Eforce[1]
            @test_throws ArgumentError H._restore_checkpoint_state!(dest, dps, corrupt)
            @test isequal(entry, _composition_snapshot(dest, dps))
        end
        force_identity = dest.inner.pressure_work.Eforce
        mktempdir() do dir
            path, archive = joinpath(dir, "transport.bin"), joinpath(dir, "transport-archive.bin")
            save_checkpoint(path, st, ps)
            @test load_checkpoint!(dest, dps, path) === dest
            @test dest.inner.pressure_work.Eforce === force_identity
            @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
            for _ = 1:2
                step!(st, ps, st.dt)
                step!(dest, dps, dest.dt)
                @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
            end
            @test randn(dest.rng, 7) == randn(st.rng, 7)
            archive_run(archive, st, ps; rng_seed = 712)
            archived = load_archive(archive).state
            @test isequal(archived.inner.Eforce, st.inner.pressure_work.Eforce)
            @test H._restore_checkpoint_state!(dest, dps, archived) === dest
            @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
        end
    end

    @testset "evolving pressure and advanced RNG continue together $R" for R in (
        Xoshiro,
        MersenneTwister,
    )

        make() = _composition_fixture(Val(1), Float64, R; closure = EvolvingElectrons())
        st, ps = make()
        st.inner.fields.pe .= [0.02, 0.022, 0.02, 0.018]
        init!(st, ps)
        randn(st.rng, 3)
        step!(st, ps, st.dt)
        saved = deepcopy(H._checkpoint_state(st, ps))
        dest, dps = make()
        init!(dest, dps)
        entry = _composition_snapshot(dest, dps)
        for pressure in (fill(-1.0, 4), 2 .* saved.inner.pe)
            # Detach this leaf: the writer intentionally shares the two pe arrays.
            corrupt = merge(saved, (fields = merge(saved.fields, (pe = pressure,)),))
            @test_throws ArgumentError H._restore_checkpoint_state!(dest, dps, corrupt)
            @test isequal(entry, _composition_snapshot(dest, dps))
        end
        missing = NamedTuple{Tuple(k for k in keys(saved.inner) if k !== :pe)}(saved.inner)
        @test_throws ArgumentError H._restore_checkpoint_state!(
            dest,
            dps,
            merge(saved, (inner = missing,)),
        )
        @test isequal(entry, _composition_snapshot(dest, dps))
        pressure_identity = dest.inner.fields.pe
        mktempdir() do dir
            path = joinpath(dir, "evolving.bin")
            save_checkpoint(path, st, ps)
            @test load_checkpoint!(dest, dps, path) === dest
            @test dest.inner.fields.pe === pressure_identity
            @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
            for _ = 1:2
                step!(st, ps, st.dt)
                step!(dest, dps, dest.dt)
                @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
            end
            @test randn(dest.rng, 7) == randn(st.rng, 7)
            archive = joinpath(dir, "evolving-archive.bin")
            archive_run(archive, st, ps; rng_seed = 712)
            state = load_archive(archive).state
            @test state.inner.pe == st.inner.fields.pe == state.fields.pe
        end
    end

    @testset "malformed payload rejects before publication $R" for R in (Xoshiro, MersenneTwister)
        st, ps = _composition_fixture(Val(1), Float64, R)
        init!(st, ps)
        step!(st, ps, st.dt)
        saved = deepcopy(H._checkpoint_state(st, ps))
        dest, dps = _composition_fixture(Val(1), Float64, R)
        init!(dest, dps)
        entry = _composition_snapshot(dest, dps)
        corruptions = Any[
            saved.inner,
            NamedTuple{Tuple(k for k in keys(saved) if k !== :rng)}(saved),
            merge(saved, (schema = 2,)),
            merge(saved, (phase = :collision,)),
            merge(saved, (initialized = false,)),
            merge(saved, (source_count = saved.source_count + 1,)),
            merge(saved, (source_time = nextfloat(saved.source_time),)),
            merge(saved, (dt = 2saved.dt,)),
            merge(saved, (coulomb_log = 2saved.coulomb_log,)),
            merge(saved, (epsilon0 = 2saved.epsilon0,)),
            merge(saved, (weight_to_number = 2saved.weight_to_number,)),
            merge(saved, (u_floor = 2saved.u_floor,)),
            merge(saved, (runtime = merge(saved.runtime, (julia = v"0.0.0",)),)),
            merge(saved, (floor_count = -1,)),
            merge(saved, (inner = merge(saved.inner, (last_dt = 2saved.dt,)),)),
            merge(saved, (inner = merge(saved.inner, (species = reverse(saved.inner.species),)),)),
        ]
        other = R === Xoshiro ? MersenneTwister(1) : Xoshiro(1)
        push!(corruptions, merge(saved, (rng = other, rng_algorithm = typeof(other))))
        corrupt = deepcopy(saved)
        corrupt.fields.E[1][1] = NaN
        push!(corruptions, corrupt)
        corrupt = deepcopy(saved)
        push!(
            corruptions,
            merge(corrupt, (inner = merge(corrupt.inner, (E = ntuple(_ -> zeros(4), 3),)),)),
        )
        if R === MersenneTwister
            for change in (
                r -> resize!(r.state.val, 1),
                r -> resize!(r.vals, 1),
                r -> resize!(r.ints, 1),
                r -> (r.idxF = -1),
                r -> (r.idxI = -1),
                r -> (r.idxF = 0; r.vals[1] = NaN),
            )
                corrupt = deepcopy(saved)
                change(corrupt.rng)
                push!(corruptions, corrupt)
            end
            # Serialization changes these objects before load sees them: Ptr
            # becomes NULL and external arrays detach. Check the actual live
            # incoming-graph boundary, then independently verify that transform.
            for kind in (:pointer, :live_array)
                corrupt = deepcopy(saved)
                corrupt.rng.seed = kind === :pointer ? Ptr{UInt8}(1) : dest.inner.fields.B[1]
                @test_throws ArgumentError H._validate_checkpoint_state(corrupt, dest)
                @test_throws ArgumentError H._restore_checkpoint_state!(dest, dps, corrupt)
                @test isequal(entry, _composition_snapshot(dest, dps))
                buffer = IOBuffer()
                serialize(buffer, corrupt)
                seekstart(buffer)
                transformed = deserialize(buffer)
                if kind === :pointer
                    @test transformed.rng.seed === Ptr{UInt8}(C_NULL)
                else
                    @test transformed.rng.seed !== dest.inner.fields.B[1]
                    @test transformed.rng.seed == dest.inner.fields.B[1]
                end
                @test H._validate_checkpoint_state(transformed, dest) === nothing
            end
        else
            corrupt = deepcopy(saved)
            corrupt.rng.s0 = corrupt.rng.s1 = corrupt.rng.s2 = corrupt.rng.s3 = 0
            push!(corruptions, corrupt)
        end
        mktempdir() do dir
            path = joinpath(dir, "state.bin")
            for corrupt in corruptions
                serialize(path, corrupt)
                @test_throws Exception load_checkpoint!(dest, dps, path)
                @test isequal(entry, _composition_snapshot(dest, dps))
                @test !dest.busy[]
            end
            write(path, "not a serialized checkpoint")
            @test_throws Exception load_checkpoint!(dest, dps, path)
            @test isequal(entry, _composition_snapshot(dest, dps))
            save_checkpoint(path, st, ps)
            bytes = read(path)
            st.phase = :inflight
            @test_throws ArgumentError save_checkpoint(path, st, ps)
            @test read(path) == bytes
            @test !st.busy[]
            st.phase = :ready
            @test_throws Exception save_checkpoint(dir, st, ps)
            @test read(path) == bytes
            @test !st.busy[]
            @test load_checkpoint!(dest, dps, path) === dest
            @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
        end
    end

    @testset "RNG detachment includes mutable seed metadata" begin
        st, ps = _composition_fixture(Val(1), Float64, MersenneTwister)
        rng = MersenneTwister(UInt32[1, 3, 5])
        c = CollisionalHybridStepper(
            st.inner,
            ps;
            dt = st.dt,
            rng,
            coulomb_log = st.coulomb_log,
            epsilon0 = st.epsilon0,
            weight_to_number = st.weight_to_number,
        )
        state = _composition_rng_state(c.rng)
        rng.seed[1] += 1
        @test isequal(state, _composition_rng_state(c.rng))
        init!(c, ps)
        step!(c, ps, c.dt)
        saved = deepcopy(H._checkpoint_state(c, ps))
        dest, dps = _composition_fixture(Val(1), Float64, MersenneTwister)
        H._restore_checkpoint_state!(dest, dps, saved)
        state = _composition_rng_state(dest.rng)
        saved.rng.seed[1] += 1
        fill!(saved.rng.vals, 0)
        @test isequal(state, _composition_rng_state(dest.rng))
        @test rand(dest.rng, UInt64, 25) == rand(c.rng, UInt64, 25)
        @test randn(dest.rng, 25) == randn(c.rng, 25)
    end

    @testset "initialization failure restores full public state" for T in (Float32, Float64)
        st, ps = _composition_fixture(Val(1), T, Xoshiro)
        for p in ps
            fill!(p.x[1], T(0.5))
            fill!(p.weight, floatmax(T))
        end
        before = _composition_snapshot(st, ps)
        @test_throws Exception init!(st, ps)
        @test isequal(before, _composition_snapshot(st, ps))
        @test !st.busy[]
        @test st.phase === :uninitialized
    end

    @testset "repeated output retains no memory" begin
        # Per-call bytes of save/load/archive depend on identity-hash table growth
        # for freshly allocated payload objects, so identical states need not
        # allocate identical bytes. The output contract is bounded retention: one
        # retained payload per cycle (about 1 MB here) would exceed 1 MiB 30-fold.
        make() = _composition_fixture(
            Val(2),
            Float64,
            Xoshiro;
            ns = 4,
            count = 4096,
            temperature = 0.0,
            magnetic = 0.0,
        )
        st, ps = make()
        dest, dps = make()
        init!(st, ps)
        step!(st, ps, st.dt)
        init!(dest, dps)
        entry = _composition_snapshot(st, ps)
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "state.bin"), joinpath(dir, "archive.bin")
            function cycle()
                save_checkpoint(checkpoint, st, ps)
                load_checkpoint!(dest, dps, checkpoint)
                archive_run(archive, st, ps; rng_seed = 712)
            end
            for _ = 1:8
                cycle()
            end
            GC.gc(true)
            GC.gc(true)
            live = Base.gc_live_bytes()
            for _ = 1:32
                cycle()
            end
            GC.gc(true)
            GC.gc(true)
            @test Base.gc_live_bytes() - live <= 2^20
        end
        @test isequal(_composition_snapshot(st, ps), entry)
        @test isequal(_composition_snapshot(dest, dps), entry)
    end
end
