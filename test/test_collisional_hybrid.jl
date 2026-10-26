using Test, Random, Serialization
using HybridPlasmaPIC
isdefined(@__MODULE__, :_composition_fixture) || include("collisional_fixtures.jl")

@testset "Collisional hybrid first-order owner" begin
    H = HybridPlasmaPIC
    @testset "canonical sequence D=$D T=$T RNG=$R logs=$logs" for D = 1:3,
        T in (Float32, Float64),
        R in (Xoshiro, MersenneTwister),
        logs in (0.0, 0.05)

        st, ps = _composition_fixture(Val(D), T, R; logs = T(logs))
        reference, rps = _composition_fixture(Val(D), T, R; logs = T(logs))
        entry = _composition_snapshot(st, ps)
        @test step!(st, ps, zero(T)) === st
        @test isequal(_composition_snapshot(st, ps), entry)
        @test_throws ArgumentError step!(st, ps, st.dt)
        init!(st, ps)
        init!(reference, rps)
        original_rng = _composition_rng_state(st.rng)
        for _ = 1:3
            @test step!(st, ps, st.dt) === st
            _composition_component!(reference, rps)
            @test isequal(_composition_physical(st, ps), _composition_physical(reference, rps))
            @test st.source_count == st.inner.step[]
            @test isequal(st.source_time, st.inner.time[])
        end
        iszero(logs) && @test isequal(_composition_rng_state(st.rng), original_rng)
        rng_before = _composition_rng_state(st.rng)
        init!(st, ps)
        @test isequal(_composition_rng_state(st.rng), rng_before)
        @test st.source_count == 0
        @test st.source_time == zero(T)
    end

    @testset "constructor ownership and rejection" for R in (Xoshiro, MersenneTwister)
        st, ps = _composition_fixture(Val(1), Float64, R)
        rng = R(44)
        randn(rng, 7)
        rand(rng, UInt64, 3)
        before = _composition_rng_state(rng)
        logs = [0.1 0.2; 0.2 0.3]
        c = CollisionalHybridStepper(
            st.inner,
            ps;
            dt = 0.01,
            rng,
            coulomb_log = logs,
            epsilon0 = 1.0,
            weight_to_number = 1.0,
        )
        @test isequal(_composition_rng_state(rng), before)
        @test isequal(_composition_rng_state(c.rng), before)
        @test c.species !== ps
        logs[1, 1] = 99
        @test c.coulomb_log[1, 1] == 0.1
        rand(rng)
        @test !isequal(_composition_rng_state(rng), _composition_rng_state(c.rng))
        init!(c, ps)
        before = _composition_snapshot(c, ps)
        for bad_dt in (-1.0, Inf, NaN, 2c.dt, BigFloat(2)^(-2000))
            @test_throws ArgumentError step!(c, ps, bad_dt)
            @test isequal(_composition_snapshot(c, ps), before)
        end
        @test_throws ArgumentError step!(c, ps, c.dt; NB = 0)
        @test_throws ArgumentError step!(c, reverse(ps), c.dt)
        @test_throws ArgumentError step!(c, deepcopy(ps), c.dt)
        c.busy[] = true
        @test_throws ArgumentError step!(c, ps, c.dt)
        @test_throws ArgumentError init!(c, ps)
        @test step!(c, ps, 0) === c
        c.busy[] = false
        @test isequal(_composition_snapshot(c, ps), before)
        c.source_count += 1
        @test_throws ArgumentError step!(c, ps, c.dt)
        c.source_count -= 1
        ps[1].q *= 2
        @test_throws ArgumentError step!(c, ps, c.dt)
        ps[1].q /= 2
        @test isequal(_composition_snapshot(c, ps), before)
        c.coulomb_log[1, 1] = 0.2
        @test_throws ArgumentError step!(c, ps, c.dt)
        c.coulomb_log[1, 1] = 0.1
        @test isequal(_composition_snapshot(c, ps), before)
    end

    @testset "other supported particle shapes" for D = 1:3,
        T in (Float32, Float64),
        shape in (NGP(), TSC())

        st, ps = _composition_fixture(Val(D), T, Xoshiro; shape)
        reference, rps = _composition_fixture(Val(D), T, Xoshiro; shape)
        init!(st, ps)
        init!(reference, rps)
        for _ = 1:3
            step!(st, ps, st.dt)
            _composition_component!(reference, rps)
            @test isequal(_composition_physical(st, ps), _composition_physical(reference, rps))
        end
    end

    @testset "constructor model and ownership boundary" begin
        st, ps = _composition_fixture(Val(1), Float64, Xoshiro)
        make(inner; dt = st.dt, rng = Xoshiro(1)) = CollisionalHybridStepper(
            inner,
            ps;
            dt,
            rng,
            coulomb_log = 0.05,
            epsilon0 = 1.0,
            weight_to_number = 1.0,
        )
        for model in (
            HybridModel(IsothermalElectrons(0.02); η = eps()),
            HybridModel(IsothermalElectrons(0.02); ηH = eps()),
            HybridModel(IsothermalElectrons(0.02); de2 = eps()),
            HybridModel(PolytropicElectrons(1.0, 1.0, 5 / 3)),
            HybridModel(EvolvingElectrons()),
        )
            inner = HybridStepper(st.inner.g, model, CIC(), sum(length(p.weight) for p in ps))
            @test_throws ArgumentError make(inner)
        end
        for dt in (0.0, -1.0, Inf, NaN, BigFloat(2)^(-2000))
            @test_throws ArgumentError make(st.inner; dt)
        end
        @test_throws ArgumentError make(st.inner; rng = RandomDevice())
        init!(st, ps)
        step!(st, ps, st.dt)
        before = _composition_snapshot(st, ps)
        @test_throws ArgumentError make(st.inner)
        @test isequal(before, _composition_snapshot(st, ps))
    end

    @testset "advanced RNG and complete checkpoint" for R in (Xoshiro, MersenneTwister),
        boundary = 0:3

        st, ps = _composition_fixture(Val(2), Float64, R)
        init!(st, ps)
        # Both scalar and bulk paths can leave native caches partially consumed.
        randn(st.rng, 11)
        rand(st.rng, UInt8, 5)
        for _ = 1:boundary
            step!(st, ps, st.dt)
        end
        mktempdir() do dir
            path = joinpath(dir, "state.bin")
            @test save_checkpoint(path, st, ps) == path
            dest, dps = _composition_fixture(Val(2), Float64, R)
            arrays = [(p.x..., p.v..., p.weight, p.id, p.tag) for p in dps]
            fields = H._hybrid_field_state(dest.inner.fields)
            @test load_checkpoint!(dest, dps, path) === dest
            @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
            for (p, identities) in zip(dps, arrays)
                @test all(
                    a === b for (a, b) in zip((p.x..., p.v..., p.weight, p.id, p.tag), identities)
                )
            end
            @test all(
                a === b for
                (a, b) in zip(values(fields), values(H._hybrid_field_state(dest.inner.fields)))
            )
            for _ = 1:3
                step!(st, ps, st.dt)
                step!(dest, dps, dest.dt)
                @test isequal(_composition_snapshot(dest, dps), _composition_snapshot(st, ps))
            end
            @test rand(st.rng, UInt64, 19) == rand(dest.rng, UInt64, 19)
            @test randn(st.rng, 17) == randn(dest.rng, 17)
            archive = joinpath(dir, "archive.bin")
            @test archive_run(archive, st, ps; rng_seed = 712) == archive
            receipt = load_archive(archive)
            @test receipt.state.rng_algorithm === R
            @test receipt.state.source_count == st.source_count
        end
    end
    @testset "fixed grid and independent public storage" begin
        st, ps = _composition_fixture(Val(1), Float64, Xoshiro; count = 4, logs = 0.0)
        init!(st, ps)
        original = st.inner.g
        st.inner.g = FourierGrid((4,), (4.0,))
        before = _composition_snapshot(st, ps)
        @test_throws ArgumentError step!(st, ps, st.dt)
        @test isequal(before, _composition_snapshot(st, ps))
        st.inner.g = original
        oldB = st.inner.fields.B
        st.inner.fields.B = (ps[1].x[1], oldB[2], oldB[3])
        before = _composition_snapshot(st, ps)
        @test_throws ArgumentError init!(st, ps)
        @test isequal(before, _composition_snapshot(st, ps))
        st.inner.fields.B = (oldB[1], oldB[1], oldB[3])
        @test_throws ArgumentError init!(st, ps)
        st.inner.fields.B = oldB
        @test step!(st, ps, st.dt) === st
    end

    @testset "retained diagonal-load field failure rolls back" for T in (Float32, Float64)
        st, ps = _composition_fixture(Val(3), T, Xoshiro; count = 16, logs = zero(T))
        for (a, p) in enumerate(ps), j = 1:16, d = 1:3
            p.x[d][j] = T(mod(j - 1, 4) + 0.25 + 0.1a) * T(0.5)
        end
        init!(st, ps)
        failed = false
        for _ = 1:3
            before = _composition_snapshot(st, ps)
            try
                step!(st, ps, st.dt)
            catch error
                @test error isa ArgumentError || error isa ErrorException
                @test isequal(before, _composition_snapshot(st, ps))
                @test !st.busy[]
                @test_throws typeof(error) step!(st, ps, st.dt)
                @test isequal(before, _composition_snapshot(st, ps))
                failed = true
                break
            end
        end
        @test failed
    end

    @testset "changing live storage and values reject before publication" begin
        mutations = (
            ps -> pop!(ps[1].v[2]),
            ps -> push!(ps[1].id, UInt64(99)),
            ps -> (ps[1].weight[1] = -1),
            ps -> (ps[1].weight[1] = Inf),
            ps -> (ps[1].x[1][1] = NaN),
            ps -> (ps[1].v[3][1] = Inf),
            ps -> (ps[1].m = 0),
            ps -> (ps[1].q = -1),
            ps -> (ps[2].v = ps[1].v),
        )
        for mutate! in mutations
            st, ps = _composition_fixture(Val(1), Float64, Xoshiro)
            init!(st, ps)
            mutate!(ps)
            before = _composition_snapshot(st, ps)
            @test_throws Union{ArgumentError,DimensionMismatch} step!(st, ps, st.dt)
            @test isequal(before, _composition_snapshot(st, ps))
            @test !st.busy[]
        end
    end

    @testset "sampled failure restores the full stream and retries identically" for T in (
            Float32,
            Float64,
        ),
        R in (Xoshiro, MersenneTwister)

        rejected = 0
        for seed = 1:8
            st, ps = _composition_sampled_failure_fixture(T, R, seed)
            init!(st, ps)
            before = _composition_snapshot(st, ps)
            failure = try
                step!(st, ps, st.dt)
                nothing
            catch error
                error
            end
            if failure === nothing
                @test all(all(isfinite, v) for p in ps for v in p.v)
                continue
            end
            rejected += 1
            @test failure isa ArgumentError
            @test isequal(before, _composition_snapshot(st, ps))
            @test !st.busy[]
            again = try
                step!(st, ps, st.dt)
                nothing
            catch error
                error
            end
            @test typeof(again) === typeof(failure)
            @test sprint(showerror, again) == sprint(showerror, failure)
            @test isequal(before, _composition_snapshot(st, ps))
            direct, dps = _composition_sampled_failure_fixture(T, R, seed)
            init!(direct, dps)
            rng_before = _composition_rng_state(direct.rng)
            cause = try
                _composition_component!(direct, dps)
                nothing
            catch error
                error
            end
            @test typeof(cause) === typeof(failure)
            @test sprint(showerror, cause) == sprint(showerror, failure)
            @test direct.inner.step[] == 1
            @test !isequal(rng_before, _composition_rng_state(direct.rng))
        end
        @test rejected > 0
    end

    @testset "zero logarithms reproduce a bare HybridStepper D=$D T=$T" for D = 1:3,
        T in (Float32, Float64)

        for logs in (zero(T), zeros(T, 2, 2))
            st, ps = _composition_fixture(Val(D), T, Xoshiro; logs)
            # Independent grid and stepper, never wrapped, with the same initial state.
            _, plain_ps = _composition_fixture(Val(D), T, Xoshiro; logs)
            g = FourierGrid(st.grid.n, st.grid.L)
            plain = HybridStepper(g, st.inner.model, st.inner.shape, sum(nparticles, plain_ps))
            foreach(copyto!, plain.fields.B, st.inner.fields.B)
            init!(st, ps)
            init!(plain, plain_ps)
            rng = _composition_rng_state(st.rng)
            for _ = 1:4
                @test step!(st, ps, st.dt) === st
                step!(plain, plain_ps, st.dt)
                @test step!(st, ps, zero(T)) === st
                step!(plain, plain_ps, zero(T))
                @test all(
                    isequal(H._checkpoint_particle_state(p), H._checkpoint_particle_state(q)) for
                    (p, q) in zip(ps, plain_ps)
                )
                @test isequal(
                    H._hybrid_field_state(st.inner.fields),
                    H._hybrid_field_state(plain.fields),
                )
                @test st.inner.fields.floor_count[] == plain.fields.floor_count[]
                @test (st.inner.time[], st.inner.step[], st.inner.last_dt[]) ===
                      (plain.time[], plain.step[], plain.last_dt[])
            end
            @test isequal(_composition_rng_state(st.rng), rng)
        end
    end

    @testset "fixed-step contract after stepping T=$T" for T in (Float32, Float64)
        st, ps = _composition_fixture(Val(2), T, Xoshiro)
        init!(st, ps)
        step!(st, ps, st.dt)
        step!(st, ps, st.dt)
        before = _composition_snapshot(st, ps)
        for dt in (zero(T), -zero(T), 0, 0.0)
            @test step!(st, ps, dt) === st
            @test isequal(_composition_snapshot(st, ps), before)
        end
        for dt in (nextfloat(st.dt), prevfloat(st.dt), 2st.dt, st.dt / 2)
            @test_throws ArgumentError step!(st, ps, dt)
            @test isequal(_composition_snapshot(st, ps), before)
            @test !st.busy[]
        end
        # A wider real equal to the bound step at grid precision is that step.
        @test step!(st, ps, BigFloat(st.dt)) === st
        @test st.source_count == st.inner.step[] == 3
    end

    @testset "restart continues the uninterrupted run R=$R" for R in (Xoshiro, MersenneTwister)
        logs = [0.05 0.02 0.01; 0.02 0.04 0.03; 0.01 0.03 0.06]
        make() = _composition_fixture(Val(2), Float64, R; ns = 3, logs)
        run, rps = make()
        init!(run, rps)
        mktempdir() do dir
            path = joinpath(dir, "boundary-4.bin")
            for k = 1:10
                step!(run, rps, run.dt)
                k == 4 && save_checkpoint(path, run, rps)
            end
            resumed, sps = make()
            @test load_checkpoint!(resumed, sps, path) === resumed
            @test resumed.source_count == resumed.inner.step[] == 4
            for _ = 1:6
                step!(resumed, sps, resumed.dt)
            end
            @test isequal(_composition_snapshot(resumed, sps), _composition_snapshot(run, rps))
            @test rand(resumed.rng, UInt64, 9) == rand(run.rng, UInt64, 9)
            # The seed alone is not a continuation state: same payload, initial stream.
            seed_only, qps = make()
            load_checkpoint!(seed_only, qps, path)
            copy!(seed_only.rng, R(712))
            for _ = 1:6
                step!(seed_only, qps, seed_only.dt)
            end
            @test !isequal(
                [H._checkpoint_particle_state(p) for p in qps],
                [H._checkpoint_particle_state(p) for p in sps],
            )
        end
    end

end
