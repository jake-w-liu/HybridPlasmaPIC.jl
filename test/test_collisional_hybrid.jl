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
        @test st.transaction.Eforce === nothing
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
            HybridModel(CGLElectrons(0.6, 0.3, 1.0, 1.0)),
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

    @testset "ideal evolving pressure uses the canonical split $T $R" for T in (Float32, Float64),
        R in (Xoshiro, MersenneTwister)

        make(logs) = _composition_fixture(Val(1), T, R; logs, closure = EvolvingElectrons(T(5 / 3)))
        st, ps = make(T(0.05))
        reference, rps = make(T(0.05))
        for run in (st, reference)
            run.inner.fields.pe .= T[0.02, 0.022, 0.02, 0.018]
        end
        initial_pressure = copy(st.inner.fields.pe)
        init!(st, ps)
        init!(reference, rps)
        @test st.inner.fields.pe == initial_pressure
        for _ = 1:2
            step!(st, ps, st.dt)
            step!(reference.inner, rps, reference.dt)
            advanced_pressure = copy(reference.inner.fields.pe)
            collide_coulomb!(
                rps,
                reference.grid,
                reference.dt;
                coulomb_log = reference.coulomb_log,
                epsilon0 = reference.epsilon0,
                weight_to_number = reference.weight_to_number,
                u_floor = reference.u_floor,
                rng = reference.rng,
                workspace = reference.workspace,
            )
            @test reference.inner.fields.pe == advanced_pressure
            H._hybrid_integer_moments!(reference.inner, rps)
            ohms_law!(reference.inner.fields, reference.inner.model, reference.grid)
            H._recenter_carried_E!(reference.inner, rps, reference.dt)
            @test reference.inner.fields.pe == advanced_pressure
            @test isequal(_composition_physical(st, ps), _composition_physical(reference, rps))
        end
        @test st.inner.fields.pe != initial_pressure

        inactive, ips = make(zero(T))
        bare, bps = make(zero(T))
        init!(inactive, ips)
        init!(bare.inner, bps)
        rng = _composition_rng_state(inactive.rng)
        for _ = 1:2
            step!(inactive, ips, inactive.dt)
            step!(bare.inner, bps, bare.dt)
            @test isequal(_composition_physical(inactive, ips), _composition_physical(bare, bps))
        end
        @test isequal(_composition_rng_state(inactive.rng), rng)
        inactive.inner.fields.pe[1] = -one(T)
        entry = _composition_snapshot(inactive, ips)
        @test step!(inactive, ips, zero(T)) === inactive
        @test_throws ArgumentError init!(inactive, ips)
        @test_throws ArgumentError step!(inactive, ips, inactive.dt)
        @test_throws ArgumentError H._checkpoint_state(inactive, ips)
        @test isequal(entry, _composition_snapshot(inactive, ips))
        @test !inactive.busy[]

        cold, cps = _composition_fixture(
            Val(1),
            T,
            R;
            temperature = zero(T),
            closure = EvolvingElectrons(T(5 / 3)),
        )
        init!(cold, cps)
        step!(cold, cps, cold.dt)
        @test all(iszero, cold.inner.fields.pe)
        @test cold.inner.pressure_work.Eforce === nothing
        cold.inner.fields.pe = cold.inner.pressure_work.ptmp
        entry = _composition_snapshot(cold, cps)
        @test_throws ArgumentError init!(cold, cps)
        @test isequal(entry, _composition_snapshot(cold, cps))
    end

    @testset "evolving pressure and sampled failure restore the stream" begin
        rejected = 0
        for seed = 1:8
            st, ps = _composition_sampled_failure_fixture(
                Float64,
                Xoshiro,
                seed;
                closure = EvolvingElectrons(),
            )
            init!(st, ps)
            entry = _composition_snapshot(st, ps)
            failure = try
                step!(st, ps, st.dt)
                nothing
            catch error
                error
            end
            failure === nothing && continue
            rejected += 1
            @test failure isa ArgumentError
            @test isequal(entry, _composition_snapshot(st, ps))
            @test !st.busy[]
            @test_throws typeof(failure) step!(st, ps, st.dt)
            @test isequal(entry, _composition_snapshot(st, ps))
            direct, dps = _composition_sampled_failure_fixture(
                Float64,
                Xoshiro,
                seed;
                closure = EvolvingElectrons(),
            )
            init!(direct, dps)
            rng = _composition_rng_state(direct.rng)
            @test_throws typeof(failure) _composition_component!(direct, dps)
            @test !isequal(rng, _composition_rng_state(direct.rng))
        end
        @test rejected > 0
    end

    @testset "rejected pressure stage restores the full composed state" begin
        # Retain the canonical rejection after an accepted earlier subcycle:
        # pressure has advanced before the second subcycle rejects the step.
        g = FourierGrid((32,), (2pi,))
        p = ParticleSet{1,Float64}(128)
        load_lattice_1d!(p, 0.0, g.L[1])
        set_density_weight!(p, 1.0, g)
        @. p.v[1] = sin(p.x[1])
        inner = HybridStepper(g, HybridModel(EvolvingElectrons(3.0)), CIC(), nparticles(p))
        fill!(inner.fields.pe, 1.0)
        st = CollisionalHybridStepper(
            inner,
            [p];
            dt = 0.8,
            rng = Xoshiro(712),
            coulomb_log = 0.05,
            epsilon0 = 1.0,
            weight_to_number = 1.0,
        )
        init!(st, [p])
        entry = _composition_snapshot(st, [p])
        pressure_identity = st.inner.fields.pe
        failure = try
            step!(st, [p], st.dt; NB = 2)
            nothing
        catch error
            error
        end
        @test failure isa ArgumentError
        @test occursin("subcycle 2", sprint(showerror, failure))
        @test isequal(entry, _composition_snapshot(st, [p]))
        @test st.inner.fields.pe === pressure_identity
        @test !st.busy[]
    end

    @testset "transport keeps the carried force and heat schedule $T etaH=$ηH" for T in (
            Float32,
            Float64,
        ),
        ηH in (zero(T), T(0.01))

        st, ps = _composition_transport_fixture(T, Xoshiro; ηH)
        reference, rps = _composition_transport_fixture(T, Xoshiro; ηH)
        init!(st, ps)
        init!(reference, rps)
        @test st.transaction.Eforce !== st.inner.pressure_work.Eforce
        if !iszero(ηH)
            # The initial Bz is a resolved first Fourier mode in a box of length 2.
            # Its current and full Laplacian obey LJ = -pi^2 J independently of heat.
            initial = st.inner.fields
            @test maximum(abs, initial.lapJ[2]) > zero(T)
            @test initial.lapJ[2] ≈ -T(pi)^2 .* initial.J[2] rtol = 32eps(T)
        end
        for _ = 1:2
            step!(st, ps, st.dt)
            pressure = copy(reference.inner.fields.pe)
            _composition_component!(reference, rps)
            @test reference.inner.fields.pe != pressure
            @test isequal(_composition_physical(st, ps), _composition_physical(reference, rps))
        end
        f, force = st.inner.fields, st.inner.pressure_work.Eforce
        ideal = ntuple(_ -> similar(f.n), 3)
        for I in eachindex(f.n)
            ux, uy, uz = (a[I] for a in f.ui)
            bx, by, bz = (a[I] for a in f.B)
            jx, jy, jz = (a[I] for a in f.J)
            ideal[1][I] =
                -(uy * bz - uz * by) + (jy * bz - jz * by) * f.ninv[I] - f.gradp[1][I] * f.ninv[I]
            ideal[2][I] = -(uz * bx - ux * bz) + (jz * bx - jx * bz) * f.ninv[I]
            ideal[3][I] = -(ux * by - uy * bx) + (jx * by - jy * bx) * f.ninv[I]
        end
        for c = 1:3
            @test force[c] ≈ ideal[c] rtol = 8eps(T)
            expected = force[c] .+ T(st.model.η) .* f.J[c]
            iszero(ηH) || (expected .-= ηH .* f.lapJ[c])
            @test isequal(f.E[c], expected)
        end
        @test any(c -> !isequal(f.E[c], force[c]), 1:3)
        saved_force = st.transaction.Eforce
        st.transaction.Eforce = st.inner.pressure_work.Eforce
        entry = _composition_snapshot(st, ps)
        @test_throws ArgumentError step!(st, ps, st.dt)
        @test isequal(entry, _composition_snapshot(st, ps))
        st.transaction.Eforce = saved_force

        unsupported = st.inner.model
        for (bad_etaH, bad_de2) in ((-0.01, 0.0), (Float64(ηH), 0.01))
            st.inner.model = HybridModel(
                unsupported.closure,
                unsupported.η,
                bad_etaH,
                unsupported.nfloor,
                bad_de2,
            )
            @test_throws ArgumentError step!(st, ps, st.dt)
            @test isequal(entry, _composition_snapshot(st, ps))
        end
        st.inner.model = unsupported

        inactive, ips = _composition_transport_fixture(T, Xoshiro; ηH, logs = zero(T))
        bare, bps = _composition_transport_fixture(T, Xoshiro; ηH, logs = zero(T))
        init!(inactive, ips)
        init!(bare.inner, bps)
        step!(inactive, ips, inactive.dt)
        step!(bare.inner, bps, bare.dt)
        @test isequal(_composition_physical(inactive, ips), _composition_physical(bare, bps))

        zero_eta, zps = _composition_transport_fixture(T, Xoshiro; η = zero(T), ηH)
        if iszero(ηH)
            evolving, epspecies = _composition_transport_fixture(
                T,
                Xoshiro;
                η = zero(T),
                closure = EvolvingElectrons(T(5 / 3)),
            )
            init!(zero_eta, zps)
            init!(evolving, epspecies)
            @test evolving.transaction.Eforce === nothing
            for _ = 1:2
                step!(zero_eta, zps, zero_eta.dt)
                step!(evolving, epspecies, evolving.dt)
                a, b =
                    _composition_snapshot(zero_eta, zps), _composition_snapshot(evolving, epspecies)
                @test isequal(a.fields, b.fields)
                @test isequal(a.particles, b.particles)
                @test isequal((a.floor, a.clock, a.rng), (b.floor, b.clock, b.rng))
                @test isequal(a.Eforce, b.fields.E)
            end
        else
            pure_hyper, hps = _composition_transport_fixture(T, Xoshiro; η = zero(T), ηH)
            init!(zero_eta, zps)
            init!(pure_hyper, hps)
            step!(zero_eta, zps, zero_eta.dt)
            _composition_component!(pure_hyper, hps)
            @test isequal(
                _composition_physical(zero_eta, zps),
                _composition_physical(pure_hyper, hps),
            )
        end
    end

    @testset "conducting pressure uses the canonical collision split $T etaH=$ηH kappa=$κ" for T in
                                                                                               (
            Float32,
            Float64,
        ),
        ηH in (zero(T), T(0.01)),
        κ in (zero(T), T(0.07))

        make(logs) = _composition_transport_fixture(
            T,
            Xoshiro;
            ηH,
            logs,
            closure = ConductingElectrons(T(5 / 3); κ),
        )
        st, ps = make(T(0.05))
        reference, rps = make(T(0.05))
        init!(st, ps)
        init!(reference, rps)
        for _ = 1:2
            step!(st, ps, st.dt)
            _composition_component!(reference, rps)
            @test isequal(_composition_physical(st, ps), _composition_physical(reference, rps))
        end
        inactive, ips = make(zero(T))
        bare, bps = make(zero(T))
        init!(inactive, ips)
        init!(bare.inner, bps)
        step!(inactive, ips, inactive.dt)
        step!(bare.inner, bps, bare.dt)
        @test isequal(_composition_physical(inactive, ips), _composition_physical(bare, bps))
        # The outer force rollback cannot alias the live inner carried force.
        snapshot = st.transaction.Eforce
        st.transaction.Eforce = (st.inner.pressure_work.Eforce[1], snapshot[2], snapshot[3])
        before = _composition_snapshot(st, ps)
        @test_throws ArgumentError step!(st, ps, st.dt)
        @test isequal(_composition_snapshot(st, ps), before)
        st.transaction.Eforce = snapshot
    end

    @testset "conduction guard restores the composed state $T" for T in (Float32, Float64)
        st, ps = _composition_transport_fixture(
            T,
            Xoshiro;
            closure = ConductingElectrons(T(2); κ = floatmax(T)),
        )
        init!(st, ps)
        before = _composition_snapshot(st, ps)
        for _ = 1:2
            @test_throws ArgumentError step!(st, ps, st.dt)
            @test isequal(_composition_snapshot(st, ps), before)
        end
    end

    @testset "outer transport transaction survives a completed inner step etaH=$ηH" for ηH in
                                                                                        (0.0, 0.01)

        st, ps = _composition_transport_fixture(Float64, Xoshiro; ηH)
        init!(st, ps)
        entry = _composition_snapshot(st, ps)
        force_identity = st.inner.pressure_work.Eforce
        injected = ErrorException("test-owned exception after completed inner step")
        H._save_collisional_transaction!(st)
        caught = try
            try
                H._step_hybrid!(st.inner, ps, st.dt, 1, st)
                advanced = _composition_snapshot(st, ps)
                @test st.inner.step[] == 1
                @test !isequal(entry.Eforce, advanced.Eforce)
                @test !isequal(entry.fields.pe, advanced.fields.pe)
                @test !isequal(entry.rng, advanced.rng)
                throw(injected)
            finally
                H._restore_collisional_transaction!(st)
            end
        catch error
            error
        end
        @test caught === injected
        @test isequal(entry, _composition_snapshot(st, ps))
        @test st.inner.pressure_work.Eforce === force_identity
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
