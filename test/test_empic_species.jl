module EMPICSpeciesTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC

function fixture(::Type{T}, D; n_sub = 1, shape = CIC()) where {T}
    g = FourierGrid(ntuple(_ -> 6, D), ntuple(_ -> one(T), D))
    species = [ParticleSet{D,T}(k + 2; q = q, m = k + 1) for (k, q) in enumerate((-2, 1, 3, 0))]
    # Exact physical charge: -2*2 + 1*1 + 3*1 = 0. The fourth species is neutral.
    for (k, ps) in enumerate(species)
        for d = 1:D
            ps.x[d] .= [T(mod(p / (nparticles(ps) + 1) + 0.07k * d, 1)) for p = 1:nparticles(ps)]
        end
        fill!(ps.weight, zero(T))
        ps.weight[1] = k == 1 ? T(2) : one(T)
        if k != 3
            for c = 1:3
                fill!(ps.v[c], T(0.01c * k))
            end
        end
    end
    state = EMPIC(g, species; mobile = Bool[true, true, false, true], n_sub, shape)
    return state, species
end

function snapshot(state, species)
    io = IOBuffer()
    serialize(
        io,
        (
            fields = H._empic_transaction_fields(state),
            phase = (state.time[], state.step[], state.last_dt[], state.last_sub_dt[]),
            species = [H._checkpoint_particle_state(ps) for ps in species],
            initialized = state.collection.initialized[],
            owners = objectid.(state.collection.owners),
        ),
    )
    return take!(io)
end

function physical(state, species)
    return deepcopy((
        fields = H._empic_checkpoint_fields(state),
        species = [H._checkpoint_particle_state(ps) for ps in species],
    ))
end

@testset "EM collection charge, summed midpoint current and fixed populations" begin
    for T in (Float32, Float64), D = 1:3, shape in (NGP(), CIC(), TSC())
        state, species = fixture(T, D; shape)
        @test_throws ArgumentError step_empic!(state, species, T(0.01))
        init_empic!(state, species)
        fixed = deepcopy(H._checkpoint_particle_state(species[3]))
        for dt in T[0.01, 0, 0.015]
            before = snapshot(state, species)
            step_empic!(state, species, dt)
            iszero(dt) && @test snapshot(state, species) == before
            @test H._checkpoint_particle_state(species[3]) == fixed
            iszero(dt) && continue
            @test charge_conservation_residual(state, dt) < 1024eps(T)
            # Independent public scalar deposition at each saved midpoint sums
            # q*w*v before one continuity correction; fixed species carry no J.
            for c = 1:3
                expected, scratch = zeros(T, state.g.n), zeros(T, state.g.n)
                for (i, ps) in enumerate(species)
                    state.collection.mobile[i] || continue
                    midpoint = deepcopy(ps)
                    foreach(copyto!, midpoint.x, state.collection.midpoints[i])
                    vals = ps.q .* ps.weight .* ps.v[c]
                    deposit_scalar!(scratch, midpoint, vals, state.g, shape)
                    expected .+= scratch ./ prod(state.g.dx)
                end
                @test isapprox(state.Jraw[c], expected; atol = 128eps(T), rtol = 128eps(T))
            end
        end
        diagnostic = electromagnetic_diagnostics(state, species)
        @test diagnostic.net_charge == 0
        @test diagnostic.kinetic == kinetic_energy(species)
        @test diagnostic.momentum == total_momentum(species)
        @test diagnostic.field == em_field_energy(state)
        @test diagnostic.total == diagnostic.kinetic + diagnostic.field
        @test diagnostic.velocity_time_offsets ==
              T[-state.last_sub_dt[]/2, -state.last_sub_dt[]/2, 0, -state.last_sub_dt[]/2]
        @test diagnostic.magnetic_time_offset == -state.last_sub_dt[] / 2
    end
end

@testset "EM collection current combines all original species before range rounding" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (one(T),))
        species = [ParticleSet{1,T}(1; q = q) for q in (1, -1, 1)]
        for (i, ps) in enumerate(species)
            ps.x[1][1] = T(0.25)
            ps.weight[1] = i == 3 ? one(T) : floatmax(T)
            ps.v[1][1] = i == 3 ? T(0.25) : T(2)
        end
        state = EMPIC(g, species; neutrality = :neutralize, shape = NGP())
        init_empic!(state, species)
        work = H._empic_collection_prepare!(state, species)
        for (ps, mid) in zip(work.staged, work.midpoints)
            foreach(copyto!, mid, ps.x)
        end
        H._empic_collection_current!(work.stage, work)
        # (2*floatmax - 2*floatmax + 1/4) / (1/4) = 1.
        @test work.stage.Jraw[1] == T[0, 1, 0, 0]
        @test all(iszero, work.stage.Jraw[2])
        @test all(iszero, work.stage.Jraw[3])
    end
end

@testset "EM collections reuse legacy two-species and fast-step evolution" begin
    for T in (Float32, Float64), D = 1:3
        state, species = fixture(T, D)
        init_empic!(state, species)
        fine_species = deepcopy(species)
        fine = EMPIC(
            FourierGrid(state.g.n, state.g.L),
            fine_species;
            mobile = Bool[true, true, false, true],
        )
        init_empic!(fine, fine_species)
        for (dt, ns) in ((T(0.03), 3), (T(0.02), 2), (T(0.012), 4))
            state.n_sub = ns
            step_empic!(state, species, dt)
            for _ = 1:ns
                step_empic!(fine, fine_species, dt / ns)
                @test charge_conservation_residual(fine, dt / ns) < 1024eps(T)
            end
            @test physical(state, species) == physical(fine, fine_species)
            @test state.time[] == fine.time[]
            @test state.last_sub_dt[] == fine.last_sub_dt[]
        end
        g = FourierGrid(ntuple(_ -> 6, D), ntuple(_ -> one(T), D))
        e, ions = ParticleSet{D,T}(3; q = -1), ParticleSet{D,T}(3; q = 1, m = 7)
        for ps in (e, ions)
            foreach(a -> (a .= T[0.1, 0.4, 0.8]), ps.x)
            fill!(ps.weight, T(0.25))
            foreach(a -> fill!(a, T(0.01)), ps.v)
        end
        legacy = EMPIC(g, 3; mobile = true, mi = T(7))
        copied = deepcopy([e, ions])
        collection = EMPIC(FourierGrid(g.n, g.L), copied)
        init_empic!(legacy, e, ions)
        init_empic!(collection, copied)
        for dt in T[0.01, 0.02, 0, 0.007]
            step_empic!(legacy, e, ions, dt)
            step_empic!(collection, copied, dt)
            @test physical(legacy, [e, ions]) == physical(collection, copied)
        end
    end
end

@testset "EM collection shared old field and late arithmetic rejection" begin
    for T in (Float32, Float64)
        state, species = fixture(T, 2)
        foreach(ps -> fill!(ps.weight, zero(T)), species)
        init_empic!(state, species)
        initial = deepcopy(species)
        force = T(0.04)
        fill!(state.E[1], force)
        dt = T(0.02)
        step_empic!(state, species, dt)
        for (i, ps) in enumerate(species)
            expected =
                initial[i].v[1] .+
                (state.collection.mobile[i] ? (ps.q / ps.m) * (dt / 2) * force : zero(T))
            @test isapprox(ps.v[1], expected; atol = 32eps(T), rtol = 32eps(T))
        end
        # Earlier populations can already be changed in private priming when the
        # last original-factor kick is unrepresentable. Nothing public commits.
        state, species = fixture(T, 1)
        foreach(ps -> fill!(ps.weight, zero(T)), species)
        species[end].q, species[end].m = floatmax(T), floatmin(T)
        init_empic!(state, species)
        fill!(state.E[1], one(T))
        before = snapshot(state, species)
        @test_throws ArgumentError step_empic!(state, species, one(T))
        @test snapshot(state, species) == before
    end
end

@testset "EM collection policy and preflight preserve complete physical state" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (one(T),))
        empty = ParticleSet{1,T}[]
        state = EMPIC(g, empty; n_sub = 2)
        init_empic!(state, empty)
        step_empic!(state, empty, T(0.02))
        @test electromagnetic_diagnostics(state, empty).total == 0
        ps = ParticleSet{1,T}(1; q = 2, m = 3)
        ps.weight[1] = one(T)
        ps.x[1][1] = T(0.25)
        @test_throws ArgumentError EMPIC(g, [ps])
        neutralized = EMPIC(g, [ps]; neutrality = :neutralize)
        init_empic!(neutralized, [ps])
        @test electromagnetic_diagnostics(neutralized, [ps]).neutralizing_charge == -2
        @test electromagnetic_diagnostics(neutralized, [ps]).net_charge == 0
        @test EMPIC(g, [ps]; relativistic = true, neutrality = :neutralize).relativistic
        @test_throws ArgumentError EMPIC(g, [ps]; mobile = Bool[])
        ps.v[1][1] = one(T)
        @test_throws ArgumentError EMPIC(g, [ps]; mobile = Bool[false], neutrality = :neutralize)
        @test_throws ArgumentError EMPIC(g, [ps, ps]; neutrality = :neutralize)
        state, species = fixture(T, 1)
        init_empic!(state, species)
        before = snapshot(state, species)
        @test_throws ArgumentError step_empic!(state, reverse(species), T(0.01))
        @test snapshot(state, species) == before
        @test_throws ArgumentError step_empic!(state, species[1], T(0.01))
        @test snapshot(state, species) == before
        saved = species[end].m
        species[end].m = -one(T)
        invalid = snapshot(state, species)
        @test_throws ArgumentError step_empic!(state, species, T(0.01))
        @test snapshot(state, species) == invalid
        @test_throws ArgumentError init_empic!(state, species)
        @test snapshot(state, species) == invalid
        species[end].m = saved
        old = state.E
        state.E = (state.g.kfull[1], state.E[2], state.E[3])
        invalid = snapshot(state, species)
        metadata = copy(state.g.kfull[1])
        @test_throws ArgumentError step_empic!(state, species, T(0.01))
        @test snapshot(state, species) == invalid
        @test state.g.kfull[1] == metadata
        state.E = old
        @test snapshot(state, species) == before
        original_grid = state.g
        state.g = FourierGrid(state.g.n, (T(2),))
        invalid = snapshot(state, species)
        @test_throws ArgumentError step_empic!(state, species, T(0.01))
        @test snapshot(state, species) == invalid
        state.g = original_grid
        # Zero weights isolate original q,m priming/kicks from self-fields.
        extreme = ParticleSet{1,T}(1; q = floatmax(T), m = floatmin(T))
        extreme.weight[1] = zero(T)
        extreme.x[1][1] = T(0.25)
        range_state = EMPIC(g, [extreme])
        init_empic!(range_state, [extreme])
        fill!(range_state.E[1], floatmin(T))
        dt = 2floatmin(T)
        step_empic!(range_state, [extreme], dt)
        @test extreme.v[1][1] == floatmax(T) * (dt / 2)
        @test all(isfinite, extreme.v[1])
        # A neutral population can have unrepresentable number-density scratch
        # while its original q-weighted physical source is exactly zero.
        neutral = ParticleSet{1,T}(1; q = 0)
        neutral.x[1][1], neutral.v[1][1] = T(0.25), T(0.1)
        neutral.weight[1] = floatmax(T)
        neutral_state = EMPIC(g, [neutral]; shape = NGP())
        init_empic!(neutral_state, [neutral])
        step_empic!(neutral_state, [neutral], T(0.01))
        @test any(isinf, neutral_state.ne)
        @test all(iszero, neutral_state.rho_np1)
        @test all(a -> all(iszero, a), neutral_state.J)
    end
end

@testset "EM collection checkpoint/archive fast phase and source detachment" begin
    for T in (Float32, Float64)
        state, species = fixture(T, 2; n_sub = 3)
        init_empic!(state, species)
        step_empic!(state, species, T(0.03))
        mktempdir() do dir
            path = joinpath(dir, "em.chk")
            save_checkpoint(path, state, species)
            # Collection and legacy files must reject the opposite destination
            # through public IO, before changing either physical state.
            legacy = EMPIC(FourierGrid(state.g.n, state.g.L), 0)
            electrons = ParticleSet{2,T}(0; q = -one(T))
            legacy_path = joinpath(dir, "legacy.chk")
            save_checkpoint(legacy_path, legacy, electrons)
            collection_before = snapshot(state, species)
            legacy_before = deepcopy(H._checkpoint_state(legacy, electrons))
            @test_throws ArgumentError load_checkpoint!(state, species, legacy_path)
            @test snapshot(state, species) == collection_before
            @test_throws ArgumentError load_checkpoint!(legacy, electrons, path)
            @test H._checkpoint_state(legacy, electrons) == legacy_before
            destination = [ParticleSet{2,T}(0) for _ = 1:4]
            restored = EMPIC(FourierGrid(state.g.n, state.g.L), destination;
                mobile = Bool[true, true, false, true], c = T(2))
            load_checkpoint!(restored, destination, path)
            @test physical(restored, destination) == physical(state, species)
            for (dt, ns) in ((T(0.02), 2), (zero(T), 4), (T(0.012), 3))
                state.n_sub = restored.n_sub = ns
                step_empic!(state, species, dt)
                step_empic!(restored, destination, dt)
                @test physical(restored, destination) == physical(state, species)
                @test (
                    restored.time[],
                    restored.step[],
                    restored.last_dt[],
                    restored.last_sub_dt[],
                ) == (state.time[], state.step[], state.last_dt[], state.last_sub_dt[])
            end
            archive = joinpath(dir, "em.run")
            archive_run(archive, state, species; rng_seed = 19)
            archived = load_archive(archive).state
            H._restore_checkpoint_state!(restored, destination, archived)
            @test physical(restored, destination) == physical(state, species)
            before = snapshot(restored, destination)
            malformed = merge(deepcopy(archived), (last_sub_dt = zero(T),))
            @test_throws ArgumentError H._restore_checkpoint_state!(
                restored,
                destination,
                malformed,
            )
            @test snapshot(restored, destination) == before
            aliased = merge(
                deepcopy(archived),
                (fields = merge(archived.fields, (B = archived.fields.E,)),),
            )
            @test_throws ArgumentError H._restore_checkpoint_state!(restored, destination, aliased)
            @test snapshot(restored, destination) == before
            # All source fields are independent, but are permuted live arrays.
            # Detaching before the first copy preserves both incoming tuples.
            swapped = merge(
                deepcopy(archived),
                (fields = merge(archived.fields, (E = restored.B, B = restored.E)),),
            )
            expected_E, expected_B = deepcopy(restored.B), deepcopy(restored.E)
            H._restore_checkpoint_state!(restored, destination, swapped)
            @test restored.E == expected_E
            @test restored.B == expected_B
        end
    end
end

@testset "EM collection strided destinations retain storage on restart" begin
    state, species = fixture(Float64, 1)
    init_empic!(state, species)
    step_empic!(state, species, 0.01)
    function strided(a)
        storage = zeros(eltype(a), 2length(a) + 1)
        result = view(storage, 2:2:2length(a))
        copyto!(result, a)
        return result
    end
    views = ParticleSet{1,Float64}[
        ParticleSet{1,Float64}(
            map(strided, ps.x),
            map(strided, ps.v),
            strided(ps.weight),
            strided(ps.id),
            strided(ps.tag),
            ps.q,
            ps.m,
        ) for ps in species
    ]
    restored =
        EMPIC(FourierGrid(state.g.n, state.g.L), views; mobile = Bool[true, true, false, true])
    identities = [ps.x[1] for ps in views]
    mktempdir() do dir
        path = joinpath(dir, "views.chk")
        save_checkpoint(path, state, species)
        load_checkpoint!(restored, views, path)
        @test all(ps.x[1] === a for (ps, a) in zip(views, identities))
        step_empic!(state, species, 0.02)
        step_empic!(restored, views, 0.02)
        @test physical(restored, views) == physical(state, species)
        payload = deepcopy(H._checkpoint_state(state, species))
        records = Any[payload.species...]
        last = records[end]
        records[end] = merge(
            last,
            (
                x = map(a -> a[1:(end-1)], last.x),
                v = map(a -> a[1:(end-1)], last.v),
                weight = last.weight[1:(end-1)],
                id = last.id[1:(end-1)],
                tag = last.tag[1:(end-1)],
            ),
        )
        malformed = merge(payload, (species = records,))
        before = snapshot(restored, views)
        @test_throws ArgumentError H._restore_checkpoint_state!(restored, views, malformed)
        @test snapshot(restored, views) == before
    end
end
end # module
