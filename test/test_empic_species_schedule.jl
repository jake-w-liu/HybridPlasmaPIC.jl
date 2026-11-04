module EMPICSpeciesScheduleTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC
Q(x) = Rational{BigInt}(x)

function schedule_fixture(
    ::Type{T},
    D;
    strides,
    n_sub = 6,
    shape = NGP(),
    relativistic = false,
) where {T}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
    species = [
        ParticleSet{D,T}(2; q = -2, m = 3),
        ParticleSet{D,T}(3; q = 1, m = 5),
        ParticleSet{D,T}(2; q = 3, m = 2),
    ]
    for (k, ps) in enumerate(species)
        for d = 1:D
            ps.x[d] .= [T(mod(p / (nparticles(ps) + 2) + 0.05k * d, 1)) for p = 1:nparticles(ps)]
        end
        fill!(ps.weight, zero(T))
        if k < 3
            for c = 1:3
                ps.v[c] .= [T(0.01c * k + 0.002p) for p = 1:nparticles(ps)]
            end
        end
    end
    es = EMPIC(
        g,
        species;
        mobile = Bool[true, true, false],
        species_stride = strides,
        n_sub,
        shape,
        relativistic,
    )
    return es, species
end

function snapshot(es, species)
    schedule =
        es.collection.schedule === nothing ? nothing :
        (
            es.n_sub,
            copy(es.collection.schedule.stride),
            copy(es.collection.schedule.last_species_dt),
        )
    io = IOBuffer()
    serialize(
        io,
        (
            fields = H._empic_transaction_fields(es),
            phase = (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]),
            species = [H._checkpoint_particle_state(ps) for ps in species],
            initialized = es.collection.initialized[],
            owners = objectid.(es.collection.owners),
            schedule = schedule,
        ),
    )
    return take!(io)
end

physical(es, species) = deepcopy((
    fields = H._empic_checkpoint_fields(es),
    species = [H._checkpoint_particle_state(ps) for ps in species],
))

cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

# Independent rational three-component Boris push on constant fields; no
# production gather, mover or derivative enters the reference.
function boris_reference(v, alpha, duration, E, B)
    h = alpha * duration / 2
    u = v .+ h .* E
    t = h .* B
    s = (2 / (1 + sum(abs2, t))) .* t
    vp = u .+ cross3(u, t)
    return u .+ cross3(vp, s) .+ h .* E
end

# Zero-weight populations freeze the uniform seeded fields (J ≡ 0, ρ ≡ 0, the
# curls vanish), so every mobile species follows an exact rational schedule:
# a -H_s/2 backward-Euler prime at step zero, due kicks at fast indices
# (index-1) % r == 0 with the previous/new half-joined species interval, and a
# fast drift at every substep.
@testset "Generic dividing events versus independent rational oracle" begin
    # Precision coverage at D = 1; the schedule code is type- and dimension-
    # generic, so Float64 alone carries the D = 2, 3 sweeps.
    for T in (Float32, Float64), D = 1:3
        (T === Float32 && D > 1) && continue
        es, species = schedule_fixture(T, D; strides = [2, 3, 0], n_sub = 6)
        init_empic!(es, species)
        Ev, Bv = (T(1 / 8), T(-1 / 16), T(1 / 32)), (T(1 / 16), T(1 / 8), T(-1 / 64))
        for c = 1:3
            fill!(es.E[c], Ev[c])
            fill!(es.B[c], Bv[c])
        end
        E, B = Q.(Ev), Q.(Bv)
        vstate = [[Tuple(Q(ps.v[c][p]) for c = 1:3) for p = 1:nparticles(ps)] for ps in species]
        xstate = [[Tuple(Q(ps.x[d][p]) for d = 1:D) for p = 1:nparticles(ps)] for ps in species]
        mobile = es.collection.mobile
        old_H = zeros(Rational{BigInt}, length(species))
        for (dt, M, strides) in ((T(1 / 32), 6, [2, 3, 0]), (T(1 / 64), 4, [4, 2, 0]))
            set_empic_schedule!(es; n_sub = M, species_stride = strides)
            h = Q(dt) / M
            Hs = [r == 0 ? zero(Rational{BigInt}) : Q(dt) / (M ÷ r) for r in strides]
            if es.step[] == 0
                for (i, ps) in enumerate(species)
                    mobile[i] || continue
                    alpha = Q(ps.q) / Q(ps.m)
                    for p = 1:nparticles(ps)
                        v0 = vstate[i][p]
                        vstate[i][p] = v0 .- Q(Hs[i]) / 2 .* alpha .* (E .+ cross3(v0, B))
                    end
                end
                old_H .= Hs
            end
            for j = 0:(M-1), (i, ps) in enumerate(species)
                mobile[i] || continue
                alpha = Q(ps.q) / Q(ps.m)
                if j % strides[i] == 0
                    duration = (old_H[i] + Hs[i]) / 2
                    for p = 1:nparticles(ps)
                        vstate[i][p] = boris_reference(vstate[i][p], alpha, duration, E, B)
                    end
                    old_H[i] = Hs[i]
                end
                for p = 1:nparticles(ps), d = 1:D
                    xstate[i][p] =
                        Base.setindex(xstate[i][p], xstate[i][p][d] + h * vstate[i][p][d], d)
                end
            end
            step_empic!(es, species, dt)
            for (i, ps) in enumerate(species)
                mobile[i] || continue
                for p = 1:nparticles(ps), d = 1:D
                    @test isapprox(
                        mod(ps.x[d][p], one(T)),
                        T(mod(xstate[i][p][d], Q(1)));
                        rtol = 128eps(T),
                    )
                end
                for p = 1:nparticles(ps), c = 1:3
                    @test isapprox(ps.v[c][p], T(vstate[i][p][c]); rtol = 128eps(T))
                end
            end
            fixed = species[.!mobile]
            @test all(ps -> all(a -> all(iszero, a), ps.v), fixed)
            @test es.collection.schedule.last_species_dt ==
                  [mobile[i] ? T(Hs[i]) : zero(T) for i in eachindex(species)]
            d = electromagnetic_diagnostics(es, species)
            @test d.species_stride == strides
            @test d.last_species_dt == es.collection.schedule.last_species_dt
            @test d.velocity_time_offsets ==
                  [mobile[i] ? -T(Hs[i]) / 2 : zero(T) for i in eachindex(species)]
            @test d.magnetic_time_offset == -T(dt) / M / 2
            before = snapshot(es, species)
            step_empic!(es, species, zero(T))
            @test snapshot(es, species) == before
        end
        # A fixed population neither drifts nor records a species interval.
        @test all(iszero, es.collection.schedule.last_species_dt[.!mobile])
    end
end

@testset "All-fast opt-in is bitwise identical to the unscheduled collection" begin
    for T in (Float32, Float64), D in (1, 2)
        g = FourierGrid(ntuple(_ -> 6, D), ntuple(_ -> one(T), D))
        mobile = Bool[true, true, false, true]
        function populations()
            species = [ParticleSet{D,T}(k + 2; q, m = k + 1) for (k, q) in enumerate((-2, 1, 3, 0))]
            for (k, ps) in enumerate(species)
                for d = 1:D
                    ps.x[d] .=
                        [T(mod(p / (nparticles(ps) + 1) + 0.07k * d, 1)) for p = 1:nparticles(ps)]
                end
                fill!(ps.weight, zero(T))
                ps.weight[1] = k == 1 ? T(2) : one(T)
                if k != 3
                    for c = 1:3
                        fill!(ps.v[c], T(0.01c * k))
                    end
                end
            end
            return species
        end
        scheduled_species, plain_species = populations(), populations()
        scheduled = EMPIC(g, scheduled_species; mobile, n_sub = 4, species_stride = [1, 1, 0, 1])
        plain = EMPIC(FourierGrid(g.n, g.L), plain_species; mobile, n_sub = 4)
        for state in (scheduled, plain)
            for c = 1:3
                fill!(state.E[c], T(0.02c))
                fill!(state.B[c], T(-0.01c))
            end
        end
        init_empic!(scheduled, scheduled_species)
        init_empic!(plain, plain_species)
        @test H._checkpoint_state(scheduled, scheduled_species).schema_version === 2
        @test H._checkpoint_state(plain, plain_species).schema_version === 1
        for dt in (T(0.02), T(0.012), T(0.03))
            step_empic!(scheduled, scheduled_species, dt)
            step_empic!(plain, plain_species, dt)
            @test physical(scheduled, scheduled_species) == physical(plain, plain_species)
            @test (
                scheduled.time[],
                scheduled.step[],
                scheduled.last_dt[],
                scheduled.last_sub_dt[],
            ) == (plain.time[], plain.step[], plain.last_dt[], plain.last_sub_dt[])
            ds, dp = electromagnetic_diagnostics(scheduled, scheduled_species),
            electromagnetic_diagnostics(plain, plain_species)
            @test all(
                getproperty(ds, k) == getproperty(dp, k) for
                k in keys(dp) if !(k in (:species_stride, :last_species_dt))
            )
            @test ds.last_species_dt == T[
                scheduled.last_sub_dt[],
                scheduled.last_sub_dt[],
                zero(T),
                scheduled.last_sub_dt[],
            ]
        end
    end
end

@testset "Schedule setter validation and phase guards" begin
    for T in (Float32, Float64)
        es, species = schedule_fixture(T, 1; strides = [1, 1, 0], n_sub = 4)
        proposal = [2, 1, 0]
        set_empic_schedule!(es; species_stride = proposal)
        proposal[1] = 4
        @test es.collection.schedule.stride == [2, 1, 0]
        for bad in ([1, 3, 0], [1, 2, 1], Bool[1, 1, 0], Int[], [1, -1, 0])
            before = snapshot(es, species)
            @test_throws ArgumentError set_empic_schedule!(es; species_stride = bad)
            @test snapshot(es, species) == before
        end
        @test_throws ArgumentError set_empic_schedule!(es; n_sub = 0)
        init_empic!(es, species)
        step_empic!(es, species, T(1 / 32))
        history = copy(es.collection.schedule.last_species_dt)
        set_empic_schedule!(es; n_sub = 2, species_stride = [2, 1, 0])
        @test es.collection.schedule.last_species_dt == history
        es.n_sub = 3
        before = snapshot(es, species)
        @test_throws ArgumentError step_empic!(es, species, T(1 / 64))
        @test snapshot(es, species) == before
        set_empic_schedule!(es; n_sub = 2, species_stride = [2, 1, 0])
        es.collection.schedule.last_species_dt[1] = prevfloat(es.last_sub_dt[])
        before = snapshot(es, species)
        @test_throws ArgumentError step_empic!(es, species, T(1 / 64))
        @test snapshot(es, species) == before
        es.collection.schedule.last_species_dt .= history
        # Aliased schedule metadata rejects before materialization.
        es.collection.schedule.previous_species_dt = es.collection.schedule.last_species_dt
        before = snapshot(es, species)
        @test_throws ArgumentError step_empic!(es, species, T(1 / 64))
        @test snapshot(es, species) == before
        neutral_ps = ParticleSet{1,T}(1; q = 1)
        fill!(neutral_ps.weight, zero(T))
        unscheduled = EMPIC(FourierGrid((4,), (one(T),)), [neutral_ps]; mobile = [true])
        @test_throws ArgumentError set_empic_schedule!(unscheduled; species_stride = [1])
        for bad in ([3, 1, 0], [1, 2, 3], Int[])
            @test_throws ArgumentError EMPIC(
                FourierGrid((4,), (one(T),)),
                species;
                mobile = Bool[true, true, false],
                species_stride = bad,
                n_sub = 4,
            )
        end
        empty = ParticleSet{1,T}[]
        vacuum = EMPIC(FourierGrid((1,), (one(T),)), empty; species_stride = Int[], n_sub = 3)
        init_empic!(vacuum, empty)
        step_empic!(vacuum, empty, T(1 / 32))
        @test isempty(electromagnetic_diagnostics(vacuum, empty).last_species_dt)
    end
end

@testset "Schema-2 checkpoint/archive round-trip and malformed rejection" begin
    for T in (Float32, Float64)
        es, species = schedule_fixture(T, 2; strides = [2, 1, 0], n_sub = 4)
        init_empic!(es, species)
        step_empic!(es, species, T(0.02))
        saved = deepcopy(H._checkpoint_state(es, species))
        @test saved.schema_version === 2 && saved.schedule_policy === 1
        @test saved.species_stride == [2, 1, 0]
        @test saved.last_species_dt == es.collection.schedule.last_species_dt
        mktempdir() do dir
            path = joinpath(dir, "em.chk")
            save_checkpoint(path, es, species)
            destination = [ParticleSet{2,T}(0) for _ = 1:3]
            restored = EMPIC(
                FourierGrid(es.g.n, es.g.L),
                destination;
                mobile = Bool[true, true, false],
                species_stride = [1, 1, 0],
                n_sub = 1,
                shape = NGP(),
            )
            load_checkpoint!(restored, destination, path)
            @test physical(restored, destination) == physical(es, species)
            @test restored.collection.schedule.stride == [2, 1, 0]
            @test restored.n_sub == 4
            for dt in (T(0.015), zero(T))
                step_empic!(es, species, dt)
                step_empic!(restored, destination, dt)
                @test physical(restored, destination) == physical(es, species)
                @test restored.collection.schedule.last_species_dt ==
                      es.collection.schedule.last_species_dt
            end
            archive = joinpath(dir, "em.run")
            archive_run(archive, es, species; rng_seed = 19)
            archived = load_archive(archive).state
            @test archived.schema_version === 2 && archived.schedule_policy === 1
            H._restore_checkpoint_state!(restored, destination, archived)
            @test physical(restored, destination) == physical(es, species)
        end
        destination = [ParticleSet{2,T}(0) for _ = 1:3]
        target = EMPIC(
            FourierGrid(es.g.n, es.g.L),
            destination;
            mobile = Bool[true, true, false],
            species_stride = [1, 1, 0],
            n_sub = 1,
            shape = NGP(),
        )
        H._restore_checkpoint_state!(target, destination, saved)
        function scheduled_before()
            return deepcopy((
                physical(target, destination),
                target.collection.schedule.stride,
                target.collection.schedule.last_species_dt,
            ))
        end
        for changes in (
            (schedule_policy = true,),
            (schedule_policy = 2,),
            (species_stride = [1, 3, 0],),
            (species_stride = [2, 1, 2],),
            (last_species_dt = fill(zero(T), 3),),
            (schema_version = 1,),
        )
            before = scheduled_before()
            @test_throws ArgumentError H._restore_checkpoint_state!(
                target,
                destination,
                merge(deepcopy(saved), changes),
            )
            @test scheduled_before() == before
        end
        unscheduled = EMPIC(
            FourierGrid(es.g.n, es.g.L),
            [ParticleSet{2,T}(0) for _ = 1:3];
            mobile = Bool[true, true, false],
            shape = NGP(),
        )
        @test_throws ArgumentError H._restore_checkpoint_state!(
            unscheduled,
            [ParticleSet{2,T}(0) for _ = 1:3],
            saved,
        )
    end
end

@testset "Schema-1 payload migrates an opted-in destination to all-fast rows" begin
    for T in (Float32, Float64)
        for stepped in (false, true)
            grid = FourierGrid((6,), (one(T),))
            # Neutral total charge: -1*3 + 1*3 + 0*2 = 0.
            source_species = [
                ParticleSet{1,T}(n; q, m = k + 1) for
                (k, (n, q)) in enumerate(((3, -1), (3, 1), (2, 0)))
            ]
            for (k, ps) in enumerate(source_species)
                ps.x[1] .= [T(mod(p / 5 + 0.11k, 1)) for p = 1:nparticles(ps)]
                fill!(ps.weight, one(T))
                k == 3 || fill!(ps.v[1], T(0.03k))
            end
            source = EMPIC(grid, source_species; mobile = Bool[true, true, false], n_sub = 4)
            init_empic!(source, source_species)
            stepped && step_empic!(source, source_species, T(0.02))
            saved = deepcopy(H._checkpoint_state(source, source_species))
            @test saved.schema_version === 1
            destination = [ParticleSet{1,T}(0) for _ = 1:3]
            scheduled = EMPIC(
                FourierGrid(grid.n, grid.L),
                destination;
                mobile = Bool[true, true, false],
                species_stride = [4, 2, 0],
                n_sub = 8,
            )
            H._restore_checkpoint_state!(scheduled, destination, saved)
            @test scheduled.collection.schedule.stride == [1, 1, 0]
            expected_history =
                stepped ? T[saved.last_sub_dt, saved.last_sub_dt, zero(T)] : zeros(T, 3)
            @test scheduled.collection.schedule.last_species_dt == expected_history
            plain_destination = [ParticleSet{1,T}(0) for _ = 1:3]
            plain = EMPIC(
                FourierGrid(grid.n, grid.L),
                plain_destination;
                mobile = Bool[true, true, false],
                n_sub = 4,
            )
            H._restore_checkpoint_state!(plain, plain_destination, saved)
            step_empic!(scheduled, destination, T(0.015))
            step_empic!(plain, plain_destination, T(0.015))
            @test physical(scheduled, destination) == physical(plain, plain_destination)
            @test scheduled.collection.schedule.last_species_dt ==
                  T[scheduled.last_sub_dt[], scheduled.last_sub_dt[], zero(T)]
        end
    end
end

# Uniform cold E/J mode: omega^2=sum(n*q^2/m). A one-cell grid removes spatial
# conditioning; dyadic coarse/fast intervals isolate temporal error.
@testset "Scheduled cold-mode physical refinement" begin
    for mass in (1.0, 7.0)
        errors = Float64[]
        for blocks in (4, 8, 16)
            species =
                [ParticleSet{1,Float64}(1; q = -1), ParticleSet{1,Float64}(1; q = 1, m = mass)]
            for ps in species
                ps.x[1][1], ps.weight[1] = 0.25, 1.0
            end
            es = EMPIC(FourierGrid((1,), (1.0,)), species; n_sub = 4, species_stride = [2, 4])
            init_empic!(es, species)
            es.E[1][1] = 1 / 16
            dt = 1 / (8blocks)
            for _ = 1:blocks
                step_empic!(es, species, dt)
            end
            omega = sqrt(1 + 1 / mass)
            error = abs(es.E[1][1] - (1 / 16) * cos(omega * es.time[]))
            for (i, ps) in enumerate(species)
                phase = es.time[] - es.collection.schedule.last_species_dt[i] / 2
                expected = (ps.q / ps.m) * (1 / 16) * sin(omega * phase) / omega
                error = max(error, abs(ps.v[1][1] - expected))
            end
            @test error > 128eps(Float64)
            push!(errors, error)
        end
        @info "Generic EMPIC dividing schedule cold refinement" mass errors observed_orders =
            log2.(errors[1:2] ./ errors[2:3])
        @test errors[3] < errors[2] < errors[1]
    end
end

@testset "Relativistic scheduled failure rollback" begin
    for T in (Float32, Float64)
        species = [
            ParticleSet{1,T}(1; m) for
            m in (one(T), reinterpret(T, T === Float32 ? UInt32(1) : UInt64(1)))
        ]
        for ps in species
            ps.x[1][1], ps.weight[1] = zero(T), zero(T)
        end
        es = EMPIC(
            FourierGrid((4,), (T(4),)),
            species;
            relativistic = true,
            c = one(T),
            n_sub = 2,
            species_stride = [1, 2],
        )
        init_empic!(es, species)
        fill!(es.E[1], one(T))
        before = snapshot(es, species)
        @test_throws ArgumentError step_empic!(es, species, one(T))
        @test snapshot(es, species) == before
        species[2].m = one(T)
        set_empic_schedule!(es; n_sub = 4, species_stride = [2, 4])
        step_empic!(es, species, T(1 / 8))
        @test es.collection.schedule.last_species_dt == T[1/16, 1/8]
    end
end
end # module
