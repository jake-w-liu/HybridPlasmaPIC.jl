using HybridPlasmaPIC, Test, Serialization
import MPI

function evolving_mpi_fixture(::Type{T}, collection, ctx) where {T}
    g = FourierGrid((16,), (T(2π),))
    full = ParticleSet{1,T}(32)
    load_lattice_1d!(full, zero(T), g.L[1])
    set_density_weight!(full, one(T), g)
    full.v[1] .= T(0.2)
    full.id .= UInt64.(1:32)
    full.tag .= UInt32(19)
    indices = ctx.mpi_rank == 0 ? collect(1:32) : Int[]
    local_ps = ParticleSet{1,T}(
        map(a -> a[indices], full.x),
        map(a -> a[indices], full.v),
        full.weight[indices],
        full.id[indices],
        full.tag[indices],
        full.q,
        full.m,
    )
    populations = collection ? [local_ps, ParticleSet{1,T}(0; q = T(2), m = T(4))] : local_ps
    originals = collection ? [full, ParticleSet{1,T}(0; q = T(2), m = T(4))] : full
    st = HybridStepper(g, HybridModel(EvolvingElectrons(T(5 / 3))), CIC(), 0)
    fill!(st.fields.B[1], T(0.5))
    for i = 1:g.n[1]
        st.fields.pe[i] = T(0.2) * (one(T) + T(0.1) * sin(T(i - 1) * g.dx[1]))
    end
    return st, populations, originals
end

function evolving_mpi_public_state(st, particles)
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

function evolving_mpi_compare(st, particles, reference, originals, ctx, migrated)
    T = eltype(st.fields.pe)
    # Parallel sums can change roundoff; the bound resolves the 0.02 pressure mode.
    tolerance = T(256) * eps(T)
    for name in (:n, :ui, :E, :B, :J, :pe, :gradp, :ninv)
        actual, expected = getproperty(st.fields, name), getproperty(reference.fields, name)
        aa = actual isa Tuple ? actual : (actual,)
        ee = expected isa Tuple ? expected : (expected,)
        @test all(isapprox(a, e; rtol = tolerance, atol = tolerance) for (a, e) in zip(aa, ee))
    end
    local_pop = particles isa ParticleSet ? (particles,) : particles
    full_pop = originals isa ParticleSet ? (originals,) : originals
    for (ps, full) in zip(local_pop, full_pop)
        @test MPI.Allreduce(nparticles(ps), +, ctx.comm) == nparticles(full)
        @test MPI.Allreduce(sum(ps.id; init = UInt64(0)), +, ctx.comm) ==
              sum(full.id; init = UInt64(0))
        for (p, id) in enumerate(ps.id)
            index = Int(id)
            @test ps.weight[p] == full.weight[index]
            @test ps.tag[p] == full.tag[index]
            @test isapprox(ps.x[1][p], full.x[1][index]; rtol = tolerance, atol = tolerance)
            @test all(
                isapprox(ps.v[c][p], full.v[c][index]; rtol = tolerance, atol = tolerance) for
                c = 1:3
            )
            migrated && @test rank_of_position((ps.x[1][p],), st.g, ctx.layout) == ctx.logical_rank
        end
    end
end

@testset "evolving pressure MPI equivalence and restart" begin
    ensure_mpi_initialized!()
    world = MPI.Comm_size(MPI.COMM_WORLD)
    ctx = create_cartesian_communicator(LogicalRankLayout((world,)); reorder = false)
    try
        for T in (Float32, Float64), collection in (false, true), migrate in (false, true)
            st, particles, originals = evolving_mpi_fixture(T, collection, ctx)
            reference = deepcopy(st)
            init!(reference, originals)
            initial = copy(st.fields.pe)
            mpi_init!(st, particles, ctx)
            @test st.fields.pe == initial
            for (dt, NB) in ((T(0.001), 1), (T(0.002), 2))
                mpi_step!(st, particles, ctx, dt; NB, migrate_particles = migrate)
                step!(reference, originals, dt; NB)
                evolving_mpi_compare(st, particles, reference, originals, ctx, migrate)
            end
            @test st.fields.pe != initial
            dir = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
            try
                path = save_mpi_checkpoint(dir, st, particles, ctx)
                manifest = deserialize(path)
                rankpath = joinpath(dir, manifest.generation, manifest.rank_files[ctx.logical_rank])
                @test deserialize(rankpath).pe == st.fields.pe
                target, restored, _ = evolving_mpi_fixture(T, collection, ctx)
                @test load_mpi_checkpoint!(target, restored, dir, ctx) === target
                @test isequal(
                    HybridPlasmaPIC._checkpoint_state(target, restored),
                    HybridPlasmaPIC._checkpoint_state(st, particles),
                )
                mpi_step!(target, restored, ctx, T(0.0015); NB = 3, migrate_particles = migrate)
                mpi_step!(st, particles, ctx, T(0.0015); NB = 3, migrate_particles = migrate)
                @test isequal(
                    evolving_mpi_public_state(target, restored),
                    evolving_mpi_public_state(st, particles),
                )
                # One corrupt rank must prevent restore on every rank.
                if ctx.mpi_rank == world - 1
                    saved = deserialize(rankpath)
                    serialize(rankpath, merge(saved, (pe = fill(T(NaN), size(saved.pe)),)))
                end
                MPI.Barrier(ctx.comm)
                before = evolving_mpi_public_state(target, restored)
                @test_throws Union{ArgumentError,ErrorException} load_mpi_checkpoint!(
                    target,
                    restored,
                    dir,
                    ctx,
                )
                @test isequal(evolving_mpi_public_state(target, restored), before)
            finally
                MPI.Barrier(ctx.comm)
                ctx.mpi_rank == 0 && rm(dir; recursive = true)
                MPI.Barrier(ctx.comm)
            end
        end
    finally
        free_mpi_communicator!(ctx)
    end
end

@testset "evolving pressure MPI collective rejection and rollback" begin
    ensure_mpi_initialized!()
    world = MPI.Comm_size(MPI.COMM_WORLD)
    ctx = create_cartesian_communicator(LogicalRankLayout((world,)); reorder = false)
    try
        for T in (Float32, Float64), collection in (false, true)
            for invalid in (:negative, :nan, :shape, :alias, :model, :xlength, :vlength)
                st, particles, _ = evolving_mpi_fixture(T, collection, ctx)
                # Rank zero owns particles; shortening an empty rank cannot
                # exercise the transaction's particle-copy boundary.
                invalid_rank = invalid in (:xlength, :vlength) ? 0 : world - 1
                if ctx.mpi_rank == invalid_rank
                    if invalid == :shape
                        st.fields.pe = zeros(T, 15)
                    elseif invalid == :alias
                        st.fields.pe = st.fields.B[1]
                    elseif invalid == :model
                        st.model = HybridModel(st.model.closure, 0.1, 0.0, st.model.nfloor, 0.0)
                    elseif invalid in (:xlength, :vlength)
                        population = collection ? first(particles) : particles
                        component = invalid == :xlength ? population.x : population.v
                        pop!(component[1])
                    else
                        st.fields.pe[1] = invalid == :negative ? -one(T) : T(NaN)
                    end
                end
                before = evolving_mpi_public_state(st, particles)
                local_shape_error =
                    invalid in (:shape, :xlength, :vlength) && ctx.mpi_rank == invalid_rank
                error_type = local_shape_error ? DimensionMismatch : ArgumentError
                @test_throws error_type mpi_init!(st, particles, ctx)
                @test isequal(evolving_mpi_public_state(st, particles), before)
                @test_throws error_type mpi_step!(st, particles, ctx, T(0.001))
                @test isequal(evolving_mpi_public_state(st, particles), before)
                @test mpi_step!(st, particles, ctx, zero(T)) === st
                @test isequal(evolving_mpi_public_state(st, particles), before)
            end
            st, particles, _ = evolving_mpi_fixture(T, collection, ctx)
            if ctx.mpi_rank == world - 1
                # The third-harmonic pressure derivative exceeds floatmax at x=0.
                # Every stored pressure is finite and nonnegative at entry.
                for i = 1:st.g.n[1]
                    phase = T(3) * T(i - 1) * st.g.dx[1]
                    st.fields.pe[i] = (floatmax(T) / T(2)) * (one(T) + sin(phase))
                end
            end
            before = evolving_mpi_public_state(st, particles)
            @test_throws ArgumentError mpi_init!(st, particles, ctx)
            @test isequal(evolving_mpi_public_state(st, particles), before)
            st, particles, _ = evolving_mpi_fixture(T, collection, ctx)
            mpi_init!(st, particles, ctx)
            if ctx.mpi_rank == world - 1
                fill!(st.fields.pe, T(1e-8))
                st.fields.pe[1] = one(T)
            end
            before = evolving_mpi_public_state(st, particles)
            pressure_identity = st.fields.pe
            field_identities = map(name -> getfield(st.fields, name), fieldnames(typeof(st.fields)))
            expected, untouched = deepcopy((st, particles))
            failure = try
                mpi_step!(st, particles, ctx, T(0.01); NB = 2)
                nothing
            catch err
                err
            end
            @test failure isa ArgumentError
            @test occursin(
                ctx.mpi_rank == world - 1 ? "stage" : "another MPI rank",
                sprint(showerror, failure),
            )
            @test st.fields.pe === pressure_identity
            @test all(
                a === getfield(st.fields, name) for
                (a, name) in zip(field_identities, fieldnames(typeof(st.fields)))
            )
            @test isequal(evolving_mpi_public_state(st, particles), before)
            mpi_step!(st, particles, ctx, T(1e-10); NB = 2, migrate_particles = false)
            mpi_step!(expected, untouched, ctx, T(1e-10); NB = 2, migrate_particles = false)
            @test isequal(
                evolving_mpi_public_state(st, particles),
                evolving_mpi_public_state(expected, untouched),
            )
        end
    finally
        free_mpi_communicator!(ctx)
    end
end
