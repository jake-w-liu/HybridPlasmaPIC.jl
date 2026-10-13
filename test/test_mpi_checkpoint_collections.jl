using HybridPlasmaPIC, Test, Serialization
import MPI

function mpi_restart_species(g::FourierGrid{D,T}, ctx) where {D,T}
    return map(1:3) do s
        count = s == 3 || isodd(ctx.mpi_rank + s) ? 0 : 2
        ps = ParticleSet{D,T}(count; q = T(s), m = T(s^2))
        for d = 1:D
            ps.x[d] .= T(d == 1 ? (ctx.mpi_rank + 0.25) / ctx.mpi_size : 0.25) * g.L[d]
        end
        if count > 0
            # The last marker crosses its owning rank boundary on the first step.
            ps.x[1][end] = T((ctx.mpi_rank + 0.999999) / ctx.mpi_size) * g.L[1]
        end
        for c = 1:3
            ps.v[c] .= T(0.01c)
        end
        ps.weight .= T(0.2 / ctx.mpi_size)
        ps.id .= UInt64(100s + 10ctx.mpi_rank) .+ UInt64.(1:count)
        ps.tag .= UInt32(s)
        ps
    end
end

mpi_restart_views(ps::ParticleSet{D,T}) where {D,T} = ParticleSet{D,T}(
    map(a -> view(a, :), ps.x),
    map(a -> view(a, :), ps.v),
    view(ps.weight, :),
    view(ps.id, :),
    view(ps.tag, :),
    ps.q,
    ps.m,
)
mpi_restart_state(st, particles) = deepcopy(HybridPlasmaPIC._checkpoint_state(st, particles))
function mpi_restart_bytes(st, particles)
    io = IOBuffer()
    serialize(
        io,
        (
            HybridPlasmaPIC._checkpoint_state(st, particles),
            st.Ep,
            st.Bp,
            st.vpred,
            st.xmid,
            st.work,
            st.fields,
            st.fn,
            st.fui,
            st.Escr,
            st.k1,
            st.k2,
            st.k3,
            st.k4,
            st.Btmp,
        ),
    )
    return take!(io)
end
function mpi_restart_error(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

@testset "MPI ordered collection checkpoint continuation" begin
    ensure_mpi_initialized!()
    world = MPI.Comm_size(MPI.COMM_WORLD)
    for T in (Float32, Float64), D = 1:3
        ctx = create_cartesian_communicator(
            LogicalRankLayout(ntuple(d -> d == 1 ? world : 1, D));
            reorder = false,
        )
        dir = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
        try
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
            model = HybridModel(IsothermalElectrons(T(0.02)); nfloor = T(0.1))
            for shape in (NGP(), CIC(), TSC())
                source = HybridStepper(g, model, shape, 0)
                species = mpi_restart_species(g, ctx)
                fill!(source.fields.B[3], one(T))
                mpi_init!(source, species, ctx)
                counts = nparticles.(species)
                identifiers = [sum(ps.id; init = UInt64(0)) for ps in species]
                global_counts = MPI.Allreduce(counts, +, ctx.comm)
                global_ids = MPI.Allreduce(identifiers, +, ctx.comm)
                mpi_step!(source, species, ctx, T(0.001); NB = 2)
                @test MPI.Allreduce(nparticles.(species), +, ctx.comm) == global_counts
                @test MPI.Allreduce(
                    [sum(ps.id; init = UInt64(0)) for ps in species],
                    +,
                    ctx.comm,
                ) == global_ids
                if world > 1
                    @test MPI.Allreduce(Int(counts != nparticles.(species)), max, ctx.comm) == 1
                end
                path = save_mpi_checkpoint(dir, source, species, ctx)
                manifest = deserialize(path)
                @test manifest.schema == 2
                @test isdir(joinpath(dir, manifest.generation))
                @test readdir(dir) == sort([basename(path), manifest.generation])
                for views in (false, true)
                    expected, expected_species = deepcopy((source, species))
                    target = HybridStepper(g, model, shape, 1)
                    destination =
                        views ? map(mpi_restart_views, deepcopy(species)) :
                        [ParticleSet{D,T}(1) for _ in species]
                    identities =
                        [(ps.x..., ps.v..., ps.weight, ps.id, ps.tag) for ps in destination]
                    @test load_mpi_checkpoint!(target, destination, dir, ctx) === target
                    @test isequal(
                        mpi_restart_state(target, destination),
                        mpi_restart_state(source, species),
                    )
                    @test length(target.work) == sum(nparticles, destination)
                    for (dt, NB) in ((T(0.0007), 3), (T(0.0012), 1))
                        mpi_step!(
                            expected,
                            expected_species,
                            ctx,
                            dt;
                            NB,
                            migrate_particles = !views,
                        )
                        mpi_step!(target, destination, ctx, dt; NB, migrate_particles = !views)
                        @test isequal(
                            mpi_restart_state(target, destination),
                            mpi_restart_state(expected, expected_species),
                        )
                    end
                    if views
                        @test all(
                            a === b for (ps, arrays) in zip(destination, identities) for
                            (a, b) in zip((ps.x..., ps.v..., ps.weight, ps.id, ps.tag), arrays)
                        )
                    end
                end
                # Global emptiness retains the ordered representation and zero workspace.
                empty_species = ParticleSet{D,T}[]
                empty_st = HybridStepper(g, model, shape, 0)
                mpi_init!(empty_st, empty_species, ctx)
                save_mpi_checkpoint(dir, empty_st, empty_species, ctx)
                restored = HybridStepper(g, model, shape, 2)
                @test load_mpi_checkpoint!(restored, ParticleSet{D,T}[], dir, ctx) === restored
                @test isempty(restored.work)
                @test !ispath(joinpath(dir, manifest.generation))
            end
        finally
            MPI.Barrier(ctx.comm)
            ctx.mpi_rank == 0 && rm(dir; recursive = true)
            MPI.Barrier(ctx.comm)
            free_mpi_communicator!(ctx)
        end
    end
end

@testset "MPI collection restore rejects a later invalid population collectively" begin
    ensure_mpi_initialized!()
    world = MPI.Comm_size(MPI.COMM_WORLD)
    ctx = create_cartesian_communicator(LogicalRankLayout((world,)); reorder = false)
    dir = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
    try
        g = FourierGrid((4,), (2π,))
        model = HybridModel(IsothermalElectrons(0.02))
        st = HybridStepper(g, model, CIC(), 0)
        species = mpi_restart_species(g, ctx)
        path = save_mpi_checkpoint(dir, st, species, ctx)
        manifest = deserialize(path)
        source_file = joinpath(dir, manifest.generation, manifest.rank_files[ctx.logical_rank])
        for failure in (
            :count,
            :alias,
            :later_resize,
            :later_payload,
            :clock,
            :species_order,
            :manifest,
            :manifest_agreement,
        )
            mktempdir() do local_dir
                cp(path, joinpath(local_dir, basename(path)))
                mkdir(joinpath(local_dir, manifest.generation))
                local_file =
                    joinpath(local_dir, manifest.generation, manifest.rank_files[ctx.logical_rank])
                cp(source_file, local_file)
                target = HybridStepper(g, model, CIC(), 1)
                destination = ParticleSet{1,Float64}[ParticleSet{1,Float64}(1) for _ in species]
                failing = ctx.mpi_rank == world - 1
                if failing
                    if failure === :count
                        pop!(destination)
                    elseif failure === :alias
                        destination[3] = destination[2]
                    elseif failure === :later_resize
                        destination[3] = mpi_restart_views(destination[3])
                    elseif failure === :later_payload
                        saved = deserialize(local_file)
                        saved.species[3] =
                            merge(saved.species[3], (v = (zeros(1), zeros(0), zeros(0)),))
                        serialize(local_file, saved)
                    elseif failure === :clock
                        saved = deserialize(local_file)
                        serialize(local_file, merge(saved, (last_dt = NaN,)))
                    elseif failure === :species_order
                        saved = deserialize(local_file)
                        saved.species[3] = merge(saved.species[3], (q = 4.0,))
                        serialize(local_file, saved)
                    elseif failure === :manifest_agreement
                        alternative = "_hybrid_mpi_alternate"
                        Base.Filesystem.rename(
                            joinpath(local_dir, manifest.generation),
                            joinpath(local_dir, alternative),
                        )
                        serialize(
                            joinpath(local_dir, basename(path)),
                            merge(manifest, (generation = alternative,)),
                        )
                    else
                        serialize(
                            joinpath(local_dir, basename(path)),
                            merge(manifest, (generation = "..",)),
                        )
                    end
                end
                before = mpi_restart_bytes(target, destination)
                err = mpi_restart_error(
                    () -> load_mpi_checkpoint!(target, destination, local_dir, ctx),
                )
                # Species agreement differs only when there is another rank.
                if failure in (:species_order, :manifest_agreement) && world == 1
                    @test err === nothing
                else
                    @test err isa Exception
                    @test mpi_restart_bytes(target, destination) == before
                end
            end
        end
        if world > 1
            before_manifest = read(path)
            before_rank = read(source_file)
            for failure in (:count, :charge, :layout)
                inputs = deepcopy(species)
                if ctx.mpi_rank == world - 1
                    failure === :count && pop!(inputs)
                    failure === :charge && (inputs[end].q = 4.0)
                    failure === :layout && (inputs = inputs[1])
                end
                err = mpi_restart_error(() -> save_mpi_checkpoint(dir, st, inputs, ctx))
                @test err isa ArgumentError
                @test read(path) == before_manifest
                @test read(source_file) == before_rank
                @test readdir(dir) == sort([basename(path), manifest.generation])
            end
        end
        target = HybridStepper(g, model, CIC(), 0)
        destination = [ParticleSet{1,Float64}(0) for _ in species]
        @test load_mpi_checkpoint!(target, destination, dir, ctx) === target
        @test isequal(mpi_restart_state(target, destination), mpi_restart_state(st, species))
    finally
        MPI.Barrier(ctx.comm)
        ctx.mpi_rank == 0 && rm(dir; recursive = true)
        MPI.Barrier(ctx.comm)
        free_mpi_communicator!(ctx)
    end
end

mutable struct MPIFailingCheckpointClosure <: HybridPlasmaPIC.ElectronClosure
    fail::Bool
    obstruction_directory::String
end
MPIFailingCheckpointClosure(fail::Bool) = MPIFailingCheckpointClosure(fail, "")
function Serialization.serialize(
    s::Serialization.AbstractSerializer,
    c::MPIFailingCheckpointClosure,
)
    if c.fail
        if !isempty(c.obstruction_directory)
            dir = c.obstruction_directory
            active = deserialize(joinpath(dir, "mpi_checkpoint_manifest.ser")).generation
            pending = only(
                filter(name -> startswith(name, "_hybrid_mpi_") && name != active, readdir(dir)),
            )
            write(joinpath(dir, pending, "unrelated.txt"), "retain failed-generation obstruction")
        end
        error("injected rank serializer failure")
    end
    invoke(Serialization.serialize, Tuple{Serialization.AbstractSerializer,Any}, s, c)
end

@testset "MPI checkpoint publication preserves the prior generation on failure" begin
    ensure_mpi_initialized!()
    world = MPI.Comm_size(MPI.COMM_WORLD)
    ctx = create_cartesian_communicator(LogicalRankLayout((world,)); reorder = false)
    dir = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
    try
        g = FourierGrid((4,), (2π,))
        closure = MPIFailingCheckpointClosure(false)
        st = HybridStepper(g, HybridModel(closure), CIC(), 1)
        ps = ParticleSet{1,Float64}(1)
        ps.x[1] .= 0.2
        path = save_mpi_checkpoint(dir, st, ps, ctx)
        manifest = deserialize(path)
        saved_manifest = read(path)
        saved_rank = read(joinpath(dir, manifest.generation, manifest.rank_files[ctx.logical_rank]))
        for failure in (:clock, :serializer, :cleanup, :directory, :phase)
            local_dir = dir
            failing = ctx.mpi_rank == world - 1
            if failing
                failure === :clock && (st.time[] = NaN)
                failure in (:serializer, :cleanup) && (closure.fail = true)
                failure === :cleanup && (closure.obstruction_directory = dir)
                failure === :directory && (local_dir = joinpath(dir, "missing"))
                failure === :phase && (st.time[] = 1.0)
            end
            err = mpi_restart_error(() -> save_mpi_checkpoint(local_dir, st, ps, ctx))
            st.time[] = 0.0
            closure.fail = false
            closure.obstruction_directory = ""
            if world == 1 && failure in (:directory, :phase)
                @test err === nothing # A different valid directory/origin is valid on COMM_SELF.
                failure === :directory && rm(local_dir; recursive = true)
                save_mpi_checkpoint(dir, st, ps, ctx)
                manifest = deserialize(path)
                saved_manifest = read(path)
                saved_rank = read(joinpath(dir, manifest.generation, only(manifest.rank_files)))
            else
                @test err isa Exception
                @test read(path) == saved_manifest
                @test read(
                    joinpath(dir, manifest.generation, manifest.rank_files[ctx.logical_rank]),
                ) == saved_rank
                if failure === :cleanup
                    @test err isa CompositeException
                    @test length(err.exceptions) == 2
                    @test all(e -> e isa CapturedException, err.exceptions)
                    @test occursin(
                        failing ? "injected rank serializer failure" : "rank write",
                        sprint(showerror, err.exceptions[1]),
                    )
                    pending = only(setdiff(readdir(dir), [basename(path), manifest.generation]))
                    @test readdir(joinpath(dir, pending)) == ["unrelated.txt"]
                    @test read(joinpath(dir, pending, "unrelated.txt"), String) ==
                          "retain failed-generation obstruction"
                    MPI.Barrier(ctx.comm)
                    if ctx.mpi_rank == 0
                        rm(joinpath(dir, pending, "unrelated.txt"))
                        rm(joinpath(dir, pending))
                    end
                    MPI.Barrier(ctx.comm)
                end
                @test readdir(dir) == sort([basename(path), manifest.generation])
            end
        end
        # A generation path redirected outside the checkpoint directory is not
        # owned storage. Replacement must not follow it when removing old files.
        external = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
        generation_path = joinpath(dir, manifest.generation)
        external_generation = joinpath(external, "saved")
        if ctx.mpi_rank == 0
            Base.Filesystem.rename(generation_path, external_generation)
            symlink(external_generation, generation_path; dir_target = true)
        end
        MPI.Barrier(ctx.comm)
        err = mpi_restart_error(() -> save_mpi_checkpoint(dir, st, ps, ctx))
        @test err isa (ctx.mpi_rank == 0 ? ArgumentError : ErrorException)
        @test occursin(
            ctx.mpi_rank == 0 ? "symbolic link" : "generation setup",
            sprint(showerror, err),
        )
        @test islink(generation_path)
        @test read(path) == saved_manifest
        @test read(joinpath(external_generation, manifest.rank_files[ctx.logical_rank])) ==
              saved_rank
        @test readdir(dir) == sort([basename(path), manifest.generation])
        MPI.Barrier(ctx.comm)
        if ctx.mpi_rank == 0
            rm(generation_path)
            Base.Filesystem.rename(external_generation, generation_path)
            rm(external)
        end
        MPI.Barrier(ctx.comm)
        save_mpi_checkpoint(dir, st, ps, ctx)
        manifest = deserialize(path)
        @test !ispath(generation_path)
        # Genuine schema 1 layout: no generation field and flat logical-rank files.
        legacy = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
        try
            record = merge(
                HybridPlasmaPIC._checkpoint_state(st, ps),
                (
                    mpi_schema = 1,
                    logical_rank = ctx.logical_rank,
                    mpi_rank = ctx.mpi_rank,
                    mpi_size = ctx.mpi_size,
                    layout_ranks = ctx.layout.ranks,
                    layout_periodic = ctx.layout.periodic,
                ),
            )
            serialize(joinpath(legacy, manifest.rank_files[ctx.logical_rank]), record)
            MPI.Barrier(ctx.comm)
            if ctx.mpi_rank == 0
                old = (
                    schema = 1,
                    D = 1,
                    T = Float64,
                    ncell = g.n,
                    L = g.L,
                    ranks = ctx.layout.ranks,
                    periodic = ctx.layout.periodic,
                    nranks = ctx.mpi_size,
                    time = 0.0,
                    step = 0,
                    rank_files = manifest.rank_files,
                )
                serialize(joinpath(legacy, basename(path)), old)
            end
            MPI.Barrier(ctx.comm)
            target = HybridStepper(g, HybridModel(MPIFailingCheckpointClosure(false)), CIC(), 0)
            destination = ParticleSet{1,Float64}(0)
            @test load_mpi_checkpoint!(target, destination, legacy, ctx) === target
            @test destination.x == ps.x
            save_mpi_checkpoint(legacy, st, ps, ctx)
            @test !isfile(joinpath(legacy, manifest.rank_files[ctx.logical_rank]))
            @test load_mpi_checkpoint!(target, destination, legacy, ctx) === target
        finally
            MPI.Barrier(ctx.comm)
            ctx.mpi_rank == 0 && rm(legacy; recursive = true)
            MPI.Barrier(ctx.comm)
        end
        # A post-publication cleanup error is explicit; unrelated files survive
        # and the newly published checkpoint remains loadable.
        old_generation = manifest.generation
        sentinel = joinpath(dir, old_generation, "unrelated.txt")
        ctx.mpi_rank == 0 && write(sentinel, "retain this file")
        MPI.Barrier(ctx.comm)
        st.time[] = 1.0
        cleanup_error = mpi_restart_error(() -> save_mpi_checkpoint(dir, st, ps, ctx))
        @test cleanup_error isa CompositeException
        @test length(cleanup_error.exceptions) == 2
        @test occursin("was published", sprint(showerror, cleanup_error.exceptions[1]))
        @test read(sentinel, String) == "retain this file"
        manifest = deserialize(path)
        @test manifest.generation != old_generation
        target = HybridStepper(g, HybridModel(MPIFailingCheckpointClosure(false)), CIC(), 0)
        destination = ParticleSet{1,Float64}(0)
        @test load_mpi_checkpoint!(target, destination, dir, ctx) === target
        @test target.time[] == 1.0
        MPI.Barrier(ctx.comm)
        if ctx.mpi_rank == 0
            rm(sentinel)
            rm(joinpath(dir, old_generation))
        end
        MPI.Barrier(ctx.comm)
        # Repeated successful saves leave exactly the active generation and manifest.
        for origin in (2.0, 3.0)
            st.time[] = origin
            previous = manifest.generation
            save_mpi_checkpoint(dir, st, ps, ctx)
            manifest = deserialize(path)
            @test !ispath(joinpath(dir, previous))
            @test readdir(dir) == sort([basename(path), manifest.generation])
        end
    finally
        MPI.Barrier(ctx.comm)
        ctx.mpi_rank == 0 && rm(dir; recursive = true)
        MPI.Barrier(ctx.comm)
        free_mpi_communicator!(ctx)
    end
end
