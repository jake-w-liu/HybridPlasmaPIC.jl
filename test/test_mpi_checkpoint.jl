using HybridPlasmaPIC, Test
import MPI
import Serialization

function mpi_checkpoint_values(st, ps)
    fields = map(fieldnames(typeof(st.fields))) do name
        value = getfield(st.fields, name)
        value isa Ref ? value[] : value
    end
    return deepcopy((
        st.model,
        st.shape,
        fields,
        st.fn,
        st.fui,
        st.Escr,
        st.k1,
        st.k2,
        st.k3,
        st.k4,
        st.Btmp,
        st.Ep,
        st.Bp,
        st.vpred,
        st.xmid,
        st.work,
        st.time[],
        st.step[],
        st.last_dt[],
        ps.x,
        ps.v,
        ps.weight,
        ps.id,
        ps.tag,
        ps.q,
        ps.m,
    ))
end

mpi_checkpoint_buffers(st, ps) = (
    st.Ep...,
    st.Bp...,
    st.vpred...,
    st.xmid...,
    st.work,
    ps.x...,
    ps.v...,
    ps.weight,
    ps.id,
    ps.tag,
)

function mpi_checkpoint_persistent(st, ps)
    return deepcopy((
        st.model,
        st.shape,
        st.fields.B,
        st.fields.E,
        st.time[],
        st.step[],
        st.last_dt[],
        ps.x,
        ps.v,
        ps.weight,
        ps.id,
        ps.tag,
        ps.q,
        ps.m,
    ))
end

function mpi_checkpoint_views(ps::ParticleSet{D,T}) where {D,T}
    return ParticleSet{D,T}(
        map(a -> view(a, :), ps.x),
        map(a -> view(a, :), ps.v),
        view(ps.weight, :),
        view(ps.id, :),
        view(ps.tag, :),
        ps.q,
        ps.m,
    )
end

@testset "MPI checkpoint restore validates all ranks before mutation" begin
    ensure_mpi_initialized!()
    size = MPI.Comm_size(MPI.COMM_WORLD)
    for T in (Float32, Float64), D = 1:3
        layout =
            LogicalRankLayout(ntuple(d -> d == 1 ? size : 1, D); periodic = ntuple(_ -> true, D))
        ctx = create_cartesian_communicator(layout; reorder = false)
        directory = MPI.bcast(ctx.mpi_rank == 0 ? mktempdir() : "", 0, ctx.comm)
        try
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
            model = HybridModel(IsothermalElectrons(T(0.01)); nfloor = T(0.1))
            source = HybridStepper(g, model, CIC(), 2)
            particles = ParticleSet{D,T}(2)
            for d = 1:D
                particles.x[d] .= T.((0.2, 0.6)) .* g.L[d]
            end
            for c = 1:3
                particles.v[c] .= T(0.01c) .* T.((1, -1))
            end
            particles.weight .= T(0.25 / size)
            particles.id .= UInt64(2ctx.mpi_rank) .+ UInt64.((1, 2))
            particles.tag .= UInt32(ctx.mpi_rank + 1)
            fill!(source.fields.B[3], one(T))
            mpi_init!(source, particles, ctx)
            mpi_step!(source, particles, ctx, T(0.001); NB = 2, migrate_particles = false)
            manifest_path = save_mpi_checkpoint(directory, source, particles, ctx)
            manifest = Serialization.deserialize(manifest_path)
            filename = manifest.rank_files[ctx.logical_rank]
            generation = joinpath(directory, manifest.generation)
            saved = Serialization.deserialize(joinpath(generation, filename))

            # Separate local copies let exactly one rank observe each failure,
            # including a manifest-read failure; the healthy ranks must not restore.
            for failure in (
                :destination,
                :manifest_missing,
                :manifest_schema,
                :rank_missing,
                :rank_truncated,
                :rank_identity,
                :rank_fields,
            )
                mktempdir() do local_directory
                    local_manifest = joinpath(local_directory, basename(manifest_path))
                    local_generation = joinpath(local_directory, manifest.generation)
                    mkdir(local_generation)
                    local_rank = joinpath(local_generation, filename)
                    cp(manifest_path, local_manifest)
                    cp(joinpath(generation, filename), local_rank)
                    target = HybridStepper(g, model, CIC(), 1)
                    destination = ParticleSet{D,T}(1; q = T(2), m = T(3))
                    fill!(target.fields.B[2], T(7))
                    target.time[] = T(4)
                    target.step[] = 5
                    target.last_dt[] = T(0.002)
                    failing = ctx.mpi_rank == size - 1
                    if failing
                        if failure === :destination
                            destination = mpi_checkpoint_views(destination)
                        elseif failure === :manifest_missing
                            rm(local_manifest)
                        elseif failure === :manifest_schema
                            Serialization.serialize(local_manifest, merge(manifest, (schema = -1,)))
                        elseif failure === :rank_missing
                            rm(local_rank)
                        elseif failure === :rank_truncated
                            write(local_rank, UInt8[])
                        elseif failure === :rank_identity
                            Serialization.serialize(local_rank, merge(saved, (logical_rank = 0,)))
                        else
                            Serialization.serialize(
                                local_rank,
                                merge(saved, (B = (zeros(T, 1), saved.B[2], saved.B[3]),)),
                            )
                        end
                    end
                    before = mpi_checkpoint_values(target, destination)
                    buffers = mpi_checkpoint_buffers(target, destination)
                    error = try
                        load_mpi_checkpoint!(target, destination, local_directory, ctx)
                        nothing
                    catch err
                        err
                    end
                    @test error isa Exception
                    if failing
                        expected =
                            failure in (:manifest_missing, :rank_missing) ? SystemError :
                            failure === :rank_truncated ? EOFError : ArgumentError
                        @test error isa expected
                    else
                        @test error isa ErrorException
                        @test occursin(
                            "restore preflight failed on another rank",
                            sprint(showerror, error),
                        )
                    end
                    @test isequal(mpi_checkpoint_values(target, destination), before)
                    @test all(
                        a === b for
                        (a, b) in zip(buffers, mpi_checkpoint_buffers(target, destination))
                    )
                end
            end

            # A failed restore must leave the communicator usable. Check both
            # resizing destinations and equal-sized views, then variable-step continuation.
            for use_views in (false, true)
                expected_st, expected_ps = deepcopy((source, particles))
                target = HybridStepper(g, model, CIC(), use_views ? 2 : 1)
                destination = ParticleSet{D,T}(use_views ? 2 : 1)
                use_views && (destination = mpi_checkpoint_views(destination))
                @test load_mpi_checkpoint!(target, destination, directory, ctx) === target
                @test isequal(
                    mpi_checkpoint_persistent(target, destination),
                    mpi_checkpoint_persistent(source, particles),
                )
                for (dt, NB) in ((T(0.0007), 3), (T(0.0012), 1))
                    mpi_step!(expected_st, expected_ps, ctx, dt; NB, migrate_particles = false)
                    mpi_step!(target, destination, ctx, dt; NB, migrate_particles = false)
                    @test isequal(
                        mpi_checkpoint_persistent(target, destination),
                        mpi_checkpoint_persistent(expected_st, expected_ps),
                    )
                end
            end
        finally
            MPI.Barrier(ctx.comm)
            ctx.mpi_rank == 0 && rm(directory; recursive = true)
            MPI.Barrier(ctx.comm)
            free_mpi_communicator!(ctx)
        end
    end
end
