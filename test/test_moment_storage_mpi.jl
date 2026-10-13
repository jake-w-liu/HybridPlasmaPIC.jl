using HybridPlasmaPIC, Test
import MPI

@testset "MPI moment buffers validate collectively before deposition" begin
    ensure_mpi_initialized!()
    ranks = MPI.Comm_size(MPI.COMM_WORLD)
    for T in (Float32, Float64), D = 1:3
        layout =
            LogicalRankLayout(ntuple(d -> d == 1 ? ranks : 1, D); periodic = ntuple(_ -> true, D))
        ctx = create_cartesian_communicator(layout; reorder = false)
        try
            g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> one(T), D))
            count = prod(g.n)
            for collection in (false, true),
                failure in (
                    :particle_work,
                    :output_work,
                    :outputs,
                    :output_particle,
                    :late_shape,
                    :work_length,
                    :position,
                )

                ps = ParticleSet{D,T}(count)
                foreach(a -> fill!(a, T(0.2)), ps.x)
                foreach(a -> fill!(a, T(0.3)), ps.v)
                ps.weight .= T(0.4 / ranks)
                population = collection ? [ps] : ps
                f = HybridFields{D,T}(g.n)
                fill!(f.n, T(7))
                foreach(a -> fill!(a, T(8)), f.ui)
                work = fill(T(9), count)
                if ctx.mpi_rank == ranks - 1
                    if failure === :particle_work
                        work = ps.weight
                    elseif failure === :output_work
                        work = vec(f.n)
                    elseif failure === :outputs
                        f.ui = (f.ui[1], f.ui[2], f.n)
                    elseif failure === :output_particle
                        f.ui = (f.ui[1], f.ui[2], reshape(ps.weight, g.n))
                    elseif failure === :late_shape
                        f.ui = (f.ui[1], f.ui[2], zeros(T, ntuple(_ -> 2, D)))
                    elseif failure === :work_length
                        resize!(work, count - 1)
                    else
                        ps.x[D][end] = T(NaN)
                    end
                end
                arrays = (f.n, f.ui..., work, ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
                before = map(copy, arrays)
                @test_throws Union{ArgumentError,DimensionMismatch} mpi_compute_moments!(
                    f,
                    population,
                    g,
                    CIC(),
                    T(0.01),
                    ctx;
                    work,
                )
                @test all(isequal(a, b) for (a, b) in zip(arrays, before))
            end
            # A failed collective leaves the communicator usable for valid
            # scalar and collection calls, including globally empty particles.
            for collection in (false, true), empty in (false, true)
                ps = ParticleSet{D,T}(empty ? 0 : count)
                foreach(a -> fill!(a, T(0.2)), ps.x)
                for c = 1:3
                    fill!(ps.v[c], T(c))
                end
                ps.weight .= T(0.4 / ranks)
                f = HybridFields{D,T}(g.n)
                population = collection ? [ps] : ps
                mpi_compute_moments!(
                    f,
                    population,
                    g,
                    CIC(),
                    T(0.01),
                    ctx;
                    work = zeros(T, nparticles(ps)),
                )
                total_weight = empty ? zero(T) : T(0.4 / ranks) * T(ranks * count)
                tolerance = 256eps(T)
                @test sum(f.n) * prod(g.dx) ≈ total_weight rtol = tolerance atol = tolerance
                for c = 1:3
                    @test sum(f.ui[c] .* max.(f.n, T(0.01))) * prod(g.dx) ≈ T(c) * total_weight rtol =
                        tolerance atol = tolerance
                end
            end
        finally
            free_mpi_communicator!(ctx)
        end
    end
end
