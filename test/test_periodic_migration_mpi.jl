module PeriodicMigrationMPITests
using HybridPlasmaPIC, Test
import MPI

ensure_mpi_initialized!()
const COMM = MPI.COMM_WORLD
const SIZE = MPI.Comm_size(COMM)

@testset "MPI periodic seam routing size=$SIZE" begin
    for T in (Float32, Float64), D = 1:3, periodic in (false, true)
        layout = LogicalRankLayout(
            ntuple(d -> d == 1 ? SIZE : 1, D);
            periodic = ntuple(d -> d == 1 ? periodic : true, D),
        )
        ctx = create_cartesian_communicator(layout; comm = COMM)
        try
            grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
            for x in (-eps(T) / T(4), -nextfloat(zero(T)))
                ps = ParticleSet{D,T}(1)
                for d = 1:D
                    ps.x[d][1] = d == 1 ? x : T(0.25)
                end
                ps.id[1] = UInt64(ctx.logical_rank)
                ps.tag[1] = UInt32(ctx.logical_rank + 10)
                ps.weight[1] = T(ctx.logical_rank) / T(2)
                for c = 1:3
                    ps.v[c][1] = T(ctx.logical_rank + c)
                end
                result = mpi_migrate_particles!(ps, grid, ctx)
                @test result.moved == (periodic ? SIZE - 1 : 0)
                @test result.lost == (periodic ? 0 : SIZE)
                @test result.sent == (periodic && ctx.logical_rank != 1 ? 1 : 0)
                @test result.received == (periodic && ctx.logical_rank == 1 ? SIZE - 1 : 0)
                @test nparticles(ps) == (periodic && ctx.logical_rank == 1 ? SIZE : 0)
                if periodic && ctx.logical_rank == 1
                    @test sort(ps.id) == UInt64.(1:SIZE)
                    for p in eachindex(ps.id)
                        id = Int(ps.id[p])
                        @test ps.x[1][p] == zero(T)
                        @test all(d -> ps.x[d][p] == T(0.25), 2:D)
                        @test ps.weight[p] == T(id) / T(2)
                        @test ps.tag[p] == UInt32(id + 10)
                        @test all(c -> ps.v[c][p] == T(id + c), 1:3)
                    end
                end
                @test mpi_migrate_particles!(ps, grid, ctx) ==
                      (moved = 0, lost = 0, sent = 0, received = 0)
            end
        finally
            free_mpi_communicator!(ctx)
        end
    end
end
end
