using HybridPlasmaPIC, Test
import MPI

HybridPlasmaPIC.ensure_mpi_initialized!()

const WORLD = MPI.COMM_WORLD
const WORLD_SIZE = MPI.Comm_size(WORLD)
const WORLD_RANK = MPI.Comm_rank(WORLD)

function _multirank_layout(n::Integer)
    n == 2 && return LogicalRankLayout((2,); periodic = (false,))
    n == 4 && return LogicalRankLayout((2, 2); periodic = (true, false))
    n == 8 && return LogicalRankLayout((2, 2, 2); periodic = (true, true, true))
    return nothing
end

function _cart_rank_for_logical(ctx, logical_rank)
    logical_rank === nothing && return MPI.PROC_NULL
    coords = rank_coords(ctx.layout, logical_rank)
    coords0 = [coords[d] - 1 for d = 1:length(coords)]
    return MPI.Cart_rank(ctx.comm, coords0)
end

function _particle_partition(rank0::Integer, nranks::Integer; total::Integer = 24)
    gids = [gid for gid = 1:total if mod(gid - 1, nranks) == rank0]
    ps = ParticleSet{1,Float64}(length(gids))
    for (i, gid) in pairs(gids)
        gidf = Float64(gid)
        ps.weight[i] = 1.0 + gidf / 32.0
        ps.v[1][i] = -0.25 + gidf / 16.0
        ps.v[2][i] = 0.5 - gidf / 48.0
        ps.v[3][i] = -0.125 + gidf / 64.0
        ps.id[i] = UInt64(gid)
    end
    return ps
end

function _particle_diagnostics(ps::ParticleSet{D,Float64}) where {D}
    px, py, pz = total_momentum(ps)
    kinetic = kinetic_energy(ps)
    return (number = sum(ps.weight), momentum = (px, py, pz), kinetic = kinetic)
end

function _budget_reference_state(::Val{D}; total::Integer = 120) where {D}
    g = FourierGrid(ntuple(d -> 16 + 2d, D), ntuple(d -> 9.0 + d, D))
    ps = ParticleSet{D,Float64}(total; q = 1.0, m = 1.25)
    denom = total + 17
    for p = 1:total
        for d = 1:D
            raw = mod(p * (2d + 1) + 7d, denom)
            ps.x[d][p] = ((raw + 0.25) / denom) * g.L[d]
        end
        pf = Float64(p)
        ps.v[1][p] = (isodd(p) ? 0.83 : -0.47) + 0.001 * pf
        ps.v[2][p] = (mod(p, 3) == 0 ? -0.61 : 0.28) + 0.0005 * pf
        ps.v[3][p] = (mod(p, 4) < 2 ? 0.39 : -0.52) - 0.0003 * pf
        ps.weight[p] = 0.75 + 0.01 * mod(p, 11)
        ps.id[p] = UInt64(10_000 + p)
        ps.tag[p] = UInt32(20_000 + p)
    end
    return g, ps
end

function _copy_particle_subset(ps::ParticleSet{D,T}, idx::AbstractVector{<:Integer}) where {D,T}
    out = ParticleSet{D,T}(length(idx); q = ps.q, m = ps.m)
    for d = 1:D
        out.x[d] .= ps.x[d][idx]
    end
    for c = 1:3
        out.v[c] .= ps.v[c][idx]
    end
    out.weight .= ps.weight[idx]
    out.id .= ps.id[idx]
    out.tag .= ps.tag[idx]
    return out
end

function _rank_budget_particles(
    serial::ParticleSet{D,T},
    g::FourierGrid{D,T},
    layout::LogicalRankLayout{D},
    rank::Integer,
) where {D,T}
    idx = Int[]
    for p = 1:nparticles(serial)
        pos = ntuple(d -> serial.x[d][p], D)
        rank_of_position(pos, g, layout) == rank && push!(idx, p)
    end
    return _copy_particle_subset(serial, idx)
end

function _advance_ballistic_periodic!(ps::ParticleSet{D,T}, g::FourierGrid{D,T}, dt) where {D,T}
    dtT = T(dt)
    @inbounds for p = 1:nparticles(ps), d = 1:D
        ps.x[d][p] += dtT * ps.v[d][p]
    end
    apply_periodic!(ps, ntuple(_ -> zero(T), D), g.L)
    return ps
end

_agreement_grid(::Val{1}) = FourierGrid((16,), (2π,))
_agreement_grid(::Val{2}) = FourierGrid((8, 6), (2π, 2π))
_agreement_grid(::Val{3}) = FourierGrid((6, 6, 6), (2π, 2π, 2π))

function _field_coupled_reference_state(::Val{D}) where {D}
    g = _agreement_grid(Val(D))
    ps = ParticleSet{D,Float64}(prod(g.n); q = 1.0, m = 1.0)
    ΔV = prod(g.dx)
    for (p, I) in enumerate(CartesianIndices(g.n))
        coords = ntuple(d -> (I[d] - 0.5) * g.dx[d], D)
        for d = 1:D
            ps.x[d][p] = coords[d]
        end
        x = coords[1]
        y = D >= 2 ? coords[2] : 0.0
        z = D >= 3 ? coords[3] : 0.0
        ps.v[1][p] = 0.03 * sin(x) + 0.01 * cos(y)
        ps.v[2][p] = 0.02 * cos(x) + 0.005 * sin(z)
        ps.v[3][p] = -0.015 * sin(x + y)
        ps.weight[p] = ΔV
        ps.id[p] = UInt64(30_000 + p)
        ps.tag[p] = UInt32(40_000 + p)
    end
    return g, ps
end

function _set_field_coupled_B!(st::HybridStepper{D,T}) where {D,T}
    g = st.g
    for c = 1:3
        fill!(st.fields.B[c], zero(T))
    end
    for I in CartesianIndices(g.n)
        x = (I[1] - 1) * g.dx[1]
        st.fields.B[2][I] = T(0.01) * sin(T(2) * x)
        st.fields.B[3][I] = one(T) + T(0.02) * cos(x)
    end
    return st
end

function _max_field_error(a::HybridFields{D,T}, b::HybridFields{D,T}) where {D,T}
    err = maximum(abs, a.n .- b.n)
    err = max(err, maximum(abs, a.pe .- b.pe))
    err = max(err, maximum(abs, a.ninv .- b.ninv))
    for d = 1:D
        err = max(err, maximum(abs, a.gradp[d] .- b.gradp[d]))
    end
    for c = 1:3
        err = max(err, maximum(abs, a.ui[c] .- b.ui[c]))
        err = max(err, maximum(abs, a.B[c] .- b.B[c]))
        err = max(err, maximum(abs, a.E[c] .- b.E[c]))
        err = max(err, maximum(abs, a.J[c] .- b.J[c]))
    end
    return err
end

_face_side_code(offset) = offset < 0 ? 1 : 2
_face_particle_id(source, axis, offset) = UInt64(100 * source + 10 * axis + _face_side_code(offset))

function _set_particle_fields!(ps::ParticleSet{D,Float64}, p, id, source, axis, offset) where {D}
    side = _face_side_code(offset)
    ps.id[p] = id
    ps.tag[p] = UInt32(1000 + id)
    ps.weight[p] = 1.0 + source / 10.0 + axis / 100.0 + side / 1000.0
    ps.v[1][p] = source
    ps.v[2][p] = axis
    ps.v[3][p] = side
    return ps
end

function _expected_face_ids(layout::LogicalRankLayout{D}, dest) where {D}
    ids = UInt64[]
    for source = 1:nranks(layout)
        coords = rank_coords(layout, source)
        for axis = 1:D, offset in (-1, 1)
            target_coords = ntuple(d -> d == axis ? coords[d] + offset : coords[d], D)
            rank_index(layout, target_coords) == dest || continue
            push!(ids, _face_particle_id(source, axis, offset))
        end
    end
    return sort!(ids)
end

function _decode_face_particle_id(id::UInt64)
    n = Int(id)
    source = n ÷ 100
    axis = (n % 100) ÷ 10
    side = n % 10
    return source, axis, side
end

function _assert_face_particle_fields(ps::ParticleSet{D,Float64}) where {D}
    for p = 1:nparticles(ps)
        source, axis, side = _decode_face_particle_id(ps.id[p])
        @test ps.tag[p] == UInt32(1000 + ps.id[p])
        @test ps.weight[p] ≈ 1.0 + source / 10.0 + axis / 100.0 + side / 1000.0
        @test ps.v[1][p] ≈ source
        @test ps.v[2][p] ≈ axis
        @test ps.v[3][p] ≈ side
    end
end

_halo_field(rank) =
    [100.0 * rank + 1.0, 10.0 * rank, 10.0 * rank + 1.0, 10.0 * rank + 2.0, 100.0 * rank + 2.0]
_halo_moment(rank) =
    [200.0 * rank + 1.0, 20.0 * rank, 20.0 * rank + 1.0, 20.0 * rank + 2.0, 200.0 * rank + 2.0]

function _halo_matrix(rank; scale)
    A = Matrix{Float64}(undef, 5, 3)
    for j = 1:3, i = 1:5
        A[i, j] = scale * rank + 10.0 * i + j
    end
    return A
end

@testset "multi-rank MPI launcher guard" begin
    if WORLD_SIZE == 1
        @test WORLD_RANK == 0
        @info "test_mpi_multirank.jl validates real multi-rank MPI only under mpiexec -n 2, -n 4, or -n 8"
    else
        @test WORLD_SIZE in (2, 4, 8)
    end
end

if WORLD_SIZE in (2, 4, 8)
    layout = _multirank_layout(WORLD_SIZE)
    ctx = create_cartesian_communicator(layout; comm = WORLD, reorder = false)
    try
        @testset "real MPI Cartesian mapping size=$WORLD_SIZE" begin
            @test ctx.mpi_size == WORLD_SIZE
            @test ctx.mpi_rank == MPI.Comm_rank(ctx.comm)
            @test MPI.Comm_size(ctx.comm) == WORLD_SIZE
            @test rank_index(ctx.layout, ctx.coords) == ctx.logical_rank
            @test rank_coords(ctx.layout, ctx.logical_rank) == ctx.coords
            @test MPI.Cart_coords(ctx.comm, ctx.mpi_rank) ==
                  [ctx.coords[d] - 1 for d = 1:length(ctx.coords)]
            @test MPI.Cart_rank(ctx.comm, [ctx.coords[d] - 1 for d = 1:length(ctx.coords)]) == ctx.mpi_rank

            gathered_logical = MPI.Allgather(ctx.logical_rank, ctx.comm)
            gathered_mpi = MPI.Allgather(ctx.mpi_rank, ctx.comm)
            @test sort(gathered_logical) == collect(1:WORLD_SIZE)
            @test sort(gathered_mpi) == collect(0:(WORLD_SIZE-1))

            for axis = 1:length(ctx.coords)
                source, dest = MPI.Cart_shift(ctx.comm, axis - 1, 1)
                @test source == _cart_rank_for_logical(ctx, mpi_cartesian_neighbor(ctx, axis, -1))
                @test dest == _cart_rank_for_logical(ctx, mpi_cartesian_neighbor(ctx, axis, 1))
            end
        end

        @testset "rank bounds classify local midpoint size=$WORLD_SIZE" begin
            D = length(ctx.coords)
            g = FourierGrid(ntuple(d -> 4 * WORLD_SIZE + d, D), ntuple(d -> 10.0 + d, D))
            bounds = rank_bounds(g, ctx.layout, ctx.logical_rank)
            midpoint = ntuple(d -> 0.5 * (bounds.lo[d] + bounds.hi[d]), D)
            @test rank_of_position(midpoint, g, ctx.layout) == ctx.logical_rank

            gathered_classification =
                MPI.Allgather(rank_of_position(midpoint, g, ctx.layout), ctx.comm)
            @test sort(gathered_classification) == collect(1:WORLD_SIZE)
        end

        @testset "real MPI diagnostic allreduce size=$WORLD_SIZE" begin
            local_value = (
                rank = ctx.logical_rank,
                vector = Float64[ctx.logical_rank, ctx.logical_rank^2, (-1)^ctx.logical_rank],
                nested = (lo = ctx.logical_rank, hi = -ctx.logical_rank),
            )
            expected_locals = [
                (rank = r, vector = Float64[r, r^2, (-1)^r], nested = (lo = r, hi = -r)) for
                r = 1:WORLD_SIZE
            ]
            status = GPUAwareMPIStatus(false, false, false, :test, "host-only multi-rank test")

            @test mpi_allreduce_diagnostics(local_value, ctx; op = :sum, gpu_status = status) ==
                  reduce_diagnostics(expected_locals; op = :sum)
            @test mpi_allreduce_diagnostics(local_value, ctx; op = :min, gpu_status = status) ==
                  reduce_diagnostics(expected_locals; op = :min)
            @test mpi_allreduce_diagnostics(local_value, ctx; op = :max, gpu_status = status) ==
                  reduce_diagnostics(expected_locals; op = :max)
        end

        @testset "destination-routed MPI byte transport size=$WORLD_SIZE" begin
            dest = mod(ctx.mpi_rank + 1, WORLD_SIZE)
            source = mod(ctx.mpi_rank - 1, WORLD_SIZE)
            send_chunks = [UInt8[] for _ = 1:WORLD_SIZE]
            send_chunks[dest+1] = UInt8[0x40+UInt8(ctx.mpi_rank)]
            recv_chunks = HybridPlasmaPIC._mpi_alltoallv_bytes(send_chunks, ctx.comm)
            for rank0 = 0:(WORLD_SIZE-1)
                if rank0 == source
                    @test recv_chunks[rank0+1] == UInt8[0x40+UInt8(rank0)]
                else
                    @test isempty(recv_chunks[rank0+1])
                end
            end
        end

        @testset "real MPI slab field halo exchange size=$WORLD_SIZE" begin
            slab_layout = LogicalRankLayout((WORLD_SIZE,); periodic = (false,))
            slab_ctx = create_cartesian_communicator(slab_layout; comm = WORLD, reorder = false)
            try
                reference = [(_halo_field(r), _halo_field(r) .+ 1000.0) for r = 1:WORLD_SIZE]
                local_fields = (
                    copy(reference[slab_ctx.logical_rank][1]),
                    copy(reference[slab_ctx.logical_rank][2]),
                )
                expected = [(copy(a), copy(b)) for (a, b) in reference]
                expected_stats =
                    exchange_field_halos!(expected, slab_layout; halo = 1, fill_value = -99.0)

                stats =
                    mpi_exchange_field_halos!(local_fields, slab_ctx; halo = 1, fill_value = -99.0)
                @test stats == expected_stats
                @test local_fields[1] == expected[slab_ctx.logical_rank][1]
                @test local_fields[2] == expected[slab_ctx.logical_rank][2]
            finally
                free_mpi_communicator!(slab_ctx)
            end
        end

        @testset "real MPI slab ghost moment exchange size=$WORLD_SIZE" begin
            slab_layout = LogicalRankLayout((WORLD_SIZE,); periodic = (false,))
            slab_ctx = create_cartesian_communicator(slab_layout; comm = WORLD, reorder = false)
            try
                reference = [(_halo_moment(r), _halo_moment(r) .+ 1000.0) for r = 1:WORLD_SIZE]
                local_moments = (
                    copy(reference[slab_ctx.logical_rank][1]),
                    copy(reference[slab_ctx.logical_rank][2]),
                )
                expected = [(copy(a), copy(b)) for (a, b) in reference]
                expected_stats = exchange_ghost_moments!(expected, slab_layout; halo = 1)

                stats = mpi_exchange_ghost_moments!(local_moments, slab_ctx; halo = 1)
                @test stats == expected_stats
                @test local_moments[1] == expected[slab_ctx.logical_rank][1]
                @test local_moments[2] == expected[slab_ctx.logical_rank][2]
            finally
                free_mpi_communicator!(slab_ctx)
            end
        end

        @testset "real MPI periodic slab halo exchange size=$WORLD_SIZE" begin
            slab_layout = LogicalRankLayout((WORLD_SIZE,); periodic = (true,))
            slab_ctx = create_cartesian_communicator(slab_layout; comm = WORLD, reorder = false)
            try
                field_ref = [_halo_field(r) for r = 1:WORLD_SIZE]
                local_field = copy(field_ref[slab_ctx.logical_rank])
                expected_field = [copy(A) for A in field_ref]
                expected_field_stats = exchange_field_halos!(expected_field, slab_layout; halo = 1)

                field_stats = mpi_exchange_field_halos!(local_field, slab_ctx; halo = 1)
                @test field_stats == expected_field_stats
                @test local_field == expected_field[slab_ctx.logical_rank]

                moment_ref = [_halo_moment(r) for r = 1:WORLD_SIZE]
                local_moment = copy(moment_ref[slab_ctx.logical_rank])
                expected_moment = [copy(A) for A in moment_ref]
                expected_moment_stats =
                    exchange_ghost_moments!(expected_moment, slab_layout; halo = 1)

                moment_stats = mpi_exchange_ghost_moments!(local_moment, slab_ctx; halo = 1)
                @test moment_stats == expected_moment_stats
                @test local_moment == expected_moment[slab_ctx.logical_rank]
            finally
                free_mpi_communicator!(slab_ctx)
            end
        end

        @testset "real MPI 2D slab halo exchange size=$WORLD_SIZE" begin
            slab_layout = LogicalRankLayout((WORLD_SIZE, 1); periodic = (false, true))
            slab_ctx = create_cartesian_communicator(slab_layout; comm = WORLD, reorder = false)
            try
                field_ref = [_halo_matrix(r; scale = 1000.0) for r = 1:WORLD_SIZE]
                local_field = copy(field_ref[slab_ctx.logical_rank])
                expected_field = [copy(A) for A in field_ref]
                expected_field_stats =
                    exchange_field_halos!(expected_field, slab_layout; halo = 1, fill_value = -77.0)

                field_stats =
                    mpi_exchange_field_halos!(local_field, slab_ctx; halo = 1, fill_value = -77.0)
                @test field_stats == expected_field_stats
                @test local_field == expected_field[slab_ctx.logical_rank]

                moment_ref = [_halo_matrix(r; scale = 2000.0) for r = 1:WORLD_SIZE]
                local_moment = copy(moment_ref[slab_ctx.logical_rank])
                expected_moment = [copy(A) for A in moment_ref]
                expected_moment_stats =
                    exchange_ghost_moments!(expected_moment, slab_layout; halo = 1)

                moment_stats = mpi_exchange_ghost_moments!(local_moment, slab_ctx; halo = 1)
                @test moment_stats == expected_moment_stats
                @test local_moment == expected_moment[slab_ctx.logical_rank]
            finally
                free_mpi_communicator!(slab_ctx)
            end
        end

        @testset "real MPI slab halo count mismatch rejects collectively size=$WORLD_SIZE" begin
            slab_layout = LogicalRankLayout((WORLD_SIZE, 1); periodic = (false, true))
            slab_ctx = create_cartesian_communicator(slab_layout; comm = WORLD, reorder = false)
            try
                ntrans = slab_ctx.logical_rank == 1 ? 2 : 3
                bad_field = fill(Float64(slab_ctx.logical_rank), 5, ntrans)
                @test_throws DimensionMismatch mpi_exchange_field_halos!(
                    bad_field,
                    slab_ctx;
                    halo = 1,
                )

                bad_moment = fill(10.0 * slab_ctx.logical_rank, 5, ntrans)
                @test_throws DimensionMismatch mpi_exchange_ghost_moments!(
                    bad_moment,
                    slab_ctx;
                    halo = 1,
                )
            finally
                free_mpi_communicator!(slab_ctx)
            end
        end

        @testset "rank-count invariant synthetic particle diagnostics size=$WORLD_SIZE" begin
            local_ps = _particle_partition(WORLD_RANK, WORLD_SIZE)
            serial_ps = _particle_partition(0, 1)
            status = GPUAwareMPIStatus(false, false, false, :test, "host-only multi-rank test")

            reduced = mpi_allreduce_diagnostics(
                _particle_diagnostics(local_ps),
                ctx;
                op = :sum,
                gpu_status = status,
            )
            expected = _particle_diagnostics(serial_ps)
            @test reduced.number ≈ expected.number
            @test reduced.momentum[1] ≈ expected.momentum[1]
            @test reduced.momentum[2] ≈ expected.momentum[2]
            @test reduced.momentum[3] ≈ expected.momentum[3]
            @test reduced.kinetic ≈ expected.kinetic
        end

        @testset "time-advanced MPI particle budgets invariant size=$WORLD_SIZE" begin
            D = length(ctx.coords)
            budget_layout = LogicalRankLayout(ctx.layout.ranks; periodic = ntuple(_ -> true, D))
            budget_ctx = create_cartesian_communicator(budget_layout; comm = WORLD, reorder = false)
            try
                g, serial0 = _budget_reference_state(Val(D))
                ps = _rank_budget_particles(serial0, g, budget_ctx.layout, budget_ctx.logical_rank)
                _, serial_ref = _budget_reference_state(Val(D))
                dt = 0.37
                nsteps = 9
                moved_total = 0
                for _ = 1:nsteps
                    _advance_ballistic_periodic!(ps, g, dt)
                    stats = mpi_migrate_particles!(ps, g, budget_ctx)
                    @test stats.lost == 0
                    moved_total += stats.moved
                    for p = 1:nparticles(ps)
                        pos = ntuple(d -> ps.x[d][p], D)
                        @test rank_of_position(pos, g, budget_ctx.layout) == budget_ctx.logical_rank
                    end
                end

                for _ = 1:nsteps
                    _advance_ballistic_periodic!(serial_ref, g, dt)
                end

                status = GPUAwareMPIStatus(false, false, false, :test, "host-only multi-rank test")
                reduced = mpi_allreduce_diagnostics(
                    _particle_diagnostics(ps),
                    budget_ctx;
                    op = :sum,
                    gpu_status = status,
                )
                expected = _particle_diagnostics(serial_ref)
                global_count = MPI.Allreduce(nparticles(ps), +, budget_ctx.comm)

                @test moved_total > 0
                @test global_count == nparticles(serial_ref)
                @test reduced.number ≈ expected.number rtol = 32eps(Float64)
                @test reduced.momentum[1] ≈ expected.momentum[1] rtol = 64eps(Float64) atol =
                    64eps(Float64)
                @test reduced.momentum[2] ≈ expected.momentum[2] rtol = 64eps(Float64) atol =
                    64eps(Float64)
                @test reduced.momentum[3] ≈ expected.momentum[3] rtol = 64eps(Float64) atol =
                    64eps(Float64)
                @test reduced.kinetic ≈ expected.kinetic rtol = 64eps(Float64)
            finally
                free_mpi_communicator!(budget_ctx)
            end
        end

        @testset "field-coupled MPI agreement with serial HybridStepper size=$WORLD_SIZE" begin
            D = length(ctx.coords)
            g, serial0 = _field_coupled_reference_state(Val(D))
            serial = _copy_particle_subset(serial0, collect(1:nparticles(serial0)))
            local_ps = _rank_budget_particles(serial0, g, ctx.layout, ctx.logical_rank)
            model = HybridModel(IsothermalElectrons(0.03); nfloor = 1e-5)
            shape = CIC()
            serial_st = HybridStepper(g, model, shape, nparticles(serial))
            mpi_st = HybridStepper(g, model, shape, nparticles(local_ps))
            _set_field_coupled_B!(serial_st)
            _set_field_coupled_B!(mpi_st)

            status = GPUAwareMPIStatus(false, false, false, :test, "host-only multi-rank test")
            init!(serial_st, serial)
            mpi_init!(mpi_st, local_ps, ctx; gpu_status = status)
            @test _max_field_error(mpi_st.fields, serial_st.fields) < 2.0e-12

            dt = 0.025
            for _ = 1:3
                step!(serial_st, serial, dt; NB = 2)
                mpi_step!(mpi_st, local_ps, ctx, dt; NB = 2, gpu_status = status)
                @test _max_field_error(mpi_st.fields, serial_st.fields) < 2.0e-10

                reduced = mpi_allreduce_diagnostics(
                    _particle_diagnostics(local_ps),
                    ctx;
                    op = :sum,
                    gpu_status = status,
                )
                expected = _particle_diagnostics(serial)
                @test reduced.number ≈ expected.number rtol = 32eps(Float64)
                @test reduced.momentum[1] ≈ expected.momentum[1] rtol = 64eps(Float64) atol =
                    64eps(Float64)
                @test reduced.momentum[2] ≈ expected.momentum[2] rtol = 64eps(Float64) atol =
                    64eps(Float64)
                @test reduced.momentum[3] ≈ expected.momentum[3] rtol = 64eps(Float64) atol =
                    64eps(Float64)
                @test reduced.kinetic ≈ expected.kinetic rtol = 64eps(Float64)
            end
        end

        @testset "distributed MPI checkpoint restart bitmatch size=$WORLD_SIZE" begin
            D = length(ctx.coords)
            g, serial0 = _field_coupled_reference_state(Val(D))
            local_ps = _rank_budget_particles(serial0, g, ctx.layout, ctx.logical_rank)
            model = HybridModel(IsothermalElectrons(0.03); nfloor = 1e-5)
            shape = CIC()
            mpi_st = HybridStepper(g, model, shape, nparticles(local_ps))
            _set_field_coupled_B!(mpi_st)
            status = GPUAwareMPIStatus(false, false, false, :test, "host-only multi-rank test")
            mpi_init!(mpi_st, local_ps, ctx; gpu_status = status)

            dt = 0.025
            for _ = 1:2
                mpi_step!(mpi_st, local_ps, ctx, dt; NB = 2, gpu_status = status)
            end

            dir = joinpath(tempdir(), "HybridPlasmaPIC_mpi_checkpoint_size$(WORLD_SIZE)")
            ctx.mpi_rank == 0 && rm(dir; recursive = true, force = true)
            MPI.Barrier(ctx.comm)
            manifest_path = save_mpi_checkpoint(dir, mpi_st, local_ps, ctx)
            @test isfile(manifest_path)

            restart_ps = ParticleSet{D,Float64}(1; q = local_ps.q, m = local_ps.m)
            restart_st = HybridStepper(g, model, shape, 1)
            load_mpi_checkpoint!(restart_st, restart_ps, dir, ctx)
            @test restart_st.step[] == mpi_st.step[]
            @test restart_st.time[] == mpi_st.time[]
            @test nparticles(restart_ps) == nparticles(local_ps)
            @test length(restart_st.work) == nparticles(restart_ps)
            @test all(length(restart_st.Ep[c]) == nparticles(restart_ps) for c = 1:3)
            @test all(length(restart_st.Bp[c]) == nparticles(restart_ps) for c = 1:3)
            @test all(length(restart_st.xmid[d]) == nparticles(restart_ps) for d = 1:D)
            for d = 1:D
                @test restart_ps.x[d] == local_ps.x[d]
            end
            for c = 1:3
                @test restart_ps.v[c] == local_ps.v[c]
                @test restart_st.fields.B[c] == mpi_st.fields.B[c]
                @test restart_st.fields.E[c] == mpi_st.fields.E[c]
            end
            @test restart_ps.weight == local_ps.weight
            @test restart_ps.id == local_ps.id
            @test restart_ps.tag == local_ps.tag

            for _ = 1:3
                mpi_step!(mpi_st, local_ps, ctx, dt; NB = 2, gpu_status = status)
                mpi_step!(restart_st, restart_ps, ctx, dt; NB = 2, gpu_status = status)
            end
            for d = 1:D
                @test restart_ps.x[d] == local_ps.x[d]
            end
            for c = 1:3
                @test restart_ps.v[c] == local_ps.v[c]
                @test restart_st.fields.B[c] == mpi_st.fields.B[c]
                @test restart_st.fields.E[c] == mpi_st.fields.E[c]
            end
            @test restart_ps.weight == local_ps.weight
            @test restart_ps.id == local_ps.id
            @test restart_ps.tag == local_ps.tag
            @test restart_st.step[] == mpi_st.step[]
            @test restart_st.time[] == mpi_st.time[]

            MPI.Barrier(ctx.comm)
            ctx.mpi_rank == 0 && rm(dir; recursive = true, force = true)
            MPI.Barrier(ctx.comm)
        end

        if WORLD_SIZE == 2
            @testset "real MPI particle migration nonperiodic slab size=2" begin
                g = FourierGrid((8,), (8.0,))
                ps = ParticleSet{1,Float64}(2)
                if ctx.logical_rank == 1
                    ps.x[1] .= [-0.05, 4.25]
                    ids = UInt64[
                        _face_particle_id(ctx.logical_rank, 1, -1),
                        _face_particle_id(ctx.logical_rank, 1, 1),
                    ]
                else
                    ps.x[1] .= [3.75, 8.05]
                    ids = UInt64[
                        _face_particle_id(ctx.logical_rank, 1, -1),
                        _face_particle_id(ctx.logical_rank, 1, 1),
                    ]
                end
                for p = 1:2
                    _set_particle_fields!(ps, p, ids[p], ctx.logical_rank, 1, p == 1 ? -1 : 1)
                end

                stats = mpi_migrate_particles!(ps, g, ctx)
                @test stats.moved == 2
                @test stats.lost == 2
                @test stats.sent == 1
                @test stats.received == 1
                @test nparticles(ps) == 1
                expected_id =
                    ctx.logical_rank == 1 ? UInt64[_face_particle_id(2, 1, -1)] :
                    UInt64[_face_particle_id(1, 1, 1)]
                @test sort(ps.id) == expected_id
                @test rank_of_position((ps.x[1][1],), g, ctx.layout) == ctx.logical_rank
                _assert_face_particle_fields(ps)
            end
        end

        if WORLD_SIZE == 8
            @testset "real MPI particle migration across every 3D face size=8" begin
                g = FourierGrid((16, 18, 20), (8.0, 9.0, 10.0))
                bounds = rank_bounds(g, ctx.layout, ctx.logical_rank)
                ps = ParticleSet{3,Float64}(6)
                p = 0
                for axis = 1:3, offset in (-1, 1)
                    p += 1
                    x = ntuple(
                        d -> begin
                            if d == axis
                                width = (bounds.hi[d] - bounds.lo[d])
                                offset < 0 ? bounds.lo[d] - 0.125 * width :
                                bounds.hi[d] + 0.125 * width
                            else
                                0.5 * (bounds.lo[d] + bounds.hi[d])
                            end
                        end,
                        3,
                    )
                    for d = 1:3
                        ps.x[d][p] = x[d]
                    end
                    id = _face_particle_id(ctx.logical_rank, axis, offset)
                    _set_particle_fields!(ps, p, id, ctx.logical_rank, axis, offset)
                end

                stats = mpi_migrate_particles!(ps, g, ctx)
                @test stats.moved == 8 * 6
                @test stats.lost == 0
                @test stats.sent == 6
                @test stats.received == 6
                @test sort(ps.id) == _expected_face_ids(ctx.layout, ctx.logical_rank)
                for p = 1:nparticles(ps)
                    pos = ntuple(d -> ps.x[d][p], 3)
                    @test rank_of_position(pos, g, ctx.layout) == ctx.logical_rank
                    @test all(d -> 0.0 <= pos[d] < g.L[d], 1:3)
                end
                _assert_face_particle_fields(ps)
            end
        end
    finally
        free_mpi_communicator!(ctx)
    end
end
