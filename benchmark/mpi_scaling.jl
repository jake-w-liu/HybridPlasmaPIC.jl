# mpi_scaling.jl — reproducible real-MPI scaling harness.
#
# Run through MPI.jl's configured launcher, for example:
#
#   julia --project=. -e 'import MPI; run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=. benchmark/mpi_scaling.jl`)'
#
# This is a measurement harness, not a pass/fail correctness test. It prints one
# TSV row per metric from rank 0. Use cluster runs to collect production scaling
# data; local laptop runs are smoke tests for the measurement path.

using HybridPlasmaPIC
using Printf
import MPI

const _DEFAULTS = (
    particles_per_rank = 4096,
    cells = 64,
    transverse = 64,
    halo = 1,
    warmup = 1,
    reps = 5,
    dt = 0.02,
    nb = 2,
)

function _usage()
    return """
    Usage: benchmark/mpi_scaling.jl [options]

    Options:
      --particles-per-rank N   rank-local particles for mpi_step! ($( _DEFAULTS.particles_per_rank ))
      --cells N                cells along the decomposed/local axis ($( _DEFAULTS.cells ))
      --transverse N           transverse cells for 2D halo benchmarks ($( _DEFAULTS.transverse ))
      --halo N                 ghost-cell width for halo benchmarks ($( _DEFAULTS.halo ))
      --warmup N               warmup repetitions before timing ($( _DEFAULTS.warmup ))
      --reps N                 measured repetitions ($( _DEFAULTS.reps ))
      --dt X                   timestep for mpi_step! ($( _DEFAULTS.dt ))
      --nb N                   magnetic subcycles for mpi_step! ($( _DEFAULTS.nb ))
      --help                   print this message
    """
end

function _parse_positive_int(value::AbstractString, key::Symbol)
    n = tryparse(Int, value)
    n !== nothing && n >= 1 || throw(ArgumentError("$key must be a positive integer, got $value"))
    return n
end

function _parse_nonnegative_int(value::AbstractString, key::Symbol)
    n = tryparse(Int, value)
    n !== nothing && n >= 0 ||
        throw(ArgumentError("$key must be a nonnegative integer, got $value"))
    return n
end

function _parse_positive_float(value::AbstractString, key::Symbol)
    x = tryparse(Float64, value)
    x !== nothing && isfinite(x) && x > 0 ||
        throw(ArgumentError("$key must be a positive finite number, got $value"))
    return x
end

function _parse_options(args)
    opts = Dict{Symbol,Any}(pairs(_DEFAULTS))
    i = 1
    while i <= length(args)
        arg = args[i]
        arg == "--help" && return nothing
        startswith(arg, "--") || throw(ArgumentError("unexpected argument $arg"))
        i < length(args) || throw(ArgumentError("$arg requires a value"))
        key = Symbol(replace(arg[3:end], "-" => "_"))
        haskey(opts, key) || throw(ArgumentError("unknown option $arg"))
        raw = args[i+1]
        if key === :dt
            opts[key] = _parse_positive_float(raw, key)
        elseif key === :warmup
            opts[key] = _parse_nonnegative_int(raw, key)
        else
            opts[key] = _parse_positive_int(raw, key)
        end
        i += 2
    end
    opts[:cells] >= 3opts[:halo] ||
        throw(ArgumentError("--cells must be at least 3*--halo for slab halo exchange"))
    return NamedTuple{keys(_DEFAULTS)}(Tuple(opts[k] for k in keys(_DEFAULTS)))
end

function _timed_max(f, comm::MPI.Comm)
    MPI.Barrier(comm)
    local_elapsed = @elapsed begin
        f()
        MPI.Barrier(comm)
    end
    return MPI.Allreduce(local_elapsed, max, comm)
end

function _fill_rank_particles!(
    ps::ParticleSet{1,Float64},
    g::FourierGrid{1,Float64},
    layout::LogicalRankLayout{1},
    logical_rank::Integer,
)
    bounds = rank_bounds(g, layout, logical_rank)
    n = nparticles(ps)
    width = bounds.hi[1] - bounds.lo[1]
    @inbounds for p = 1:n
        frac = (p - 0.5) / n
        ps.x[1][p] = bounds.lo[1] + width * frac
        phase = 0.01 * (p + 97 * logical_rank)
        ps.v[1][p] = 0.05 * sin(phase)
        ps.v[2][p] = 0.03 * cos(2phase)
        ps.v[3][p] = -0.02 * sin(3phase)
        ps.weight[p] = g.L[1] / (n * nranks(layout))
        gid = UInt64(logical_rank - 1) * UInt64(n) + UInt64(p)
        ps.id[p] = gid
        ps.tag[p] = UInt32(gid % UInt64(typemax(UInt32)))
    end
    return ps
end

function _bench_mpi_step(opts, comm::MPI.Comm)
    nr = MPI.Comm_size(comm)
    layout = LogicalRankLayout((nr,); periodic = (true,))
    ctx = create_cartesian_communicator(layout; comm, reorder = false)
    try
        g = FourierGrid((opts.cells,), (2π,))
        ps = ParticleSet{1,Float64}(opts.particles_per_rank)
        _fill_rank_particles!(ps, g, layout, ctx.logical_rank)
        model = HybridModel(IsothermalElectrons(0.03); nfloor = 1e-5)
        st = HybridStepper(g, model, CIC(), nparticles(ps))
        fill!(st.fields.B[3], 1.0)
        status =
            GPUAwareMPIStatus(false, false, false, :benchmark, "host-only MPI scaling benchmark")
        mpi_init!(st, ps, ctx; gpu_status = status)

        for _ = 1:opts.warmup
            mpi_step!(st, ps, ctx, opts.dt; NB = opts.nb, gpu_status = status)
        end

        seconds = _timed_max(comm) do
            for _ = 1:opts.reps
                mpi_step!(st, ps, ctx, opts.dt; NB = opts.nb, gpu_status = status)
            end
        end
        work = opts.reps * opts.particles_per_rank * nr
        return (; seconds, work, throughput = work / seconds)
    finally
        free_mpi_communicator!(ctx)
    end
end

function _fill_halo_array!(A::AbstractArray{Float64}, logical_rank::Integer, scale::Float64)
    @inbounds for I in CartesianIndices(A)
        A[I] = scale * logical_rank + sum(Tuple(I))
    end
    return A
end

function _bench_field_halo(opts, comm::MPI.Comm)
    nr = MPI.Comm_size(comm)
    layout = LogicalRankLayout((nr, 1); periodic = (true, true))
    ctx = create_cartesian_communicator(layout; comm, reorder = false)
    try
        A = Array{Float64}(undef, opts.cells, opts.transverse)
        _fill_halo_array!(A, ctx.logical_rank, 1000.0)
        for _ = 1:opts.warmup
            mpi_exchange_field_halos!(A, ctx; halo = opts.halo)
        end
        stats_ref = mpi_exchange_field_halos!(A, ctx; halo = opts.halo)
        seconds = _timed_max(comm) do
            for _ = 1:opts.reps
                mpi_exchange_field_halos!(A, ctx; halo = opts.halo)
            end
        end
        work = opts.reps * stats_ref.exchanged
        return (; seconds, work, throughput = work / seconds)
    finally
        free_mpi_communicator!(ctx)
    end
end

function _bench_ghost_moments(opts, comm::MPI.Comm)
    nr = MPI.Comm_size(comm)
    layout = LogicalRankLayout((nr, 1); periodic = (true, true))
    ctx = create_cartesian_communicator(layout; comm, reorder = false)
    try
        A = Array{Float64}(undef, opts.cells, opts.transverse)
        _fill_halo_array!(A, ctx.logical_rank, 2000.0)
        for _ = 1:opts.warmup
            _fill_halo_array!(A, ctx.logical_rank, 2000.0)
            mpi_exchange_ghost_moments!(A, ctx; halo = opts.halo)
        end
        _fill_halo_array!(A, ctx.logical_rank, 2000.0)
        stats_ref = mpi_exchange_ghost_moments!(A, ctx; halo = opts.halo)
        seconds = _timed_max(comm) do
            for _ = 1:opts.reps
                _fill_halo_array!(A, ctx.logical_rank, 2000.0)
                mpi_exchange_ghost_moments!(A, ctx; halo = opts.halo)
            end
        end
        work = opts.reps * stats_ref.exchanged
        return (; seconds, work, throughput = work / seconds)
    finally
        free_mpi_communicator!(ctx)
    end
end

function _print_metric(name, nranks, result)
    @printf(
        "%s\t%d\t%.9f\t%d\t%.6e\n",
        name,
        nranks,
        result.seconds,
        result.work,
        result.throughput
    )
end

function main(args = ARGS)
    ensure_mpi_initialized!()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    opts = _parse_options(args)
    if opts === nothing
        rank == 0 && print(_usage())
        return nothing
    end

    if rank == 0
        println("# HybridPlasmaPIC MPI scaling harness")
        println("# options=$(opts)")
        println("metric\tnranks\tseconds\twork\tthroughput_per_s")
    end

    step_result = _bench_mpi_step(opts, comm)
    field_result = _bench_field_halo(opts, comm)
    moment_result = _bench_ghost_moments(opts, comm)

    if rank == 0
        _print_metric("mpi_step_particle_steps", nranks, step_result)
        _print_metric("field_halo_scalar_values", nranks, field_result)
        _print_metric("ghost_moment_scalar_values", nranks, moment_result)
    end
    return nothing
end

main()
