# mpi.jl (§5.3 parallel/)
#
# Thin MPI.jl transport hooks over the deterministic logical-rank reference
# kernels in domain_decomposition.jl. These wrappers are intentionally small:
# they create real MPI Cartesian communicators, expose rank/layout provenance,
# and provide structured collective diagnostics, real slab halo exchanges, and
# destination-routed particle migration. Scaling still requires cluster tests.

import MPI
import Serialization

"""
    GPUAwareMPIStatus

Result of querying whether MPI can accept GPU-resident buffers directly.
`enabled` is true when at least one supported GPU transport is reported.
`cuda` and `rocm` mirror MPI.jl's capability checks. `source` and `reason`
record how the decision was made for run metadata and debugging.
"""
struct GPUAwareMPIStatus
    enabled::Bool
    cuda::Bool
    rocm::Bool
    source::Symbol
    reason::String
end

"""
    ensure_mpi_initialized!()

Initialize MPI if it has not already been initialized. Throws if MPI has already
been finalized in this process.
"""
function ensure_mpi_initialized!()
    MPI.Finalized() && error("MPI has already been finalized")
    MPI.Initialized() || MPI.Init()
    return nothing
end

mpi_initialized() = MPI.Initialized() && !MPI.Finalized()

function _mpi_bool_capability(f)
    try
        return Bool(f()), nothing
    catch err
        return false, sprint(showerror, err)
    end
end

"""
    gpu_aware_mpi_status(; initialize=true) -> GPUAwareMPIStatus

Query MPI.jl's CUDA/ROCm-aware transport capability. The implementation uses
`MPI.has_cuda()` and `MPI.has_rocm()`; those functions can also be controlled
through MPI.jl's documented `JULIA_MPI_HAS_CUDA` and `JULIA_MPI_HAS_ROCM`
environment overrides. MPI is initialized by default because some MPI
implementations require initialization before capability queries are valid.
"""
function gpu_aware_mpi_status(; initialize::Bool = true)
    initialize && ensure_mpi_initialized!()
    cuda, cuda_error = _mpi_bool_capability(MPI.has_cuda)
    rocm, rocm_error = _mpi_bool_capability(MPI.has_rocm)
    enabled = cuda || rocm
    details = String[]
    push!(details, "cuda=$(cuda)")
    push!(details, "rocm=$(rocm)")
    cuda_error === nothing || push!(details, "cuda_error=$(cuda_error)")
    rocm_error === nothing || push!(details, "rocm_error=$(rocm_error)")
    reason = enabled ? "MPI reports GPU-aware transport: " : "MPI reports no GPU-aware transport: "
    return GPUAwareMPIStatus(enabled, cuda, rocm, :mpi, reason * join(details, "; "))
end

function mpi_comm_size(comm::MPI.Comm = MPI.COMM_WORLD)
    ensure_mpi_initialized!()
    return MPI.Comm_size(comm)
end

function mpi_comm_rank(comm::MPI.Comm = MPI.COMM_WORLD)
    ensure_mpi_initialized!()
    return MPI.Comm_rank(comm)
end

"""
    mpi_dims_create(nranks, D) -> NTuple{D,Int}

Use `MPI_Dims_create` to choose a Cartesian process grid for `nranks` ranks and
`D` dimensions.
"""
function mpi_dims_create(nranks::Integer, D::Integer)
    nranks > 0 || throw(ArgumentError("nranks must be positive, got $nranks"))
    D > 0 || throw(ArgumentError("D must be positive, got $D"))
    ensure_mpi_initialized!()
    dims = MPI.Dims_create(Int(nranks), fill(0, Int(D)))
    return Tuple(Int.(dims))
end

"""
    MPICartesianCommunicator

Real MPI Cartesian communicator plus the matching [`LogicalRankLayout`](@ref).
`mpi_rank` is MPI's zero-based rank in `comm`; `coords` and `logical_rank` use
HybridPlasmaPIC's one-based logical-rank convention.
"""
struct MPICartesianCommunicator{D}
    comm::MPI.Comm
    layout::LogicalRankLayout{D}
    mpi_rank::Int
    mpi_size::Int
    coords::NTuple{D,Int}
    logical_rank::Int
end

"""
    create_cartesian_communicator(layout; comm=MPI.COMM_WORLD, reorder=false)

Create a real MPI Cartesian communicator whose dimensions and periodicity match
`layout`. The communicator size must equal `nranks(layout)`.
"""
function create_cartesian_communicator(
    layout::LogicalRankLayout{D};
    comm::MPI.Comm = MPI.COMM_WORLD,
    reorder::Bool = false,
) where {D}
    ensure_mpi_initialized!()
    ncomm = MPI.Comm_size(comm)
    ncomm == nranks(layout) || throw(
        ArgumentError("MPI communicator size $ncomm must equal nranks(layout)=$(nranks(layout))"),
    )

    dims = Int[layout.ranks...]
    periodic = Bool[layout.periodic...]
    cart = MPI.Cart_create(comm, dims; periodic, reorder)
    mpi_rank = MPI.Comm_rank(cart)
    coords0 = MPI.Cart_coords(cart, mpi_rank)
    coords = ntuple(d -> Int(coords0[d]) + 1, D)
    logical_rank = rank_index(layout, coords)
    logical_rank === nothing &&
        error("MPI Cartesian coordinates $coords are outside layout $layout")
    return MPICartesianCommunicator{D}(cart, layout, mpi_rank, ncomm, coords, logical_rank)
end

function free_mpi_communicator!(ctx::MPICartesianCommunicator)
    mpi_initialized() && MPI.free(ctx.comm)
    return nothing
end

"""
    mpi_cartesian_neighbor(ctx, axis, offset) -> Union{Int,Nothing}

Return the one-based logical-rank neighbor of `ctx` at signed `offset` along
one-based `axis`. Periodic wrapping and nonperiodic exterior handling follow
the stored [`LogicalRankLayout`](@ref).
"""
function mpi_cartesian_neighbor(
    ctx::MPICartesianCommunicator{D},
    axis::Integer,
    offset::Integer,
) where {D}
    ax = Int(axis)
    1 <= ax <= D || throw(ArgumentError("axis must be in 1:$D, got $axis"))
    off = Int(offset)
    coords = ntuple(d -> d == ax ? ctx.coords[d] + off : ctx.coords[d], D)
    return rank_index(ctx.layout, coords)
end

function _tuple_desc(xs)
    return string("(", join(xs, ","), ")")
end

"""
    mpi_rank_layout_description(ctx) -> String

Provenance string suitable for `RunMetadata.rank_layout`.
"""
function mpi_rank_layout_description(ctx::MPICartesianCommunicator)
    coords0 = ntuple(d -> ctx.coords[d] - 1, length(ctx.coords))
    return string(
        "mpi;ranks=",
        ctx.mpi_size,
        ";dims=",
        _tuple_desc(ctx.layout.ranks),
        ";coords=",
        _tuple_desc(coords0),
        ";periodic=",
        _tuple_desc(ctx.layout.periodic),
        ";mpi_rank=",
        ctx.mpi_rank,
        ";mpi_size=",
        ctx.mpi_size,
    )
end

const MPI_CHECKPOINT_SCHEMA_VERSION = 1
const _MPI_CHECKPOINT_MANIFEST = "mpi_checkpoint_manifest.ser"

_mpi_checkpoint_manifest_path(dir::AbstractString) = joinpath(dir, _MPI_CHECKPOINT_MANIFEST)

function _mpi_checkpoint_rank_filename(logical_rank::Integer)
    logical_rank >= 1 || throw(ArgumentError("logical rank must be positive, got $logical_rank"))
    return string("rank_", lpad(string(Int(logical_rank)), 6, '0'), ".ser")
end

_mpi_checkpoint_rank_path(dir::AbstractString, logical_rank::Integer) =
    joinpath(dir, _mpi_checkpoint_rank_filename(logical_rank))

function _mpi_checkpoint_manifest(
    st::HybridStepper{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    rank_files = Tuple(_mpi_checkpoint_rank_filename(r) for r = 1:nranks(ctx.layout))
    return (
        schema = MPI_CHECKPOINT_SCHEMA_VERSION,
        D = D,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        ranks = ctx.layout.ranks,
        periodic = ctx.layout.periodic,
        nranks = nranks(ctx.layout),
        time = st.time[],
        step = st.step[],
        rank_files = rank_files,
    )
end

function _validate_mpi_checkpoint_filename(file)
    file isa AbstractString ||
        throw(ArgumentError("MPI checkpoint rank file names must be strings"))
    s = String(file)
    !isempty(s) || throw(ArgumentError("MPI checkpoint rank file name must be nonempty"))
    basename(s) == s || throw(ArgumentError("MPI checkpoint rank file must be relative: $s"))
    return s
end

function _validate_mpi_checkpoint_manifest(
    manifest,
    st::HybridStepper{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    manifest isa NamedTuple ||
        throw(ArgumentError("MPI checkpoint manifest is not a valid container"))
    required = (:schema, :D, :T, :ncell, :L, :ranks, :periodic, :nranks, :time, :step, :rank_files)
    missing = Symbol[k for k in required if !hasproperty(manifest, k)]
    if !isempty(missing)
        throw(
            ArgumentError(
                "MPI checkpoint manifest is missing fields: $(join(string.(missing), ", "))",
            ),
        )
    end
    manifest.schema == MPI_CHECKPOINT_SCHEMA_VERSION || throw(
        ArgumentError(
            "MPI checkpoint schema $(manifest.schema) ≠ $(MPI_CHECKPOINT_SCHEMA_VERSION)",
        ),
    )
    manifest.D == D || throw(ArgumentError("MPI checkpoint dimension $(manifest.D) ≠ $D"))
    manifest.T == T || throw(ArgumentError("MPI checkpoint eltype $(manifest.T) ≠ $T"))
    manifest.ncell == st.g.n ||
        throw(ArgumentError("MPI checkpoint grid $(manifest.ncell) ≠ $(st.g.n)"))
    manifest.L == st.g.L ||
        throw(ArgumentError("MPI checkpoint box lengths $(manifest.L) ≠ $(st.g.L)"))
    manifest.ranks == ctx.layout.ranks ||
        throw(ArgumentError("MPI checkpoint rank grid $(manifest.ranks) ≠ $(ctx.layout.ranks)"))
    manifest.periodic == ctx.layout.periodic || throw(
        ArgumentError("MPI checkpoint periodicity $(manifest.periodic) ≠ $(ctx.layout.periodic)"),
    )
    manifest.nranks == nranks(ctx.layout) ||
        throw(ArgumentError("MPI checkpoint rank count $(manifest.nranks) ≠ $(nranks(ctx.layout))"))
    length(manifest.rank_files) == nranks(ctx.layout) || throw(
        ArgumentError(
            "MPI checkpoint manifest has $(length(manifest.rank_files)) rank files, expected $(nranks(ctx.layout))",
        ),
    )
    foreach(_validate_mpi_checkpoint_filename, manifest.rank_files)
    return nothing
end

function _mpi_checkpoint_state(
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    base = _checkpoint_state(st, ps)
    _validate_checkpoint_state(base, st)
    return merge(
        base,
        (
            mpi_schema = MPI_CHECKPOINT_SCHEMA_VERSION,
            logical_rank = ctx.logical_rank,
            mpi_rank = ctx.mpi_rank,
            mpi_size = ctx.mpi_size,
            layout_ranks = ctx.layout.ranks,
            layout_periodic = ctx.layout.periodic,
        ),
    )
end

function _validate_mpi_rank_checkpoint_state(
    state,
    manifest,
    st::HybridStepper{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    _validate_checkpoint_state(state, st)
    required = (:mpi_schema, :logical_rank, :mpi_rank, :mpi_size, :layout_ranks, :layout_periodic)
    missing = Symbol[k for k in required if !hasproperty(state, k)]
    if !isempty(missing)
        throw(
            ArgumentError("MPI rank checkpoint is missing fields: $(join(string.(missing), ", "))"),
        )
    end
    state.mpi_schema == MPI_CHECKPOINT_SCHEMA_VERSION || throw(
        ArgumentError(
            "MPI rank checkpoint schema $(state.mpi_schema) ≠ $(MPI_CHECKPOINT_SCHEMA_VERSION)",
        ),
    )
    state.logical_rank == ctx.logical_rank || throw(
        ArgumentError(
            "MPI rank checkpoint logical rank $(state.logical_rank) ≠ $(ctx.logical_rank)",
        ),
    )
    state.mpi_size == ctx.mpi_size ||
        throw(ArgumentError("MPI rank checkpoint size $(state.mpi_size) ≠ $(ctx.mpi_size)"))
    state.layout_ranks == ctx.layout.ranks || throw(
        ArgumentError("MPI rank checkpoint layout $(state.layout_ranks) ≠ $(ctx.layout.ranks)"),
    )
    state.layout_periodic == ctx.layout.periodic || throw(
        ArgumentError(
            "MPI rank checkpoint periodicity $(state.layout_periodic) ≠ $(ctx.layout.periodic)",
        ),
    )
    state.time == manifest.time ||
        throw(ArgumentError("MPI rank checkpoint time $(state.time) ≠ manifest $(manifest.time)"))
    state.step == manifest.step ||
        throw(ArgumentError("MPI rank checkpoint step $(state.step) ≠ manifest $(manifest.step)"))
    return nothing
end

"""
    save_mpi_checkpoint(dir, stepper, ps, ctx) -> manifest_path

Collectively write a restartable MPI checkpoint directory. Each rank writes its
rank-local particles plus replicated fields to a logical-rank file, and MPI rank
0 writes a manifest describing the grid, layout, time, step, and rank-file map.
All ranks return the shared manifest path after the files are complete.
"""
function save_mpi_checkpoint(
    dir::AbstractString,
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    ensure_mpi_initialized!()
    ctx.mpi_rank == 0 && mkpath(dir)
    MPI.Barrier(ctx.comm)

    rank_path = _mpi_checkpoint_rank_path(dir, ctx.logical_rank)
    Serialization.serialize(rank_path, _mpi_checkpoint_state(st, ps, ctx))
    MPI.Barrier(ctx.comm)

    manifest_path = _mpi_checkpoint_manifest_path(dir)
    if ctx.mpi_rank == 0
        manifest = _mpi_checkpoint_manifest(st, ctx)
        Serialization.serialize(manifest_path, manifest)
    end
    MPI.Barrier(ctx.comm)
    return manifest_path
end

"""
    load_mpi_checkpoint!(stepper, ps, dir, ctx) -> stepper

Collectively restore the rank-local state written by [`save_mpi_checkpoint`](@ref).
The manifest and rank-local file are validated against `ctx`, `stepper`, grid
identity, time, step, and logical-rank identity before the destination state is
mutated. After loading, [`mpi_step!`](@ref) can continue the run.
"""
function load_mpi_checkpoint!(
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    dir::AbstractString,
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    ensure_mpi_initialized!()
    manifest = Serialization.deserialize(_mpi_checkpoint_manifest_path(dir))
    _validate_mpi_checkpoint_manifest(manifest, st, ctx)
    rank_file = _validate_mpi_checkpoint_filename(manifest.rank_files[ctx.logical_rank])
    state = Serialization.deserialize(joinpath(dir, rank_file))
    _validate_mpi_rank_checkpoint_state(state, manifest, st, ctx)
    return _restore_checkpoint_state!(st, ps, state)
end

function _mpi_reduction_op(op::Symbol)
    op === :sum && return +
    op === :min && return min
    op === :max && return max
    throw(ArgumentError("MPI diagnostic reduction op must be :sum, :min, or :max, got $op"))
end

function _gpu_array_backend(value)
    type_name = string(typeof(value))
    occursin("CUDA.CuArray", type_name) && return :cuda
    occursin("AMDGPU.ROCArray", type_name) && return :rocm
    occursin("Metal.MtlArray", type_name) && return :metal
    return nothing
end

"""
    mpi_buffer_uses_host_staging(buffer; status=gpu_aware_mpi_status()) -> Bool

Return whether `buffer` should be copied through host memory before passing it
to MPI. Plain Julia `Array`s are direct host buffers. Known CUDA/ROCm arrays are
direct only when MPI reports matching GPU-aware support. Metal and unknown
`AbstractArray` implementations are staged conservatively.
"""
function mpi_buffer_uses_host_staging(
    buffer::AbstractArray;
    status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
)
    buffer isa Array && return false
    backend = _gpu_array_backend(buffer)
    backend === :cuda && return !status.cuda
    backend === :rocm && return !status.rocm
    return true
end

host_staging_buffer(buffer::Array) = buffer
host_staging_buffer(buffer::AbstractArray) = Array(buffer)

function copy_from_host_staging!(dest::AbstractArray, staging::AbstractArray)
    axes(dest) == axes(staging) || throw(
        DimensionMismatch(
            "destination axes $(axes(dest)) do not match staging axes $(axes(staging))",
        ),
    )
    copyto!(dest, staging)
    return dest
end

"""
    MPIBufferPlan

Buffer plan returned by [`prepare_mpi_buffer`](@ref). `buffer` is the object that
can be passed to MPI. `original` is the caller-provided array. `used_host_staging`
records whether `buffer` is a host copy, and `copy_back` records whether
[`finish_mpi_buffer!`](@ref) should write staged receive data back to `original`.
"""
struct MPIBufferPlan{B,O}
    buffer::B
    original::O
    used_host_staging::Bool
    copy_back::Bool
end

function _validate_mpi_buffer_intent(intent::Symbol)
    intent in (:send, :recv, :sendrecv, :inplace) || throw(
        ArgumentError(
            "MPI buffer intent must be :send, :recv, :sendrecv, or :inplace, got $intent",
        ),
    )
    return intent
end

"""
    prepare_mpi_buffer(buffer; status=gpu_aware_mpi_status(), intent=:send)

Prepare an `AbstractArray` for MPI transport. If direct transport is not safe,
the returned plan contains a host `Array` staging buffer. Use
[`finish_mpi_buffer!`](@ref) after receive, sendrecv, or in-place operations to
copy staged data back into the original array.
"""
function prepare_mpi_buffer(
    buffer::AbstractArray;
    status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
    intent::Symbol = :send,
)
    intent = _validate_mpi_buffer_intent(intent)
    if mpi_buffer_uses_host_staging(buffer; status)
        staging = host_staging_buffer(buffer)
        return MPIBufferPlan(staging, buffer, true, intent !== :send)
    end
    return MPIBufferPlan(buffer, buffer, false, false)
end

function finish_mpi_buffer!(plan::MPIBufferPlan)
    plan.used_host_staging && plan.copy_back && copy_from_host_staging!(plan.original, plan.buffer)
    return plan.original
end

function _mpi_allreduce_value(
    value::NamedTuple{names},
    op,
    comm::MPI.Comm,
    status::GPUAwareMPIStatus,
) where {names}
    reduced = ntuple(
        i -> _mpi_allreduce_value(getfield(value, names[i]), op, comm, status),
        Val(length(names)),
    )
    return NamedTuple{names}(reduced)
end

function _mpi_allreduce_value(value::Tuple, op, comm::MPI.Comm, status::GPUAwareMPIStatus)
    return ntuple(i -> _mpi_allreduce_value(value[i], op, comm, status), length(value))
end

function _mpi_allreduce_array(value::AbstractArray, op, comm::MPI.Comm)
    reduced = similar(value)
    MPI.Allreduce!(MPI.RBuffer(value, reduced), op, comm)
    return reduced
end

function _mpi_allreduce_value(value::AbstractArray, op, comm::MPI.Comm, status::GPUAwareMPIStatus)
    plan = prepare_mpi_buffer(value; status, intent = :send)
    return _mpi_allreduce_array(plan.buffer, op, comm)
end

_mpi_allreduce_value(value::Number, op, comm::MPI.Comm, status::GPUAwareMPIStatus) =
    MPI.Allreduce(value, op, comm)

function _mpi_allreduce_value(value, op, comm::MPI.Comm, status::GPUAwareMPIStatus)
    throw(ArgumentError("unsupported MPI diagnostic reduction leaf type $(typeof(value))"))
end

"""
    mpi_allreduce_diagnostics(local_value, ctx; op=:sum, gpu_status=gpu_aware_mpi_status())

Allreduce a rank-local diagnostic value over `ctx.comm`. Supports the same leaf
types as [`reduce_diagnostics`](@ref): numbers, arrays, tuples, and named
tuples. This is real MPI transport; one-rank tests exercise it through
`MPI.COMM_SELF`, and multi-rank tests exercise two, four, and eight local ranks.
Array leaves are routed through the host-staging fallback when GPU-aware MPI is
not available.
"""
function mpi_allreduce_diagnostics(
    local_value,
    ctx::MPICartesianCommunicator;
    op::Symbol = :sum,
    gpu_status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
)
    ensure_mpi_initialized!()
    return _mpi_allreduce_value(local_value, _mpi_reduction_op(op), ctx.comm, gpu_status)
end

function _mpi_moments!(
    nout::Array{T,D},
    uout::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor,
    ctx::MPICartesianCommunicator{D};
    work::Union{Nothing,Vector{T}} = nothing,
    gpu_status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
) where {D,T}
    nf = _require_finite_positive_real("nfloor", nfloor, T)
    density!(nout, ps, g, shape)
    momentum!(uout, ps, g, shape; work)

    copyto!(nout, mpi_allreduce_diagnostics(nout, ctx; op = :sum, gpu_status))
    for c = 1:3
        copyto!(uout[c], mpi_allreduce_diagnostics(uout[c], ctx; op = :sum, gpu_status))
    end

    @inbounds for I in eachindex(nout)
        inv = one(T) / max(nout[I], nf)
        uout[1][I] *= inv
        uout[2][I] *= inv
        uout[3][I] *= inv
    end
    return nout
end

"""
    mpi_compute_moments!(f, ps, g, shape, nfloor, ctx; work=..., gpu_status=...)

Compute globally replicated hybrid moments from rank-local particles. Each rank
deposits its local density and momentum density on the full grid, then real MPI
`Allreduce` sums the arrays so every rank obtains the same global density and
ion bulk velocity as the serial particle set, up to floating-point reduction
order. This is a correctness reference for MPI agreement; scalable domain-local
moment exchange for slab ghost zones is provided by
[`mpi_exchange_ghost_moments!`](@ref), while cluster scaling data remains a
separate gate.
"""
function mpi_compute_moments!(
    f::HybridFields{D,T},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor,
    ctx::MPICartesianCommunicator{D};
    work::Union{Nothing,Vector{T}} = nothing,
    gpu_status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
) where {D,T}
    _mpi_moments!(f.n, f.ui, ps, g, shape, nfloor, ctx; work, gpu_status)
    return f
end

"""
    mpi_init!(stepper, ps, ctx; gpu_status=...)

MPI reference initialization matching [`init!`](@ref), but with moments summed
over rank-local particle subsets by [`mpi_compute_moments!`](@ref).
"""
function mpi_init!(
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    ctx::MPICartesianCommunicator{D};
    gpu_status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
) where {D,T}
    _resize_hybrid_particle_workspaces!(st, nparticles(ps))
    mpi_compute_moments!(
        st.fields,
        ps,
        st.g,
        st.shape,
        st.model.nfloor,
        ctx;
        work = st.work,
        gpu_status,
    )
    ohms_law!(st.fields, st.model, st.g)
    st.time[] = zero(T)                 # (re)start at t=0; step==0 triggers leapfrog priming
    st.step[] = 0
    return st
end

"""
    mpi_step!(stepper, ps, ctx, dt; NB=1, gpu_status=..., migrate_particles=true)

Advance a replicated-field MPI reference step. Particle push, Ohm's law, and
Faraday use the same kernels as [`step!`](@ref); moment deposition is summed
over all ranks with MPI. If `migrate_particles` is true, the rank-local particle
set is migrated through [`mpi_migrate_particles!`](@ref) after the field update.
"""
function mpi_step!(
    st::HybridStepper{D,T},
    ps::ParticleSet{D,T},
    ctx::MPICartesianCommunicator{D},
    dt::Real;
    NB::Integer = 1,
    gpu_status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
    migrate_particles::Bool = true,
) where {D,T}
    dtT = _validated_step_dt(T, dt, NB; min_NB = 1, name = "mpi_step!")
    _resize_hybrid_particle_workspaces!(st, nparticles(ps))

    g = st.g
    h = dtT / 2
    nf = T(st.model.nfloor)
    lo = ntuple(_ -> zero(T), D)
    hi = g.L

    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)
    qm = ps.q / ps.m
    # prime the leapfrog once: loaded v is physical v^0 → v^{-1/2} for 2nd-order accuracy.
    st.step[] == 0 && _prime_leapfrog!(ps.v, st.Ep, st.Bp, qm, h, nparticles(ps))
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        nx, ny, nz = boris_kick(
            vx[p],
            vy[p],
            vz[p],
            st.Ep[1][p],
            st.Ep[2][p],
            st.Ep[3][p],
            st.Bp[1][p],
            st.Bp[2][p],
            st.Bp[3][p],
            qm,
            dtT,
        )
        vx[p] = nx
        vy[p] = ny
        vz[p] = nz
        for d = 1:D
            st.xmid[d][p] = ps.x[d][p] + h * ps.v[d][p]
            ps.x[d][p] += dtT * ps.v[d][p]
        end
    end
    apply_periodic!(ps, lo, hi)

    psmid = ParticleSet{D,T}(st.xmid, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    apply_periodic!(psmid, lo, hi)
    _mpi_moments!(st.fn, st.fui, psmid, g, st.shape, nf, ctx; work = st.work, gpu_status)

    f = st.fields
    # scalar closures freeze pe/∇pe/1/n; the anisotropic (CGL) closure freezes only 1/n and
    # recomputes ∇·P_e per subcycle stage (in _rk4_B!→_bfield_rhs!) — mirror HybridStepper.
    if is_anisotropic(st.model.closure)
        _ohm_ninv!(f.ninv, st.fn, nf, f.floor_count)
    else
        _ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, st.model.closure, nf, f.floor_count, g)
    end

    hb = dtT / NB
    for _ = 1:NB
        _rk4_B!(st, hb)
    end

    mpi_compute_moments!(st.fields, ps, g, st.shape, nf, ctx; work = st.work, gpu_status)
    ohms_law!(st.fields, st.model, st.g)
    # re-center u_i to integer level n+1 (predictor half-kick from the fresh E^{n+1}, B^{n+1}) so
    # the carried E is 2nd-order accurate, matching HybridStepper (see _recenter_carried_E!). The
    # re-deposit uses the MPI moment routine; particles are still pre-migration (as for step 4).
    gather_vector!(st.Ep, st.fields.E, ps, g, st.shape)
    gather_vector!(st.Bp, st.fields.B, ps, g, st.shape)
    _predict_half_kick!(st.vpred, ps.v, st.Ep, st.Bp, ps.q / ps.m, h, nparticles(ps))
    pspred = ParticleSet{D,T}(ps.x, st.vpred, ps.weight, ps.id, ps.tag, ps.q, ps.m)
    mpi_compute_moments!(st.fields, pspred, g, st.shape, nf, ctx; work = st.work, gpu_status)
    ohms_law!(st.fields, st.model, st.g)

    migrate_particles && mpi_migrate_particles!(ps, g, ctx)

    st.time[] += dtT
    st.step[] += 1
    return st
end

function _serialize_mpi_payload(payload)
    io = IOBuffer()
    Serialization.serialize(io, payload)
    return take!(io)
end

function _deserialize_mpi_payload(::Type{T}, bytes::AbstractVector{UInt8}) where {T}
    payload = Serialization.deserialize(IOBuffer(bytes))
    payload isa T ||
        throw(ArgumentError("MPI particle payload has type $(typeof(payload)); expected $T"))
    return payload
end

function _check_mpi_count(n::Integer, label::AbstractString)
    0 <= n <= typemax(Cint) || throw(ArgumentError("$label exceeds Cint count capacity"))
    return Cint(n)
end

_check_mpi_byte_count(n::Integer) = _check_mpi_count(n, "MPI byte payload")

function _mpi_alltoallv_bytes(send_chunks::AbstractVector{<:AbstractVector{UInt8}}, comm::MPI.Comm)
    ncomm = MPI.Comm_size(comm)
    length(send_chunks) == ncomm || throw(
        ArgumentError(
            "send_chunks length $(length(send_chunks)) must equal communicator size $ncomm",
        ),
    )

    send_counts = Cint[_check_mpi_byte_count(length(chunk)) for chunk in send_chunks]
    recv_counts = MPI.Alltoall(MPI.UBuffer(send_counts, 1), comm)
    send_total = sum(Int, send_counts)
    recv_total = sum(Int, recv_counts)

    sendbuf = Vector{UInt8}(undef, send_total)
    offset = 1
    for chunk in send_chunks
        n = length(chunk)
        if n > 0
            copyto!(sendbuf, offset, chunk, firstindex(chunk), n)
            offset += n
        end
    end

    recvbuf = Vector{UInt8}(undef, recv_total)
    MPI.Alltoallv!(MPI.VBuffer(sendbuf, send_counts), MPI.VBuffer(recvbuf, recv_counts), comm)

    chunks = Vector{Vector{UInt8}}(undef, ncomm)
    offset = 1
    for i = 1:ncomm
        n = Int(recv_counts[i])
        chunks[i] = n == 0 ? UInt8[] : recvbuf[offset:(offset+n-1)]
        offset += n
    end
    return chunks
end

function _mpi_rank_for_logical(ctx::MPICartesianCommunicator{D}, logical_rank::Integer) where {D}
    coords = rank_coords(ctx.layout, logical_rank)
    coords0 = [coords[d] - 1 for d = 1:D]
    return MPI.Cart_rank(ctx.comm, coords0)
end

function _validate_mpi_particle_migration_inputs!(
    ps::ParticleSet{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T<:Real}
    local_bad_position = _has_nonfinite_particle_positions(ps)
    any_bad_position = MPI.Allreduce(local_bad_position ? 1 : 0, max, ctx.comm)
    if any_bad_position != 0
        local_bad_position && _validate_migration_positions!(ps)
        throw(ArgumentError("non-finite particle position detected on another MPI rank"))
    end

    q::Float64 = Float64(ps.q)
    m::Float64 = Float64(ps.m)
    q_min::Float64 = MPI.Allreduce(q, min, ctx.comm)
    q_max::Float64 = MPI.Allreduce(q, max, ctx.comm)
    m_min::Float64 = MPI.Allreduce(m, min, ctx.comm)
    m_max::Float64 = MPI.Allreduce(m, max, ctx.comm)
    species_mismatch = q_min != q_max || m_min != m_max
    species_mismatch && throw(
        ArgumentError("MPI particle migration requires identical charge and mass on all ranks"),
    )
    return nothing
end

function _validate_mpi_particle_migration_inputs!(
    ps::ParticleSet{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    throw(ArgumentError("MPI particle migration requires real particle scalar type, got $T"))
end

const _MPI_FIELD_HALO_TAG_UPPER = 0x4841
const _MPI_FIELD_HALO_TAG_LOWER = 0x4842
const _MPI_MOMENT_HALO_TAG_UPPER = 0x4843
const _MPI_MOMENT_HALO_TAG_LOWER = 0x4844
const _MPI_COUNT_TAG_OFFSET = 0x100

function _validate_mpi_slab_array!(
    A::AbstractArray{T,D},
    ctx::MPICartesianCommunicator{D},
    halo::Integer,
) where {T,D}
    ensure_mpi_initialized!()
    MPI.Comm_size(ctx.comm) == nranks(ctx.layout) ||
        throw(ArgumentError("MPI communicator size must equal nranks(ctx.layout)"))
    MPI.Comm_size(ctx.comm) == ctx.mpi_size ||
        throw(ArgumentError("ctx.mpi_size does not match MPI communicator size"))
    h = Int(halo)
    h >= 1 || throw(ArgumentError("halo must be >= 1, got $halo"))
    axis = _slab_axis(ctx.layout)
    axis == 0 && return h, axis
    size(A, axis) >= 3h ||
        throw(DimensionMismatch("axis $axis size $(size(A, axis)) is too small for halo=$halo"))
    return h, axis
end

function _mpi_neighbor_rank(
    ctx::MPICartesianCommunicator{D},
    axis::Integer,
    offset::Integer,
) where {D}
    logical =
        rank_index(ctx.layout, ntuple(d -> d == axis ? ctx.coords[d] + offset : ctx.coords[d], D))
    logical === nothing && return MPI.PROC_NULL
    return _mpi_rank_for_logical(ctx, logical)
end

function _slab_face_ranges(A::AbstractArray, axis::Integer, halo::Integer)
    h = Int(halo)
    nloc = size(A, axis)
    return (;
        lower_ghost = _face_ranges(A, axis, 1:h),
        lower_owned = _face_ranges(A, axis, (h+1):(2h)),
        upper_owned = _face_ranges(A, axis, (nloc-2h+1):(nloc-h)),
        upper_ghost = _face_ranges(A, axis, (nloc-h+1):nloc),
    )
end

function _face_payload(A::AbstractArray{T}, face) where {T}
    return collect(vec(view(A, face...)))
end

function _copy_face_payload!(A::AbstractArray, face, payload::AbstractVector)
    dst = view(A, face...)
    length(dst) == length(payload) || throw(
        DimensionMismatch("received face has $(length(payload)) values; expected $(length(dst))"),
    )
    copyto!(vec(dst), payload)
    return A
end

function _mpi_sendrecv_checked!(
    sendbuf::AbstractVector{T},
    dest::Integer,
    recvbuf::AbstractVector{T},
    source::Integer,
    comm::MPI.Comm,
    tag::Integer,
    label::AbstractString,
) where {T}
    send_count = Cint[_check_mpi_count(length(sendbuf), "$label send count")]
    recv_count = Cint[0]
    MPI.Sendrecv!(
        send_count,
        Int(dest),
        Int(tag + _MPI_COUNT_TAG_OFFSET),
        recv_count,
        Int(source),
        Int(tag + _MPI_COUNT_TAG_OFFSET),
        comm,
    )
    local_mismatch = source != MPI.PROC_NULL && Int(recv_count[1]) != length(recvbuf)
    any_mismatch = MPI.Allreduce(local_mismatch ? 1 : 0, max, comm)
    if any_mismatch != 0
        local_mismatch || throw(
            DimensionMismatch(
                "$label count mismatch detected on another MPI rank; aborting payload exchange",
            ),
        )
        throw(
            DimensionMismatch(
                "$label received $(Int(recv_count[1])) values; expected $(length(recvbuf))",
            ),
        )
    end
    MPI.Sendrecv!(sendbuf, Int(dest), Int(tag), recvbuf, Int(source), Int(tag), comm)
    return recvbuf
end

function _allreduce_exchange_stats(local_a::Integer, local_b::Integer, comm::MPI.Comm)
    return (;
        exchanged = Int(MPI.Allreduce(Int(local_a), +, comm)),
        filled = Int(MPI.Allreduce(Int(local_b), +, comm)),
    )
end

"""
    mpi_exchange_field_halos!(A, ctx; halo=1, fill_value=zero(eltype(A)))

Exchange slab field halos for this rank's local array through real MPI
point-to-point transport. The layout must have at most one decomposed axis.
Lower and upper ghost cells receive neighboring owned boundary cells; exterior
nonperiodic ghosts are filled with `fill_value`. Returns global counts
`(; exchanged, filled)` reduced across `ctx.comm`.
"""
function mpi_exchange_field_halos!(
    A::AbstractArray{T,D},
    ctx::MPICartesianCommunicator{D};
    halo::Integer = 1,
    fill_value = zero(T),
) where {T,D}
    h, axis = _validate_mpi_slab_array!(A, ctx, halo)
    axis == 0 && return (; exchanged = 0, filled = 0)

    faces = _slab_face_ranges(A, axis, h)
    lower = _mpi_neighbor_rank(ctx, axis, -1)
    upper = _mpi_neighbor_rank(ctx, axis, 1)

    local_exchanged = 0
    local_filled = 0

    lower_recv = Vector{T}(undef, length(view(A, faces.lower_ghost...)))
    _mpi_sendrecv_checked!(
        _face_payload(A, faces.upper_owned),
        upper,
        lower_recv,
        lower,
        ctx.comm,
        _MPI_FIELD_HALO_TAG_UPPER,
        "field lower halo",
    )
    if lower == MPI.PROC_NULL
        view(A, faces.lower_ghost...) .= fill_value
        local_filled += length(lower_recv)
    else
        _copy_face_payload!(A, faces.lower_ghost, lower_recv)
        local_exchanged += length(lower_recv)
    end

    upper_recv = Vector{T}(undef, length(view(A, faces.upper_ghost...)))
    _mpi_sendrecv_checked!(
        _face_payload(A, faces.lower_owned),
        lower,
        upper_recv,
        upper,
        ctx.comm,
        _MPI_FIELD_HALO_TAG_LOWER,
        "field upper halo",
    )
    if upper == MPI.PROC_NULL
        view(A, faces.upper_ghost...) .= fill_value
        local_filled += length(upper_recv)
    else
        _copy_face_payload!(A, faces.upper_ghost, upper_recv)
        local_exchanged += length(upper_recv)
    end

    return _allreduce_exchange_stats(local_exchanged, local_filled, ctx.comm)
end

function mpi_exchange_field_halos!(
    fields::NTuple{N,<:AbstractArray{T,D}},
    ctx::MPICartesianCommunicator{D};
    halo::Integer = 1,
    fill_value = zero(T),
) where {N,T,D}
    exchanged = 0
    filled = 0
    for A in fields
        stats = mpi_exchange_field_halos!(A, ctx; halo, fill_value)
        exchanged += stats.exchanged
        filled += stats.filled
    end
    return (; exchanged, filled)
end

"""
    mpi_exchange_ghost_moments!(A, ctx; halo=1, clear_ghosts=true)

Accumulate this rank's slab ghost-zone moment contributions into neighboring
rank interiors through real MPI point-to-point transport. Contributions crossing
nonperiodic exterior boundaries are dropped. When `clear_ghosts` is true, local
ghost zones are zeroed after their values have been copied into send buffers.
Returns global counts `(; exchanged, dropped)` reduced across `ctx.comm`.
"""
function mpi_exchange_ghost_moments!(
    A::AbstractArray{T,D},
    ctx::MPICartesianCommunicator{D};
    halo::Integer = 1,
    clear_ghosts::Bool = true,
) where {T,D}
    h, axis = _validate_mpi_slab_array!(A, ctx, halo)
    axis == 0 && return (; exchanged = 0, dropped = 0)

    faces = _slab_face_ranges(A, axis, h)
    lower = _mpi_neighbor_rank(ctx, axis, -1)
    upper = _mpi_neighbor_rank(ctx, axis, 1)

    local_exchanged = 0
    local_dropped = 0

    lower_contrib = Vector{T}(undef, length(view(A, faces.lower_owned...)))
    _mpi_sendrecv_checked!(
        _face_payload(A, faces.upper_ghost),
        upper,
        lower_contrib,
        lower,
        ctx.comm,
        _MPI_MOMENT_HALO_TAG_UPPER,
        "moment lower interior",
    )
    if lower == MPI.PROC_NULL
        local_dropped += length(view(A, faces.lower_ghost...))
    else
        view(A, faces.lower_owned...) .+=
            reshape(lower_contrib, size(view(A, faces.lower_owned...)))
        local_exchanged += length(lower_contrib)
    end

    upper_contrib = Vector{T}(undef, length(view(A, faces.upper_owned...)))
    _mpi_sendrecv_checked!(
        _face_payload(A, faces.lower_ghost),
        lower,
        upper_contrib,
        upper,
        ctx.comm,
        _MPI_MOMENT_HALO_TAG_LOWER,
        "moment upper interior",
    )
    if upper == MPI.PROC_NULL
        local_dropped += length(view(A, faces.upper_ghost...))
    else
        view(A, faces.upper_owned...) .+=
            reshape(upper_contrib, size(view(A, faces.upper_owned...)))
        local_exchanged += length(upper_contrib)
    end

    if clear_ghosts
        view(A, faces.lower_ghost...) .= zero(T)
        view(A, faces.upper_ghost...) .= zero(T)
    end

    stats = _allreduce_exchange_stats(local_exchanged, local_dropped, ctx.comm)
    return (; exchanged = stats.exchanged, dropped = stats.filled)
end

function mpi_exchange_ghost_moments!(
    moments::NTuple{N,<:AbstractArray{T,D}},
    ctx::MPICartesianCommunicator{D};
    halo::Integer = 1,
    clear_ghosts::Bool = true,
) where {N,T,D}
    exchanged = 0
    dropped = 0
    for A in moments
        stats = mpi_exchange_ghost_moments!(A, ctx; halo, clear_ghosts)
        exchanged += stats.exchanged
        dropped += stats.dropped
    end
    return (; exchanged, dropped)
end

"""
    mpi_migrate_particles!(ps, g, ctx) -> (; moved, lost, sent, received)

Migrate this rank's particles through destination-routed real MPI transport
according to `ctx`'s Cartesian layout. Periodic coordinates are wrapped before
destination classification; particles outside a nonperiodic global boundary are
removed. Each rank serializes only the particles needed by each destination and
exchanges variable-size byte payloads with `MPI.Alltoallv!`, so inbound data is
received only by its owning rank.
"""
function mpi_migrate_particles!(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    ensure_mpi_initialized!()
    MPI.Comm_size(ctx.comm) == nranks(ctx.layout) ||
        throw(ArgumentError("MPI communicator size must equal nranks(ctx.layout)"))
    _validate_mpi_particle_migration_inputs!(ps, ctx)

    _wrap_for_layout!(ps, g, ctx.layout)
    keep = Int[]
    by_dest = Dict{Int,Vector{Int}}()
    lost_local = 0

    @inbounds for p = 1:nparticles(ps)
        pos = ntuple(d -> ps.x[d][p], D)
        dest = rank_of_position(pos, g, ctx.layout)
        if dest === nothing
            lost_local += 1
        elseif dest == ctx.logical_rank
            push!(keep, p)
        else
            push!(get!(by_dest, dest, Int[]), p)
        end
    end

    send_chunks = [UInt8[] for _ = 1:ctx.mpi_size]
    sent_local = 0
    for dest in sort!(collect(keys(by_dest)))
        subset = _subset_particles(ps, by_dest[dest])
        mpi_dest = _mpi_rank_for_logical(ctx, dest)
        send_chunks[mpi_dest+1] = _serialize_mpi_payload(subset)
        sent_local += nparticles(subset)
    end

    _replace_particles!(ps, keep)

    payload_type = ParticleSet{D,T}
    chunks = _mpi_alltoallv_bytes(send_chunks, ctx.comm)
    received_local = 0
    for chunk in chunks
        if !isempty(chunk)
            incoming = _deserialize_mpi_payload(payload_type, chunk)
            append_particles!(ps, incoming)
            received_local += nparticles(incoming)
        end
    end

    moved_global = MPI.Allreduce(sent_local, +, ctx.comm)
    lost_global = MPI.Allreduce(lost_local, +, ctx.comm)
    return (;
        moved = Int(moved_global),
        lost = Int(lost_global),
        sent = sent_local,
        received = received_local,
    )
end
