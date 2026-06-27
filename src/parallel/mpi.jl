# mpi.jl (§5.3 parallel/)
#
# Thin MPI.jl transport hooks over the deterministic logical-rank reference
# kernels in domain_decomposition.jl. These wrappers are intentionally small:
# they create real MPI Cartesian communicators, expose rank/layout provenance,
# and provide structured collective diagnostics plus correctness-first particle
# migration. Time-advanced rank invariance and scaling still require
# mpiexec/cluster tests.

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
`MPI.COMM_SELF`, while multi-rank invariance remains a separate gate. Array
leaves are routed through the host-staging fallback when GPU-aware MPI is not
available.
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

function _resize_hybrid_particle_workspaces!(st::HybridStepper{D,T}, n::Integer) where {D,T}
    n >= 0 || throw(ArgumentError("particle workspace length must be nonnegative, got $n"))
    N = Int(n)
    length(st.work) == N &&
        all(length(st.Ep[c]) == N for c = 1:3) &&
        all(length(st.Bp[c]) == N for c = 1:3) &&
        all(length(st.xmid[d]) == N for d = 1:D) &&
        return st

    st.Ep = ntuple(_ -> Vector{T}(undef, N), 3)
    st.Bp = ntuple(_ -> Vector{T}(undef, N), 3)
    st.xmid = ntuple(_ -> Vector{T}(undef, N), D)
    st.work = Vector{T}(undef, N)
    return st
end

function _mpi_moments!(
    nout::Array{T,D},
    uout::NTuple{3,<:Array{T,D}},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor,
    ctx::MPICartesianCommunicator{D};
    work::Vector{T} = Vector{T}(undef, nparticles(ps)),
    gpu_status::GPUAwareMPIStatus = gpu_aware_mpi_status(),
) where {D,T}
    density!(nout, ps, g, shape)
    momentum!(uout, ps, g, shape; work)

    copyto!(nout, mpi_allreduce_diagnostics(nout, ctx; op = :sum, gpu_status))
    for c = 1:3
        copyto!(uout[c], mpi_allreduce_diagnostics(uout[c], ctx; op = :sum, gpu_status))
    end

    nf = T(nfloor)
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
moment exchange remains a separate gate.
"""
function mpi_compute_moments!(
    f::HybridFields{D,T},
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    shape::ShapeFunction,
    nfloor,
    ctx::MPICartesianCommunicator{D};
    work::Vector{T} = Vector{T}(undef, nparticles(ps)),
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
    _ohm_prep!(f.pe, f.gradp, f.ninv, st.fn, st.model.closure, nf, f.floor_count, g)

    hb = dtT / NB
    for _ = 1:NB
        _rk4_B!(st, hb)
    end

    mpi_compute_moments!(st.fields, ps, g, st.shape, nf, ctx; work = st.work, gpu_status)
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

function _mpi_allgather_bytes(bytes::Vector{UInt8}, comm::MPI.Comm)
    length(bytes) <= typemax(Cint) ||
        throw(ArgumentError("MPI byte payload exceeds Cint count capacity"))
    counts = Cint.(MPI.Allgather(Cint(length(bytes)), comm))
    total = sum(Int, counts)
    recv = Vector{UInt8}(undef, total)
    MPI.Allgatherv!(bytes, MPI.VBuffer(recv, counts), comm)

    chunks = Vector{Vector{UInt8}}(undef, length(counts))
    offset = 1
    for i = 1:length(counts)
        n = Int(counts[i])
        chunks[i] = n == 0 ? UInt8[] : recv[offset:(offset+n-1)]
        offset += n
    end
    return chunks
end

"""
    mpi_migrate_particles!(ps, g, ctx) -> (; moved, lost, sent, received)

Migrate this rank's particles through real MPI transport according to `ctx`'s
Cartesian layout. Periodic coordinates are wrapped before destination
classification; particles outside a nonperiodic global boundary are removed.

This is a correctness-first collective transport: each rank serializes its
variable-size outbound particle subsets and exchanges them with `Allgatherv`.
It preserves the same destination rules as [`migrate_particles!`](@ref), but it
is not the scalable production neighbor-exchange path. Scaling remains a
separate Phase-9 gate.
"""
function mpi_migrate_particles!(
    ps::ParticleSet{D,T},
    g::FourierGrid{D,T},
    ctx::MPICartesianCommunicator{D},
) where {D,T}
    ensure_mpi_initialized!()
    MPI.Comm_size(ctx.comm) == nranks(ctx.layout) ||
        throw(ArgumentError("MPI communicator size must equal nranks(ctx.layout)"))

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

    outbound = Tuple{Int,ParticleSet{D,T}}[]
    sent_local = 0
    for dest in sort!(collect(keys(by_dest)))
        subset = _subset_particles(ps, by_dest[dest])
        push!(outbound, (dest, subset))
        sent_local += nparticles(subset)
    end

    _replace_particles!(ps, keep)

    payload_type = Vector{Tuple{Int,ParticleSet{D,T}}}
    chunks = _mpi_allgather_bytes(_serialize_mpi_payload(outbound), ctx.comm)
    received_local = 0
    for chunk in chunks
        rank_payload = _deserialize_mpi_payload(payload_type, chunk)
        for (dest, incoming) in rank_payload
            if dest == ctx.logical_rank
                append_particles!(ps, incoming)
                received_local += nparticles(incoming)
            end
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
