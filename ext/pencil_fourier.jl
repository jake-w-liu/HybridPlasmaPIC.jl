const SpectralOperators = HybridPlasmaPIC.SpectralOperators

struct DistributedFourierGrid{T,F,A,B,R}
    n::NTuple{3,Int}
    L::NTuple{3,T}
    dx::NTuple{3,T}
    midx::NTuple{3,Vector{Int}}
    kfull::NTuple{3,Vector{T}}
    ik::NTuple{3,Vector{Complex{T}}}
    kvec::NTuple{3,Vector{T}}
    gain::NTuple{3,T}
    count::Int
    fft::F
    physical::A
    spectral::B
    staging::R
end

"""
    distributed_fourier_grid(n, L; comm=MPI.COMM_WORLD)

Collectively construct a periodic 3D Fourier workspace with rank-owned fields.
Every rank must supply the same three grid sizes and positive finite lengths.
Lengths have one common Float32 or Float64 type, which selects field precision.
The workspace belongs to one caller at a time; its metadata vectors are read-only.
The communicator remains caller-owned. This workspace provides field derivatives;
it does not provide particle halos, distributed stepping or checkpoint storage.
"""
function HybridPlasmaPIC.distributed_fourier_grid(n, L; comm = MPI.COMM_WORLD)
    HybridPlasmaPIC.ensure_mpi_initialized!()
    nn, lengths = _agree_configuration(comm, "Distributed Fourier grid") do
        nn = _checked_fft_shape(n)
        (L isa NTuple{3,Float32} || L isa NTuple{3,Float64}) || throw(
            ArgumentError(
                "distributed Fourier lengths must be a three-tuple of Float32 or Float64",
            ),
        )
        all(x -> isfinite(x) && x > zero(x), L) ||
            throw(ArgumentError("distributed Fourier lengths must be positive and finite"))
        nn, L
    end
    metadata, count = HybridPlasmaPIC._mpi_collective_stage(
        comm,
        "Distributed Fourier grid",
        "metadata construction",
    ) do
        count = Base.checked_mul(Base.checked_mul(nn[1], nn[2]), nn[3])
        SpectralOperators._fourier_metadata(nn, lengths), count
    end
    T = eltype(lengths)
    fft = HybridPlasmaPIC.distributed_fft_plan(nn; comm, real_type = T)
    return HybridPlasmaPIC._mpi_collective_stage(
        comm,
        "Distributed Fourier grid",
        "workspace allocation",
    ) do
        physical = HybridPlasmaPIC.distributed_fft_input(fft)
        spectral = HybridPlasmaPIC.distributed_fft_output(fft)
        staging = similar(physical, T)
        gain = ntuple(d -> maximum(abs, metadata.kvec[d]), 3)
        DistributedFourierGrid(
            metadata.n,
            metadata.L,
            metadata.dx,
            metadata.midx,
            metadata.kfull,
            metadata.ik,
            metadata.kvec,
            gain,
            count,
            fft,
            physical,
            spectral,
            staging,
        )
    end
end

"""
    distributed_field(grid)

Allocate a zeroed real CPU PencilArray on `grid`'s exact physical pencil.
Use `PencilArrays.global_view(field)` to access global logical cell indices.
Allocation is local; every collective derivative still requires all grid ranks.
"""
function HybridPlasmaPIC.distributed_field(g::DistributedFourierGrid)
    return fill!(similar(g.staging), zero(eltype(g.staging)))
end

function _require_physical_field(name, a, g::DistributedFourierGrid{T}) where {T}
    a isa PencilArrays.PencilArray{T,3} ||
        throw(ArgumentError("$name must be a three-dimensional PencilArray with element type $T"))
    parent(a) isa Array{T,3} || throw(ArgumentError("$name must use CPU Array storage"))
    PencilArrays.pencil(a) === PencilArrays.pencil(g.physical) ||
        throw(ArgumentError("$name must use this grid's physical pencil"))
    axes(a) == axes(g.physical) ||
        throw(DimensionMismatch("$name axes do not match the grid's owned cells"))
    for workspace in (g.physical, g.spectral, g.staging)
        Base.mightalias(a, workspace) &&
            throw(ArgumentError("$name must not alias private Fourier workspace storage"))
    end
    return nothing
end

"""
    deriv!(out, field, grid::DistributedFourierGrid, axis)

Collectively differentiate a real rank-owned field along the agreed global axis.
Input and output may alias each other. Both must use the grid's physical pencil,
precision and CPU Array storage, and must not alias its private workspaces.
Recoverable local failures are reported on every rank before publishing output.
MPI process/transport failures and fatal resource failures are excluded.
"""
function SpectralOperators.deriv!(out, field, g::DistributedFourierGrid{T}, axis) where {T}
    comm = g.fft.comm
    j::Int, local_magnitude::Union{Float32,Float64} =
        HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "input validation") do
            axis isa Integer && 1 <= axis <= 3 ||
                throw(ArgumentError("distributed derivative axis must be an integer in 1:3"))
            _require_physical_field("input", field, g)
            _require_physical_field("output", out, g)
            magnitude = maximum(abs, parent(field))
            isfinite(magnitude) ||
                throw(ArgumentError("spectral operator input values must be finite"))
            Int(axis), magnitude
        end
    agreed_axis = MPI.Bcast(j, 0, comm)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "axis agreement") do
        j == agreed_axis ||
            throw(ArgumentError("distributed derivative axis differs between MPI ranks"))
    end
    magnitude = MPI.Allreduce(local_magnitude, MPI.MAX, comm)
    safe = SpectralOperators._fft_safe_magnitude(T, g.count, g.gain[j])
    scale = SpectralOperators._fft_input_scale(magnitude, safe)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "input scaling") do
        SpectralOperators._copy_fft_input!(parent(g.physical), parent(field), scale)
    end
    HybridPlasmaPIC.distributed_fft_forward!(g.spectral, g.fft, g.physical)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "mode multiplication") do
        modes = PencilArrays.global_view(g.spectral)
        SpectralOperators._apply_ik_store!(modes, modes, g.ik[j], j)
    end
    HybridPlasmaPIC.distributed_fft_inverse!(g.physical, g.fft, g.spectral)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed derivative", "output preparation") do
        SpectralOperators._store_fft_output!(parent(g.staging), parent(g.physical), scale)
    end
    copyto!(parent(out), parent(g.staging))
    return out
end

function SpectralOperators.deriv(field, g::DistributedFourierGrid, axis)
    out = HybridPlasmaPIC.distributed_field(g)
    return SpectralOperators.deriv!(out, field, g, axis)
end
