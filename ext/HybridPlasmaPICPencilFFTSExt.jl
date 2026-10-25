module HybridPlasmaPICPencilFFTSExt

import HybridPlasmaPIC
import MPI
import PencilArrays
import PencilFFTs
import LinearAlgebra: ldiv!, mul!

struct DistributedPencilFFTPlan{P,C}
    plan::P
    n::NTuple{3,Int}
    comm::C
end

HybridPlasmaPIC.extension_dependency_module(::Val{:pencilfft}) = (PencilFFTs, PencilArrays)

function _checked_fft_shape(n::Tuple)
    length(n) == 3 ||
        throw(ArgumentError("distributed FFT requires exactly three grid sizes, got $(length(n))"))
    all(x -> x isa Integer, n) ||
        throw(ArgumentError("distributed FFT grid sizes must be integers, got $n"))
    return ntuple(d -> HybridPlasmaPIC._require_positive_intlike("n[$d]", n[d]), 3)
end

_checked_fft_shape(n) = throw(ArgumentError("distributed FFT grid sizes must be a tuple"))

function _checked_fft_configuration(n, transform, periodic, comm)
    periodic isa NTuple{3,Bool} ||
        throw(ArgumentError("distributed FFT periodicity must be a three-tuple of Bool values"))
    all(periodic) || throw(
        ArgumentError("distributed FFT support is currently limited to fully periodic 3-D domains"),
    )
    nn = _checked_fft_shape(n)
    transforms = if transform isa PencilFFTs.Transforms.AbstractTransform
        PencilFFTs.Transforms.expand_dims(transform, Val(3))
    elseif transform isa NTuple{3,PencilFFTs.Transforms.AbstractTransform}
        transform
    else
        throw(ArgumentError("distributed FFT requires one transform or a three-transform tuple"))
    end
    # Validate local transform metadata before any topology collective. BigInt
    # lengths avoid wrapping a real-transform output size before checking it.
    output_sizes = ntuple(3) do d
        value = PencilFFTs.Transforms.length_output(transforms[d], BigInt(nn[d]))
        HybridPlasmaPIC._require_positive_intlike("transform output size[$d]", value)
    end
    parameters = PencilFFTs.GlobalFFTParams(nn, transforms)
    data_type = PencilFFTs.input_data_type(parameters)
    for transform in transforms
        data_type = PencilFFTs.Transforms.eltype_output(transform, data_type)
    end
    # PencilFFTs moves the decomposed axes from (2,3) to (1,3), then (1,2).
    # Its backend does not support empty partitions during those transposes.
    ranks = MPI.Dims_create(MPI.Comm_size(comm), (Cint(0), Cint(0)))
    ranks[1] <= min(nn[2], output_sizes[1]) && ranks[2] <= min(nn[3], output_sizes[2]) || throw(
        ArgumentError("distributed FFT process grid $ranks produces an empty pencil partition"),
    )
    return nn, transforms
end

"""
    distributed_fft_plan(n; comm=MPI.COMM_WORLD, transform=Transforms.FFT(), periodic=(true,true,true))

Collectively construct a fully periodic three-dimensional distributed transform.
All ranks in `comm` must call with equal grid sizes and equal transforms after
expanding the single-transform shorthand. Transform configuration must support
serialization and value equality. Empty pencil partitions are rejected.

Recoverable local validation and FFT-plan errors are reported collectively;
the originating rank retains its exception and peers receive an explicit error.
Communicator or MPI transport failures are outside this contract. Plans and data
remain owned by PencilFFTs/PencilArrays; this wrapper does not gather global fields.
"""
function HybridPlasmaPIC.distributed_fft_plan(
    n;
    comm = MPI.COMM_WORLD,
    transform = PencilFFTs.Transforms.FFT(),
    periodic = (true, true, true),
)
    HybridPlasmaPIC.ensure_mpi_initialized!()
    configuration, encoded = HybridPlasmaPIC._mpi_collective_stage(
        comm,
        "Distributed FFT",
        "configuration validation",
    ) do
        value = _checked_fft_configuration(n, transform, periodic, comm)
        bytes = HybridPlasmaPIC._serialize_mpi_payload(value)
        length(bytes) <= typemax(Cint) || throw(
            ArgumentError("distributed FFT configuration exceeds the MPI message-size limit"),
        )
        value, bytes
    end
    # The payload has already been serialized and validated collectively.
    count = MPI.Bcast(Cint(length(encoded)), 0, comm)
    reference = MPI.Comm_rank(comm) == 0 ? encoded : Vector{UInt8}(undef, count)
    MPI.Bcast!(reference, 0, comm)
    HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed FFT", "configuration agreement") do
        agreed = HybridPlasmaPIC._deserialize_mpi_payload(Tuple, reference)
        isequal(configuration, agreed) ||
            throw(ArgumentError("distributed FFT configuration differs between MPI ranks"))
    end
    nn, transforms = configuration
    pencil = PencilFFTs.Pencil(nn, comm)
    plan = HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed FFT", "plan construction") do
        PencilFFTs.PencilFFTPlan(pencil, transforms)
    end
    return DistributedPencilFFTPlan{typeof(plan),typeof(comm)}(plan, nn, comm)
end

HybridPlasmaPIC.distributed_fft_input(plan::DistributedPencilFFTPlan) =
    PencilFFTs.allocate_input(plan.plan)

HybridPlasmaPIC.distributed_fft_output(plan::DistributedPencilFFTPlan) =
    PencilFFTs.allocate_output(plan.plan)

function HybridPlasmaPIC.distributed_fft_forward!(output, plan::DistributedPencilFFTPlan, input)
    mul!(output, plan.plan, input)
    return output
end

function HybridPlasmaPIC.distributed_fft_inverse!(input, plan::DistributedPencilFFTPlan, output)
    ldiv!(input, plan.plan, output)
    return input
end

function HybridPlasmaPIC.distributed_fft_roundtrip_error(
    plan::DistributedPencilFFTPlan,
    input;
    output = HybridPlasmaPIC.distributed_fft_output(plan),
)
    original = copy(parent(input))
    HybridPlasmaPIC.distributed_fft_forward!(output, plan, input)
    HybridPlasmaPIC.distributed_fft_inverse!(input, plan, output)
    local_error = isempty(original) ? 0.0 : maximum(abs, parent(input) .- original)
    return MPI.Allreduce(local_error, MPI.MAX, plan.comm)
end

end # module
