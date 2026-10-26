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

function _checked_fft_configuration(n, transform, periodic, real_type, comm)
    (real_type === Float32 || real_type === Float64) ||
        throw(ArgumentError("distributed FFT real_type must be Float32 or Float64"))
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
    parameters = PencilFFTs.GlobalFFTParams(nn, transforms, real_type)
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
    return nn, transforms, real_type
end

@inline function _agree_configuration(f, comm, context)
    configuration, encoded =
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "configuration validation") do
            value = f()
            bytes = HybridPlasmaPIC._serialize_mpi_payload(value)
            length(bytes) <= typemax(Cint) || throw(
                ArgumentError(
                    "$(lowercase(context)) configuration exceeds the MPI message-size limit",
                ),
            )
            value, bytes
        end
    count = MPI.Bcast(Cint(length(encoded)), 0, comm)
    reference = MPI.Comm_rank(comm) == 0 ? encoded : Vector{UInt8}(undef, count)
    MPI.Bcast!(reference, 0, comm)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "configuration agreement") do
        agreed = HybridPlasmaPIC._deserialize_mpi_payload(Tuple, reference)
        isequal(configuration, agreed) ||
            throw(ArgumentError("$(lowercase(context)) configuration differs between MPI ranks"))
    end
    return configuration
end

"""
    distributed_fft_plan(n; comm=MPI.COMM_WORLD, transform=Transforms.FFT(), periodic=(true,true,true), real_type=Float64)

Collectively construct a fully periodic three-dimensional distributed transform.
All ranks in `comm` must call with equal grid sizes, equal `real_type` (Float32 or
Float64), and equal transforms after expanding the single-transform shorthand.
Transform configuration must support
serialization and value equality. Empty pencil partitions are rejected.

Recoverable local validation and FFT-plan errors are reported collectively;
the originating rank retains its exception and peers receive an explicit error.
Communicator or MPI transport failures are outside this contract. Plans and data
remain owned by PencilFFTs/PencilArrays; this wrapper does not gather global fields.
"""
@inline function HybridPlasmaPIC.distributed_fft_plan(
    n;
    comm = MPI.COMM_WORLD,
    transform = PencilFFTs.Transforms.FFT(),
    periodic = (true, true, true),
    real_type = Float64,
)
    HybridPlasmaPIC.ensure_mpi_initialized!()
    configuration = _agree_configuration(comm, "Distributed FFT") do
        _checked_fft_configuration(n, transform, periodic, real_type, comm)
    end
    nn, transforms, precision = configuration
    pencil = PencilFFTs.Pencil(nn, comm)
    plan = HybridPlasmaPIC._mpi_collective_stage(comm, "Distributed FFT", "plan construction") do
        PencilFFTs.PencilFFTPlan(pencil, transforms, precision)
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

"""
    distributed_fft_roundtrip_error(plan, input; output=...)

Apply a forward and inverse transform and return the communicator-wide maximum
absolute change in the input. This diagnostic allocates a copy of the input.
All ranks must provide compatible storage for the plan. In-place transforms
reuse `input` as the default output; an explicit output must be that same object.
"""
function HybridPlasmaPIC.distributed_fft_roundtrip_error(
    plan::DistributedPencilFFTPlan,
    input;
    output = PencilFFTs.is_inplace(plan.plan) ? input :
             HybridPlasmaPIC.distributed_fft_output(plan),
)
    input_view = PencilFFTs.is_inplace(plan.plan) ? first(input) : input
    original = copy(parent(input_view))
    HybridPlasmaPIC.distributed_fft_forward!(output, plan, input)
    HybridPlasmaPIC.distributed_fft_inverse!(input, plan, output)
    local_error = isempty(original) ? 0.0 : maximum(abs, parent(input_view) .- original)
    return MPI.Allreduce(local_error, MPI.MAX, plan.comm)
end

include("pencil_fourier.jl")

end # module
