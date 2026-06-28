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

function _checked_fft_shape(n::NTuple{3,<:Integer})
    nn = ntuple(d -> Int(n[d]), 3)
    all(>(0), nn) || throw(ArgumentError("distributed FFT grid sizes must be positive, got $nn"))
    return nn
end

function _checked_fft_shape(n::Tuple)
    length(n) == 3 ||
        throw(ArgumentError("distributed FFT requires exactly three grid sizes, got $(length(n))"))
    all(x -> x isa Integer, n) ||
        throw(ArgumentError("distributed FFT grid sizes must be integers, got $n"))
    return _checked_fft_shape((n[1], n[2], n[3]))
end

function HybridPlasmaPIC.distributed_fft_plan(
    n::Tuple;
    comm = MPI.COMM_WORLD,
    transform = PencilFFTs.Transforms.FFT(),
    periodic::NTuple{3,Bool} = (true, true, true),
)
    all(periodic) || throw(
        ArgumentError("distributed FFT support is currently limited to fully periodic 3-D domains"),
    )
    nn = _checked_fft_shape(n)
    HybridPlasmaPIC.ensure_mpi_initialized!()
    pencil = PencilFFTs.Pencil(nn, comm)
    plan = PencilFFTs.PencilFFTPlan(pencil, transform)
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
