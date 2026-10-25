module DistributedFFTDiagnosticTests
using HybridPlasmaPIC, Test
import MPI, PencilArrays, PencilFFTs, FFTW
const H = HybridPlasmaPIC
const F = PencilFFTs.Transforms
H.ensure_mpi_initialized!()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const SIZE = MPI.Comm_size(COMM)

@testset "Distributed roundtrip diagnostic uses the input view" begin
    unequal_rank_errors = false
    configurations =
        (F.FFT(), F.RFFT(), F.R2R(FFTW.REDFT10), F.FFT!(), F.RFFT!(), F.R2R!(FFTW.REDFT10))
    for n in ((8, 6, 4), (9, 7, 5)), transform in configurations
        plan = distributed_fft_plan(n; comm = COMM, transform)
        inplace = PencilFFTs.is_inplace(plan.plan)
        input = distributed_fft_input(plan)
        output = inplace ? input : distributed_fft_output(plan)
        input_view = inplace ? first(input) : input
        global_input = PencilArrays.global_view(input_view)
        for I in CartesianIndices(global_input)
            value = (RANK + 1) * (sin(0.2I[1]) + cos(0.3I[2]) + 0.7sin(0.4I[3]))
            global_input[I] =
                eltype(global_input) <: Complex ? complex(value, 0.1I[3] - 0.2I[2]) : value
        end
        original = copy(parent(input_view))
        distributed_fft_forward!(output, plan, input)
        distributed_fft_inverse!(input, plan, output)
        restored = copy(parent(input_view))
        local_error = maximum(abs, restored .- original)
        errors = MPI.Allgather(local_error, COMM)
        expected = maximum(errors)
        unequal_rank_errors |= minimum(errors) < expected
        @test expected <= 256eps(Float64) * max(1.0, maximum(abs, original))
        for explicit in (false, true)
            copyto!(parent(input_view), original)
            measured =
                explicit ? distributed_fft_roundtrip_error(plan, input; output) :
                distributed_fft_roundtrip_error(plan, input)
            @test isequal(measured, expected)
            @test isequal(parent(input_view), restored)
        end
        copyto!(parent(input_view), original)
        incompatible = inplace ? distributed_fft_output(plan) : input
        @test_throws ArgumentError distributed_fft_roundtrip_error(
            plan,
            input;
            output = incompatible,
        )
        @test isequal(parent(input_view), original)
    end
    @test SIZE == 1 || unequal_rank_errors
end
end
