using HybridPlasmaPIC, FFTW, MPI, Test
import PencilArrays
import PencilFFTs

@testset "PencilFFTs distributed FFT extension" begin
    ensure_mpi_initialized!()
    @test extension_loaded(Val(:pencilfft))
    @test extension_dependency_module(Val(:pencilfft)) == (PencilFFTs, PencilArrays)

    @test_throws ArgumentError distributed_fft_plan((8, 6); comm = MPI.COMM_SELF)
    @test_throws ArgumentError distributed_fft_plan((8, 6, 0); comm = MPI.COMM_SELF)
    @test_throws ArgumentError distributed_fft_plan(
        (8, 6, 4);
        comm = MPI.COMM_SELF,
        periodic = (true, false, true),
    )

    plan = distributed_fft_plan((8, 6, 4); comm = MPI.COMM_SELF)
    input = distributed_fft_input(plan)
    output = distributed_fft_output(plan)
    @test size(input) == (8, 6, 4)
    @test size(output) == (8, 6, 4)

    local_input = parent(input)
    for I in CartesianIndices(local_input)
        i, j, k = Tuple(I)
        local_input[I] = complex(sin(0.2i) + cos(0.3j), 0.1k - 0.2j)
    end
    reference = copy(Array(input))

    distributed_fft_forward!(output, plan, input)
    @test PencilArrays.gather(output) ≈ fft(reference)
    distributed_fft_inverse!(input, plan, output)
    @test Array(input) ≈ reference

    @test distributed_fft_roundtrip_error(plan, input) < 1e-12
end
