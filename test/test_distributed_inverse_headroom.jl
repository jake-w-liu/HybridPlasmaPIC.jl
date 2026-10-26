using HybridPlasmaPIC, PencilArrays, PencilFFTs, MPI, Test

# Independent physical-space reference, without transforms or shared metadata.
function distributed_cutoff_response(n)
    kernel = setprecision(BigFloat, 256) do
        cutoff = (n - 1) ÷ 3
        [Float64(sum(cospi(BigFloat(2m) * j / n) for m = -cutoff:cutoff) / n) for j = 0:n-1]
    end
    signs = sign.(kernel)
    response = [sum(kernel[mod(x - j, n)+1] * signs[j+1] for j = 0:n-1) for x = 0:n-1]
    return signs, response
end

@testset "Distributed inverse FFT range uses complete operation bound" begin
    ensure_mpi_initialized!()
    comm = MPI.COMM_WORLD
    for T in (Float32, Float64), n in ((32, 32, 32), (32, 24, 16))
        axes_reference = map(distributed_cutoff_response, n)
        signs, response = map(first, axes_reference), map(last, axes_reference)
        amplitude = floatmax(T) / T(100)
        retained = prod(2 * ((size - 1) ÷ 3) + 1 for size in n)
        @test Rational{BigInt}(amplitude)^2 * retained < Rational{BigInt}(floatmax(T))^2
        grid = distributed_fourier_grid(n, (T(2), T(3), T(5)); comm)
        field = distributed_field(grid)
        view = PencilArrays.global_view(field)
        for shift in (0, 5)
            local_reference_scale = 0.0
            for I in CartesianIndices(view)
                view[I] = amplitude * T(prod(signs[d][mod(I[d] - 1 - shift, n[d])+1] for d = 1:3))
                reference = prod(response[d][mod(I[d] - 1 - shift, n[d])+1] for d = 1:3)
                local_reference_scale = max(local_reference_scale, abs(reference))
            end
            @test dealias_two_thirds!(field, grid) === field
            @test MPI.Allreduce(count(!isfinite, parent(field)), +, comm) == 0
            local_error = 0.0
            for I in CartesianIndices(view)
                expected = prod(response[d][mod(I[d] - 1 - shift, n[d])+1] for d = 1:3)
                local_error = max(local_error, abs(Float64(view[I] / amplitude) - expected))
            end
            error = MPI.Allreduce(local_error, MPI.MAX, comm)
            reference_scale = MPI.Allreduce(local_reference_scale, MPI.MAX, comm)
            @test error <= 128eps(T) * max(reference_scale, 1.0)
        end
    end
end
