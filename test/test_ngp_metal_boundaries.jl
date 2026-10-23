module NGPMetalBoundaryTests
using HybridPlasmaPIC, Test

# The package's ordinary CI has no Metal dependency/device. The native Metal
# acceptance command requires a functional device and runs these same kernels.
if Sys.isapple() && Base.find_package("Metal") !== nothing
    import Metal
    if Metal.functional()
        prepare_gpu_backend!(Val(:metal))
        @testset "Metal shape boundaries share CPU coupling policy" begin
            for D = 1:3, shape in (NGP(), CIC(), TSC())
                n = ntuple(_ -> 4, D)
                g = FourierGrid(n, ntuple(_ -> 4.0f0, D))
                ps = ParticleSet{D,Float32}(1)
                field = reshape(Float32.(1:prod(n)), n)
                device_field = Metal.MtlArray(field)
                host_out = zeros(Float32, 1)
                device_out = Metal.MtlArray(host_out)
                host_deposit = zeros(Float32, n)
                device_deposit = Metal.MtlArray(host_deposit)
                values = Float32[1]
                device_values = Metal.MtlArray(values)
                for center = 0:3, neighbor in (-1, 0, 1)
                    tie = Float32(center) + 0.5f0
                    s = neighbor < 0 ? prevfloat(tie) : neighbor > 0 ? nextfloat(tie) : tie
                    for x in ps.x
                        x[1] = s
                    end
                    gpu = copy_particles_to_backend(Val(:metal), ps)
                    gather_scalar!(host_out, field, ps, g, shape)
                    gather_scalar!(device_out, device_field, gpu, g, shape)
                    deposit_scalar!(host_deposit, ps, values, g, shape)
                    deposit_scalar!(device_deposit, gpu, device_values, g, shape)
                    Metal.synchronize()
                    gathered, deposited = Array(device_out), Array(device_deposit)
                    @test gathered ≈ host_out rtol = 8eps(Float32) atol = 0
                    @test deposited ≈ host_deposit rtol = 8eps(Float32) atol = 0
                    @test Array(device_field) == field
                    @test Array(device_values) == values
                    if shape isa NGP
                        node = floor(Int, Rational{BigInt}(s) + 1 // 2)
                        index = CartesianIndex(ntuple(_ -> mod(node, 4) + 1, D))
                        expected = zeros(Float32, n)
                        expected[index] = 1.0f0
                        @test deposited == expected
                        @test gathered == Float32[field[index]]
                        # Exact ties must translate by one node on every axis.
                        if neighbor == 0
                            gpu.x[1] .+= 1.0f0
                            deposit_scalar!(device_deposit, gpu, device_values, g, shape)
                            Metal.synchronize()
                            shift = ntuple(d -> d == 1 ? 1 : 0, D)
                            @test Array(device_deposit) == circshift(deposited, shift)
                        end
                    end
                end
            end
        end
    else
        @test_skip "NGP Metal boundaries require a functional Metal device"
    end
end
end
