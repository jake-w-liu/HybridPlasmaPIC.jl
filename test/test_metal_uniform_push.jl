module MetalUniformPushTests
using HybridPlasmaPIC, Test

function population(D, count; q = 1.5f0, m = 2.0f0)
    ps = ParticleSet{D,Float32}(count; q, m)
    for p = 1:count
        for d = 1:D
            ps.x[d][p] = Float32(mod(p, 11)) / 8 + Float32(d) / 16
        end
        for c = 1:3
            ps.v[c][p] = Float32(mod(p + 3c, 7) - 3) / 8
        end
        ps.weight[p] = Float32(mod(p, 3) + 1) / 4
        ps.tag[p] = UInt32(mod(p, 2))
    end
    return ps
end
snapshot(ps) =
    (map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(ps.id), copy(ps.tag), ps.q, ps.m)

if Sys.isapple() && Base.find_package("Metal") !== nothing
    import Metal
    if Metal.functional()
        prepare_gpu_backend!(Val(:metal))
        @testset "Metal uniform push covers dimensions and thread groups" begin
            E = (0.125f0, -0.0625f0, 0.03125f0)
            B = (0.25f0, -0.125f0, 0.5f0)
            for D = 1:3, count in (0, 1, 255, 256, 257, 4096)
                cpu = population(D, count)
                gpu = copy_particles_to_backend(Val(:metal), cpu)
                before = snapshot(cpu)
                for steps in (1, 16)
                    for _ = 1:steps
                        push_uniform!(cpu, E, B, 0.03125f0)
                        @test push_uniform!(gpu, E, B, 0.03125f0) === gpu
                    end
                    Metal.synchronize()
                    actual = copy_particles_to_host(gpu)
                    for d = 1:D
                        @test actual.x[d] ≈ cpu.x[d] rtol = 32eps(Float32) atol = 16eps(Float32)
                    end
                    for c = 1:3
                        @test actual.v[c] ≈ cpu.v[c] rtol = 32eps(Float32) atol = 16eps(Float32)
                    end
                    @test isequal(
                        (actual.weight, actual.id, actual.tag, actual.q, actual.m),
                        before[3:7],
                    )
                end
            end
        end

        @testset "Metal free drift and electric acceleration have independent oracles" begin
            for D = 1:3,
                count in (0, 1, 257),
                electric in (false, true),
                dt in (-0.125f0, 0.0f0, 0.125f0)

                cpu = population(D, count; q = 2.0f0, m = 4.0f0)
                gpu = copy_particles_to_backend(Val(:metal), cpu)
                E = electric ? (0.125f0, -0.25f0, 0.5f0) : (0.0f0, 0.0f0, 0.0f0)
                # Every quantity is dyadic and small enough that these complete
                # electric/free-drift expressions are exactly represented.
                expected_v = ntuple(c -> cpu.v[c] .+ 0.5f0 * dt * E[c], 3)
                expected_x = ntuple(d -> cpu.x[d] .+ dt .* expected_v[d], D)
                push_uniform!(gpu, E, (0.0f0, 0.0f0, 0.0f0), dt)
                Metal.synchronize()
                actual = copy_particles_to_host(gpu)
                @test actual.v == expected_v
                @test actual.x == expected_x
            end
        end

        @testset "Metal magnetic rotation retains speed and the discrete gyro angle" begin
            for D = 1:3
                cpu = population(D, 257; q = 1.0f0, m = 1.0f0)
                fill!(cpu.v[1], 1.0f0)
                fill!(cpu.v[2], 0.0f0)
                fill!(cpu.v[3], 0.0f0)
                gpu = copy_particles_to_backend(Val(:metal), cpu)
                dt = 0.25f0
                angle = -2atan(Float64(dt) / 2)
                for step = 1:8
                    push_uniform!(gpu, (0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 1.0f0), dt)
                    Metal.synchronize()
                    actual = copy_particles_to_host(gpu)
                    @test all(
                        v -> isapprox(v, cos(step * angle); atol = 32eps(Float32)),
                        actual.v[1],
                    )
                    @test all(
                        v -> isapprox(v, sin(step * angle); atol = 32eps(Float32)),
                        actual.v[2],
                    )
                    @test all(iszero, actual.v[3])
                    @test all(abs.(actual.v[1] .^ 2 .+ actual.v[2] .^ 2 .- 1) .<= 64eps(Float32))
                end
            end
        end

        @testset "Metal uniform parameter errors precede device mutation" begin
            zero_field = (0.0f0, 0.0f0, 0.0f0)
            invalid = (
                (Inf32, 1.0f0, zero_field, zero_field, 0.25f0),
                (NaN32, 1.0f0, zero_field, zero_field, 0.25f0),
                (1.0f0, 0.0f0, zero_field, zero_field, 0.25f0),
                (1.0f0, -1.0f0, zero_field, zero_field, 0.25f0),
                (1.0f0, Inf32, zero_field, zero_field, 0.25f0),
                (floatmax(Float32), floatmin(Float32), zero_field, zero_field, 0.25f0),
                (1.0f0, 1.0f0, zero_field, zero_field, Inf32),
                (1.0f0, 1.0f0, zero_field, zero_field, NaN32),
                (1.0f0, 1.0f0, zero_field, zero_field, 1.0e40),
                (1.0f0, 1.0f0, (Inf32, 0.0f0, 0.0f0), zero_field, 0.25f0),
                (1.0f0, 1.0f0, (0.0f0, NaN32, 0.0f0), zero_field, 0.25f0),
                (1.0f0, 1.0f0, zero_field, (0.0f0, 0.0f0, Inf32), 0.25f0),
                (1.0f0, 1.0f0, (0.0, 0.0, 1.0e40), zero_field, 0.25f0),
            )
            for count in (0, 2), (q, m, E, B, dt) in invalid
                cpu = population(2, count; q, m)
                gpu = copy_particles_to_backend(Val(:metal), cpu)
                before = snapshot(cpu)
                @test_throws ArgumentError push_uniform!(gpu, E, B, dt)
                Metal.synchronize()
                @test isequal(snapshot(copy_particles_to_host(gpu)), before)
            end
        end
    else
        @test_skip "Metal uniform-push regression requires a functional Metal device"
    end
end
end
