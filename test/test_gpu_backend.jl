using HybridPlasmaPIC, Test

function _fill_backend_test_particles!(::Type{T} = Float64) where {T<:AbstractFloat}
    ps = ParticleSet{2,T}(4; q = T(2), m = T(3))
    ps.x[1] .= T[0.1, 0.2, 0.3, 0.4]
    ps.x[2] .= T[1.1, 1.2, 1.3, 1.4]
    ps.v[1] .= T[2.1, 2.2, 2.3, 2.4]
    ps.v[2] .= T[3.1, 3.2, 3.3, 3.4]
    ps.v[3] .= T[4.1, 4.2, 4.3, 4.4]
    ps.weight .= T[0.5, 0.6, 0.7, 0.8]
    ps.id .= UInt64[11, 12, 13, 14]
    ps.tag .= UInt32[21, 22, 23, 24]
    return ps
end

function _assert_same_particles(a, b)
    @test nparticles(a) == nparticles(b)
    @test a.q == b.q
    @test a.m == b.m
    for d = 1:length(a.x)
        @test collect(a.x[d]) == collect(b.x[d])
    end
    for c = 1:3
        @test collect(a.v[c]) == collect(b.v[c])
    end
    @test collect(a.weight) == collect(b.weight)
    @test collect(a.id) == collect(b.id)
    @test collect(a.tag) == collect(b.tag)
end

function _assert_particles_approx(a, b; rtol, atol)
    @test nparticles(a) == nparticles(b)
    @test a.q == b.q
    @test a.m == b.m
    for d = 1:length(a.x)
        @test collect(a.x[d]) ≈ collect(b.x[d]) rtol = rtol atol = atol
    end
    for c = 1:3
        @test collect(a.v[c]) ≈ collect(b.v[c]) rtol = rtol atol = atol
    end
    @test collect(a.weight) == collect(b.weight)
    @test collect(a.id) == collect(b.id)
    @test collect(a.tag) == collect(b.tag)
end

@testset "backend particle storage helpers" begin
    ps = _fill_backend_test_particles!()
    @test particle_storage_backend(ps) == :cpu
    @test all(
        particle_array_backend(A) == :cpu for A in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
    )
    cpu_memory = backend_memory_status(Val(:cpu))
    @test cpu_memory isa BackendMemoryStatus
    @test cpu_memory.backend == :cpu
    @test cpu_memory.device_available
    @test !cpu_memory.pool_supported
    @test cpu_memory.total_bytes !== nothing
    @test cpu_memory.free_bytes !== nothing
    @test cpu_memory.used_bytes !== nothing
    @test cpu_memory.total_bytes >= cpu_memory.free_bytes
    @test memory_pressure(cpu_memory) !== nothing
    @test 0.0 <= memory_pressure(cpu_memory) <= 1.0
    @test backend_memory_status(ps).backend == :cpu
    @test !reclaim_backend_memory!(Val(:cpu))
    @test_throws ArgumentError BackendMemoryStatus(:cpu, true, false; total_bytes = -1)

    host = copy_particles_to_backend(Val(:cpu), ps)
    @test host isa ParticleSet{2,Float64}
    @test particle_storage_backend(host) == :cpu
    _assert_same_particles(host, ps)
    @test host !== ps
    @test host.x[1] !== ps.x[1]
    host.x[1][1] = -1.0
    @test ps.x[1][1] == 0.1

    roundtrip = copy_particles_to_host(host)
    @test particle_storage_backend(roundtrip) == :cpu
    _assert_same_particles(roundtrip, host)

    @test_throws DimensionMismatch ParticleSet{1,Float64}(
        ([1.0, 2.0],),
        ([0.0], [0.0], [0.0]),
        [1.0],
        UInt64[1],
        UInt32[0],
        1.0,
        1.0,
    )
    @test_throws ErrorException prepare_gpu_backend!(Val(:cuda))
    @test_throws ErrorException copy_particles_to_backend(Val(:cuda), ps)
    @test_throws ErrorException copy_particles_to_backend(Val(:metal), ps)
    @test_throws ErrorException backend_memory_status(Val(:cuda))
    @test_throws ErrorException reclaim_backend_memory!(Val(:metal))
end

@testset "optional CUDA particle storage" begin
    if Base.find_package("CUDA") !== nothing
        @eval import CUDA
        cuda_memory = backend_memory_status(Val(:cuda))
        @test cuda_memory.backend == :cuda
        @test cuda_memory.pool_supported ==
              (isdefined(CUDA, :pool_status) || isdefined(CUDA, :reclaim))
        pressure = memory_pressure(cuda_memory)
        @test pressure === nothing || 0.0 <= pressure <= 1.0
        @test reclaim_backend_memory!(Val(:cuda)) isa Bool
        if CUDA.functional()
            ps = _fill_backend_test_particles!()
            dev = copy_particles_to_backend(Val(:cuda), ps)
            @test extension_loaded(Val(:cuda))
            @test particle_storage_backend(dev) == :cuda
            @test backend_memory_status(dev).backend == :cuda
            @test all(A isa CUDA.CuArray for A in (dev.x..., dev.v..., dev.weight, dev.id, dev.tag))
            _assert_same_particles(copy_particles_to_host(dev), ps)
        else
            @test_skip "CUDA package is available but no functional CUDA device was found"
        end
    end
end

@testset "optional Metal particle storage" begin
    if Sys.isapple() && Base.find_package("Metal") !== nothing
        @eval import Metal
        metal_memory = backend_memory_status(Val(:metal))
        @test metal_memory.backend == :metal
        @test metal_memory.pool_supported == false
        @test metal_memory.note == "Metal.jl does not expose memory-pool counters"
        @test memory_pressure(metal_memory) === nothing
        @test !reclaim_backend_memory!(Val(:metal))
        metal_functional = !isdefined(Metal, :functional) || getproperty(Metal, :functional)()
        if metal_functional
            ps = _fill_backend_test_particles!(Float32)
            dev = copy_particles_to_backend(Val(:metal), ps)
            @test extension_loaded(Val(:metal))
            @test particle_storage_backend(dev) == :metal
            @test backend_memory_status(dev).backend == :metal
            @test all(
                A isa Metal.MtlArray for A in (dev.x..., dev.v..., dev.weight, dev.id, dev.tag)
            )
            _assert_same_particles(copy_particles_to_host(dev), ps)
            cpu = _fill_backend_test_particles!(Float32)
            gpu = copy_particles_to_backend(Val(:metal), _fill_backend_test_particles!(Float32))
            E = (0.25f0, -0.125f0, 0.0625f0)
            B = (0.05f0, -0.15f0, 0.3f0)
            dt = 0.02f0
            for _ = 1:25
                push_uniform!(cpu, E, B, dt)
                push_uniform!(gpu, E, B, dt)
            end
            _assert_particles_approx(copy_particles_to_host(gpu), cpu; rtol = 2.0f-5, atol = 2.0f-6)

            cpu = _fill_backend_test_particles!(Float32)
            gpu = copy_particles_to_backend(Val(:metal), _fill_backend_test_particles!(Float32))
            Eg = (
                Float32[0.22, -0.18, 0.11, -0.07],
                Float32[-0.08, 0.16, -0.12, 0.05],
                Float32[0.03, -0.04, 0.06, -0.02],
            )
            Bg = (
                Float32[0.04, -0.02, 0.03, -0.01],
                Float32[-0.06, 0.05, -0.03, 0.02],
                Float32[0.31, 0.28, 0.35, 0.26],
            )
            Eg_dev = ntuple(c -> Metal.MtlArray(Eg[c]), 3)
            Bg_dev = ntuple(c -> Metal.MtlArray(Bg[c]), 3)
            xmid_cpu = ntuple(_ -> zeros(Float32, nparticles(cpu)), 2)
            xmid_gpu = ntuple(_ -> Metal.MtlArray(zeros(Float32, nparticles(gpu))), 2)
            for _ = 1:17
                push_gathered!(cpu, Eg, Bg, dt; xmid = xmid_cpu)
                push_gathered!(gpu, Eg_dev, Bg_dev, dt; xmid = xmid_gpu)
            end
            _assert_particles_approx(copy_particles_to_host(gpu), cpu; rtol = 2.0f-5, atol = 2.0f-6)
            for d = 1:2
                @test collect(xmid_gpu[d]) ≈ xmid_cpu[d] rtol = 2.0f-5 atol = 2.0f-6
            end

            for D = 1:3, shape in (NGP(), CIC(), TSC())
                n = ntuple(d -> 6 + d, D)
                L = ntuple(d -> 1.0f0 + 0.25f0 * Float32(d), D)
                g = FourierGrid(n, L)
                Np = 7
                ps_cpu = ParticleSet{D,Float32}(Np)
                for d = 1:D, p = 1:Np
                    ps_cpu.x[d][p] = mod(0.09f0 * Float32(p) + 0.13f0 * Float32(d), L[d])
                end
                field_cpu = Array{Float32,D}(undef, n)
                for I in CartesianIndices(field_cpu)
                    idx = Tuple(I)
                    field_cpu[I] =
                        0.25f0 + sum((0.11f0 * Float32(d)) * Float32(idx[d] - 1) for d = 1:D)
                end
                out_cpu = zeros(Float32, Np)
                gather_scalar!(out_cpu, field_cpu, ps_cpu, g, shape)

                ps_gpu = copy_particles_to_backend(Val(:metal), ps_cpu)
                field_gpu = Metal.MtlArray(field_cpu)
                out_gpu = Metal.MtlArray(zeros(Float32, Np))
                gather_scalar!(out_gpu, field_gpu, ps_gpu, g, shape)
                Metal.synchronize()
                @test collect(out_gpu) ≈ out_cpu rtol = 2.0f-5 atol = 2.0f-6

                vals_cpu = Float32[0.2, -0.3, 0.4, 0.7, -0.1, 0.5, -0.6]
                dep_cpu = Array{Float32,D}(undef, n)
                deposit_scalar!(dep_cpu, ps_cpu, vals_cpu, g, shape)
                vals_gpu = Metal.MtlArray(vals_cpu)
                dep_gpu = Metal.MtlArray(zeros(Float32, n))
                deposit_scalar!(dep_gpu, ps_gpu, vals_gpu, g, shape)
                Metal.synchronize()
                @test Array(dep_gpu) ≈ dep_cpu rtol = 2.0f-5 atol = 2.0f-6

                if D == 2 && shape isa CIC
                    fields_cpu = ntuple(c -> field_cpu .+ Float32(c) / 10, 3)
                    outs_cpu = ntuple(_ -> zeros(Float32, Np), 3)
                    gather_vector!(outs_cpu, fields_cpu, ps_cpu, g, shape)
                    fields_gpu = ntuple(c -> Metal.MtlArray(fields_cpu[c]), 3)
                    outs_gpu = ntuple(_ -> Metal.MtlArray(zeros(Float32, Np)), 3)
                    gather_vector!(outs_gpu, fields_gpu, ps_gpu, g, shape)
                    Metal.synchronize()
                    for c = 1:3
                        @test collect(outs_gpu[c]) ≈ outs_cpu[c] rtol = 2.0f-5 atol = 2.0f-6
                    end
                end
            end

            if isdefined(Metal, :allowscalar)
                @test_throws Exception dev.x[1][1]
            end
            @test_throws ArgumentError copy_particles_to_backend(
                Val(:metal),
                _fill_backend_test_particles!(Float64),
            )
        else
            @test_skip "Metal package is available but no functional Metal device was found"
        end
    end
end
