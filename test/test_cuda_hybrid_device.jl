# Actual NVIDIA acceptance only. Missing hardware is a failure, never a skip.
include("test_cuda_hybrid_host.jl")
module CUDAHybridDeviceTests
using Test
using HybridPlasmaPIC
import CUDA

dev(a::Vector{T}) where {T} = CUDA.CuArray{T,1,CUDA.DeviceMemory}(a)
device_particle(ps) = ParticleSet{1,Float64}(
    map(dev, ps.x),
    map(dev, ps.v),
    dev(ps.weight),
    dev(ps.id),
    dev(ps.tag),
    ps.q,
    ps.m,
)
field_leaves(f) = (f.n, f.ui..., f.B..., f.E..., f.J..., f.pe, f.gradp..., f.ninv)
public_roles(st, species) = (
    st.fields,
    field_leaves(st.fields)...,
    species...,
    (a for ps in species for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))...,
)

function loaded_species(n, copies)
    species = ParticleSet{1,Float64}[]
    for (s, (q, m, density)) in enumerate(((1.0, 1.0, 1.0), (2.0, 4.0, 0.5), (3.0, 5.0, 0.0)))
        count = iszero(density) ? 0 : copies * n
        ps = ParticleSet{1,Float64}(count; q, m)
        for p = 1:count
            x = (p - 0.5) / copies
            ps.x[1][p] = x
            ps.v[1][p] = sinpi(2x / n) / 1024
            # Keep all nonempty rows inside one velocity-scale binade.
            ps.v[2][p] = (isone(s) ? 3.0 : -3.0) / 64 + sinpi(2x / n) / 1024
            ps.v[3][p] = s * cospi(2x / n) / 2048
        end
        fill!(ps.weight, density / copies)
        push!(species, ps)
    end
    return species
end

function physical_close(actual, fields, particles)
    # Declared small-fixture acceptance envelope; atomic/cuFFT bits may differ.
    for (a, b) in zip(field_leaves(actual.fields), field_leaves(fields))
        @test isapprox(a, b; atol = 2e-10, rtol = 0.0)
    end
    @test length(actual.species) == length(particles)
    for (a, b) in zip(actual.species, particles)
        @test isequal((a.q, a.m), (b.q, b.m))
        for (x, y) in zip((a.x..., a.v...), (b.x..., b.v...))
            @test isapprox(x, y; atol = 2e-10, rtol = 0.0)
        end
        @test isequal((a.weight, a.id, a.tag), (b.weight, b.id, b.tag))
    end
end

@testset "Actual resident CUDA lifecycle" begin
    @test CUDA.functional()
    @test CUDA.cuFFT.functional()
    CUDA.allowscalar(false)
    ext = Base.get_extension(HybridPlasmaPIC, :HybridPlasmaPICCUDAExt)
    @testset "empty allocation role identity" begin
        a, b = dev(Float64[]), dev(Float64[])
        alias = reinterpret(UInt64, a)
        @test !ext._ch_alias(a, b)
        @test ext._ch_alias(a, alias)
        @test ext._ch_independent(Any[a, b]) === nothing
        @test_throws ArgumentError ext._ch_independent(Any[a, alias])
    end
    for n in (8, 9)
        @testset "CPU discrete map and continuation N=$n" begin
            h, nb = 1 / 128, 2
            model = HybridModel(IsothermalElectrons(1 / 16); nfloor = 1e-8)
            host = loaded_species(n, 4)
            particles = map(device_particle, host)
            x = Float64.(0:(n-1))
            seed = (fill(0.5, n), cospi.(2 .* x ./ n) ./ 16, sinpi.(4 .* x ./ n) ./ 32)
            cpu = HybridStepper(FourierGrid((n,), (Float64(n),)), model, CIC(), 0)
            foreach(copyto!, cpu.fields.B, seed)
            st = cuda_hybrid_stepper(
                (n,),
                (Float64(n),),
                model,
                particles;
                dt = h,
                NB = nb,
                B = map(dev, seed),
            )
            init!(cpu, host)
            init!(st)
            initial = cuda_hybrid_host_snapshot(st)
            physical_close(initial, cpu.fields, host)
            @test initial.velocity_phase === :loaded_integer && !initial.primed
            for k = 1:2
                step!(cpu, host, h; NB = nb)
                step!(st)
                actual = cuda_hybrid_host_snapshot(st)
                physical_close(actual, cpu.fields, host)
                @test actual.time == cpu.time[] && actual.step == cpu.step[]
                @test actual.primed && actual.velocity_offset == -h / 2
                @test actual.velocity_phase === :integer_fields_half_velocity
                @test isequal(actual.fields.B[1], seed[1])
            end
            saved = cuda_hybrid_host_snapshot(st)
            @test maximum(abs, saved.fields.B[2] - seed[2]) > 1e-8
            roles = public_roles(st, particles)
            @test step!(st; dt = 0.0) === st
            @test step!(st; dt = -0.0) === st
            @test_throws ArgumentError step!(st; dt = -h)
            @test_throws ArgumentError step!(st; dt = h / 2)
            @test_throws ArgumentError step!(st; dt = nextfloat(0.0))
            @test isequal(cuda_hybrid_host_snapshot(st), saved)
            @test all(a === b for (a, b) in zip(public_roles(st, particles), roles))
            diagnostics = cuda_hybrid_host_diagnostics(st)
            @test isfinite(diagnostics.energy.total)
            @test diagnostics.time == saved.time && diagnostics.velocity_offset == -h / 2

            mktempdir() do directory
                raw, archive = joinpath(directory, "cuda.chk"), joinpath(directory, "cuda.run")
                save_checkpoint(raw, st)
                archive_run(archive, st; rng_seed = 0)
                step!(st)
                continuation = cuda_hybrid_host_snapshot(st)
                for (restore!, path) in ((load_checkpoint!, raw), (load_archive!, archive))
                    changed = map(device_particle, loaded_species(n, 2))
                    target = cuda_hybrid_stepper(
                        (n,),
                        (Float64(n),),
                        model,
                        changed;
                        dt = h,
                        NB = nb,
                        B = map(dev, seed),
                    )
                    init!(target)
                    step!(target)
                    @test length(changed[1].weight) != length(saved.species[1].weight)
                    restore!(target, path)
                    @test isequal(cuda_hybrid_host_snapshot(target), saved)
                    step!(target)
                    resumed = cuda_hybrid_host_snapshot(target)
                    physical_close(resumed, continuation.fields, continuation.species)
                    @test (resumed.time, resumed.step, resumed.primed) ==
                          (continuation.time, continuation.step, continuation.primed)
                end
            end
        end
    end

    @testset "real private push then density rejection preserves public state" begin
        h = 1 / 8
        host = ParticleSet{1,Float64}(4)
        host.x[1] .= Float64.(0:3)
        host.v[1] .= -host.x[1] ./ h
        particles = [device_particle(host)]
        seed = (fill(-0.0, 4), zeros(4), zeros(4))
        st = cuda_hybrid_stepper(
            (4,),
            (4.0,),
            HybridModel(IsothermalElectrons(0.0)),
            particles;
            dt = h,
            B = map(dev, seed),
        )
        init!(st)
        before, roles = cuda_hybrid_host_snapshot(st), public_roles(st, particles)
        @test all(x -> isequal(x, -0.0), before.fields.B[1])
        @test_throws r"CUDA hybrid density/ion velocity rejected category 5" step!(st)
        @test all(iszero, Array(st.motion[1].x[1]))
        @test iszero(Array(st.work.frozen_n)[4])
        @test isequal(cuda_hybrid_host_snapshot(st), before)
        @test all(a === b for (a, b) in zip(public_roles(st, particles), roles))
    end
end
end
