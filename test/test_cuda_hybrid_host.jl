# Run in a fresh process: the ordinary extension tests require CUDA unloaded.
module CUDAHybridHostTests
using Test
using HybridPlasmaPIC
import CUDA

@testset "Resident CUDA extension load and pure host admission" begin
    ext = Base.get_extension(HybridPlasmaPIC, :HybridPlasmaPICCUDAExt)
    @test ext !== nothing
    @test extension_dependency_module(Val(:cuda)) === CUDA
    versions = map(Base.pkgversion, (CUDA, CUDA.cuFFT, HybridPlasmaPIC.SpectralOperators))
    @test all(v -> v isa VersionNumber, versions)
    resident = ext._ResidentCUDAHybrid{Nothing,Nothing,Nothing}
    @test isconcretetype(resident) && resident <: CUDAHybridStepper
    for operation in (init!, step!, cuda_hybrid_host_snapshot, cuda_hybrid_host_diagnostics)
        @test hasmethod(operation, Tuple{resident})
    end
    @test hasmethod(save_checkpoint, Tuple{String,resident})
    @test hasmethod(load_checkpoint!, Tuple{resident,String})
    @test hasmethod(archive_run, Tuple{String,resident})
    @test hasmethod(load_archive!, Tuple{resident,String})

    # No context, device, stream, plan, CuArray or functional() call is made here.
    @test HybridPlasmaPIC._boris_accuracy_mode()
    @test ext._ch_supported(0.0) && ext._ch_supported(-0.0)
    @test ext._ch_supported(floatmin(Float64))
    @test !ext._ch_supported(prevfloat(floatmin(Float64)))
    @test !ext._ch_supported(Inf) && !ext._ch_supported(NaN)
    h = 1 / 128
    @test ext._ch_fixed_step(h, 2) == (h, 2)
    @test isequal(ext._ch_duration(-0.0, h), -0.0)
    @test_throws ArgumentError ext._ch_duration(-h, h)
    @test_throws ArgumentError ext._ch_duration(h / 2, h)
    @test_throws ArgumentError ext._ch_duration(nextfloat(0.0), h)
    @test_throws ArgumentError ext._ch_positive("q", -1.0)
    @test ext._ch_phase(true, false, 0.0, 0, 0.0, h) === nothing
    @test ext._ch_phase(true, true, h, 1, h, h) === nothing
    @test_throws ArgumentError ext._ch_phase(true, true, 0.0, 0, 0.0, h)
    @test_throws ArgumentError ext._ch_phase(true, false, h, 1, h, h)

    model = HybridModel(IsothermalElectrons(1 / 16); nfloor = 1e-8)
    factory_types = Tuple{Tuple{Int},Tuple{Float64},typeof(model),Vector{ParticleSet{1,Float64}}}
    @test which(cuda_hybrid_stepper, factory_types).module === ext
    @test ext._ch_model(model) === nothing
    @test_throws ArgumentError ext._ch_model(
        HybridModel(IsothermalElectrons(0.0); η = nextfloat(0.0)),
    )
    g = (n = (4,), L = (4.0,), dx = (1.0,), gain = 1.0)
    scalars = [(1.0, 1.0), (2.0, 4.0)]
    configuration = ext._ch_configuration(g, model, h, 2, scalars)
    @test (configuration.dt, configuration.NB, configuration.scalars) == (h, 2, Tuple(scalars))
    @test configuration.shape === :CIC && configuration.boundary === :periodic_constant_Bx

    # Exact original-operand witness: a physical-scale FMA residual underflows.
    q, m = floatmin(Float64), 3 / 4
    factor = HybridPlasmaPIC._BorisSpeciesFactors(q, m, q / m)
    Q = Rational{BigInt}
    @test Q(factor.qm) * Q(m) != Q(q)
    @test !HybridPlasmaPIC._rel_prime_exact_coefficient(factor)
    @test HybridPlasmaPIC._boris_relative_coefficient(factor, Float64) > 2.0^-53
    exact = HybridPlasmaPIC._BorisSpeciesFactors(m, m, 1.0)
    @test HybridPlasmaPIC._rel_prime_exact_coefficient(exact)

    # Independently derived maxima/minimum for the actual device rollback fixture.
    unit = HybridPlasmaPIC._BorisSpeciesFactors(1.0, 1.0, 1.0)
    @test HybridPlasmaPIC._euler_batch_certificate(1 / 16, unit, 1 / 16, 24.0, 0.0, 0.0, 8.0)
    admitted, vbound =
        HybridPlasmaPIC._boris_batch_admission(1 / 16, unit, 1 / 8, 24.0, 0.0, 0.0, 8.0)
    @test admitted && HybridPlasmaPIC._coordinate_batch_certified(3.0, 1 / 8, vbound)
end
end
