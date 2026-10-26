include("local_coulomb/angle_range.jl")
include("local_coulomb/component_range.jl")
include("local_coulomb/conditioned_rates.jl")
include("local_coulomb/widened_updates.jl")
include("local_coulomb/compensated_updates.jl")
include("local_coulomb/completion_cancellation.jl")

module LocalCoulombNormalizationRangeTests
using HybridPlasmaPIC, Test
@testset "Local collision normalization in an isolated process" begin
    project = Base.active_project()
    @test project !== nothing
    fixture = joinpath(@__DIR__, "local_coulomb", "normalization_fixture.jl")
    package = dirname(@__DIR__)
    command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(project)) $fixture $package`
    log = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = log, stderr = log))
    transcript = String(take!(log))
    print(transcript)
    success(process) || @error "Local collision normalization fixture failed" transcript
    @test success(process)
    @test occursin("Local Coulomb normalization settings restored", transcript)
end
end
