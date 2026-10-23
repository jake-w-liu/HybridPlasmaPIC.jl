module NormalizationRangeRegressionTests
using HybridPlasmaPIC, Test
include("normalization/range.jl")
@testset "Normalization MPFR range in an isolated process" begin
    project = Base.active_project()
    @test project !== nothing
    fixture = joinpath(@__DIR__, "normalization", "mpfr_fixture.jl")
    package = dirname(@__DIR__)
    command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(project)) $fixture $package`
    log = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = log, stderr = log))
    transcript = String(take!(log))
    print(transcript)
    success(process) || @error "Normalization range fixture failed" transcript
    @test success(process)
    @test occursin("Normalization numerical settings restored", transcript)
end
end
include("normalization/metal.jl")
