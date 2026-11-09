# Private CI launcher; scientific tests remain in test/.
using Pkg

length(ARGS) == 3 || error("Usage: ci_distributed_fft.jl MODE ENVIRONMENT RANKS")
mode, environment, rank_text = ARGS
mode in ("prepare", "core", "metadata", "models") ||
    error("Unsupported distributed FFT mode: $mode")
ranks = tryparse(Int, rank_text)
ranks in (2, 4, 8) || error("Unsupported MPI rank count: $rank_text")
environment = abspath(environment)
repository = realpath(joinpath(@__DIR__, ".."))

if mode == "prepare"
    Pkg.activate(environment)
    Pkg.develop(path = repository)
    Pkg.add(["MPI", "FFTW", "PencilArrays", "PencilFFTs"])
else
    isfile(joinpath(environment, "Project.toml")) || error("Missing prepared project")
    isfile(joinpath(environment, "Manifest.toml")) || error("Missing prepared manifest")
    Pkg.activate(environment)
    import MPI

    launcher = MPI.mpiexec()
    launcher === nothing && error("MPI launcher unavailable")
    filename = if mode == "core"
        "test_distributed_fft_core.jl"
    elseif mode == "metadata"
        "test_distributed_metadata_alias.jl"
    else
        "test_distributed_models_core.jl"
    end
    script = joinpath(repository, "test", filename)
    command = `$launcher -n $ranks $(Base.julia_cmd()) --project=$environment $script`
    println("Distributed FFT MPI command: ", command.exec)
    run(command)
end
