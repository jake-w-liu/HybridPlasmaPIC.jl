#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_17_distributed_fft_roundtrip(artifact_dir::AbstractString)
    id = "17_distributed_fft_roundtrip"
    try
        @eval import MPI
        @eval import PencilArrays
        @eval import PencilFFTs
    catch err
        artifact = joinpath(artifact_dir, "17_distributed_fft_roundtrip.csv")
        rows = (("distributed_fft_dependencies_available_error", 1.0, 0.0, "absolute", 1.0, 0.0),)
        _write_metric_csv(artifact, rows)
        return _metric_rows_to_results(
            id = id,
            category = "parallel_fft",
            reference_kind = "external_library",
            reference = "PencilFFTs/PencilArrays extension compared with FFTW",
            rows = rows,
            artifact = artifact,
            notes = "PencilFFTs/PencilArrays dependencies are not loadable: $(typeof(err))",
        )
    end

    return Base.invokelatest(_case_17_distributed_fft_roundtrip_loaded, artifact_dir)
end

function _case_17_distributed_fft_roundtrip_loaded(artifact_dir::AbstractString)
    id = "17_distributed_fft_roundtrip"
    plan = distributed_fft_plan((8, 6, 4); comm = MPI.COMM_SELF)
    input = distributed_fft_input(plan)
    output = distributed_fft_output(plan)
    local_input = parent(input)
    for index in CartesianIndices(local_input)
        i, j, k = Tuple(index)
        local_input[index] = complex(sin(0.2i) + cos(0.3j), 0.1k - 0.2j)
    end
    reference = copy(Array(input))
    distributed_fft_forward!(output, plan, input)
    fft_error = maximum(abs, PencilArrays.gather(output) .- fft(reference))
    distributed_fft_inverse!(input, plan, output)
    roundtrip_error = maximum(abs, Array(input) .- reference)

    artifact = joinpath(artifact_dir, "17_distributed_fft_roundtrip.csv")
    rows = (
        ("forward_fft_max_abs_error", fft_error, 0.0, "absolute", fft_error, 1e-11),
        ("roundtrip_max_abs_error", roundtrip_error, 0.0, "absolute", roundtrip_error, 1e-11),
    )
    _write_metric_csv(artifact, rows)
    return _metric_rows_to_results(
        id = id,
        category = "parallel_fft",
        reference_kind = "external_library",
        reference = "PencilFFTs/PencilArrays extension compared with FFTW on MPI.COMM_SELF and inverse roundtrip",
        rows = rows,
        artifact = artifact,
    )
end


VALIDATION_CASE = ValidationCase(
    id = "17_distributed_fft_roundtrip",
    default = true,
    description = "Distributed FFT extension against FFTW and inverse roundtrip.",
    runner = case_17_distributed_fft_roundtrip,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(
        _run_single_case_main(
            VALIDATION_CASE,
            ARGS;
            default_artifact_dir = joinpath(@__DIR__, "artifacts"),
        ),
    )
end
