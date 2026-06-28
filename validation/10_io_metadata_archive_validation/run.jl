#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_10_io_metadata_archive_validation(artifact_dir::AbstractString)
    id = "10_io_metadata_archive_validation"
    paths = String[]
    try
        field_path = tempname()
        push!(paths, field_path)
        A = reshape([sin(0.2i) + cos(0.3j) for i = 1:4, j = 1:5], 4, 5)
        write_field(field_path, A)
        B = read_field(field_path)
        field_error = maximum(abs, B .- A)

        async_path = tempname()
        push!(paths, async_path)
        async_state = (a = fill(1.0, 10_000), nested = ([3.0, 4.0],), step = 9)
        task = async_save(async_path, async_state)
        fill!(async_state.a, 2.0)
        async_state.nested[1][1] = -7.0
        wait(task)
        async_loaded = deserialize(async_path)
        async_error =
            all(==(1.0), async_loaded.a) && async_loaded.nested[1] == [3.0, 4.0] &&
            async_loaded.step == 9 ? 0.0 : 1.0

        run_path = tempname()
        push!(paths, run_path)
        meta = capture_metadata(; rng_seed = 123, timestamp = "2026-06-28T00:00:00Z")
        state = (a = [1.0, 2.0, 3.0], label = "validation", step = 7)
        save_run(run_path, state, meta)
        loaded = load_run(run_path)
        metadata_type_error = meta isa RunMetadata && loaded.meta isa RunMetadata ? 0.0 : 1.0
        save_run_error =
            loaded.schema == CHECKPOINT_SCHEMA_VERSION && loaded.meta.rng_seed == 123 &&
            loaded.state == state ? 0.0 : 1.0

        archive_path = tempname()
        push!(paths, archive_path)
        ps, st = _checkpoint_validation_run(55)
        archive_run(archive_path, st, ps; rng_seed = 55, diagnostic_desc = "validation")
        archive = load_archive(archive_path)
        archive_error = max(
            abs(archive.meta.rng_seed - 55),
            abs(archive.state.step - st.step[]),
            maximum(abs, archive.state.x[1] .- ps.x[1]),
            maximum(abs, archive.state.B[2] .- st.fields.B[2]),
        )

        sample = sample_particles(ps, 7)
        sample_expected = collect(1:7:nparticles(ps))
        sample_error =
            sample.index == sample_expected &&
            sample.x[1] == ps.x[1][sample_expected] &&
            sample.v[1] == ps.v[1][sample_expected] ? 0.0 : 1.0

        operator_error = operators_match(FourierGrid((8, 6), (2π, 2π))) ? 0.0 : 1.0

        artifact = joinpath(artifact_dir, "10_io_metadata_archive_validation.csv")
        rows = (
            ("write_read_field_max_abs_error", field_error, 0.0, "absolute", field_error, 0.0),
            ("async_save_snapshot_error", async_error, 0.0, "absolute", async_error, 0.0),
            ("run_metadata_type_contract_error", metadata_type_error, 0.0, "absolute", metadata_type_error, 0.0),
            ("save_run_load_run_roundtrip_error", save_run_error, 0.0, "absolute", save_run_error, 0.0),
            ("archive_run_load_archive_max_abs_error", archive_error, 0.0, "absolute", archive_error, 0.0),
            ("sample_particles_stride_error", sample_error, 0.0, "absolute", sample_error, 0.0),
            ("operators_match_boolean_error", operator_error, 0.0, "absolute", operator_error, 0.0),
        )
        _write_metric_csv(artifact, rows)
        return _metric_rows_to_results(
            id = id,
            category = "io_metadata_archive",
            reference_kind = "roundtrip",
            reference = "exact Serialization/raw-field roundtrips, archive checksum load, sampled particle indices, operator parity",
            rows = rows,
            artifact = artifact,
        )
    finally
        for path in paths
            rm(path; force = true)
        end
    end
end


VALIDATION_CASE = ValidationCase(
    id = "10_io_metadata_archive_validation",
    default = true,
    description = "Metadata, raw field IO, async save, archive, and sampled dump roundtrips.",
    runner = case_10_io_metadata_archive_validation,
)

if abspath(PROGRAM_FILE) == @__FILE__
    exit(_run_single_case_main(VALIDATION_CASE, ARGS; default_artifact_dir = joinpath(@__DIR__, "artifacts")))
end
