#!/usr/bin/env julia

if !isdefined(@__MODULE__, :ValidationCase)
    include(joinpath(@__DIR__, "..", "common.jl"))
end

function case_09_normalization_migration_io(artifact_dir::AbstractString)
    id = "09_normalization_migration_io"
    units = PlasmaUnits(n0 = 1e6, B0 = 5e-9, mi = 1.6726e-27)
    kinds = (:length, :time, :velocity, :magnetic, :electric, :density, :current, :pressure)
    values = [-3.7, 0.0, 1.0, 2.5e3]
    unit_roundtrip = maximum(
        abs(to_normalized(to_SI(x, kind, units), kind, units) - x) / (abs(x) + 1.0) for
        kind in kinds for x in values
    )
    unit_identity_error = abs(inertial_length(units) - alfven_speed(units) / gyrofrequency(units))

    g = FourierGrid((8,), (8.0,))
    layout = LogicalRankLayout((4,); periodic = (true,))
    rank_error = abs(rank_of_position((8.2,), g, layout) - 1)
    bounds = rank_bounds(g, layout, 3)
    bounds_error = max(abs(bounds.lo[1] - 4.0), abs(bounds.hi[1] - 6.0))
    rank_index_error = nranks(layout) == 4 && rank_index(layout, (3,)) == 3 ? 0.0 : 1.0

    ranks = [ParticleSet{1,Float64}(1) for _ = 1:4]
    xs = (2.6, 5.1, 7.8, 10.2)
    ids = UInt64[11, 21, 31, 41]
    for r = 1:4
        ranks[r].x[1][1] = xs[r]
        ranks[r].id[1] = ids[r]
    end
    migration = migrate_particles!(ranks, g, layout)
    all_ids = sort(vcat((collect(rank.id) for rank in ranks)...))
    migration_id_error = all_ids == sort(ids) ? 0.0 : 1.0
    migration_count_error = abs(migration.moved - 4) + abs(migration.lost)
    append_dest = ParticleSet{1,Float64}(1)
    append_src = ParticleSet{1,Float64}(2)
    append_dest.id .= UInt64[7]
    append_src.id .= UInt64[8, 9]
    append_particles!(append_dest, append_src)
    append_error = nparticles(append_dest) == 3 && append_dest.id == UInt64[7, 8, 9] ? 0.0 : 1.0

    dt = 0.02
    ps_ref, st_ref = _checkpoint_validation_run(42)
    for _ = 1:8
        step!(st_ref, ps_ref, dt; NB = 2)
    end
    ps_a, st_a = _checkpoint_validation_run(42)
    for _ = 1:4
        step!(st_a, ps_a, dt; NB = 2)
    end
    path = tempname()
    try
        save_checkpoint(path, st_a, ps_a)
        ps_b, st_b = _checkpoint_validation_run(1)
        load_checkpoint!(st_b, ps_b, path)
        for _ = 1:4
            step!(st_b, ps_b, dt; NB = 2)
        end
        particle_error = maximum(abs, ps_b.x[1] .- ps_ref.x[1])
        velocity_error = maximum(maximum(abs, ps_b.v[c] .- ps_ref.v[c]) for c = 1:3)
        field_error = maximum(maximum(abs, st_b.fields.B[c] .- st_ref.fields.B[c]) for c = 1:3)
        step_error = abs(st_b.step[] - st_ref.step[])
        checkpoint_error = max(particle_error, velocity_error, field_error, step_error)
        artifact = joinpath(artifact_dir, "09_normalization_migration_io.csv")
        rows = (
            (
                "unit_roundtrip_max_relative_error",
                unit_roundtrip,
                0.0,
                "relative",
                unit_roundtrip,
                1e-12,
            ),
            (
                "normalization_scalar_identity_error",
                unit_identity_error,
                0.0,
                "absolute",
                unit_identity_error,
                1e-12,
            ),
            ("rank_periodic_wrap_abs_error", rank_error, 0.0, "absolute", rank_error, 0.0),
            (
                "rank_index_nranks_contract_error",
                rank_index_error,
                0.0,
                "absolute",
                rank_index_error,
                0.0,
            ),
            ("rank_bounds_max_abs_error", bounds_error, 0.0, "absolute", bounds_error, 0.0),
            (
                "migration_id_preservation_error",
                migration_id_error,
                0.0,
                "absolute",
                migration_id_error,
                0.0,
            ),
            (
                "migration_count_error",
                migration_count_error,
                0.0,
                "absolute",
                migration_count_error,
                0.0,
            ),
            (
                "append_particles_id_contract_error",
                append_error,
                0.0,
                "absolute",
                append_error,
                0.0,
            ),
            (
                "checkpoint_restart_max_abs_error",
                checkpoint_error,
                0.0,
                "absolute",
                checkpoint_error,
                0.0,
            ),
        )
        _write_metric_csv(artifact, rows)
        return _metric_rows_to_results(
            id = id,
            category = "normalization_parallel_io",
            reference_kind = "analytic",
            reference = "unit inverse maps, logical-rank geometry, migration invariants, bitwise checkpoint restart",
            rows = rows,
            artifact = artifact,
        )
    finally
        rm(path; force = true)
    end
end


VALIDATION_CASE = ValidationCase(
    id = "09_normalization_migration_io",
    default = true,
    description = "Unit normalization, logical migration, and checkpoint restart invariants.",
    runner = case_09_normalization_migration_io,
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
