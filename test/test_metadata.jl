using HybridPlasmaPIC, Test
using Serialization

@testset "RunMetadata / capture_metadata" begin
    meta = capture_metadata(;
        rng_seed = 12345,
        normalization = "Omega_ci",
        filter_desc = "binomial(2)",
        boundary_desc = "periodic",
        diagnostic_desc = "energy+spectra",
    )
    @test meta isa RunMetadata
    # Populated, type-correct fields.
    @test meta.julia_version == string(VERSION)
    @test !isempty(meta.git_commit)            # "unknown" or a real SHA, never empty
    @test meta.rng_seed == 12345
    @test meta.normalization == "Omega_ci"
    @test meta.filter_desc == "binomial(2)"
    @test meta.boundary_desc == "periodic"
    @test meta.diagnostic_desc == "energy+spectra"
    @test !isempty(meta.timestamp)             # auto-filled
    @test meta.rank_layout == "serial;ranks=1;dims=(1);coords=(0);periodic=false"
    @test occursin("ranks=1", meta.rank_layout)
    # Hashes: either "absent" or a parseable hash string. In this repo Project.toml
    # exists, so project_hash must be a numeric-string hash, not "absent".
    @test meta.project_hash != "absent"
    @test tryparse(UInt64, meta.project_hash) !== nothing
    @test meta.manifest_hash == "absent" || tryparse(UInt64, meta.manifest_hash) !== nothing

    # Defaults branch.
    meta2 = capture_metadata(; rng_seed = -7)
    @test meta2.normalization == "Omega_ci"
    @test meta2.filter_desc == "" && meta2.boundary_desc == "" && meta2.diagnostic_desc == ""
    @test meta2.rng_seed == -7

    # Explicit timestamp is preserved verbatim.
    meta3 = capture_metadata(; rng_seed = 1, timestamp = "2026-06-25T00:00:00Z")
    @test meta3.timestamp == "2026-06-25T00:00:00Z"

    # Explicit distributed-memory rank layout is stored verbatim for provenance.
    layout = "mpi;ranks=4;dims=(2,2);coords=(1,0);periodic=(true,false)"
    meta4 = capture_metadata(; rng_seed = 2, rank_layout = layout)
    @test meta4.rank_layout == layout
end

@testset "schema constant" begin
    @test CHECKPOINT_SCHEMA_VERSION == 2
end

@testset "save_run / load_run round-trip" begin
    meta = capture_metadata(; rng_seed = 999, filter_desc = "none")
    state = (
        B = [randn(4, 4) for _ = 1:3],
        E = (zeros(4, 4), ones(4, 4), fill(2.0, 4, 4)),
        time = 3.5,
        step = 42,
        label = "shock_run",
    )

    mktemp() do path, io
        close(io)
        ret = save_run(path, state, meta)
        @test ret == path

        loaded = load_run(path)
        @test loaded.schema == CHECKPOINT_SCHEMA_VERSION
        # meta exact round-trip (compare every field).
        for f in fieldnames(RunMetadata)
            @test getfield(loaded.meta, f) == getfield(meta, f)
        end
        # state exact round-trip.
        @test loaded.state.B == state.B
        @test loaded.state.E == state.E
        @test loaded.state.time == state.time
        @test loaded.state.step == state.step
        @test loaded.state.label == state.label
        # hash equality is what the integrity check relies on.
        @test hash(loaded.state) == hash(state)
    end
end

@testset "corruption detection (bad checksum)" begin
    meta = capture_metadata(; rng_seed = 5)
    state = (a = [1.0, 2.0, 3.0], n = 7)

    mktemp() do path, io
        close(io)
        # Valid save then load succeeds.
        save_run(path, state, meta)
        @test load_run(path).state == state

        # Hand-craft a container with a deliberately wrong checksum.
        bad = (
            schema = CHECKPOINT_SCHEMA_VERSION,
            meta = meta,
            state = state,
            checksum = hash(state) + 0x1,
        )   # wrong on purpose
        serialize(path, bad)
        @test_throws ErrorException load_run(path)
    end
end

@testset "schema mismatch detection" begin
    meta = capture_metadata(; rng_seed = 5)
    state = (a = 1, b = 2)
    mktemp() do path, io
        close(io)
        wrong = (
            schema = CHECKPOINT_SCHEMA_VERSION + 1,
            meta = meta,
            state = state,
            checksum = hash(state),
        )   # checksum correct, schema wrong
        serialize(path, wrong)
        @test_throws ErrorException load_run(path)
    end
end

@testset "malformed container detection" begin
    mktemp() do path, io
        close(io)
        serialize(path, "not a container")
        @test_throws ErrorException load_run(path)
    end
end

@testset "malformed metadata detection" begin
    state = (a = 1, b = 2)
    mktemp() do path, io
        close(io)
        wrong = (
            schema = CHECKPOINT_SCHEMA_VERSION,
            meta = (not = "metadata",),
            state = state,
            checksum = hash(state),
        )
        serialize(path, wrong)
        @test_throws ErrorException load_run(path)
    end
end
