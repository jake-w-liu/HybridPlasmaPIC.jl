using HybridPlasmaPIC, Test
using Dates
using Serialization

_sha256_hex(s) = s isa AbstractString && occursin(r"^[0-9a-f]{64}$", s)

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
    @test endswith(meta.timestamp, "Z")
    @test meta.rank_layout == "serial;ranks=1;dims=(1);coords=(0);periodic=false"
    @test occursin("ranks=1", meta.rank_layout)
    # Hashes: either "absent" or a stable SHA-256 hex digest. In this repo
    # Project.toml exists, so project_hash must be a digest, not "absent".
    @test meta.project_hash != "absent"
    @test _sha256_hex(meta.project_hash)
    @test meta.manifest_hash == "absent" || _sha256_hex(meta.manifest_hash)

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

    @test_throws ArgumentError capture_metadata(; rng_seed = big(typemax(Int)) + 1)
    @test_throws ArgumentError capture_metadata(; rng_seed = big(typemin(Int)) - 1)
end

@testset "capture_metadata auto timestamp is UTC" begin
    before = Dates.unix2datetime(time())
    meta = capture_metadata(; rng_seed = 3)
    after = Dates.unix2datetime(time())
    parsed = DateTime(chop(meta.timestamp; tail = 1), dateformat"yyyy-mm-ddTHH:MM:SS")
    @test endswith(meta.timestamp, "Z")
    @test before - Second(1) <= parsed <= after + Second(1)
end

@testset "schema constant" begin
    @test CHECKPOINT_SCHEMA_VERSION == 3
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
        container = deserialize(path)
        @test container.checksum isa String
        @test _sha256_hex(container.checksum)

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
    end
end

@testset "save_run / load_run checksum survives process boundary" begin
    project_dir = pkgdir(HybridPlasmaPIC)
    external_state = (a = [1.0, 2.0, 3.0], label = "external", step = 11)

    mktemp() do path, io
        close(io)
        writer = """
        using HybridPlasmaPIC
        meta = capture_metadata(; rng_seed = 17, timestamp = "2026-06-28T00:00:00Z")
        state = (a = [1.0, 2.0, 3.0], label = "external", step = 11)
        save_run(ARGS[1], state, meta)
        """
        run(`$(Base.julia_cmd()) --project=$project_dir -e $writer $path`)

        loaded = load_run(path)
        @test loaded.state == external_state
        @test loaded.meta.rng_seed == 17

        current_state = (a = [4.0, 5.0], label = "current", step = 12)
        save_run(path, current_state, capture_metadata(; rng_seed = 18))
        reader = """
        using HybridPlasmaPIC
        loaded = load_run(ARGS[1])
        loaded.state == (a = [4.0, 5.0], label = "current", step = 12) ||
            error("external load_run saw the wrong state")
        loaded.meta.rng_seed == 18 || error("external load_run saw the wrong metadata")
        """
        run(`$(Base.julia_cmd()) --project=$project_dir -e $reader $path`)
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
        container = deserialize(path)
        bad = (
            schema = CHECKPOINT_SCHEMA_VERSION,
            meta = meta,
            state = state,
            checksum = container.checksum == repeat("0", 64) ? repeat("1", 64) : repeat("0", 64),
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
            checksum = repeat("0", 64),
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
            checksum = repeat("0", 64),
        )
        serialize(path, wrong)
        @test_throws ErrorException load_run(path)
    end
end
