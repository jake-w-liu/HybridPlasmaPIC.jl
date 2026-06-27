using HybridPlasmaPIC, Test
using Serialization

@testset "write_field / read_field round-trip" begin
    # 1D, 2D, 3D Float64 arrays must round-trip *exactly* (bit-for-bit).
    a1 = randn(17)
    a2 = randn(5, 8)
    a3 = randn(3, 4, 5)
    for A in (a1, a2, a3)
        mktemp() do path, io
            close(io)
            ret = write_field(path, A)
            @test ret == path
            B = read_field(path)
            @test B isa Array{Float64,ndims(A)}
            @test size(B) == size(A)
            @test B == A                      # exact, not approx
        end
    end

    # Edge cases: empty array and 0-dim-ish small arrays.
    for A in (Float64[], zeros(0, 3), reshape([42.0], 1, 1, 1))
        mktemp() do path, io
            close(io)
            write_field(path, A)
            B = read_field(path)
            @test size(B) == size(A)
            @test B == A
        end
    end

    # Non-Float64 input is converted to Float64 and round-trips by value.
    Ai = reshape(collect(1:12), 3, 4)         # Int matrix
    mktemp() do path, io
        close(io)
        write_field(path, Ai)
        B = read_field(path)
        @test B == Float64.(Ai)
        @test eltype(B) === Float64
    end

    # A non-contiguous view must serialize by value, not by layout assumption.
    base = randn(6, 6)
    v = @view base[2:5, 1:3]
    mktemp() do path, io
        close(io)
        write_field(path, v)
        B = read_field(path)
        @test size(B) == size(v)
        @test B == Array(v)
    end

    # Special Float64 values (Inf/-Inf/NaN/-0.0) survive the byte round-trip.
    sp = [Inf, -Inf, 0.0, -0.0, 1e-300, 1e300]
    mktemp() do path, io
        close(io)
        write_field(path, sp)
        B = read_field(path)
        @test B[1] == Inf && B[2] == -Inf
        @test B[3] === 0.0 && signbit(B[4])    # -0.0 sign bit preserved
        @test B[5] == 1e-300 && B[6] == 1e300
    end
    # NaN separately (NaN != NaN).
    mktemp() do path, io
        close(io)
        write_field(path, [NaN, 1.0])
        B = read_field(path)
        @test isnan(B[1]) && B[2] == 1.0
    end
end

@testset "read_field rejects malformed files" begin
    mktemp() do path, io
        close(io)
        write(path, "not a field dump at all, definitely > 8 bytes")
        @test_throws ErrorException read_field(path)
    end
    # Truncated: write a valid header claiming more data than present.
    mktemp() do path, io
        # Valid magic + version + ndims=1 + size=100, but no data follows.
        write(io, b"HPSTDFLD")
        write(io, UInt8(1))
        write(io, Int64(1))
        write(io, Int64(100))
        close(io)
        @test_throws Exception read_field(path)
    end
end

@testset "async_save returns a Task and stages the write" begin
    state =
        (B = [randn(4, 4) for _ = 1:2], E = (zeros(3), ones(3)), t = 1.25, step = 9, tag = "async")
    mktemp() do path, io
        close(io)
        task = async_save(path, state)
        @test task isa Task
        wait(task)                            # block until the write finishes
        @test istaskdone(task)
        loaded = deserialize(path)
        @test loaded == state
        @test loaded.B == state.B
        @test loaded.E == state.E
        @test loaded.t == state.t
        @test loaded.step == state.step
        @test loaded.tag == state.tag
    end
end

@testset "capture_metadata backend/hardware fields" begin
    meta = capture_metadata(; rng_seed = 1)
    @test meta.backend == "CPU"
    @test !isempty(meta.backend)
    @test !isempty(meta.hardware)
    # hardware encodes the logical thread count via the documented " x<N>" tail.
    @test occursin(" x$(Sys.CPU_THREADS)", meta.hardware)
    @test meta.rank_layout == "serial;ranks=1;dims=(1);coords=(0);periodic=false"

    # New fields must survive the save_run/load_run checksum round-trip.
    state = (a = [1.0, 2.0], n = 3)
    mktemp() do path, io
        close(io)
        save_run(path, state, meta)
        loaded = load_run(path)
        @test loaded.meta.backend == meta.backend
        @test loaded.meta.hardware == meta.hardware
        for f in fieldnames(RunMetadata)
            @test getfield(loaded.meta, f) == getfield(meta, f)
        end
    end
end
