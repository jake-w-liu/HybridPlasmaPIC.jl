using HybridPlasmaPIC, Test, Serialization

struct PublicationFailure <: Exception
    stage::Symbol
end

struct PublicationValue
    values::Vector{Int}
    fail::Bool
    result::Symbol
end

function Serialization.serialize(s::Serialization.AbstractSerializer, value::PublicationValue)
    Serialization.serialize_any(s, value)
    value.fail && throw(PublicationFailure(:serialize))
    return value.result
end

struct PublicationElectrons <: HybridPlasmaPIC.ElectronClosure
    temperature::Float64
    fail::Bool
end

function Serialization.serialize(s::Serialization.AbstractSerializer, value::PublicationElectrons)
    Serialization.serialize_any(s, value)
    value.fail && throw(PublicationFailure(:serialize))
    return nothing
end

mutable struct PublicationDeepcopyFailure end
Base.deepcopy_internal(::PublicationDeepcopyFailure, ::IdDict) =
    throw(PublicationFailure(:deepcopy))

mutable struct PublicationFaultIO{F,G} <: IO
    inner::IOStream
    remaining::Int
    write_error::Union{Nothing,Exception}
    flush_error::Union{Nothing,Exception}
    close_error::Union{Nothing,Exception}
    before_close::F
    after_close::G
end

function PublicationFaultIO(
    io;
    remaining = typemax(Int),
    write_error = nothing,
    flush_error = nothing,
    close_error = nothing,
    before_close = () -> nothing,
    after_close = () -> nothing,
)
    return PublicationFaultIO(
        io,
        remaining,
        write_error,
        flush_error,
        close_error,
        before_close,
        after_close,
    )
end

Base.stat(io::PublicationFaultIO) = stat(io.inner)
Base.isopen(io::PublicationFaultIO) = isopen(io.inner)
Base.iswritable(io::PublicationFaultIO) = iswritable(io.inner)
Base.position(io::PublicationFaultIO) = position(io.inner)
Base.write(io::PublicationFaultIO, x::UInt8) = write(io, UInt8[x])
function Base.unsafe_write(io::PublicationFaultIO, pointer::Ptr{UInt8}, count::UInt)
    if io.write_error !== nothing && count > io.remaining
        written = unsafe_write(io.inner, pointer, UInt(io.remaining))
        io.remaining -= Int(written)
        throw(io.write_error)
    end
    written = unsafe_write(io.inner, pointer, count)
    io.remaining -= Int(written)
    return written
end
function Base.flush(io::PublicationFaultIO)
    flush(io.inner)
    io.flush_error === nothing || throw(io.flush_error)
    return nothing
end
function Base.close(io::PublicationFaultIO)
    io.before_close()
    close(io.inner)
    io.after_close()
    io.close_error === nothing || throw(io.close_error)
    return nothing
end

function publication_error(f)
    try
        f()
    catch err
        return err
    end
    return nothing
end

function publication_run()
    grid = FourierGrid((4,), (2π,))
    state = HybridStepper(grid, HybridModel(PublicationElectrons(0.1, false)), CIC(), 2)
    particles = ParticleSet{1,Float64}(2; q = 1.0, m = 1.0)
    particles.x[1] .= [0.1, 0.2]
    foreach(a -> fill!(a, 0.0), particles.v)
    fill!(particles.weight, 1.0)
    return state, particles
end

@testset "Serial publication retains previous checkpoints on serializer errors" begin
    state, particles = publication_run()
    mktempdir() do dir
        path = joinpath(dir, "checkpoint")
        @test save_checkpoint(path, state, particles) == path
        before, entries = read(path), readdir(dir)
        state.model = HybridModel(PublicationElectrons(0.1, true))
        err = publication_error(() -> save_checkpoint(path, state, particles))
        @test err isa PublicationFailure
        @test err.stage === :serialize
        @test read(path) == before
        @test readdir(dir) == entries
        state.model = HybridModel(PublicationElectrons(0.2, false))
        @test save_checkpoint(path, state, particles) == path
        @test deserialize(path).model.closure.temperature == 0.2
        @test readdir(dir) == entries
    end
end

@testset "Atomic writer preserves bytes, return values and owned resources" begin
    for stage in (:write, :flush, :close, :close_before)
        mktempdir() do dir
            target = joinpath(dir, "checkpoint")
            serialize(target, :previous)
            before = read(target)
            unrelated = joinpath(dir, "unrelated")
            write(unrelated, "retain")
            temporary, raw = mktemp(dir; cleanup = false)
            failure = PublicationFailure(stage)
            close_calls = Ref(0)
            before_close = () -> begin
                close_calls[] += 1
                stage === :close_before && close_calls[] == 1 && throw(failure)
            end
            io = PublicationFaultIO(
                raw;
                remaining = 2,
                write_error = stage === :write ? failure : nothing,
                flush_error = stage === :flush ? failure : nothing,
                close_error = stage === :close ? failure : nothing,
                before_close,
            )
            err = publication_error() do
                HybridPlasmaPIC._publish_serialized(target, temporary, io, 0o600) do stream
                    write(stream, UInt8[1, 2, 3, 4])
                end
            end
            @test err === failure
            @test close_calls[] == (stage === :close_before ? 2 : 1)
            @test !isopen(raw)
            @test !ispath(temporary)
            @test read(target) == before
            @test read(unrelated, String) == "retain"
            @test sort(readdir(dir)) == ["checkpoint", "unrelated"]
            @test HybridPlasmaPIC._atomic_serialize(
                target,
                PublicationValue([4, 5], false, :saved),
            ) === :saved
            @test deserialize(target).values == [4, 5]
        end
    end
    # A close that still fails on cleanup remains observable, and the test owns
    # the final raw-stream close. The production owner retries exactly once.
    mktempdir() do dir
        target = joinpath(dir, "checkpoint")
        serialize(target, :previous)
        before = read(target)
        temporary, raw = mktemp(dir; cleanup = false)
        failures = (PublicationFailure(:close), PublicationFailure(:cleanup_close))
        calls = Ref(0)
        io = PublicationFaultIO(raw; before_close = () -> begin
            calls[] += 1
            throw(failures[calls[]])
        end)
        try
            err = publication_error() do
                HybridPlasmaPIC._publish_serialized(target, temporary, io, 0o600) do stream
                    write(stream, "new")
                end
            end
            @test err isa CompositeException
            @test length(err.exceptions) == 2
            @test err.exceptions[1].ex === failures[1]
            @test err.exceptions[2].ex === failures[2]
            @test calls[] == 2
            @test isopen(raw)
            @test !ispath(temporary)
            @test read(target) == before
        finally
            close(raw)
        end
    end
end

@testset "Atomic publication handles missing staging and chmod errors" begin
    for failure in (:missing, :mode)
        mktempdir() do dir
            target = joinpath(dir, "checkpoint")
            serialize(target, :previous)
            before = read(target)
            temporary, raw = mktemp(dir; cleanup = false)
            io = PublicationFaultIO(
                raw;
                after_close = () ->
                    (failure === :missing && ispath(temporary) ? rm(temporary) : nothing),
            )
            written = Ref(false)
            # A mode outside Cint's domain fails in Base.chmod before its syscall.
            mode = failure === :mode ? big(1) << 100 : 0o600
            err = publication_error() do
                HybridPlasmaPIC._publish_serialized(target, temporary, io, mode) do stream
                    written[] = true
                    write(stream, "new")
                end
            end
            @test err isa (failure === :missing ? ArgumentError : InexactError)
            @test written[] == (failure === :missing)
            @test !isopen(raw)
            @test !ispath(temporary)
            @test read(target) == before
            @test readdir(dir) == ["checkpoint"]
        end
    end
end

@testset "Atomic publication preserves replacement staging entries" begin
    for replacement in (:symlink, :dangling_link, :file, :directory), fail_write in (false, true)
        mktempdir() do dir
            target = joinpath(dir, "checkpoint")
            serialize(target, :previous)
            before = read(target)
            external = joinpath(dir, "external")
            write(external, "retain")
            chmod(external, 0o640)
            temporary, raw = mktemp(dir; cleanup = false)
            failure = PublicationFailure(:write)
            displaced = joinpath(dir, "displaced-stage")
            after_close = () -> begin
                # Keep the old inode alive so this fixture is independent of
                # filesystem inode-reuse policies after unlinking a closed file.
                Base.Filesystem.rename(temporary, displaced)
                if replacement === :symlink
                    symlink(external, temporary)
                elseif replacement === :dangling_link
                    symlink(joinpath(dir, "absent"), temporary)
                elseif replacement === :file
                    write(temporary, "unrelated")
                else
                    mkdir(temporary)
                    write(joinpath(temporary, "unrelated"), "retain")
                end
            end
            io = PublicationFaultIO(
                raw;
                remaining = 2,
                write_error = fail_write ? failure : nothing,
                after_close,
            )
            err = publication_error() do
                HybridPlasmaPIC._publish_serialized(target, temporary, io, 0o600) do stream
                    write(stream, UInt8[1, 2, 3, 4])
                end
            end
            @test err isa CompositeException
            @test length(err.exceptions) == 2
            @test fail_write ? err.exceptions[1].ex === failure :
                  err.exceptions[1].ex isa ArgumentError
            @test err.exceptions[2].ex isa ArgumentError
            @test !isopen(raw)
            @test read(displaced) == (fail_write ? UInt8[1, 2] : UInt8[1, 2, 3, 4])
            @test read(target) == before
            @test read(external, String) == "retain"
            @test filemode(external) & 0o777 == 0o640
            if replacement in (:symlink, :dangling_link)
                @test islink(temporary)
                @test readlink(temporary) ==
                      (replacement === :symlink ? external : joinpath(dir, "absent"))
            elseif replacement === :file
                @test read(temporary, String) == "unrelated"
            else
                @test read(joinpath(temporary, "unrelated"), String) == "retain"
            end
        end
    end
    # Real rename and unlink failures: the directory remains searchable but is
    # no longer writable. Restore its mode before the fixture removes its files.
    if !Sys.iswindows() && ccall(:geteuid, Cuint, ()) != 0
        mktempdir() do dir
            target = joinpath(dir, "checkpoint")
            serialize(target, :previous)
            before = read(target)
            temporary, raw = mktemp(dir; cleanup = false)
            io = PublicationFaultIO(raw; after_close = () -> chmod(dir, 0o500))
            try
                err = publication_error() do
                    HybridPlasmaPIC._publish_serialized(target, temporary, io, 0o600) do stream
                        write(stream, "new")
                    end
                end
                @test err isa CompositeException
                @test length(err.exceptions) == 2
                @test err.exceptions[1].ex isa Base.IOError
                @test occursin("rename", sprint(showerror, err.exceptions[1].ex))
                @test err.exceptions[2].ex isa Base.IOError
                @test isfile(temporary)
                @test !isopen(raw)
                @test read(target) == before
            finally
                chmod(dir, 0o700)
            end
        end
    end
end

@testset "Atomic serial callers preserve encoding, snapshots and task errors" begin
    mktempdir() do dir
        target = joinpath(dir, "checkpoint")
        state, particles = publication_run()
        buffer = IOBuffer()
        serialize(buffer, HybridPlasmaPIC._checkpoint_state(state, particles))
        @test save_checkpoint(target, state, particles) == target
        @test read(target) == take!(buffer)
        meta = capture_metadata(; rng_seed = 7)
        value = (x = [1.0, 2.0], step = 3)
        payload = IOBuffer()
        serialize(payload, value)
        bytes = take!(payload)
        container = (
            schema = CHECKPOINT_SCHEMA_VERSION,
            meta = meta,
            state = bytes,
            checksum = bytes2hex(HybridPlasmaPIC.SHA.sha256(bytes)),
        )
        buffer = IOBuffer()
        serialize(buffer, container)
        @test save_run(target, value, meta) == target
        @test read(target) == take!(buffer)
        @test load_run(target).state == value
        before, entries = read(target), readdir(dir)
        @test_throws PublicationFailure save_run(target, PublicationValue([1], true, :unused), meta)
        @test read(target) == before
        @test readdir(dir) == entries
        @test_throws PublicationFailure async_save(target, PublicationDeepcopyFailure())
        @test read(target) == before
        task = async_save(target, PublicationValue([1], true, :unused))
        @test task isa Task
        @test_throws TaskFailedException wait(task)
        @test istaskfailed(task)
        @test read(target) == before
        @test readdir(dir) == entries
        value = PublicationValue([1, 2, 3], false, :saved)
        task = async_save(target, value)
        value.values .= 9
        @test fetch(task) === :saved
        @test deserialize(target).values == [1, 2, 3]
        @test readdir(dir) == entries
    end
end

@testset "Serial callers retain path, permission and validation contracts" begin
    for caller in (:checkpoint, :archive, :async)
        state, particles = publication_run()
        metadata = capture_metadata(; rng_seed = 12)
        write_value = if caller === :checkpoint
            path -> save_checkpoint(path, state, particles)
        elseif caller === :archive
            path -> save_run(path, (value = [1, 2],), metadata)
        else
            path -> fetch(async_save(path, (value = [1, 2],)))
        end
        mktempdir() do dir
            target = joinpath(dir, "target")
            write_value(target)
            @test filemode(target) & 0o777 == 0o600
            before = read(target)
            chmod(target, 0o751)
            write_value(target)
            @test filemode(target) & 0o777 == 0o751
            @test read(target) == before
            link = joinpath(dir, "link")
            symlink("target", link)
            write_value(link)
            @test islink(link)
            @test readlink(link) == "target"
            @test read(target) == before
            @test filemode(target) & 0o777 == 0o751
            created = joinpath(dir, "created")
            dangling = joinpath(dir, "dangling")
            symlink("created", dangling)
            write_value(dangling)
            @test islink(dangling)
            @test read(created) == before
            @test filemode(created) & 0o777 == 0o600
            hardlink = joinpath(dir, "hardlink")
            Base.Filesystem.hardlink(target, hardlink)
            write(target, "old hardlink bytes")
            write_value(target)
            @test read(target) == before
            @test read(hardlink, String) == "old hardlink bytes"
            @test !Base.samefile(target, hardlink)
            if !Sys.iswindows() && ccall(:geteuid, Cuint, ()) != 0
                chmod(target, 0o440)
                entries = readdir(dir)
                err = publication_error(() -> write_value(target))
                @test err isa (caller === :async ? TaskFailedException : SystemError)
                @test read(target) == before
                @test filemode(target) & 0o777 == 0o440
                @test readdir(dir) == entries
            end
            entries = readdir(dir)
            @test publication_error(() -> write_value(joinpath(dir, "missing", "target"))) !==
                  nothing
            @test publication_error(() -> write_value(dir)) !== nothing
            @test readdir(dir) == entries
        end
        mktempdir() do dir
            mkpath(joinpath(dir, "left"))
            mkpath(joinpath(dir, "right", "nested"))
            symlink(joinpath(dir, "right", "nested"), joinpath(dir, "left", "directory"))
            symlink("directory/../new", joinpath(dir, "left", "dangling"))
            symlink("dangling", joinpath(dir, "left", "chain"))
            target = joinpath(dir, "right", "new")
            for path in (
                joinpath(dir, "left", "directory", "..", "new"),
                joinpath(dir, "left", "dangling"),
                joinpath(dir, "left", "chain"),
            )
                write_value(path)
                @test isfile(target)
                @test !ispath(joinpath(dir, "left", "new"))
                @test Base.samefile(path, target)
                @test sort(readdir(joinpath(dir, "right"))) == ["nested", "new"]
                rm(target)
            end
            cd(joinpath(dir, "left")) do
                write_value("chain")
                @test Base.samefile("chain", target)
            end
            @test islink(joinpath(dir, "left", "dangling"))
            @test islink(joinpath(dir, "left", "chain"))
            symlink("b", joinpath(dir, "a"))
            symlink("a", joinpath(dir, "b"))
            err = publication_error(() -> write_value(joinpath(dir, "a")))
            @test err isa (caller === :async ? TaskFailedException : ArgumentError)
            @test islink(joinpath(dir, "a")) && islink(joinpath(dir, "b"))
        end
    end
    state, particles = publication_run()
    mktempdir() do dir
        target = joinpath(dir, "checkpoint")
        save_checkpoint(target, state, particles)
        before, entries = read(target), readdir(dir)
        state.time[] = NaN
        @test_throws ArgumentError save_checkpoint(target, state, particles)
        @test read(target) == before
        @test readdir(dir) == entries
        missing = joinpath(dir, "new")
        state.time[] = 0.0
        state.model = HybridModel(PublicationElectrons(0.1, true))
        @test_throws PublicationFailure save_checkpoint(missing, state, particles)
        @test !ispath(missing)
        @test readdir(dir) == entries
        state.model = HybridModel(PublicationElectrons(0.1, false))
        @test save_checkpoint(missing, state, particles) == missing
        @test read(missing) == before
    end
end
