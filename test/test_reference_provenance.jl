module ExternalReferenceProvenanceTests
using Test
include(joinpath(@__DIR__, "..", "validation", "reference_provenance.jl"))
const RP = ReferenceProvenance

function make_attempt(dir)
    root = joinpath(dir, "attempt with spaces")
    src = joinpath(root, "source")
    mkpath(src)
    git = Sys.which("git")
    quiet(cmd) = run(pipeline(cmd; stdout = devnull, stderr = devnull))
    quiet(`$git init $src`)
    write(joinpath(src, "input.dat"), "independent test input\n")
    quiet(`$git -C $src add input.dat`)
    quiet(
        `$git -C $src -c user.name=ReferenceTest -c user.email=reference@example.invalid -c commit.gpgsign=false commit -m fixture`,
    )
    revision = strip(read(`$git -C $src rev-parse HEAD`, String))
    identity = Dict(
        "source" => Dict("revision" => revision, "output" => "output.dat"),
        "tool" => "fixture",
    )
    write(joinpath(root, "output.dat"), "independent test output\n")
    record = Dict(
        "schema" => 1,
        "identity" => identity,
        "output" => "output.dat",
        "files" => RP._inventory(root),
    )
    RP._write_record(joinpath(root, "provenance.toml"), record)
    return root, identity, record
end

@testset "Pinned reference cache rejects altered artifacts" begin
    mktempdir() do dir
        root, identity, record = make_attempt(dir)
        result = RP._verify_attempt(root, identity)
        @test result.output == joinpath(root, "output.dat")
        @test read(result.output, String) == "independent test output\n"
        @test RP._read_record(result.provenance) == record
        original = read(result.provenance)
        for rel in ("output.dat", "source/input.dat")
            path = joinpath(root, rel)
            bytes = read(path)
            write(path, reverse(bytes))
            @test_throws RP.ReferenceFailure RP._verify_attempt(root, identity)
            write(path, bytes)
            @test RP._verify_attempt(root, identity) == result
            rm(path)
            @test_throws RP.ReferenceFailure RP._verify_attempt(root, identity)
            write(path, bytes)
        end
        write(joinpath(root, "extra.dat"), "unexpected artifact")
        @test_throws RP.ReferenceFailure RP._verify_attempt(root, identity)
        rm(joinpath(root, "extra.dat"))
        other = deepcopy(identity)
        other["tool"] = "changed compiler"
        @test_throws RP.ReferenceFailure RP._verify_attempt(root, other)
        for mutation in (
            Dict("schema" => 99),
            Dict("output" => "../outside"),
            Dict("output" => "source/input.dat"),
            Dict("files" => Dict()),
        )
            RP._write_record(result.provenance, merge(record, mutation))
            @test_throws RP.ReferenceFailure RP._verify_attempt(root, identity)
        end
        write(result.provenance, "invalid = [")
        @test_throws RP.ReferenceFailure RP._verify_attempt(root, identity)
        write(result.provenance, original)
        for path in ("", "/outside", "../outside", "source/../../outside", "./output.dat")
            @test_throws RP.ReferenceFailure RP._inside(root, path)
        end
        outside = joinpath(dir, "outside")
        write(outside, "external content")
        symlink(outside, joinpath(root, "linked.dat"))
        @test_throws RP.ReferenceFailure RP._inventory(root)
        rm(joinpath(root, "linked.dat"))
        mkdir(joinpath(dir, "external directory"))
        symlink(joinpath(dir, "external directory"), joinpath(root, "linked directory"))
        @test_throws RP.ReferenceFailure RP._inventory(root)
        rm(joinpath(root, "linked directory"))
        @test RP._verify_attempt(root, identity) == result
        @test length(readdir(root)) == 3
        git = Sys.which("git")
        source = joinpath(root, "source")
        run(
            pipeline(
                `$git -C $source -c user.name=ReferenceTest -c user.email=reference@example.invalid -c commit.gpgsign=false commit --allow-empty -m changed-revision`;
                stdout = devnull,
                stderr = devnull,
            ),
        )
        @test_throws RP.ReferenceFailure RP._verify_attempt(root, identity)
    end
end

@testset "Receipt publication preserves existing state on failure" begin
    mktempdir() do dir
        destination = joinpath(dir, "existing-directory")
        mkdir(destination)
        write(joinpath(destination, "retained"), "owned test content")
        @test_throws Base.IOError RP._write_record(destination, Dict("schema" => 1))
        @test read(joinpath(destination, "retained"), String) == "owned test content"
        @test readdir(dir) == ["existing-directory"]
        file = joinpath(dir, "current.toml")
        RP._write_record(file, Dict("generation" => 1))
        @test_throws ErrorException RP._write_record(file, Dict("invalid" => nothing))
        @test RP._read_record(file) == Dict("generation" => 1)
        RP._write_record(file, Dict("generation" => 2))
        @test RP._read_record(file) == Dict("generation" => 2)
        @test sort(readdir(dir)) == ["current.toml", "existing-directory"]
    end
end

@testset "External command failure retains its cause and log" begin
    mktempdir() do dir
        log = joinpath(dir, "failed.log")
        err = open(log, "w") do io
            try
                RP._run_logged(
                    `$(Base.julia_cmd()) --startup-file=no -e 'println(stderr,"reference failure witness"); exit(7)'`,
                    dir,
                    Dict(ENV),
                    io,
                    log,
                )
                nothing
            catch failure
                failure
            end
        end
        @test err isa RP.ReferenceFailure
        @test occursin("7", sprint(showerror, err))
        @test occursin(log, sprint(showerror, err))
        @test occursin("reference failure witness", read(log, String))
        @test !isfile(joinpath(dir, "provenance.toml"))
        script = joinpath(dir, "failing-tool")
        write(script, "#!/bin/sh\nprintf '%s\\n' 'incompatible runtime witness' >&2\nexit 23\n")
        chmod(script, 0o700)
        err = try
            RP._tool(script)
            nothing
        catch failure
            failure
        end
        @test err isa RP.ReferenceFailure
        @test occursin("incompatible runtime witness", sprint(showerror, err))
        @test occursin("23", sprint(showerror, err))
    end
end

@testset "Reference writer locks and missing tools fail before building" begin
    validation = abspath(joinpath(@__DIR__, "..", "validation"))
    @test_throws ArgumentError RP._request("unknown", validation)
    withenv("FC" => "unavailable-reference-compiler-4928") do
        @test_throws RP.ReferenceFailure RP._request("nhds", validation)
    end
    # Only preflight/lock/pointer paths run; no compiler or network is needed.
    withenv("FC" => Sys.which("git")) do
        mktempdir() do dir
            lock = joinpath(dir, ".writer-lock")
            mkdir(lock)
            write(joinpath(lock, "owner.toml"), "pid = 1\n")
            @test_throws RP.ReferenceFailure RP.ensure_reference("nhds"; cache_root = dir)
            @test read(joinpath(lock, "owner.toml"), String) == "pid = 1\n"
            rm(lock; recursive = true)
            RP._write_record(joinpath(dir, "current.toml"), Dict("attempt" => "../escape"))
            @test_throws RP.ReferenceFailure RP.ensure_reference("nhds"; cache_root = dir)
            @test !isdir(lock)
            @test readdir(dir) == ["current.toml"]
            identity = RP._request("nhds", validation)
            RP._write_record(
                joinpath(dir, "current.toml"),
                Dict(
                    "schema" => 1,
                    "attempt" => "../escape",
                    "identity_sha256" => RP._identity_digest(identity),
                ),
            )
            @test_throws RP.ReferenceFailure RP.ensure_reference("nhds"; cache_root = dir)
            @test !isdir(lock)
        end
    end
end
end # module
