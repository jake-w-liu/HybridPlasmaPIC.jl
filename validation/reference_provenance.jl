module ReferenceProvenance

using SHA
using TOML

export ensure_reference, ReferenceFailure

struct ReferenceFailure <: Exception
    message::String
end
Base.showerror(io::IO, e::ReferenceFailure) = print(io, e.message)

_digest(path) = open(sha256, path) |> bytes2hex
_toml(data) = sprint(io -> TOML.print(io, data; sorted = true))
_identity_digest(data) = bytes2hex(sha256(_toml(data)))

function _inside(root, relative)
    parts = splitpath(relative)
    (isempty(parts) || isabspath(relative) || any(p -> p in (".", ".."), parts)) &&
        throw(ReferenceFailure("Invalid relative artifact path: $(repr(relative))"))
    path = joinpath(root, relative)
    # Existing symlinks must not redirect verification outside this attempt.
    if ispath(path) &&
       !startswith(realpath(path), realpath(root) * string(Base.Filesystem.path_separator))
        throw(ReferenceFailure("Artifact escapes its attempt directory: $relative"))
    end
    return path
end

function _inventory(root)
    files = Dict{String,Any}()
    for (dir, dirs, names) in walkdir(root)
        if dir == joinpath(root, "source")
            filter!(!=(".git"), dirs)
        end
        for name in dirs
            islink(joinpath(dir, name)) && throw(
                ReferenceFailure(
                    "Artifact directory is a symlink: $(relpath(joinpath(dir, name), root))",
                ),
            )
        end
        for name in names
            relative = relpath(joinpath(dir, name), root)
            relative == "provenance.toml" && continue
            path = _inside(root, relative)
            isfile(path) || throw(ReferenceFailure("Artifact is not a regular file: $relative"))
            files[relative] = Dict(
                "sha256" => _digest(path),
                "bytes" => filesize(path),
                "mode" => Int(filemode(path) & 0o777),
                "link" => islink(path) ? readlink(path) : "",
            )
        end
    end
    return files
end

function _read_record(path)
    isfile(path) || throw(ReferenceFailure("Missing provenance record: $path"))
    try
        return TOML.parsefile(path)
    catch err
        err isa TOML.ParserError || rethrow()
        throw(ReferenceFailure("Invalid provenance record $path: $(sprint(showerror, err))"))
    end
end

function _write_record(path, data)
    # A unique sibling temp file permits atomic publication on the same filesystem.
    tmp, io = mktemp(dirname(path))
    try
        write(io, _toml(data))
        close(io)
        # Match the package publication owner: no delete/copy fallback on error.
        Base.Filesystem.rename(tmp, path)
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp)
    end
    return path
end

function _verify_attempt(root, identity; git = "git")
    record = _read_record(joinpath(root, "provenance.toml"))
    get(record, "schema", nothing) == 1 ||
        throw(ReferenceFailure("Unsupported provenance schema in $root"))
    get(record, "identity", nothing) == identity ||
        throw(ReferenceFailure("Reference identity changed in $root"))
    files = _inventory(root)
    get(record, "files", nothing) == files || throw(
        ReferenceFailure(
            "Reference artifacts changed in $root; preserve this attempt and rebuild with a fresh cache directory",
        ),
    )
    source = joinpath(root, "source")
    revision = strip(read(`$git -C $source rev-parse HEAD`, String))
    revision == identity["source"]["revision"] ||
        throw(ReferenceFailure("Source revision changed in $root"))
    output = get(record, "output", "")
    output == identity["source"]["output"] ||
        throw(ReferenceFailure("Receipt selects the wrong output in $root"))
    haskey(files, output) ||
        throw(ReferenceFailure("Output is absent from the artifact manifest in $root"))
    files[output]["bytes"] > 0 || throw(ReferenceFailure("Empty reference output in $root"))
    return (output = _inside(root, output), provenance = joinpath(root, "provenance.toml"))
end

function _tool(name; version = true)
    path = Sys.which(name)
    isnothing(path) && throw(
        ReferenceFailure(
            "Required external-reference tool '$name' is unavailable; install it explicitly before running this case",
        ),
    )
    path = abspath(path)
    info = Dict("path" => path, "realpath" => realpath(path), "sha256" => _digest(path))
    if version
        info["version"] = mktemp() do _, io
            try
                run(pipeline(setenv(`$path --version`, _build_environment()); stdout = io, stderr = io))
            catch err
                err isa ProcessFailedException || rethrow()
                seekstart(io)
                detail = read(io, String)
                throw(
                    ReferenceFailure(
                        "Cannot identify tool $path: $(sprint(showerror, err))\n$detail",
                    ),
                )
            end
            seekstart(io)
            read(io, String)
        end
    end
    return info
end

function _build_environment()
    # The build receives only these variables, so unrelated shell flags cannot
    # silently alter a cached reference. The effective environment is recorded.
    env = Dict(
        "PATH" => ENV["PATH"],
        "HOME" => homedir(),
        "LC_ALL" => "C",
        "LANG" => "C",
        "GIT_CONFIG_NOSYSTEM" => "1",
        "GIT_CONFIG_GLOBAL" => "/dev/null",
    )
    for key in ("TMPDIR", "SDKROOT", "SYSTEMROOT")
        haskey(ENV, key) && (env[key] = ENV[key])
    end
    return env
end

function _request(id, validation_dir)
    pins = TOML.parsefile(joinpath(validation_dir, "reference_sources.toml"))
    pins["schema"] == 1 || throw(ArgumentError("Unsupported reference lock schema"))
    id in ("nhds", "hybridvpic") || throw(ArgumentError("Unknown reference: $id"))
    source = pins[id]
    occursin(r"^[0-9a-f]{40}$", source["revision"]) ||
        throw(ArgumentError("Reference revision must be a full Git commit"))
    tools = Dict("git" => _tool("git"), "sh" => _tool("sh"; version = false))
    inputs = Dict{String,String}()
    if id == "nhds"
        tools["fc"] = _tool(get(ENV, "FC", "gfortran"))
        tools["make"] = _tool("make")
    else
        tools["cc"] = _tool(get(ENV, "MPICC", "mpicc"))
        tools["cxx"] = _tool(get(ENV, "MPICXX", "mpicxx"))
        tools["cmake"] = _tool("cmake")
        tools["mpi"] = _tool(get(ENV, "MPIEXEC", "mpirun"))
        inputs["shock_perp.cxx"] =
            _digest(joinpath(validation_dir, source["case_directory"], "shock_perp.cxx"))
    end
    return Dict(
        "reference" => id,
        "source" => source,
        "tools" => tools,
        "inputs" => inputs,
        "environment" => _build_environment(),
        "builder_sha256" => _digest(@__FILE__),
        "machine" => Sys.MACHINE,
        "kernel" => string(Sys.KERNEL),
        "cpu" => Sys.CPU_NAME,
    )
end

function _run_logged(cmd, dir, env, io, logpath)
    println(io, "\nDirectory: ", dir, "\nCommand: ", repr(cmd))
    flush(io)
    try
        run(pipeline(setenv(Cmd(cmd; dir = dir), env); stdout = io, stderr = io))
    catch err
        err isa ProcessFailedException || rethrow()
        throw(
            ReferenceFailure(
                "External reference command failed: $(sprint(showerror, err)). Full log: $logpath",
            ),
        )
    end
end

function _build_attempt(root, identity, validation_dir)
    source = identity["source"]
    tools = identity["tools"]
    env = copy(identity["environment"])
    src = joinpath(root, "source")
    build = joinpath(root, "build")
    rundir = joinpath(root, "run")
    mkdir(src)
    mkdir(build)
    mkdir(rundir)
    if identity["reference"] == "nhds"
        # Autoconf treats FC as command text. A fixed launcher preserves executable
        # paths containing spaces or shell syntax without reinterpreting the path.
        bindir = joinpath(build, "reference-bin")
        mkdir(bindir)
        launcher = joinpath(bindir, "reference-fortran")
        write(launcher, "#!/bin/sh\nexec \"\$REFERENCE_FORTRAN_COMPILER\" \"\$@\"\n")
        chmod(launcher, 0o700)
        env["PATH"] = bindir * string(Sys.iswindows() ? ';' : ':') * env["PATH"]
        env["REFERENCE_FORTRAN_COMPILER"] = tools["fc"]["path"]
    end
    git = tools["git"]["path"]
    logpath = joinpath(root, "build.log")
    output = source["output"]
    open(logpath, "w") do io
        println(io, "Effective build environment:\n", _toml(env))
        execute(cmd, dir = root) = _run_logged(cmd, dir, env, io, logpath)
        execute(`$git init $src`)
        execute(`$git -C $src fetch --depth 1 $(source["url"]) $(source["revision"])`)
        execute(`$git -C $src checkout --detach FETCH_HEAD`)
        revision = strip(read(`$git -C $src rev-parse HEAD`, String))
        revision == source["revision"] ||
            throw(ReferenceFailure("Fetched source does not match pinned revision; see $logpath"))
        license = _inside(src, source["license_path"])
        _digest(license) == source["license_sha256"] ||
            throw(ReferenceFailure("Pinned license checksum mismatch; see $logpath"))
        if identity["reference"] == "nhds"
            # Upstream configure permits spaces in cwd, but rejects them in srcdir.
            configure = joinpath("..", "source", "configure")
            sh = tools["sh"]["path"]
            fc = "reference-fortran"
            make = tools["make"]["path"]
            execute(`$sh $configure FC=$fc FCFLAGS=`, build)
            # Upstream Fortran module dependencies require serial make.
            execute(`$make -j1`, build)
            cp(joinpath(src, "whistler.in"), joinpath(rundir, "whistler.in"))
            execute(`$(joinpath(build, "src", "NHDS")) whistler.in`, rundir)
        else
            deck = joinpath(rundir, "shock_perp.cxx")
            cp(joinpath(validation_dir, source["case_directory"], "shock_perp.cxx"), deck)
            cp(
                joinpath(src, "examples", "shock", "injection.cxx"),
                joinpath(rundir, "injection.cxx"),
            )
            cmake = tools["cmake"]["path"]
            cc = tools["cc"]["path"]
            cxx = tools["cxx"]["path"]
            execute(
                `$cmake -S $src -B $build -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release -DENABLE_INTEGRATED_TESTS=OFF -DCMAKE_C_COMPILER=$cc -DCMAKE_CXX_COMPILER=$cxx -DCMAKE_CXX_FLAGS=-fno-strict-aliasing -DUSER_DECKS=ON -DUSER_DECK=$deck`,
            )
            execute(`$cmake --build $build --target shock_perp --parallel 4`)
            mpi = tools["mpi"]["path"]
            execute(`$mpi -n 1 $(joinpath(build, "shock_perp"))`, rundir)
        end
    end
    open(logpath, "a") do io
        _run_logged(`$git -C $src diff --no-ext-diff --exit-code HEAD --`, root, env, io, logpath)
    end
    path = _inside(root, output)
    isfile(path) && filesize(path) > 0 ||
        throw(ReferenceFailure("External run produced no output: $path; see $logpath"))
    # Detect a deck/tool change during the build before publishing its receipt.
    _request(identity["reference"], validation_dir) == identity || throw(
        ReferenceFailure(
            "Reference inputs or tools changed during the build; retained attempt: $root",
        ),
    )
    record = Dict(
        "schema" => 1,
        "identity" => identity,
        "output" => output,
        "julia_version" => string(VERSION),
        "files" => _inventory(root),
    )
    _write_record(joinpath(root, "provenance.toml"), record)
    return _verify_attempt(root, identity; git)
end

"""Build or verify a pinned external reference, returning output and provenance paths.

Legacy caches are untouched. A changed configuration gets a new attempt; a
corrupt current attempt raises `ReferenceFailure` and is retained for inspection.
"""
function ensure_reference(id::AbstractString; cache_root = nothing, validation_dir = @__DIR__)
    identity = _request(String(id), validation_dir)
    base =
        isnothing(cache_root) ?
        joinpath(validation_dir, identity["source"]["case_directory"], "reference_cache") :
        abspath(cache_root)
    mkpath(base)
    lockdir = joinpath(base, ".writer-lock")
    try
        mkdir(lockdir)
    catch err
        err isa Base.IOError && isdir(lockdir) || rethrow()
        throw(
            ReferenceFailure(
                "Reference writer lock exists: $lockdir. Verify its owner before removing a stale lock",
            ),
        )
    end
    try
        _write_record(
            joinpath(lockdir, "owner.toml"),
            Dict("pid" => getpid(), "reference" => String(id)),
        )
        current = joinpath(base, "current.toml")
        key = _identity_digest(identity)
        if isfile(current)
            pointer = _read_record(current)
            get(pointer, "schema", nothing) == 1 &&
                get(pointer, "identity_sha256", nothing) isa String &&
                get(pointer, "attempt", nothing) isa String ||
                throw(ReferenceFailure("Invalid cache pointer: $current"))
            if get(pointer, "identity_sha256", nothing) == key
                root = _inside(base, get(pointer, "attempt", ""))
                return _verify_attempt(root, identity; git = identity["tools"]["git"]["path"])
            end
        end
        root = mktempdir(base; prefix = "attempt-", cleanup = false)
        _write_record(joinpath(root, "request.toml"), identity)
        result = _build_attempt(root, identity, validation_dir)
        _write_record(
            current,
            Dict("schema" => 1, "identity_sha256" => key, "attempt" => basename(root)),
        )
        return result
    finally
        # This invocation created this lock; no other process-owned state is removed.
        rm(lockdir; recursive = true)
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 2) ||
        error("Usage: reference_provenance.jl nhds|hybridvpic [cache-directory]")
    result = ReferenceProvenance.ensure_reference(
        ARGS[1];
        cache_root = length(ARGS) == 2 ? ARGS[2] : nothing,
    )
    println("Verified output: ", result.output)
    println("Provenance: ", result.provenance)
end
