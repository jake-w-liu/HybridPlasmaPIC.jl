# metadata.jl — run provenance + checkpoint schema/checksum (purely additive).
#
# Captures the information needed to reproduce and audit a run (git commit,
# Julia version, environment hashes, RNG seed, normalization/filter/boundary/
# diagnostic descriptions, timestamp), and wraps Serialization-based save/load
# with a schema-version tag and an integrity digest over the stored state.
#
# Serialization is already imported by checkpoint.jl (module-wide); no new dep.

using Dates
using SHA

struct _SHA256WriteSink <: IO
    ctx::SHA.SHA256_CTX
end

Base.iswritable(::_SHA256WriteSink) = true

function Base.write(s::_SHA256WriteSink, x::UInt8)
    SHA.update!(s.ctx, UInt8[x])
    return 1
end

function Base.unsafe_write(s::_SHA256WriteSink, p::Ptr{UInt8}, n::UInt)
    bytes = unsafe_wrap(Vector{UInt8}, p, Int(n))
    SHA.update!(s.ctx, bytes)
    return n
end

"""
    RunMetadata

Provenance record for a simulation run. All fields are `String` except `rng_seed`.

Fields:
- `git_commit`   — `git rev-parse HEAD`, or `"unknown"` if unavailable.
- `julia_version`— `string(VERSION)`.
- `project_hash` — SHA-256 digest of `Project.toml` contents, or `"absent"`.
- `manifest_hash`— SHA-256 digest of `Manifest.toml` contents, or `"absent"`.
- `rng_seed`     — the RNG seed used to make the run reproducible.
- `normalization`— units convention (default `"Omega_ci"`).
- `filter_desc`, `boundary_desc`, `diagnostic_desc` — free-form descriptions.
- `timestamp`    — ISO-ish timestamp string (caller-supplied or auto).
- `backend`      — compute backend label (e.g. `"CPU"`).
- `hardware`     — hardware description: CPU model + logical thread count.
- `rank_layout`  — distributed-memory layout description; serial runs record
  `ranks=1`.
"""
struct RunMetadata
    git_commit::String
    julia_version::String
    project_hash::String
    manifest_hash::String
    rng_seed::Int
    normalization::String
    filter_desc::String
    boundary_desc::String
    diagnostic_desc::String
    timestamp::String
    backend::String
    hardware::String
    rank_layout::String
end

# Stable content digest of a file's bytes, or "absent" if the file does not exist.
function _file_content_hash(path::AbstractString)
    isfile(path) || return "absent"
    return bytes2hex(SHA.sha256(read(path)))
end

function _state_checksum(state)
    ctx = SHA.SHA256_CTX()
    serialize(_SHA256WriteSink(ctx), state)
    return bytes2hex(SHA.digest!(ctx))
end

# Locate the package root = nearest ancestor directory containing Project.toml.
# Robust to this file's depth in the src/ tree (it lives at src/io/), so the
# reorganization into §5.3 subdirectories cannot break the Project/Manifest hash.
function _pkg_root()
    d = @__DIR__
    for _ = 1:8
        isfile(joinpath(d, "Project.toml")) && return d
        nd = dirname(d)
        nd == d && break          # reached the filesystem root
        d = nd
    end
    return dirname(@__DIR__)        # fallback (matches the original flat-src behavior)
end

_serial_rank_layout() = "serial;ranks=1;dims=(1);coords=(0);periodic=false"

function _validated_rng_seed(rng_seed::Integer)
    typemin(Int) <= rng_seed <= typemax(Int) ||
        throw(ArgumentError("rng_seed must fit in Int, got $rng_seed"))
    return Int(rng_seed)
end

"""
    capture_metadata(; rng_seed, normalization="Omega_ci", filter_desc="",
                       boundary_desc="", diagnostic_desc="", timestamp="",
                       rank_layout="")

Build a [`RunMetadata`] for the current process. `git_commit` is read from
`git rev-parse HEAD` (falls back to `"unknown"` outside a repo / no git).
`project_hash`/`manifest_hash` are SHA-256 digests of this package's
`Project.toml`/`Manifest.toml` if present (else `"absent"`). An empty
`timestamp` is auto-filled with the current UTC time as an ISO-ish string ending
in `Z`.

`rank_layout` records the MPI/distributed-memory topology used by the run. Leave
it empty for the serial default, which is recorded as `ranks=1` rather than
omitted.
"""
function capture_metadata(;
    rng_seed::Integer,
    normalization::AbstractString = "Omega_ci",
    filter_desc::AbstractString = "",
    boundary_desc::AbstractString = "",
    diagnostic_desc::AbstractString = "",
    timestamp::AbstractString = "",
    rank_layout::AbstractString = "",
)
    root = _pkg_root()
    git_commit = try
        readchomp(pipeline(`git -C $root rev-parse HEAD`; stderr = devnull))
    catch
        "unknown"
    end
    # Guard against an empty success (e.g. git returns nothing useful).
    isempty(git_commit) && (git_commit = "unknown")

    project_hash = _file_content_hash(joinpath(root, "Project.toml"))
    manifest_hash = _file_content_hash(joinpath(root, "Manifest.toml"))

    ts =
        isempty(timestamp) ?
        string(Dates.format(Dates.unix2datetime(time()), dateformat"yyyy-mm-ddTHH:MM:SS"), 'Z') :
        String(timestamp)

    # Backend/hardware provenance. This solver runs on the CPU; the hardware
    # string combines the CPU model with the logical-thread count. Guard the
    # cpu_info() lookup so a runtime that returns an empty vector (rare, but not
    # impossible on exotic platforms) still yields a nonempty hardware string.
    backend = "CPU"
    cpu_model = try
        info = Sys.cpu_info()
        isempty(info) ? "unknown CPU" : strip(info[1].model)
    catch
        "unknown CPU"
    end
    isempty(cpu_model) && (cpu_model = "unknown CPU")
    hardware = string(cpu_model, " x", Sys.CPU_THREADS)
    layout = isempty(rank_layout) ? _serial_rank_layout() : String(rank_layout)

    return RunMetadata(
        String(git_commit),
        string(VERSION),
        project_hash,
        manifest_hash,
        _validated_rng_seed(rng_seed),
        String(normalization),
        String(filter_desc),
        String(boundary_desc),
        String(diagnostic_desc),
        ts,
        backend,
        hardware,
        layout,
    )
end

"Schema version of the checkpoint container written by [`save_run`](@ref)."
const CHECKPOINT_SCHEMA_VERSION = 3

"""
    save_run(path, state, meta::RunMetadata)

Serialize a tagged container to `path`: a `NamedTuple`
`(schema, meta, state, checksum)` where `schema == CHECKPOINT_SCHEMA_VERSION`
and `checksum` is a SHA-256 digest of the serialized `state`. `state` is any
serializable value (typically a `NamedTuple` of simulation arrays/scalars).
Returns `path`.
"""
function save_run(path::AbstractString, state, meta::RunMetadata)
    container =
        (schema = CHECKPOINT_SCHEMA_VERSION, meta = meta, state = state, checksum = _state_checksum(state))
    serialize(path, container)
    return path
end

"""
    load_run(path) -> (; schema, meta, state)

Deserialize a container written by [`save_run`](@ref). Verifies that
`schema == CHECKPOINT_SCHEMA_VERSION` and that the recomputed SHA-256 state
digest matches the stored `checksum`; throws an `ErrorException` on either
mismatch (version skew or corruption).
"""
function load_run(path::AbstractString)
    container = deserialize(path)
    # Validate that the deserialized object has the expected shape.
    if !(container isa NamedTuple) ||
       !all(k -> hasproperty(container, k), (:schema, :meta, :state, :checksum))
        error("load_run: file $(path) is not a valid run container")
    end
    if container.schema != CHECKPOINT_SCHEMA_VERSION
        error(
            "load_run: schema version $(container.schema) ≠ expected " *
            "$(CHECKPOINT_SCHEMA_VERSION)",
        )
    end
    if !(container.meta isa RunMetadata)
        error("load_run: metadata in $(path) is not a RunMetadata record")
    end
    if !(container.checksum isa AbstractString)
        error("load_run: checksum in $(path) is not a SHA-256 digest string")
    end
    checksum = _state_checksum(container.state)
    if checksum != container.checksum
        error(
            "load_run: checksum mismatch — state is corrupted " *
            "(stored $(container.checksum), recomputed $(checksum))",
        )
    end
    return (schema = container.schema, meta = container.meta, state = container.state)
end
