# io_extra.jl — portable field dumps + asynchronous output staging (additive).
#
# Two independent I/O helpers:
#
#   write_field / read_field : a tiny, self-describing, cross-language-readable
#       binary container for dense Float64 arrays — no HDF5 / external dep. The
#       layout is fixed little-endian-on-host native order (Julia `write` uses
#       host byte order; the format records ndims+shape so any reader knowing
#       the (documented) host endianness can reconstruct the array). Round-trips
#       any Float64 AbstractArray exactly.
#
#   async_save : stage a checkpoint write off the calling thread via
#       `Threads.@spawn` + Serialization, returning the Task so the caller can
#       overlap I/O with computation and `wait` on it later.
#
# Serialization is already imported module-wide (checkpoint.jl); no new dep.

# On-disk magic + format version for write_field (8-byte tag, then a UInt8
# version). Lets a reader sanity-check the stream before trusting the sizes.
const _FIELD_MAGIC = b"HPSTDFLD"      # 8 bytes
const _FIELD_FORMAT_VERSION = UInt8(1)

"""
    write_field(path, A::AbstractArray) -> path

Write `A` to `path` as a portable, self-describing binary dump. The element
type is `Float64`: a non-`Float64` array is converted (`Float64.(A)`) before
writing, so the stored bytes are always 8-byte IEEE-754 doubles.

Layout (host byte order, which for all supported platforms is little-endian):

| bytes                | meaning                                   |
|----------------------|-------------------------------------------|
| 8                    | magic `"HPSTDFLD"`                         |
| 1  (`UInt8`)         | format version (`1`)                      |
| 8  (`Int64`)         | `ndims(A)` = `N`                           |
| 8·N (`Int64` each)   | `size(A)` per dimension                   |
| 8·length(A) (`Float64`) | array data in column-major order       |

The data is written in Julia's native column-major iteration order via the
`Array` itself, so `read_field` reconstructs the array exactly. Returns `path`.
"""
function write_field(path::AbstractString, A::AbstractArray)
    data = eltype(A) === Float64 ? A : Float64.(A)
    # Materialize to a dense Array so the on-disk byte order is column-major and
    # contiguous regardless of the input array type (views, ranges, etc.).
    arr = data isa Array{Float64} ? data : convert(Array{Float64}, data)
    open(path, "w") do io
        write(io, _FIELD_MAGIC)
        write(io, _FIELD_FORMAT_VERSION)
        write(io, Int64(ndims(arr)))
        for d in size(arr)
            write(io, Int64(d))
        end
        # Bulk write of the contiguous Float64 buffer.
        write(io, arr)
    end
    return path
end

"""
    read_field(path) -> Array{Float64}

Read a binary dump written by [`write_field`](@ref) and return the reconstructed
`Array{Float64}` (with the original number of dimensions and shape). Throws an
`ErrorException` if the magic tag or format version does not match, or if the
file is truncated relative to the declared shape.
"""
function read_field(path::AbstractString)
    open(path, "r") do io
        file_bytes = filesize(io)
        # Check the fixed header before typed reads, and bound allocations by
        # actual file bytes rather than trusting a possibly corrupt shape.
        if file_bytes < length(_FIELD_MAGIC) + sizeof(UInt8) + sizeof(Int64)
            error("read_field: truncated header in $(path)")
        end
        magic = read(io, length(_FIELD_MAGIC))
        if magic != _FIELD_MAGIC
            error("read_field: bad magic in $(path) (not a write_field dump)")
        end
        ver = read(io, UInt8)
        if ver != _FIELD_FORMAT_VERSION
            error("read_field: unsupported format version $(ver) in $(path)")
        end
        N = read(io, Int64)
        if N < 0
            error("read_field: corrupt ndims=$(N) in $(path)")
        end
        if N > (file_bytes - position(io)) ÷ sizeof(Int64)
            error("read_field: truncated dimensions in $(path)")
        end
        szs = Vector{Int}(undef, N)
        for i = 1:N
            szs[i] = Int(read(io, Int64))
            if szs[i] < 0
                error("read_field: corrupt size $(szs[i]) (dim $(i)) in $(path)")
            end
        end
        available = (file_bytes - position(io)) ÷ sizeof(Float64)
        # Empty fields need no payload. For nonempty shapes, dividing before
        # multiplying rejects truncated or overflowing extents without overflow.
        n = any(iszero, szs) ? 0 : 1
        if !iszero(n)
            for size in szs
                if n > available ÷ size
                    error("read_field: truncated data in $(path)")
                end
                n *= size
            end
        end
        if n > available
            error("read_field: truncated data in $(path)")
        end
        A = Array{Float64}(undef, szs...)
        read!(io, A)                          # fills column-major, matching write
        return A
    end
end

"""
    async_save(path, state) -> Task

Asynchronously serialize `state` to `path` on a `Threads.@spawn` task, returning
the `Task` immediately so the caller can continue computing and overlap the I/O.
`wait(task)` blocks until the write has completed (and rethrows any error raised
inside the task). After `wait`, `Serialization.deserialize(path)` returns a value
`==` to `state`.

This is deliberately lighter than [`save_run`](@ref): it stages a raw
Serialization dump with no schema/checksum wrapper, intended for transient
output staging where the goal is to get bytes off the hot path quickly.

`async_save` snapshots `state` before spawning the write, so mutations made by
the caller after `async_save(path, state)` returns do not race or leak into the
serialized file.

Publish atomically after writing and closing the snapshot. A failed write leaves
the previous file intact and fails the returned task. File permissions, symlinks,
hardlinks and directory-access requirements follow [`save_checkpoint`](@ref).
"""
function async_save(path::AbstractString, state)
    p = String(path)
    payload = deepcopy(state)
    return Threads.@spawn _atomic_serialize(p, payload)
end
