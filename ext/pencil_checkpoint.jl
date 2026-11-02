# Shared file/checksum/generation mechanics for the two actual pencil models.
# Model schemas, phase/configuration and save/load orchestration remain explicit.

function _pencil_checkpoint_generation_name(name, prefix, label)
    name = HybridPlasmaPIC._validate_mpi_checkpoint_filename(name)
    startswith(name, prefix) && length(name) > length(prefix) ||
        throw(ArgumentError("invalid $label checkpoint generation"))
    return name
end

function _pencil_checkpoint_generation_path(dir, name, prefix, label)
    location = joinpath(dir, _pencil_checkpoint_generation_name(name, prefix, label))
    !islink(location) && isdir(location) ||
        throw(ArgumentError("$label checkpoint generation must be an ordinary directory"))
    return location
end

function _pencil_checkpoint_digest(value, label)
    value isa String &&
        ncodeunits(value) == 64 &&
        all(c -> '0' <= c <= '9' || 'a' <= c <= 'f', value) ||
        throw(ArgumentError("$label checksum must be a lowercase SHA-256 digest"))
    return value
end

function _pencil_checkpoint_file_digest(path, label)
    !islink(path) && isfile(path) ||
        throw(ArgumentError("$label rank file must be an ordinary file"))
    return open(path, "r") do io
        bytes2hex(HybridPlasmaPIC.SHA.sha256(io))
    end
end

function _pencil_checkpoint_broadcast_generation(name, comm, rank, prefix, label)
    context = "$label checkpoint"
    encoded = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation encoding") do
        bytes = HybridPlasmaPIC._serialize_mpi_payload(name)
        HybridPlasmaPIC._check_mpi_count(length(bytes), "$label generation name")
        bytes
    end
    count = MPI.Bcast(Cint(length(encoded)), 0, comm)
    bytes = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation receive allocation") do
        rank == 0 ? encoded : Vector{UInt8}(undef, count)
    end
    MPI.Bcast!(bytes, 0, comm)
    return HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation decoding") do
        decoded = HybridPlasmaPIC._deserialize_mpi_payload(String, bytes)
        _pencil_checkpoint_generation_name(decoded, prefix, label)
    end
end

function _pencil_checkpoint_remove_generation(dir, generation, nranks, prefix, label)
    location = _pencil_checkpoint_generation_path(dir, generation, prefix, label)
    files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    for name in files
        path = joinpath(location, name)
        islink(path) && throw(ArgumentError("refusing to clean a linked $label rank file"))
        !ispath(path) ||
            isfile(path) ||
            throw(ArgumentError("refusing to clean a non-file $label rank entry"))
    end
    for name in files
        rm(joinpath(location, name); force = true)
    end
    # Nonrecursive removal preserves unknown entries and reports their presence.
    rm(location)
    return nothing
end

function _pencil_checkpoint_rank_state(path, digest, label)
    !islink(path) && isfile(path) ||
        throw(ArgumentError("$label rank file must be an ordinary file"))
    return open(path, "r") do io
        bytes2hex(HybridPlasmaPIC.SHA.sha256(io)) == digest ||
            throw(ArgumentError("$label rank checksum mismatch"))
        seekstart(io)
        HybridPlasmaPIC.Serialization.deserialize(io)
    end
end
