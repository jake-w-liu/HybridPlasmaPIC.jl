# Resolve physical parents before following a leaf link. Lexically normalizing
# `..` first would change paths that cross a directory symlink. The final leaf
# may not exist yet, including when it is the referent of a dangling link.
function _serial_publication_target(path::AbstractString)
    seen = Set{String}()
    while true
        directory = dirname(path)
        directory = realpath(isempty(directory) ? "." : directory)
        target = joinpath(directory, basename(path))
        islink(target) || return target
        target in seen && throw(ArgumentError("checkpoint destination has a symbolic link cycle"))
        push!(seen, target)
        linked = readlink(target)
        path = isabspath(linked) ? linked : joinpath(directory, linked)
    end
end

# Check the entry itself: following a replacement symlink would authorize an
# unrelated file. Missing and replaced entries are never publication candidates.
function _require_publication_file(temporary, identity)
    current = lstat(temporary)
    isfile(current) && current.device == identity.device && current.inode == identity.inode ||
        throw(ArgumentError("serialized staging file is missing or has been replaced"))
    return nothing
end

# This function owns the file-backed stream and the original staging inode.
# The caller must exclude concurrent directory mutation; identity checks also
# prevent serializer/IO callbacks from substituting another entry for that file.
function _publish_serialized(write_contents, target, temporary, io::IO, mode)
    identity = nothing
    try
        identity = stat(io)
        _require_publication_file(temporary, identity)
        # Set permissions before invoking serializer callbacks, which may change
        # filesystem entries. Base provides pathname-based chmod, not fchmod.
        chmod(temporary, mode)
        result = write_contents(io)
        flush(io)
        close(io)
        _require_publication_file(temporary, identity)
        # Same-directory rename publishes the complete file without a delete or
        # copy fallback. No fallible operation follows the publication point.
        Base.Filesystem.rename(temporary, target)
        return result
    catch original
        errors = CapturedException[CapturedException(original, catch_backtrace())]
        try
            isopen(io) && close(io)
        catch cleanup_error
            push!(errors, CapturedException(cleanup_error, catch_backtrace()))
        end
        try
            # A replacement at our former path is not owned by this call. Keep
            # it intact and report the ownership loss instead of removing it.
            current = lstat(temporary)
            if ispath(current)
                identity === nothing &&
                    throw(ArgumentError("serialized staging file identity is unavailable"))
                _require_publication_file(temporary, identity)
                rm(temporary)
            end
        catch cleanup_error
            push!(errors, CapturedException(cleanup_error, catch_backtrace()))
        end
        length(errors) == 1 && rethrow()
        throw(CompositeException(errors))
    end
end

function _atomic_serialize(path::AbstractString, value)
    target = _serial_publication_target(path)
    mode = if ispath(target)
        isfile(target) || throw(ArgumentError("checkpoint destination must be a regular file"))
        # Rename only checks directory permissions. Retain the existing writer's
        # file-access requirement, without truncating or changing its contents.
        open(target; write = true, create = false, truncate = false) do io
            filemode(io) & 0o777
        end
    else
        0o600
    end
    temporary, io = mktemp(dirname(target); cleanup = false)
    return _publish_serialized(target, temporary, io, mode) do stream
        Serialization.serialize(stream, value)
    end
end
