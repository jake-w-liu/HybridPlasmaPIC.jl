const _PENCIL_MAXWELL_MANIFEST = "pencil_maxwell_manifest.ser"
const _PENCIL_MAXWELL_GENERATION = "_pencil_maxwell_"

function _pencil_maxwell_checkpoint_config(model::DistributedMaxwellVacuum{T}, nranks) where {T}
    return (
        schema = 1,
        model = :distributed_maxwell_vacuum,
        integrator_policy = :synchronized_spectral_kdk_exact_affine_v1,
        derivative_policy = :represented_odd_nyquist_zero_curl_v1,
        layout_policy = :matching_physical_pencil_v1,
        T = T,
        n = model.g.n,
        L = model.g.L,
        dx = model.g.dx,
        c = model.c,
        nranks = nranks,
    )
end

function _pencil_maxwell_checkpoint_agreement(model::DistributedMaxwellVacuum{T}, dir, nranks) where {T}
    return _agree_configuration(model.g.fft.comm, "Pencil Maxwell checkpoint") do
        _distributed_maxwell_mode(T)
        HybridPlasmaPIC._require_finite_positive_real("c", model.c, T)
        dir isa AbstractString || throw(ArgumentError("checkpoint directory must be a string"))
        !isempty(dir) || throw(ArgumentError("checkpoint directory must be nonempty"))
        String(dir), _pencil_maxwell_checkpoint_config(model, nranks)
    end
end

function _pencil_maxwell_checkpoint_keys(record, required)
    HybridPlasmaPIC._validate_checkpoint_container(record, required)
    length(keys(record)) == length(required) ||
        throw(ArgumentError("Pencil Maxwell checkpoint contains unsupported fields"))
    return nothing
end

function _pencil_maxwell_checkpoint_common(record, model::DistributedMaxwellVacuum{T}, nranks) where {T}
    expected = _pencil_maxwell_checkpoint_config(model, nranks)
    for name in keys(expected)
        value, reference = getproperty(record, name), getproperty(expected, name)
        typeof(value) === typeof(reference) && isequal(value, reference) ||
            throw(ArgumentError("Pencil Maxwell checkpoint $name does not match the model"))
    end
    record.time isa T && isfinite(record.time) && record.time >= zero(T) ||
        throw(ArgumentError("Pencil Maxwell checkpoint time must be finite, nonnegative and $T"))
    record.step isa Int && record.step >= 0 ||
        throw(ArgumentError("Pencil Maxwell checkpoint step must be a nonnegative Int"))
    return nothing
end

function _pencil_maxwell_checkpoint_layout(model::DistributedMaxwellVacuum)
    field = model.fields.E[1]
    ranges = PencilArrays.range_local(field, PencilArrays.LogicalOrder())
    return (
        logical_ranges = ntuple(d -> (Int(first(ranges[d])), Int(last(ranges[d]))), 3),
        permutation = PencilArrays.permutation(field) * (1, 2, 3),
        parent_size = size(parent(field)),
    )
end

function _pencil_maxwell_generation_name(name)
    name = HybridPlasmaPIC._validate_mpi_checkpoint_filename(name)
    startswith(name, _PENCIL_MAXWELL_GENERATION) &&
        length(name) > length(_PENCIL_MAXWELL_GENERATION) ||
        throw(ArgumentError("invalid Pencil Maxwell checkpoint generation"))
    return name
end

function _pencil_maxwell_generation_path(dir, name)
    location = joinpath(dir, _pencil_maxwell_generation_name(name))
    !islink(location) && isdir(location) ||
        throw(ArgumentError("Pencil Maxwell checkpoint generation must be an ordinary directory"))
    return location
end

function _pencil_maxwell_digest(value)
    value isa String && ncodeunits(value) == 64 &&
        all(c -> '0' <= c <= '9' || 'a' <= c <= 'f', value) ||
        throw(ArgumentError("Pencil Maxwell checksum must be a lowercase SHA-256 digest"))
    return value
end

function _pencil_maxwell_file_digest(path)
    !islink(path) && isfile(path) ||
        throw(ArgumentError("Pencil Maxwell rank file must be an ordinary file"))
    return open(path, "r") do io
        bytes2hex(HybridPlasmaPIC.SHA.sha256(io))
    end
end

function _pencil_maxwell_manifest(path, model, nranks)
    envelope = HybridPlasmaPIC.Serialization.deserialize(path)
    _pencil_maxwell_checkpoint_keys(envelope, (:schema, :model, :payload, :sha256))
    envelope.schema === 1 && envelope.model === :pencil_maxwell_manifest ||
        throw(ArgumentError("unsupported Pencil Maxwell manifest envelope"))
    envelope.payload isa Vector{UInt8} ||
        throw(ArgumentError("Pencil Maxwell manifest payload must be serialized bytes"))
    digest = _pencil_maxwell_digest(envelope.sha256)
    bytes2hex(HybridPlasmaPIC.SHA.sha256(envelope.payload)) == digest ||
        throw(ArgumentError("Pencil Maxwell manifest checksum mismatch"))
    record = HybridPlasmaPIC._deserialize_mpi_payload(NamedTuple, envelope.payload)
    required = (
        keys(_pencil_maxwell_checkpoint_config(model, nranks))...,
        :time,
        :step,
        :generation,
        :rank_files,
        :rank_sha256,
    )
    _pencil_maxwell_checkpoint_keys(record, required)
    _pencil_maxwell_checkpoint_common(record, model, nranks)
    record.generation isa String || throw(ArgumentError("checkpoint generation must be a String"))
    _pencil_maxwell_generation_name(record.generation)
    expected_files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    record.rank_files isa Tuple && isequal(record.rank_files, expected_files) ||
        throw(ArgumentError("Pencil Maxwell manifest requires canonical ordered rank filenames"))
    record.rank_sha256 isa Tuple && length(record.rank_sha256) == nranks ||
        throw(ArgumentError("Pencil Maxwell manifest requires one digest per rank"))
    foreach(_pencil_maxwell_digest, record.rank_sha256)
    return record
end

function _pencil_maxwell_envelope(manifest)
    payload = HybridPlasmaPIC._serialize_mpi_payload(manifest)
    HybridPlasmaPIC._check_mpi_count(length(payload), "Pencil Maxwell manifest payload")
    return (
        schema = 1,
        model = :pencil_maxwell_manifest,
        payload = payload,
        sha256 = bytes2hex(HybridPlasmaPIC.SHA.sha256(payload)),
    )
end

# The fixed generation-name broadcast cannot use unguarded MPI.bcast: encoding,
# receiver allocation and decoding must agree before the next collective.
function _pencil_maxwell_broadcast_generation(name, comm, rank)
    context = "Pencil Maxwell checkpoint"
    encoded = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation encoding") do
        bytes = HybridPlasmaPIC._serialize_mpi_payload(name)
        HybridPlasmaPIC._check_mpi_count(length(bytes), "Pencil Maxwell generation name")
        bytes
    end
    count = MPI.Bcast(Cint(length(encoded)), 0, comm)
    bytes = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation receive allocation") do
        rank == 0 ? encoded : Vector{UInt8}(undef, count)
    end
    MPI.Bcast!(bytes, 0, comm)
    return HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation decoding") do
        decoded = HybridPlasmaPIC._deserialize_mpi_payload(String, bytes)
        _pencil_maxwell_generation_name(decoded)
    end
end

function _pencil_maxwell_remove_generation(dir, generation, nranks)
    location = _pencil_maxwell_generation_path(dir, generation)
    files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    for name in files
        path = joinpath(location, name)
        islink(path) && throw(ArgumentError("refusing to clean a linked Pencil Maxwell rank file"))
        !ispath(path) || isfile(path) ||
            throw(ArgumentError("refusing to clean a non-file Pencil Maxwell rank entry"))
    end
    for name in files
        rm(joinpath(location, name); force = true)
    end
    # Nonrecursive removal preserves unknown entries and reports their presence.
    rm(location)
    return nothing
end

"""
    save_mpi_checkpoint(dir, model::DistributedMaxwellVacuum) -> manifest_path

Collectively checkpoint matching-layout periodic physical-pencil Maxwell state.
Each rank writes only six detached owned-cell arrays. A distinct schema-1 manifest
records configuration/phase, canonical rank files and their streamed SHA-256
digests; its serialized payload is checksummed before atomic publication. No
plan, communicator or scratch is serialized. Every rank file closes before the
manifest replaces the prior checkpoint. Pre-publication failure preserves the
old manifest; superseded cleanup failure reports the new checkpoint already usable.

The directory is shared and has exclusive save/load access. Generation symlinks
and noncanonical rank names reject; cleanup never removes unknown files or recurses.
Path resolution retains symlink/`..` semantics. MPI/fatal resource exclusions,
trusted Julia Serialization files and no crash-durability guarantee apply.
"""
function HybridPlasmaPIC.save_mpi_checkpoint(dir, model::DistributedMaxwellVacuum)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil Maxwell checkpoint"
    raw_dir, _ = _pencil_maxwell_checkpoint_agreement(model, dir, nranks)
    _distributed_maxwell_state_agreement(model, context)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "save preflight") do
        _distributed_maxwell_storage(model)
        _distributed_maxwell_finite_fields(model.fields)
    end
    state = HybridPlasmaPIC._mpi_collective_stage(comm, context, "snapshot preparation") do
        merge(
            _pencil_maxwell_checkpoint_config(model, nranks),
            (time = model.time[], step = model.step[], rank = rank),
            _pencil_maxwell_checkpoint_layout(model),
            (fields = (E = map(a -> copy(parent(a)), model.fields.E), B = map(a -> copy(parent(a)), model.fields.B)),),
        )
    end
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "directory setup") do
        rank == 0 && mkpath(raw_dir)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest_path = joinpath(directory, _PENCIL_MAXWELL_MANIFEST)
    previous, owned_generation = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation setup") do
        if rank == 0
            old = ispath(manifest_path) ? _pencil_maxwell_manifest(manifest_path, model, nranks) : nothing
            if old !== nothing
                location = _pencil_maxwell_generation_path(directory, old.generation)
                for r = 1:nranks
                    _pencil_maxwell_file_digest(joinpath(location, old.rank_files[r])) == old.rank_sha256[r] ||
                        throw(ArgumentError("superseded Pencil Maxwell checkpoint rank checksum mismatch"))
                end
            end
            old, basename(mktempdir(directory; prefix = _PENCIL_MAXWELL_GENERATION, cleanup = false))
        else
            nothing, ""
        end
    end
    try
        generation = _pencil_maxwell_broadcast_generation(owned_generation, comm, rank)
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank write") do
            location = _pencil_maxwell_generation_path(directory, generation)
            path = joinpath(location, HybridPlasmaPIC._mpi_checkpoint_rank_filename(rank + 1))
            !ispath(path) && !islink(path) ||
                throw(ArgumentError("new Pencil Maxwell rank file already exists"))
            HybridPlasmaPIC._atomic_serialize(path, state)
        end
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "manifest publication") do
            if rank == 0
                location = _pencil_maxwell_generation_path(directory, owned_generation)
                files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
                digests = map(name -> _pencil_maxwell_file_digest(joinpath(location, name)), files)
                manifest = merge(
                    _pencil_maxwell_checkpoint_config(model, nranks),
                    (time = state.time, step = state.step, generation = owned_generation, rank_files = files, rank_sha256 = digests),
                )
                HybridPlasmaPIC._atomic_serialize(manifest_path, _pencil_maxwell_envelope(manifest))
            end
        end
    catch original
        (original isa InterruptException || original isa OutOfMemoryError || original isa StackOverflowError) && rethrow()
        trace = catch_backtrace()
        try
            HybridPlasmaPIC._mpi_collective_stage(comm, context, "failed generation cleanup") do
                rank == 0 && _pencil_maxwell_remove_generation(directory, owned_generation, nranks)
            end
        catch cleanup_error
            (cleanup_error isa InterruptException || cleanup_error isa OutOfMemoryError || cleanup_error isa StackOverflowError) && rethrow()
            throw(CompositeException([CapturedException(original, trace), CapturedException(cleanup_error, catch_backtrace())]))
        end
        rethrow()
    end
    try
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "superseded generation cleanup") do
            rank == 0 && previous !== nothing &&
                _pencil_maxwell_remove_generation(directory, previous.generation, nranks)
        end
    catch cleanup_error
        (cleanup_error isa InterruptException || cleanup_error isa OutOfMemoryError || cleanup_error isa StackOverflowError) && rethrow()
        throw(CompositeException([
            ErrorException("Pencil Maxwell checkpoint was published, but superseded-file cleanup failed"),
            CapturedException(cleanup_error, catch_backtrace()),
        ]))
    end
    return manifest_path
end

function _pencil_maxwell_rank_state(path, digest)
    !islink(path) && isfile(path) ||
        throw(ArgumentError("Pencil Maxwell rank file must be an ordinary file"))
    return open(path, "r") do io
        bytes2hex(HybridPlasmaPIC.SHA.sha256(io)) == digest ||
            throw(ArgumentError("Pencil Maxwell rank checksum mismatch"))
        seekstart(io)
        HybridPlasmaPIC.Serialization.deserialize(io)
    end
end

function _pencil_maxwell_validate_rank_state(state, manifest, model::DistributedMaxwellVacuum{T}, rank, nranks) where {T}
    layout = _pencil_maxwell_checkpoint_layout(model)
    required = (keys(_pencil_maxwell_checkpoint_config(model, nranks))..., :time, :step, :rank, keys(layout)..., :fields)
    _pencil_maxwell_checkpoint_keys(state, required)
    _pencil_maxwell_checkpoint_common(state, model, nranks)
    state.rank isa Int && state.rank == rank || throw(ArgumentError("Pencil Maxwell checkpoint rank mismatch"))
    isequal(state.time, manifest.time) && state.step == manifest.step ||
        throw(ArgumentError("Pencil Maxwell rank phase differs from the manifest"))
    for name in keys(layout)
        value, expected = getproperty(state, name), getproperty(layout, name)
        typeof(value) === typeof(expected) && isequal(value, expected) ||
            throw(ArgumentError("Pencil Maxwell checkpoint $name differs from the current layout"))
    end
    _pencil_maxwell_checkpoint_keys(state.fields, (:E, :B))
    for components in (state.fields.E, state.fields.B)
        components isa NTuple{3,Array{T,3}} ||
            throw(ArgumentError("Pencil Maxwell checkpoint fields must be three-tuples of Array{$T,3}"))
        for array in components
            size(array) == layout.parent_size ||
                throw(DimensionMismatch("Pencil Maxwell checkpoint parent-array shape mismatch"))
        end
    end
    _distributed_maxwell_finite_fields(state.fields)
    incoming = (state.fields.E..., state.fields.B...)
    destination = (model.fields.E..., model.fields.B..., model._candidate.E..., model._candidate.B..., model._curl...)
    for i in eachindex(incoming)
        for j = 1:(i-1)
            Base.mightalias(incoming[i], incoming[j]) &&
                throw(ArgumentError("Pencil Maxwell incoming field arrays must be independent"))
        end
        for array in destination
            Base.mightalias(incoming[i], parent(array)) &&
                throw(ArgumentError("Pencil Maxwell incoming fields must be detached from destination storage"))
        end
    end
    return nothing
end

"""
    load_mpi_checkpoint!(model::DistributedMaxwellVacuum, dir) -> model

Collectively restore the distinct matching-rank/layout checkpoint. Checksummed
manifest/configuration agreement and each rank-file digest precede decoding the
detached fields. Required policies, phase, logical ranges, memory permutation,
ordinary array shape/type, finite values and input/destination independence all
validate on every rank before six plain copies and clock stores. Stale destination
field/time values are replaceable; destination storage and immutable configuration
must remain valid. Physical arrays and clock Ref identities persist.

No plan/scratch reconstruction, repartition, schema migration or cross-model load.
The shared-filesystem/exclusive-access and MPI/fatal-resource exclusions of save
apply. Runtime/restart/compiler/resource acceptance remains deferred.
"""
function HybridPlasmaPIC.load_mpi_checkpoint!(model::DistributedMaxwellVacuum, dir)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil Maxwell checkpoint"
    raw_dir, _ = _pencil_maxwell_checkpoint_agreement(model, dir, nranks)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "restore destination storage") do
        _distributed_maxwell_storage(model)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest = first(_agree_configuration(comm, context) do
        record = _pencil_maxwell_manifest(joinpath(directory, _PENCIL_MAXWELL_MANIFEST), model, nranks)
        _pencil_maxwell_generation_path(directory, record.generation)
        (record,)
    end)
    saved = HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank restore preflight") do
        location = _pencil_maxwell_generation_path(directory, manifest.generation)
        path = joinpath(location, manifest.rank_files[rank + 1])
        state = _pencil_maxwell_rank_state(path, manifest.rank_sha256[rank + 1])
        _pencil_maxwell_validate_rank_state(state, manifest, model, rank, nranks)
        state
    end
    for c = 1:3
        copyto!(parent(model.fields.E[c]), saved.fields.E[c])
        copyto!(parent(model.fields.B[c]), saved.fields.B[c])
    end
    model.time[] = saved.time
    model.step[] = saved.step
    return model
end
