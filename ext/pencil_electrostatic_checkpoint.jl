const _PENCIL_ES_MANIFEST = "pencil_electrostatic_manifest.ser"
const _PENCIL_ES_GENERATION = "_pencil_electrostatic_"
const _PENCIL_ES_CHECKPOINT_PHASE = (:initialized, :primed, :last_dt, :time, :step)

function _pencil_es_checkpoint_config(model::DistributedElectrostaticPIC, nranks)
    return merge(
        (
            schema = 1,
            model = :distributed_electrostatic_pic,
            integrator_policy = :shared_cic_leapfrog_v1,
            charge_policy = :exact_original_cic_route_round_once_v1,
            derivative_policy = :represented_odd_nyquist_zero_poisson_v1,
            layout_policy = :matching_physical_pencil_cic_anchor_v1,
        ),
        _pencil_es_bound_config(model),
        (nranks = nranks,),
    )
end

function _pencil_es_checkpoint_agreement(model::DistributedElectrostaticPIC{T}, dir, nranks) where {T}
    return _agree_configuration(model.g.fft.comm, "Pencil electrostatic checkpoint") do
        _pencil_es_mode(T)
        dir isa AbstractString || throw(ArgumentError("checkpoint directory must be a string"))
        !isempty(dir) || throw(ArgumentError("checkpoint directory must be nonempty"))
        String(dir), _pencil_es_checkpoint_config(model, nranks)
    end
end

function _pencil_es_checkpoint_keys(record, required)
    HybridPlasmaPIC._validate_checkpoint_container(record, required)
    length(keys(record)) == length(required) ||
        throw(ArgumentError("Pencil electrostatic checkpoint contains unsupported fields"))
    return nothing
end

function _pencil_es_checkpoint_common(record, model::DistributedElectrostaticPIC{T}, nranks) where {T}
    expected = _pencil_es_checkpoint_config(model, nranks)
    for name in keys(expected)
        value, reference = getproperty(record, name), getproperty(expected, name)
        typeof(value) === typeof(reference) && isequal(value, reference) ||
            throw(ArgumentError("Pencil electrostatic checkpoint $name does not match the model"))
    end
    record.initialized === true && record.primed isa Bool &&
        record.last_dt isa T && record.time isa T && record.step isa Int ||
        throw(ArgumentError("Pencil electrostatic checkpoint requires typed initialized phase and clock"))
    _pencil_pic_phase_values(record.initialized, record.primed, record.last_dt, record.time, record.step)
    return nothing
end

function _pencil_es_checkpoint_layout(model::DistributedElectrostaticPIC)
    field = model.rho
    ranges = PencilArrays.range_local(field, PencilArrays.LogicalOrder())
    return (
        logical_ranges = ntuple(d -> (Int(first(ranges[d])), Int(last(ranges[d]))), 3),
        permutation = PencilArrays.permutation(field) * (1, 2, 3),
        parent_size = size(parent(field)),
    )
end

_pencil_es_generation_name(name) =
    _pencil_checkpoint_generation_name(name, _PENCIL_ES_GENERATION, "Pencil electrostatic")
_pencil_es_generation_path(dir, name) =
    _pencil_checkpoint_generation_path(dir, name, _PENCIL_ES_GENERATION, "Pencil electrostatic")
_pencil_es_digest(value) = _pencil_checkpoint_digest(value, "Pencil electrostatic")
_pencil_es_file_digest(path) = _pencil_checkpoint_file_digest(path, "Pencil electrostatic")
_pencil_es_broadcast_generation(name, comm, rank) =
    _pencil_checkpoint_broadcast_generation(name, comm, rank, _PENCIL_ES_GENERATION, "Pencil electrostatic")
_pencil_es_remove_generation(dir, generation, nranks) =
    _pencil_checkpoint_remove_generation(dir, generation, nranks, _PENCIL_ES_GENERATION, "Pencil electrostatic")
_pencil_es_rank_state(path, digest) = _pencil_checkpoint_rank_state(path, digest, "Pencil electrostatic")

function _pencil_es_manifest(path, model, nranks)
    envelope = HybridPlasmaPIC.Serialization.deserialize(path)
    _pencil_es_checkpoint_keys(envelope, (:schema, :model, :payload, :sha256))
    envelope.schema === 1 && envelope.model === :pencil_electrostatic_manifest ||
        throw(ArgumentError("unsupported Pencil electrostatic manifest envelope"))
    envelope.payload isa Vector{UInt8} ||
        throw(ArgumentError("Pencil electrostatic manifest payload must be serialized bytes"))
    digest = _pencil_es_digest(envelope.sha256)
    bytes2hex(HybridPlasmaPIC.SHA.sha256(envelope.payload)) == digest ||
        throw(ArgumentError("Pencil electrostatic manifest checksum mismatch"))
    record = HybridPlasmaPIC._deserialize_mpi_payload(NamedTuple, envelope.payload)
    required = (
        keys(_pencil_es_checkpoint_config(model, nranks))...,
        _PENCIL_ES_CHECKPOINT_PHASE...,
        :generation,
        :rank_files,
        :rank_sha256,
    )
    _pencil_es_checkpoint_keys(record, required)
    _pencil_es_checkpoint_common(record, model, nranks)
    record.generation isa String || throw(ArgumentError("checkpoint generation must be a String"))
    _pencil_es_generation_name(record.generation)
    expected_files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    record.rank_files isa Tuple && isequal(record.rank_files, expected_files) ||
        throw(ArgumentError("Pencil electrostatic manifest requires canonical ordered rank filenames"))
    record.rank_sha256 isa Tuple && length(record.rank_sha256) == nranks ||
        throw(ArgumentError("Pencil electrostatic manifest requires one digest per rank"))
    foreach(_pencil_es_digest, record.rank_sha256)
    return record
end

_pencil_es_envelope(manifest) =
    _pencil_pic_prepared_envelope(manifest, :pencil_electrostatic_manifest, "Pencil electrostatic")

function _pencil_pic_prepared_envelope(manifest, tag, label)
    payload = HybridPlasmaPIC._serialize_mpi_payload(manifest)
    HybridPlasmaPIC._check_mpi_count(length(payload), "$label manifest payload")
    return (
        schema = 1,
        model = tag,
        payload = payload,
        sha256 = bytes2hex(HybridPlasmaPIC.SHA.sha256(payload)),
    )
end

function _pencil_es_checkpoint_snapshot(model, rank, nranks)
    merge(
        _pencil_es_checkpoint_config(model, nranks),
        (
            initialized = model.initialized[], primed = model.primed[], last_dt = model.last_dt[],
            time = model.time[], step = model.step[], rank = rank,
        ),
        _pencil_es_checkpoint_layout(model),
        (
            fields = (E = map(a -> copy(parent(a)), model.E), rho = copy(parent(model.rho))),
            species = deepcopy([HybridPlasmaPIC._checkpoint_particle_state(ps) for ps in model.owners]),
        ),
    )
end

"""
    save_mpi_checkpoint(dir, model::DistributedElectrostaticPIC) -> manifest_path

Collectively checkpoint an initialized matching-layout periodic ES collection.
Each rank detaches its owned E/rho and ordered local particle leaves/scalars.
Configuration, actual ownership and common leapfrog/clock phase are explicit;
derived charge ledgers follow from the saved physical state, not extra scratch.
Plans, communicators, lookup/ghost/staged storage are omitted. Each rank file
closes before the checksummed manifest atomically publishes a unique generation.

An existing checkpoint directory must match the current cached bound q/m and
configuration. After a scalar-rebinding refresh, use another directory. Earlier
generation/configuration/read/write failures preserve the prior manifest;
superseded cleanup failure reports the new checkpoint already usable. Shared
filesystem, exclusive access, trusted Serialization, MPI/fatal-resource and
no-crash-durability limits match the physical-pencil Maxwell checkpoint.
"""
function HybridPlasmaPIC.save_mpi_checkpoint(dir, model::DistributedElectrostaticPIC)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil electrostatic checkpoint"
    raw_dir, common = _pencil_es_checkpoint_agreement(model, dir, nranks)
    _pencil_es_agreement(model, model.owners, context)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "save particle ownership") do
        _pencil_owned_particles(model.owners, model.g, model._layout)
    end
    _pencil_pic_charges(model, model.owners, context)
    state, phase = HybridPlasmaPIC._mpi_collective_stage(comm, context, "snapshot preparation") do
        state = _pencil_es_checkpoint_snapshot(model, rank, nranks)
        phase = (initialized = state.initialized, primed = state.primed, last_dt = state.last_dt,
                 time = state.time, step = state.step)
        state, phase
    end
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "directory setup") do
        rank == 0 && mkpath(raw_dir)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest_path = joinpath(directory, _PENCIL_ES_MANIFEST)
    previous, owned_generation = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation setup") do
        if rank == 0
            old = ispath(manifest_path) ? _pencil_es_manifest(manifest_path, model, nranks) : nothing
            if old !== nothing
                location = _pencil_es_generation_path(directory, old.generation)
                for r = 1:nranks
                    _pencil_es_file_digest(joinpath(location, old.rank_files[r])) == old.rank_sha256[r] ||
                        throw(ArgumentError("superseded Pencil electrostatic checkpoint rank checksum mismatch"))
                end
            end
            old, basename(mktempdir(directory; prefix = _PENCIL_ES_GENERATION, cleanup = false))
        else
            nothing, ""
        end
    end
    return _pencil_pic_write_prepared(directory, manifest_path, previous, owned_generation,
        state, common, phase,
        (prefix = _PENCIL_ES_GENERATION, tag = :pencil_electrostatic_manifest, label = "Pencil electrostatic"),
        comm, rank, nranks, context)
end

# Prepared ES/composition generations use explicit caller-owned schemas/policies.
# Bare EM persistence retains its separate established orchestration.
function _pencil_pic_write_prepared(directory, manifest_path, previous, owned_generation,
    state, common, phase, policy, comm, rank, nranks, context)
    try
        generation = _pencil_checkpoint_broadcast_generation(owned_generation, comm, rank, policy.prefix, policy.label)
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank write") do
            location = _pencil_checkpoint_generation_path(directory, generation, policy.prefix, policy.label)
            path = joinpath(location, HybridPlasmaPIC._mpi_checkpoint_rank_filename(rank + 1))
            !ispath(path) && !islink(path) ||
                throw(ArgumentError("new $(policy.label) rank file already exists"))
            HybridPlasmaPIC._atomic_serialize(path, state)
        end
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "manifest publication") do
            if rank == 0
                location = _pencil_checkpoint_generation_path(directory, owned_generation, policy.prefix, policy.label)
                files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
                digests = map(name -> _pencil_checkpoint_file_digest(joinpath(location, name), policy.label), files)
                manifest = merge(
                    common,
                    (;
                        phase..., generation = owned_generation,
                        rank_files = files, rank_sha256 = digests,
                    ),
                )
                HybridPlasmaPIC._atomic_serialize(manifest_path, _pencil_pic_prepared_envelope(manifest, policy.tag, policy.label))
            end
        end
    catch original
        (original isa InterruptException || original isa OutOfMemoryError || original isa StackOverflowError) && rethrow()
        trace = catch_backtrace()
        try
            HybridPlasmaPIC._mpi_collective_stage(comm, context, "failed generation cleanup") do
                rank == 0 && _pencil_checkpoint_remove_generation(directory, owned_generation, nranks, policy.prefix, policy.label)
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
                _pencil_checkpoint_remove_generation(directory, previous.generation, nranks, policy.prefix, policy.label)
        end
    catch cleanup_error
        (cleanup_error isa InterruptException || cleanup_error isa OutOfMemoryError || cleanup_error isa StackOverflowError) && rethrow()
        throw(CompositeException([
            ErrorException("$(policy.label) checkpoint was published, but superseded-file cleanup failed"),
            CapturedException(cleanup_error, catch_backtrace()),
        ]))
    end
    return manifest_path
end

function _pencil_es_prepare_rank_state(state, manifest, model::DistributedElectrostaticPIC{T,G,R,P}, rank, nranks) where {T,G,R,P}
    layout = _pencil_es_checkpoint_layout(model)
    required = (
        keys(_pencil_es_checkpoint_config(model, nranks))..., _PENCIL_ES_CHECKPOINT_PHASE...,
        :rank, keys(layout)..., :fields, :species,
    )
    _pencil_es_checkpoint_keys(state, required)
    _pencil_es_checkpoint_common(state, model, nranks)
    state.rank isa Int && state.rank == rank || throw(ArgumentError("Pencil electrostatic checkpoint rank mismatch"))
    for name in _PENCIL_ES_CHECKPOINT_PHASE
        isequal(getproperty(state, name), getproperty(manifest, name)) ||
            throw(ArgumentError("Pencil electrostatic rank phase differs from the manifest"))
    end
    for name in keys(layout)
        value, expected = getproperty(state, name), getproperty(layout, name)
        typeof(value) === typeof(expected) && isequal(value, expected) ||
            throw(ArgumentError("Pencil electrostatic checkpoint $name differs from the current layout"))
    end
    _pencil_es_checkpoint_keys(state.fields, (:E, :rho))
    state.fields.E isa NTuple{3,Array{T,3}} && state.fields.rho isa Array{T,3} ||
        throw(ArgumentError("Pencil electrostatic checkpoint requires ordinary Array E/rho fields"))
    fields = (state.fields.E..., state.fields.rho)
    destination = (model.E..., model.rho, model._candidate_E..., model._candidate_rho)
    for (i, array) in enumerate(fields)
        size(array) == layout.parent_size ||
            throw(DimensionMismatch("Pencil electrostatic checkpoint parent-array shape mismatch"))
        all(isfinite, array) || throw(ArgumentError("Pencil electrostatic checkpoint fields must be finite"))
        for j = 1:(i-1)
            HybridPlasmaPIC._particle_arrays_mightalias(array, fields[j]) &&
                throw(ArgumentError("Pencil electrostatic incoming fields must be independent"))
        end
        for other in destination
            HybridPlasmaPIC._particle_arrays_mightalias(array, parent(other)) &&
                throw(ArgumentError("Pencil electrostatic incoming fields must be detached from destination storage"))
        end
        _aliases_workspace(array, model.g) &&
            throw(ArgumentError("Pencil electrostatic incoming fields must be detached from Fourier storage"))
    end
    state.species isa Vector && length(state.species) == length(model.owners) ||
        throw(ArgumentError("Pencil electrostatic checkpoint requires matching ordered population records"))
    populations = P[]
    for (s, record) in enumerate(state.species)
        _pencil_es_checkpoint_keys(record, HybridPlasmaPIC._CHECKPOINT_PARTICLE_FIELDS)
        ps = HybridPlasmaPIC._checkpoint_species(record, Val(3), T)
        ps isa P || throw(ArgumentError("Pencil electrostatic checkpoint requires ordinary Vector particle leaves"))
        isequal((ps.q, ps.m), model._species_scalars[s]) ||
            throw(ArgumentError("Pencil electrostatic checkpoint particle q/m differs from the cached model binding"))
        push!(populations, ps)
    end
    _pencil_es_species(model, populations)
    _pencil_owned_particles(populations, model.g, model._layout)
    for ps in populations, live in model.owners
        HybridPlasmaPIC._particle_sets_mightalias(ps, live) &&
            throw(ArgumentError("Pencil electrostatic incoming particles must be detached from destination storage"))
    end
    for array in fields, species in (populations, model.owners, model._staged), ps in species
        HybridPlasmaPIC._mightalias_particle_storage(array, ps) &&
            throw(ArgumentError("Pencil electrostatic incoming fields must not alias particle storage"))
    end
    # File decoding detaches fields; complete same-type indexed copies prepare
    # every replacement leaf even when local counts change or a population empties.
    prepared = P[HybridPlasmaPIC._prepare_particle_replacement(ps, eachindex(ps.weight)) for ps in populations]
    return (state = state, particles = prepared)
end

"""
    load_mpi_checkpoint!(model::DistributedElectrostaticPIC, dir) -> model

Collectively restore matching rank count, actual pencil ownership/permutation and
cached model policy/q/m. Destination physical values and all clock/phase values
may be stale; valid ordinary storage/configuration is required. Source checksum,
typed phase, particles, fields, ownership and original physical neutrality pass
on all ranks before publication. Local counts may change through fully prepared
same-type leaves; ParticleSet owners, field arrays and phase Ref identities persist.
E/rho restore verbatim, without migration or a Poisson solve. Derived lookup/ghost
and particle/FFT scratch is omitted and reused by later normal staging.

No current-phase agreement, public resizing or late validation occurs on restore.
Shared-filesystem/exclusive-access, trusted-Serialization and MPI/fatal-resource
exclusions match save. Runtime/restart/liveness/resource acceptance is deferred.
"""
function HybridPlasmaPIC.load_mpi_checkpoint!(model::DistributedElectrostaticPIC, dir)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil electrostatic checkpoint"
    raw_dir, _ = _pencil_es_checkpoint_agreement(model, dir, nranks)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "restore destination storage") do
        _pencil_es_storage(model)
        _pencil_es_species(model, model.owners; bound = true, values = false, factors = false)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest = first(_agree_configuration(comm, context) do
        record = _pencil_es_manifest(joinpath(directory, _PENCIL_ES_MANIFEST), model, nranks)
        _pencil_es_generation_path(directory, record.generation)
        (record,)
    end)
    saved = HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank restore preparation") do
        location = _pencil_es_generation_path(directory, manifest.generation)
        path = joinpath(location, manifest.rank_files[rank + 1])
        state = _pencil_es_rank_state(path, manifest.rank_sha256[rank + 1])
        _pencil_es_prepare_rank_state(state, manifest, model, rank, nranks)
    end
    # This collective checks source values only. Destination clocks/phase and
    # particle values are still replaceable, and no field reconstruction occurs.
    _pencil_pic_charges(model, saved.particles, context)
    return _pencil_es_restore_prepared!(model, saved)
end

function _pencil_es_restore_prepared!(model, saved)
    for (ps, replacement) in zip(model.owners, saved.particles)
        HybridPlasmaPIC._publish_particle_replacement!(ps, replacement)
        ps.q = replacement.q
        ps.m = replacement.m
    end
    for c = 1:3
        copyto!(parent(model.E[c]), saved.state.fields.E[c])
    end
    copyto!(parent(model.rho), saved.state.fields.rho)
    model.initialized[] = saved.state.initialized
    model.primed[] = saved.state.primed
    model.last_dt[] = saved.state.last_dt
    model.time[] = saved.state.time
    model.step[] = saved.state.step
    return model
end
