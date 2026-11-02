const _PENCIL_EMPIC_MANIFEST = "pencil_empic_manifest.ser"
const _PENCIL_EMPIC_GENERATION = "_pencil_empic_"
const _PENCIL_EMPIC_CHECKPOINT_PHASE = (:initialized, :primed, :last_dt, :time, :step)

function _pencil_empic_checkpoint_config(model::DistributedEMPIC, nranks)
    return merge(
        (
            schema = 1,
            model = :distributed_empic,
            integrator_policy = :canonical_boris_carried_b_single_rate_v1,
            current_policy = :midpoint_cic_continuity_dc_raw_other_null_zero_v1,
            charge_policy = :exact_original_cic_route_round_once_v1,
            derivative_policy = :represented_odd_nyquist_zero_v1,
            layout_policy = :matching_physical_pencil_cic_anchor_v1,
            field_phase = :integer_e_x_initial_full_then_carried_b_v_v1,
            shape = :cic,
            mobility = :all_mobile,
            relativistic = false,
            n_sub = 1,
        ),
        _pencil_empic_bound_config(model),
        (nranks = nranks,),
    )
end

function _pencil_empic_checkpoint_agreement(model::DistributedEMPIC{T}, dir, nranks) where {T}
    return _agree_configuration(model.g.fft.comm, "Pencil electromagnetic checkpoint") do
        _pencil_empic_mode(T)
        dir isa AbstractString || throw(ArgumentError("checkpoint directory must be a string"))
        !isempty(dir) || throw(ArgumentError("checkpoint directory must be nonempty"))
        String(dir), _pencil_empic_checkpoint_config(model, nranks)
    end
end

function _pencil_empic_checkpoint_keys(record, required)
    HybridPlasmaPIC._validate_checkpoint_container(record, required)
    length(keys(record)) == length(required) ||
        throw(ArgumentError("Pencil electromagnetic checkpoint contains unsupported fields"))
    return nothing
end

function _pencil_empic_checkpoint_common(record, model::DistributedEMPIC{T}, nranks) where {T}
    expected = _pencil_empic_checkpoint_config(model, nranks)
    for name in keys(expected)
        value, reference = getproperty(record, name), getproperty(expected, name)
        typeof(value) === typeof(reference) && isequal(value, reference) ||
            throw(ArgumentError("Pencil electromagnetic checkpoint $name does not match the model"))
    end
    record.initialized === true && record.primed isa Bool &&
        record.last_dt isa T && record.time isa T && record.step isa Int ||
        throw(ArgumentError("Pencil electromagnetic checkpoint requires typed initialized phase and clock"))
    _pencil_empic_phase_values(record.initialized, record.primed, record.last_dt, record.time, record.step, model._cfl)
    return nothing
end

function _pencil_empic_checkpoint_layout(model::DistributedEMPIC)
    field = model.rho_np1
    ranges = PencilArrays.range_local(field, PencilArrays.LogicalOrder())
    return (
        logical_ranges = ntuple(d -> (Int(first(ranges[d])), Int(last(ranges[d]))), 3),
        permutation = PencilArrays.permutation(field) * (1, 2, 3),
        parent_size = size(parent(field)),
    )
end

_pencil_empic_generation_name(name) =
    _pencil_checkpoint_generation_name(name, _PENCIL_EMPIC_GENERATION, "Pencil electromagnetic")
_pencil_empic_generation_path(dir, name) =
    _pencil_checkpoint_generation_path(dir, name, _PENCIL_EMPIC_GENERATION, "Pencil electromagnetic")
_pencil_empic_digest(value) = _pencil_checkpoint_digest(value, "Pencil electromagnetic")
_pencil_empic_file_digest(path) = _pencil_checkpoint_file_digest(path, "Pencil electromagnetic")
_pencil_empic_broadcast_generation(name, comm, rank) =
    _pencil_checkpoint_broadcast_generation(name, comm, rank, _PENCIL_EMPIC_GENERATION, "Pencil electromagnetic")
_pencil_empic_remove_generation(dir, generation, nranks) =
    _pencil_checkpoint_remove_generation(dir, generation, nranks, _PENCIL_EMPIC_GENERATION, "Pencil electromagnetic")
_pencil_empic_rank_state(path, digest) = _pencil_checkpoint_rank_state(path, digest, "Pencil electromagnetic")

function _pencil_empic_manifest(path, model, nranks)
    envelope = HybridPlasmaPIC.Serialization.deserialize(path)
    _pencil_empic_checkpoint_keys(envelope, (:schema, :model, :payload, :sha256))
    envelope.schema === 1 && envelope.model === :pencil_empic_manifest ||
        throw(ArgumentError("unsupported Pencil electromagnetic manifest envelope"))
    envelope.payload isa Vector{UInt8} ||
        throw(ArgumentError("Pencil electromagnetic manifest payload must be serialized bytes"))
    digest = _pencil_empic_digest(envelope.sha256)
    bytes2hex(HybridPlasmaPIC.SHA.sha256(envelope.payload)) == digest ||
        throw(ArgumentError("Pencil electromagnetic manifest checksum mismatch"))
    record = HybridPlasmaPIC._deserialize_mpi_payload(NamedTuple, envelope.payload)
    required = (
        keys(_pencil_empic_checkpoint_config(model, nranks))...,
        _PENCIL_EMPIC_CHECKPOINT_PHASE...,
        :generation,
        :rank_files,
        :rank_sha256,
    )
    _pencil_empic_checkpoint_keys(record, required)
    _pencil_empic_checkpoint_common(record, model, nranks)
    record.generation isa String || throw(ArgumentError("checkpoint generation must be a String"))
    _pencil_empic_generation_name(record.generation)
    expected_files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    record.rank_files isa Tuple && isequal(record.rank_files, expected_files) ||
        throw(ArgumentError("Pencil electromagnetic manifest requires canonical ordered rank filenames"))
    record.rank_sha256 isa Tuple && length(record.rank_sha256) == nranks ||
        throw(ArgumentError("Pencil electromagnetic manifest requires one digest per rank"))
    foreach(_pencil_empic_digest, record.rank_sha256)
    return record
end

function _pencil_empic_envelope(manifest)
    payload = HybridPlasmaPIC._serialize_mpi_payload(manifest)
    HybridPlasmaPIC._check_mpi_count(length(payload), "Pencil electromagnetic manifest payload")
    return (
        schema = 1,
        model = :pencil_empic_manifest,
        payload = payload,
        sha256 = bytes2hex(HybridPlasmaPIC.SHA.sha256(payload)),
    )
end

function _pencil_empic_checkpoint_snapshot(model, rank, nranks)
    merge(
        _pencil_empic_checkpoint_config(model, nranks),
        (
            initialized = model.initialized[], primed = model.primed[], last_dt = model.last_dt[],
            time = model.time[], step = model.step[], rank = rank,
        ),
        _pencil_empic_checkpoint_layout(model),
        (
            fields = (
                E = map(a -> copy(parent(a)), model.E),
                B = map(a -> copy(parent(a)), model.B),
                J = map(a -> copy(parent(a)), model.J),
                rho_n = copy(parent(model.rho_n)), rho_np1 = copy(parent(model.rho_np1)),
            ),
            species = deepcopy([HybridPlasmaPIC._checkpoint_particle_state(ps) for ps in model.owners]),
        ),
    )
end

"""
    save_mpi_checkpoint(dir, model::DistributedEMPIC) -> manifest_path

Collectively save an initialized matching-layout single-rate EM PIC boundary.
Each rank detaches E/B/J/rho_n/rho_np1 and complete ordered local particles.
Exact cached configuration, actual layout and sole magnetic/velocity timestep
history persist. Initial current is stored zero with identical density endpoints;
advanced J retains the last corrected midpoint current. No field reconstruction,
count-cache array, plans, communicator or particle/FFT scratch is serialized.

An existing directory must match current cached q/m/c and policy. All canonical
rank files close before the checksummed manifest atomically publishes a unique
generation. Earlier failures preserve its previous manifest; superseded cleanup
failure reports the new checkpoint usable. Shared/exclusive-access, trusted
Serialization, MPI/fatal-resource and no-crash-durability limits match ES save.
"""
function HybridPlasmaPIC.save_mpi_checkpoint(dir, model::DistributedEMPIC)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil electromagnetic checkpoint"
    raw_dir, _ = _pencil_empic_checkpoint_agreement(model, dir, nranks)
    _pencil_empic_agreement(model, model.owners, context)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "save particle ownership") do
        _pencil_owned_particles(model.owners, model.g, model._layout)
    end
    _pencil_pic_charges(model, model.owners, context)
    state = HybridPlasmaPIC._mpi_collective_stage(comm, context, "snapshot preparation") do
        _pencil_empic_checkpoint_snapshot(model, rank, nranks)
    end
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "directory setup") do
        rank == 0 && mkpath(raw_dir)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest_path = joinpath(directory, _PENCIL_EMPIC_MANIFEST)
    previous, owned_generation = HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation setup") do
        if rank == 0
            old = ispath(manifest_path) ? _pencil_empic_manifest(manifest_path, model, nranks) : nothing
            if old !== nothing
                location = _pencil_empic_generation_path(directory, old.generation)
                for r = 1:nranks
                    _pencil_empic_file_digest(joinpath(location, old.rank_files[r])) == old.rank_sha256[r] ||
                        throw(ArgumentError("superseded Pencil electromagnetic checkpoint rank checksum mismatch"))
                end
            end
            old, basename(mktempdir(directory; prefix = _PENCIL_EMPIC_GENERATION, cleanup = false))
        else
            nothing, ""
        end
    end
    try
        generation = _pencil_empic_broadcast_generation(owned_generation, comm, rank)
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank write") do
            location = _pencil_empic_generation_path(directory, generation)
            path = joinpath(location, HybridPlasmaPIC._mpi_checkpoint_rank_filename(rank + 1))
            !ispath(path) && !islink(path) ||
                throw(ArgumentError("new Pencil electromagnetic rank file already exists"))
            HybridPlasmaPIC._atomic_serialize(path, state)
        end
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "manifest publication") do
            if rank == 0
                location = _pencil_empic_generation_path(directory, owned_generation)
                files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
                digests = map(name -> _pencil_empic_file_digest(joinpath(location, name)), files)
                manifest = merge(
                    _pencil_empic_checkpoint_config(model, nranks),
                    (
                        initialized = state.initialized, primed = state.primed, last_dt = state.last_dt,
                        time = state.time, step = state.step, generation = owned_generation,
                        rank_files = files, rank_sha256 = digests,
                    ),
                )
                HybridPlasmaPIC._atomic_serialize(manifest_path, _pencil_empic_envelope(manifest))
            end
        end
    catch original
        (original isa InterruptException || original isa OutOfMemoryError || original isa StackOverflowError) && rethrow()
        trace = catch_backtrace()
        try
            HybridPlasmaPIC._mpi_collective_stage(comm, context, "failed generation cleanup") do
                rank == 0 && _pencil_empic_remove_generation(directory, owned_generation, nranks)
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
                _pencil_empic_remove_generation(directory, previous.generation, nranks)
        end
    catch cleanup_error
        (cleanup_error isa InterruptException || cleanup_error isa OutOfMemoryError || cleanup_error isa StackOverflowError) && rethrow()
        throw(CompositeException([
            ErrorException("Pencil electromagnetic checkpoint was published, but superseded-file cleanup failed"),
            CapturedException(cleanup_error, catch_backtrace()),
        ]))
    end
    return manifest_path
end

function _pencil_empic_prepare_rank_state(state, manifest, model::DistributedEMPIC{T,G,R,P}, rank, nranks) where {T,G,R,P}
    layout = _pencil_empic_checkpoint_layout(model)
    required = (
        keys(_pencil_empic_checkpoint_config(model, nranks))..., _PENCIL_EMPIC_CHECKPOINT_PHASE...,
        :rank, keys(layout)..., :fields, :species,
    )
    _pencil_empic_checkpoint_keys(state, required)
    _pencil_empic_checkpoint_common(state, model, nranks)
    state.rank isa Int && state.rank == rank || throw(ArgumentError("Pencil electromagnetic checkpoint rank mismatch"))
    for name in _PENCIL_EMPIC_CHECKPOINT_PHASE
        isequal(getproperty(state, name), getproperty(manifest, name)) ||
            throw(ArgumentError("Pencil electromagnetic rank phase differs from the manifest"))
    end
    for name in keys(layout)
        value, expected = getproperty(state, name), getproperty(layout, name)
        typeof(value) === typeof(expected) && isequal(value, expected) ||
            throw(ArgumentError("Pencil electromagnetic checkpoint $name differs from the current layout"))
    end
    _pencil_empic_checkpoint_keys(state.fields, (:E, :B, :J, :rho_n, :rho_np1))
    state.fields.E isa NTuple{3,Array{T,3}} && state.fields.B isa NTuple{3,Array{T,3}} &&
        state.fields.J isa NTuple{3,Array{T,3}} &&
        state.fields.rho_n isa Array{T,3} && state.fields.rho_np1 isa Array{T,3} ||
        throw(ArgumentError("Pencil electromagnetic checkpoint requires eleven ordinary Array fields"))
    fields = _pencil_empic_field_arrays(state.fields)
    destination = (
        _pencil_empic_arrays(model)..., model._source_spectral,
        _pencil_empic_metadata(model)..., _pencil_empic_midpoints(model)...,
    )
    for (i, array) in enumerate(fields)
        size(array) == layout.parent_size ||
            throw(DimensionMismatch("Pencil electromagnetic checkpoint parent-array shape mismatch"))
        all(isfinite, array) || throw(ArgumentError("Pencil electromagnetic checkpoint fields must be finite"))
        for j = 1:(i-1)
            HybridPlasmaPIC._particle_arrays_mightalias(array, fields[j]) &&
                throw(ArgumentError("Pencil electromagnetic incoming fields must be independent"))
        end
        for other in destination
            HybridPlasmaPIC._particle_arrays_mightalias(array, parent(other)) &&
                throw(ArgumentError("Pencil electromagnetic incoming fields must be detached from destination storage"))
        end
        _aliases_workspace(array, model.g) &&
            throw(ArgumentError("Pencil electromagnetic incoming fields must be detached from Fourier storage"))
    end
    _pencil_empic_phase_fields(state.fields.J, state.fields.rho_n, state.fields.rho_np1, state.step)
    state.species isa Vector && length(state.species) == length(model.owners) ||
        throw(ArgumentError("Pencil electromagnetic checkpoint requires matching ordered population records"))
    populations = P[]
    for (s, record) in enumerate(state.species)
        _pencil_empic_checkpoint_keys(record, HybridPlasmaPIC._CHECKPOINT_PARTICLE_FIELDS)
        ps = HybridPlasmaPIC._checkpoint_species(record, Val(3), T)
        ps isa P || throw(ArgumentError("Pencil electromagnetic checkpoint requires ordinary Vector particle leaves"))
        isequal((ps.q, ps.m), model._species_scalars[s]) ||
            throw(ArgumentError("Pencil electromagnetic checkpoint particle q/m differs from the cached model binding"))
        push!(populations, ps)
    end
    _pencil_empic_species(model, populations; bound = false, values = true, counts = false)
    _pencil_owned_particles(populations, model.g, model._layout)
    for ps in populations, live in model.owners
        HybridPlasmaPIC._particle_sets_mightalias(ps, live) &&
            throw(ArgumentError("Pencil electromagnetic incoming particles must be detached from destination storage"))
    end
    for array in fields, species in (populations, model.owners, model._staged), ps in species
        HybridPlasmaPIC._mightalias_particle_storage(array, ps) &&
            throw(ArgumentError("Pencil electromagnetic incoming fields must not alias particle storage"))
    end
    # File decoding detaches fields; complete same-type indexed copies prepare
    # every replacement leaf even when local counts change or a population empties.
    prepared = P[HybridPlasmaPIC._prepare_particle_replacement(ps, eachindex(ps.weight)) for ps in populations]
    scalars = Tuple{T,T}[(ps.q, ps.m) for ps in prepared]
    counts = Int[length(ps.weight) for ps in prepared]
    # Storage preflight fixed both cache lengths to this ordered species count.
    return (state = state, particles = prepared, scalars = scalars, counts = counts)
end

"""
    load_mpi_checkpoint!(model::DistributedEMPIC, dir) -> model

Collectively restore exact cached configuration and actual rank/layout match.
Old physical values, q/m, phase/clocks, local counts and count-cache entries may
be stale; representation and immutable policy must remain valid. Source checksums,
phase/current history, finite independent fields/particles and ownership validate,
and every same-P leaf plus derived count/scalar row prepares on all ranks before
source physical neutrality and publication. Counts may change, including zero.

Publish saved q/m/leaves and all eleven physical fields verbatim, then prepared
cache rows and typed phase/clock Refs. Owner/field/Ref/cache Vector identities
persist. No ready destination phase/count check, public resizing, priming, migration,
deposition, projection or current regeneration occurs. Derived scratch is omitted.
The shared-filesystem/exclusive-access and trusted-Serialization/MPI/fatal-resource
premises match save; runtime/restart/liveness/resource acceptance is deferred.
"""
function HybridPlasmaPIC.load_mpi_checkpoint!(model::DistributedEMPIC, dir)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil electromagnetic checkpoint"
    raw_dir, _ = _pencil_empic_checkpoint_agreement(model, dir, nranks)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "restore destination storage") do
        _pencil_empic_storage(model)
        _pencil_empic_species(model, model.owners; bound = true, values = false, factors = false, counts = false)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest = first(_agree_configuration(comm, context) do
        record = _pencil_empic_manifest(joinpath(directory, _PENCIL_EMPIC_MANIFEST), model, nranks)
        _pencil_empic_generation_path(directory, record.generation)
        (record,)
    end)
    saved = HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank restore preparation") do
        location = _pencil_empic_generation_path(directory, manifest.generation)
        path = joinpath(location, manifest.rank_files[rank + 1])
        state = _pencil_empic_rank_state(path, manifest.rank_sha256[rank + 1])
        _pencil_empic_prepare_rank_state(state, manifest, model, rank, nranks)
    end
    # This collective checks source values only. Destination clocks/phase and
    # particle values are still replaceable, and no field reconstruction occurs.
    _pencil_pic_charges(model, saved.particles, context)
    return _pencil_empic_restore_prepared!(model, saved)
end

function _pencil_empic_restore_prepared!(model, saved)
    for (ps, replacement) in zip(model.owners, saved.particles)
        HybridPlasmaPIC._publish_particle_replacement!(ps, replacement)
        ps.q = replacement.q
        ps.m = replacement.m
    end
    for c = 1:3
        copyto!(parent(model.E[c]), saved.state.fields.E[c])
        copyto!(parent(model.B[c]), saved.state.fields.B[c])
        copyto!(parent(model.J[c]), saved.state.fields.J[c])
    end
    copyto!(parent(model.rho_n), saved.state.fields.rho_n)
    copyto!(parent(model.rho_np1), saved.state.fields.rho_np1)
    copyto!(model._species_scalars, saved.scalars)
    copyto!(model._population_counts, saved.counts)
    model.initialized[] = saved.state.initialized
    model.primed[] = saved.state.primed
    model.last_dt[] = saved.state.last_dt
    model.time[] = saved.state.time
    model.step[] = saved.state.step
    return model
end
