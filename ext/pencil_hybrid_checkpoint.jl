# Complete carried isothermal hybrid state; field and phase reconstruction is absent.
const _PENCIL_HYBRID_MANIFEST = "pencil_hybrid_manifest.ser"
const _PENCIL_HYBRID_GENERATION = "_pencil_hybrid_"
const _PENCIL_HYBRID_CHECKPOINT_PHASE = (:initialized, :primed, :last_dt, :time, :step)
const _PENCIL_HYBRID_FIELD_NAMES = (:B, :E, :n, :ui, :pe, :J, :gradp, :ninv)

function _pencil_hybrid_checkpoint_config(model::DistributedHybrid, nranks)
    return merge(
        (
            schema = 1,
            model = :distributed_isothermal_hybrid,
            integrator_policy = :canonical_boris_frozen_midpoint_rk4_endpoint_predictor_v1,
            moment_policy = :exact_original_charge_weighted_cic_route_round_once_v1,
            derivative_policy = :represented_odd_nyquist_zero_v1,
            faraday_policy = :differentiated_ohm_source_risk_rejected_v1,
            layout_policy = :matching_physical_pencil_cic_anchor_v1,
            field_phase = :integer_b_x_predictor_moments_carried_e_half_v_v1,
            shape = :cic,
            mobility = :all_mobile,
            relativistic = false,
        ),
        _pencil_hybrid_bound_config(model),
        (nranks = nranks,),
    )
end

function _pencil_hybrid_checkpoint_agreement(model::DistributedHybrid{T}, dir, nranks) where {T}
    return _agree_configuration(model.g.fft.comm, "Pencil hybrid checkpoint") do
        _pencil_hybrid_mode(T)
        dir isa AbstractString && !isempty(dir) ||
            throw(ArgumentError("checkpoint directory must be a nonempty string"))
        String(dir), _pencil_hybrid_checkpoint_config(model, nranks)
    end
end

function _pencil_hybrid_checkpoint_keys(record, required)
    HybridPlasmaPIC._validate_checkpoint_container(record, required)
    length(keys(record)) == length(required) ||
        throw(ArgumentError("Pencil hybrid checkpoint contains unsupported fields"))
    return nothing
end

function _pencil_hybrid_checkpoint_common(record, model::DistributedHybrid{T}, nranks) where {T}
    expected = _pencil_hybrid_checkpoint_config(model, nranks)
    for name in keys(expected)
        value, reference = getproperty(record, name), getproperty(expected, name)
        typeof(value) === typeof(reference) && isequal(value, reference) ||
            throw(ArgumentError("Pencil hybrid checkpoint $name does not match the model"))
    end
    record.initialized === true &&
        record.primed isa Bool &&
        record.last_dt isa T &&
        record.time isa T &&
        record.step isa Int ||
        throw(ArgumentError("Pencil hybrid checkpoint requires typed initialized phase and clock"))
    _pencil_hybrid_phase_values(
        record.initialized,
        record.primed,
        record.last_dt,
        record.time,
        record.step,
    )
    return nothing
end

function _pencil_hybrid_checkpoint_layout(model::DistributedHybrid)
    field = model.fields.n
    ranges = PencilArrays.range_local(field, PencilArrays.LogicalOrder())
    return (
        logical_ranges = ntuple(d -> (Int(first(ranges[d])), Int(last(ranges[d]))), 3),
        permutation = PencilArrays.permutation(field) * (1, 2, 3),
        parent_size = size(parent(field)),
    )
end

_pencil_hybrid_generation_name(name) =
    _pencil_checkpoint_generation_name(name, _PENCIL_HYBRID_GENERATION, "Pencil hybrid")
_pencil_hybrid_generation_path(dir, name) =
    _pencil_checkpoint_generation_path(dir, name, _PENCIL_HYBRID_GENERATION, "Pencil hybrid")
_pencil_hybrid_digest(value) = _pencil_checkpoint_digest(value, "Pencil hybrid")
_pencil_hybrid_file_digest(path) = _pencil_checkpoint_file_digest(path, "Pencil hybrid")
_pencil_hybrid_broadcast_generation(name, comm, rank) = _pencil_checkpoint_broadcast_generation(
    name,
    comm,
    rank,
    _PENCIL_HYBRID_GENERATION,
    "Pencil hybrid",
)
_pencil_hybrid_remove_generation(dir, generation, nranks) = _pencil_checkpoint_remove_generation(
    dir,
    generation,
    nranks,
    _PENCIL_HYBRID_GENERATION,
    "Pencil hybrid",
)
_pencil_hybrid_rank_state(path, digest) =
    _pencil_checkpoint_rank_state(path, digest, "Pencil hybrid")

function _pencil_hybrid_manifest(path, model, nranks)
    envelope = HybridPlasmaPIC.Serialization.deserialize(path)
    _pencil_hybrid_checkpoint_keys(envelope, (:schema, :model, :payload, :sha256))
    envelope.schema === 1 && envelope.model === :pencil_hybrid_manifest ||
        throw(ArgumentError("unsupported Pencil hybrid manifest envelope"))
    envelope.payload isa Vector{UInt8} ||
        throw(ArgumentError("Pencil hybrid manifest payload must be serialized bytes"))
    digest = _pencil_hybrid_digest(envelope.sha256)
    bytes2hex(HybridPlasmaPIC.SHA.sha256(envelope.payload)) == digest ||
        throw(ArgumentError("Pencil hybrid manifest checksum mismatch"))
    record = HybridPlasmaPIC._deserialize_mpi_payload(NamedTuple, envelope.payload)
    _pencil_hybrid_checkpoint_keys(
        record,
        (
            keys(_pencil_hybrid_checkpoint_config(model, nranks))...,
            _PENCIL_HYBRID_CHECKPOINT_PHASE...,
            :generation,
            :rank_files,
            :rank_sha256,
        ),
    )
    _pencil_hybrid_checkpoint_common(record, model, nranks)
    record.generation isa String || throw(ArgumentError("checkpoint generation must be a String"))
    _pencil_hybrid_generation_name(record.generation)
    files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    record.rank_files isa Tuple && isequal(record.rank_files, files) ||
        throw(ArgumentError("Pencil hybrid manifest requires canonical ordered rank filenames"))
    record.rank_sha256 isa Tuple && length(record.rank_sha256) == nranks ||
        throw(ArgumentError("Pencil hybrid manifest requires one digest per rank"))
    foreach(_pencil_hybrid_digest, record.rank_sha256)
    return record
end

function _pencil_hybrid_envelope(manifest)
    payload = HybridPlasmaPIC._serialize_mpi_payload(manifest)
    HybridPlasmaPIC._check_mpi_count(length(payload), "Pencil hybrid manifest payload")
    return (
        schema = 1,
        model = :pencil_hybrid_manifest,
        payload = payload,
        sha256 = bytes2hex(HybridPlasmaPIC.SHA.sha256(payload)),
    )
end

# Preserve the canonical named physical inventory, copying every owned parent.
_pencil_hybrid_checkpoint_fields(fields) =
    map(fields) do value
        value isa Tuple ? map(a -> copy(parent(a)), value) : copy(parent(value))
    end

"""
    save_mpi_checkpoint(dir, model::DistributedHybrid) -> manifest_path

Collectively save all eighteen carried isothermal hybrid fields, complete ordered
local particles and actual initialized/primed/last_dt/time/step. Matching cached
physics/qm, layout, counts and finite density>=nfloor are required. Predictor ui/E
and stored half-time velocities persist verbatim; no field reconstruction occurs.

Each rank file closes before a checksummed manifest atomically publishes its unique
generation. Earlier failures preserve the previous manifest; superseded cleanup
failure reports the new checkpoint usable. Shared filesystem, exclusive access,
trusted Serialization and recoverable-local/MPI/fatal-resource limits apply.
"""
function HybridPlasmaPIC.save_mpi_checkpoint(dir, model::DistributedHybrid)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil hybrid checkpoint"
    raw_dir, common = _pencil_hybrid_checkpoint_agreement(model, dir, nranks)
    _pencil_hybrid_agreement(model, model.owners, context)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "save particle ownership") do
        _pencil_owned_particles(model.owners, model.g, model._layout)
    end
    state, phase = HybridPlasmaPIC._mpi_collective_stage(comm, context, "snapshot preparation") do
        phase = (
            initialized = model.initialized[],
            primed = model.primed[],
            last_dt = model.last_dt[],
            time = model.time[],
            step = model.step[],
        )
        state = merge(
            common,
            phase,
            (rank = rank,),
            _pencil_hybrid_checkpoint_layout(model),
            (
                fields = _pencil_hybrid_checkpoint_fields(model.fields),
                species = deepcopy([
                    HybridPlasmaPIC._checkpoint_particle_state(ps) for ps in model.owners
                ]),
            ),
        )
        state, phase
    end
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "directory setup") do
        rank == 0 && mkpath(raw_dir)
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest_path = joinpath(directory, _PENCIL_HYBRID_MANIFEST)
    previous, owned_generation =
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation setup") do
            if rank == 0
                old =
                    ispath(manifest_path) ? _pencil_hybrid_manifest(manifest_path, model, nranks) :
                    nothing
                if old !== nothing
                    location = _pencil_hybrid_generation_path(directory, old.generation)
                    for r = 1:nranks
                        _pencil_hybrid_file_digest(joinpath(location, old.rank_files[r])) ==
                        old.rank_sha256[r] || throw(
                            ArgumentError(
                                "superseded Pencil hybrid checkpoint rank checksum mismatch",
                            ),
                        )
                    end
                end
                old,
                basename(mktempdir(directory; prefix = _PENCIL_HYBRID_GENERATION, cleanup = false))
            else
                nothing, ""
            end
        end
    try
        generation = _pencil_hybrid_broadcast_generation(owned_generation, comm, rank)
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank write") do
            location = _pencil_hybrid_generation_path(directory, generation)
            path = joinpath(location, HybridPlasmaPIC._mpi_checkpoint_rank_filename(rank + 1))
            !ispath(path) && !islink(path) ||
                throw(ArgumentError("new Pencil hybrid rank file already exists"))
            HybridPlasmaPIC._atomic_serialize(path, state)
        end
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "manifest publication") do
            if rank == 0
                location = _pencil_hybrid_generation_path(directory, owned_generation)
                files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
                digests = map(name -> _pencil_hybrid_file_digest(joinpath(location, name)), files)
                manifest = merge(
                    common,
                    phase,
                    (generation = owned_generation, rank_files = files, rank_sha256 = digests),
                )
                HybridPlasmaPIC._atomic_serialize(manifest_path, _pencil_hybrid_envelope(manifest))
            end
        end
    catch original
        (
            original isa InterruptException ||
            original isa OutOfMemoryError ||
            original isa StackOverflowError
        ) && rethrow()
        trace = catch_backtrace()
        try
            HybridPlasmaPIC._mpi_collective_stage(comm, context, "failed generation cleanup") do
                rank == 0 && _pencil_hybrid_remove_generation(directory, owned_generation, nranks)
            end
        catch cleanup_error
            (
                cleanup_error isa InterruptException ||
                cleanup_error isa OutOfMemoryError ||
                cleanup_error isa StackOverflowError
            ) && rethrow()
            throw(
                CompositeException([
                    CapturedException(original, trace),
                    CapturedException(cleanup_error, catch_backtrace()),
                ]),
            )
        end
        rethrow()
    end
    try
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "superseded generation cleanup") do
            rank == 0 &&
                previous !== nothing &&
                _pencil_hybrid_remove_generation(directory, previous.generation, nranks)
        end
    catch cleanup_error
        (
            cleanup_error isa InterruptException ||
            cleanup_error isa OutOfMemoryError ||
            cleanup_error isa StackOverflowError
        ) && rethrow()
        throw(
            CompositeException([
                ErrorException(
                    "Pencil hybrid checkpoint was published, but superseded-file cleanup failed",
                ),
                CapturedException(cleanup_error, catch_backtrace()),
            ]),
        )
    end
    return manifest_path
end

function _pencil_hybrid_prepare_rank_state(
    state,
    manifest,
    model::DistributedHybrid{T,G,R,P},
    rank,
    nranks,
) where {T,G,R,P}
    layout = _pencil_hybrid_checkpoint_layout(model)
    _pencil_hybrid_checkpoint_keys(
        state,
        (
            keys(_pencil_hybrid_checkpoint_config(model, nranks))...,
            _PENCIL_HYBRID_CHECKPOINT_PHASE...,
            :rank,
            keys(layout)...,
            :fields,
            :species,
        ),
    )
    _pencil_hybrid_checkpoint_common(state, model, nranks)
    state.rank isa Int && state.rank == rank ||
        throw(ArgumentError("Pencil hybrid checkpoint rank mismatch"))
    for name in _PENCIL_HYBRID_CHECKPOINT_PHASE
        isequal(getproperty(state, name), getproperty(manifest, name)) ||
            throw(ArgumentError("Pencil hybrid rank phase differs from manifest"))
    end
    for name in keys(layout)
        value, expected = getproperty(state, name), getproperty(layout, name)
        typeof(value) === typeof(expected) && isequal(value, expected) ||
            throw(ArgumentError("Pencil hybrid checkpoint $name differs from current layout"))
    end
    _pencil_hybrid_checkpoint_keys(state.fields, _PENCIL_HYBRID_FIELD_NAMES)
    for name in (:B, :E, :ui, :J, :gradp)
        getproperty(state.fields, name) isa NTuple{3,Array{T,3}} || throw(
            ArgumentError(
                "Pencil hybrid checkpoint requires ordinary three-component Array fields",
            ),
        )
    end
    for name in (:n, :pe, :ninv)
        getproperty(state.fields, name) isa Array{T,3} ||
            throw(ArgumentError("Pencil hybrid checkpoint requires ordinary scalar Array fields"))
    end
    fields = _pencil_hybrid_field_arrays(state.fields)
    destination = (
        _pencil_hybrid_arrays(model)...,
        _pencil_hybrid_metadata(model)...,
        _pencil_hybrid_particle_work(model)...,
    )
    for (i, array) in enumerate(fields)
        size(array) == layout.parent_size ||
            throw(DimensionMismatch("Pencil hybrid incoming field shape mismatch"))
        for j = 1:(i-1)
            HybridPlasmaPIC._particle_arrays_mightalias(array, fields[j]) &&
                throw(ArgumentError("Pencil hybrid incoming fields must be independent"))
        end
        for other in destination
            HybridPlasmaPIC._particle_arrays_mightalias(array, parent(other)) && throw(
                ArgumentError("Pencil hybrid incoming fields must be detached from model storage"),
            )
        end
        _aliases_workspace(array, model.g) && throw(
            ArgumentError("Pencil hybrid incoming fields must be detached from Fourier storage"),
        )
    end
    _pencil_hybrid_finite_fields(state.fields, model.nfloor)
    state.species isa Vector && length(state.species) == length(model.owners) || throw(
        ArgumentError("Pencil hybrid checkpoint requires matching ordered population records"),
    )
    populations = P[]
    for (s, record) in enumerate(state.species)
        _pencil_hybrid_checkpoint_keys(record, HybridPlasmaPIC._CHECKPOINT_PARTICLE_FIELDS)
        ps = HybridPlasmaPIC._checkpoint_species(record, Val(3), T)
        ps isa P || throw(
            ArgumentError("Pencil hybrid checkpoint requires ordinary Vector particle leaves"),
        )
        isequal((ps.q, ps.m), model._species_scalars[s]) ||
            throw(ArgumentError("Pencil hybrid incoming q/m differs from cached binding"))
        push!(populations, ps)
    end
    _pencil_hybrid_species(
        model,
        populations;
        bound = false,
        values = true,
        factors = true,
        counts = false,
    )
    _pencil_owned_particles(populations, model.g, model._layout)
    for ps in populations, live in model.owners
        HybridPlasmaPIC._particle_sets_mightalias(ps, live) && throw(
            ArgumentError("Pencil hybrid incoming particles must be detached from destination"),
        )
    end
    for array in fields, species in (populations, model.owners, model._staged), ps in species
        HybridPlasmaPIC._mightalias_particle_storage(array, ps) &&
            throw(ArgumentError("Pencil hybrid incoming fields must not alias particle storage"))
    end
    prepared = P[
        HybridPlasmaPIC._prepare_particle_replacement(ps, eachindex(ps.weight)) for
        ps in populations
    ]
    scalars = Tuple{T,T}[(ps.q, ps.m) for ps in prepared]
    counts = Int[length(ps.weight) for ps in prepared]
    targets = map(parent, _pencil_hybrid_field_arrays(model.fields))
    return (
        state = state,
        particles = prepared,
        scalars = scalars,
        counts = counts,
        fields = fields,
        targets = targets,
    )
end

"""
    load_mpi_checkpoint!(model::DistributedHybrid, dir) -> model

Restore matching rank/layout, exact typed cached physics/qm and complete carried
fields/particles/phase. Old destination field/qm/count/cache/clock/phase values may
be stale; ordinary representation and immutable configuration remain required.
Incoming fields/particles detach and validate, and replacement leaves/count rows
prepare on every rank before any public store. No reconstruction, deposition,
Poisson, curl, Ohm, predictor, recentering, priming or field projection occurs.
"""
function HybridPlasmaPIC.load_mpi_checkpoint!(model::DistributedHybrid, dir)
    comm = model.g.fft.comm
    rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
    context = "Pencil hybrid checkpoint"
    raw_dir, _ = _pencil_hybrid_checkpoint_agreement(model, dir, nranks)
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "restore destination storage") do
        _pencil_hybrid_storage(model)
        _pencil_hybrid_species(
            model,
            model.owners;
            bound = true,
            values = false,
            factors = false,
            counts = false,
        )
    end
    directory = first(_agree_configuration(comm, context) do
        (realpath(raw_dir),)
    end)
    manifest = first(
        _agree_configuration(comm, context) do
            record = _pencil_hybrid_manifest(
                joinpath(directory, _PENCIL_HYBRID_MANIFEST),
                model,
                nranks,
            )
            _pencil_hybrid_generation_path(directory, record.generation)
            (record,)
        end,
    )
    saved = HybridPlasmaPIC._mpi_collective_stage(comm, context, "rank restore preparation") do
        location = _pencil_hybrid_generation_path(directory, manifest.generation)
        path = joinpath(location, manifest.rank_files[rank+1])
        incoming = _pencil_hybrid_rank_state(path, manifest.rank_sha256[rank+1])
        _pencil_hybrid_prepare_rank_state(incoming, manifest, model, rank, nranks)
        detached = deepcopy(incoming)
        _pencil_hybrid_prepare_rank_state(detached, manifest, model, rank, nranks)
    end
    # This explicit final allrank gate has no reconstruction or neutrality policy.
    HybridPlasmaPIC._mpi_collective_stage(comm, context, "complete publication preflight") do
        _pencil_hybrid_storage(model)
        _pencil_hybrid_species(
            model,
            model.owners;
            bound = true,
            values = false,
            factors = false,
            counts = false,
        )
    end
    for s in eachindex(model.owners)
        replacement = saved.particles[s]
        HybridPlasmaPIC._publish_particle_replacement!(model.owners[s], replacement)
        model.owners[s].q = replacement.q
        model.owners[s].m = replacement.m
    end
    for c = 1:18
        copyto!(saved.targets[c], saved.fields[c])
    end
    copyto!(model._species_scalars, saved.scalars)
    copyto!(model._population_counts, saved.counts)
    model.initialized[] = saved.state.initialized
    model.primed[] = saved.state.primed
    model.last_dt[] = saved.state.last_dt
    model.time[] = saved.state.time
    model.step[] = saved.state.step
    model._floor_count[] = 0
    return model
end
