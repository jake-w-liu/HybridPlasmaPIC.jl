# One generation contains complete actual-pencil physical state and full rank RNG.
const _PENCIL_CES_MANIFEST = "pencil_collisional_electrostatic_manifest.ser"
const _PENCIL_CES_GENERATION = "_pencil_collisional_electrostatic_"
const _PENCIL_CES_POLICY = (
    prefix = _PENCIL_CES_GENERATION,
    tag = :pencil_collisional_electrostatic_manifest,
    label = "Pencil collisional electrostatic",
)
const _PENCIL_CES_PHASE = (_PENCIL_ES_CHECKPOINT_PHASE..., :phase, :source_time, :source_count)

_pencil_ces_checkpoint_config(st, nranks) = _pencil_collision_checkpoint_config(
    st,
    _pencil_es_checkpoint_config(st.inner, nranks),
    :distributed_collisional_espic,
)

function _pencil_ces_checkpoint_agreement(st, dir, nranks, context)
    return _agree_configuration(st.g.fft.comm, context) do
        _pencil_es_mode(eltype(st.g.L))
        dir isa AbstractString && !isempty(dir) ||
            throw(ArgumentError("checkpoint directory must be a nonempty string"))
        (String(dir), _pencil_ces_checkpoint_config(st, nranks))
    end
end

function _pencil_ces_checkpoint_common(
    record,
    st::DistributedCollisionalElectrostaticPIC{T},
    nranks,
) where {T}
    return _pencil_collision_checkpoint_common(
        record,
        st,
        _pencil_ces_checkpoint_config(st, nranks),
        T,
    )
end

function _pencil_ces_manifest(path, st, nranks)
    envelope = HybridPlasmaPIC.Serialization.deserialize(path)
    _pencil_es_checkpoint_keys(envelope, (:schema, :model, :payload, :sha256))
    envelope.schema === 1 && envelope.model === _PENCIL_CES_POLICY.tag ||
        throw(ArgumentError("unsupported Pencil collisional ES manifest envelope"))
    envelope.payload isa Vector{UInt8} ||
        throw(ArgumentError("checkpoint payload must be serialized bytes"))
    digest = _pencil_checkpoint_digest(envelope.sha256, _PENCIL_CES_POLICY.label)
    bytes2hex(HybridPlasmaPIC.SHA.sha256(envelope.payload)) == digest ||
        throw(ArgumentError("Pencil collisional ES manifest checksum mismatch"))
    record = HybridPlasmaPIC._deserialize_mpi_payload(NamedTuple, envelope.payload)
    _pencil_es_checkpoint_keys(
        record,
        (
            keys(_pencil_ces_checkpoint_config(st, nranks))...,
            _PENCIL_CES_PHASE...,
            :generation,
            :rank_files,
            :rank_sha256,
        ),
    )
    _pencil_ces_checkpoint_common(record, st, nranks)
    record.generation isa String || throw(ArgumentError("checkpoint generation must be a String"))
    _pencil_checkpoint_generation_name(
        record.generation,
        _PENCIL_CES_GENERATION,
        _PENCIL_CES_POLICY.label,
    )
    files = ntuple(r -> HybridPlasmaPIC._mpi_checkpoint_rank_filename(r), nranks)
    record.rank_files isa Tuple && isequal(record.rank_files, files) ||
        throw(ArgumentError("Pencil collisional ES requires canonical ordered rank filenames"))
    record.rank_sha256 isa Tuple && length(record.rank_sha256) == nranks ||
        throw(ArgumentError("Pencil collisional ES requires one rank checksum per rank"))
    for digest in record.rank_sha256
        _pencil_checkpoint_digest(digest, _PENCIL_CES_POLICY.label)
    end
    return record
end

"""
    save_mpi_checkpoint(dir, st::DistributedCollisionalElectrostaticPIC)

Collectively publish one checksummed generation containing complete E/rho,
ordered rank-owned particles, actual clock/half phase and the full per-rank stock
RNG. Requires a completed initialized boundary and matching fixed configuration.
Shared filesystem/exclusive access and trusted Serialization are required. Earlier
failures preserve the prior manifest; superseded cleanup can report a usable new
checkpoint. No crash durability, changed-layout or cross-runtime stream migration.
"""
function HybridPlasmaPIC.save_mpi_checkpoint(dir, st::DistributedCollisionalElectrostaticPIC)
    context = "Pencil collisional electrostatic checkpoint"
    return _pencil_collision_operation(st, context) do
        comm = st.g.fft.comm
        rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
        raw_dir, common = _pencil_ces_checkpoint_agreement(st, dir, nranks, context)
        _pencil_ces_agreement(st, st.species, context)
        HybridPlasmaPIC._mpi_collective_stage(comm, context, "save ownership") do
            _pencil_owned_particles(st.species, st.g, st.inner._layout)
            _pencil_collision_cell_admission(st.inner, st.species)
        end
        _pencil_pic_charges(st.inner, st.species, context)
        state, phase =
            HybridPlasmaPIC._mpi_collective_stage(comm, context, "complete snapshot preparation") do
                phase = _pencil_collision_checkpoint_phase(st)
                state = merge(
                    common,
                    phase,
                    (
                        rank = rank,
                        inner = _pencil_es_checkpoint_snapshot(st.inner, rank, nranks),
                        rng = deepcopy(st.rng),
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
        manifest_path = joinpath(directory, _PENCIL_CES_MANIFEST)
        previous, owned_generation =
            HybridPlasmaPIC._mpi_collective_stage(comm, context, "generation setup") do
                if rank == 0
                    old =
                        ispath(manifest_path) ? _pencil_ces_manifest(manifest_path, st, nranks) :
                        nothing
                    if old !== nothing
                        location = _pencil_checkpoint_generation_path(
                            directory,
                            old.generation,
                            _PENCIL_CES_GENERATION,
                            _PENCIL_CES_POLICY.label,
                        )
                        for r = 1:nranks
                            _pencil_checkpoint_file_digest(
                                joinpath(location, old.rank_files[r]),
                                _PENCIL_CES_POLICY.label,
                            ) == old.rank_sha256[r] || throw(
                                ArgumentError(
                                    "superseded Pencil collisional ES rank checksum mismatch",
                                ),
                            )
                        end
                    end
                    old,
                    basename(
                        mktempdir(directory; prefix = _PENCIL_CES_GENERATION, cleanup = false),
                    )
                else
                    nothing, ""
                end
            end
        _pencil_pic_write_prepared(
            directory,
            manifest_path,
            previous,
            owned_generation,
            state,
            common,
            phase,
            _PENCIL_CES_POLICY,
            comm,
            rank,
            nranks,
            context,
        )
    end
end

function _pencil_ces_prepare_rank(state, manifest, st, rank, nranks)
    _pencil_es_checkpoint_keys(
        state,
        (
            keys(_pencil_ces_checkpoint_config(st, nranks))...,
            _PENCIL_CES_PHASE...,
            :rank,
            :inner,
            :rng,
        ),
    )
    _pencil_ces_checkpoint_common(state, st, nranks)
    state.rank isa Int && state.rank == rank ||
        throw(ArgumentError("Pencil collisional ES checkpoint rank mismatch"))
    for name in _PENCIL_CES_PHASE
        isequal(getproperty(state, name), getproperty(manifest, name)) ||
            throw(ArgumentError("Pencil collisional ES rank phase differs from manifest"))
    end
    _pencil_es_checkpoint_keys(
        state.inner,
        (
            keys(_pencil_es_checkpoint_config(st.inner, nranks))...,
            _PENCIL_ES_CHECKPOINT_PHASE...,
            :rank,
            keys(_pencil_es_checkpoint_layout(st.inner))...,
            :fields,
            :species,
        ),
    )
    for name in _PENCIL_ES_CHECKPOINT_PHASE
        isequal(getproperty(state.inner, name), getproperty(state, name)) ||
            throw(ArgumentError("Pencil collisional ES nested physical phase differs"))
    end
    _pencil_collision_incoming_configuration(st, state.inner, state.coulomb_log, state.rng)
    inner = _pencil_es_prepare_rank_state(state.inner, state.inner, st.inner, rank, nranks)
    _pencil_collision_cell_admission(st.inner, inner.particles)
    _pencil_ces_roles(st, inner.particles)
    return (
        inner = inner,
        rng = state.rng,
        source_time = state.source_time,
        source_count = state.source_count,
    )
end

"""
    load_mpi_checkpoint!(st::DistributedCollisionalElectrostaticPIC, dir)

Collectively restore one complete matching-layout physical/RNG generation.
Destination physical, clock, phase and old RNG values may be stale; valid fixed
configuration and stock writable storage are required. Incoming physical fields
restore verbatim without migration/Poisson; counts may change via prepared leaves.
All validation/detachment and collective gates precede stock copies and ordinary
physical publication. No independent inner restore or RNG sidecar is published.
"""
function HybridPlasmaPIC.load_mpi_checkpoint!(st::DistributedCollisionalElectrostaticPIC, dir)
    context = "Pencil collisional electrostatic checkpoint"
    return _pencil_collision_operation(st, context) do
        comm = st.g.fft.comm
        rank, nranks = Int(MPI.Comm_rank(comm)), Int(MPI.Comm_size(comm))
        raw_dir, _ = _pencil_ces_checkpoint_agreement(st, dir, nranks, context)
        _pencil_ces_agreement(st, st.species, context; ready = false, values = false)
        directory = first(_agree_configuration(comm, context) do
            (realpath(raw_dir),)
        end)
        manifest = first(
            _agree_configuration(comm, context) do
                record =
                    _pencil_ces_manifest(joinpath(directory, _PENCIL_CES_MANIFEST), st, nranks)
                _pencil_checkpoint_generation_path(
                    directory,
                    record.generation,
                    _PENCIL_CES_GENERATION,
                    _PENCIL_CES_POLICY.label,
                )
                (record,)
            end,
        )
        saved = HybridPlasmaPIC._mpi_collective_stage(
            comm,
            context,
            "complete rank restore preparation",
        ) do
            location = _pencil_checkpoint_generation_path(
                directory,
                manifest.generation,
                _PENCIL_CES_GENERATION,
                _PENCIL_CES_POLICY.label,
            )
            path = joinpath(location, manifest.rank_files[rank+1])
            state = _pencil_checkpoint_rank_state(
                path,
                manifest.rank_sha256[rank+1],
                _PENCIL_CES_POLICY.label,
            )
            _pencil_ces_prepare_rank(state, manifest, st, rank, nranks)
            detached = deepcopy(state)
            _pencil_ces_prepare_rank(detached, manifest, st, rank, nranks)
        end
        _pencil_pic_charges(st.inner, saved.inner.particles, context)
        # Every fixed stock buffer size and ordinary replacement is prepared;
        # no validation, allocation, resize, normalization or arithmetic follows.
        copy!(st.rng, saved.rng)
        copy!(st.candidate_rng, saved.rng)
        _pencil_es_restore_prepared!(st.inner, saved.inner)
        st.source_time = saved.source_time
        st.source_count = saved.source_count
        st.phase = :ready
        st
    end
end
