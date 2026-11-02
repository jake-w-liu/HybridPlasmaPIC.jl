# Complete static-composite physical state. Geometry, source and factor policies
# are fixed; no field/phase/history is silently regenerated during publication.
const _STATIC_COMPOSITE_ES_CHECKPOINT_FIELDS = (
    :integrator,
    :schema_version,
    :configuration,
    :particle_layout,
    :species,
    :E,
    :phi,
    :initialized,
    :phase,
    :time,
    :step,
    :last_dt,
    :right_plate_charge,
    :last_impulse,
    :last_report,
)

function _sce_checkpoint_keys(s, required)
    _validate_checkpoint_container(s, required)
    length(keys(s)) == length(required) ||
        throw(ArgumentError("static composite checkpoint contains unsupported fields"))
    return nothing
end

function _sce_restore_storage(st, species)
    _sce_mode()
    _sce_storage(st, species)
    _sce_configuration(st)
    # Only the current owner/component identities are bindings here. Old values,
    # factors, counts, phase, fields and report will be replaced by saved state.
    _sce_binding(st, species; values = false, counts = false)
    return nothing
end

function _sce_checkpoint_data(s, st::StaticCompositeElectrostaticPIC{T}) where {T}
    _sce_restore_storage(st, st._binding[].owners)
    _sce_checkpoint_keys(s, _STATIC_COMPOSITE_ES_CHECKPOINT_FIELDS)
    s.integrator === :static_composite_es &&
        s.schema_version === 1 &&
        s.particle_layout === :collection &&
        s.initialized === true &&
        s.phase === :synchronous_integer ||
        throw(ArgumentError("unsupported static composite checkpoint family/schema/phase"))
    configuration = _sce_configuration(st)
    typeof(s.configuration) === typeof(configuration) && isequal(s.configuration, configuration) ||
        throw(
            ArgumentError(
                "static composite checkpoint hierarchy/source/factor policies must match",
            ),
        )
    s.species isa Vector && length(s.species) == length(st._factors) ||
        throw(ArgumentError("static composite checkpoint requires its ordered species count"))
    populations = _ElectrostaticHostParticle{1,T}[]
    for record in s.species
        _sce_checkpoint_keys(record, _CHECKPOINT_PARTICLE_FIELDS)
        ps = _checkpoint_species(record, Val(1), T)
        ps isa _ElectrostaticHostParticle{1,T} ||
            throw(ArgumentError("static composite checkpoint particles must use ordinary Vectors"))
        push!(populations, ps)
    end
    _sce_population_storage(populations, T)
    _sce_values(st.mesh, populations, st._factors)
    total = 0
    for ps in populations
        total = Base.checked_add(total, length(ps.weight))
    end
    Base.checked_add(total, 1)
    # Admit all scalar/Q/history types before inspecting record-owned arrays or
    # evaluating the original-source field map.
    _sce_phase_values(
        st,
        s.initialized,
        s.time,
        s.step,
        s.last_dt,
        s.last_impulse,
        s.last_report,
        total,
    )
    s.E isa Vector{T} &&
        s.phi isa Vector{T} &&
        length(s.E) == length(st.E) &&
        length(s.phi) == length(st.phi) ||
        throw(ArgumentError("static composite checkpoint fields require ordinary matching Vectors"))
    incoming = Any[s.E, s.phi, s.last_impulse]
    for ps in populations
        append!(incoming, _sce_particle_arrays(ps))
    end
    s.last_report isa StaticCompositeESReport && push!(incoming, s.last_report.interface_crossings)
    writable, readonly = _sce_arrays(st, st._binding[].owners)
    destination = vcat(writable, readonly)
    for i in eachindex(incoming)
        for j = 1:(i-1)
            _particle_arrays_mightalias(incoming[i], incoming[j]) && throw(
                ArgumentError(
                    "static composite checkpoint physical/history arrays must be independent",
                ),
            )
        end
        for a in destination
            _particle_arrays_mightalias(incoming[i], a) && throw(
                ArgumentError(
                    "static composite checkpoint must be detached from live/private storage",
                ),
            )
        end
    end
    _sce_state_values(
        st,
        populations,
        s.E,
        s.phi,
        s.initialized,
        s.time,
        s.step,
        s.last_dt,
        s.right_plate_charge,
        s.last_impulse,
        s.last_report,
    )
    return populations
end

function _validate_checkpoint_state(s, st::StaticCompositeElectrostaticPIC)
    _sce_checkpoint_data(s, st)
    return nothing
end

function _checkpoint_state(st::StaticCompositeElectrostaticPIC, species::AbstractVector)
    _sce_ready(st, species)
    state = (
        integrator = :static_composite_es,
        schema_version = 1,
        configuration = _sce_configuration(st),
        particle_layout = :collection,
        species = [_checkpoint_particle_state(ps) for ps in species],
        E = st.E,
        phi = st.phi,
        initialized = true,
        phase = :synchronous_integer,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        right_plate_charge = st.right_plate_charge[],
        last_impulse = st.last_impulse,
        last_report = st.last_report[],
    )
    # Closed concrete scalar records allow legitimate repeated BigInt values;
    # all writable array roles are detached as one graph before validation.
    saved = deepcopy(state)
    _validate_checkpoint_state(saved, st)
    return saved
end

function _restore_checkpoint_state!(
    st::StaticCompositeElectrostaticPIC{T},
    species::AbstractVector,
    s,
) where {T}
    _sce_restore_storage(st, species)
    # Original links are checked before copying can hide a retained live alias.
    _sce_checkpoint_data(s, st)
    saved = deepcopy(s)
    incoming = _sce_checkpoint_data(saved, st)
    binding = _sce_make_binding(species, incoming, st.mesh, T; components = incoming)
    _sce_validate_bins(binding.bins, incoming, st.mesh)
    # Incoming ordinary particle types are exactly the destination type; every
    # replacement leaf and binding/report allocation exists before this point.
    for (destination, prepared) in zip(species, incoming)
        typeof(destination) === typeof(prepared) ||
            throw(ArgumentError("static composite checkpoint replacement storage types must match"))
    end
    for (destination, prepared) in zip(species, incoming)
        _publish_particle_replacement!(destination, prepared)
        destination.q, destination.m = prepared.q, prepared.m
    end
    copyto!(st.E, saved.E)
    copyto!(st.phi, saved.phi)
    copyto!(st.last_impulse, saved.last_impulse)
    st._binding[] = binding
    st.right_plate_charge[] = saved.right_plate_charge
    st.last_report[] = saved.last_report
    st.time[], st.step[], st.last_dt[] = saved.time, saved.step, saved.last_dt
    st.initialized[] = true
    return st
end

"""
    save_checkpoint(path, state::StaticCompositeElectrostaticPIC, species)

Atomically save complete distinct schema1 `:static_composite_es`: immutable
hierarchy/source/factor configuration, particles, both physical fields, clock,
exact integrated current/derived plate and full signed report. Ready boundaries
only. Raw serialization has no separate checksum; use `archive_run` for one.
"""
function save_checkpoint(
    path::AbstractString,
    st::StaticCompositeElectrostaticPIC,
    species::AbstractVector,
)
    return _save_checkpoint(path, st, _checkpoint_state(st, species))
end

"""
    load_checkpoint!(state::StaticCompositeElectrostaticPIC, species, path)

Validate the same hierarchy and immutable factors, detach and prepare all leaves
and caches, then restore the saved physical fields/history verbatim. Old bound
particle values/counts, fields and clock may be stale; storage/owner identities
must remain valid. No reinitialization, projection or half-step priming occurs.
"""
function load_checkpoint!(
    st::StaticCompositeElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    return _restore_checkpoint_state!(st, species, deserialize(path))
end

function archive_run(
    path::AbstractString,
    st::StaticCompositeElectrostaticPIC,
    species::AbstractVector;
    boundary_desc::AbstractString = "closed isolated composite 1D",
    kwargs...,
)
    return _archive_run(path, st, _checkpoint_state(st, species); boundary_desc, kwargs...)
end

"""
    load_archive!(state::StaticCompositeElectrostaticPIC, species, path)

Verify the existing run-container schema and byte checksum, then use the same
prepared restore as raw checkpoints. `load_archive(path)` inspects saved metadata.
"""
function load_archive!(
    st::StaticCompositeElectrostaticPIC,
    species::AbstractVector,
    path::AbstractString,
)
    loaded = load_run(path)
    return _restore_checkpoint_state!(st, species, loaded.state)
end
