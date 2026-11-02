# Capture has a distinct source/lifecycle schema. Shared closed-Yee owners
# retain configuration bits, whole-payload detachment and file publication.
const _YA_CHECKPOINT_POLICIES = (
    gauss_policy = :original_source_exact_rounded,
    plate_policy = :isolated_dynamic_free_plates,
    plate_policy_version = 1,
    force_policy = :original_geometry_boris_nearest_even,
    force_policy_version = 1,
    integrator_policy = :event_local_longitudinal_kdk,
    integrator_policy_version = 1,
    duration_policy = :exact_rational_chunks,
    duration_policy_version = 1,
    endpoint_policy = :four_ulp_constrained_endpoint,
    endpoint_policy_version = 1,
    event_policy = :one_isolated_capture_per_external_step,
    event_policy_version = 1,
    coverage_policy = :external_start_true_hit_carried_capture,
    coverage_policy_version = 1,
    capture_policy = :center_full_rigid_cloud_to_isolated_plate,
    capture_policy_version = 1,
    passivity_policy = :certified_nonnegative_received_energy,
    passivity_policy_version = 1,
    storage_policy = :ordinary_vector_stable_survivors,
    storage_policy_version = 1,
    capture_record_policy = :last_step_local_elapsed,
    capture_record_policy_version = 1,
    velocity_phase = :synchronized,
    field_phase = :synchronized,
    phase_policy = :strict_completed_interior,
    phase_policy_version = 1,
    timestep_policy = :fixed_after_first_success,
    timestep_policy_version = 1,
    arithmetic_policy = :nearest_gradual,
    arithmetic_policy_version = 1,
)
const _YA_CHECKPOINT_FIELDS = (
    :integrator,
    :schema_version,
    :D,
    :T,
    :ncell,
    :L,
    :dx,
    :periodic,
    :shape,
    :shape_width_type,
    :shape_width_bits,
    :background_charge,
    :background_charge_type,
    :background_charge_bits,
    :left_plate_charge,
    :left_plate_charge_type,
    :left_plate_charge_bits,
    keys(_YA_CHECKPOINT_POLICIES)...,
    :max_root_refinements,
    :particle_layout,
    :species,
    :fields,
    :phase,
    :initialized,
    :time,
    :step,
    :last_dt,
    :last_impulse,
    :last_work,
    :last_energy_defect,
    :last_physical_energy_change,
    :last_wall_count,
    :cumulative_wall_count,
    :last_wall_impulse,
    :cumulative_wall_impulse,
    :last_electric_impulse,
    :last_max_snap,
    :last_snap_energy_bound,
    :last_predictor_time_error,
    :last_physical_time_error,
    :last_source_invariant_error,
    :last_position_error_bound,
    :last_velocity_error_bound,
    :right_plate_charge,
    :free_left_plate_charge,
    :initial_right_plate_charge,
    :initial_counts,
    :initial_charges,
    :last_wall_charge,
    :cumulative_wall_charge,
    :last_wall_kinetic,
    :cumulative_wall_kinetic,
    :last_wall_work,
    :cumulative_wall_work,
    :last_wall_energy,
    :cumulative_wall_energy,
    :last_transport_work,
    :last_passivity_error,
    :last_field_energy_defect,
    :last_total_energy_defect,
    :last_capture,
)

function _ya_validate_checkpoint_configuration(s, st::YeeAbsorbingLongitudinalPIC)
    _yee_closed_validate_checkpoint_configuration(s, st)
    all(name -> getproperty(s, name) === getproperty(_YA_CHECKPOINT_POLICIES, name),
        keys(_YA_CHECKPOINT_POLICIES)) ||
        throw(ArgumentError("unsupported Yee absorbing checkpoint policy or version"))
    s.max_root_refinements isa Int && s.max_root_refinements > 0 &&
        s.max_root_refinements === st.max_root_refinements ||
        throw(ArgumentError("checkpoint absorbing root budget does not match destination"))
    return nothing
end

function _ya_validate_checkpoint_metadata_arrays(s, populations)
    counts, charges = s.initial_counts, s.initial_charges
    counts isa Vector{Int} && charges isa Vector{Rational{BigInt}} ||
        throw(ArgumentError("checkpoint initial source records must use ordinary typed Vectors"))
    metadata = (counts, charges)
    any(_particle_requires_special_stores, metadata) &&
        throw(ArgumentError("checkpoint initial source records require ordinary writes"))
    _particle_arrays_mightalias(counts, charges) &&
        throw(ArgumentError("checkpoint initial source records must be independent"))
    fields = (s.fields.E..., s.fields.B..., s.last_impulse...)
    for a in metadata, b in fields
        _particle_arrays_mightalias(a, b) &&
            throw(ArgumentError("checkpoint initial metadata must not alias fields or impulse"))
    end
    for a in metadata, ps in populations
        _mightalias_particle_storage(a, ps) &&
            throw(ArgumentError("checkpoint initial source records must not alias particles"))
    end
    return nothing
end

function _validate_checkpoint_state(s, st::YeeAbsorbingLongitudinalPIC{T}) where {T}
    _validate_checkpoint_container(s, (:integrator, :schema_version))
    s.integrator === :yee_absorbing_longitudinal_pic && s.schema_version === 1 ||
        throw(ArgumentError("unsupported Yee absorbing checkpoint integrator or schema"))
    _validate_checkpoint_container(s, _YA_CHECKPOINT_FIELDS)
    _ya_validate_checkpoint_configuration(s, st)
    s.initialized === true && s.phase === :ready ||
        throw(ArgumentError("Yee absorbing checkpoint must be an initialized completed boundary"))
    populations = _yl_checkpoint_populations(s, T)
    _yl_validate_populations(st, populations)
    fields = _yee_closed_checkpoint_arrays(s, st, populations)
    _yl_validate_longitudinal_fields(st, fields)
    _ya_validate_checkpoint_metadata_arrays(s, populations)
    # Dynamic plates and initial/live/loss relations belong to the incoming
    # history owner; the destination's old free-left value must not enter them.
    _ya_validate_history_values(s, st, populations)
    return nothing
end

function _ya_checkpoint_state(st::YeeAbsorbingLongitudinalPIC{T}, species, layout) where {T}
    populations = _yl_species(species)
    _yl_validate_populations(st, populations)
    _ya_validate_phase(st, populations)
    st.initialized[] ||
        throw(ArgumentError("Yee absorbing checkpoints require initialized state"))
    _ya_validate_bound(st, populations)
    _yl_validate_longitudinal_fields(st, st.fields)
    config = (
        integrator = :yee_absorbing_longitudinal_pic,
        schema_version = 1,
        D = 1,
        T = T,
        ncell = st.g.n,
        L = st.g.L,
        dx = st.g.dx,
        periodic = st.g.periodic,
        shape = st.shape,
        shape_width_type = typeof(st.shape.half_width),
        shape_width_bits = _yl_checkpoint_bits(st.shape.half_width),
        background_charge = st.background_charge,
        background_charge_type = typeof(st.background_charge),
        background_charge_bits = _yl_checkpoint_bits(st.background_charge),
        left_plate_charge = st.left_plate_charge,
        left_plate_charge_type = typeof(st.left_plate_charge),
        left_plate_charge_bits = _yl_checkpoint_bits(st.left_plate_charge),
        max_root_refinements = st.max_root_refinements,
        particle_layout = layout,
        species = [_checkpoint_particle_state(ps) for ps in populations],
        fields = _yl_checkpoint_fields(st),
        phase = :ready,
    )
    return merge(config, _YA_CHECKPOINT_POLICIES, _ya_history_values(st))
end

_checkpoint_state(st::YeeAbsorbingLongitudinalPIC, ps::ParticleSet) =
    _ya_checkpoint_state(st, ps, :scalar)
_checkpoint_state(st::YeeAbsorbingLongitudinalPIC, species::AbstractVector) =
    _ya_checkpoint_state(st, species, :collection)

_yee_closed_make_stage(st::YeeAbsorbingLongitudinalPIC, populations, sources, T) =
    _ya_make_stages(populations, sources, T, st.g.n[1])

function _ya_restore_checkpoint!(st::YeeAbsorbingLongitudinalPIC, species, s, layout)
    _ya_validate_particle_vectors(st, _yl_species(species))
    populations, saved, prepared = _yee_closed_checkpoint_restore_data(st, species, s, layout)
    _yee_closed_restore_checkpoint_values!(st, populations, saved, prepared.transport)
    st._binding[] = prepared.binding
    _ya_publish_extra_history_refs!(st, saved)
    return st
end

_restore_checkpoint_state!(st::YeeAbsorbingLongitudinalPIC, ps::ParticleSet, s) =
    _ya_restore_checkpoint!(st, ps, s, :scalar)
_restore_checkpoint_state!(st::YeeAbsorbingLongitudinalPIC, species::AbstractVector, s) =
    _ya_restore_checkpoint!(st, species, s, :collection)

function save_checkpoint(path::AbstractString, st::YeeAbsorbingLongitudinalPIC, species)
    return _save_checkpoint(path, st, _checkpoint_state(st, species))
end

function load_checkpoint!(st::YeeAbsorbingLongitudinalPIC, species, path::AbstractString)
    return _restore_checkpoint_state!(st, species, deserialize(path))
end

function archive_run(
    path::AbstractString,
    st::YeeAbsorbingLongitudinalPIC,
    species;
    boundary_desc::AbstractString = "closed PEC with isolated passive capture",
    kwargs...,
)
    return _archive_run(
        path,
        st,
        _checkpoint_state(st, species);
        boundary_desc = boundary_desc,
        kwargs...,
    )
end
