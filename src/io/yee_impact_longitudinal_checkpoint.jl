# The impact map has its own completed-state schema. Physical configuration,
# whole-payload detachment and file publication keep their closed-Yee owners.
const _YI_CHECKPOINT_POLICIES = (
    gauss_policy = :original_source_exact_rounded,
    plate_policy = :isolated_left_derived_right,
    force_policy = :original_geometry_boris_nearest_even,
    force_policy_version = 1,
    integrator_policy = :event_local_longitudinal_kdk,
    integrator_policy_version = 1,
    duration_policy = :exact_rational_chunks,
    duration_policy_version = 1,
    endpoint_policy = :four_ulp_constrained_endpoint,
    endpoint_policy_version = 1,
    event_policy = :one_cluster_per_external_step,
    event_policy_version = 1,
    coverage_policy = :external_start_carried_original_flow,
    coverage_policy_version = 1,
    reflection_policy = :elastic_center_wall,
    reflection_policy_version = 1,
    velocity_phase = :synchronized,
    field_phase = :synchronized,
    phase_policy = :strict_completed_interior,
    phase_policy_version = 1,
    timestep_policy = :fixed_after_first_success,
    timestep_policy_version = 1,
    arithmetic_policy = :nearest_gradual,
    arithmetic_policy_version = 1,
    wall_transfer_policy = :fixed_elastic_zero_charge_work,
    wall_transfer_policy_version = 1,
)
const _YI_CHECKPOINT_FIELDS = (
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
    keys(_YI_CHECKPOINT_POLICIES)...,
    :max_root_refinements,
    :particle_layout,
    :species,
    :fields,
    :phase,
    :right_plate_charge,
    :wall_charge_transfer,
    :wall_work_transfer,
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
)

function _yi_validate_checkpoint_configuration(s, st::YeeImpactLongitudinalPIC)
    _yee_closed_validate_checkpoint_configuration(s, st)
    all(
        name -> getproperty(s, name) === getproperty(_YI_CHECKPOINT_POLICIES, name),
        keys(_YI_CHECKPOINT_POLICIES),
    ) || throw(ArgumentError("unsupported Yee impact checkpoint policy or version"))
    s.max_root_refinements isa Int &&
        s.max_root_refinements > 0 &&
        s.max_root_refinements === st.max_root_refinements ||
        throw(ArgumentError("checkpoint impact root budget does not match destination"))
    for name in (:wall_charge_transfer, :wall_work_transfer)
        values = getproperty(s, name)
        values isa NTuple{2,Rational{BigInt}} &&
            all(x -> denominator(x) > 0 && iszero(x), values) ||
            throw(ArgumentError("checkpoint elastic wall charge/work transfers must be exact zero"))
    end
    return nothing
end

function _validate_checkpoint_state(s, st::YeeImpactLongitudinalPIC{T}) where {T}
    _validate_checkpoint_container(s, (:integrator, :schema_version))
    s.integrator === :yee_impact_longitudinal_pic && s.schema_version === 1 ||
        throw(ArgumentError("unsupported Yee impact checkpoint integrator or schema"))
    _validate_checkpoint_container(s, _YI_CHECKPOINT_FIELDS)
    _yi_validate_checkpoint_configuration(s, st)
    s.initialized === true && s.phase === :ready ||
        throw(ArgumentError("Yee impact checkpoint must be an initialized completed boundary"))
    populations = _yl_checkpoint_populations(s, T)
    _yl_validate_populations(st, populations)
    count = 0
    for ps in populations
        count = Base.checked_add(count, length(ps.weight))
    end
    _yi_validate_history_values(s, st.g, count)
    fields = _yee_closed_checkpoint_arrays(s, st, populations)
    _yl_validate_longitudinal_fields(st, fields)
    s.right_plate_charge isa Rational{BigInt} &&
        denominator(s.right_plate_charge) > 0 &&
        s.right_plate_charge == _yl_right_plate(st, populations) ||
        throw(ArgumentError("checkpoint derived right plate charge is inconsistent with sources"))
    return nothing
end

function _yi_checkpoint_state(st::YeeImpactLongitudinalPIC{T}, species, layout) where {T}
    populations = _yl_species(species)
    _yl_validate_populations(st, populations)
    _yi_validate_phase(st)
    st.initialized[] || throw(ArgumentError("Yee impact checkpoints require initialized state"))
    _yi_validate_bound(st, populations)
    _yl_validate_longitudinal_fields(st, st.fields)
    config = (
        integrator = :yee_impact_longitudinal_pic,
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
        right_plate_charge = st.right_plate_charge[],
        wall_charge_transfer = ntuple(_ -> zero(Rational{BigInt}), 2),
        wall_work_transfer = ntuple(_ -> zero(Rational{BigInt}), 2),
    )
    return merge(config, _YI_CHECKPOINT_POLICIES, _yi_history_values(st))
end

_checkpoint_state(st::YeeImpactLongitudinalPIC, ps::ParticleSet) =
    _yi_checkpoint_state(st, ps, :scalar)
_checkpoint_state(st::YeeImpactLongitudinalPIC, species::AbstractVector) =
    _yi_checkpoint_state(st, species, :collection)

_yee_closed_make_stage(st::YeeImpactLongitudinalPIC, populations, sources, T) =
    _yi_make_stage(populations, sources, T, st.g.n[1])

function _yi_restore_checkpoint!(st::YeeImpactLongitudinalPIC, species, s, layout)
    populations, saved, stage = _yee_closed_checkpoint_restore_data(st, species, s, layout)
    _yee_closed_restore_checkpoint_values!(st, populations, saved, stage)
    _yi_publish_extra_history_refs!(st, saved)
    return st
end

_restore_checkpoint_state!(st::YeeImpactLongitudinalPIC, ps::ParticleSet, s) =
    _yi_restore_checkpoint!(st, ps, s, :scalar)
_restore_checkpoint_state!(st::YeeImpactLongitudinalPIC, species::AbstractVector, s) =
    _yi_restore_checkpoint!(st, species, s, :collection)

function save_checkpoint(path::AbstractString, st::YeeImpactLongitudinalPIC, species)
    return _save_checkpoint(path, st, _checkpoint_state(st, species))
end

function load_checkpoint!(st::YeeImpactLongitudinalPIC, species, path::AbstractString)
    return _restore_checkpoint_state!(st, species, deserialize(path))
end

function archive_run(
    path::AbstractString,
    st::YeeImpactLongitudinalPIC,
    species;
    boundary_desc::AbstractString = "closed PEC",
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
