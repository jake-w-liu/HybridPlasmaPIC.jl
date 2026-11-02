# Synchronized full-vector cloud PIC has a distinct schema. Closed-source
# configuration/storage/restore mechanics are shared with longitudinal IO;
# file publication, archive metadata and checksums keep their existing owners.
const _YE_CHECKPOINT_POLICY_VERSION = 1
const _YE_CHECKPOINT_POLICY_FIELDS = (
    :gather_policy_version,
    :current_policy_version,
    :pec_policy_version,
    :event_policy_version,
    :integrator_policy_version,
    :arithmetic_policy_version,
)
const _YE_CHECKPOINT_FIELDS = (
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
    :c,
    :c_type,
    :c_bits,
    :gauss_policy,
    :plate_policy,
    :gather_policy,
    :current_policy,
    :pec_policy,
    :event_policy,
    :integrator_policy,
    :velocity_phase,
    :field_phase,
    :timestep_policy,
    :arithmetic_policy,
    :force_policy,
    :force_policy_version,
    _YE_CHECKPOINT_POLICY_FIELDS...,
    :particle_layout,
    :species,
    :fields,
    :initialized,
    :phase,
    :time,
    :step,
    :last_dt,
    :right_plate_charge,
    :last_impulse,
    :last_work,
    :last_modified_energy_defect,
    :last_physical_energy_change,
)

function _ye_validate_checkpoint_configuration(s, st::YeeElectromagneticPIC)
    _yee_closed_validate_checkpoint_configuration(s, st)
    _yl_validate_checkpoint_ieee(:c, s.c, s.c_type, s.c_bits, st.c)
    s.gauss_policy === :original_source_exact_rounded &&
        s.plate_policy === :isolated_left_derived_right &&
        s.gather_policy === :matched_cloud_edge_nodal &&
        s.current_policy === :original_source_path_exact_rounded &&
        s.pec_policy === :stored_zero_tangent_E_Bx &&
        s.event_policy === :strict_event_free &&
        s.integrator_policy === :synchronized_boris_kdk_yee &&
        s.velocity_phase === :synchronized &&
        s.field_phase === :synchronized &&
        s.timestep_policy === :fixed_after_first_success &&
        s.arithmetic_policy === :nearest_gradual || throw(
        ArgumentError("checkpoint Yee electromagnetic source or stepping policy is unsupported"),
    )
    all(
        name -> getproperty(s, name) === _YE_CHECKPOINT_POLICY_VERSION,
        _YE_CHECKPOINT_POLICY_FIELDS,
    ) || throw(ArgumentError("unsupported Yee electromagnetic checkpoint policy version"))
    return nothing
end

function _validate_checkpoint_state(s, st::YeeElectromagneticPIC{T}) where {T}
    _yee_closed_validate_checkpoint_envelope(s, :yee_electromagnetic_pic)
    _validate_checkpoint_container(s, _YE_CHECKPOINT_FIELDS)
    _ye_validate_checkpoint_configuration(s, st)
    s.initialized === true && s.phase === :ready || throw(
        ArgumentError("Yee electromagnetic checkpoint must be an initialized completed boundary"),
    )
    s.time isa T &&
        s.step isa Int &&
        s.last_dt isa T &&
        s.last_work isa T &&
        s.last_modified_energy_defect isa T &&
        s.last_physical_energy_change isa T ||
        throw(ArgumentError("checkpoint clock and diagnostics must have destination precision"))
    # Re-derive the owner bound from immutable grid/c, not a saved cache.
    _ye_validate_phase_values(
        s.initialized,
        s.time,
        s.step,
        s.last_dt,
        s.last_work,
        s.last_modified_energy_defect,
        s.last_physical_energy_change,
        s.last_impulse,
        st.g,
        _yee_cfl(st.g, st.c),
    )
    populations = _yl_checkpoint_populations(s, T)
    _yl_validate_populations(st, populations)
    fields = _yee_closed_checkpoint_arrays(s, st, populations)
    _ye_validate_fields(st, fields)
    s.right_plate_charge isa Rational{BigInt} &&
        s.right_plate_charge == _yl_right_plate(st, populations) ||
        throw(ArgumentError("checkpoint derived right plate charge is inconsistent with sources"))
    return nothing
end

function _ye_checkpoint_state(st::YeeElectromagneticPIC{T}, species, layout) where {T}
    populations = _yl_species(species)
    _yl_validate_populations(st, populations)
    _ye_validate_phase(st)
    st.initialized[] ||
        throw(ArgumentError("Yee electromagnetic checkpoints require initialized state"))
    _yl_validate_bound(st, populations)
    _ye_validate_fields(st, st.fields)
    return (
        integrator = :yee_electromagnetic_pic,
        schema_version = 2,
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
        c = st.c,
        c_type = typeof(st.c),
        c_bits = _yl_checkpoint_bits(st.c),
        gauss_policy = :original_source_exact_rounded,
        plate_policy = :isolated_left_derived_right,
        gather_policy = :matched_cloud_edge_nodal,
        current_policy = :original_source_path_exact_rounded,
        pec_policy = :stored_zero_tangent_E_Bx,
        event_policy = :strict_event_free,
        integrator_policy = :synchronized_boris_kdk_yee,
        velocity_phase = :synchronized,
        field_phase = :synchronized,
        timestep_policy = :fixed_after_first_success,
        arithmetic_policy = :nearest_gradual,
        force_policy = :original_geometry_boris_nearest_even,
        force_policy_version = 1,
        gather_policy_version = _YE_CHECKPOINT_POLICY_VERSION,
        current_policy_version = _YE_CHECKPOINT_POLICY_VERSION,
        pec_policy_version = _YE_CHECKPOINT_POLICY_VERSION,
        event_policy_version = _YE_CHECKPOINT_POLICY_VERSION,
        integrator_policy_version = _YE_CHECKPOINT_POLICY_VERSION,
        arithmetic_policy_version = _YE_CHECKPOINT_POLICY_VERSION,
        particle_layout = layout,
        species = [_checkpoint_particle_state(ps) for ps in populations],
        fields = _yl_checkpoint_fields(st),
        initialized = true,
        phase = :ready,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
        right_plate_charge = st.right_plate_charge[],
        last_impulse = st.last_impulse,
        last_work = st.last_work[],
        last_modified_energy_defect = st.last_modified_energy_defect[],
        last_physical_energy_change = st.last_physical_energy_change[],
    )
end

_checkpoint_state(st::YeeElectromagneticPIC, ps::ParticleSet) =
    _ye_checkpoint_state(st, ps, :scalar)
_checkpoint_state(st::YeeElectromagneticPIC, species::AbstractVector) =
    _ye_checkpoint_state(st, species, :collection)

function _ye_restore_checkpoint!(st::YeeElectromagneticPIC, species, s, layout)
    populations, saved, stage = _yee_closed_checkpoint_restore_data(st, species, s, layout)
    _yee_closed_restore_checkpoint_values!(st, populations, saved, stage)
    st.last_modified_energy_defect[] = saved.last_modified_energy_defect
    st.last_physical_energy_change[] = saved.last_physical_energy_change
    return st
end

_restore_checkpoint_state!(st::YeeElectromagneticPIC, ps::ParticleSet, s) =
    _ye_restore_checkpoint!(st, ps, s, :scalar)
_restore_checkpoint_state!(st::YeeElectromagneticPIC, species::AbstractVector, s) =
    _ye_restore_checkpoint!(st, species, s, :collection)

function save_checkpoint(path::AbstractString, st::YeeElectromagneticPIC, species)
    return _save_checkpoint(path, st, _checkpoint_state(st, species))
end

function load_checkpoint!(st::YeeElectromagneticPIC, species, path::AbstractString)
    return _restore_checkpoint_state!(st, species, deserialize(path))
end

function archive_run(
    path::AbstractString,
    st::YeeElectromagneticPIC,
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
