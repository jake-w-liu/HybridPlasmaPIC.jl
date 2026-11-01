# Versioned synchronized longitudinal state. Atomic file publication and archive
# metadata remain owned by checkpoint.jl and restart.jl respectively.
const _YL_CHECKPOINT_FIELDS = (
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
    :gauss_policy,
    :plate_policy,
    :event_policy,
    :integrator_policy,
    :velocity_phase,
    :timestep_policy,
    :arithmetic_policy,
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
    :last_energy_defect,
)

_yl_checkpoint_bits(x::Float32) = reinterpret(UInt32, x)
_yl_checkpoint_bits(x::Float64) = reinterpret(UInt64, x)
_yl_checkpoint_fields(st::YeeLongitudinalPIC) = (E = st.fields.E, B = st.fields.B)

function _yl_validate_checkpoint_ieee(name, value, scalar_type, bits, expected)
    scalar_type === typeof(expected) && typeof(value) === typeof(expected) ||
        throw(ArgumentError("checkpoint $name precision does not match destination"))
    expected_bits = _yl_checkpoint_bits(expected)
    typeof(bits) === typeof(expected_bits) &&
        bits === expected_bits && _yl_checkpoint_bits(value) === bits ||
        throw(ArgumentError("checkpoint $name bits do not match destination"))
    return nothing
end

function _yl_validate_checkpoint_configuration(s, st::YeeLongitudinalPIC{T}) where {T}
    s.integrator === :yee_longitudinal_pic && s.schema_version === 1 ||
        throw(ArgumentError("unsupported Yee longitudinal checkpoint schema"))
    s.D === 1 && s.T === T ||
        throw(ArgumentError("checkpoint dimension and precision must match destination"))
    for name in (:ncell, :L, :dx, :periodic)
        value, expected = getproperty(s, name), getfield(st.g, name === :ncell ? :n : name)
        typeof(value) === typeof(expected) && isequal(value, expected) ||
            throw(ArgumentError("checkpoint grid $name does not match destination"))
    end
    typeof(s.shape) === typeof(st.shape) ||
        throw(ArgumentError("checkpoint cloud shape precision does not match destination"))
    _yl_validate_checkpoint_ieee(
        :shape_width,
        s.shape.half_width,
        s.shape_width_type,
        s.shape_width_bits,
        st.shape.half_width,
    )
    for name in (:background_charge, :left_plate_charge)
        _yl_validate_checkpoint_ieee(
            name,
            getproperty(s, name),
            getproperty(s, Symbol(name, :_type)),
            getproperty(s, Symbol(name, :_bits)),
            getfield(st, name),
        )
    end
    s.gauss_policy === :original_source_exact_rounded &&
        s.plate_policy === :isolated_left_derived_right &&
        s.event_policy === :strict_event_free &&
        s.integrator_policy === :velocity_verlet &&
        s.velocity_phase === :synchronized &&
        s.timestep_policy === :fixed_after_first_success &&
        s.arithmetic_policy === :nearest_gradual ||
        throw(ArgumentError("checkpoint Yee longitudinal source or stepping policy is unsupported"))
    (s.particle_layout === :scalar || s.particle_layout === :collection) ||
        throw(ArgumentError("checkpoint particle layout must be scalar or collection"))
    return nothing
end

function _yl_checkpoint_populations(s, ::Type{T}) where {T}
    s.species isa AbstractVector && axes(s.species, 1) == Base.OneTo(length(s.species)) ||
        throw(ArgumentError("checkpoint species must be an ordered one-based vector"))
    s.particle_layout === :scalar && length(s.species) != 1 &&
        throw(ArgumentError("scalar Yee longitudinal checkpoint requires exactly one population"))
    return ParticleSet{1,T}[_checkpoint_species(record, Val(1), T) for record in s.species]
end

function _yl_validate_checkpoint_arrays(s, st::YeeLongitudinalPIC{T}, populations) where {T}
    _validate_checkpoint_fields(s.fields, _yl_checkpoint_fields(st))
    s.fields.E isa NTuple{3,Vector{T}} && s.fields.B isa NTuple{3,Vector{T}} ||
        throw(ArgumentError("checkpoint Yee fields must use ordinary vectors at grid precision"))
    fields = YeeFields{1,T}(s.fields.E, s.fields.B)
    _yl_validate_longitudinal_fields(st, fields)
    _yee_validate_impulse(s.last_impulse, st.g)
    all(_yl_stored_zero, s.last_impulse[2]) && all(_yl_stored_zero, s.last_impulse[3]) ||
        throw(ArgumentError("checkpoint Yee longitudinal impulse must have zero transverse data"))
    s.step == 0 && !all(a -> all(iszero, a), s.last_impulse) &&
        throw(ArgumentError("checkpoint initial impulse must be zero"))
    arrays = (fields.E..., fields.B..., s.last_impulse...)
    _yee_check_independent_arrays(arrays)
    for a in arrays, ps in populations
        _mightalias_particle_storage(a, ps) &&
            throw(ArgumentError("checkpoint particle, field and impulse arrays must be independent"))
    end
    return nothing
end

function _validate_checkpoint_state(s, st::YeeLongitudinalPIC{T}) where {T}
    _validate_checkpoint_container(s, _YL_CHECKPOINT_FIELDS)
    _yl_validate_checkpoint_configuration(s, st)
    s.initialized === true && s.phase === :ready ||
        throw(ArgumentError("Yee longitudinal checkpoint must be an initialized completed boundary"))
    s.time isa T && s.step isa Int && s.last_dt isa T &&
        s.last_work isa T && s.last_energy_defect isa T ||
        throw(ArgumentError("checkpoint clock and diagnostics must have destination precision"))
    _yl_validate_phase_values(
        s.initialized,
        s.time,
        s.step,
        s.last_dt,
        s.last_work,
        s.last_energy_defect,
    )
    populations = _yl_checkpoint_populations(s, T)
    _yl_validate_populations(st, populations)
    _yl_validate_checkpoint_arrays(s, st, populations)
    s.right_plate_charge isa Rational{BigInt} &&
        s.right_plate_charge == _yl_right_plate(st, populations) ||
        throw(ArgumentError("checkpoint derived right plate charge is inconsistent with sources"))
    return nothing
end

function _yl_checkpoint_state(st::YeeLongitudinalPIC{T}, species, layout) where {T}
    populations = _yl_species(species)
    _yl_validate_populations(st, populations)
    _yl_validate_phase(st)
    st.initialized[] ||
        throw(ArgumentError("Yee longitudinal checkpoints require initialized state"))
    _yl_validate_bound(st, populations)
    _yl_validate_longitudinal_fields(st, st.fields)
    return (
        integrator = :yee_longitudinal_pic,
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
        gauss_policy = :original_source_exact_rounded,
        plate_policy = :isolated_left_derived_right,
        event_policy = :strict_event_free,
        integrator_policy = :velocity_verlet,
        velocity_phase = :synchronized,
        timestep_policy = :fixed_after_first_success,
        arithmetic_policy = :nearest_gradual,
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
        last_energy_defect = st.last_energy_defect[],
    )
end

_checkpoint_state(st::YeeLongitudinalPIC, ps::ParticleSet) =
    _yl_checkpoint_state(st, ps, :scalar)
_checkpoint_state(st::YeeLongitudinalPIC, species::AbstractVector) =
    _yl_checkpoint_state(st, species, :collection)

function _yl_restore_checkpoint!(st::YeeLongitudinalPIC{T}, species, s, layout) where {T}
    _validate_checkpoint_state(s, st)
    s.particle_layout === layout ||
        throw(ArgumentError("Yee longitudinal scalar and collection checkpoint layouts cannot mix"))
    populations = _yl_species(species)
    _yl_validate_storage(st, populations)
    length(populations) == length(s.species) ||
        throw(ArgumentError("checkpoint population count does not match destination"))
    incoming = _yl_checkpoint_populations(s, s.T)
    incoming_arrays = (s.fields.E..., s.fields.B..., s.last_impulse...)
    for (ps, record) in zip(populations, s.species)
        _validate_checkpoint_destination(ps, length(record.weight))
        for source in incoming
            _particle_sets_mightalias(ps, source) && throw(
                ArgumentError("checkpoint source and destination particles must be independent"),
            )
        end
        for a in incoming_arrays
            _mightalias_particle_storage(a, ps) &&
                throw(ArgumentError("checkpoint fields must not alias destination particles"))
        end
        for a in _yl_particle_arrays(ps)
            length(a) == length(record.weight) || a isa Vector ||
                throw(ArgumentError("only Vector checkpoint destinations may resize"))
        end
    end
    # Detach all incoming leaves together: a source field may be a permutation
    # of live field storage, and later copies must still read the saved values.
    saved = deepcopy(s)
    _validate_checkpoint_state(saved, st)
    saved_populations = _yl_checkpoint_populations(saved, saved.T)
    stage = _yl_make_stage(populations, saved_populations, T)
    # A Vector can grow only after every destination has reserved its capacity.
    # Matched strided views retain their existing shape and storage throughout.
    for (ps, record) in zip(populations, saved.species)
        count = length(record.weight)
        for a in _yl_particle_arrays(ps)
            length(a) == count || sizehint!(a, count)
        end
    end
    for (ps, record) in zip(populations, saved.species)
        _restore_checkpoint_particles!(ps, record)
    end
    _restore_checkpoint_fields!(_yl_checkpoint_fields(st), saved.fields)
    foreach(copyto!, st.last_impulse, saved.last_impulse)
    st.time[] = saved.time
    st.step[] = saved.step
    st.last_dt[] = saved.last_dt
    st.right_plate_charge[] = saved.right_plate_charge
    st.last_work[] = saved.last_work
    st.last_energy_defect[] = saved.last_energy_defect
    st._stage[] = stage
    st.initialized[] = true
    return st
end

_restore_checkpoint_state!(st::YeeLongitudinalPIC, ps::ParticleSet, s) =
    _yl_restore_checkpoint!(st, ps, s, :scalar)
_restore_checkpoint_state!(st::YeeLongitudinalPIC, species::AbstractVector, s) =
    _yl_restore_checkpoint!(st, species, s, :collection)

function save_checkpoint(path::AbstractString, st::YeeLongitudinalPIC, species)
    return _save_checkpoint(path, st, _checkpoint_state(st, species))
end

function load_checkpoint!(st::YeeLongitudinalPIC, species, path::AbstractString)
    return _restore_checkpoint_state!(st, species, deserialize(path))
end

function archive_run(
    path::AbstractString,
    st::YeeLongitudinalPIC,
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
