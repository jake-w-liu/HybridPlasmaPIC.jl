struct _YeeCaptureRecord{T}
    species::Int
    id::UInt64
    tag::UInt32
    q::T
    m::T
    weight::T
    velocity::T
    side::Int8
    elapsed::_YeeImpactQ
    external_step::Int
end
struct _YeeAbsorptionStages{T}
    transport::_YeeImpactStage{T}
    binding::_YeeLongitudinalStage{T}
end

"""
    YeeAbsorbingLongitudinalPIC(g; shape::YeeTriangularCloud,
        background_charge=zero(T), left_plate_charge=zero(T),
        max_root_refinements=precision(T))

Nonrelativistic CPU Float32/64 longitudinal PIC with zero or one certified
isolated passive center capture per fixed external step. The complete rigid
cloud transfers to its isolated plate, including particle-bound endpoint tails.
The incoming final half-kick precedes removal. Signed death-current work plus
incoming kinetic energy must pass the original-trajectory passivity certificate.
This ideal law is not a material sheath, emission, circuit or radiation model.

Original geometry/factors and rational local durations retain their canonical
owners. Public live centers are strict interior; grazing, terminal/second capture,
unresolved ordering/budget or passivity reject atomically. Particle leaves must
be ordinary Vectors: stable survivor compaction preserves objects/array identities.
Private scratch may change on error; source binding changes only on success.
Initial left-plate configuration stays in `left_plate_charge`; current isolated
left charge is `free_left_plate_charge[]`. Exact plate/wall records are separate
from particle-bound cloud tails. External mutation/reinitialization is explicit.
One instance must not be stepped or mutated concurrently. Nearest/gradual mode
and ordinary Julia writes are required; custom setter/OS protection changes are
outside this transaction scope. Local elapsed capture history is not an absolute
wall-clock timestamp. No exact discrete energy or whole-run trajectory claim.
"""
struct YeeAbsorbingLongitudinalPIC{T,S<:YeeTriangularCloud,B,P}
    g::YeeGrid{1,T}
    shape::S
    background_charge::B
    left_plate_charge::P
    max_root_refinements::Int
    fields::YeeFields{1,T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    initialized::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    right_plate_charge::Base.RefValue{_YeeImpactQ}
    last_impulse::NTuple{3,Vector{T}}
    last_work::Base.RefValue{T}
    last_energy_defect::Base.RefValue{T}
    last_physical_energy_change::Base.RefValue{T}
    last_wall_count::NTuple{2,Base.RefValue{Int}}
    cumulative_wall_count::NTuple{2,Base.RefValue{Int}}
    last_wall_impulse::NTuple{2,Base.RefValue{_YeeImpactQ}}
    cumulative_wall_impulse::NTuple{2,Base.RefValue{_YeeImpactQ}}
    last_electric_impulse::Base.RefValue{_YeeImpactQ}
    last_max_snap::Base.RefValue{_YeeImpactQ}
    last_snap_energy_bound::Base.RefValue{_YeeImpactQ}
    last_predictor_time_error::Base.RefValue{_YeeImpactQ}
    last_physical_time_error::Base.RefValue{_YeeImpactQ}
    last_source_invariant_error::Base.RefValue{_YeeImpactQ}
    last_position_error_bound::Base.RefValue{_YeeImpactQ}
    last_velocity_error_bound::Base.RefValue{_YeeImpactQ}
    _candidate::YeeFields{1,T}
    _event_fields::YeeFields{1,T}
    _impulse::NTuple{3,Vector{T}}
    _charge_workspace::YeeChargeWorkspace{T}
    _stage::Base.RefValue{_YeeImpactStage{T}}
    free_left_plate_charge::Base.RefValue{_YeeImpactQ}
    initial_right_plate_charge::Base.RefValue{_YeeImpactQ}
    last_wall_charge::NTuple{2,Base.RefValue{_YeeImpactQ}}
    cumulative_wall_charge::NTuple{2,Base.RefValue{_YeeImpactQ}}
    last_wall_kinetic::NTuple{2,Base.RefValue{_YeeImpactQ}}
    cumulative_wall_kinetic::NTuple{2,Base.RefValue{_YeeImpactQ}}
    last_wall_work::NTuple{2,Base.RefValue{_YeeImpactQ}}
    cumulative_wall_work::NTuple{2,Base.RefValue{_YeeImpactQ}}
    last_wall_energy::NTuple{2,Base.RefValue{_YeeImpactQ}}
    cumulative_wall_energy::NTuple{2,Base.RefValue{_YeeImpactQ}}
    last_transport_work::Base.RefValue{_YeeImpactQ}
    last_passivity_error::Base.RefValue{_YeeImpactQ}
    last_field_energy_defect::Base.RefValue{_YeeImpactQ}
    last_total_energy_defect::Base.RefValue{T}
    last_capture::Base.RefValue{Union{Nothing,_YeeCaptureRecord{T}}}
    _initial_counts::Base.RefValue{Vector{Int}}
    _initial_charges::Base.RefValue{Vector{_YeeImpactQ}}
    _binding::Base.RefValue{_YeeLongitudinalStage{T}}
    _death_fields::YeeFields{1,T}
end

function _ya_make_stages(populations, sources, ::Type{T}, n::Int) where {T}
    return _YeeAbsorptionStages{T}(
        _yi_make_stage(populations, sources, T, n),
        _yl_make_stage(populations, sources, T),
    )
end
function YeeAbsorbingLongitudinalPIC(
    g::YeeGrid{1,T};
    shape::YeeTriangularCloud,
    background_charge = zero(T),
    left_plate_charge = zero(T),
    max_root_refinements = precision(T),
) where {T}
    g.periodic[1] && throw(ArgumentError("absorbing longitudinal PIC requires closed PEC"))
    cap = _require_positive_intlike("max_root_refinements", max_root_refinements)
    b = _yl_source_scalar("background_charge", background_charge, T)
    sigma = _yl_source_scalar("left_plate_charge", left_plate_charge, T)
    fields = YeeFields(g)
    impulse() = ntuple(c -> zeros(T, length(fields.E[c])), 3)
    qref() = Ref(zero(_YeeImpactQ))
    qpair() = ntuple(_ -> qref(), 2)
    empty = ParticleSet{1,T}[]
    prepared = _ya_make_stages(empty, empty, T, g.n[1])
    stage, binding = prepared.transport, prepared.binding
    right = _yee_closed_right_plate(sigma, b, g.L[1], empty)
    return YeeAbsorbingLongitudinalPIC{T,typeof(shape),typeof(b),typeof(sigma)}(
        g, # g
        shape, # shape
        b, # background_charge
        sigma, # left_plate_charge
        cap, # max_root_refinements
        fields, # fields
        Ref(zero(T)), # time
        Ref(0), # step
        Ref(false), # initialized
        Ref(zero(T)), # last_dt
        Ref(right), # right_plate_charge
        impulse(), # last_impulse
        Ref(zero(T)), # last_work
        Ref(zero(T)), # last_energy_defect
        Ref(zero(T)), # last_physical_energy_change
        ntuple(_ -> Ref(0), 2), # last_wall_count
        ntuple(_ -> Ref(0), 2), # cumulative_wall_count
        qpair(), # last_wall_impulse
        qpair(), # cumulative_wall_impulse
        qref(), # last_electric_impulse
        qref(), # last_max_snap
        qref(), # last_snap_energy_bound
        qref(), # last_predictor_time_error
        qref(), # last_physical_time_error
        qref(), # last_source_invariant_error
        qref(), # last_position_error_bound
        qref(), # last_velocity_error_bound
        YeeFields(g), # _candidate
        YeeFields(g), # _event_fields
        impulse(), # _impulse
        YeeChargeWorkspace(g), # _charge_workspace
        Ref(stage), # _stage
        Ref(_exact_rational(sigma)), # free_left_plate_charge
        Ref(right), # initial_right_plate_charge
        qpair(), # last_wall_charge
        qpair(), # cumulative_wall_charge
        qpair(), # last_wall_kinetic
        qpair(), # cumulative_wall_kinetic
        qpair(), # last_wall_work
        qpair(), # cumulative_wall_work
        qpair(), # last_wall_energy
        qpair(), # cumulative_wall_energy
        qref(), # last_transport_work
        qref(), # last_passivity_error
        qref(), # last_field_energy_defect
        Ref(zero(T)), # last_total_energy_defect
        Ref{Union{Nothing,_YeeCaptureRecord{T}}}(nothing), # last_capture
        Ref(Int[]), # _initial_counts
        Ref(_YeeImpactQ[]), # _initial_charges
        Ref(binding), # _binding
        YeeFields(g), # _death_fields
    )
end
YeeAbsorbingLongitudinalPIC(g::YeeGrid; kwargs...) =
    throw(ArgumentError("absorbing longitudinal PIC requires one spatial dimension"))

_yl_validate_longitudinal_fields(st::YeeAbsorbingLongitudinalPIC, fields::YeeFields) =
    _yee_closed_validate_longitudinal_fields(st, fields)
function _ya_validate_particle_vectors(st::YeeAbsorbingLongitudinalPIC{T}, species) where {T}
    all(ps -> ps isa _ElectrostaticHostParticle{1,T}, species) ||
        throw(ArgumentError("absorbing longitudinal particles require ordinary Vector leaves"))
    return nothing
end
function _yl_validate_populations(st::YeeAbsorbingLongitudinalPIC, species)
    _ya_validate_particle_vectors(st, species)
    return _yee_closed_validate_longitudinal_populations(st, species)
end
_yl_right_plate(st::YeeAbsorbingLongitudinalPIC, species) =
    _yee_closed_right_plate(st.free_left_plate_charge[], st.background_charge, st.g.L[1], species)
_yee_closed_plate_values(st::YeeAbsorbingLongitudinalPIC) =
    (st.free_left_plate_charge[], st.right_plate_charge[])
_ya_validate_bound(st, species) = _yee_closed_validate_binding(st, species, st._binding[])

function _yee_closed_references(st::YeeAbsorbingLongitudinalPIC)
    common = (
        st.time,
        st.step,
        st.initialized,
        st.last_dt,
        st.right_plate_charge,
        st.last_work,
        st.last_energy_defect,
        st.last_physical_energy_change,
        st.last_wall_count...,
        st.cumulative_wall_count...,
        st.last_wall_impulse...,
        st.cumulative_wall_impulse...,
        st.last_electric_impulse,
        st.last_max_snap,
        st.last_snap_energy_bound,
        st.last_predictor_time_error,
        st.last_physical_time_error,
        st.last_source_invariant_error,
        st.last_position_error_bound,
        st.last_velocity_error_bound,
        st._stage,
        st.free_left_plate_charge,
        st.initial_right_plate_charge,
        st.last_transport_work,
        st.last_passivity_error,
        st.last_field_energy_defect,
        st.last_total_energy_defect,
        st.last_capture,
        st._initial_counts,
        st._initial_charges,
        st._binding,
    )
    return (
        common...,
        st.last_wall_charge...,
        st.cumulative_wall_charge...,
        st.last_wall_kinetic...,
        st.cumulative_wall_kinetic...,
        st.last_wall_work...,
        st.cumulative_wall_work...,
        st.last_wall_energy...,
        st.cumulative_wall_energy...,
    )
end
function _yl_storage_arrays(st::YeeAbsorbingLongitudinalPIC)
    # The shared inventory owns transport leaves. The binding's owner/component
    # references intentionally point at public objects; only its frozen copies
    # are independent writable leaves and belong in this alias inventory.
    _yee_check_shape(st._death_fields.E, st.g, true)
    _yee_check_shape(st._death_fields.B, st.g, false)
    arrays = _yee_impact_storage_arrays(st)
    append!(
        arrays,
        (st._death_fields.E..., st._death_fields.B..., st._initial_counts[], st._initial_charges[]),
    )
    for ps in st._binding[].candidates
        _yl_check_particle_storage(ps)
        append!(arrays, _yl_particle_arrays(ps))
    end
    return arrays
end

function _ya_history_values(st)
    return merge(
        _yi_history_values(st),
        (
            right_plate_charge = st.right_plate_charge[],
            free_left_plate_charge = st.free_left_plate_charge[],
            initial_right_plate_charge = st.initial_right_plate_charge[],
            initial_counts = st._initial_counts[],
            initial_charges = st._initial_charges[],
            last_wall_charge = map(x -> x[], st.last_wall_charge),
            cumulative_wall_charge = map(x -> x[], st.cumulative_wall_charge),
            last_wall_kinetic = map(x -> x[], st.last_wall_kinetic),
            cumulative_wall_kinetic = map(x -> x[], st.cumulative_wall_kinetic),
            last_wall_work = map(x -> x[], st.last_wall_work),
            cumulative_wall_work = map(x -> x[], st.cumulative_wall_work),
            last_wall_energy = map(x -> x[], st.last_wall_energy),
            cumulative_wall_energy = map(x -> x[], st.cumulative_wall_energy),
            last_transport_work = st.last_transport_work[],
            last_passivity_error = st.last_passivity_error[],
            last_field_energy_defect = st.last_field_energy_defect[],
            last_total_energy_defect = st.last_total_energy_defect[],
            last_capture = st.last_capture[],
        ),
    )
end

function _ya_validate_history_values(s, st::YeeAbsorbingLongitudinalPIC{T}, species) where {T}
    Q, z = _exact_rational, zero(_YeeImpactQ)
    _yl_validate_phase_values(
        s.initialized,
        s.time,
        s.step,
        s.last_dt,
        s.last_work,
        s.last_energy_defect,
    )
    all(
        x -> x isa T && isfinite(x),
        (s.time, s.last_physical_energy_change, s.last_total_energy_defect),
    ) || throw(ArgumentError("invalid absorbing energy/clock history"))
    _yee_validate_impulse(s.last_impulse, st.g)
    all(a -> all(_yl_stored_zero, a), (s.last_impulse[2], s.last_impulse[3])) ||
        throw(ArgumentError("absorbing transverse impulse must be stored zero"))
    qhistory = _yi_q_history_values(s)
    for (name, value) in zip(
        (
            _YI_Q_HISTORY...,
            :right_plate_charge,
            :free_left_plate_charge,
            :initial_right_plate_charge,
            :last_transport_work,
            :last_passivity_error,
            :last_field_energy_defect,
        ),
        (
            qhistory...,
            s.right_plate_charge,
            s.free_left_plate_charge,
            s.initial_right_plate_charge,
            s.last_transport_work,
            s.last_passivity_error,
            s.last_field_energy_defect,
        ),
    )
        _yi_finite_q(value) || throw(ArgumentError("invalid absorbing rational history $name"))
    end
    all(value -> value >= 0, (Base.tail(qhistory)..., s.last_passivity_error)) ||
        throw(ArgumentError("absorbing error bounds must be nonnegative"))
    s.initial_counts isa Vector{Int} &&
        s.initial_charges isa Vector{_YeeImpactQ} &&
        length(s.initial_counts) == length(s.initial_charges) == length(species) ||
        throw(ArgumentError("invalid absorbing initial source records"))
    initial_count, live_count, initial_charge, live_charge = BigInt(0), BigInt(0), z, z
    for (k, ps) in enumerate(species)
        count, charge = length(ps.weight), sum(w -> Q(ps.q) * Q(w), ps.weight; init = z)
        s.initial_counts[k] >= count && _yi_finite_q(s.initial_charges[k]) ||
            throw(ArgumentError("absorbing initial records precede live sources"))
        lost_charge = s.initial_charges[k] - charge
        (iszero(Q(ps.q)) ? iszero(lost_charge) : lost_charge * Q(ps.q) >= 0) &&
            (s.initial_counts[k] != count || iszero(lost_charge)) ||
            throw(ArgumentError("absorbing per-species initial charge is inconsistent"))
        initial_count += s.initial_counts[k]
        live_count += count
        initial_charge += s.initial_charges[k]
        live_charge += charge
    end
    for name in (:last_wall_count, :cumulative_wall_count)
        values = getproperty(s, name)
        values isa NTuple{2,Int} && all(>=(0), values) ||
            throw(ArgumentError("invalid absorbing wall count"))
    end
    sum(BigInt, s.last_wall_count) <= 1 &&
        sum(BigInt, s.cumulative_wall_count) == initial_count - live_count &&
        sum(BigInt, s.cumulative_wall_count) <= s.step ||
        throw(ArgumentError("absorbing loss count is inconsistent with initial/live sources"))
    for name in (
        :last_wall_impulse,
        :cumulative_wall_impulse,
        :last_wall_charge,
        :cumulative_wall_charge,
        :last_wall_kinetic,
        :cumulative_wall_kinetic,
        :last_wall_work,
        :cumulative_wall_work,
        :last_wall_energy,
        :cumulative_wall_energy,
    )
        values = getproperty(s, name)
        values isa NTuple{2,_YeeImpactQ} && all(_yi_finite_q, values) ||
            throw(ArgumentError("invalid absorbing wall ledger $name"))
    end
    for side = 1:2
        0 <= s.last_wall_count[side] <= s.cumulative_wall_count[side] ||
            throw(ArgumentError("absorbing last/cumulative counts disagree"))
        sign = side == 1 ? -1 : 1
        sign * s.last_wall_impulse[side] >= 0 &&
            sign * s.cumulative_wall_impulse[side] >= abs(s.last_wall_impulse[side]) ||
            throw(ArgumentError("absorbing wall momentum has an invalid sign/history"))
        for (lasts, cumulatives) in (
            (s.last_wall_kinetic, s.cumulative_wall_kinetic),
            (s.last_wall_energy, s.cumulative_wall_energy),
        )
            last, cumulative = lasts[side], cumulatives[side]
            0 <= last <= cumulative ||
                throw(ArgumentError("absorbing passive energy history is invalid"))
        end
        s.last_wall_energy[side] == s.last_wall_kinetic[side] + s.last_wall_work[side] &&
            s.cumulative_wall_energy[side] ==
            s.cumulative_wall_kinetic[side] + s.cumulative_wall_work[side] ||
            throw(ArgumentError("absorbing received energy double-counts or omits death work"))
        for (counts, impulse, charge, kinetic, work, energy) in (
            (
                s.last_wall_count,
                s.last_wall_impulse,
                s.last_wall_charge,
                s.last_wall_kinetic,
                s.last_wall_work,
                s.last_wall_energy,
            ),
            (
                s.cumulative_wall_count,
                s.cumulative_wall_impulse,
                s.cumulative_wall_charge,
                s.cumulative_wall_kinetic,
                s.cumulative_wall_work,
                s.cumulative_wall_energy,
            ),
        )
            counts[side] != 0 ||
                all(
                    iszero,
                    (impulse[side], charge[side], kinetic[side], work[side], energy[side]),
                ) ||
                throw(ArgumentError("zero absorbing count has nonzero wall ledger"))
        end
    end
    sum(s.cumulative_wall_charge) == initial_charge - live_charge &&
        s.free_left_plate_charge == Q(st.left_plate_charge) + s.cumulative_wall_charge[1] &&
        s.initial_right_plate_charge ==
        -Q(st.left_plate_charge) - Q(st.background_charge) * Q(st.g.L[1]) - initial_charge &&
        s.right_plate_charge == s.initial_right_plate_charge + s.cumulative_wall_charge[2] &&
        s.right_plate_charge ==
        _yee_closed_right_plate(s.free_left_plate_charge, st.background_charge, st.g.L[1], species) ||
        throw(ArgumentError("absorbing dynamic plates or signed captured charge are inconsistent"))
    isequal(s.last_work, _from_exact(T, s.last_transport_work + sum(s.last_wall_work))) ||
        throw(ArgumentError("absorbing total work omits or duplicates death work"))
    if sum(s.last_wall_count) == 0
        s.last_capture === nothing &&
            iszero(s.last_passivity_error) &&
            all(
                iszero,
                (
                    s.last_max_snap,
                    s.last_snap_energy_bound,
                    s.last_predictor_time_error,
                    s.last_physical_time_error,
                ),
            ) || throw(ArgumentError("no-capture step has capture records or event errors"))
    else
        r = s.last_capture
        r isa _YeeCaptureRecord{T} || throw(ArgumentError("invalid absorbing capture record type"))
        1 <= r.species <= length(species) &&
            r.side in (1, 2) &&
            r.external_step == s.step &&
            _yi_finite_q(r.elapsed) &&
            0 < r.elapsed < Q(s.last_dt) &&
            all(isfinite, (r.q, r.m, r.weight, r.velocity)) &&
            Q(r.m) > 0 &&
            Q(r.weight) >= 0 ||
            throw(ArgumentError("invalid absorbing local capture phase/source record"))
        side = Int(r.side)
        ps = species[r.species]
        isequal(r.q, ps.q) &&
            isequal(r.m, ps.m) &&
            s.last_wall_count[side] == 1 &&
            (side == 1 ? Q(r.velocity) < 0 : Q(r.velocity) > 0) &&
            s.last_wall_charge[side] == Q(r.q) * Q(r.weight) &&
            s.last_wall_impulse[side] == Q(r.m) * Q(r.weight) * Q(r.velocity) &&
            s.last_wall_kinetic[side] == Q(r.m) * Q(r.weight) * Q(r.velocity)^2 / 2 &&
            s.last_wall_energy[side] >= s.last_passivity_error ||
            throw(ArgumentError("absorbing capture record does not agree with passive wall ledger"))
    end
    if s.step == 0
        all(
                _yl_stored_zero,
                (
                    s.time,
                    s.last_dt,
                    s.last_work,
                    s.last_energy_defect,
                    s.last_physical_energy_change,
                    s.last_total_energy_defect,
                ),
            ) &&
            all(a -> all(_yl_stored_zero, a), s.last_impulse) &&
            all(
                iszero,
                (
                    qhistory...,
                    s.last_transport_work,
                    s.last_passivity_error,
                    s.last_field_energy_defect,
                ),
            ) || throw(ArgumentError("unadvanced absorbing histories must be zero"))
    end
    return nothing
end
_ya_validate_phase(st, species) = _ya_validate_history_values(_ya_history_values(st), st, species)

function _ya_publish_extra_history_refs!(st, s)
    _yi_publish_extra_history_refs!(st, s)
    st.free_left_plate_charge[] = s.free_left_plate_charge
    st.initial_right_plate_charge[] = s.initial_right_plate_charge
    st._initial_counts[], st._initial_charges[] = s.initial_counts, s.initial_charges
    st.last_wall_charge[1][], st.last_wall_charge[2][] = s.last_wall_charge
    st.cumulative_wall_charge[1][], st.cumulative_wall_charge[2][] = s.cumulative_wall_charge
    st.last_wall_kinetic[1][], st.last_wall_kinetic[2][] = s.last_wall_kinetic
    st.cumulative_wall_kinetic[1][], st.cumulative_wall_kinetic[2][] = s.cumulative_wall_kinetic
    st.last_wall_work[1][], st.last_wall_work[2][] = s.last_wall_work
    st.cumulative_wall_work[1][], st.cumulative_wall_work[2][] = s.cumulative_wall_work
    st.last_wall_energy[1][], st.last_wall_energy[2][] = s.last_wall_energy
    st.cumulative_wall_energy[1][], st.cumulative_wall_energy[2][] = s.cumulative_wall_energy
    st.last_transport_work[] = s.last_transport_work
    st.last_passivity_error[] = s.last_passivity_error
    st.last_field_energy_defect[] = s.last_field_energy_defect
    st.last_total_energy_defect[] = s.last_total_energy_defect
    st.last_capture[] = s.last_capture
    return nothing
end

function _ya_initial_history(st::YeeAbsorbingLongitudinalPIC{T}, species, right) where {T}
    qz, z = zero(_YeeImpactQ), zero(T)
    pair = (qz, qz)
    counts = Int[length(ps.weight) for ps in species]
    charges = _YeeImpactQ[
        sum(w -> _exact_rational(ps.q) * _exact_rational(w), ps.weight; init = qz) for ps in species
    ]
    return (
        initialized = true,
        time = z,
        step = 0,
        last_dt = z,
        last_impulse = st._impulse,
        last_work = z,
        last_energy_defect = z,
        last_physical_energy_change = z,
        last_wall_count = (0, 0),
        cumulative_wall_count = (0, 0),
        last_wall_impulse = pair,
        cumulative_wall_impulse = pair,
        last_electric_impulse = qz,
        last_max_snap = qz,
        last_snap_energy_bound = qz,
        last_predictor_time_error = qz,
        last_physical_time_error = qz,
        last_source_invariant_error = qz,
        last_position_error_bound = qz,
        last_velocity_error_bound = qz,
        right_plate_charge = right,
        free_left_plate_charge = _exact_rational(st.left_plate_charge),
        initial_right_plate_charge = right,
        initial_counts = counts,
        initial_charges = charges,
        last_wall_charge = pair,
        cumulative_wall_charge = pair,
        last_wall_kinetic = pair,
        cumulative_wall_kinetic = pair,
        last_wall_work = pair,
        cumulative_wall_work = pair,
        last_wall_energy = pair,
        cumulative_wall_energy = pair,
        last_transport_work = qz,
        last_passivity_error = qz,
        last_field_energy_defect = qz,
        last_total_energy_defect = z,
        last_capture = nothing,
    )
end
function init!(st::YeeAbsorbingLongitudinalPIC{T}, particles) where {T}
    species = _yl_species(particles)
    _yl_validate_populations(st, species)
    _yl_validate_longitudinal_fields(st, st.fields)
    prepared = _ya_make_stages(species, species, T, st.g.n[1])
    right = _yee_closed_right_plate(st.left_plate_charge, st.background_charge, st.g.L[1], species)
    _yl_initialize_fields!(st, species)
    foreach(a -> fill!(a, zero(T)), st._impulse)
    history = _ya_initial_history(st, species, right)
    _ya_validate_history_values(history, st, species)
    copyto!(st.fields.E[1], st._candidate.E[1])
    foreach(copyto!, st.last_impulse, st._impulse)
    st._stage[], st._binding[] = prepared.transport, prepared.binding
    st.right_plate_charge[] = right
    st.time[], st.step[], st.last_dt[], st.last_work[] = zero(T), 0, zero(T), zero(T)
    _ya_publish_extra_history_refs!(st, history)
    st.initialized[] = true
    return st
end

function _ya_reset_transport!(st, species)
    stage = st._stage[]
    _pic_prepare_species_stages!(stage.candidates, species)
    _pic_prepare_species_stages!(stage.event_candidates, species)
    resize!(stage.offsets, length(species) + 1)
    stage.offsets[1] = 0
    for (k, ps) in enumerate(species)
        stage.offsets[k+1] = Base.checked_add(stage.offsets[k], length(ps.weight))
    end
    count = last(stage.offsets)
    resize!(stage.factors, length(species))
    for a in (
        stage.acceleration...,
        stage.gather...,
        stage.bounds...,
        stage.residual...,
        stage.carried...,
        stage.terminal...,
    )
        resize!(a, count)
        fill!(a, zero(_YeeImpactQ))
    end
    resize!(stage.selected_side, count)
    fill!(stage.selected_side, 0)
    empty!(stage.roots)
    empty!(stage.contacts)
    resize!(stage.rho, st.g.n[1])
    resize!(stage.net_impulse, st.g.n[1])
    fill!(stage.rho, zero(_YeeImpactQ))
    fill!(stage.net_impulse, zero(_YeeImpactQ))
    return stage
end
function _yi_cluster(st::YeeAbsorbingLongitudinalPIC{T}, species, geom) where {T}
    cluster = _yee_longitudinal_cluster_candidate(st, T, species, geom)
    return cluster === nothing || cluster.group != 1 ? nothing : cluster
end
_yi_local_chunk!(
    st::YeeAbsorbingLongitudinalPIC{T},
    source,
    oldfields,
    destination,
    newfields,
    duration,
    first_acceleration,
    next_acceleration,
    constrained,
) where {T} = _yee_longitudinal_local_chunk!(
    st,
    T,
    source,
    oldfields,
    destination,
    newfields,
    duration,
    first_acceleration,
    next_acceleration,
    constrained,
)

# Align at the unique true PREDEATH hit t in [lower,upper]. The unchanged
# scalar X has a two-case proof: for t>=s use prehit motion; for t<s shift only
# the predictor polynomial. No smooth physical trajectory is extended past t.
function _ya_incoming_seed!(st, species, cluster, geom, captured)
    stage = st._stage[]
    xerr, verr, aerr = stage.carried
    dx, kappa = stage.residual
    A, V, J = stage.bounds
    for i in eachindex(xerr)
        xerr[i] =
            i == captured ? zero(geom.tau) :
            abs(dx[i]) + J[i] * cluster.duration^3 / 6 + V[i] * cluster.delta
    end
    weighted = _yi_weighted_error(stage, species, xerr)
    rho = maximum(abs, stage.rho)
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        aerr[i] = _yi_force_error(stage, i, k, weighted, xerr, rho, geom)
        verr[i] =
            abs(kappa[i]) +
            cluster.delta * abs(stage.acceleration[1][i] + stage.acceleration[2][i]) / 2 +
            cluster.upper * aerr[i] / 2 +
            J[i] * cluster.upper^2 / 4
    end
    return rho + geom.M * weighted
end

function _ya_death!(st::YeeAbsorbingLongitudinalPIC{T}, charge, side, geom) where {T}
    Q, stage = _exact_rational, st._stage[]
    raw = _yee_cloud_death_sums(charge, side, st.g, st.shape)
    _yee_charge_store_exact_candidate!(st._charge_workspace, raw, st.g.n[1])
    copyto!(st._impulse[1], 1, st._charge_workspace.candidate, 1, st.g.n[1])
    for c = 2:3
        fill!(st._impulse[c], zero(T))
        fill!(st._death_fields.E[c], zero(T))
    end
    foreach(a -> fill!(a, zero(T)), st._death_fields.B)
    idealwork, magnitude = zero(geom.tau), zero(geom.tau)
    for e in eachindex(st._death_fields.E[1])
        before = Q(st._event_fields.E[1][e])
        st._death_fields.E[1][e] = _yee_source_row_exact(before, st._impulse[1], e, T)
        stage.rho[e] += Q(st._death_fields.E[1][e]) - before + raw[e]
        stage.net_impulse[e] += Q(st._impulse[1][e])
        idealwork += raw[e] * (before - raw[e] / 2)
        magnitude += abs(raw[e])
    end
    _yl_validate_longitudinal_fields(st, st._death_fields)
    work =
        _yee_source_work_rational(st._event_fields.E, st._death_fields.E, st._impulse, st.g.dx[1])
    return (work = work, ideal = geom.delta * idealwork, magnitude = magnitude)
end

function _ya_compact!(st, species_index, slot, captured, geom)
    stage = st._stage[]
    ps = stage.event_candidates[species_index]
    for a in _yl_particle_arrays(ps)
        deleteat!(a, slot)
    end
    for a in (
        stage.acceleration...,
        stage.gather...,
        stage.bounds...,
        stage.residual...,
        stage.carried...,
        stage.terminal...,
        stage.selected_side,
    )
        deleteat!(a, captured)
    end
    stage.offsets[1] = 0
    for (k, population) in enumerate(stage.event_candidates)
        stage.offsets[k+1] = Base.checked_add(stage.offsets[k], length(population.weight))
    end
    fill!(stage.selected_side, 0)
    _pic_prepare_species_stages!(stage.candidates, stage.event_candidates)
    weighted = _yi_weighted_error(stage, stage.event_candidates, stage.carried[1])
    rho = maximum(abs, stage.rho)
    for (k, population) in enumerate(stage.event_candidates), p in eachindex(population.weight)
        i = stage.offsets[k] + p
        field = _yee_cloud_field_rational(
            st._death_fields.E[1],
            population.x[1][p],
            nothing,
            st.g,
            st.shape,
        )
        stage.gather[2][i] = field
        stage.acceleration[2][i] = stage.factors[k] * field
        stage.carried[3][i] = _yi_force_error(stage, i, k, weighted, stage.carried[1], rho, geom)
    end
    return stage.event_candidates
end

function step!(st::YeeAbsorbingLongitudinalPIC{T}, particles, dt::Real) where {T}
    _yl_require_mode()
    species = _yl_species(particles)
    tau = _require_finite_positive_real("dt", dt, T)
    st.initialized[] || throw(ArgumentError("call init! before absorbing stepping"))
    _yl_validate_populations(st, species)
    _ya_validate_phase(st, species)
    _ya_validate_bound(st, species)
    _yl_validate_longitudinal_fields(st, st.fields)
    st.step[] == 0 ||
        isequal(tau, st.last_dt[]) ||
        throw(ArgumentError("absorbing dt is fixed after first success"))
    next_time, next_step = _validated_step_clock(st.time[], st.step[], tau)
    next_time > st.time[] || throw(ArgumentError("absorbing time must advance finitely"))
    oldkinetic = _kinetic_energy_checked(species, T, Val(false))
    isfinite(oldkinetic) || throw(ArgumentError("absorbing old kinetic diagnostic is nonfinite"))
    oldK = _kinetic_energy_rational(species)
    oldU = _field_energy_rational(st.fields.E, st.fields.B, one(_YeeImpactQ), st.g.dx)
    stage = _ya_reset_transport!(st, species)
    geom = _yi_global_bounds!(st, species, _exact_rational(tau))
    _yi_first_brackets!(st, species, geom)
    Q, qz = _exact_rational, zero(_YeeImpactQ)
    pair = (qz, qz)
    counts, walls, charges, kinetic, deathworks, received = (0, 0), pair, pair, pair, pair, pair
    record, cluster = nothing, nothing
    passivity = qz
    left, right = st.free_left_plate_charge[], st.right_plate_charge[]
    if isempty(stage.roots)
        for i in eachindex(stage.selected_side), side = 1:2
            c = _yi_coefficients(st, species, i, side, stage.acceleration[1])
            _yi_minimum(c, qz, geom.tau) > stage.bounds[3][i] * geom.tau^3 / 6 ||
                throw(ArgumentError("absorbing physical no-contact certificate is inconclusive"))
        end
        local_step = _yi_local_chunk!(
            st,
            species,
            st.fields,
            stage.candidates,
            st._candidate,
            geom.tau,
            1,
            2,
            false,
        )
        transport, electric = local_step.work, local_step.electric
        terminal_sources = species
    else
        cluster = _yi_certify_cluster!(st, species, geom)
        captured = findfirst(!=(0), stage.selected_side)
        captured === nothing && throw(ArgumentError("absorbing contact selection is empty"))
        side = Int(stage.selected_side[captured])
        source, slot, k = _yi_particle(stage, species, captured)
        local_step = _yi_local_chunk!(
            st,
            species,
            st.fields,
            stage.event_candidates,
            st._event_fields,
            cluster.duration,
            1,
            2,
            true,
        )
        incoming = stage.event_candidates[k].v[1][slot]
        (side == 1 ? -Q(incoming) : Q(incoming)) >= stage.terminal[2][captured] / 2 ||
            throw(ArgumentError("completed capture kick is not sufficiently incoming"))
        eE = _ya_incoming_seed!(st, species, cluster, geom, captured)
        ev = stage.carried[2][captured]
        mass = Q(source.m) * Q(source.weight[slot])
        charge = Q(source.q) * Q(source.weight[slot])
        momentum, Kcap = mass * Q(incoming), mass * Q(incoming)^2 / 2
        death = _ya_death!(st, charge, side, geom)
        energy = Kcap + death.work
        energy >= 0 || throw(ArgumentError("capture has negative represented received wall energy"))
        passivity =
            mass * (abs(Q(incoming)) * ev + ev^2 / 2) +
            geom.delta * death.magnitude * eE +
            abs(death.ideal - death.work)
        energy >= passivity ||
            throw(ArgumentError("physical passive-capture certificate is inconclusive"))
        record = _YeeCaptureRecord{T}(
            k,
            source.id[slot],
            source.tag[slot],
            source.q,
            source.m,
            source.weight[slot],
            incoming,
            Int8(side),
            cluster.duration,
            next_step,
        )
        counts = Base.setindex(counts, 1, side)
        walls, charges = Base.setindex(walls, momentum, side), Base.setindex(charges, charge, side)
        kinetic, deathworks =
            Base.setindex(kinetic, Kcap, side), Base.setindex(deathworks, death.work, side)
        received = Base.setindex(received, energy, side)
        side == 1 ? (left += charge) : (right += charge)
        terminal_sources = _ya_compact!(st, k, slot, captured, geom)
        # All selected flags are now zero: every survivor uses strict both-wall
        # carried clearance; no reflection/departing-wall exception applies.
        remainder = _yi_clear_remainder!(st, cluster, geom)
        final_step = _yi_local_chunk!(
            st,
            terminal_sources,
            st._death_fields,
            stage.candidates,
            st._candidate,
            remainder,
            2,
            3,
            false,
        )
        transport = local_step.work + final_step.work
        electric = local_step.electric + final_step.electric
    end
    _yi_terminal!(st, terminal_sources, cluster, geom)
    newkinetic = kinetic_energy(stage.candidates)
    newK = _kinetic_energy_rational(stage.candidates)
    newU = _field_energy_rational(st._candidate.E, st._candidate.B, one(_YeeImpactQ), st.g.dx)
    work = transport + sum(deathworks)
    defect = _from_exact(T, newK - oldK + sum(kinetic) - transport)
    physical = _from_exact(T, (newK - oldK) + (newU - oldU))
    totaldefect = _from_exact(T, (newK - oldK) + (newU - oldU) + sum(received))
    lastwork = _from_exact(T, work)
    all(isfinite, (newkinetic, defect, physical, totaldefect, lastwork)) ||
        throw(ArgumentError("absorbing required energy/work diagnostics do not round finite"))
    for i in eachindex(stage.net_impulse)
        value = _from_exact(T, stage.net_impulse[i])
        isfinite(value) || throw(ArgumentError("net absorbing impulse does not round finite"))
        st._impulse[1][i] = value
    end
    foreach(a -> fill!(a, zero(T)), (st._impulse[2], st._impulse[3]))
    cumulative_count =
        ntuple(side -> Base.checked_add(st.cumulative_wall_count[side][], counts[side]), 2)
    cumulative_wall = ntuple(side -> st.cumulative_wall_impulse[side][] + walls[side], 2)
    cumulative_charge = ntuple(side -> st.cumulative_wall_charge[side][] + charges[side], 2)
    cumulative_kinetic = ntuple(side -> st.cumulative_wall_kinetic[side][] + kinetic[side], 2)
    cumulative_work = ntuple(side -> st.cumulative_wall_work[side][] + deathworks[side], 2)
    cumulative_energy = ntuple(side -> st.cumulative_wall_energy[side][] + received[side], 2)
    history = (
        initialized = true,
        time = next_time,
        step = next_step,
        last_dt = tau,
        last_impulse = st._impulse,
        last_work = lastwork,
        last_energy_defect = defect,
        last_physical_energy_change = physical,
        last_wall_count = counts,
        cumulative_wall_count = cumulative_count,
        last_wall_impulse = walls,
        cumulative_wall_impulse = cumulative_wall,
        last_electric_impulse = electric,
        last_max_snap = cluster === nothing ? qz : cluster.max_snap,
        last_snap_energy_bound = cluster === nothing ? qz :
                                 geom.delta * (
            geom.field * cluster.charge_snap / geom.h +
            geom.M * cluster.charge_snap^2 / (2 * geom.h)
        ),
        last_predictor_time_error = cluster === nothing ? qz : cluster.predictor_error,
        last_physical_time_error = cluster === nothing ? qz : cluster.delta,
        last_source_invariant_error = maximum(abs, stage.rho),
        last_position_error_bound = maximum(stage.terminal[1]; init = qz),
        last_velocity_error_bound = maximum(stage.terminal[2]; init = qz),
        right_plate_charge = right,
        free_left_plate_charge = left,
        initial_right_plate_charge = st.initial_right_plate_charge[],
        initial_counts = st._initial_counts[],
        initial_charges = st._initial_charges[],
        last_wall_charge = charges,
        cumulative_wall_charge = cumulative_charge,
        last_wall_kinetic = kinetic,
        cumulative_wall_kinetic = cumulative_kinetic,
        last_wall_work = deathworks,
        cumulative_wall_work = cumulative_work,
        last_wall_energy = received,
        cumulative_wall_energy = cumulative_energy,
        last_transport_work = transport,
        last_passivity_error = passivity,
        last_field_energy_defect = newU - oldU + work,
        last_total_energy_defect = totaldefect,
        last_capture = record,
    )
    _ya_validate_history_values(history, st, stage.candidates)
    # A no-capture step changes x/v only: its authoritative source snapshot is
    # still valid, so do not rematerialize it. A loss needs a new frozen binding.
    binding = cluster === nothing ? st._binding[] : _yl_make_stage(species, stage.candidates, T)
    owners = _ElectrostaticHostParticle{1,T}[ps for ps in species]
    # All validation, allocation, counts and conversions precede public stores.
    # Equal logical population count and shrink-only leaves make this existing
    # materializer an ordinary typed copy/shrink operation, with no grow branch.
    _pic_prepare_species_stages!(owners, stage.candidates)
    copyto!(st.fields.E[1], st._candidate.E[1])
    foreach(copyto!, st.last_impulse, st._impulse)
    st.time[], st.step[], st.last_dt[], st.last_work[] = next_time, next_step, tau, lastwork
    st.right_plate_charge[], st._binding[] = right, binding
    _ya_publish_extra_history_refs!(st, history)
    return st
end
