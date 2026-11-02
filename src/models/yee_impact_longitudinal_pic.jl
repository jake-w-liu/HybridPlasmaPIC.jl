const _YeeImpactQ = Rational{BigInt}

struct _YeeImpactBracket
    particle::Int
    side::Int8
    lo::_YeeImpactQ
    hi::_YeeImpactQ
    exact::Bool
    refinements::Int
end
struct _YeeImpactContact
    bracket::Int
    lo::_YeeImpactQ
    hi::_YeeImpactQ
    incoming::_YeeImpactQ
end
struct _YeeImpactStage{T}
    owners::Vector{ParticleSet{1,T}}
    components::Vector{Tuple}
    candidates::Vector{_ElectrostaticHostParticle{1,T}}
    event_candidates::Vector{_ElectrostaticHostParticle{1,T}}
    offsets::Vector{Int}
    factors::Vector{_YeeImpactQ}
    acceleration::NTuple{3,Vector{_YeeImpactQ}}
    gather::NTuple{3,Vector{_YeeImpactQ}}
    bounds::NTuple{3,Vector{_YeeImpactQ}}
    residual::NTuple{2,Vector{_YeeImpactQ}}
    carried::NTuple{3,Vector{_YeeImpactQ}}
    terminal::NTuple{2,Vector{_YeeImpactQ}}
    selected_side::Vector{Int8}
    roots::Vector{_YeeImpactBracket}
    contacts::Vector{_YeeImpactContact}
    rho::Vector{_YeeImpactQ}
    net_impulse::Vector{_YeeImpactQ}
end

"""
    YeeImpactLongitudinalPIC(g; shape::YeeTriangularCloud,
        background_charge=zero(T), left_plate_charge=zero(T),
        max_root_refinements=precision(T))

Nonrelativistic CPU Float32/64 longitudinal cloud PIC with elastic center-wall
reflection. At most one certified contact cluster and two synchronized local
KDK chunks are admitted per external step. All completed public centers are
strict interior; initial wall contact, grazing, terminal contact, unresolved
ordering, a second cluster or exhausted refinement reject atomically.

Internal durations, gathers and factors remain rational through the canonical
Boris map. Only final half-kick velocities, coordinates, current and fields
round. Actual constrained endpoints have a four-ulp operation-scale snap gate;
carried physical certificates retain event-time/state/source uncertainty from
this step's original stored input. They are not a whole-run accuracy theorem.
The root-refinement cap is a finite work budget, not a success tolerance.

Call `init!`, then fixed-positive-dt `step!`. Original source configuration,
counts/provenance and storage identities remain fixed until init/restore.
Nearest rounding, gradual underflow and structurally zero transverse inputs
are required. `last_energy_defect` is once-rounded complete rational DeltaK
minus accepted local source work; it is generally nonzero and differs from
`last_physical_energy_change`. Wall charge/work are zero by the fixed elastic
law; exact represented wall momentum is delivered to each wall.

Public array/Ref identities are retained. Private scratch may change on error.
Storage must support ordinary Julia writes and alias detection; special-store
ancestry rejects. OS protection changes, arbitrary custom setter failures and
concurrent stepping/mutation are outside the transaction guarantee.
"""
struct YeeImpactLongitudinalPIC{T,S<:YeeTriangularCloud,B,P}
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
end

function _yi_make_stage(populations, sources, ::Type{T}, n::Int) where {T}
    final = _yl_make_stage(populations, sources, T)
    event = _ElectrostaticHostParticle{1,T}[]
    _pic_prepare_species_stages!(event, sources)
    offsets = Int[0]
    for ps in sources
        push!(offsets, Base.checked_add(last(offsets), length(ps.weight)))
    end
    count = last(offsets)
    capacity = Base.checked_mul(2, count)
    qvector() = fill(zero(_YeeImpactQ), count)
    roots, contacts = _YeeImpactBracket[], _YeeImpactContact[]
    sizehint!(roots, capacity)
    sizehint!(contacts, capacity)
    return _YeeImpactStage{T}(final.owners, final.components, final.candidates, event,
        offsets, fill(zero(_YeeImpactQ), length(sources)), ntuple(_ -> qvector(), 3), ntuple(_ -> qvector(), 3),
        ntuple(_ -> qvector(), 3), ntuple(_ -> qvector(), 2), ntuple(_ -> qvector(), 3),
        ntuple(_ -> qvector(), 2), zeros(Int8, count), roots, contacts,
        fill(zero(_YeeImpactQ), n), fill(zero(_YeeImpactQ), n))
end

function YeeImpactLongitudinalPIC(g::YeeGrid{1,T}; shape::YeeTriangularCloud,
    background_charge = zero(T), left_plate_charge = zero(T),
    max_root_refinements = precision(T)) where {T}
    g.periodic[1] && throw(ArgumentError("impact longitudinal PIC requires closed PEC"))
    cap = _require_positive_intlike("max_root_refinements", max_root_refinements)
    b = _yl_source_scalar("background_charge", background_charge, T)
    sigma = _yl_source_scalar("left_plate_charge", left_plate_charge, T)
    fields, final, event = YeeFields(g), YeeFields(g), YeeFields(g)
    impulse() = ntuple(c -> zeros(T, length(fields.E[c])), 3)
    qref() = Ref(zero(_YeeImpactQ))
    stage = _yi_make_stage(ParticleSet{1,T}[], ParticleSet{1,T}[], T, g.n[1])
    return YeeImpactLongitudinalPIC{T,typeof(shape),typeof(b),typeof(sigma)}(
        g, shape, b, sigma, cap, fields, Ref(zero(T)), Ref(0), Ref(false), Ref(zero(T)),
        Ref(-_exact_rational(sigma) - _exact_rational(b) * _exact_rational(g.L[1])),
        impulse(), Ref(zero(T)), Ref(zero(T)), Ref(zero(T)),
        ntuple(_ -> Ref(0), 2), ntuple(_ -> Ref(0), 2), ntuple(_ -> qref(), 2),
        ntuple(_ -> qref(), 2), qref(), qref(), qref(), qref(), qref(), qref(), qref(), qref(),
        final, event, impulse(), YeeChargeWorkspace(g), Ref(stage))
end
YeeImpactLongitudinalPIC(g::YeeGrid; kwargs...) =
    throw(ArgumentError("impact longitudinal PIC requires one spatial dimension"))
_yl_validate_longitudinal_fields(st::YeeImpactLongitudinalPIC, fields::YeeFields) =
    _yee_closed_validate_longitudinal_fields(st, fields)
_yl_validate_populations(st::YeeImpactLongitudinalPIC, species) =
    _yee_closed_validate_longitudinal_populations(st, species)

function _yee_impact_storage_arrays(st)
    _yee_check_shape(st._event_fields.E, st.g, true)
    _yee_check_shape(st._event_fields.B, st.g, false)
    stage = st._stage[]
    arrays = Any[st.fields.E..., st.fields.B..., st._candidate.E..., st._candidate.B...,
        st._event_fields.E..., st._event_fields.B..., st.last_impulse..., st._impulse...,
        st._charge_workspace.sum, st._charge_workspace.candidate,
        stage.offsets, stage.factors, stage.acceleration..., stage.gather..., stage.bounds...,
        stage.residual..., stage.carried..., stage.terminal..., stage.selected_side,
        stage.roots, stage.contacts, stage.rho, stage.net_impulse]
    for candidates in (stage.candidates, stage.event_candidates), ps in candidates
        _yl_check_particle_storage(ps)
        append!(arrays, _yl_particle_arrays(ps))
    end
    return arrays
end
_yl_storage_arrays(st::YeeImpactLongitudinalPIC) = _yee_impact_storage_arrays(st)
_yee_closed_references(st::YeeImpactLongitudinalPIC) =
    (st.time, st.step, st.initialized, st.last_dt, st.right_plate_charge, st.last_work,
    st.last_energy_defect, st.last_physical_energy_change, st.last_wall_count...,
    st.cumulative_wall_count..., st.last_wall_impulse..., st.cumulative_wall_impulse...,
    st.last_electric_impulse, st.last_max_snap, st.last_snap_energy_bound,
    st.last_predictor_time_error, st.last_physical_time_error, st.last_source_invariant_error,
    st.last_position_error_bound, st.last_velocity_error_bound, st._stage)

function _yi_validate_bound(st, species)
    _yl_validate_bound(st, species)
    stage = st._stage[]
    length(stage.event_candidates) == length(species) == length(stage.factors) ||
        throw(ArgumentError("impact private population count changed"))
    length(stage.offsets) == length(species) + 1 && first(stage.offsets) == 0 ||
        throw(ArgumentError("impact private population offsets changed"))
    count = 0
    for (k, ps) in enumerate(species)
        _yee_closed_frozen_sources_match(ps, stage.event_candidates[k]) ||
            throw(ArgumentError("impact event candidate sources changed"))
        length(stage.event_candidates[k].weight) == length(ps.weight) ||
            throw(ArgumentError("impact event candidate count changed"))
        count = Base.checked_add(count, length(ps.weight))
        stage.offsets[k+1] == count || throw(ArgumentError("impact offsets changed"))
    end
    all(a -> length(a) == count, (stage.acceleration..., stage.gather..., stage.bounds..., stage.residual...,
        stage.carried..., stage.terminal..., stage.selected_side)) ||
        throw(DimensionMismatch("impact certificate buffer shape changed"))
    length(stage.rho) == st.g.n[1] == length(stage.net_impulse) ||
        throw(DimensionMismatch("impact field certificate shape changed"))
    return nothing
end

function _yi_history_values(st)
    return (initialized = st.initialized[], time = st.time[], step = st.step[],
        last_dt = st.last_dt[], last_impulse = st.last_impulse, last_work = st.last_work[],
        last_energy_defect = st.last_energy_defect[],
        last_physical_energy_change = st.last_physical_energy_change[],
        last_wall_count = map(x -> x[], st.last_wall_count),
        cumulative_wall_count = map(x -> x[], st.cumulative_wall_count),
        last_wall_impulse = map(x -> x[], st.last_wall_impulse),
        cumulative_wall_impulse = map(x -> x[], st.cumulative_wall_impulse),
        last_electric_impulse = st.last_electric_impulse[], last_max_snap = st.last_max_snap[],
        last_snap_energy_bound = st.last_snap_energy_bound[],
        last_predictor_time_error = st.last_predictor_time_error[],
        last_physical_time_error = st.last_physical_time_error[],
        last_source_invariant_error = st.last_source_invariant_error[],
        last_position_error_bound = st.last_position_error_bound[],
        last_velocity_error_bound = st.last_velocity_error_bound[])
end
const _YI_Q_HISTORY = (:last_electric_impulse, :last_max_snap, :last_snap_energy_bound,
    :last_predictor_time_error, :last_physical_time_error, :last_source_invariant_error,
    :last_position_error_bound, :last_velocity_error_bound)
_yi_finite_q(x) = x isa _YeeImpactQ && denominator(x) > 0
function _yi_validate_history_values(s, g::YeeGrid{1,T}, count::Int) where {T}
    _yl_validate_phase_values(s.initialized, s.time, s.step, s.last_dt, s.last_work, s.last_energy_defect)
    s.time isa T && s.last_physical_energy_change isa T && isfinite(s.last_physical_energy_change) ||
        throw(ArgumentError("invalid impact energy-change history"))
    _yee_validate_impulse(s.last_impulse, g)
    all(a -> all(_yl_stored_zero, a), (s.last_impulse[2], s.last_impulse[3])) ||
        throw(ArgumentError("impact transverse impulse history must be stored zero"))
    for name in _YI_Q_HISTORY
        value = getproperty(s, name)
        _yi_finite_q(value) && (name === :last_electric_impulse || value >= 0) ||
            throw(ArgumentError("invalid impact rational history $name"))
    end
    for name in (:last_wall_count, :cumulative_wall_count)
        values = getproperty(s, name)
        values isa NTuple{2,Int} && all(>=(0), values) ||
            throw(ArgumentError("invalid impact wall counts"))
    end
    for name in (:last_wall_impulse, :cumulative_wall_impulse)
        values = getproperty(s, name)
        values isa NTuple{2,_YeeImpactQ} && all(_yi_finite_q, values) && values[1] <= 0 <= values[2] ||
            throw(ArgumentError("invalid impact wall momentum"))
    end
    for side = 1:2
        0 <= s.last_wall_count[side] <= s.cumulative_wall_count[side] &&
            abs(s.last_wall_impulse[side]) <= abs(s.cumulative_wall_impulse[side]) ||
            throw(ArgumentError("inconsistent impact cumulative history"))
        (s.last_wall_count[side] != 0 || iszero(s.last_wall_impulse[side])) &&
            (s.cumulative_wall_count[side] != 0 || iszero(s.cumulative_wall_impulse[side])) ||
            throw(ArgumentError("zero impact count has nonzero momentum"))
    end
    sum(BigInt, s.last_wall_count) <= count &&
        sum(BigInt, s.cumulative_wall_count) <= BigInt(s.step) * count ||
        throw(ArgumentError("impact counts exceed completed step scope"))
    if sum(BigInt, s.last_wall_count) == 0
        all(name -> iszero(getproperty(s, name)), (:last_max_snap, :last_snap_energy_bound,
            :last_predictor_time_error, :last_physical_time_error)) ||
            throw(ArgumentError("no-event impact history has contact residuals"))
    end
    if s.step == 0
        all(_yl_stored_zero, (s.time, s.last_dt, s.last_work, s.last_energy_defect,
            s.last_physical_energy_change)) && all(a -> all(_yl_stored_zero, a), s.last_impulse) &&
            all(name -> iszero(getproperty(s, name)), _YI_Q_HISTORY) &&
            all(iszero, (s.last_wall_count..., s.cumulative_wall_count...,
                s.last_wall_impulse..., s.cumulative_wall_impulse...)) ||
            throw(ArgumentError("unadvanced impact histories must be zero"))
    end
    return nothing
end
function _yi_validate_phase(st)
    count = 0
    for ps in st._stage[].candidates
        count = Base.checked_add(count, length(ps.weight))
    end
    return _yi_validate_history_values(_yi_history_values(st), st.g, count)
end

function _yi_publish_extra_history_refs!(st, s)
    st.last_energy_defect[] = s.last_energy_defect
    st.last_physical_energy_change[] = s.last_physical_energy_change
    for side = 1:2
        st.last_wall_count[side][] = s.last_wall_count[side]
        st.cumulative_wall_count[side][] = s.cumulative_wall_count[side]
        st.last_wall_impulse[side][] = s.last_wall_impulse[side]
        st.cumulative_wall_impulse[side][] = s.cumulative_wall_impulse[side]
    end
    st.last_electric_impulse[] = s.last_electric_impulse
    st.last_max_snap[] = s.last_max_snap
    st.last_snap_energy_bound[] = s.last_snap_energy_bound
    st.last_predictor_time_error[] = s.last_predictor_time_error
    st.last_physical_time_error[] = s.last_physical_time_error
    st.last_source_invariant_error[] = s.last_source_invariant_error
    st.last_position_error_bound[] = s.last_position_error_bound
    st.last_velocity_error_bound[] = s.last_velocity_error_bound
    return nothing
end

function init!(st::YeeImpactLongitudinalPIC{T}, particles) where {T}
    species = _yl_species(particles)
    _yl_validate_populations(st, species)
    _yl_validate_longitudinal_fields(st, st.fields)
    stage = _yi_make_stage(species, species, T, st.g.n[1])
    right = _yl_right_plate(st, species)
    _yl_initialize_fields!(st, species)
    qz, z = zero(_YeeImpactQ), zero(T)
    history = (last_energy_defect = z, last_physical_energy_change = z,
        last_wall_count = (0, 0), cumulative_wall_count = (0, 0),
        last_wall_impulse = (qz, qz), cumulative_wall_impulse = (qz, qz),
        last_electric_impulse = qz, last_max_snap = qz, last_snap_energy_bound = qz,
        last_predictor_time_error = qz, last_physical_time_error = qz,
        last_source_invariant_error = qz, last_position_error_bound = qz,
        last_velocity_error_bound = qz)
    copyto!(st.fields.E[1], st._candidate.E[1])
    foreach(a -> fill!(a, z), st.last_impulse)
    st._stage[] = stage
    st.right_plate_charge[] = right
    st.time[], st.step[], st.last_dt[], st.last_work[] = z, 0, z, z
    _yi_publish_extra_history_refs!(st, history)
    st.initialized[] = true
    return st
end

function _yi_global_bounds!(st, species, tau)
    Q, stage = _exact_rational, st._stage[]
    L, delta = Q(st.g.L[1]), Q(st.g.dx[1])
    h, a = L / st.g.n[1], Q(st.shape.half_width)
    charge = _yee_longitudinal_rest_bound(species, L, delta, tau, h, a)
    field = maximum(x -> abs(Q(x)), st.fields.E[1]) + charge
    M = min(inv(h), inv(a))
    A, V, J = stage.bounds
    a0 = stage.acceleration[1]
    sumQV = zero(tau)
    for (k, ps) in enumerate(species)
        alpha = _boris_exact_factor(_validated_boris_factor(ps))
        stage.factors[k] = alpha
        for p in eachindex(ps.weight)
            i = stage.offsets[k] + p
            A[i] = abs(alpha) * delta * field / h
            V[i] = abs(Q(ps.v[1][p])) + A[i] * tau
            stage.gather[1][i] = _yee_cloud_field_rational(st.fields.E[1], ps.x[1][p], nothing, st.g, st.shape)
            a0[i] = alpha * stage.gather[1][i]
            sumQV += abs(Q(ps.q) * Q(ps.weight[p])) * V[i]
        end
    end
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        J[i] = abs(stage.factors[k]) * delta * (M / h * sumQV + 2field / (a * h) * V[i])
    end
    return (L = L, delta = delta, h = h, a = a, M = M, field = field, tau = tau)
end
function _yi_particle(stage, species, i)
    k = searchsortedlast(stage.offsets, i - 1)
    return species[k], i - stage.offsets[k], k
end
function _yi_coefficients(st, species, i, side, acceleration)
    ps, p, _ = _yi_particle(st._stage[], species, i)
    x, v = _exact_rational(ps.x[1][p]), _exact_rational(ps.v[1][p])
    return side == 1 ? (x, v, acceleration[i]) :
        (_exact_rational(st.g.L[1]) - x, -v, -acceleration[i])
end
_yi_polynomial(c, t) = c[1] + t * c[2] + t^2 * c[3] / 2
_yi_slope(c, t) = c[2] + t * c[3]
function _yi_minimum(c, lo, hi)
    value = min(_yi_polynomial(c, lo), _yi_polynomial(c, hi))
    if !iszero(c[3])
        stationary = -c[2] / c[3]
        lo < stationary < hi && (value = min(value, _yi_polynomial(c, stationary)))
    end
    return value
end

function _yi_first_root(c, horizon, i, side)
    z = zero(horizon)
    c[1] > 0 || throw(ArgumentError("impact initial center is not strict interior"))
    stationary = iszero(c[3]) ? z : -c[2] / c[3]
    points = 0 < stationary < horizon ? (z, stationary, horizon) : (z, horizon)
    for j = 1:length(points)-1
        lo, hi = points[j], points[j+1]
        before, after = _yi_polynomial(c, lo), _yi_polynomial(c, hi)
        if iszero(after)
            slope = _yi_slope(c, hi)
            iszero(slope) && throw(ArgumentError("grazing quadratic contact is unsupported"))
            if slope < 0
                hi < horizon || throw(ArgumentError("terminal quadratic contact is unsupported"))
                return _YeeImpactBracket(i, Int8(side), hi, hi, true, 0)
            end
        elseif before > 0 && after < 0
            return _YeeImpactBracket(i, Int8(side), lo, hi, false, 0)
        end
    end
    return nothing
end
function _yi_first_brackets!(st, species, geom)
    stage = st._stage[]
    empty!(stage.roots)
    for i in eachindex(stage.selected_side), side = 1:2
        c = _yi_coefficients(st, species, i, side, stage.acceleration[1])
        root = _yi_first_root(c, geom.tau, i, side)
        root === nothing || push!(stage.roots, root)
    end
    return nothing
end
function _yi_refine!(st, species)
    stage = st._stage[]
    for k in eachindex(stage.roots)
        r = stage.roots[k]
        r.exact && continue
        r.refinements < st.max_root_refinements || continue
        c = _yi_coefficients(st, species, r.particle, r.side, stage.acceleration[1])
        middle = (r.lo + r.hi) / 2
        value = _yi_polynomial(c, middle)
        lo, hi = value > 0 ? (middle, r.hi) : value < 0 ? (r.lo, middle) : (middle, middle)
        stage.roots[k] = _YeeImpactBracket(r.particle, r.side, lo, hi,
            iszero(value), r.refinements + 1)
    end
    return nothing
end

function _yi_physical_contact(st, species, rootindex, geom)
    stage = st._stage[]
    r = stage.roots[rootindex]
    c = _yi_coefficients(st, species, r.particle, r.side, stage.acceleration[1])
    J = stage.bounds[3][r.particle]
    H = J * geom.tau^3 / 6
    nu = -max(_yi_slope(c, r.lo), _yi_slope(c, r.hi))
    nu > 0 || return nothing
    width = 2H / nu
    lo, hi = r.lo - width, r.hi + width
    0 < lo <= hi < geom.tau || return nothing
    maximum_slope = max(_yi_slope(c, lo), _yi_slope(c, hi))
    maximum_slope <= -nu / 2 || return nothing
    exact = r.exact && iszero(H)
    if exact
        _yi_polynomial(c, lo) == 0 && _yi_minimum(c, zero(lo), lo) >= 0 || return nothing
    else
        _yi_polynomial(c, lo) > H && _yi_polynomial(c, hi) < -H &&
            _yi_minimum(c, zero(lo), lo) > H || return nothing
    end
    opposite = _yi_coefficients(st, species, r.particle, 3 - r.side, stage.acceleration[1])
    _yi_minimum(opposite, zero(lo), hi) > H || return nothing
    incoming = -maximum_slope - J * hi^2 / 2
    incoming > 0 || return nothing
    return _YeeImpactContact(rootindex, lo, hi, incoming)
end
function _yi_contacts!(st, species, geom)
    stage = st._stage[]
    empty!(stage.contacts)
    # Only a particle's first predictor branch may claim its first physical hit.
    cursor = 1
    for i in eachindex(stage.selected_side)
        first = cursor
        first > length(stage.roots) && continue
        stage.roots[first].particle == i || continue
        last = first
        while last < length(stage.roots) && stage.roots[last+1].particle == i
            last += 1
        end
        cursor = last + 1
        indexes = first:last
        chosen = first
        for k in indexes
            stage.roots[k].hi < stage.roots[chosen].hi && (chosen = k)
        end
        r = stage.roots[chosen]
        for k in indexes
            k == chosen && continue
            r.hi < stage.roots[k].lo || return false
        end
        contact = _yi_physical_contact(st, species, chosen, geom)
        contact === nothing && return false
        push!(stage.contacts, contact)
    end
    sort!(stage.contacts; by = c -> c.lo)
    return true
end

function _yee_longitudinal_local_chunk!(st, ::Type{T}, source, oldfields,
    destination, newfields, duration, first_acceleration, next_acceleration, constrained) where {T}
    stage, Q = st._stage[], _exact_rational
    half = duration / 2
    dx, kappa = stage.residual
    for (k, (ps, candidate)) in enumerate(zip(source, destination)), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        u = _yl_kick(ps.v[1][p], stage.gather[first_acceleration][i], stage.factors[k], half)
        free = _particle_coordinate_native(Q(ps.x[1][p]), Q(u), duration)
        side = constrained ? stage.selected_side[i] : 0
        y = side == 1 ? zero(T) : side == 2 ? st.g.L[1] : _from_exact(T, free)
        isfinite(y) && (side != 0 || 0 < Q(y) < Q(st.g.L[1])) ||
            throw(ArgumentError("impact local drift does not produce an admitted endpoint"))
        candidate.x[1][p], candidate.v[1][p] = y, u
        candidate.v[2][p], candidate.v[3][p] = zero(T), zero(T)
        dx[i] = Q(y) - (Q(ps.x[1][p]) + duration * Q(ps.v[1][p]) +
            duration^2 * stage.acceleration[first_acceleration][i] / 2)
    end
    old, new = [ps.x[1] for ps in source], [ps.x[1] for ps in destination]
    _yee_validate_charge(st._impulse[1], source, old, new, st.g,
        st._charge_workspace, Val(:impulse), false)
    raw = _yee_charge_exact_sums(source, old, new, st.g, Val(:impulse), st.shape)
    _yee_charge_store_exact_candidate!(st._charge_workspace, raw, st.g.n[1])
    copyto!(st._impulse[1], 1, st._charge_workspace.candidate, 1, st.g.n[1])
    for c = 2:3
        fill!(st._impulse[c], zero(T))
        fill!(newfields.E[c], zero(T))
    end
    foreach(a -> fill!(a, zero(T)), newfields.B)
    for i in eachindex(newfields.E[1])
        newfields.E[1][i] = _yee_source_row_exact(Q(oldfields.E[1][i]), st._impulse[1], i, T)
        stage.rho[i] += Q(newfields.E[1][i]) - Q(oldfields.E[1][i]) + raw[i]
        stage.net_impulse[i] += Q(st._impulse[1][i])
    end
    _yl_validate_longitudinal_fields(st, newfields)
    work = _yee_source_work_rational(oldfields.E, newfields.E, st._impulse, st.g.dx[1])
    electric = zero(duration)
    for (k, (ps, candidate)) in enumerate(zip(source, destination)), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        field = _yee_cloud_field_rational(newfields.E[1], candidate.x[1][p], nothing, st.g, st.shape)
        stage.gather[next_acceleration][i] = field
        stage.acceleration[next_acceleration][i] = stage.factors[k] * field
        candidate.v[1][p] = _yl_kick(candidate.v[1][p], field, stage.factors[k], half)
        kappa[i] = Q(candidate.v[1][p]) - Q(ps.v[1][p]) -
            half * (stage.acceleration[first_acceleration][i] + stage.acceleration[next_acceleration][i])
        electric += Q(ps.m) * Q(ps.weight[p]) * (Q(candidate.v[1][p]) - Q(ps.v[1][p]))
    end
    return (work = work, electric = electric)
end
_yi_local_chunk!(st::YeeImpactLongitudinalPIC{T}, source, oldfields,
    destination, newfields, duration, first_acceleration, next_acceleration, constrained) where {T} =
    _yee_longitudinal_local_chunk!(st, T, source, oldfields, destination, newfields,
        duration, first_acceleration, next_acceleration, constrained)

function _yi_weighted_error(stage, species, errors)
    total = zero(_YeeImpactQ)
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        total += abs(_exact_rational(ps.q) * _exact_rational(ps.weight[p])) * errors[stage.offsets[k]+p]
    end
    return total
end
function _yi_force_error(stage, i, k, weighted, errors, rho, geom)
    return abs(stage.factors[k]) * geom.delta *
        (rho / geom.h + geom.M / geom.h * weighted +
            2 * geom.field / (geom.a * geom.h) * errors[i])
end
function _yi_carry!(st, species, cluster, geom)
    stage = st._stage[]
    A, V, J = stage.bounds
    xerr, verr, aerr = stage.carried
    dx, kappa = stage.residual
    for i in eachindex(xerr)
        selected = stage.selected_side[i] != 0
        if selected
            incoming = stage.terminal[2][i]
            A[i] * cluster.width < incoming / 2 && V[i] * cluster.width < geom.L ||
                throw(ArgumentError("physical contact hull permits another wall event"))
        end
        xerr[i] = selected ? V[i] * cluster.width :
            abs(dx[i]) + J[i] * cluster.duration^3 / 6 + V[i] * cluster.delta
    end
    weighted = _yi_weighted_error(stage, species, xerr)
    rho = maximum(abs, stage.rho)
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        aerr[i] = _yi_force_error(stage, i, k, weighted, xerr, rho, geom)
        verr[i] = abs(kappa[i]) + cluster.delta * abs(stage.acceleration[1][i] +
            stage.acceleration[2][i]) / 2 + cluster.upper * aerr[i] / 2 + J[i] * cluster.upper^2 / 4
        stage.selected_side[i] != 0 && (verr[i] += A[i] * cluster.width)
    end
    return nothing
end

function _yi_clear_remainder!(st, cluster, geom)
    stage = st._stage[]
    particles = stage.event_candidates
    duration, horizon = geom.tau - cluster.duration, geom.tau - cluster.lower
    A, _, J = stage.bounds
    xerr, verr, aerr = stage.carried
    for i in eachindex(stage.selected_side)
        selected = stage.selected_side[i]
        H = xerr[i] + horizon * verr[i] + horizon^2 * aerr[i] / 2 + J[i] * horizon^3 / 6
        for side = 1:2
            c = _yi_coefficients(st, particles, i, side, stage.acceleration[2])
            if selected == side
                # Exempt only the certified inward departure at this same wall.
                iszero(c[1]) && c[2] > 0 && c[2] + duration * c[3] / 2 > 0 ||
                    throw(ArgumentError("numerical remainder has same-wall return/contact"))
                inward = max(stage.terminal[2][i] - A[i] * cluster.width, c[2] - verr[i])
                inward > 0 && inward + (c[3] - aerr[i]) * horizon / 2 - J[i] * horizon^2 / 6 > 0 ||
                    throw(ArgumentError("carried physical same-wall return cannot be excluded"))
            else
                _yi_minimum(c, zero(duration), duration) > 0 ||
                    throw(ArgumentError("numerical remainder contains another wall event"))
                _yi_minimum(c, zero(duration), horizon) > H ||
                    throw(ArgumentError("carried physical remainder clearance is inconclusive"))
            end
        end
    end
    return duration
end

function _yi_terminal!(st, species, cluster, geom)
    stage = st._stage[]
    X, V = stage.terminal
    J = stage.bounds[3]
    dx, kappa = stage.residual
    event = cluster !== nothing
    horizon = event ? geom.tau - cluster.lower : geom.tau
    lower = event ? geom.tau - cluster.upper : geom.tau
    first = event ? 2 : 1
    next = event ? 3 : 2
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        if event
            xerr, verr, aerr = stage.carried
            v = _exact_rational(stage.event_candidates[k].v[1][p])
            acceleration = stage.acceleration[2][i]
            Z = max(abs(v + acceleration * lower), abs(v + acceleration * horizon))
            X[i] = xerr[i] + horizon * verr[i] + horizon^2 * aerr[i] / 2 +
                J[i] * horizon^3 / 6 + cluster.delta * Z + abs(dx[i])
        else
            X[i] = abs(dx[i]) + J[i] * horizon^3 / 6
        end
    end
    weighted = _yi_weighted_error(stage, species, X)
    rho = maximum(abs, stage.rho)
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        force = _yi_force_error(stage, i, k, weighted, X, rho, geom)
        V[i] = abs(kappa[i]) + horizon * force / 2 + J[i] * horizon^2 / 4
        if event
            V[i] += stage.carried[2][i] + horizon * stage.carried[3][i] / 2 +
                cluster.delta * abs(stage.acceleration[first][i] + stage.acceleration[next][i]) / 2
        end
    end
    return nothing
end

function step!(st::YeeImpactLongitudinalPIC{T}, particles, dt::Real) where {T}
    _yl_require_mode()
    species = _yl_species(particles)
    tau = _require_finite_positive_real("dt", dt, T)
    st.initialized[] || throw(ArgumentError("call init! before impact stepping"))
    _yi_validate_phase(st)
    st.step[] == 0 || isequal(tau, st.last_dt[]) ||
        throw(ArgumentError("impact dt is fixed after first success"))
    next_time, next_step = _validated_step_clock(st.time[], st.step[], tau)
    next_time > st.time[] || throw(ArgumentError("impact time must advance finitely"))
    _yl_validate_populations(st, species)
    _yi_validate_bound(st, species)
    _yl_validate_longitudinal_fields(st, st.fields)
    oldkinetic = _kinetic_energy_checked(species, T, Val(false))
    isfinite(oldkinetic) || throw(ArgumentError("impact old kinetic diagnostic is nonfinite"))
    oldK = _kinetic_energy_rational(species)
    oldU = _field_energy_rational(st.fields.E, st.fields.B, one(_YeeImpactQ), st.g.dx)
    geom = _yi_global_bounds!(st, species, _exact_rational(tau))
    stage = st._stage[]
    fill!(stage.rho, zero(_YeeImpactQ))
    fill!(stage.net_impulse, zero(_YeeImpactQ))
    fill!(stage.selected_side, 0)
    _yi_first_brackets!(st, species, geom)
    cluster = nothing
    counts, walls = (0, 0), (zero(_YeeImpactQ), zero(_YeeImpactQ))
    if isempty(stage.roots)
        for i in eachindex(stage.selected_side), side = 1:2
            c = _yi_coefficients(st, species, i, side, stage.acceleration[1])
            _yi_minimum(c, zero(geom.tau), geom.tau) > stage.bounds[3][i] * geom.tau^3 / 6 ||
                throw(ArgumentError("physical no-contact certificate is inconclusive"))
        end
        local_step = _yi_local_chunk!(st, species, st.fields, stage.candidates,
            st._candidate, geom.tau, 1, 2, false)
        work, electric = local_step.work, local_step.electric
    else
        cluster = _yi_certify_cluster!(st, species, geom)
        local_step = _yi_local_chunk!(st, species, st.fields, stage.event_candidates,
            st._event_fields, cluster.duration, 1, 2, true)
        _yi_carry!(st, species, cluster, geom)
        for (k, ps) in enumerate(stage.event_candidates), p in eachindex(ps.weight)
            i = stage.offsets[k] + p
            side = stage.selected_side[i]
            side == 0 && continue
            incoming = side == 1 ? -_exact_rational(ps.v[1][p]) : _exact_rational(ps.v[1][p])
            incoming >= stage.terminal[2][i] / 2 ||
                throw(ArgumentError("completed impact kick is not sufficiently incoming"))
            counts = Base.setindex(counts, Base.checked_add(counts[side], 1), side)
            impulse = 2 * _exact_rational(ps.m) * _exact_rational(ps.weight[p]) * _exact_rational(ps.v[1][p])
            walls = Base.setindex(walls, walls[side] + impulse, side)
            ps.v[1][p] = -ps.v[1][p]
        end
        remainder = _yi_clear_remainder!(st, cluster, geom)
        final_step = _yi_local_chunk!(st, stage.event_candidates, st._event_fields,
            stage.candidates, st._candidate, remainder, 2, 3, false)
        work = local_step.work + final_step.work
        electric = local_step.electric + final_step.electric
    end
    _yi_terminal!(st, species, cluster, geom)
    newkinetic = kinetic_energy(stage.candidates)
    newK = _kinetic_energy_rational(stage.candidates)
    newU = _field_energy_rational(st._candidate.E, st._candidate.B, one(_YeeImpactQ), st.g.dx)
    lastwork = _from_exact(T, work)
    defect = _from_exact(T, newK - oldK - work)
    physical = _from_exact(T, (newK - oldK) + (newU - oldU))
    all(isfinite, (newkinetic, lastwork, defect, physical)) ||
        throw(ArgumentError("impact required energy/work diagnostics do not round finite"))
    for i in eachindex(stage.net_impulse)
        value = _from_exact(T, stage.net_impulse[i])
        isfinite(value) || throw(ArgumentError("net impact impulse does not round finite"))
        st._impulse[1][i] = value
    end
    cumulative_count = ntuple(side -> Base.checked_add(st.cumulative_wall_count[side][], counts[side]), 2)
    cumulative_wall = ntuple(side -> st.cumulative_wall_impulse[side][] + walls[side], 2)
    qz = zero(_YeeImpactQ)
    snapenergy = cluster === nothing ? qz : geom.delta *
        (geom.field * cluster.charge_snap / geom.h + geom.M * cluster.charge_snap^2 / (2 * geom.h))
    history = (initialized = true, time = next_time, step = next_step, last_dt = tau,
        last_impulse = st._impulse, last_work = lastwork, last_energy_defect = defect,
        last_physical_energy_change = physical, last_wall_count = counts,
        cumulative_wall_count = cumulative_count, last_wall_impulse = walls,
        cumulative_wall_impulse = cumulative_wall, last_electric_impulse = electric,
        last_max_snap = cluster === nothing ? qz : cluster.max_snap,
        last_snap_energy_bound = snapenergy,
        last_predictor_time_error = cluster === nothing ? qz : cluster.predictor_error,
        last_physical_time_error = cluster === nothing ? qz : cluster.delta,
        last_source_invariant_error = maximum(abs, stage.rho),
        last_position_error_bound = maximum(stage.terminal[1]; init = qz),
        last_velocity_error_bound = maximum(stage.terminal[2]; init = qz))
    _yi_validate_history_values(history, st.g, last(stage.offsets))
    # All failing work precedes the first public store.
    for (ps, candidate) in zip(species, stage.candidates)
        copyto!(ps.x[1], candidate.x[1])
        copyto!(ps.v[1], candidate.v[1])
    end
    copyto!(st.fields.E[1], st._candidate.E[1])
    foreach(copyto!, st.last_impulse, st._impulse)
    st.time[], st.step[], st.last_dt[], st.last_work[] = next_time, next_step, tau, lastwork
    _yi_publish_extra_history_refs!(st, history)
    return st
end

# Exact operation-scale ulp, including scales outside the IEEE output range.
function _yi_ulp(::Type{T}, scale::_YeeImpactQ) where {T}
    minimum = exponent(floatmin(T)) - precision(T) + 1
    if iszero(scale)
        quantum = minimum
    else
        n, d = numerator(scale), denominator(scale)
        power = ndigits(n; base = 2) - ndigits(d; base = 2)
        (power >= 0 ? n < d << power : n << -power < d) && (power -= 1)
        quantum = max(power - precision(T) + 1, minimum)
    end
    return quantum >= 0 ? (BigInt(1) << quantum) // BigInt(1) :
        BigInt(1) // (BigInt(1) << -quantum)
end
function _yee_longitudinal_cluster_candidate(st, ::Type{T}, species, geom) where {T}
    stage = st._stage[]
    isempty(stage.contacts) && return nothing
    last = 1
    upper = stage.contacts[1].hi
    while last < length(stage.contacts) && stage.contacts[last+1].lo <= upper
        last += 1
        upper = max(upper, stage.contacts[last].hi)
    end
    group = view(stage.contacts, 1:last)
    lower = minimum(c -> c.lo, group)
    s = minimum(c -> stage.roots[c.bracket].hi, group)
    lower <= s <= upper || return nothing
    fill!(stage.selected_side, 0)
    fill!(stage.residual[1], zero(s))
    fill!(stage.terminal[2], zero(s)) # Incoming margins, dead before terminal output.
    for contact in group
        r = stage.roots[contact.bracket]
        stage.selected_side[r.particle] = r.side
        stage.terminal[2][r.particle] = contact.incoming
    end
    # Every omitted particle has physical first-hit clearance through the hull.
    for i in eachindex(stage.selected_side)
        stage.selected_side[i] != 0 && continue
        H = stage.bounds[3][i] * geom.tau^3 / 6
        for side = 1:2
            c = _yi_coefficients(st, species, i, side, stage.acceleration[1])
            _yi_minimum(c, zero(s), upper) > H || return nothing
        end
    end
    eta, maxsnap, Rcharge = zero(s), zero(s), zero(s)
    for contact in group
        r = stage.roots[contact.bracket]
        ps, p, k = _yi_particle(stage, species, r.particle)
        x, v = _exact_rational(ps.x[1][p]), _exact_rational(ps.v[1][p])
        field = stage.gather[1][r.particle]
        qz = zero(field)
        out = _yee_particle_half_kick((ps.v[1][p], zero(T), zero(T)),
            (field, qz, qz), (qz, qz, qz), stage.factors[k], s / 2)
        all(isfinite, out) || return nothing
        u = _exact_rational(out[1])
        free = _particle_coordinate_native(x, u, s)
        wall = r.side == 1 ? zero(s) : geom.L
        snap = wall - free
        abs(snap) <= 4 * _yi_ulp(T, abs(x) + abs(s * u)) || return nothing
        c = _yi_coefficients(st, species, r.particle, r.side, stage.acceleration[1])
        nu = -max(_yi_slope(c, min(s, r.lo)), _yi_slope(c, max(s, r.hi)))
        nu > 0 || return nothing
        predictor = wall - (x + s * v + stage.acceleration[1][r.particle] * s^2 / 2)
        eta = max(eta, abs(predictor) / nu)
        maxsnap = max(maxsnap, abs(snap))
        stage.residual[1][r.particle] = snap
        Rcharge += abs(_exact_rational(ps.q) * _exact_rational(ps.weight[p])) * abs(snap)
    end
    A, V, _ = stage.bounds
    for (k, ps) in enumerate(species), p in eachindex(ps.weight)
        i = stage.offsets[k] + p
        x, v = _exact_rational(ps.x[1][p]), _exact_rational(ps.v[1][p])
        V[i] * eta + A[i] * eta^2 / 2 <= 4 * _yi_ulp(T, abs(x) + geom.tau * V[i]) &&
            2A[i] * eta <= 4 * _yi_ulp(T, abs(v) + geom.tau * A[i]) || return nothing
        snap = abs(stage.residual[1][i])
        force = abs(stage.factors[k]) * geom.delta *
            (geom.M / geom.h * Rcharge + 2 * geom.field / (geom.a * geom.h) * snap)
        s / 2 * force <= 4 * _yi_ulp(T, abs(v) + geom.tau * A[i]) || return nothing
    end
    return (duration = s, lower = lower, upper = upper,
        delta = max(s - lower, upper - s), width = upper - lower,
        group = last, predictor_error = eta, max_snap = maxsnap, charge_snap = Rcharge)
end
_yi_cluster(st::YeeImpactLongitudinalPIC{T}, species, geom) where {T} =
    _yee_longitudinal_cluster_candidate(st, T, species, geom)

function _yi_certify_cluster!(st, species, geom)
    for iteration = 0:st.max_root_refinements
        if _yi_contacts!(st, species, geom)
            cluster = _yi_cluster(st, species, geom)
            cluster === nothing || return cluster
        end
        iteration == st.max_root_refinements && break
        _yi_refine!(st, species)
    end
    throw(ArgumentError("impact root/cluster certificate is unresolved within its finite budget"))
end
