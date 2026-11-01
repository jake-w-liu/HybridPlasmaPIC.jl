"""
    YeeElectromagneticPIC(g; shape::YeeTriangularCloud, c=1,
                         background_charge=zero(T), left_plate_charge=zero(T))

Synchronized fixed-step, event-free nonrelativistic CPU Float32/64 1D3V
electromagnetic cloud
PIC on a closed PEC grid. All x/v/E/B are at integer time. Set admissible
transverse E/B before `init!`; Ex is initialized from the existing cloud Gauss
and isolated-plate convention. Bx must be stored zero; tangential E endpoints
must be stored zero. Both transverse polarizations and velocities evolve.

Each half-kick retains complete original-geometry GE/GB, original q/m and the
original resolved full-dt half through canonical Boris. Only its final velocity
components round to grid precision. A gathered component need not fit grid
precision when the completed velocity is finite. Ordinary native velocity bits
are not preserved; other stage quantities still use their existing IEEE owners.

Explicitly select a fixed physical cloud width. Nearest rounding and gradual
underflow, frozen source/provenance/storage identities, strict CFL/rest guards
and a conservative full-speed physical-energy exit tube are required. These
are not general nonlinear stability or accuracy certificates. Candidate drift
endpoints must also be strictly interior. No reflection, clipping or wrapping
is supplied. Every required candidate/diagnostic is completed before public
publication; private scratch can change on rejection.

`last_modified_energy_defect` is rounded-diagnostic DeltaK minus field source
work. `last_physical_energy_change` is the complete original-input change of
physical K+U at the two stored states. Stepping requires finite old/new rounded
kinetic diagnostics, source work and stored histories. Neither defect is forced
to zero. This model
has a distinct checkpoint schema and one clock; all arrays retain identity.
Concurrent mutation or stepping of one instance is unsupported.
"""
struct YeeElectromagneticPIC{T,S<:YeeTriangularCloud,B,P}
    g::YeeGrid{1,T}
    shape::S
    background_charge::B
    left_plate_charge::P
    c::T
    fields::YeeFields{1,T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    initialized::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    right_plate_charge::Base.RefValue{Rational{BigInt}}
    last_impulse::NTuple{3,Vector{T}}
    last_work::Base.RefValue{T}
    last_modified_energy_defect::Base.RefValue{T}
    last_physical_energy_change::Base.RefValue{T}
    _candidate::YeeFields{1,T}
    _impulse::NTuple{3,Vector{T}}
    _charge_workspace::YeeChargeWorkspace{T}
    _stage::Base.RefValue{_YeeLongitudinalStage{T}}
    _cfl::T
end
function YeeElectromagneticPIC(g::YeeGrid{1,T}; shape::YeeTriangularCloud,
    c = 1, background_charge = zero(T), left_plate_charge = zero(T)) where {T}
    g.periodic[1] && throw(ArgumentError("electromagnetic cloud PIC requires closed PEC"))
    b = _yl_source_scalar("background_charge", background_charge, T)
    sigma = _yl_source_scalar("left_plate_charge", left_plate_charge, T)
    speed = _require_finite_positive_real("c", c, T)
    fields, candidate = YeeFields(g), YeeFields(g)
    impulse() = ntuple(c -> zeros(T, length(fields.E[c])), 3)
    stage = _YeeLongitudinalStage{T}(ParticleSet{1,T}[], Tuple[], _ElectrostaticHostParticle{1,T}[])
    return YeeElectromagneticPIC{T,typeof(shape),typeof(b),typeof(sigma)}(
        g, shape, b, sigma, speed, fields, Ref(zero(T)), Ref(0), Ref(false), Ref(zero(T)),
        Ref(-_exact_rational(sigma) - _exact_rational(b) * _exact_rational(g.L[1])),
        impulse(), Ref(zero(T)), Ref(zero(T)), Ref(zero(T)), candidate, impulse(),
        YeeChargeWorkspace(g), Ref(stage), _yee_cfl(g, speed))
end
YeeElectromagneticPIC(g::YeeGrid; kwargs...) =
    throw(ArgumentError("electromagnetic cloud PIC requires one spatial dimension"))
maxwell_cfl(st::YeeElectromagneticPIC) = st._cfl
_yee_closed_references(st::YeeElectromagneticPIC) =
    (st.time, st.last_dt, st.last_work, st.last_modified_energy_defect, st.last_physical_energy_change)
_yl_validate_populations(st::YeeElectromagneticPIC, species) =
    _yee_closed_validate_populations(st, species)

function _ye_validate_fields(st::YeeElectromagneticPIC, fields::YeeFields)
    _yee_validate_fields(fields, st.g)
    for a in (fields.E[2], fields.E[3])
        _yl_stored_zero(first(a)) && _yl_stored_zero(last(a)) ||
            throw(ArgumentError("electromagnetic PEC requires stored zero tangential E endpoints"))
    end
    all(_yl_stored_zero, fields.B[1]) ||
        throw(ArgumentError("one-dimensional divergence-free PEC requires stored zero Bx"))
    return nothing
end
function _ye_validate_phase_values(initialized, time, step, last_dt, work, modified,
    physical, impulse, g, cfl)
    _yl_validate_phase_values(initialized, time, step, last_dt, work, modified)
    physical isa typeof(time) && isfinite(physical) && (step != 0 || iszero(physical)) ||
        throw(ArgumentError("invalid physical energy-change history"))
    step == 0 || last_dt < cfl || throw(ArgumentError("saved fixed dt must remain below Maxwell CFL"))
    _yee_validate_impulse(impulse, g)
    for c = 2:3
        _yl_stored_zero(first(impulse[c])) && _yl_stored_zero(last(impulse[c])) ||
            throw(ArgumentError("electromagnetic current history requires zero PEC endpoints"))
    end
    step == 0 && !all(a -> all(_yl_stored_zero, a), impulse) &&
        throw(ArgumentError("initial electromagnetic impulse must be zero"))
    return nothing
end
_ye_validate_phase(st::YeeElectromagneticPIC) = _ye_validate_phase_values(st.initialized[],
    st.time[], st.step[], st.last_dt[], st.last_work[], st.last_modified_energy_defect[],
    st.last_physical_energy_change[], st.last_impulse, st.g, maxwell_cfl(st))

function init!(st::YeeElectromagneticPIC{T}, particles) where {T}
    species = _yl_species(particles)
    _yl_validate_populations(st, species)
    _ye_validate_fields(st, st.fields)
    stage = _yl_make_stage(species, species, T)
    right = _yl_right_plate(st, species)
    _yl_initialize_fields!(st, species)
    copyto!(st.fields.E[1], st._candidate.E[1])
    for a in st.last_impulse
        fill!(a, zero(T))
    end
    st._stage[] = stage
    st.right_plate_charge[] = right
    st.time[], st.step[], st.last_dt[] = zero(T), 0, zero(T)
    st.last_work[], st.last_modified_energy_defect[], st.last_physical_energy_change[] = zero(T), zero(T), zero(T)
    st.initialized[] = true
    return st
end

_ye_physical_energy(st, species, fields) = _kinetic_energy_rational(species) +
    _field_energy_rational(fields.E, fields.B, _exact_rational(st.c)^2, st.g.dx)

# Exact dyadic upper enclosure of sqrt(s), with no finite-float bound required.
function _ye_sqrt_upper(s::Rational{BigInt})
    iszero(s) && return s
    n, d = numerator(s), denominator(s)
    k = ndigits(n; base = 2) - ndigits(d; base = 2)
    (k >= 0 ? n < d << k : n << -k < d) && (k -= 1)
    power = cld(k + 1, 2)
    return power >= 0 ? (BigInt(1) << power) // BigInt(1) : BigInt(1) // (BigInt(1) << -power)
end
function _ye_admit_step(st::YeeElectromagneticPIC, species, dt, energy)
    Q = _exact_rational
    tau, delta, L = Q(dt), Q(st.g.dx[1]), Q(st.g.L[1])
    h, a = L / st.g.n[1], Q(st.shape.half_width)
    stiffness = sum((Q(ps.q)^2 * Q(w) / Q(ps.m) for ps in species for w in ps.weight);
        init = zero(tau))
    tau^2 * delta * min(inv(h), inv(a)) / h * stiffness < 4 ||
        throw(ArgumentError("dt violates longitudinal rest-mode bound"))
    tau^2 * (4 * Q(st.c)^2 / delta^2 + stiffness / delta) < 4 ||
        throw(ArgumentError("dt violates combined transverse Maxwell/plasma rest-mode bound"))
    fieldbound = max(one(h), delta / h) * _ye_sqrt_upper(2 * energy / delta)
    for ps in species
        acceleration = abs(Q(ps.q) / Q(ps.m)) * fieldbound
        for p in eachindex(ps.weight)
            x = Q(ps.x[1][p])
            speed = sum(abs(Q(ps.v[c][p])) for c = 1:3)
            tau * speed + tau^2 * acceleration / 2 < min(x, L - x) ||
                throw(ArgumentError("dt fails electromagnetic full-speed physical-energy exit tube"))
        end
    end
    return nothing
end
function _ye_gather(fields, x, g, shape)
    E = (_yee_cloud_field_rational(fields.E[1], x, nothing, g, shape),
        _yee_cloud_nodal_rational(fields.E[2], x, g, shape),
        _yee_cloud_nodal_rational(fields.E[3], x, g, shape))
    B = (zero(E[1]), _yee_cloud_field_rational(fields.B[2], x, nothing, g, shape),
        _yee_cloud_field_rational(fields.B[3], x, nothing, g, shape))
    return E, B
end

function step!(st::YeeElectromagneticPIC{T}, particles, dt::Real) where {T}
    _yl_require_mode()
    species = _yl_species(particles)
    tau = _require_finite_positive_real("dt", dt, T)
    st.initialized[] || throw(ArgumentError("call init! before electromagnetic stepping"))
    _ye_validate_phase(st)
    st.step[] == 0 || isequal(tau, st.last_dt[]) ||
        throw(ArgumentError("electromagnetic dt is fixed after first success"))
    tau < maxwell_cfl(st) || throw(ArgumentError("dt must be strictly below Maxwell CFL"))
    next_time, next_step = _validated_step_clock(st.time[], st.step[], tau)
    next_time > st.time[] || throw(ArgumentError("electromagnetic time must advance finitely"))
    _yl_validate_populations(st, species)
    _yl_validate_bound(st, species)
    _ye_validate_fields(st, st.fields)
    oldenergy = _ye_physical_energy(st, species, st.fields)
    _ye_admit_step(st, species, tau, oldenergy)
    oldkinetic = _kinetic_energy_checked(species, T, Val(false))
    stage, g = st._stage[], st.g
    half = _exact_rational(tau) / 2
    factors = [_boris_exact_factor(_validated_boris_factor(ps)) for ps in species]
    for (ps, candidate, factor) in zip(species, stage.candidates, factors), p in eachindex(ps.weight)
        x = ps.x[1][p]
        E, B = _ye_gather(st.fields, x, g, st.shape)
        u = _yee_particle_half_kick(ntuple(c -> ps.v[c][p], 3), E, B, factor, half)
        all(isfinite, u) || throw(ArgumentError("electromagnetic first kick does not round finite"))
        y = _particle_coordinate(x, u[1], tau)
        isfinite(y) && 0 < _exact_rational(y) < _exact_rational(g.L[1]) ||
            throw(ArgumentError("electromagnetic drift leaves strict interior"))
        candidate.x[1][p] = y
        for c = 1:3
            candidate.v[c][p] = u[c]
        end
    end
    old, new = [ps.x[1] for ps in species], [ps.x[1] for ps in stage.candidates]
    yee_longitudinal_impulse!(st._impulse[1], species, old, new, g, st._charge_workspace; shape = st.shape)
    _yee_transverse_impulse!(st._impulse, species, stage.candidates, g, st.shape, tau)
    _yee_candidate_step!(st._candidate, st.fields, g, tau, st.c, st._impulse)
    _ye_validate_fields(st, st._candidate)
    for (candidate, factor) in zip(stage.candidates, factors), p in eachindex(candidate.weight)
        E, B = _ye_gather(st._candidate, candidate.x[1][p], g, st.shape)
        v = _yee_particle_half_kick(ntuple(c -> candidate.v[c][p], 3), E, B, factor, half)
        all(isfinite, v) || throw(ArgumentError("electromagnetic final kick does not round finite"))
        for c = 1:3
            candidate.v[c][p] = v[c]
        end
    end
    newkinetic = kinetic_energy(stage.candidates)
    work = yee_source_work(st.fields.E, st._candidate.E, st._impulse, g)
    all(isfinite, (oldkinetic, newkinetic, work)) ||
        throw(ArgumentError("electromagnetic work/kinetic diagnostics must round finite"))
    defect = _from_exact(T, _exact_rational(newkinetic) - _exact_rational(oldkinetic) - _exact_rational(work))
    physical = _from_exact(T, _ye_physical_energy(st, stage.candidates, st._candidate) - oldenergy)
    all(isfinite, (defect, physical)) || throw(ArgumentError("electromagnetic energy histories must round finite"))
    # Every required candidate, diagnostic and clock is complete before publication.
    for (ps, candidate) in zip(species, stage.candidates)
        copyto!(ps.x[1], candidate.x[1])
        foreach(copyto!, ps.v, candidate.v)
    end
    foreach(copyto!, st.fields.E, st._candidate.E)
    foreach(copyto!, st.fields.B, st._candidate.B)
    foreach(copyto!, st.last_impulse, st._impulse)
    st.time[], st.step[], st.last_dt[] = next_time, next_step, tau
    st.last_work[], st.last_modified_energy_defect[], st.last_physical_energy_change[] = work, defect, physical
    return st
end
