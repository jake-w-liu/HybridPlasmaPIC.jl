# One synchronized, event-free longitudinal model. Static cloud primitives own
# all shape algebra; the Boris, coordinate and sourced Yee row owners own stages.
struct _YeeLongitudinalStage{T}
    owners::Vector{ParticleSet{1,T}}
    components::Vector{Tuple}
    candidates::Vector{_ElectrostaticHostParticle{1,T}}
end

"""
    YeeLongitudinalPIC(g; shape::YeeTriangularCloud,
                       background_charge=zero(T), left_plate_charge=zero(T))

Synchronized fixed-step CPU Float32/64 longitudinal electrostatic cloud PIC on a
closed 1D Yee grid. Call `init!(state, particles)` then `step!(state, particles, dt)`;
particles may be one `ParticleSet` or an ordered one-based vector. Initialization,
stepping and restart require nearest hardware rounding and gradual underflow;
other contexts reject before staging. Static cloud APIs retain their own policy. Explicitly
select a physical cloud width, held fixed under mesh refinement. This model has
no transverse dynamics or wall events: transverse fields/velocities and any
step failing its strict conservative impact-exclusion tube are rejected.

The original Float32/64 shape width, background density and left free-plate
charge retain their type and bits. Other real background/plate inputs resolve
once to grid precision. `right_plate_charge[]` is the exact derived free charge,
separate from cloud endpoint tails. Initialization builds Gauss-consistent Ex
with longitudinal DC and resets the clock/history. Source, provenance, particle
objects and component identities are fixed until explicit reinitialization.

Each half-kick retains the complete original-geometry gather, original q/m and
original resolved full-dt half through canonical Boris, rounding only its final
velocity to grid precision. A gathered field need not fit grid precision when
the completed kick is finite. Ordinary native velocity bits are not preserved.

Velocity-Verlet advances integer-time x/v/Ex. The first accepted positive dt
fixes the resolved timestep. The rest-mode guard is not a general nonlinear
stability or accuracy certificate. `last_energy_defect[]` is measured
`Knew-Kold-last_work[]`, not an exactly conserved energy; kinetic energies use
the existing rounded diagnostic owner. Physical arrays retain their identities.
A rejected operation leaves physical state and diagnostics unchanged; private
scratch may change. One instance must not be stepped concurrently.
"""
struct YeeLongitudinalPIC{T,S<:YeeTriangularCloud,B,P}
    g::YeeGrid{1,T}
    shape::S
    background_charge::B
    left_plate_charge::P
    fields::YeeFields{1,T}
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    initialized::Base.RefValue{Bool}
    last_dt::Base.RefValue{T}
    right_plate_charge::Base.RefValue{Rational{BigInt}}
    last_impulse::NTuple{3,Vector{T}}
    last_work::Base.RefValue{T}
    last_energy_defect::Base.RefValue{T}
    _candidate::YeeFields{1,T}
    _impulse::NTuple{3,Vector{T}}
    _charge_workspace::YeeChargeWorkspace{T}
    _stage::Base.RefValue{_YeeLongitudinalStage{T}}
end

function _yl_source_scalar(name, value, ::Type{T}) where {T}
    if value isa Union{Float32,Float64}
        isfinite(value) || throw(ArgumentError("$name must be finite"))
        return value
    end
    return _require_finite_real(name, value, T)
end

function YeeLongitudinalPIC(
    g::YeeGrid{1,T}; shape::YeeTriangularCloud,
    background_charge = zero(T), left_plate_charge = zero(T),
) where {T}
    g.periodic[1] && throw(ArgumentError("longitudinal PIC requires a closed PEC grid"))
    b = _yl_source_scalar("background_charge", background_charge, T)
    sigma = _yl_source_scalar("left_plate_charge", left_plate_charge, T)
    fields, candidate = YeeFields(g), YeeFields(g)
    impulse() = ntuple(c -> zeros(T, length(fields.E[c])), 3)
    stage = _YeeLongitudinalStage{T}(ParticleSet{1,T}[], Tuple[], _ElectrostaticHostParticle{1,T}[])
    return YeeLongitudinalPIC{T,typeof(shape),typeof(b),typeof(sigma)}(
        g, shape, b, sigma, fields, Ref(zero(T)), Ref(0), Ref(false), Ref(zero(T)),
        Ref(-_exact_rational(sigma) - _exact_rational(b) * _exact_rational(g.L[1])),
        impulse(), Ref(zero(T)), Ref(zero(T)), candidate, impulse(),
        YeeChargeWorkspace(g), Ref(stage),
    )
end

YeeLongitudinalPIC(g::YeeGrid; kwargs...) =
    throw(ArgumentError("longitudinal PIC requires one spatial dimension"))

_yl_species(ps::ParticleSet) = (ps,)
function _yl_species(species::AbstractVector)
    Base.require_one_based_indexing(species)
    return species
end
_yl_species(species) =
    throw(ArgumentError("closed Yee PIC particles must be a ParticleSet or ordered vector"))
_yl_particle_arrays(ps::ParticleSet) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)

# Structural zero admission must inspect stored bits, including when the caller
# has enabled flush-to-zero before entering an API that will reject that mode.
_yl_stored_zero(x::Float32) = iszero(reinterpret(UInt32, x) & 0x7fffffff)
_yl_stored_zero(x::Float64) = iszero(reinterpret(UInt64, x) & 0x7fffffffffffffff)
function _yl_require_mode()
    _boris_accuracy_mode() ||
        throw(ArgumentError("closed Yee PIC requires nearest rounding and gradual underflow"))
    return nothing
end

function _yl_make_stage(species, sources, ::Type{T}) where {T}
    length(species) == length(sources) || throw(ArgumentError("population counts must match"))
    owners = ParticleSet{1,T}[ps for ps in species]
    components = Tuple[_yl_particle_arrays(ps) for ps in species]
    # Candidate storage belongs to the model, independently of public views.
    candidates = _ElectrostaticHostParticle{1,T}[]
    _pic_prepare_species_stages!(candidates, sources)
    return _YeeLongitudinalStage{T}(owners, components, candidates)
end

function _yl_check_particle_storage(ps::ParticleSet)
    _check_particle_storage_axes(ps)
    all(a -> a isa StridedVector && _particle_batch_injective(a), _yl_particle_arrays(ps)) ||
        throw(ArgumentError("closed Yee PIC particles require injective CPU strided storage"))
    return nothing
end

function _yl_storage_arrays(st)
    arrays = Any[
        st.fields.E..., st.fields.B..., st._candidate.E..., st._candidate.B...,
        st.last_impulse..., st._impulse..., st._charge_workspace.sum,
        st._charge_workspace.candidate,
    ]
    for ps in st._stage[].candidates
        _yl_check_particle_storage(ps)
        append!(arrays, _yl_particle_arrays(ps))
    end
    return arrays
end

function _yl_validate_longitudinal_fields(st::YeeLongitudinalPIC, fields::YeeFields)
    _yee_validate_fields(fields, st.g)
    all(a -> all(_yl_stored_zero, a), (fields.E[2], fields.E[3], fields.B...)) ||
        throw(ArgumentError("longitudinal PIC requires zero transverse E and all B"))
    return nothing
end

function _yl_validate_storage(st, species)
    T = eltype(st.fields.E[1])
    _yee_check_shape(st.fields.E, st.g, true)
    _yee_check_shape(st.fields.B, st.g, false)
    _yee_check_shape(st._candidate.E, st.g, true)
    _yee_check_shape(st._candidate.B, st.g, false)
    _yee_check_shape(st.last_impulse, st.g, true)
    _yee_check_shape(st._impulse, st.g, true)
    work = st._charge_workspace
    work.n == st.g.n[1] && work.L === st.g.L[1] ||
        throw(ArgumentError("closed Yee PIC charge workspace does not match grid"))
    length(work.sum) == st.g.n[1] + 1 && length(work.candidate) == st.g.n[1] + 1 ||
        throw(DimensionMismatch("closed Yee PIC charge workspace shape changed"))
    clocks = _yee_closed_references(st)
    for j in eachindex(clocks), i = 1:j-1
        clocks[i] === clocks[j] &&
            throw(ArgumentError("closed Yee PIC clock and diagnostic references must be independent"))
    end
    arrays = _yl_storage_arrays(st)
    for ps in species
        ps isa ParticleSet{1,T} ||
            throw(ArgumentError("closed Yee PIC particles must match grid precision and dimension"))
        _yl_check_particle_storage(ps)
        append!(arrays, _yl_particle_arrays(ps))
    end
    for j in eachindex(arrays), i = 1:j-1
        _particle_arrays_mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("closed Yee PIC public, stage and particle arrays must be independent"))
    end
    return nothing
end

_yee_closed_references(st::YeeLongitudinalPIC) =
    (st.time, st.last_dt, st.last_work, st.last_energy_defect)

function _yee_closed_validate_populations(st, species)
    _yl_require_mode()
    _yl_validate_storage(st, species)
    L = _exact_rational(st.g.L[1])
    for ps in species
        _validated_charge_mass(ps)
        all(isfinite, ps.weight) && all(w -> _exact_rational(w) >= 0, ps.weight) ||
            throw(ArgumentError("closed Yee weights must be finite and nonnegative"))
        all(a -> all(isfinite, a), ps.v) ||
            throw(ArgumentError("closed Yee velocities must be finite"))
        all(x -> isfinite(x) && 0 < _exact_rational(x) < L, ps.x[1]) ||
            throw(ArgumentError("closed Yee centers must lie strictly inside (0,L)"))
    end
    return nothing
end
function _yl_validate_populations(st::YeeLongitudinalPIC, species)
    _yee_closed_validate_populations(st, species)
    for ps in species
        all(_yl_stored_zero, ps.v[2]) && all(_yl_stored_zero, ps.v[3]) ||
            throw(ArgumentError("longitudinal velocities must have zero transverse components"))
    end
    return nothing
end

function _yl_validate_phase_values(initialized, time::T, step, last_dt, work, defect) where {T}
    initialized isa Bool && step isa Int && step >= 0 && last_dt isa T &&
        work isa T && defect isa T && all(isfinite, (time, last_dt, work, defect)) ||
        throw(ArgumentError("invalid closed Yee PIC phase or diagnostics"))
    if step == 0
        iszero(time) && iszero(last_dt) && iszero(work) && iszero(defect) ||
            throw(ArgumentError("unadvanced closed Yee PIC clock/history must be zero"))
    else
        initialized && time > zero(T) && last_dt > zero(T) ||
            throw(ArgumentError("advanced closed Yee PIC phase must be initialized with positive dt/time"))
    end
    return nothing
end
function _yl_validate_phase(st::YeeLongitudinalPIC)
    _yl_validate_phase_values(st.initialized[], st.time[], st.step[], st.last_dt[],
        st.last_work[], st.last_energy_defect[])
    _yee_validate_impulse(st.last_impulse, st.g)
    all(_yl_stored_zero, st.last_impulse[2]) && all(_yl_stored_zero, st.last_impulse[3]) ||
        throw(ArgumentError("longitudinal impulse history must have zero transverse data"))
    st.step[] == 0 && !all(a -> all(iszero, a), st.last_impulse) &&
        throw(ArgumentError("unadvanced longitudinal impulse must be zero"))
    return nothing
end

function _yl_right_plate(st, species)
    total = _exact_rational(st.left_plate_charge) +
        _exact_rational(st.background_charge) * _exact_rational(st.g.L[1])
    for ps in species, w in ps.weight
        total += _exact_rational(ps.q) * _exact_rational(w)
    end
    return -total
end

function _yl_validate_bound(st, species)
    stage = st._stage[]
    length(species) == length(stage.owners) == length(stage.candidates) == length(stage.components) ||
        throw(ArgumentError("population count changed; call init! explicitly"))
    for (k, ps) in enumerate(species)
        ps === stage.owners[k] || throw(ArgumentError("population identity/order changed; call init!"))
        all(t -> t[1] === t[2], zip(_yl_particle_arrays(ps), stage.components[k])) ||
            throw(ArgumentError("particle component identity changed; call init!"))
        frozen = stage.candidates[k]
        isequal(ps.q, frozen.q) && isequal(ps.m, frozen.m) &&
            isequal(ps.weight, frozen.weight) && isequal(ps.id, frozen.id) && isequal(ps.tag, frozen.tag) ||
            throw(ArgumentError("particle source/provenance changed; call init!"))
    end
    st.right_plate_charge[] == _yl_right_plate(st, species) ||
        throw(ArgumentError("derived closed Yee PIC free-plate charge changed"))
    return nothing
end

function _yl_initialize_fields!(st, species)
    T = eltype(st.fields.E[1])
    g, shape = st.g, st.shape
    h = _exact_rational(g.L[1]) / g.n[1]
    b, sigma = _exact_rational(st.background_charge), _exact_rational(st.left_plate_charge)
    charge = fill(zero(Rational{BigInt}), g.n[1] + 1)
    _yee_charge_kernel!(charge, species, nothing, nothing, g, Val(true), Val(:charge), shape)
    cumulative = zero(h)
    for e = 0:g.n[1]-1
        # The exact nodal prefix telescopes to sum(q*w*(1-A_e)); do not
        # convert nodal charges, which can overflow even when Ex is finite.
        cumulative += charge[e+1]
        total = sigma + b * (e + 1//2) * h + cumulative
        value = _from_exact(T, total)
        isfinite(value) || throw(ArgumentError("initial closed Yee PIC Ex does not round finite"))
        st._candidate.E[1][e+1] = value
    end
    return nothing
end

function init!(st::YeeLongitudinalPIC{T}, particles) where {T}
    _yl_require_mode()
    species = _yl_species(particles)
    _yl_validate_populations(st, species)
    _yl_validate_longitudinal_fields(st, st.fields)
    # All new bindings and storage are allocated before any public write.
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
    st.last_work[], st.last_energy_defect[] = zero(T), zero(T)
    st.initialized[] = true
    return st
end

# Exact original-operand admission, before staging. This bound excludes a first
# physical exit; it deliberately rejects some safe trajectories as well.
function _yl_admit_step(st::YeeLongitudinalPIC, species, dt)
    g = st.g
    L, delta, tau = _exact_rational(g.L[1]), _exact_rational(g.dx[1]), _exact_rational(dt)
    h, a = L / g.n[1], _exact_rational(st.shape.half_width)
    charge = zero(L)
    stiffness = zero(L)
    for ps in species
        q, m = _exact_rational(ps.q), _exact_rational(ps.m)
        for w in ps.weight
            weight = _exact_rational(w)
            charge += abs(q * weight)
            stiffness += q^2 * weight / m
        end
    end
    tau^2 * delta * min(inv(h), inv(a)) / h * stiffness < 4 ||
        throw(ArgumentError("dt violates the strict longitudinal rest-mode bound"))
    fieldbound = maximum(x -> abs(_exact_rational(x)), st.fields.E[1]) + charge
    for ps in species
        acceleration = abs(_exact_rational(ps.q) / _exact_rational(ps.m)) * (delta / h) * fieldbound
        for p in eachindex(ps.weight)
            x, v = _exact_rational(ps.x[1][p]), _exact_rational(ps.v[1][p])
            abs(v) * tau + acceleration * tau^2 / 2 < min(x, L - x) ||
                throw(ArgumentError("dt fails the strict longitudinal impact-exclusion tube"))
        end
    end
    return nothing
end

# Shared original-force owner. The caller prepares the original q/m and full-dt
# half once; gathered fields remain rational through the canonical Boris map.
function _yee_particle_half_kick(v::NTuple{3,T}, E::NTuple{3,Rational{BigInt}},
    B::NTuple{3,Rational{BigInt}}, factor::Rational{BigInt}, half::Rational{BigInt}) where {T}
    exact, _ = _boris_kick_parts(map(_exact_rational, v)..., E..., B..., factor, half)
    return map(value -> _from_exact(T, value), exact)
end
function _yl_kick(v::T, field::Rational{BigInt}, factor, half) where {T}
    z = zero(T)
    qz = zero(field)
    out = _yee_particle_half_kick((v, z, z), (field, qz, qz), (qz, qz, qz), factor, half)
    all(isfinite, out) && iszero(out[2]) && iszero(out[3]) ||
        throw(ArgumentError("longitudinal kick does not produce finite longitudinal velocity"))
    return out[1]
end

function step!(st::YeeLongitudinalPIC{T}, particles, dt::Real) where {T}
    _yl_require_mode()
    species = _yl_species(particles)
    tau = _require_finite_positive_real("dt", dt, T)
    st.initialized[] || throw(ArgumentError("call init! before longitudinal stepping"))
    _yl_validate_phase(st)
    st.step[] == 0 || isequal(tau, st.last_dt[]) ||
        throw(ArgumentError("longitudinal dt is fixed after the first accepted step"))
    next_time, next_step = _validated_step_clock(st.time[], st.step[], tau)
    next_time > st.time[] || throw(ArgumentError("longitudinal time must advance finitely"))
    _yl_validate_populations(st, species)
    _yl_validate_bound(st, species)
    _yl_validate_longitudinal_fields(st, st.fields)
    _yl_admit_step(st, species, tau)
    stage, g = st._stage[], st.g
    oldkinetic = _kinetic_energy_checked(species, T, Val(false))
    half = _exact_rational(tau) / 2
    factors = [_boris_exact_factor(_validated_boris_factor(ps)) for ps in species]
    for (ps, candidate, factor) in zip(species, stage.candidates, factors)
        for p in eachindex(ps.weight)
            x = ps.x[1][p]
            field = _yee_cloud_field_rational(st.fields.E[1], x, nothing, g, st.shape)
            u = _yl_kick(ps.v[1][p], field, factor, half)
            y = _particle_coordinate(x, u, tau)
            isfinite(y) && 0 < _exact_rational(y) < _exact_rational(g.L[1]) ||
                throw(ArgumentError("longitudinal drift leaves the strict interior"))
            candidate.x[1][p], candidate.v[1][p] = y, u
        end
    end
    old = [ps.x[1] for ps in species]
    new = [ps.x[1] for ps in stage.candidates]
    yee_longitudinal_impulse!(st._impulse[1], species, old, new, g,
        st._charge_workspace; shape = st.shape)
    # Only Ex is advanced. Structurally zero curl reuses the canonical source
    # row arithmetic without an extra Maxwell clock, transverse stage or CFL.
    for c = 2:3
        fill!(st._impulse[c], zero(T))
    end
    for a in st._candidate.B
        fill!(a, zero(T))
    end
    _yee_advance!(st._candidate.E, st.fields.E, st._candidate.B, g, false,
        tau, one(T), st._impulse)
    _yl_validate_longitudinal_fields(st, st._candidate)
    for (ps, candidate, factor) in zip(species, stage.candidates, factors)
        for p in eachindex(ps.weight)
            field = _yee_cloud_field_rational(st._candidate.E[1], candidate.x[1][p], nothing, g, st.shape)
            candidate.v[1][p] = _yl_kick(candidate.v[1][p], field, factor, half)
        end
    end
    newkinetic = kinetic_energy(stage.candidates)
    work = yee_source_work(st.fields.E, st._candidate.E, st._impulse, g)
    all(isfinite, (oldkinetic, newkinetic, work)) ||
        throw(ArgumentError("longitudinal step diagnostics do not round finite"))
    defect = _from_exact(T, _exact_rational(newkinetic) - _exact_rational(oldkinetic) - _exact_rational(work))
    isfinite(defect) || throw(ArgumentError("longitudinal energy defect does not round finite"))
    # Every potentially failing operation precedes the first public store.
    for (ps, candidate) in zip(species, stage.candidates)
        copyto!(ps.x[1], candidate.x[1])
        copyto!(ps.v[1], candidate.v[1])
    end
    copyto!(st.fields.E[1], st._candidate.E[1])
    for c = 1:3
        copyto!(st.last_impulse[c], st._impulse[c])
    end
    st.last_work[], st.last_energy_defect[] = work, defect
    st.time[], st.step[], st.last_dt[] = next_time, next_step, tau
    return st
end
