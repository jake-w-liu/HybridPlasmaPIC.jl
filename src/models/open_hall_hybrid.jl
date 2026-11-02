# Synchronous mixed-SBP/Fourier Hall family. Its numerical reservoir map is
# deliberately separate from periodic leapfrog and reduced PerpShock methods.
const _OpenHallQ = Rational{BigInt}

"One fixed inward Maxwellian numerical reservoir at one open x face."
struct OpenHallReservoir{T<:Union{Float32,Float64}}
    density::T
    drift::T
    spread::T
    tangent_drift::NTuple{2,T}
    tangent_spread::NTuple{2,T}
    weight::T
    tag::UInt32
end

function _open_hall_source(source::OpenHallReservoir{T}) where {T}
    all(isfinite, (source.density, source.drift, source.spread,
        source.tangent_drift..., source.tangent_spread..., source.weight)) ||
        throw(ArgumentError("open Hall source parameters must be finite"))
    source.density >= 0 && source.spread >= 0 && all(>=(zero(T)), source.tangent_spread) &&
        source.weight > 0 || throw(ArgumentError("open Hall source densities/spreads must be nonnegative and weight positive"))
    flux = flux_per_density(source.drift, source.spread)
    isfinite(flux) && flux >= 0 || throw(ArgumentError("open Hall source flux crosses its finite represented boundary"))
    return flux
end

function OpenHallReservoir(density::T, drift::T, spread::T;
    tangent_drift = (zero(T), zero(T)), tangent_spread = (zero(T), zero(T)),
    weight::T = one(T), tag::UInt32 = UInt32(0)) where {T<:Union{Float32,Float64}}
    tangent_drift isa NTuple{2,T} && tangent_spread isa NTuple{2,T} ||
        throw(ArgumentError("open Hall source records require one Float32/64 precision"))
    source = OpenHallReservoir{T}(density, drift, spread, tangent_drift, tangent_spread, weight, tag)
    _open_hall_source(source)
    return source
end

# Fixed-size signed aggregates, not an append-only event history.
struct _OpenHallExchange
    count::BigInt
    mass::_OpenHallQ
    charge::_OpenHallQ
    momentum::NTuple{3,_OpenHallQ}
    kinetic::_OpenHallQ
end
_OpenHallExchange() = _OpenHallExchange(BigInt(0), zero(_OpenHallQ), zero(_OpenHallQ),
    ntuple(_ -> zero(_OpenHallQ), 3), zero(_OpenHallQ))
function _open_hall_exchange_add(a::_OpenHallExchange, b::_OpenHallExchange, sign = 1)
    return _OpenHallExchange(a.count + sign * b.count, a.mass + sign * b.mass,
        a.charge + sign * b.charge,
        ntuple(c -> a.momentum[c] + sign * b.momentum[c], 3), a.kinetic + sign * b.kinetic)
end

struct _OpenHallFieldLedger
    delta_magnetic::_OpenHallQ
    raw_power::_OpenHallQ
    current_work::_OpenHallQ
    outward_poynting::_OpenHallQ
    raw_identity_defect::_OpenHallQ
    sat_work::_OpenHallQ
    assembly_work::_OpenHallQ
    rate_projection_work::_OpenHallQ
    rk_defect::_OpenHallQ
    endpoint_projection_work::_OpenHallQ
    stage_projection_energy::_OpenHallQ
end
_OpenHallFieldLedger() = _OpenHallFieldLedger(ntuple(_ -> zero(_OpenHallQ), 11)...)
_open_hall_field_add(a::_OpenHallFieldLedger, b::_OpenHallFieldLedger) =
    _OpenHallFieldLedger(ntuple(i -> getfield(a, i) + getfield(b, i), 11)...)

struct _OpenHallCumulative
    duration::_OpenHallQ
    field::_OpenHallFieldLedger
    kick_energy::_OpenHallQ
    kick_work::_OpenHallQ
    kick_momentum::NTuple{3,_OpenHallQ}
    event_energy::_OpenHallQ
    electron_enthalpy::NTuple{2,_OpenHallQ}
    energy_residual::NTuple{2,_OpenHallQ}
end
_OpenHallCumulative() = _OpenHallCumulative(zero(_OpenHallQ), _OpenHallFieldLedger(),
    zero(_OpenHallQ), zero(_OpenHallQ), ntuple(_ -> zero(_OpenHallQ), 3),
    zero(_OpenHallQ), (zero(_OpenHallQ), zero(_OpenHallQ)),
    (zero(_OpenHallQ), zero(_OpenHallQ)))

"Detached scalar record of one accepted OpenHallHybrid interval; no borrowed work arrays."
struct OpenHallStepReport{T}
    step::Int
    time::T
    dt::T
    level::Int
    segments::Int
    rhs_evaluations::Int
    field_projections::Int
    rate_projections::Int
    uniform_requests::Int
    births::Int
    exits::Int
    error::_OpenHallQ
    tolerance::_OpenHallQ
    worst_constraint::_OpenHallQ
    worst_equation::_OpenHallQ
    worst_work::_OpenHallQ
    field::_OpenHallFieldLedger
    kick_energy::_OpenHallQ
    kick_work::_OpenHallQ
    kick_momentum::NTuple{3,_OpenHallQ}
    event_energy::_OpenHallQ
    energy_before::NTuple{2,_OpenHallQ}
    energy_after::NTuple{2,_OpenHallQ}
    energy_change::NTuple{2,_OpenHallQ}
    electron_enthalpy::NTuple{2,_OpenHallQ}
    energy_residual::NTuple{2,_OpenHallQ}
    grid_charge_before::_OpenHallQ
    grid_charge_after::_OpenHallQ
    physical_charge_before::_OpenHallQ
    physical_charge_after::_OpenHallQ
    clock_rounding::_OpenHallQ
end

struct _OpenHallControls{T,P}
    field_rtol::T
    field_atol::T
    max_levels::Int
    max_births::Int
    max_particles::Int
    max_draws::Int
    log_bits::Int
    projection::P
end

mutable struct _OpenHallWorkspace{D,T,H,P}
    endpoint::HybridFields{D,T}
    force::HybridFields{D,T}
    hall::H
    projection::P
    base::NTuple{3,Array{T,D}}
    coarse::NTuple{3,Array{T,D}}
    half::NTuple{3,Array{T,D}}
    fine::NTuple{3,Array{T,D}}
    stage::NTuple{3,Array{T,D}}
    raw::NTuple{3,Array{T,D}}
    combined::NTuple{3,Array{T,D}}
    sat::NTuple{3,Array{T,D}}
    rates::NTuple{4,NTuple{3,Array{T,D}}}
    moment::NTuple{4,Array{_OpenHallQ,D}}
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    rhs_evaluations::Int
    field_projections::Int
    rate_projections::Int
    worst_constraint::_OpenHallQ
    worst_equation::_OpenHallQ
    worst_work::_OpenHallQ
end

mutable struct OpenHallHybrid{D,T,O,M,R,C,W,K} <: AbstractModel
    const op::O
    const model::M
    const shape::CIC
    const fields::HybridFields{D,T}
    const owners::Vector{_ElectrostaticHostParticle{D,T}}
    const _species_scalars::Vector{Tuple{T,T}}
    const _population_counts::Vector{Int}
    const namespaces::Vector{UInt64}
    const sources::Vector{NTuple{2,OpenHallReservoir{T}}}
    const fluxes::Vector{NTuple{2,T}}
    const targets::NTuple{2,NTuple{2,T}}
    const sat_speeds::NTuple{2,T}
    const controls::C
    const _configuration::K
    const rng::R
    const rng_candidate::R
    const alpha::Vector{_OpenHallQ}
    const next_index::Vector{UInt64}
    const initial_next_index::Vector{UInt64}
    const initial_totals::Vector{_OpenHallExchange}
    const births::Vector{_OpenHallExchange}
    const exits::Vector{_OpenHallExchange}
    const last_births::Vector{_OpenHallExchange}
    const last_exits::Vector{_OpenHallExchange}
    source_elapsed::_OpenHallQ
    cumulative::_OpenHallCumulative
    const initialized::Base.RefValue{Bool}
    const time::Base.RefValue{T}
    const step::Base.RefValue{Int}
    const last_dt::Base.RefValue{T}
    last_report::Union{Nothing,OpenHallStepReport{T}}
    const work::W
    busy::Bool
end

"Bounded private field trials failed; no physical/source/RNG state was published."
struct OpenHallControllerError <: Exception
    level::Int
    reason::Symbol
    detail::Any
end
Base.showerror(io::IO, e::OpenHallControllerError) =
    print(io, "open Hall field controller rejected level ", e.level, " (", e.reason, "): ", e.detail)
struct _OpenHallRetry <: Exception
    reason::Symbol
    detail::Any
end

function _open_hall_limit(name, value)
    value isa Integer && 0 <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be a nonnegative Int-sized integer"))
    return Int(value)
end

function _open_hall_controls(::Type{T}; max_levels, max_births, max_particles, max_draws,
    field_rtol = sqrt(eps(T)), field_atol = zero(T), projection_rtol = 100eps(T),
    projection_atol = zero(T), state_rtol = projection_rtol, state_atol = projection_atol,
    work_rtol = 100eps(T), work_atol = zero(T), maxiter = 200,
    log_bits = precision(T) + 16) where {T}
    r = _require_finite_nonnegative_real("field_rtol", field_rtol, T)
    a = _require_finite_nonnegative_real("field_atol", field_atol, T)
    r < 1 && (r > 0 || a > 0) && (iszero(field_rtol) || r > 0) &&
        (iszero(field_atol) || a > 0) || throw(ArgumentError("open Hall field tolerances must survive conversion and define a positive gate"))
    levels = _open_hall_limit("max_levels", max_levels)
    # All cost counters/level shifts must fit their concrete Int storage.
    levels < 8sizeof(Int) - 2 || throw(ArgumentError("open Hall max_levels exceeds checked segment-count range"))
    trials = Base.checked_sub(Base.checked_mul(2, Int(1) << levels), 1)
    Base.checked_mul(12, trials)
    bits = _open_hall_limit("log_bits", log_bits)
    precision(T) + 16 <= bits <= 4096 || throw(ArgumentError("open Hall log_bits must lie between precision(T)+16 and 4096"))
    projection = _mixed_projection_controls(T, projection_rtol, projection_atol,
        work_rtol, work_atol, state_rtol, state_atol, maxiter)
    for (old, converted) in ((projection_rtol, projection.rtol), (projection_atol, projection.atol),
        (state_rtol, projection.state_rtol), (state_atol, projection.state_atol),
        (work_rtol, projection.work_rtol), (work_atol, projection.work_atol))
        iszero(old) || converted != 0 || throw(ArgumentError("open Hall projection tolerance disappears at grid precision"))
    end
    return _OpenHallControls(r, a, levels, _open_hall_limit("max_births", max_births),
        _open_hall_limit("max_particles", max_particles), _open_hall_limit("max_draws", max_draws), bits, projection)
end

function _open_hall_physics(model::HybridModel, ::Type{T}) where {T}
    model.closure isa IsothermalElectrons{T} || throw(ArgumentError("open Hall requires grid-precision isothermal electrons"))
    _require_finite_nonnegative_real("open Hall Te", model.closure.Te, T)
    all(_particle_stored_zero, (model.η, model.ηH, model.de2)) || throw(ArgumentError("open Hall requires original eta=etaH=de2=0"))
    return _require_finite_positive_real("open Hall nfloor", model.nfloor, T)
end

_open_hall_field_arrays(f::HybridFields) =
    (f.n, f.ui..., f.B..., f.E..., f.J..., f.pe, f.gradp..., f.ninv,
        f.lapJ..., f.pforce..., f.cgl_dp)
_open_hall_particle_arrays(species) =
    Tuple(a for ps in species for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
function _open_hall_work_arrays(w::_OpenHallWorkspace)
    return (_open_hall_field_arrays(w.endpoint)..., _open_hall_field_arrays(w.force)...,
        _mixed_hall_arrays(w.hall)..., _mixed_projection_arrays(w.projection)...,
        w.base..., w.coarse..., w.half..., w.fine..., w.stage..., w.raw...,
        w.combined..., w.sat..., (a for v in w.rates for a in v)..., w.moment..., w.Ep..., w.Bp...)
end
function _open_hall_owned_arrays(st::OpenHallHybrid)
    return (_open_hall_field_arrays(st.fields)..., st.owners, st._species_scalars,
        st._population_counts, st.namespaces, st.sources, st.fluxes, st.alpha,
        st.next_index, st.initial_next_index, st.initial_totals, st.births, st.exits,
        st.last_births, st.last_exits,
        _open_hall_particle_arrays(st.owners)..., _open_hall_work_arrays(st.work)...,
        _mixed_operator_arrays(st.op)..., _collisional_rng_buffers(st.rng)...,
        _collisional_rng_buffers(st.rng_candidate)...)
end

function _open_hall_configuration(op, model, scalars, namespaces, sources, targets, speeds, controls)
    return (family = :open_hall_hybrid, dimension = length(op.n), precision = eltype(op.temp),
        n = op.n, L = op.L, dx = op.dx, sbp_policy = :sbp_2_1, derivative_policy = :represented_fourier_odd_nyquist_zero,
        model = model, species_scalars = Tuple(scalars), namespaces = Tuple(namespaces),
        sources = Tuple(sources), targets = targets, sat_speeds = speeds, controls = controls,
        phase = :synchronous_integer_lie_fields_velocity_ballistic_v1,
        rk = :exact_original_duration_rk4_v1, moments = :exact_sbp_charge_weighted_cic_v1,
        gather = :complete_original_cic_v1, drift = :complete_ballistic_first_exit_v1,
        fields = :complete_public_hybrid_v1, source_policy = :bounded_stock_uniform_requests_v1,
        source_order = :species_left_right_birth, event_equality = :exit_at_endpoint,
        report_policy = :fixed_log_interval_signed_work_v1)
end
function _open_hall_bound_configuration(st::OpenHallHybrid{D,T}) where {D,T}
    _mixed_storage(st.op)
    _open_hall_physics(st.model, T)
    c = st.controls
    _open_hall_controls(T; max_levels = c.max_levels, max_births = c.max_births,
        max_particles = c.max_particles, max_draws = c.max_draws,
        field_rtol = c.field_rtol, field_atol = c.field_atol,
        projection_rtol = c.projection.rtol, projection_atol = c.projection.atol,
        state_rtol = c.projection.state_rtol, state_atol = c.projection.state_atol,
        work_rtol = c.projection.work_rtol, work_atol = c.projection.work_atol,
        maxiter = c.projection.maxiter, log_bits = c.log_bits)
    current = _open_hall_configuration(st.op, st.model, st._species_scalars,
        st.namespaces, st.sources, st.targets, st.sat_speeds, st.controls)
    isequal(current, st._configuration) || throw(ArgumentError("open Hall immutable configuration/cache binding changed"))
    for i in eachindex(st.sources)
        map(_open_hall_source, st.sources[i]) == st.fluxes[i] || throw(ArgumentError("open Hall cached source flux changed"))
    end
    return current
end

function _open_hall_storage(st::OpenHallHybrid{D,T}; rng_metadata = true) where {D,T}
    ns = length(st.owners)
    ns > 0 && ns <= 65536 || throw(ArgumentError("open Hall requires 1..65536 ordered species"))
    for a in (st._species_scalars, st._population_counts, st.namespaces, st.sources,
        st.fluxes, st.next_index, st.initial_next_index, st.initial_totals)
        axes(a, 1) == Base.OneTo(ns) || throw(DimensionMismatch("open Hall binding/cache length differs"))
    end
    for a in (st.alpha, st.births, st.exits, st.last_births, st.last_exits)
        axes(a, 1) == Base.OneTo(2ns) || throw(DimensionMismatch("open Hall face-source cache length differs"))
    end
    length(unique(st.namespaces)) == ns && all(x -> x <= _PARTICLE_ID_MAX_SPECIES, st.namespaces) ||
        throw(ArgumentError("open Hall species namespaces must be distinct 16-bit values"))
    for f in (st.fields, st.work.endpoint, st.work.force)
        for (name, value) in pairs(_hybrid_field_state(f)), a in (value isa Tuple ? value : (value,))
            a isa Array{T,D} && size(a) == (name in (:pforce, :cgl_dp) ? ntuple(_ -> 0, D) : st.op.n) ||
                throw(DimensionMismatch("open Hall field role differs from its ordinary mixed-grid storage"))
        end
    end
    _mixed_projection_storage(st.work.projection, st.op)
    st.work.hall.n == st.op.n || throw(DimensionMismatch("open Hall electric workspace shape changed"))
    for a in (_mixed_hall_arrays(st.work.hall)..., st.work.base..., st.work.coarse...,
        st.work.half..., st.work.fine..., st.work.stage..., st.work.raw...,
        st.work.combined..., st.work.sat..., (a for v in st.work.rates for a in v)...)
        _mixed_array(a, st.op, "open Hall private field")
    end
    for a in st.work.moment
        size(a) == st.op.n || throw(DimensionMismatch("open Hall exact moment work differs"))
    end
    all(a -> a isa Vector{T}, (st.work.Ep..., st.work.Bp...)) || throw(ArgumentError("open Hall gather work requires Vectors"))
    refs = (st.initialized, st.time, st.step, st.last_dt, st.fields.floor_count,
        st.work.endpoint.floor_count, st.work.force.floor_count, st.work.hall.floor_count)
    for i in eachindex(refs), j = 1:(i-1)
        refs[i] === refs[j] && throw(ArgumentError("open Hall phase/field references must be independent"))
    end
    _validate_collisional_rng_storage(st.rng)
    _validate_collisional_rng_storage(st.rng_candidate)
    st.rng === st.rng_candidate && throw(ArgumentError("open Hall live/candidate RNG identities must differ"))
    arrays = _open_hall_owned_arrays(st)
    for i in eachindex(arrays), j = 1:(i-1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(ArgumentError("open Hall public/private storage roles must be independent"))
    end
    if rng_metadata
        for metadata in (_collisional_rng_metadata(st.rng), _collisional_rng_metadata(st.rng_candidate))
            _validate_checkpoint_configuration_registered(metadata) do remember, backing, owned
                _open_hall_register_owned!(remember, backing, owned, st)
            end
        end
    end
    return nothing
end

function _open_hall_binding(st::OpenHallHybrid{D,T}, species; values = true, counts = values) where {D,T}
    _open_hall_bound_configuration(st)
    _open_hall_storage(st; rng_metadata = values)
    species isa Vector{<:ParticleSet{D,T}} && length(species) == length(st.owners) || throw(ArgumentError("open Hall requires its bound ordered ordinary population Vector"))
    for (s, ps) in enumerate(species)
        ps === st.owners[s] && ps isa _ElectrostaticHostParticle{D,T} || throw(ArgumentError("open Hall particle owner identity/order/storage changed"))
        _check_particle_storage_axes(ps)
        q, m = st._species_scalars[s]
        _validated_charge_mass(q, m, T)
        _stored_positive(q) || throw(ArgumentError("open Hall cached ion charge must be positive"))
        values && !isequal((q, m), (ps.q, ps.m)) && throw(ArgumentError("open Hall species charge/mass binding changed"))
        counts && st._population_counts[s] != nparticles(ps) && throw(ArgumentError("open Hall counts changed outside prepared model publication"))
    end
    values && _open_hall_particles(st, species)
    return nothing
end

function _open_hall_particles(st::OpenHallHybrid{D,T}, species) where {D,T}
    total = _validate_hybrid_species(species)
    total <= st.controls.max_particles || throw(ArgumentError("open Hall population limit exceeded"))
    for (s, ps) in enumerate(species)
        seen = Set{UInt64}()
        for p in eachindex(ps.weight)
            id = ps.id[p]
            id >> _PARTICLE_ID_INDEX_BITS == st.namespaces[s] &&
                (id & _PARTICLE_ID_MAX_INDEX) != 0 && !(id in seen) ||
                throw(ArgumentError("open Hall IDs require unique nonzero indices in the bound species namespace"))
            push!(seen, id)
            for d = 1:D
                x = ps.x[d][p]
                (d == 1 ? 0 <= x <= st.op.L[d] : 0 <= x < st.op.L[d]) ||
                    throw(ArgumentError("open Hall positions must lie in the closed x / periodic transverse domain"))
            end
        end
    end
    return total
end

function _open_hall_fields(st::OpenHallHybrid{D,T}, f::HybridFields{D,T}) where {D,T}
    all(a -> all(isfinite, a), _open_hall_field_arrays(f)) || throw(ArgumentError("open Hall physical fields must be finite"))
    floor = _open_hall_physics(st.model, T)
    all(>=(floor), f.n) && all(>=(zero(T)), f.pe) && all(>(zero(T)), f.ninv) &&
        f.floor_count[] == 0 && all(a -> all(_particle_stored_zero, a), f.lapJ) ||
        throw(ArgumentError("open Hall requires strict density, positive reciprocal, nonnegative pressure and inactive zero auxiliary fields"))
    for side in (1, st.op.n[1]), I in CartesianIndices(Base.tail(st.op.n))
        _particle_stored_zero(f.B[1][CartesianIndex((side, Tuple(I)...))]) ||
            throw(ArgumentError("open Hall normal magnetic faces must be stored zero"))
    end
    return nothing
end

function _open_hall_phase_values(st::OpenHallHybrid{D,T}, initialized, time, step, last_dt, report) where {D,T}
    initialized isa Bool && time isa T && step isa Int && last_dt isa T &&
        isfinite(time) && time >= 0 && step >= 0 && isfinite(last_dt) && last_dt >= 0 ||
        throw(ArgumentError("open Hall phase/clock types or values are invalid"))
    if step == 0
        _particle_stored_zero(time) && _particle_stored_zero(last_dt) && report === nothing ||
            throw(ArgumentError("unstepped open Hall requires zero clock/history and no report"))
    else
        initialized && time > 0 && last_dt > 0 && last_dt <= time || throw(ArgumentError("advanced open Hall requires initialized positive clock/history"))
        _open_hall_pair_duration(last_dt, 0)
        _open_hall_report(st, report, time, step, last_dt)
    end
    return nothing
end

function _open_hall_direct_state(st, B)
    c = st.controls.projection
    source = _exact_rational(_mixed_vector_maxabs(B))
    result = _mixed_direct_constraints(B, st.op, st.work.projection, c.state_rtol, c.state_atol, source)
    result.passed || throw(MixedHallProjectionError("open Hall physical B state constraint gate", result))
    return result
end

function _open_hall_boundary(st::OpenHallHybrid, species = st.owners; differentiated = true)
    _open_hall_binding(st, species)
    st.initialized[] || throw(ArgumentError("open Hall model is not initialized"))
    _open_hall_phase_values(st, true, st.time[], st.step[], st.last_dt[], st.last_report)
    _open_hall_fields(st, st.fields)
    _validate_collisional_rng(st.rng)
    _open_hall_source_values(st, _open_hall_source_state(st), species, st.step[]; report = st.last_report)
    differentiated && _open_hall_direct_state(st, st.fields.B)
    return nothing
end

function _acquire_open_hall!(st::OpenHallHybrid)
    st.busy && throw(ArgumentError("open Hall model permits one nonreentrant caller"))
    st.busy = true
    return nothing
end

# Own metadata is compared by value after Serialization/deepcopy, including the
# Q-valued controls; default object identity is not a restart binding policy.
Base.isequal(a::_OpenHallControls, b::_OpenHallControls) = typeof(a) === typeof(b) &&
    isequal(ntuple(i -> getfield(a, i), fieldcount(typeof(a))), ntuple(i -> getfield(b, i), fieldcount(typeof(b))))
_open_hall_exchange_equal(a::_OpenHallExchange, b::_OpenHallExchange) =
    a.count == b.count && a.mass == b.mass && a.charge == b.charge &&
    a.momentum == b.momentum && a.kinetic == b.kinetic

function _open_hall_register_owned!(remember, backing, owned, st::OpenHallHybrid)
    # Register writable roles, not the stock metadata being checked. Its graph
    # may contain arbitrary user seed objects and cannot be treated as an Array.
    owned[st] = true
    for object in (st.op, st.fields, st.work, st.work.endpoint, st.work.force,
        st.work.hall, st.work.projection, st.work.projection.solver,
        st.rng, st.rng_candidate, st.owners...)
        remember(object)
    end
    for spectral in st.op.transverse
        _visit_checkpoint_objects(remember, spectral; backing)
    end
    for a in _open_hall_owned_arrays(st)
        _visit_checkpoint_objects(remember, a; backing)
    end
    _visit_checkpoint_objects(remember, _open_hall_source_state(st); backing)
    _visit_checkpoint_objects(remember, st.last_report; backing)
    for ref in (st.initialized, st.time, st.step, st.last_dt, st.fields.floor_count,
        st.work.endpoint.floor_count, st.work.force.floor_count, st.work.hall.floor_count)
        _visit_checkpoint_objects(remember, ref; backing)
    end
    return nothing
end

function _open_hall_particle_totals(species)
    result = _OpenHallExchange()
    for ps in species, p in eachindex(ps.weight)
        Q = _exact_rational
        mass, charge = Q(ps.m) * Q(ps.weight[p]), Q(ps.q) * Q(ps.weight[p])
        velocity = ntuple(c -> Q(ps.v[c][p]), 3)
        event = _OpenHallExchange(BigInt(1), mass, charge,
            ntuple(c -> mass * velocity[c], 3), zero(_OpenHallQ))
        result = _open_hall_exchange_add(result, event)
    end
    return _OpenHallExchange(result.count, result.mass, result.charge, result.momentum,
        _kinetic_energy_rational(species))
end
# Q is a concrete scalar carrier, but Rational{BigInt}(n,0) can represent
# infinities. Validate every leaf of our known report/ledger records explicitly;
# this is not opaque metadata traversal or a serialization-alias policy.
function _open_hall_finite_scalars(value::_OpenHallQ)
    denominator(value) > 0 || throw(ArgumentError("open Hall exact scalar must be finite"))
    return nothing
end
function _open_hall_finite_scalars(value::Real)
    isfinite(value) || throw(ArgumentError("open Hall scalar record must be finite"))
    return nothing
end
function _open_hall_finite_scalars(values::Tuple)
    foreach(_open_hall_finite_scalars, values)
    return nothing
end
function _open_hall_finite_scalars(record::Union{_OpenHallExchange,_OpenHallFieldLedger,_OpenHallCumulative,OpenHallStepReport})
    for i = 1:fieldcount(typeof(record))
        _open_hall_finite_scalars(getfield(record, i))
    end
    return nothing
end

function _open_hall_exchange_valid(a::_OpenHallExchange)
    _open_hall_finite_scalars(a)
    a.count >= 0 && a.mass >= 0 && a.charge >= 0 && a.kinetic >= 0 &&
        sum(abs2, a.momentum) <= 2 * a.mass * a.kinetic ||
        throw(ArgumentError("open Hall exchange has invalid count/mass/charge/kinetic or momentum bound"))
    (a.mass != 0 || a.kinetic == 0) && (a.count != 0 ||
        (a.mass == 0 && a.charge == 0 && a.kinetic == 0 && all(iszero, a.momentum))) ||
        throw(ArgumentError("open Hall zero-count/zero-mass exchange cannot carry physical event energy"))
    return nothing
end
function _open_hall_field_ledger_valid(a::_OpenHallFieldLedger)
    _open_hall_finite_scalars(a)
    a.raw_power + a.current_work + a.outward_poynting == a.raw_identity_defect &&
        a.delta_magnetic == a.raw_power + a.sat_work + a.assembly_work +
            a.rate_projection_work + a.rk_defect + a.endpoint_projection_work ||
        throw(ArgumentError("open Hall magnetic ledger identities are inconsistent"))
    return nothing
end

function _open_hall_report(st::OpenHallHybrid{D,T}, r, time, step, h) where {D,T}
    r isa OpenHallStepReport{T} && r.step == step && isequal(r.time, time) && isequal(r.dt, h) ||
        throw(ArgumentError("open Hall report must match its completed clock/history"))
    _open_hall_finite_scalars(r)
    c = st.controls
    0 <= r.level <= c.max_levels && r.segments == Int(1) << r.level &&
        0 <= r.births <= c.max_births && 0 <= r.exits <= c.max_particles &&
        0 <= r.uniform_requests <= c.max_draws || throw(ArgumentError("open Hall report counters exceed its fixed limits"))
    _open_hall_pair_duration(h, r.level)
    maximum_calls = 12 * (2r.segments - 1)
    8r.segments <= r.rhs_evaluations <= maximum_calls &&
        8r.segments <= r.field_projections <= maximum_calls &&
        8r.segments <= r.rate_projections <= maximum_calls &&
        0 <= r.error <= r.tolerance &&
        min(r.worst_constraint, r.worst_equation, r.worst_work) >= 0 ||
        throw(ArgumentError("open Hall report has invalid field-controller metrics"))
    for interval in (r.energy_before, r.energy_after, r.energy_change, r.electron_enthalpy, r.energy_residual)
        interval[1] <= interval[2] || throw(ArgumentError("open Hall energy interval endpoints are reversed"))
    end
    _open_hall_field_ledger_valid(r.field)
    r.energy_after[1] - r.energy_before[2] <= r.energy_change[1] <= r.energy_change[2] <=
        r.energy_after[2] - r.energy_before[1] || throw(ArgumentError("open Hall merged energy-change interval is inconsistent with endpoint bounds"))
    expected = (r.energy_change[1] - r.event_energy + r.field.outward_poynting + r.electron_enthalpy[1],
        r.energy_change[2] - r.event_energy + r.field.outward_poynting + r.electron_enthalpy[2])
    r.energy_residual == expected || throw(ArgumentError("open Hall report residual differs from its signed reservoir definition"))
    return nothing
end

_open_hall_source_state(st::OpenHallHybrid) =
    (alpha = st.alpha, next_index = st.next_index, initial_next_index = st.initial_next_index,
        initial_totals = st.initial_totals, births = st.births, exits = st.exits,
        last_births = st.last_births, last_exits = st.last_exits,
        source_elapsed = st.source_elapsed, cumulative = st.cumulative)

function _open_hall_source_values(st::OpenHallHybrid{D,T}, source, species, step; report = nothing) where {D,T}
    source isa NamedTuple && keys(source) == (:alpha, :next_index, :initial_next_index,
        :initial_totals, :births, :exits, :last_births, :last_exits, :source_elapsed, :cumulative) ||
        throw(ArgumentError("open Hall source state has an unsupported layout"))
    ns = length(st.owners)
    source.alpha isa Vector{_OpenHallQ} && length(source.alpha) == 2ns &&
        source.next_index isa Vector{UInt64} && length(source.next_index) == ns &&
        source.initial_next_index isa Vector{UInt64} && length(source.initial_next_index) == ns &&
        source.initial_totals isa Vector{_OpenHallExchange} && length(source.initial_totals) == ns &&
        source.births isa Vector{_OpenHallExchange} && length(source.births) == 2ns &&
        source.exits isa Vector{_OpenHallExchange} && length(source.exits) == 2ns &&
        source.last_births isa Vector{_OpenHallExchange} && length(source.last_births) == 2ns &&
        source.last_exits isa Vector{_OpenHallExchange} && length(source.last_exits) == 2ns &&
        source.source_elapsed isa _OpenHallQ && source.source_elapsed >= 0 &&
        source.cumulative isa _OpenHallCumulative || throw(ArgumentError("open Hall source storage/types differ"))
    foreach(_open_hall_finite_scalars, source.alpha)
    _open_hall_finite_scalars(source.source_elapsed)
    _open_hall_finite_scalars(source.cumulative)
    source.cumulative.duration == source.source_elapsed || throw(ArgumentError("open Hall exact source/controller duration differs"))
    for event in (source.initial_totals..., source.births..., source.exits..., source.last_births..., source.last_exits...)
        _open_hall_exchange_valid(event)
    end
    for s = 1:ns
        initial, next = source.initial_next_index[s], source.next_index[s]
        1 <= initial <= next <= _PARTICLE_ID_MAX_INDEX + UInt64(1) ||
            throw(ArgumentError("open Hall provenance counter/exhausted sentinel is invalid"))
        source.initial_totals[s].count <= BigInt(initial) - 1 ||
            throw(ArgumentError("open Hall initial population cannot fit below its initial next index"))
        produced = BigInt(0)
        for face = 1:2
            i = 2(s - 1) + face
            a, birth, exit = source.alpha[i], source.births[i], source.exits[i]
            0 <= a < 1 || throw(ArgumentError("open Hall source fractional carry must lie in[0,1)"))
            last_birth, last_exit = source.last_births[i], source.last_exits[i]
            _open_hall_exchange_valid(_open_hall_exchange_add(birth, last_birth, -1))
            _open_hall_exchange_valid(_open_hall_exchange_add(exit, last_exit, -1))
            record = st.sources[s][face]
            rate = _exact_rational(record.density) * _exact_rational(st.fluxes[s][face]) *
                prod(_exact_rational, Base.tail(st.op.L)) / _exact_rational(record.weight)
            a + birth.count == rate * source.source_elapsed ||
                throw(ArgumentError("open Hall fractional source/count/duration identity is inconsistent"))
            for event in (birth, last_birth)
                event.mass == event.count * _exact_rational(record.weight) * _exact_rational(st._species_scalars[s][2]) &&
                    event.charge == event.count * _exact_rational(record.weight) * _exact_rational(st._species_scalars[s][1]) ||
                    throw(ArgumentError("open Hall birth ledger differs from its fixed source weight/factors"))
            end
            for event in (source.initial_totals[s], birth, exit, last_birth, last_exit)
                event.mass * _exact_rational(st._species_scalars[s][1]) ==
                    event.charge * _exact_rational(st._species_scalars[s][2]) ||
                    throw(ArgumentError("open Hall mass/charge ledger differs from cached original factors"))
            end
            produced += birth.count
        end
        BigInt(next) == BigInt(initial) + produced || throw(ArgumentError("open Hall next index differs from reserved births"))
        all(id -> (id & _PARTICLE_ID_MAX_INDEX) < next, species[s].id) ||
            throw(ArgumentError("open Hall retained ID is outside the consumed provenance range"))
    end
    initial = foldl(_open_hall_exchange_add, source.initial_totals; init = _OpenHallExchange())
    births = foldl(_open_hall_exchange_add, source.births; init = _OpenHallExchange())
    exits = foldl(_open_hall_exchange_add, source.exits; init = _OpenHallExchange())
    initial.count <= st.controls.max_particles && births.count <= BigInt(step) * st.controls.max_births &&
        exits.count <= initial.count + births.count || throw(ArgumentError("open Hall cumulative source counts exceed their fixed policy"))
    last_births = foldl(_open_hall_exchange_add, source.last_births; init = _OpenHallExchange())
    last_exits = foldl(_open_hall_exchange_add, source.last_exits; init = _OpenHallExchange())
    last_births.count <= st.controls.max_births && last_exits.count <= st.controls.max_particles ||
        throw(ArgumentError("open Hall last event counts exceed fixed per-step limits"))
    if step == 0
        last_births.count == 0 && last_exits.count == 0 ||
            throw(ArgumentError("unstepped open Hall cannot carry last-step events"))
    else
        report isa OpenHallStepReport{T} && report.births == last_births.count &&
            report.exits == last_exits.count && report.event_energy == last_births.kinetic - last_exits.kinetic ||
            throw(ArgumentError("open Hall last per-face events differ from the completed scalar report"))
    end
    current = _open_hall_particle_totals(species)
    expected = _open_hall_exchange_add(_open_hall_exchange_add(initial, births), exits, -1)
    cumulative = source.cumulative
    current.count == expected.count && current.mass == expected.mass && current.charge == expected.charge &&
        current.momentum == ntuple(c -> expected.momentum[c] + cumulative.kick_momentum[c], 3) &&
        current.kinetic == expected.kinetic + cumulative.kick_energy &&
        cumulative.event_energy == births.kinetic - exits.kinetic ||
        throw(ArgumentError("open Hall retained-record/source/kick ledger identities differ"))
    _open_hall_field_ledger_valid(cumulative.field)
    cumulative.electron_enthalpy[1] <= cumulative.electron_enthalpy[2] &&
        cumulative.energy_residual[1] <= cumulative.energy_residual[2] ||
        throw(ArgumentError("open Hall cumulative interval endpoints are invalid"))
    if step == 0
        for s = 1:ns
            _open_hall_exchange_equal(source.initial_totals[s], _open_hall_particle_totals((species[s],))) ||
                throw(ArgumentError("unstepped open Hall initial particle totals differ"))
            highest = maximum(id -> id & _PARTICLE_ID_MAX_INDEX, species[s].id; init = UInt64(0))
            source.initial_next_index[s] == highest + UInt64(1) ||
                throw(ArgumentError("unstepped open Hall initial index does not follow its actual IDs"))
        end
        source.source_elapsed == 0 && births.count == 0 && exits.count == 0 &&
            all(iszero, source.alpha) && cumulative.kick_energy == 0 && cumulative.kick_work == 0 &&
            all(iszero, cumulative.kick_momentum) && cumulative.event_energy == 0 &&
            all(iszero, cumulative.electron_enthalpy) && all(iszero, cumulative.energy_residual) &&
            all(i -> iszero(getfield(cumulative.field, i)), 1:11) ||
            throw(ArgumentError("unstepped open Hall source/ledgers must be zero"))
    else
        source.source_elapsed > 0 || throw(ArgumentError("advanced open Hall requires positive accumulated exact duration"))
    end
    return nothing
end

function _open_hall_replacements(st::OpenHallHybrid{D,T}, incoming) where {D,T}
    prepared = _ElectrostaticHostParticle{D,T}[]
    counts = Int[]
    for (s, ps) in enumerate(incoming)
        isequal((ps.q, ps.m), st._species_scalars[s]) || throw(ArgumentError("open Hall replacement original factors differ"))
        push!(prepared, _prepare_particle_replacement(ps, collect(eachindex(ps.weight))))
        push!(counts, nparticles(ps))
    end
    return prepared, counts
end
function _open_hall_reserve!(st::OpenHallHybrid, species)
    largest = maximum(nparticles, species)
    for a in (st.work.Ep..., st.work.Bp...)
        resize!(a, largest)
    end
    return nothing
end

# Inputs have passed all graph/storage/config/value checks and allocations.
function _open_hall_commit!(st::OpenHallHybrid, prepared, counts, fields, source, rng,
    initialized, time, step, last_dt, report)
    for s in eachindex(st.owners)
        _publish_particle_replacement!(st.owners[s], prepared[s])
        st.owners[s].q, st.owners[s].m = st._species_scalars[s]
    end
    _copy_hybrid_fields!(st.fields, fields)
    copyto!(st._population_counts, counts)
    for name in (:alpha, :next_index, :initial_next_index, :initial_totals, :births, :exits, :last_births, :last_exits)
        copyto!(getproperty(st, name), getproperty(source, name))
    end
    st.source_elapsed = source.source_elapsed
    st.cumulative = source.cumulative
    rng === st.rng_candidate || copy!(st.rng_candidate, rng)
    rng === st.rng || copy!(st.rng, rng)
    st.last_report = report
    st.last_dt[] = last_dt
    st.time[] = time
    st.step[] = step
    st.initialized[] = initialized
    return st
end

function _open_hall_workspace(op::MixedShockOperators{D,T}) where {D,T}
    A() = zeros(T, op.n)
    V() = (A(), A(), A())
    return _OpenHallWorkspace(HybridFields{D,T}(op.n), HybridFields{D,T}(op.n),
        MixedHallWorkspace(op), MixedHallProjectionWorkspace(op), V(), V(), V(), V(),
        V(), V(), V(), V(), ntuple(_ -> V(), 4), ntuple(_ -> zeros(_OpenHallQ, op.n), 4),
        ntuple(_ -> T[], 3), ntuple(_ -> T[], 3), 0, 0, 0,
        zero(_OpenHallQ), zero(_OpenHallQ), zero(_OpenHallQ))
end

"""
    OpenHallHybrid(op, physics, species, sources; rng, max_levels, max_births,
        max_particles, max_draws, targets, sat_speeds, B=nothing, ...)

CPU D2/3 positive-ion CIC/isothermal synchronous Hall numerical reservoir model.
The intended perpendicular space has Bx=Jx=0 at both x faces; represented
residual gates apply.
Transverse axes are periodic. This is a restrictive
numerical boundary model, not transparent Hall outflow or a stability theorem.
`species` must have ordinary Vector leaves and distinct 16/48 provenance namespaces;
`sources[s]` is a same-precision `(left,right)` pair of OpenHallReservoir records.
The supplied stock RNG is transferred to model ownership; do not draw externally.
Call init! once after filling particle values. See docs/OPEN_HALL_HYBRID.md.
"""
function OpenHallHybrid(op::MixedShockOperators{D,T}, physics::HybridModel,
    species::AbstractVector, sources::AbstractVector;
    rng, max_levels, max_births, max_particles, max_draws,
    namespaces = collect(0:(length(species)-1)),
    targets = ((zero(T), zero(T)), (zero(T), zero(T))),
    sat_speeds = (zero(T), zero(T)), B = nothing, kwargs...) where {D,T}
    _mixed_storage(op)
    _open_hall_physics(physics, T)
    1 <= length(species) <= 65536 && all(ps -> ps isa _ElectrostaticHostParticle{D,T}, species) ||
        throw(ArgumentError("open Hall requires 1..65536 same-grid-precision Vector-backed populations"))
    length(sources) == length(species) && all(s -> s isa NTuple{2,OpenHallReservoir{T}}, sources) ||
        throw(ArgumentError("open Hall requires one fixed same-precision source pair per species"))
    targets isa NTuple{2,NTuple{2,T}} && sat_speeds isa NTuple{2,T} &&
        all(isfinite, (targets[1]..., targets[2]..., sat_speeds...)) && all(>=(zero(T)), sat_speeds) ||
        throw(ArgumentError("open Hall tangent targets/speeds must be finite at grid precision; speeds nonnegative"))
    length(namespaces) == length(species) || throw(DimensionMismatch("open Hall namespace count differs"))
    ids = UInt64[_particle_species(s) for s in namespaces]
    length(unique(ids)) == length(ids) || throw(ArgumentError("open Hall requires distinct species namespaces"))
    owners = _ElectrostaticHostParticle{D,T}[ps for ps in species]
    scalars = Tuple{T,T}[(ps.q, ps.m) for ps in species]
    for (q, m) in scalars
        _validated_charge_mass(q, m, T)
        _stored_positive(q) || throw(ArgumentError("open Hall requires positive original ion charge"))
    end
    controls = _open_hall_controls(T; max_levels, max_births, max_particles, max_draws, kwargs...)
    counts = Int[nparticles(ps) for ps in species]
    sum(BigInt, counts) <= controls.max_particles || throw(ArgumentError("open Hall initial population exceeds max_particles"))
    records = NTuple{2,OpenHallReservoir{T}}[pair for pair in sources]
    fluxes = NTuple{2,T}[map(_open_hall_source, pair) for pair in records]
    rng isa Union{Random.Xoshiro,Random.MersenneTwister} || throw(ArgumentError("open Hall requires Xoshiro or MersenneTwister"))
    _validate_collisional_rng(rng)
    candidate_rng = deepcopy(rng)
    _validate_collisional_rng_pair(rng, candidate_rng, "open Hall")
    fields = HybridFields{D,T}(op.n)
    if B === nothing
        Q = _exact_rational
        for c = 1:3, I in CartesianIndices(fields.n)
            fraction = Q(op.coordinates[1][I[1]]) / Q(op.L[1])
            fields.B[c][I] = c == 1 ? zero(T) : _from_exact(T,
                (1 - fraction) * Q(targets[1][c - 1]) + fraction * Q(targets[2][c - 1]))
        end
    else
        B isa NTuple{3,Array{T,D}} || throw(ArgumentError("open Hall seed B requires matching ordinary Arrays"))
        for c = 1:3
            _mixed_array(B[c], op, "open Hall magnetic seed"; finite = true)
            copyto!(fields.B[c], B[c])
        end
    end
    # Q-valued controls contain BigInts; cache their values independently rather
    # than letting a later mutable scalar edit change both sides of the binding.
    configuration = deepcopy(_open_hall_configuration(op, physics, scalars, ids,
        records, targets, sat_speeds, controls))
    ns = length(owners)
    st = OpenHallHybrid(op, physics, CIC(), fields, owners, scalars, counts, ids,
        records, fluxes, targets, sat_speeds, controls, configuration, rng, candidate_rng,
        zeros(_OpenHallQ, 2ns), ones(UInt64, ns), ones(UInt64, ns),
        [_OpenHallExchange() for _ = 1:ns], [_OpenHallExchange() for _ = 1:2ns],
        [_OpenHallExchange() for _ = 1:2ns], [_OpenHallExchange() for _ = 1:2ns],
        [_OpenHallExchange() for _ = 1:2ns], zero(_OpenHallQ), _OpenHallCumulative(),
        Ref(false), Ref(zero(T)), Ref(0), Ref(zero(T)), nothing, _open_hall_workspace(op), false)
    _open_hall_binding(st, species; values = false, counts = false)
    _open_hall_storage(st) # Original opaque RNG metadata is checked before use.
    return st
end

function _open_hall_pair_duration(h::T, level::Int) where {T}
    delta = ldexp(h, -level)
    half = delta / T(2)
    isfinite(delta) && isfinite(half) && delta > 0 && half > 0 &&
        (BigInt(1) << level) * _exact_rational(delta) == _exact_rational(h) &&
        2 * _exact_rational(half) == _exact_rational(delta) ||
        throw(ArgumentError("open Hall dyadic field durations must preserve the complete external interval"))
    _require_resolved_dt("OpenHallHybrid paired field trial", half)
    return delta, half
end

_open_hall_magnetic(B, op) = prod(_exact_rational, op.dx) * _mixed_qdot(B, B, op) / 2
_open_hall_power(B, rate, op) = prod(_exact_rational, op.dx) * _mixed_qdot(B, rate, op)

function _open_hall_projection_metrics!(w, report; field)
    equation = report.solve.equation
    constraint = report.solve.constraints
    w.worst_equation = max(w.worst_equation, equation.residual)
    w.worst_constraint = max(w.worst_constraint, constraint.scaled_divergence, constraint.scaled_normal_current)
    work = field ? report.orthogonality_work : report.power_change
    w.worst_work = max(w.worst_work, abs(work))
    return nothing
end
function _open_hall_project_field!(st, out, B)
    c, w = st.controls.projection, st.work
    w.field_projections += 1
    result = try
        project_mixed_field!(out, B, st.op, w.projection; rtol = c.rtol, atol = c.atol,
            work_rtol = c.work_rtol, work_atol = c.work_atol, maxiter = c.maxiter)
    catch err
        err isa MixedHallProjectionError || rethrow()
        throw(_OpenHallRetry(:projection, err))
    end
    _open_hall_projection_metrics!(w, result; field = true)
    return nothing
end
function _open_hall_project_rate!(st, out, rate, B)
    c, w = st.controls.projection, st.work
    w.rate_projections += 1
    result = try
        project_mixed_rate!(out, rate, B, st.op, w.projection; rtol = c.rtol, atol = c.atol,
            state_rtol = c.state_rtol, state_atol = c.state_atol,
            work_rtol = c.work_rtol, work_atol = c.work_atol, maxiter = c.maxiter)
    catch err
        err isa MixedHallProjectionError || rethrow()
        throw(_OpenHallRetry(:projection, err))
    end
    _open_hall_projection_metrics!(w, result; field = false)
    return nothing
end

function _open_hall_rate!(st::OpenHallHybrid{D,T}, out, B, n, ui) where {D,T}
    w, op = st.work, st.op
    w.rhs_evaluations += 1
    result = mixed_hall_rhs!(w.combined, B, n, ui, st.model, op, w.hall)
    for a in w.sat
        fill!(a, zero(T))
    end
    Q = _exact_rational
    for side = 1:2, I in CartesianIndices(Base.tail(op.n)), c = 2:3
        ix = side == 1 ? 1 : op.n[1]
        index = CartesianIndex((ix, Tuple(I)...))
        value = -Q(st.sat_speeds[side]) * (Q(B[c][index]) - Q(st.targets[side][c - 1])) / Q(op.sbp.H[ix])
        w.sat[c][index] = _from_exact(T, value)
    end
    foreach(a -> _mixed_finite(a, "open Hall materialized SAT"), w.sat)
    for c = 1:3, I in CartesianIndices(w.raw[c])
        w.raw[c][I] = _from_exact(T, Q(w.combined[c][I]) + Q(w.sat[c][I]))
    end
    foreach(a -> _mixed_finite(a, "open Hall materialized combined rate"), w.raw)
    raw = result.budget
    sat = _open_hall_power(B, w.sat, op)
    assembled = _open_hall_power(B, w.raw, op)
    _open_hall_project_rate!(st, out, w.raw, B)
    projected = _open_hall_power(B, out, op)
    return (raw = raw.magnetic_power, current = raw.current_work,
        poynting = raw.outward_poynting, identity = raw.identity_residual,
        sat = sat, assembly = assembled - raw.magnetic_power - sat,
        projection = projected - assembled, total = projected)
end

function _open_hall_affine!(out, initial, h, divisor, rates, weights)
    affine = _RK4Affine(h, divisor)
    for c = 1:3, I in CartesianIndices(out[c])
        values = map(k -> k[c][I], rates)
        value = _rk4_affine_exact(affine, initial[c][I], values, weights)
        isfinite(value) || throw(_OpenHallRetry(:affine_range, (component = c, index = I)))
        out[c][I] = value
    end
    return nothing
end

function _open_hall_rk4!(st, out, input, h, n, ui)
    w, op = st.work, st.op
    k = w.rates
    p1 = _open_hall_rate!(st, k[1], input, n, ui)
    stage_energy = zero(_OpenHallQ)
    _open_hall_affine!(w.raw, input, h, 2, (k[1],), (1,))
    before = _open_hall_magnetic(w.raw, op)
    _open_hall_project_field!(st, w.stage, w.raw)
    stage_energy += _open_hall_magnetic(w.stage, op) - before
    p2 = _open_hall_rate!(st, k[2], w.stage, n, ui)
    _open_hall_affine!(w.raw, input, h, 2, (k[2],), (1,))
    before = _open_hall_magnetic(w.raw, op)
    _open_hall_project_field!(st, w.stage, w.raw)
    stage_energy += _open_hall_magnetic(w.stage, op) - before
    p3 = _open_hall_rate!(st, k[3], w.stage, n, ui)
    _open_hall_affine!(w.raw, input, h, 1, (k[3],), (1,))
    before = _open_hall_magnetic(w.raw, op)
    _open_hall_project_field!(st, w.stage, w.raw)
    stage_energy += _open_hall_magnetic(w.stage, op) - before
    p4 = _open_hall_rate!(st, k[4], w.stage, n, ui)
    _open_hall_affine!(w.raw, input, h, 6, k, (1, 2, 2, 1))
    raw_final = _open_hall_magnetic(w.raw, op)
    initial = _open_hall_magnetic(input, op)
    _open_hall_project_field!(st, out, w.raw)
    final = _open_hall_magnetic(out, op)
    coefficients = (_exact_rational(h) / 6, _exact_rational(h) / 3,
        _exact_rational(h) / 3, _exact_rational(h) / 6)
    powers = (p1, p2, p3, p4)
    integral(name) = sum(coefficients[i] * getproperty(powers[i], name) for i = 1:4)
    ledger = _OpenHallFieldLedger(final - initial, integral(:raw), integral(:current),
        integral(:poynting), integral(:identity), integral(:sat), integral(:assembly),
        integral(:projection), raw_final - initial - integral(:total), final - raw_final, stage_energy)
    _open_hall_field_ledger_valid(ledger)
    return ledger
end

function _open_hall_field_flow!(st, h)
    w = st.work
    w.rhs_evaluations = w.field_projections = w.rate_projections = 0
    w.worst_constraint = w.worst_equation = w.worst_work = zero(_OpenHallQ)
    detail = nothing
    for level = 0:st.controls.max_levels
        delta, half = try
            _open_hall_pair_duration(h, level)
        catch err
            err isa ArgumentError || rethrow()
            throw(OpenHallControllerError(level, :duration, err))
        end
        for c = 1:3
            copyto!(w.base[c], st.fields.B[c])
        end
        ledger = _OpenHallFieldLedger()
        error = tolerance = zero(_OpenHallQ)
        accepted = true
        for segment = 1:(Int(1) << level)
            trial = try
                _open_hall_rk4!(st, w.coarse, w.base, delta, st.fields.n, st.fields.ui)
                first = _open_hall_rk4!(st, w.half, w.base, half, st.fields.n, st.fields.ui)
                second = _open_hall_rk4!(st, w.fine, w.half, half, st.fields.n, st.fields.ui)
                _open_hall_field_add(first, second)
            catch err
                err isa _OpenHallRetry || rethrow()
                detail = err
                accepted = false
                break
            end
            difference = maximum(abs(_exact_rational(w.coarse[c][I]) - _exact_rational(w.fine[c][I]))
                for c = 1:3 for I in CartesianIndices(w.base[c]))
            scale = max(_exact_rational(_mixed_vector_maxabs(w.base)), _exact_rational(_mixed_vector_maxabs(w.fine)))
            threshold = _exact_rational(st.controls.field_atol) + _exact_rational(st.controls.field_rtol) * scale
            if difference >= error
                error, tolerance = difference, threshold
            end
            if difference > threshold
                detail = (segment = segment, error = difference, tolerance = threshold)
                accepted = false
                break
            end
            ledger = _open_hall_field_add(ledger, trial)
            for c = 1:3
                copyto!(w.base[c], w.fine[c])
            end
        end
        accepted && return (B = w.base, ledger = ledger, level = level,
            error = error, tolerance = tolerance)
    end
    throw(OpenHallControllerError(st.controls.max_levels, :trial_limit, detail))
end

function _open_hall_candidate_storage(st, prepared)
    arrays = _open_hall_particle_arrays(prepared)
    owned = _open_hall_owned_arrays(st)
    for i in eachindex(arrays)
        for j = 1:(i - 1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(ArgumentError("open Hall candidate particle leaves must be independent"))
        end
        for live in owned
            _particle_arrays_mightalias(arrays[i], live) && throw(ArgumentError("open Hall candidate particles must be detached from live/private storage"))
        end
    end
    return nothing
end

"Initialize one synchronous open Hall boundary without drawing or priming particles."
function init!(st::OpenHallHybrid{D,T}, species = st.owners) where {D,T}
    _acquire_open_hall!(st)
    try
        _open_hall_binding(st, species; values = true, counts = false)
        !st.initialized[] || throw(ArgumentError("open Hall init! is one-time; construct a fresh model to reset"))
        _open_hall_phase_values(st, false, st.time[], st.step[], st.last_dt[], st.last_report)
        st.step[] == 0 || throw(ArgumentError("open Hall initialization requires its unadvanced clock"))
        _validate_collisional_rng(st.rng)
        for a in st.fields.B
            _mixed_finite(a, "open Hall seed B")
        end
        _open_hall_direct_state(st, st.fields.B)
        endpoint = st.work.endpoint
        for c = 1:3
            copyto!(endpoint.B[c], st.fields.B[c])
        end
        _open_hall_moments!(st, endpoint, species)
        _open_hall_complete_fields!(st, endpoint)
        _open_hall_fields(st, endpoint)
        prepared, counts = _open_hall_replacements(st, species)
        _open_hall_reserve!(st, prepared)
        _open_hall_candidate_storage(st, prepared)
        initial = [_open_hall_particle_totals((ps,)) for ps in prepared]
        indices = UInt64[maximum(id -> id & _PARTICLE_ID_MAX_INDEX, ps.id; init = UInt64(0)) + UInt64(1)
            for ps in prepared]
        ns = length(prepared)
        source = (alpha = zeros(_OpenHallQ, 2ns), next_index = copy(indices), initial_next_index = indices,
            initial_totals = initial, births = [_OpenHallExchange() for _ = 1:2ns],
            exits = [_OpenHallExchange() for _ = 1:2ns], last_births = [_OpenHallExchange() for _ = 1:2ns],
            last_exits = [_OpenHallExchange() for _ = 1:2ns],
            source_elapsed = zero(_OpenHallQ), cumulative = _OpenHallCumulative())
        _open_hall_source_values(st, source, prepared, 0)
        _open_hall_phase_values(st, true, zero(T), 0, zero(T), nothing)
        return _open_hall_commit!(st, prepared, counts, endpoint, source, st.rng,
            true, zero(T), 0, zero(T), nothing)
    finally
        st.busy = false
    end
end

"""
    step!(model::OpenHallHybrid, species, dt)

Advance the complete synchronous first-order Lie map: frozen-moment projected
Hall field -> original-factor Boris velocities -> full-time ballistic exits and
bounded reservoir births -> strict endpoint moments/Ohm. Only a completely
prepared finite boundary publishes; source/RNG/counters roll back by staging.
True zero validates read-only state and performs no transforms or source draws.
"""
function step!(st::OpenHallHybrid{D,T}, species::AbstractVector, dt::Real) where {D,T}
    _acquire_open_hall!(st)
    try
        _open_hall_boundary(st, species; differentiated = false)
        isfinite(dt) && dt >= zero(dt) || throw(ArgumentError("open Hall dt must be originally finite and nonnegative"))
        h = _validated_nonnegative_dt(T, dt; name = "open Hall step!")
        _particle_stored_zero(dt) || !_particle_stored_zero(h) ||
            throw(ArgumentError("open Hall nonzero original dt disappears at grid precision"))
        _particle_stored_zero(h) && return st
        next_time, next_step = _validated_step_clock(st.time[], st.step[], h)
        next_time > st.time[] || throw(ArgumentError("open Hall clock must strictly advance at grid precision"))
        _open_hall_pair_duration(h, 0)
        plan = _open_hall_birth_plan(st, h)
        _open_hall_reserve!(st, species)
        _open_hall_direct_state(st, st.fields.B)
        field = _open_hall_field_flow!(st, h)
        force, endpoint = st.work.force, st.work.endpoint
        _copy_hybrid_fields!(force, st.fields)
        for c = 1:3
            copyto!(force.B[c], field.B[c])
        end
        _open_hall_complete_fields!(st, force)
        _open_hall_fields(st, force)
        transport = _open_hall_transport!(st, h, plan)
        _open_hall_candidate_storage(st, transport.particles)
        for c = 1:3
            copyto!(endpoint.B[c], field.B[c])
        end
        _open_hall_moments!(st, endpoint, transport.particles)
        _open_hall_complete_fields!(st, endpoint)
        _open_hall_fields(st, endpoint)
        _open_hall_direct_state(st, endpoint.B)
        energy_before = _open_hall_energy_interval(st, st.fields, species)
        energy_after = _open_hall_energy_interval(st, endpoint, transport.particles)
        energy_change = _open_hall_energy_change_interval(st, st.fields, endpoint, species, transport.particles)
        enthalpy = _open_hall_enthalpy_interval(st, st.fields, endpoint, h)
        residual = (energy_change[1] - transport.event_energy + field.ledger.outward_poynting + enthalpy[1],
            energy_change[2] - transport.event_energy + field.ledger.outward_poynting + enthalpy[2])
        report = OpenHallStepReport(next_step, next_time, h, field.level, Int(1) << field.level,
            st.work.rhs_evaluations, st.work.field_projections, st.work.rate_projections,
            transport.draws, transport.births, transport.exits, field.error, field.tolerance,
            st.work.worst_constraint, st.work.worst_equation, st.work.worst_work,
            field.ledger, transport.kick_energy, transport.kick_work, transport.kick_momentum,
            transport.event_energy, energy_before, energy_after, energy_change, enthalpy, residual,
            _open_hall_grid_charge(st, st.fields), _open_hall_grid_charge(st, endpoint),
            transport.old_totals.charge, transport.new_totals.charge,
            _exact_rational(next_time) - _exact_rational(st.time[]) - _exact_rational(h))
        old = st.cumulative
        cumulative = _OpenHallCumulative(old.duration + _exact_rational(h),
            _open_hall_field_add(old.field, field.ledger), old.kick_energy + transport.kick_energy,
            old.kick_work + transport.kick_work,
            ntuple(c -> old.kick_momentum[c] + transport.kick_momentum[c], 3),
            old.event_energy + transport.event_energy,
            ntuple(i -> old.electron_enthalpy[i] + enthalpy[i], 2),
            ntuple(i -> old.energy_residual[i] + residual[i], 2))
        source = (alpha = plan.alpha, next_index = plan.next_index,
            initial_next_index = copy(st.initial_next_index), initial_totals = copy(st.initial_totals),
            births = [_open_hall_exchange_add(st.births[i], transport.born[i]) for i in eachindex(st.births)],
            exits = [_open_hall_exchange_add(st.exits[i], transport.exited[i]) for i in eachindex(st.exits)],
            last_births = transport.born, last_exits = transport.exited,
            source_elapsed = st.source_elapsed + _exact_rational(h), cumulative = cumulative)
        _open_hall_source_values(st, source, transport.particles, next_step; report)
        _open_hall_phase_values(st, true, next_time, next_step, h, report)
        return _open_hall_commit!(st, transport.particles, transport.counts, endpoint, source,
            st.rng_candidate, true, next_time, next_step, h, report)
    finally
        st.busy = false
    end
end
step!(st::OpenHallHybrid, dt::Real) = step!(st, st.owners, dt)
