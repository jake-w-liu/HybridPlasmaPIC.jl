"Actual signed KDK work on the receipt's own Cartesian composite mesh."
struct AdaptiveCompositeESStepReport{D,T<:Union{Float32,Float64}}
    dt::T
    time_before::T
    time_after::T
    step_before::Int
    step_after::Int
    epoch_start_time::T
    epoch_start_step::Int
    refinements::Tuple
    kinetic_before::_CompositeQ
    kinetic_after::_CompositeQ
    field_before::_CompositeQ
    field_after::_CompositeQ
    source_work::_CompositeQ
    kick_work::_CompositeQ
    field_defect::_CompositeQ
    kick_defect::_CompositeQ
    trajectory_defect::_CompositeQ
    total_defect::_CompositeQ
    max_residual::_CompositeQ
    crossings::Vector{_AdaptiveCompositeCrossing{D}}
end
"Signed field-energy exchange of a completed-boundary mesh reconstruction."
struct AdaptiveCompositeESRegridReport{D,T<:Union{Float32,Float64}}
    time::T
    step::Int
    regrid::Int
    before_refinements::Tuple
    after_refinements::Tuple
    before_leaves::Int
    after_leaves::Int
    before_dofs::Int
    after_dofs::Int
    field_before::_CompositeQ
    field_after::_CompositeQ
    kinetic_before::_CompositeQ
    kinetic_after::_CompositeQ
    mesh_work::_CompositeQ
    before_inventory::NamedTuple
    after_inventory::NamedTuple
    max_residual::_CompositeQ
    selection::NamedTuple
end
struct _AdaptiveCompositePhysical{D,T<:Union{Float32,Float64}}
    mesh::CartesianAMRHierarchy{D,T}
    operator::_AdaptiveQ1Data
    species::Vector{_ElectrostaticHostParticle{D,T}}
    phi::Vector{T}
    bins::_CompositeParticleBins
    solve_report::AdaptiveCompositeESSolveReport{T}
    step_report::Union{Nothing,AdaptiveCompositeESStepReport{D,T}}
    transport::Union{Nothing,_PeriodicCompositeTransport{D,T}}
    background_density::Union{Nothing,_CompositeQ}
end
mutable struct AdaptiveCompositeElectrostaticPIC{D,T<:Union{Float32,Float64},S,B,F}
    mesh::CartesianAMRHierarchy{D,T}
    phi::Vector{T}
    time::T
    step::Int
    last_dt::T
    initialized::Bool
    epoch_start_time::T
    epoch_start_step::Int
    regrid_count::Int
    initial_energy::_CompositeQ
    field_defect::_CompositeQ
    kick_defect::_CompositeQ
    trajectory_defect::_CompositeQ
    physical_defect::_CompositeQ
    mesh_work::_CompositeQ
    last_regrid_report::Union{Nothing,AdaptiveCompositeESRegridReport{D,T}}
    clouds::S
    background_charge::B
    rtol::T
    atol::T
    maxiter::Int
    periodic_neutrality::Symbol
    _factors::F
    _configuration::NamedTuple
    _state::Union{Nothing,_AdaptiveCompositePhysical{D,T}}
    _history::Union{Nothing,NamedTuple}
    _owners::Vector{_ElectrostaticHostParticle{D,T}}
    _components::Vector{Tuple}
    _counts::Vector{Int}
    _busy::Threads.Atomic{Bool}
end
function _ace_configuration(st)
    st.periodic_neutrality in (:require, :neutralize) &&
        (st.mesh.periodic || st.periodic_neutrality === :require) ||
        throw(ArgumentError("adaptive composite neutrality policy does not match its boundary"))
    base = (
        D = length(st.mesh.nbase),
        T = eltype(st.phi),
        nbase = st.mesh.nbase,
        lengths = map(_sce_ieee_config, st.mesh.L),
        ratio = 2,
        ghosts = 0,
        clouds = map(c -> _sce_ieee_config(c.half_width), st.clouds),
        background = _sce_ieee_config(st.background_charge),
        factors = map(pair -> map(_sce_ieee_config, pair), st._factors),
        rtol = _sce_ieee_config(st.rtol),
        atol = _sce_ieee_config(st.atol),
        maxiter = st.maxiter,
        boundary = :grounded_closest_point_bound_tails_v1,
        field = :constrained_q1_original_gradient_v1,
        solve = :cold_scaled_jacobi_cg_original_residual_v1,
        force = :original_tensor_cloud_boris_v1,
        phase = :synchronous_global_fixed_duration_v1,
        adaptation = :multilevel_curvature_hysteresis_v1,
        arithmetic = :nearest_gradual,
        events = :strict_straight_path_interior_v1,
    )
    return st.mesh.periodic ?
           merge(
        base,
        (
            boundary = :periodic_torus_v1,
            field = :periodic_constrained_q1_original_gradient_v1,
            events = :represented_lifted_multitile_v1,
            periodic_neutrality = st.periodic_neutrality,
            gauge = :origin_master_zero_v1,
            harmonic_field = :zero_periodic_gradient_v1,
            periodization = :exact_integer_cubic_bands_v1,
        ),
    ) : base
end
function _ace_populations(species, ::Val{D}, ::Type{T}) where {D,T}
    species isa AbstractVector ||
        throw(ArgumentError("adaptive composite populations require an ordered vector"))
    Base.require_one_based_indexing(species)
    for ps in species
        ps isa _ElectrostaticHostParticle{D,T} || throw(
            ArgumentError(
                "adaptive composite particles require D-dimensional ordinary Vector storage at mesh precision",
            ),
        )
        _check_particle_storage_axes(ps)
        _require_particle_ordinary_stores(ps, "adaptive composite ES")
    end
    return nothing
end
function _ace_values(mesh::CartesianAMRHierarchy{D,T}, species, factors) where {D,T}
    _ace_populations(species, Val(D), T)
    length(species) == length(factors) ||
        throw(ArgumentError("adaptive composite ordered population count changed"))
    for (ps, original) in zip(species, factors)
        isequal(_validated_charge_mass(ps), original) ||
            throw(ArgumentError("adaptive composite original q/m differ from construction"))
        all(v -> all(isfinite, v), ps.v) ||
            throw(ArgumentError("adaptive composite velocity is nonfinite"))
        all(w -> isfinite(w) && _stored_nonnegative(w), ps.weight) ||
            throw(ArgumentError("adaptive composite weights must be finite nonnegative"))
        for p in eachindex(ps.weight)
            all(d -> isfinite(ps.x[d][p]), 1:D) ||
                throw(ArgumentError("adaptive composite position is nonfinite"))
            _camr_locate(mesh, ntuple(d -> _exact_rational(ps.x[d][p]), D))
        end
    end
    return nothing
end
function _ace_copy(species, ::Val{D}, ::Type{T}) where {D,T}
    out = _ElectrostaticHostParticle{D,T}[]
    _pic_prepare_species_stages!(out, species)
    return out
end
function _ace_species_equal(a, b)
    length(a) == length(b) || return false
    return all(zip(a, b)) do pair
        left, right = pair
        isequal((left.q, left.m), (right.q, right.m)) && all(
            p -> isequal(p[1], p[2]),
            zip(_sce_particle_arrays(left), _sce_particle_arrays(right)),
        )
    end
end
function _ace_independent(arrays)
    for i in eachindex(arrays), j = 1:(i-1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) && throw(
            ArgumentError("adaptive composite physical/private/history storage roles overlap"),
        )
    end
    return nothing
end
function _ace_particle_arrays(species)
    arrays = Any[]
    for ps in species
        append!(arrays, _sce_particle_arrays(ps))
    end
    return arrays
end
function _ace_physical_arrays(state; geometry = true)
    arrays = _ace_particle_arrays(state.species)
    push!(arrays, state.phi)
    append!(arrays, (state.bins.offsets, state.bins.species, state.bins.particle))
    state.step_report === nothing || push!(arrays, state.step_report.crossings)
    state.transport === nothing || append!(arrays, _pacc_transport_arrays(state.transport))
    if geometry
        append!(arrays, _camr_arrays(state.mesh))
        push!(arrays, state.operator.elements, state.operator.diagonal)
        append!(arrays, state.operator.elements)
    end
    return arrays
end
function _ace_live_arrays(st, species)
    arrays = _ace_particle_arrays(species)
    push!(arrays, st.phi, st._counts)
    append!(arrays, _camr_arrays(st.mesh))
    if st._state !== nothing
        append!(arrays, _ace_physical_arrays(st._state; geometry = false))
        push!(arrays, st._state.operator.elements, st._state.operator.diagonal)
        append!(arrays, st._state.operator.elements)
    end
    st._history === nothing || append!(arrays, _ace_history_arrays(st._history, st))
    return arrays
end
function _ace_storage(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    species;
    binding = true,
    counts = true,
) where {D,T}
    _sce_mode()
    _camr_validate(st.mesh)
    _ace_populations(species, Val(D), T)
    (st._state === nothing || st.mesh === st._state.mesh) ||
        throw(ArgumentError("adaptive composite physical mesh binding changed"))
    st._state === nothing || _acq_validate(st._state.operator, st.mesh)
    isequal(st._configuration, _ace_configuration(st)) ||
        throw(ArgumentError("adaptive composite immutable configuration changed"))
    st.phi isa Vector{T} && length(st.phi) == length(st.mesh.free) ||
        throw(DimensionMismatch("adaptive composite public potential storage changed"))
    length(species) == length(st._factors) ||
        throw(ArgumentError("adaptive composite population count differs"))
    if binding
        length(st._owners) == length(st._components) == length(st._counts) == length(species) ||
            throw(ArgumentError("adaptive composite binding metadata changed"))
        for (i, ps) in enumerate(species)
            ps === st._owners[i] &&
                length(st._components[i]) == D + 6 &&
                all(p -> p[1] === p[2], zip(_sce_particle_arrays(ps), st._components[i])) || throw(
                ArgumentError(
                    "adaptive composite population/component identity changed; call init!",
                ),
            )
            counts &&
                length(ps.weight) != st._counts[i] &&
                throw(ArgumentError("adaptive composite particle count changed; call init!"))
        end
    end
    _ace_independent(_ace_live_arrays(st, species))
    return nothing
end
function _ace_acquire(st)
    Threads.atomic_cas!(st._busy, false, true) == false ||
        throw(ArgumentError("adaptive composite ES already has an active operation"))
    return nothing
end
_ace_release(st) = (Threads.atomic_xchg!(st._busy, false); nothing)
_ace_source(st, mesh, species) = _acc_source(
    mesh,
    species,
    st.clouds,
    st.background_charge;
    periodic_neutrality = st.periodic_neutrality,
)
# The bound background is a periodic-source quantity. Dispatch keeps the shared
# source union from reading a field that only one inventory owns.
_ace_bound_background(source::_PeriodicCompositeSource) = source.background_density
_ace_bound_background(::_AdaptiveCompositeSource) =
    throw(ArgumentError("grounded composite source has no periodic bound background"))
function _ace_physical(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    mesh,
    species;
    operator = _acq_data(mesh),
    step_report = nothing,
    transport = nothing,
) where {D,T}
    _ace_values(mesh, species, st._factors)
    source = _ace_source(st, mesh, species)
    phi, report = _acq_solve(operator, source.total, st.rtol, st.atol, st.maxiter)
    _acq_means(mesh, phi)
    return _AdaptiveCompositePhysical{D,T}(
        mesh,
        operator,
        species,
        phi,
        _acc_bins(species, mesh),
        report,
        step_report,
        transport,
        source isa _PeriodicCompositeSource ? source.background_density : nothing,
    )
end
function _ace_physical_check(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    state::_AdaptiveCompositePhysical{D,T},
) where {D,T}
    _camr_validate(state.mesh)
    _acq_validate(state.operator, state.mesh)
    _ace_values(state.mesh, state.species, st._factors)
    length(state.phi) == length(state.mesh.free) && all(isfinite, state.phi) ||
        throw(ArgumentError("adaptive composite physical potential is invalid"))
    source = _ace_source(st, state.mesh, state.species)
    _acq_report(
        state.operator,
        state.phi,
        source.total,
        state.solve_report,
        st.rtol,
        st.atol,
        st.maxiter,
    )
    if state.mesh.periodic
        _sce_q_finite(state.background_density) &&
            state.background_density == _ace_bound_background(source) ||
            throw(ArgumentError("periodic composite bound background density differs"))
    else
        state.background_density === nothing ||
            throw(ArgumentError("grounded composite has periodic background state"))
    end
    _acq_means(state.mesh, state.phi)
    _acc_validate_bins(state.bins, state.species, state.mesh)
    _ace_independent(_ace_physical_arrays(state))
    return source
end
_ace_energy(state) =
    _kinetic_energy_rational(state.species) + _acq_energy(state.operator, state.phi)
function _ace_inventory(state, source::_AdaptiveCompositeSource)
    reaction, free = _acq_residual(state.operator, state.phi, source.total)
    wall = sum(reaction[state.mesh.boundary]; init = zero(_CompositeQ))
    tail = sum(source.tail; init = zero(_CompositeQ))
    background = sum(source.background; init = zero(_CompositeQ))
    inside = sum(source.inside; init = zero(_CompositeQ))
    full = sum(source.strata; init = zero(_CompositeQ))
    residual = sum(free; init = zero(_CompositeQ))
    full + background + wall + residual == 0 && inside + background + wall + tail + residual == 0 ||
        throw(ArgumentError("adaptive composite weak conductor charge closure failed"))
    return (
        particle = full,
        inside = inside,
        tail = tail,
        strata = Tuple(source.strata),
        background = background,
        free_conductor = wall,
        net_wall = wall + tail,
        free_residual_sum = residual,
    )
end
function _ace_inventory(state, source::_PeriodicCompositeSource)
    residual, _ = _acq_residual(state.operator, state.phi, source.total)
    particle = sum(source.particle; init = zero(_CompositeQ))
    background = sum(source.background; init = zero(_CompositeQ))
    total = sum(source.total; init = zero(_CompositeQ))
    residual_sum = sum(residual; init = zero(_CompositeQ))
    particle + background == total == 0 &&
        background == source.declared_background + source.neutralizing_charge &&
        iszero(residual_sum) ||
        throw(ArgumentError("periodic composite physical charge/null-mode closure failed"))
    return (
        particle = particle,
        declared_background = source.declared_background,
        neutralizing_charge = source.neutralizing_charge,
        background = background,
        source_total = total,
        gauge_residual = residual[state.mesh.gauge],
        full_residual_sum = residual_sum,
    )
end
function _ace_step_report(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    state,
    time,
    step,
    dt,
    epoch_time,
    epoch_step,
) where {D,T}
    report = state.step_report
    report isa AdaptiveCompositeESStepReport{D,T} ||
        throw(ArgumentError("adaptive composite last step receipt is absent"))
    isequal(report.dt, dt) &&
        report.time_after === time &&
        report.step_after == step &&
        report.step_before == step - 1 &&
        report.epoch_start_time === epoch_time &&
        report.epoch_start_step == epoch_step &&
        isequal(report.refinements, state.mesh.refinements) ||
        throw(ArgumentError("adaptive composite last step chronology differs"))
    next, count = _validated_step_clock(report.time_before, report.step_before, dt)
    isequal(next, time) &&
        count == step &&
        report.time_before < time &&
        epoch_time <= report.time_before ||
        throw(ArgumentError("adaptive composite last step interval is invalid"))
    values = (
        report.kinetic_before,
        report.kinetic_after,
        report.field_before,
        report.field_after,
        report.source_work,
        report.kick_work,
        report.field_defect,
        report.kick_defect,
        report.trajectory_defect,
        report.total_defect,
        report.max_residual,
    )
    all(_sce_q_finite, values) &&
        all(>=(zero(_CompositeQ)), values[1:4]) &&
        report.max_residual >= 0 ||
        throw(ArgumentError("adaptive composite step exact scalars are invalid"))
    du, dk = report.field_after - report.field_before, report.kinetic_after - report.kinetic_before
    du + report.source_work == report.field_defect &&
        dk - report.kick_work == report.kick_defect &&
        report.kick_work - report.source_work == report.trajectory_defect &&
        du + dk ==
        report.total_defect ==
        report.field_defect + report.kick_defect + report.trajectory_defect ||
        throw(ArgumentError("adaptive composite step work identities differ"))
    report.field_after == _acq_energy(state.operator, state.phi) &&
        report.kinetic_after == _kinetic_energy_rational(state.species) &&
        report.max_residual == state.solve_report.max_residual || throw(
        ArgumentError("adaptive composite step receipt does not match its ending physical state"),
    )
    for event in report.crossings
        1 <= event.population <= length(state.species) &&
            1 <= event.particle <= length(state.species[event.population].weight) &&
            0 <= event.parameter <= 1 &&
            _sce_q_finite(event.parameter) &&
            1 <= event.from <= length(state.mesh.addresses) &&
            1 <= event.to <= length(state.mesh.addresses) &&
            event.from != event.to &&
            any(event.axes) ||
            throw(ArgumentError("adaptive composite crossing receipt is invalid"))
    end
    if state.mesh.periodic
        _pacc_transport_check(state.mesh, state.species, state.transport, report)
    else
        state.transport === nothing ||
            throw(ArgumentError("grounded composite state has periodic transport"))
    end
    return nothing
end
_ace_phase(st) = (
    time = st.time,
    step = st.step,
    last_dt = st.last_dt,
    epoch_start_time = st.epoch_start_time,
    epoch_start_step = st.epoch_start_step,
    regrid_count = st.regrid_count,
    initial_energy = st.initial_energy,
    field_defect = st.field_defect,
    kick_defect = st.kick_defect,
    trajectory_defect = st.trajectory_defect,
    physical_defect = st.physical_defect,
    mesh_work = st.mesh_work,
    last_regrid_report = st.last_regrid_report,
)
function _ace_phase_check(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    state,
    history,
    phase,
) where {D,T}
    t, n, h = phase.time, phase.step, phase.last_dt
    et, en = phase.epoch_start_time, phase.epoch_start_step
    t isa T &&
        h isa T &&
        et isa T &&
        n isa Int &&
        en isa Int &&
        phase.regrid_count isa Int &&
        isfinite(t) &&
        isfinite(h) &&
        isfinite(et) &&
        _stored_nonnegative(t) &&
        _stored_nonnegative(h) &&
        _stored_nonnegative(et) &&
        0 <= en <= n &&
        0 <= et <= t &&
        phase.regrid_count >= 0 ||
        throw(ArgumentError("adaptive composite global/epoch phase types or bounds differ"))
    if iszero(n)
        _particle_stored_zero(t) &&
            _particle_stored_zero(h) &&
            _particle_stored_zero(et) &&
            en == 0 || throw(ArgumentError("unadvanced adaptive composite clock must be zero"))
    else
        _stored_positive(t) && _stored_positive(h) && h <= t ||
            throw(ArgumentError("advanced adaptive composite clock must retain positive duration"))
    end
    (iszero(en) ? _particle_stored_zero(et) : _stored_positive(et) && h <= et) ||
        throw(ArgumentError("adaptive composite epoch origin is inconsistent"))
    fields = (
        phase.initial_energy,
        phase.field_defect,
        phase.kick_defect,
        phase.trajectory_defect,
        phase.physical_defect,
        phase.mesh_work,
    )
    all(_sce_q_finite, fields) &&
        phase.initial_energy >= 0 &&
        phase.physical_defect == phase.field_defect + phase.kick_defect + phase.trajectory_defect &&
        _ace_energy(state) - phase.initial_energy == phase.physical_defect + phase.mesh_work ||
        throw(ArgumentError("adaptive composite cumulative energy ledger differs"))
    if iszero(n)
        all(iszero, fields[2:5]) &&
            state.step_report === nothing &&
            state.transport === nothing &&
            history === nothing ||
            throw(ArgumentError("unadvanced adaptive composite physical history is not empty"))
    elseif n == en
        isequal(t, et) && state.step_report === nothing && state.transport === nothing ||
            throw(ArgumentError("empty adaptive composite epoch has a step/time change"))
        _ace_validate_history(st, history, state, phase)
    else
        t > et && history === nothing ||
            throw(ArgumentError("advanced adaptive composite epoch has stale prior history"))
        _ace_step_report(st, state, t, n, h, et, en)
    end
    if iszero(phase.regrid_count)
        phase.last_regrid_report === nothing &&
            iszero(phase.mesh_work) &&
            en == 0 &&
            _particle_stored_zero(et) ||
            throw(ArgumentError("unregridded adaptive composite has mesh accounting/origin"))
    else
        _ace_regrid_report_check(st, state, phase)
    end
    return nothing
end
function _ace_ready(st::AdaptiveCompositeElectrostaticPIC{D,T}, species::AbstractVector) where {D,T}
    _ace_storage(st, species)
    st.initialized || throw(ArgumentError("call init! before adaptive composite evolution or save"))
    state = st._state
    state === nothing &&
        throw(ArgumentError("call init! before adaptive composite evolution or save"))
    st.mesh === state.mesh &&
        isequal(st.phi, state.phi) &&
        _ace_species_equal(species, state.species) ||
        throw(ArgumentError("adaptive composite physical/source values changed; call init!"))
    source = _ace_physical_check(st, state)
    _ace_phase_check(st, state, st._history, _ace_phase(st))
    return source
end
function _ace_prepared(st, state, public, public_phi, history)
    arrays = _ace_particle_arrays(public)
    push!(arrays, public_phi)
    append!(arrays, _ace_physical_arrays(state))
    history === nothing || append!(arrays, _ace_history_arrays(history, st))
    _ace_independent(arrays)
    return nothing
end
function _ace_publish!(st, state, phi, owners, components, counts, history, phase)
    st.mesh, st.phi, st._state = state.mesh, phi, state
    st._owners, st._components, st._counts, st._history = owners, components, counts, history
    st.time, st.step, st.last_dt = phase.time, phase.step, phase.last_dt
    st.epoch_start_time, st.epoch_start_step = phase.epoch_start_time, phase.epoch_start_step
    st.regrid_count, st.initial_energy = phase.regrid_count, phase.initial_energy
    st.field_defect, st.kick_defect, st.trajectory_defect =
        phase.field_defect, phase.kick_defect, phase.trajectory_defect
    st.physical_defect, st.mesh_work, st.last_regrid_report =
        phase.physical_defect, phase.mesh_work, phase.last_regrid_report
    st.initialized = true
    return st
end
"Coupled 2D3V/3D3V PIC with fixed physical tensor clouds and grounded or periodic adaptive Q1 fields."
function AdaptiveCompositeElectrostaticPIC(
    mesh::CartesianAMRHierarchy{D,T},
    species::AbstractVector;
    clouds::Tuple,
    background_charge = zero(T),
    rtol = sqrt(eps(T)),
    atol = zero(T),
    maxiter::Integer = 0,
    periodic_neutrality::Symbol = :require,
) where {D,T}
    _sce_mode()
    _camr_validate(mesh)
    length(clouds) == D && all(c -> c isa YeeTriangularCloud, clouds) ||
        throw(ArgumentError("adaptive composite requires one original physical cloud per axis"))
    all(c -> isfinite(c.half_width) && _stored_positive(c.half_width), clouds) ||
        throw(ArgumentError("adaptive composite original cloud widths must be finite positive"))
    all(v -> isfinite(v) && v >= 0, (rtol, atol)) || throw(
        ArgumentError("adaptive composite tolerances must be finite nonnegative before conversion"),
    )
    relative, absolute = T(rtol), T(atol)
    isfinite(relative) &&
        isfinite(absolute) &&
        0 <= relative < 1 &&
        absolute >= 0 &&
        (relative > 0 || absolute > 0) ||
        throw(ArgumentError("adaptive composite converted tolerance pair is invalid"))
    0 <= maxiter <= typemax(Int) ||
        throw(ArgumentError("adaptive composite iteration limit must be a nonnegative Int"))
    periodic_neutrality in (:require, :neutralize) &&
        (mesh.periodic || periodic_neutrality === :require) ||
        throw(ArgumentError("adaptive composite neutrality policy does not match its boundary"))
    _ace_populations(species, Val(D), T)
    factors = Tuple(_validated_charge_mass(ps) for ps in species)
    _ace_values(mesh, species, factors)
    _ace_independent(_ace_particle_arrays(species))
    owned = CartesianAMRHierarchy(mesh.nbase, mesh.L, mesh.refinements; periodic = mesh.periodic)
    b = _yl_source_scalar("background_charge", background_charge, T)
    mesh.periodic &&
        periodic_neutrality === :require &&
        !iszero(_pic_collection_charges(owned, b, species).net) &&
        throw(ArgumentError("periodic composite Poisson requires exact physical neutrality"))
    z = zero(_CompositeQ)
    st = AdaptiveCompositeElectrostaticPIC{D,T,typeof(clouds),typeof(b),typeof(factors)}(
        owned,
        zeros(T, length(owned.free)),
        zero(T),
        0,
        zero(T),
        false,
        zero(T),
        0,
        0,
        z,
        z,
        z,
        z,
        z,
        z,
        nothing,
        clouds,
        b,
        relative,
        absolute,
        Int(maxiter),
        periodic_neutrality,
        factors,
        (;),
        nothing,
        nothing,
        _ElectrostaticHostParticle{D,T}[ps for ps in species],
        Tuple[_sce_particle_arrays(ps) for ps in species],
        Int[length(ps.weight) for ps in species],
        Threads.Atomic{Bool}(false),
    )
    st._configuration = _ace_configuration(st)
    _ace_storage(st, species)
    return st
end
function init!(st::AdaptiveCompositeElectrostaticPIC{D,T}, species::AbstractVector) where {D,T}
    _ace_acquire(st)
    try
        _ace_storage(st, species; binding = false, counts = false)
        _ace_values(st.mesh, species, st._factors)
        private = _ace_copy(species, Val(D), T)
        state = _ace_physical(st, st.mesh, private)
        z = zero(_CompositeQ)
        phase = (
            time = zero(T),
            step = 0,
            last_dt = zero(T),
            epoch_start_time = zero(T),
            epoch_start_step = 0,
            regrid_count = 0,
            initial_energy = _ace_energy(state),
            field_defect = z,
            kick_defect = z,
            trajectory_defect = z,
            physical_defect = z,
            mesh_work = z,
            last_regrid_report = nothing,
        )
        phi = copy(state.phi)
        owners, components, counts = _ElectrostaticHostParticle{D,T}[ps for ps in species],
        Tuple[_sce_particle_arrays(ps) for ps in species],
        Int[length(ps.weight) for ps in species]
        _ace_phase_check(st, state, nothing, phase)
        _ace_prepared(st, state, species, phi, nothing)
        return _ace_publish!(st, state, phi, owners, components, counts, nothing, phase)
    finally
        _ace_release(st)
    end
end
function step!(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    species::AbstractVector,
    dt::Real,
) where {D,T}
    _ace_acquire(st)
    try
        tau = _validated_nonnegative_dt(T, dt; name = "adaptive composite ES")
        _particle_stored_zero(tau) &&
            !iszero(dt) &&
            throw(ArgumentError("nonzero adaptive composite duration rounds to zero"))
        oldsource = _ace_ready(st, species)
        _particle_stored_zero(tau) && return st
        iszero(st.step) ||
            isequal(tau, st.last_dt) ||
            throw(ArgumentError("adaptive composite duration is fixed after first success"))
        time, count = _validated_step_clock(st.time, st.step, tau)
        time > st.time ||
            throw(ArgumentError("adaptive composite global time must advance strictly"))
        old = st._state
        private = _ace_copy(old.species, Val(D), T)
        raw = _acq_raw_potential(old.mesh, old.phi)
        half, work = _exact_rational(tau) / 2, zero(_CompositeQ)
        factors = [_boris_exact_factor(_validated_boris_factor(ps)) for ps in old.species]
        transport = old.mesh.periodic ? _pacc_transport(old.species, Val(D), T) : nothing
        for (s, (before, after, factor)) in enumerate(zip(old.species, private, factors)),
            p in eachindex(before.weight)

            x = ntuple(d -> _exact_rational(before.x[d][p]), D)
            field = _acc_gather(old.mesh, raw, x, st.clouds)
            v0 = ntuple(c -> before.v[c][p], 3)
            vhalf =
                _yee_particle_half_kick(v0, field, ntuple(_ -> zero(_CompositeQ), 3), factor, half)
            all(isfinite, vhalf) ||
                throw(ArgumentError("adaptive composite first kick is nonfinite"))
            if old.mesh.periodic
                for d = 1:D
                    ideal =
                        _particle_coordinate(x[d], _exact_rational(vhalf[d]), _exact_rational(tau))
                    represented = _periodic_represented_lift(ideal, old.mesh.L[d])
                    after.x[d][p] = represented.wrapped
                    transport.half_velocity[s][d][p] = vhalf[d]
                    transport.winding[s][d][p] = represented.winding
                    transport.defect[s][d][p] = represented.defect
                end
            else
                exact = ntuple(d -> x[d] + _exact_rational(tau) * _exact_rational(vhalf[d]), D)
                _camr_locate(old.mesh, exact)
                for d = 1:D
                    after.x[d][p] = _particle_coordinate(before.x[d][p], vhalf[d], tau)
                end
            end
            charge = _yee_macrocharge(before, p, Val(true))
            for c = 1:3
                after.v[c][p] = vhalf[c]
                work +=
                    charge *
                    half *
                    field[c] *
                    (_exact_rational(v0[c]) + _exact_rational(vhalf[c])) / 2
            end
        end
        _ace_values(old.mesh, private, st._factors)
        endpoint = _ace_physical(st, old.mesh, private; operator = old.operator)
        rawnew = _acq_raw_potential(old.mesh, endpoint.phi)
        for (after, factor) in zip(private, factors), p in eachindex(after.weight)
            field = _acc_gather(
                old.mesh,
                rawnew,
                ntuple(d -> _exact_rational(after.x[d][p]), D),
                st.clouds,
            )
            halfv = ntuple(c -> after.v[c][p], 3)
            v1 = _yee_particle_half_kick(
                halfv,
                field,
                ntuple(_ -> zero(_CompositeQ), 3),
                factor,
                half,
            )
            all(isfinite, v1) || throw(ArgumentError("adaptive composite final kick is nonfinite"))
            charge = _yee_macrocharge(after, p, Val(true))
            for c = 1:3
                work +=
                    charge *
                    half *
                    field[c] *
                    (_exact_rational(halfv[c]) + _exact_rational(v1[c])) / 2
                after.v[c][p] = v1[c]
            end
        end
        _ace_values(old.mesh, private, st._factors)
        newsource = _ace_source(st, old.mesh, private)
        u0, u1 =
            _acq_master_potential(old.mesh, old.phi), _acq_master_potential(old.mesh, endpoint.phi)
        sourcework =
            -sum(
                i -> (u0[i] + u1[i]) * (newsource.total[i] - oldsource.total[i]) / 2,
                eachindex(u0);
                init = zero(_CompositeQ),
            )
        U0, U1, K0, K1 = _acq_energy(old.operator, old.phi),
        _acq_energy(old.operator, endpoint.phi),
        _kinetic_energy_rational(old.species),
        _kinetic_energy_rational(private)
        df, dk, dx = U1 - U0 + sourcework, K1 - K0 - work, work - sourcework
        _, r0 = _acq_residual(old.operator, old.phi, oldsource.total)
        _, r1 = _acq_residual(old.operator, endpoint.phi, newsource.total)
        df == sum(
            i -> (u0[old.mesh.free[i]] + u1[old.mesh.free[i]]) * (r1[i] - r0[i]) / 2,
            eachindex(r0);
            init = zero(_CompositeQ),
        ) || throw(ArgumentError("adaptive composite original-operator field work identity failed"))
        if old.mesh.periodic
            crossings, before_tiles, after_tiles =
                _pacc_crossings(old.mesh, transport.old_positions, private, transport.winding)
            append!(transport.tiles_before, before_tiles)
            append!(transport.tiles_after, after_tiles)
        else
            crossings = _acc_crossings(old.species, private, old.mesh)
        end
        report = AdaptiveCompositeESStepReport{D,T}(
            tau,
            st.time,
            time,
            st.step,
            count,
            st.epoch_start_time,
            st.epoch_start_step,
            old.mesh.refinements,
            K0,
            K1,
            U0,
            U1,
            sourcework,
            work,
            df,
            dk,
            dx,
            U1 - U0 + K1 - K0,
            endpoint.solve_report.max_residual,
            crossings,
        )
        state = _AdaptiveCompositePhysical{D,T}(
            old.mesh,
            old.operator,
            private,
            endpoint.phi,
            endpoint.bins,
            endpoint.solve_report,
            report,
            transport,
            endpoint.background_density,
        )
        phase = merge(
            _ace_phase(st),
            (
                time = time,
                step = count,
                last_dt = tau,
                field_defect = st.field_defect + df,
                kick_defect = st.kick_defect + dk,
                trajectory_defect = st.trajectory_defect + dx,
                physical_defect = st.physical_defect + report.total_defect,
            ),
        )
        phi = copy(state.phi)
        _ace_physical_check(st, state)
        _ace_phase_check(st, state, nothing, phase)
        _ace_prepared(st, state, species, phi, nothing)
        for (public, prepared) in zip(species, private)
            foreach(copyto!, public.x, prepared.x)
            foreach(copyto!, public.v, prepared.v)
        end
        return _ace_publish!(st, state, phi, st._owners, st._components, st._counts, nothing, phase)
    finally
        _ace_release(st)
    end
end
function _ace_regrid_report_check(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    state,
    phase,
) where {D,T}
    r = phase.last_regrid_report
    r isa AdaptiveCompositeESRegridReport{D,T} &&
        r.regrid == phase.regrid_count &&
        r.step == phase.epoch_start_step &&
        isequal(r.time, phase.epoch_start_time) &&
        isequal(r.after_refinements, state.mesh.refinements) ||
        throw(ArgumentError("adaptive composite regrid receipt chronology differs"))
    before = CartesianAMRHierarchy(
        st.mesh.nbase,
        st.mesh.L,
        r.before_refinements;
        periodic = st.mesh.periodic,
    )
    r.before_leaves == length(before.addresses) &&
        r.after_leaves == length(state.mesh.addresses) &&
        r.before_dofs == length(before.free) &&
        r.after_dofs == length(state.mesh.free) &&
        !isequal(r.before_refinements, r.after_refinements) ||
        throw(ArgumentError("adaptive composite regrid topology receipt differs"))
    all(
            _sce_q_finite,
            (
                r.field_before,
                r.field_after,
                r.kinetic_before,
                r.kinetic_after,
                r.mesh_work,
                r.max_residual,
            ),
        ) &&
        all(
            >=(zero(_CompositeQ)),
            (r.field_before, r.field_after, r.kinetic_before, r.kinetic_after, r.max_residual),
        ) &&
        r.mesh_work == r.field_after - r.field_before &&
        r.kinetic_before == r.kinetic_after ||
        throw(ArgumentError("adaptive composite mesh-work receipt differs"))
    _ace_selection_check(r.selection, T)
    if r.selection.kind === :multilevel_curvature
        r.selection.max_level >= maximum(first, state.mesh.addresses) &&
            r.selection.refine_tags <= r.before_leaves &&
            r.selection.coarsened <= length(r.before_refinements) &&
            r.selection.balance_refined <= length(r.after_refinements) ||
            throw(ArgumentError("adaptive composite selection/topology bounds differ"))
        _aca_denominators(state.mesh.nbase, r.selection.max_level)
    end
    _ace_inventory_check(r.before_inventory, D; periodic = st.mesh.periodic)
    _ace_inventory_check(r.after_inventory, D; periodic = st.mesh.periodic)
    preserved =
        st.mesh.periodic ?
        (:particle, :declared_background, :neutralizing_charge, :background, :source_total) :
        (:particle, :inside, :tail, :strata, :background)
    all(
        name ->
            isequal(getproperty(r.before_inventory, name), getproperty(r.after_inventory, name)),
        preserved,
    ) || throw(ArgumentError("adaptive composite regrid changed physical source inventory"))
    if phase.step == phase.epoch_start_step
        source = _ace_source(st, state.mesh, state.species)
        isequal(r.after_inventory, _ace_inventory(state, source)) &&
            r.field_after == _acq_energy(state.operator, state.phi) &&
            r.kinetic_after == _kinetic_energy_rational(state.species) &&
            r.max_residual == state.solve_report.max_residual || throw(
            ArgumentError(
                "adaptive composite regrid ending receipt differs from current physical state",
            ),
        )
    end
    return nothing
end
function _ace_regrid!(
    st::AdaptiveCompositeElectrostaticPIC{D,T},
    species,
    mesh,
    selection,
) where {D,T}
    mesh isa CartesianAMRHierarchy{D,T} ||
        throw(ArgumentError("adaptive composite regrid dimension/precision differs"))
    _camr_validate(mesh)
    mesh.periodic == st.mesh.periodic &&
        mesh.nbase == st.mesh.nbase &&
        isequal(mesh.L, st.mesh.L) ||
        throw(ArgumentError("adaptive composite regrid must preserve original base/domain"))
    isequal(mesh.refinements, st.mesh.refinements) && return st
    count = Base.checked_add(st.regrid_count, 1)
    owned = CartesianAMRHierarchy(mesh.nbase, mesh.L, mesh.refinements; periodic = mesh.periodic)
    old = st._state
    private = _ace_copy(old.species, Val(D), T)
    state = _ace_physical(st, owned, private)
    before, after = _ace_inventory(old, _ace_source(st, old.mesh, old.species)),
    _ace_inventory(state, _ace_source(st, owned, private))
    U0, U1, K0, K1 = _acq_energy(old.operator, old.phi),
    _acq_energy(state.operator, state.phi),
    _kinetic_energy_rational(old.species),
    _kinetic_energy_rational(private)
    report = AdaptiveCompositeESRegridReport{D,T}(
        st.time,
        st.step,
        count,
        old.mesh.refinements,
        owned.refinements,
        length(old.mesh.addresses),
        length(owned.addresses),
        length(old.mesh.free),
        length(owned.free),
        U0,
        U1,
        K0,
        K1,
        U1 - U0,
        before,
        after,
        state.solve_report.max_residual,
        selection,
    )
    history =
        old.step_report === nothing ? st._history :
        (
            epoch_start_time = st.epoch_start_time,
            epoch_start_step = st.epoch_start_step,
            global_time = st.time,
            global_step = st.step,
            last_dt = st.last_dt,
            state = _ace_physical_payload(old),
        )
    history = history === nothing ? nothing : deepcopy(history)
    phase = merge(
        _ace_phase(st),
        (
            epoch_start_time = st.time,
            epoch_start_step = st.step,
            regrid_count = count,
            mesh_work = st.mesh_work + U1 - U0,
            last_regrid_report = report,
        ),
    )
    phi = copy(state.phi)
    _ace_physical_check(st, state)
    _ace_phase_check(st, state, history, phase)
    _ace_prepared(st, state, species, phi, history)
    return _ace_publish!(st, state, phi, st._owners, st._components, st._counts, history, phase)
end
function regrid!(st::AdaptiveCompositeElectrostaticPIC, species::AbstractVector, mesh)
    _ace_acquire(st)
    try
        _ace_ready(st, species)
        return _ace_regrid!(st, species, mesh, _ace_manual_selection())
    finally
        _ace_release(st)
    end
end
function _ace_inventory_check(inventory, D; periodic = false)
    if periodic
        names = (
            :particle,
            :declared_background,
            :neutralizing_charge,
            :background,
            :source_total,
            :gauge_residual,
            :full_residual_sum,
        )
        inventory isa NamedTuple &&
            keys(inventory) == names &&
            all(_sce_q_finite, values(inventory)) &&
            inventory.background == inventory.declared_background + inventory.neutralizing_charge &&
            inventory.particle + inventory.background == inventory.source_total == 0 &&
            iszero(inventory.full_residual_sum) ||
            throw(ArgumentError("periodic composite inventory schema/closure differs"))
        return nothing
    end
    names = (
        :particle,
        :inside,
        :tail,
        :strata,
        :background,
        :free_conductor,
        :net_wall,
        :free_residual_sum,
    )
    inventory isa NamedTuple &&
        keys(inventory) == names &&
        inventory.strata isa Tuple &&
        length(inventory.strata) == 3^D ||
        throw(ArgumentError("adaptive composite inventory schema differs"))
    all(
            _sce_q_finite,
            (
                inventory.particle,
                inventory.inside,
                inventory.tail,
                inventory.strata...,
                inventory.background,
                inventory.free_conductor,
                inventory.net_wall,
                inventory.free_residual_sum,
            ),
        ) &&
        inventory.inside == inventory.strata[(3^D+1)÷2] &&
        inventory.particle ==
        inventory.inside + inventory.tail ==
        sum(inventory.strata; init = zero(_CompositeQ)) &&
        inventory.net_wall == inventory.free_conductor + inventory.tail &&
        inventory.particle +
        inventory.background +
        inventory.free_conductor +
        inventory.free_residual_sum == 0 ||
        throw(ArgumentError("adaptive composite inventory identities differ"))
    return nothing
end
