# Synchronous periodic hybrid: one concrete joint particle-velocity/B solve.
# Canonical Boris, moments, gather, Ohm and Fourier owners remain unchanged.
struct _ImplicitHybridControls{T}
    rtol::T
    atol::T
    max_newton::Int
    max_linear::Int
    memory::Int
    max_backtrack::Int
    max_jv_shrink::Int
end

struct ImplicitHybridSolveReport{T}
    step::Int
    time::T
    dt::T
    initial_residual::T
    tolerance::T
    velocity_residual::T
    magnetic_residual::T
    velocity_scale::T
    magnetic_scale::T
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
    electron_free_energy_before::T
    electron_free_energy_after::T
    pressure_work::T
    continuity_defect::T
    free_energy_defect::Union{Nothing,T}
    range_status::Symbol
end

"A bounded joint implicit solve failed before public state publication."
struct ImplicitHybridConvergenceError{T} <: Exception
    reason::Symbol
    velocity_residual::T
    magnetic_residual::T
    tolerance::T
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
    cause::Union{Nothing,CapturedException}
end

function Base.showerror(io::IO, err::ImplicitHybridConvergenceError)
    print(io, "implicit hybrid ", err.reason, " after ", err.newton_iterations,
        " Newton iterations: velocity/B residuals ", err.velocity_residual, "/",
        err.magnetic_residual, ", tolerance ", err.tolerance)
    err.cause === nothing || (print(io, "; "); showerror(io, err.cause.ex))
end

# Only these owned trial-domain checks permit bounded probe/line-search shrink.
# Exceptions from canonical physics/FFT/storage code are never broadly caught.
struct _ImplicitHybridTrialError <: Exception
    reason::Symbol
end
Base.showerror(io::IO, err::_ImplicitHybridTrialError) =
    print(io, "implicit hybrid trial is outside the admitted ", err.reason, " domain")

mutable struct _ImplicitHybridWorkspace{D,T}
    old::HybridFields{D,T}
    midpoint::HybridFields{D,T}
    endpoint::HybridFields{D,T}
    old_species::Vector{_ElectrostaticHostParticle{D,T}}
    midpoint_species::Vector{_ElectrostaticHostParticle{D,T}}
    endpoint_species::Vector{_ElectrostaticHostParticle{D,T}}
    moment::NTuple{3,Array{T,D}}
    density_scratch::Array{T,D}
    moment_scratch::NTuple{3,Array{T,D}}
    electron_velocity::NTuple{3,Array{T,D}}
    divergence::Array{T,D}
    rate::NTuple{3,Array{T,D}}
    Ep::NTuple{3,Vector{T}}
    Bp::NTuple{3,Vector{T}}
    z::Vector{T}
    residual::Vector{T}
    probe::Vector{T}
    probe_residual::Vector{T}
    trial::Vector{T}
    trial_residual::Vector{T}
    rhs::Vector{T}
    zero::Vector{T}
    action::Vector{T}
    solver::Krylov.GmresWorkspace{T,T,Vector{T}}
    h::T
    velocity_scale::T
    magnetic_scale::T
    initial_residual::T
    tolerance::T
    velocity_residual::T
    magnetic_residual::T
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
end

mutable struct ImplicitHybridStepper{D,T,G,M}
    const g::G
    const model::M
    const shape::CIC
    const fields::HybridFields{D,T}
    const owners::Vector{_ElectrostaticHostParticle{D,T}}
    const _species_scalars::Vector{Tuple{T,T}}
    const _population_counts::Vector{Int}
    const controls::_ImplicitHybridControls{T}
    const work::_ImplicitHybridWorkspace{D,T}
    const initialized::Base.RefValue{Bool}
    const time::Base.RefValue{T}
    const step::Base.RefValue{Int}
    const last_dt::Base.RefValue{T}
    last_report::Union{Nothing,ImplicitHybridSolveReport{T}}
    busy::Bool
end

function _implicit_hybrid_mode(::Type{T}) where {T}
    T in (Float32, Float64) && Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() ||
        throw(ArgumentError("implicit hybrid requires CPU Float32/64 nearest rounding and gradual underflow"))
    return nothing
end

function _implicit_hybrid_model(model::HybridModel, ::Type{T}) where {T}
    model.closure isa IsothermalElectrons{T} ||
        throw(ArgumentError("implicit hybrid requires isothermal electrons at grid precision"))
    _require_finite_nonnegative_real("implicit hybrid Te", model.closure.Te, T)
    iszero(model.η) && iszero(model.ηH) && iszero(model.de2) ||
        throw(ArgumentError("implicit hybrid requires eta = etaH = de2 = 0"))
    return _require_finite_positive_real("implicit hybrid nfloor", model.nfloor, T)
end

function _implicit_hybrid_nonnegative_limit(name, value)
    value isa Integer && 0 <= value <= typemax(Int) ||
        throw(ArgumentError("$name must be a nonnegative Int-sized integer"))
    return Int(value)
end

function _implicit_hybrid_controls(::Type{T}, rtol, atol, max_newton, max_linear,
    memory, max_backtrack, max_jv_shrink) where {T}
    r = _require_finite_nonnegative_real("nonlinear_rtol", rtol, T)
    a = _require_finite_nonnegative_real("nonlinear_atol", atol, T)
    r < one(T) && (r > zero(T) || a > zero(T)) ||
        throw(ArgumentError("implicit hybrid requires rtol < 1 and a positive tolerance"))
    (iszero(rtol) || !iszero(r)) && (iszero(atol) || !iszero(a)) ||
        throw(ArgumentError("implicit hybrid positive tolerances must survive conversion"))
    return _ImplicitHybridControls(r, a,
        _implicit_hybrid_nonnegative_limit("max_newton", max_newton),
        _require_positive_intlike("max_linear", max_linear),
        _require_positive_intlike("memory", memory),
        _implicit_hybrid_nonnegative_limit("max_backtrack", max_backtrack),
        _implicit_hybrid_nonnegative_limit("max_jv_shrink", max_jv_shrink))
end

function _implicit_hybrid_configuration(st::ImplicitHybridStepper{D,T}) where {D,T}
    nf = _implicit_hybrid_model(st.model, T)
    c = st.controls
    return (D = D, T = T, ncell = st.g.n, L = st.g.L, dx = st.g.dx,
        model = st.model, nfloor = nf, shape = :cic, mobility = :all_mobile,
        species_scalars = Tuple(st._species_scalars),
        residual_policy = :boris_map_faraday_v1,
        quadrature_policy = :lifted_single_midpoint_v1,
        phase_policy = :synchronous_integer,
        moment_policy = :canonical_charge_weighted_collection_v1,
        derivative_policy = :represented_odd_nyquist_zero_v1,
        density_policy = :strict_reject,
        field_policy = :complete_public_hybrid_fields_v1,
        norm_policy = :fixed_reference_block_infinity_v1,
        solver_policy = :newton_forward_jv_gmres_identity_armijo_v1,
        nonlinear_rtol = c.rtol, nonlinear_atol = c.atol,
        max_newton = c.max_newton, max_linear = c.max_linear, memory = c.memory,
        max_backtrack = c.max_backtrack, max_jv_shrink = c.max_jv_shrink)
end

function _implicit_hybrid_phase_values(initialized, time::T, step::Int, last_dt::T,
    report) where {T}
    initialized isa Bool && isfinite(time) && time >= zero(T) && step >= 0 ||
        throw(ArgumentError("implicit hybrid phase/clock is invalid"))
    isfinite(last_dt) && _stored_nonnegative(last_dt) ||
        throw(ArgumentError("implicit hybrid history must be finite and nonnegative"))
    if step == 0
        _particle_stored_zero(time) && _particle_stored_zero(last_dt) && report === nothing ||
            throw(ArgumentError("unstepped implicit hybrid requires zero clock/history and no solve report"))
    else
        initialized && time > zero(T) && _stored_positive(last_dt) ||
            throw(ArgumentError("advanced implicit hybrid requires initialized positive clock/history"))
        _require_resolved_dt("implicit hybrid history", last_dt)
        _implicit_hybrid_validate_report(report, T, time, step, last_dt)
    end
    return nothing
end

function _implicit_hybrid_validate_report(r, ::Type{T}, time, step, h) where {T}
    r isa ImplicitHybridSolveReport{T} && r.step == step && isequal(r.time, time) && isequal(r.dt, h) ||
        throw(ArgumentError("implicit hybrid report must match its completed clock/history"))
    all(x -> isfinite(x) && _stored_nonnegative(x),
        (r.initial_residual, r.tolerance, r.velocity_residual, r.magnetic_residual)) &&
        all(x -> isfinite(x) && x >= one(T), (r.velocity_scale, r.magnetic_scale)) &&
        max(r.velocity_residual, r.magnetic_residual) <= r.tolerance &&
        (iszero(r.initial_residual) || r.tolerance > zero(T)) ||
        throw(ArgumentError("implicit hybrid report has invalid residual metrics"))
    all(x -> x >= 0, (r.newton_iterations, r.linear_iterations, r.backtracks)) &&
        r.residual_evaluations > 0 && !isnan(r.electron_free_energy_before) &&
        !isnan(r.electron_free_energy_after) && !isnan(r.pressure_work) &&
        !isnan(r.continuity_defect) && r.continuity_defect >= zero(T) ||
        throw(ArgumentError("implicit hybrid report has invalid counters or work values"))
    finite_energy = isfinite(r.electron_free_energy_before) && isfinite(r.electron_free_energy_after)
    (finite_energy ? r.range_status === :represented_reservoir &&
        r.free_energy_defect isa T && !isnan(r.free_energy_defect) :
        r.range_status === :unresolved_free_energy_range && r.free_energy_defect === nothing) ||
        throw(ArgumentError("implicit hybrid report range status does not match its reservoirs"))
    return nothing
end

function _implicit_hybrid_report_policy(st::ImplicitHybridStepper{D,T}, r) where {D,T}
    r === nothing && return nothing
    c = st.controls
    expected = _from_exact(T, _exact_rational(c.atol) +
        _exact_rational(c.rtol) * _exact_rational(r.initial_residual))
    updates = BigInt(r.newton_iterations)
    isequal(r.tolerance, expected) && r.newton_iterations <= c.max_newton &&
        r.linear_iterations <= updates * c.max_linear &&
        r.backtracks <= updates * c.max_backtrack &&
        r.residual_evaluations >= updates + 2 ||
        throw(ArgumentError("implicit hybrid report differs from its fixed solver policy"))
    # Necessary metadata relations do not re-solve the preceding interval.
    return nothing
end

function _acquire_implicit_hybrid!(st::ImplicitHybridStepper)
    st.busy && throw(ArgumentError("implicit hybrid does not permit concurrent or reentrant operations"))
    st.busy = true
    return nothing
end

_implicit_hybrid_field_arrays(f::HybridFields) =
    (f.n, f.ui..., f.B..., f.E..., f.J..., f.pe, f.gradp..., f.ninv,
        f.lapJ..., f.pforce..., f.cgl_dp)

function _implicit_hybrid_failure(st::ImplicitHybridStepper{D,T}, reason, cause = nothing) where {D,T}
    w = st.work
    return ImplicitHybridConvergenceError{T}(reason, w.velocity_residual,
        w.magnetic_residual, w.tolerance, w.newton_iterations, w.linear_iterations,
        w.residual_evaluations, w.backtracks, cause)
end

_implicit_hybrid_particle_arrays(species) =
    Tuple(a for ps in species for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))

function _implicit_hybrid_work_arrays(w::_ImplicitHybridWorkspace)
    return (_implicit_hybrid_field_arrays(w.old)...,
        _implicit_hybrid_field_arrays(w.midpoint)...,
        _implicit_hybrid_field_arrays(w.endpoint)...,
        w.old_species, w.midpoint_species, w.endpoint_species,
        _implicit_hybrid_particle_arrays(w.old_species)...,
        _implicit_hybrid_particle_arrays(w.midpoint_species)...,
        _implicit_hybrid_particle_arrays(w.endpoint_species)...,
        w.moment..., w.density_scratch, w.moment_scratch..., w.electron_velocity...,
        w.divergence, w.rate..., w.Ep..., w.Bp..., w.z, w.residual, w.probe,
        w.probe_residual, w.trial, w.trial_residual, w.rhs, w.zero, w.action,
        _gmres_workspace_arrays(w.solver)...)
end

function _implicit_hybrid_storage(st::ImplicitHybridStepper{D,T}) where {D,T}
    g, w = st.g, st.work
    ns = length(st.owners)
    axes(st._species_scalars, 1) == Base.OneTo(ns) &&
        axes(st._population_counts, 1) == Base.OneTo(ns) ||
        throw(DimensionMismatch("implicit hybrid binding storage differs from its ordered species"))
    for f in (st.fields, w.old, w.midpoint, w.endpoint)
        for (name, value) in pairs(_hybrid_field_state(f))
            expected = name in (:pforce, :cgl_dp) ? ntuple(_ -> 0, D) : g.n
            for a in (value isa Tuple ? value : (value,))
                size(a) == expected || throw(DimensionMismatch("implicit hybrid $name storage differs from its grid role"))
            end
        end
    end
    refs = (st.initialized, st.time, st.step, st.last_dt, st.fields.floor_count,
        w.old.floor_count, w.midpoint.floor_count, w.endpoint.floor_count)
    for i in eachindex(refs), j = 1:(i-1)
        refs[i] === refs[j] && throw(ArgumentError("implicit hybrid clock and field references must be independent"))
    end
    for stages in (w.old_species, w.midpoint_species, w.endpoint_species)
        length(stages) == ns || throw(DimensionMismatch("implicit hybrid stage species count differs"))
        for ps in stages
            _check_particle_storage_axes(ps)
        end
    end
    for a in (g.cbuf, g.tbuf, g.abuf, g.sbuf)
        size(a) == g.n || throw(DimensionMismatch("implicit hybrid Fourier work dimensions differ"))
    end
    for d = 1:D, a in (g.midx[d], g.kfull[d], g.ik[d], g.kvec[d])
        axes(a, 1) == Base.OneTo(g.n[d]) ||
            throw(DimensionMismatch("implicit hybrid Fourier metadata dimensions differ"))
    end
    for a in (w.moment..., w.density_scratch, w.moment_scratch...,
        w.electron_velocity..., w.divergence, w.rate...)
        size(a) == g.n || throw(DimensionMismatch("implicit hybrid grid work dimensions differ"))
    end
    n = length(w.z)
    _validate_gmres_workspace_layout(w.solver, n)
    all(a -> length(a) == n, _gmres_workspace_vectors(w.solver)) ||
        throw(DimensionMismatch("implicit hybrid GMRES vector lengths differ"))
    all(a -> length(a) == n, (w.residual, w.probe, w.probe_residual, w.trial,
        w.trial_residual, w.rhs, w.zero, w.action)) ||
        throw(DimensionMismatch("implicit hybrid nonlinear vector dimensions differ"))
    old_count = 0
    for ps in w.old_species
        old_count = Base.checked_add(old_count, length(ps.weight))
    end
    n == _implicit_hybrid_dimension(g, old_count) &&
        all(a -> length(a) == old_count, (w.Ep..., w.Bp...)) ||
        throw(DimensionMismatch("implicit hybrid particle work dimensions differ"))
    # Writers, cached bindings and particle containers are separate roles even
    # when empty. Read-only Fourier metadata may share only with other metadata.
    arrays = (_implicit_hybrid_field_arrays(st.fields)..., st.owners,
        st._species_scalars, st._population_counts,
        _implicit_hybrid_particle_arrays(st.owners)...,
        _implicit_hybrid_work_arrays(w)..., g.cbuf, g.tbuf, g.abuf, g.sbuf)
    metadata = (g.midx..., g.kfull..., g.ik..., g.kvec...)
    for i in eachindex(arrays)
        for j = 1:(i-1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) &&
                throw(ArgumentError("implicit hybrid public/private storage roles must be independent"))
        end
        for a in metadata
            _particle_arrays_mightalias(arrays[i], a) &&
                throw(ArgumentError("implicit hybrid writers must not alias Fourier metadata"))
        end
    end
    return nothing
end

function _implicit_hybrid_binding(st::ImplicitHybridStepper{D,T}, species;
    values::Bool = true, counts::Bool = values) where {D,T}
    _implicit_hybrid_mode(T)
    _implicit_hybrid_model(st.model, T)
    c = st.controls
    _implicit_hybrid_controls(T, c.rtol, c.atol, c.max_newton, c.max_linear,
        c.memory, c.max_backtrack, c.max_jv_shrink)
    _implicit_hybrid_storage(st)
    species isa Vector{<:ParticleSet{D,T}} && length(species) == length(st.owners) ||
        throw(ArgumentError("implicit hybrid requires its fixed ordinary ordered species Vector"))
    for (s, ps) in enumerate(species)
        ps isa _ElectrostaticHostParticle{D,T} && ps === st.owners[s] ||
            throw(ArgumentError("implicit hybrid owner identity/order and Vector-backed storage are fixed"))
        _check_particle_storage_axes(ps)
        q, m = st._species_scalars[s]
        _validated_charge_mass(q, m, T)
        _stored_positive(q) || throw(ArgumentError("implicit hybrid cached ion charge must be positive"))
        values && !isequal((ps.q, ps.m), (q, m)) &&
            throw(ArgumentError("implicit hybrid charge/mass binding changed"))
        counts && st._population_counts[s] != length(ps.weight) &&
            throw(ArgumentError("implicit hybrid population counts changed"))
    end
    values && _validate_hybrid_species(species)
    return nothing
end

function _implicit_hybrid_boundary(st::ImplicitHybridStepper{D,T}) where {D,T}
    st.initialized[] || throw(ArgumentError("initialize implicit hybrid before ready use"))
    _implicit_hybrid_phase_values(st.initialized[], st.time[], st.step[], st.last_dt[], st.last_report)
    _implicit_hybrid_report_policy(st, st.last_report)
    nf = _implicit_hybrid_model(st.model, T)
    all(a -> all(isfinite, a), _implicit_hybrid_field_arrays(st.fields)) &&
        all(n -> n >= nf, st.fields.n) && st.fields.floor_count[] == 0 ||
        throw(ArgumentError("implicit hybrid public fields must be finite with inactive density floor"))
    return nothing
end

function _implicit_hybrid_dimension(g, particles::Int)
    cells = foldl(Base.checked_mul, g.n; init = 1)
    return Base.checked_mul(3, Base.checked_add(particles, cells))
end

function _implicit_hybrid_solver(::Type{T}, n, memory) where {T}
    s = Krylov.GmresWorkspace(n, n, Vector{T}; memory = min(memory, n))
    resize!(s.Δx, n)
    return s
end

function _implicit_hybrid_workspace(g::FourierGrid{D,T}, particles, controls) where {D,T}
    vector() = ntuple(_ -> zeros(T, g.n), 3)
    particle_vector() = ntuple(_ -> zeros(T, particles), 3)
    n = _implicit_hybrid_dimension(g, particles)
    return _ImplicitHybridWorkspace(HybridFields{D,T}(g.n), HybridFields{D,T}(g.n),
        HybridFields{D,T}(g.n), _ElectrostaticHostParticle{D,T}[],
        _ElectrostaticHostParticle{D,T}[], _ElectrostaticHostParticle{D,T}[],
        vector(), zeros(T, g.n), vector(), vector(), zeros(T, g.n), vector(),
        particle_vector(), particle_vector(),
        zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n),
        zeros(T, n), zeros(T, n), zeros(T, n), zeros(T, n),
        _implicit_hybrid_solver(T, n, controls.memory),
        zero(T), one(T), one(T), zero(T), zero(T), zero(T), zero(T), 0, 0, 0, 0)
end

function _implicit_hybrid_reserve!(st::ImplicitHybridStepper{D,T}, species) where {D,T}
    w = st.work
    total = 0
    for ps in species
        total = Base.checked_add(total, length(ps.weight))
    end
    n = _implicit_hybrid_dimension(st.g, total)
    for stages in (w.old_species, w.midpoint_species, w.endpoint_species)
        _pic_prepare_species_stages!(stages, species)
    end
    for a in (w.Ep..., w.Bp...)
        resize!(a, total)
    end
    for a in (w.z, w.residual, w.probe, w.probe_residual, w.trial, w.trial_residual,
        w.rhs, w.zero, w.action)
        resize!(a, n)
    end
    length(w.solver.x) == n || (w.solver = _implicit_hybrid_solver(T, n, st.controls.memory))
    return nothing
end

"""
    ImplicitHybridStepper(g, model, CIC(), species; nonlinear_rtol=sqrt(eps(T)),
        nonlinear_atol=64eps(T), max_newton=20, max_linear=200, memory=20,
        max_backtrack=12, max_jv_shrink=8)

Allocate a synchronous periodic CPU Float32/64 isothermal massless hybrid.
Every ordinary Vector-backed positive-ion population is mobile. Seed `fields.B`
and call `init!` with full-time velocities. The dimensionless solver defaults are
bounded policy, not convergence or physical-accuracy guarantees. No leapfrog
conversion, orbit averaging, transport, source or timestep fallback is provided.
See `docs/IMPLICIT_HYBRID.md` for the represented residual and work limitations.
"""
function ImplicitHybridStepper(g::FourierGrid{D,T}, model::HybridModel, shape::CIC,
    species; nonlinear_rtol = sqrt(eps(T)), nonlinear_atol = 64eps(T),
    max_newton = 20, max_linear = 200, memory = 20, max_backtrack = 12,
    max_jv_shrink = 8) where {D,T}
    _implicit_hybrid_mode(T)
    _implicit_hybrid_model(model, T)
    species isa Vector{<:ParticleSet{D,T}} &&
        all(ps -> ps isa _ElectrostaticHostParticle{D,T}, species) ||
        throw(ArgumentError("implicit hybrid requires an ordinary Vector of Vector-backed matching species"))
    total = _validate_hybrid_species(species)
    controls = _implicit_hybrid_controls(T, nonlinear_rtol, nonlinear_atol,
        max_newton, max_linear, memory, max_backtrack, max_jv_shrink)
    owners = _ElectrostaticHostParticle{D,T}[ps for ps in species]
    st = ImplicitHybridStepper(g, model, shape, HybridFields{D,T}(g.n), owners,
        Tuple{T,T}[(ps.q, ps.m) for ps in owners], Int[length(ps.weight) for ps in owners],
        controls, _implicit_hybrid_workspace(g, total, controls),
        Ref(false), Ref(zero(T)), Ref(0), Ref(zero(T)), nothing, false)
    _implicit_hybrid_reserve!(st, owners)
    _implicit_hybrid_binding(st, owners)
    return st
end

@inline function _implicit_hybrid_mean(a::T, b::T) where {T}
    value = _from_exact(T, (_exact_rational(a) + _exact_rational(b)) / 2)
    isfinite(value) || throw(_ImplicitHybridTrialError(:affine_range))
    return value
end

function _implicit_hybrid_moments!(st::ImplicitHybridStepper{D,T}, f, species;
    keep_momentum::Bool = false) where {D,T}
    w = st.work
    _zero_ion_moments!(f.n, f.ui)
    for ps in species
        _accumulate_charge_moments!(f.n, f.ui, ps, st.g, st.shape,
            w.density_scratch, w.moment_scratch, nothing)
    end
    nf = T(st.model.nfloor)
    all(n -> isfinite(n) && n >= nf, f.n) ||
        throw(_ImplicitHybridTrialError(:strict_density))
    all(a -> all(isfinite, a), f.ui) || throw(_ImplicitHybridTrialError(:momentum_range))
    keep_momentum && foreach(copyto!, w.moment, f.ui)
    _normalize_ion_moments!(f.n, f.ui, nf)
    all(a -> all(isfinite, a), f.ui) || throw(_ImplicitHybridTrialError(:velocity_range))
    return nothing
end

function _implicit_hybrid_ohm!(st::ImplicitHybridStepper, f)
    # These are private fixed-zero inputs to the inactive hyperresistive row.
    foreach(a -> fill!(a, zero(eltype(a))), f.lapJ)
    ohms_law!(f, st.model, st.g)
    all(a -> all(isfinite, a), _implicit_hybrid_field_arrays(f)) && f.floor_count[] == 0 ||
        throw(_ImplicitHybridTrialError(:field_range))
    return nothing
end

function _implicit_hybrid_residual!(out::Vector{T}, st::ImplicitHybridStepper{D,T},
    z::Vector{T}) where {D,T}
    w, g = st.work, st.g
    w.residual_evaluations = Base.checked_add(w.residual_evaluations, 1)
    # An incomplete attempted block has no finite norm; do not report a stale
    # accepted/base norm as if that invalid trial had completed.
    w.velocity_residual = T(Inf)
    w.magnetic_residual = T(Inf)
    all(isfinite, z) || throw(_ImplicitHybridTrialError(:iterate_range))
    i = 0
    for (mid, old) in zip(w.midpoint_species, w.old_species)
        for c = 1:3, p in eachindex(old.weight)
            i += 1
            mid.v[c][p] = _implicit_hybrid_mean(old.v[c][p], z[i])
        end
        for d = 1:D, p in eachindex(old.weight)
            lifted = _coordinate_eval(Val(true), old.x[d][p], mid.v[d][p], w.h / 2)
            endpoint_lift = _coordinate_eval(Val(true), old.x[d][p], mid.v[d][p], w.h)
            isfinite(lifted) && isfinite(endpoint_lift) ||
                throw(_ImplicitHybridTrialError(:lift_range))
            mid.x[d][p] = _periodic_coordinate(lifted, zero(T), g.L[d], g.L[d])
        end
    end
    for c = 1:3, I in eachindex(w.midpoint.B[c])
        i += 1
        w.midpoint.B[c][I] = _implicit_hybrid_mean(w.old.B[c][I], z[i])
    end
    _implicit_hybrid_moments!(st, w.midpoint, w.midpoint_species; keep_momentum = true)
    _implicit_hybrid_ohm!(st, w.midpoint)
    i, offset = 0, 0
    vmax = zero(T)
    Q = _exact_rational
    for (mid, old) in zip(w.midpoint_species, w.old_species)
        n = length(old.weight)
        indices = (offset + 1):(offset + n)
        Ep = ntuple(c -> view(w.Ep[c], indices), 3)
        Bp = ntuple(c -> view(w.Bp[c], indices), 3)
        gather_vector!(Ep, w.midpoint.E, mid, g, st.shape)
        gather_vector!(Bp, w.midpoint.B, mid, g, st.shape)
        all(a -> all(isfinite, a), (Ep..., Bp...)) ||
            throw(_ImplicitHybridTrialError(:gather_range))
        factor = _validated_boris_factor(old)
        # The unknown uses component-major storage within each population.
        for p in eachindex(old.weight)
            kicked = _kick_eval(Val(true), old.v[1][p], old.v[2][p], old.v[3][p],
                Ep[1][p], Ep[2][p], Ep[3][p], Bp[1][p], Bp[2][p], Bp[3][p], factor, w.h)
            all(isfinite, kicked) || throw(_ImplicitHybridTrialError(:boris_range))
            for c = 1:3
                index = i + (c - 1) * n + p
                value = _from_exact(T, (Q(z[index]) - Q(kicked[c])) / Q(w.velocity_scale))
                isfinite(value) || throw(_ImplicitHybridTrialError(:residual_range))
                out[index] = value
                vmax = max(vmax, abs(value))
            end
        end
        i += 3n
        offset += n
    end
    w.velocity_residual = vmax
    faraday_rhs!(w.rate, w.midpoint.E, g)
    all(a -> all(isfinite, a), w.rate) || throw(_ImplicitHybridTrialError(:faraday_range))
    bmax = zero(T)
    for c = 1:3, I in eachindex(w.old.B[c])
        i += 1
        value = _from_exact(T,
            (Q(z[i]) - Q(w.old.B[c][I]) - Q(w.h) * Q(w.rate[c][I])) / Q(w.magnetic_scale))
        isfinite(value) || throw(_ImplicitHybridTrialError(:residual_range))
        out[i] = value
        bmax = max(bmax, abs(value))
    end
    w.magnetic_residual = bmax
    return vmax, bmax
end

struct _ImplicitHybridJacobian{S}
    state::S
end
Base.size(A::_ImplicitHybridJacobian) = (length(A.state.work.z), length(A.state.work.z))
Base.eltype(::_ImplicitHybridJacobian{<:ImplicitHybridStepper{D,T}}) where {D,T} = T

function LinearAlgebra.mul!(out::Vector{T}, A::_ImplicitHybridJacobian,
    p::Vector{T}) where {T}
    st, w = A.state, A.state.work
    all(isfinite, p) || throw(_implicit_hybrid_failure(st, :nonfinite_jv_direction))
    all(iszero, p) && return fill!(out, zero(T))
    pnorm = BLAS.nrm2(p)
    isfinite(pnorm) && pnorm > zero(T) || throw(_implicit_hybrid_failure(st, :jv_scale))
    nv = 3length(w.Ep[1])
    zscale = one(T)
    for i in eachindex(w.z)
        scale = i <= nv ? w.velocity_scale : w.magnetic_scale
        zscale = max(zscale, abs(w.z[i] / scale))
    end
    Q = _exact_rational
    epsilon = _from_exact(T, Q(sqrt(eps(T))) * Q(zscale) / Q(pnorm))
    isfinite(epsilon) && epsilon > zero(T) || throw(_implicit_hybrid_failure(st, :jv_scale))
    cause = nothing
    for attempt = 0:st.controls.max_jv_shrink
        changed = false
        try
            for i in eachindex(w.z)
                scale = i <= nv ? w.velocity_scale : w.magnetic_scale
                value = _from_exact(T, Q(w.z[i]) + Q(epsilon) * Q(scale) * Q(p[i]))
                isfinite(value) || throw(_ImplicitHybridTrialError(:probe_range))
                changed |= value != w.z[i]
                w.probe[i] = value
            end
            changed || throw(_implicit_hybrid_failure(st, :unresolved_jv))
            _implicit_hybrid_residual!(w.probe_residual, st, w.probe)
            for i in eachindex(out)
                value = _from_exact(T, (Q(w.probe_residual[i]) - Q(w.residual[i])) / Q(epsilon))
                isfinite(value) || throw(_implicit_hybrid_failure(st, :jv_action_range))
                out[i] = value
            end
            return out
        catch err
            err isa _ImplicitHybridTrialError || rethrow()
            cause = CapturedException(err, catch_backtrace())
        end
        attempt == st.controls.max_jv_shrink && break
        epsilon /= 2
        epsilon > zero(T) || throw(_implicit_hybrid_failure(st, :unresolved_jv, cause))
    end
    throw(_implicit_hybrid_failure(st, :jv_domain, cause))
end

function _implicit_hybrid_linear_tick!(st)
    # Installed GMRES invokes its internal callback once per completed Arnoldi
    # iteration, including the final one. This is not a user callback surface.
    st.work.linear_iterations = Base.checked_add(st.work.linear_iterations, 1)
    return false
end

function _implicit_hybrid_set_metrics!(st, metrics)
    st.work.velocity_residual, st.work.magnetic_residual = metrics
    return max(metrics...)
end

function _implicit_hybrid_solve!(st::ImplicitHybridStepper{D,T}) where {D,T}
    w, controls = st.work, st.controls
    initial = try
        _implicit_hybrid_residual!(w.residual, st, w.z)
    catch err
        err isa _ImplicitHybridTrialError || rethrow()
        throw(_implicit_hybrid_failure(st, :initial_domain, CapturedException(err, catch_backtrace())))
    end
    rho = _implicit_hybrid_set_metrics!(st, initial)
    w.initial_residual = rho
    Q = _exact_rational
    w.tolerance = _from_exact(T, Q(controls.atol) + Q(controls.rtol) * Q(rho))
    isfinite(w.tolerance) && (iszero(rho) || w.tolerance > zero(T)) ||
        throw(_implicit_hybrid_failure(st, :tolerance_range))
    A = _ImplicitHybridJacobian(st)
    for iteration = 0:controls.max_newton
        rho <= w.tolerance && break
        iteration == controls.max_newton && throw(_implicit_hybrid_failure(st, :nonlinear_limit))
        for i in eachindex(w.rhs)
            w.rhs[i] = -w.residual[i]
        end
        rhsnorm = BLAS.nrm2(w.rhs)
        target = rhsnorm / 2
        estimated_target = target / 2
        isfinite(rhsnorm) && isfinite(target) && estimated_target > zero(T) ||
            throw(_implicit_hybrid_failure(st, :linear_tolerance_range))
        fill!(w.zero, zero(T))
        # Check the defined action at the explicit zero correction before
        # touching restarted history. No previous correction is retained.
        mul!(w.action, A, w.zero)
        for i in eachindex(w.action)
            w.action[i] = _from_exact(T, Q(w.action[i]) + Q(w.residual[i]))
        end
        zero_residual = BLAS.nrm2(w.action)
        isfinite(zero_residual) || throw(_implicit_hybrid_failure(st, :linear_initial_residual))
        delta = if zero_residual <= target
            w.zero
        else
            Krylov.gmres!(w.solver, A, w.rhs, w.zero; restart = true,
                reorthogonalization = true, atol = estimated_target, rtol = zero(T),
                itmax = controls.max_linear, callback = _ -> _implicit_hybrid_linear_tick!(st))
            Krylov.solution(w.solver)
        end
        all(isfinite, delta) || throw(_implicit_hybrid_failure(st, :nonfinite_linear_solution))
        mul!(w.action, A, delta)
        for i in eachindex(w.action)
            w.action[i] = _from_exact(T, Q(w.action[i]) + Q(w.residual[i]))
        end
        linear_residual = BLAS.nrm2(w.action)
        isfinite(linear_residual) && linear_residual <= target ||
            throw(_implicit_hybrid_failure(st, :linear_residual))
        accepted = false
        cause = nothing
        nv = 3length(w.Ep[1])
        alpha = one(T)
        for backtrack = 0:controls.max_backtrack
            changed = false
            try
                for i in eachindex(w.z)
                    scale = i <= nv ? w.velocity_scale : w.magnetic_scale
                    value = _from_exact(T, Q(w.z[i]) + Q(alpha) * Q(scale) * Q(delta[i]))
                    isfinite(value) || throw(_ImplicitHybridTrialError(:line_range))
                    changed |= value != w.z[i]
                    w.trial[i] = value
                end
                changed || throw(_implicit_hybrid_failure(st, :unresolved_line_search))
                metrics = _implicit_hybrid_residual!(w.trial_residual, st, w.trial)
                trial_rho = max(metrics...)
                armijo = _from_exact(T, (1 - Q(alpha) / 4) * Q(rho))
                if trial_rho <= w.tolerance || trial_rho <= armijo
                    copyto!(w.z, w.trial)
                    copyto!(w.residual, w.trial_residual)
                    rho = _implicit_hybrid_set_metrics!(st, metrics)
                    w.newton_iterations = Base.checked_add(w.newton_iterations, 1)
                    accepted = true
                    break
                end
            catch err
                err isa _ImplicitHybridTrialError || rethrow()
                cause = CapturedException(err, catch_backtrace())
            end
            w.backtracks = Base.checked_add(w.backtracks, 1)
            backtrack == controls.max_backtrack && break
            alpha /= 2
            alpha > zero(T) || throw(_implicit_hybrid_failure(st, :unresolved_line_search, cause))
        end
        accepted || throw(_implicit_hybrid_failure(st, :backtrack_limit, cause))
    end
    # Reconstruct every coupled quantity from immutable old snapshots and z.
    # GMRES's action/probe and the accepted trial cache are not acceptance data.
    final = try
        _implicit_hybrid_residual!(w.residual, st, w.z)
    catch err
        err isa _ImplicitHybridTrialError || rethrow()
        throw(_implicit_hybrid_failure(st, :final_domain, CapturedException(err, catch_backtrace())))
    end
    _implicit_hybrid_set_metrics!(st, final) <= w.tolerance ||
        throw(_implicit_hybrid_failure(st, :final_residual))
    return nothing
end

function _implicit_hybrid_stage!(st::ImplicitHybridStepper{D,T}, species, h::T) where {D,T}
    _implicit_hybrid_reserve!(st, species)
    w = st.work
    _copy_hybrid_fields!(w.old, st.fields)
    w.h = h
    w.velocity_scale = one(T)
    w.magnetic_scale = one(T)
    w.initial_residual = zero(T)
    w.tolerance = zero(T)
    w.velocity_residual = zero(T)
    w.magnetic_residual = zero(T)
    w.newton_iterations = 0
    w.linear_iterations = 0
    w.residual_evaluations = 0
    w.backtracks = 0
    i = 0
    for ps in w.old_species, c = 1:3, p in eachindex(ps.weight)
        i += 1
        value = ps.v[c][p]
        w.z[i] = value
        w.velocity_scale = max(w.velocity_scale, abs(value))
    end
    for c = 1:3, I in eachindex(w.old.B[c])
        i += 1
        value = w.old.B[c][I]
        w.z[i] = value
        w.magnetic_scale = max(w.magnetic_scale, abs(value))
    end
    return nothing
end

function _implicit_hybrid_endpoint!(st::ImplicitHybridStepper{D,T}) where {D,T}
    w, g = st.work, st.g
    i = 0
    for (endpoint, midpoint, old) in zip(w.endpoint_species, w.midpoint_species, w.old_species)
        for c = 1:3, p in eachindex(old.weight)
            i += 1
            endpoint.v[c][p] = w.z[i]
        end
        for d = 1:D, p in eachindex(old.weight)
            lifted = _coordinate_eval(Val(true), old.x[d][p], midpoint.v[d][p], w.h)
            isfinite(lifted) || throw(_ImplicitHybridTrialError(:lift_range))
            endpoint.x[d][p] = _periodic_coordinate(lifted, zero(T), g.L[d], g.L[d])
        end
    end
    for c = 1:3, I in eachindex(w.endpoint.B[c])
        i += 1
        w.endpoint.B[c][I] = w.z[i]
    end
    _implicit_hybrid_moments!(st, w.endpoint, w.endpoint_species)
    _implicit_hybrid_ohm!(st, w.endpoint)
    return nothing
end

function _implicit_hybrid_report(st::ImplicitHybridStepper{D,T}, time::T, step::Int) where {D,T}
    w, g = st.work, st.g
    before = electron_internal_energy(w.old.n, st.model.closure, g)
    after = electron_internal_energy(w.endpoint.n, st.model.closure, g)
    !isnan(before) && !isnan(after) || throw(_ImplicitHybridTrialError(:free_energy_nan))
    electron_velocity!(w.electron_velocity, w.midpoint.ui, w.midpoint.J,
        w.midpoint.n; nfloor = T(st.model.nfloor))
    all(a -> all(isfinite, a), w.electron_velocity) ||
        throw(_ImplicitHybridTrialError(:electron_velocity_range))
    # Use the actual canonical standalone divergence, not a reconstructed flux
    # or a new fused derivative with different range/rounding behavior.
    divergence!(w.divergence, w.moment, g)
    all(isfinite, w.divergence) || throw(_ImplicitHybridTrialError(:continuity_range))
    Q = _exact_rational
    volume = prod(Q.(g.dx))
    pressure_work_exact = zero(Q(w.h))
    continuity_exact = zero(Q(w.h))
    for I in eachindex(w.old.n)
        for d = 1:D
            pressure_work_exact -= Q(w.h) * volume *
                Q(w.electron_velocity[d][I]) * Q(w.midpoint.gradp[d][I])
        end
        value = Q(w.endpoint.n[I]) - Q(w.old.n[I]) + Q(w.h) * Q(w.divergence[I])
        continuity_exact = max(continuity_exact, abs(value))
    end
    pressure_work = _from_exact(T, pressure_work_exact)
    continuity_defect = _from_exact(T, continuity_exact)
    finite_energy = isfinite(before) && isfinite(after)
    defect = finite_energy ? _from_exact(T, Q(after) - Q(before) + pressure_work_exact) : nothing
    report = ImplicitHybridSolveReport{T}(step, time, w.h, w.initial_residual,
        w.tolerance, w.velocity_residual, w.magnetic_residual,
        w.velocity_scale, w.magnetic_scale, w.newton_iterations,
        w.linear_iterations, w.residual_evaluations, w.backtracks,
        before, after, pressure_work, continuity_defect, defect,
        finite_energy ? :represented_reservoir : :unresolved_free_energy_range)
    _implicit_hybrid_validate_report(report, T, time, step, w.h)
    _implicit_hybrid_report_policy(st, report)
    return report
end

function _implicit_hybrid_replacements(species::Vector{_ElectrostaticHostParticle{D,T}}) where {D,T}
    prepared = _ElectrostaticHostParticle{D,T}[
        _prepare_particle_replacement(ps, eachindex(ps.weight)) for ps in species]
    counts = Int[length(ps.weight) for ps in prepared]
    return prepared, counts
end

function _implicit_hybrid_commit!(st::ImplicitHybridStepper{D,T}, prepared,
    counts, fields, initialized::Bool, time::T, step::Int, last_dt::T,
    report::Union{Nothing,ImplicitHybridSolveReport{T}}) where {D,T}
    # Every leaf/type/count and field shape has already passed preparation.
    # Preserve public ParticleSet identities and field-array identities.
    for s in eachindex(st.owners)
        destination, source = st.owners[s], prepared[s]
        _publish_particle_replacement!(destination, source)
        destination.q = source.q
        destination.m = source.m
    end
    _copy_hybrid_fields!(st.fields, fields)
    copyto!(st._population_counts, counts)
    st.initialized[] = initialized
    st.time[] = time
    st.step[] = step
    st.last_dt[] = last_dt
    st.last_report = report
    return st
end

"""
    init!(st::ImplicitHybridStepper, species)

Initialize integer-time fields from the bound full-time particle state and the
seeded `st.fields.B`. Positions are wrapped on private copies. Strict density,
canonical Ohm completion and all replacements finish before publication. This
starts a new zero-clock history without a leapfrog half kick.
"""
function init!(st::ImplicitHybridStepper{D,T}, species::AbstractVector) where {D,T}
    _acquire_implicit_hybrid!(st)
    try
        _implicit_hybrid_binding(st, species)
        all(a -> all(isfinite, a), st.fields.B) ||
            throw(ArgumentError("implicit hybrid initial magnetic field must be finite"))
        _implicit_hybrid_stage!(st, species, zero(T))
        w = st.work
        for ps in w.endpoint_species, d = 1:D, p in eachindex(ps.weight)
            ps.x[d][p] = _periodic_coordinate(ps.x[d][p], zero(T), st.g.L[d], st.g.L[d])
        end
        foreach(copyto!, w.endpoint.B, st.fields.B)
        try
            _implicit_hybrid_moments!(st, w.endpoint, w.endpoint_species)
            _implicit_hybrid_ohm!(st, w.endpoint)
        catch err
            err isa _ImplicitHybridTrialError || rethrow()
            throw(_implicit_hybrid_failure(st, :initialization_domain, CapturedException(err, catch_backtrace())))
        end
        prepared, counts = _implicit_hybrid_replacements(w.endpoint_species)
        return _implicit_hybrid_commit!(st, prepared, counts, w.endpoint,
            true, zero(T), 0, zero(T), nothing)
    finally
        st.busy = false
    end
end

"""
    step!(st::ImplicitHybridStepper, species, dt)

Solve one synchronous joint velocity/B interval. A zero stored timestep is a
validated no-op only for an originally zero input. Nonzero inputs must survive
conversion. Positive timesteps require resolved half steps and strict clock
progress. Newton/GMRES/probe/backtracking budgets are fixed; failure never retries
with a different timestep and leaves the public state and prior report intact.
"""
function step!(st::ImplicitHybridStepper{D,T}, species::AbstractVector, dt::Real) where {D,T}
    _acquire_implicit_hybrid!(st)
    try
        _implicit_hybrid_binding(st, species)
        _implicit_hybrid_boundary(st)
        h = _validated_step_dt(T, dt, 1; min_NB = 1, name = "implicit hybrid step!")
        iszero(h) && !iszero(dt) &&
            throw(ArgumentError("implicit hybrid nonzero dt must survive conversion to grid precision"))
        iszero(h) && return st
        time, count = _validated_step_clock(st.time[], st.step[], h)
        time > st.time[] || throw(ArgumentError("implicit hybrid clock must advance at stored precision"))
        _implicit_hybrid_stage!(st, species, h)
        _implicit_hybrid_solve!(st)
        report = try
            _implicit_hybrid_endpoint!(st)
            _implicit_hybrid_report(st, time, count)
        catch err
            err isa _ImplicitHybridTrialError || rethrow()
            throw(_implicit_hybrid_failure(st, :endpoint_domain, CapturedException(err, catch_backtrace())))
        end
        prepared, counts = _implicit_hybrid_replacements(st.work.endpoint_species)
        return _implicit_hybrid_commit!(st, prepared, counts, st.work.endpoint,
            true, time, count, h, report)
    finally
        st.busy = false
    end
end

"""
    implicit_hybrid_diagnostics(st, species=st.owners)

Return synchronous integer-time energy/momentum reservoirs, clock and the last
accepted solve report. The report's represented work and continuity defects are
diagnostics of the completed interval, not an exact energy-conservation theorem.
"""
function implicit_hybrid_diagnostics(st::ImplicitHybridStepper, species = st.owners)
    _acquire_implicit_hybrid!(st)
    try
        _implicit_hybrid_binding(st, species)
        _implicit_hybrid_boundary(st)
        return (time = st.time[], step = st.step[], last_dt = st.last_dt[],
            phase = :synchronous_integer, velocity_time_offset = zero(st.time[]),
            energy = energy_budget(species, st.fields.B, st.fields.n, st.model.closure, st.g),
            momentum = momentum_budget(species, st.fields.B, st.g),
            last_solve = st.last_report)
    finally
        st.busy = false
    end
end
