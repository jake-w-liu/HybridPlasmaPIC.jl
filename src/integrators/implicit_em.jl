# Synchronous implicit Vlasov–Maxwell with actual endpoint-charge coupling.
struct _ImplicitEMFields{D,T}
    E::NTuple{3,Array{T,D}}
    B::NTuple{3,Array{T,D}}
    J::NTuple{3,Array{T,D}}
    Jraw::NTuple{3,Array{T,D}}
    rho_n::Array{T,D}
    rho_np1::Array{T,D}
end
function _ImplicitEMFields(g::FourierGrid{D,T}) where {D,T}
    vector() = ntuple(_ -> zeros(T, g.n), 3)
    return _ImplicitEMFields(vector(), vector(), vector(), vector(), zeros(T, g.n), zeros(T, g.n))
end
_implicit_em_field_state(f::_ImplicitEMFields) =
    (E = f.E, B = f.B, J = f.J, Jraw = f.Jraw, rho_n = f.rho_n, rho_np1 = f.rho_np1)
_implicit_em_field_arrays(f::_ImplicitEMFields) =
    (f.E..., f.B..., f.J..., f.Jraw..., f.rho_n, f.rho_np1)

struct _ImplicitEMBinding{T}
    weight::Vector{T}
    id::Vector{UInt64}
    tag::Vector{UInt32}
end
struct ImplicitEMSolveReport{T}
    step::Int
    time::T
    dt::T
    initial_residual::T
    tolerance::T
    velocity_residual::T
    electric_residual::T
    magnetic_residual::T
    scales::NTuple{3,T}
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
    continuity::T
    gauss::T
    magnetic_divergence::T
    constraint_scales::NTuple{2,Rational{BigInt}}
    delta_kinetic::Rational{BigInt}
    delta_field::Rational{BigInt}
    particle_work::Rational{BigInt}
    raw_work::Rational{BigInt}
    current_work::Rational{BigInt}
    defects::NTuple{4,Rational{BigInt}}
    energy_change::Rational{BigInt}
    momentum_change::NTuple{3,Rational{BigInt}}
end
struct _ImplicitEMTrialError <: Exception
    reason::Symbol
end
Base.showerror(io::IO, e::_ImplicitEMTrialError) =
    print(io, "implicit EM trial outside admitted ", e.reason, " domain")
struct ImplicitEMConvergenceError{T} <: Exception
    reason::Symbol
    residuals::NTuple{3,T}
    tolerance::T
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
    cause::Union{Nothing,CapturedException}
end
function Base.showerror(io::IO, e::ImplicitEMConvergenceError)
    print(
        io,
        "implicit EM ",
        e.reason,
        "; velocity/E/B residuals=",
        e.residuals,
        ", tolerance=",
        e.tolerance,
        ", Newton iterations=",
        e.newton_iterations,
    )
    e.cause === nothing || (print(io, "; "); showerror(io, e.cause.ex))
end
mutable struct _ImplicitEMWorkspace{D,T}
    old::_ImplicitEMFields{D,T}
    midpoint::_ImplicitEMFields{D,T}
    endpoint::_ImplicitEMFields{D,T}
    old_species::Vector{_ElectrostaticHostParticle{D,T}}
    midpoint_species::Vector{_ElectrostaticHostParticle{D,T}}
    endpoint_species::Vector{_ElectrostaticHostParticle{D,T}}
    mobile::Vector{Bool}
    density::Array{T,D}
    curlE::NTuple{3,Array{T,D}}
    curlB::NTuple{3,Array{T,D}}
    divergence::NTuple{3,Array{T,D}}
    current_work::Vector{T}
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
    electric_scale::T
    magnetic_scale::T
    charge_scale::Rational{BigInt}
    divergence_scale::Rational{BigInt}
    initial_residual::T
    tolerance::T
    velocity_residual::T
    electric_residual::T
    magnetic_residual::T
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
end
mutable struct ImplicitEMPICStepper{D,T,G,K}
    const g::G
    const shape::CIC
    const c::T
    const fields::_ImplicitEMFields{D,T}
    const owners::Vector{_ElectrostaticHostParticle{D,T}}
    const _species_scalars::Vector{Tuple{T,T}}
    const _population_counts::Vector{Int}
    const source_bindings::Vector{_ImplicitEMBinding{T}}
    const controls::_ImplicitHybridControls{T}
    const constraint_tolerances::Tuple{T,T}
    const _configuration::K
    const work::_ImplicitEMWorkspace{D,T}
    const initialized::Base.RefValue{Bool}
    const time::Base.RefValue{T}
    const step::Base.RefValue{Int}
    const last_dt::Base.RefValue{T}
    last_report::Union{Nothing,ImplicitEMSolveReport{T}}
    busy::Bool
end

function _implicit_em_mode(::Type{T}) where {T}
    T in (Float32, Float64) && Base.rounding(T) === RoundNearest && !Base.get_zero_subnormals() ||
        throw(
            ArgumentError(
                "implicit EM requires CPU Float32/64 nearest rounding and gradual underflow",
            ),
        )
    return nothing
end
function _implicit_em_constraint_controls(::Type{T}, rtol, atol) where {T}
    r = _require_finite_nonnegative_real("constraint_rtol", rtol, T)
    a = _require_finite_nonnegative_real("constraint_atol", atol, T)
    r < one(T) && (r > 0 || a > 0) && (iszero(rtol) || r > 0) && (iszero(atol) || a > 0) || throw(
        ArgumentError(
            "implicit EM constraint tolerances must survive conversion and define a positive gate",
        ),
    )
    return r, a
end
function _implicit_em_configuration(g, c, scalars, controls, constraints)
    return (
        family = :implicit_empic_midpoint,
        D = length(g.n),
        T = typeof(c),
        ncell = g.n,
        L = g.L,
        dx = g.dx,
        c = c,
        species_scalars = Tuple(scalars),
        shape = :cic,
        mobility = :all_mobile,
        background = :zero_exact_neutral,
        grid_policy = :odd_axes,
        source_policy = :closed_frozen_weight_id_tag,
        phase = :synchronous_integer,
        residual = :boris_map_maxwell_midpoint_v1,
        force = :finite_represented_gather_original_factor_v1,
        quadrature = :lifted_single_represented_midpoint_v1,
        current = :raw_midpoint_spectral_continuity_v1,
        derivative = :represented_fourier_odd_v1,
        norm = :fixed_reference_block_infinity_v1,
        solver = :newton_forward_jv_gmres_identity_armijo_v1,
        nonlinear_rtol = controls.rtol,
        nonlinear_atol = controls.atol,
        max_newton = controls.max_newton,
        max_linear = controls.max_linear,
        memory = controls.memory,
        max_backtrack = controls.max_backtrack,
        max_jv_shrink = controls.max_jv_shrink,
        constraint_rtol = constraints[1],
        constraint_atol = constraints[2],
    )
end
_implicit_em_dimension(g, particles) = Base.checked_add(
    Base.checked_mul(3, particles),
    Base.checked_mul(6, foldl(Base.checked_mul, g.n; init = 1)),
)
function _implicit_em_workspace(g::FourierGrid{D,T}, species, controls) where {D,T}
    count = foldl((n, ps) -> Base.checked_add(n, nparticles(ps)), species; init = 0)
    n = _implicit_em_dimension(g, count)
    gridvector() = ntuple(_ -> zeros(T, g.n), 3)
    particlevector() = ntuple(_ -> zeros(T, count), 3)
    vector() = zeros(T, n)
    return _ImplicitEMWorkspace(
        _ImplicitEMFields(g),
        _ImplicitEMFields(g),
        _ImplicitEMFields(g),
        _ElectrostaticHostParticle{D,T}[],
        _ElectrostaticHostParticle{D,T}[],
        _ElectrostaticHostParticle{D,T}[],
        fill(true, length(species)),
        zeros(T, g.n),
        gridvector(),
        gridvector(),
        gridvector(),
        T[],
        particlevector(),
        particlevector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        _implicit_hybrid_solver(T, n, controls.memory),
        zero(T),
        one(T),
        one(T),
        one(T),
        one(Rational{BigInt}),
        one(Rational{BigInt}),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        0,
        0,
        0,
        0,
    )
end
function _implicit_em_reserve!(st::ImplicitEMPICStepper{D,T}, species) where {D,T}
    w = st.work
    total = foldl((n, ps) -> Base.checked_add(n, nparticles(ps)), species; init = 0)
    n = _implicit_em_dimension(st.g, total)
    for stages in (w.old_species, w.midpoint_species, w.endpoint_species)
        _pic_prepare_species_stages!(stages, species)
    end
    for a in (w.Ep..., w.Bp...)
        resize!(a, total)
    end
    for a in (
        w.z,
        w.residual,
        w.probe,
        w.probe_residual,
        w.trial,
        w.trial_residual,
        w.rhs,
        w.zero,
        w.action,
    )
        resize!(a, n)
    end
    length(w.solver.x) == n || (w.solver = _implicit_hybrid_solver(T, n, st.controls.memory))
    return nothing
end
_implicit_em_binding_arrays(bindings) = Tuple(a for b in bindings for a in (b.weight, b.id, b.tag))
function _implicit_em_work_arrays(w::_ImplicitEMWorkspace)
    return (
        _implicit_em_field_arrays(w.old)...,
        _implicit_em_field_arrays(w.midpoint)...,
        _implicit_em_field_arrays(w.endpoint)...,
        w.old_species,
        w.midpoint_species,
        w.endpoint_species,
        _implicit_hybrid_particle_arrays(w.old_species)...,
        _implicit_hybrid_particle_arrays(w.midpoint_species)...,
        _implicit_hybrid_particle_arrays(w.endpoint_species)...,
        w.mobile,
        w.density,
        w.curlE...,
        w.curlB...,
        w.divergence...,
        w.current_work,
        w.Ep...,
        w.Bp...,
        w.z,
        w.residual,
        w.probe,
        w.probe_residual,
        w.trial,
        w.trial_residual,
        w.rhs,
        w.zero,
        w.action,
        _gmres_workspace_arrays(w.solver)...,
    )
end
function _implicit_em_owned_arrays(st::ImplicitEMPICStepper)
    return (
        _implicit_em_field_arrays(st.fields)...,
        st.owners,
        st._species_scalars,
        st._population_counts,
        st.source_bindings,
        _implicit_em_binding_arrays(st.source_bindings)...,
        _implicit_hybrid_particle_arrays(st.owners)...,
        _implicit_em_work_arrays(st.work)...,
        st.g.cbuf,
        st.g.tbuf,
        st.g.abuf,
        st.g.sbuf,
    )
end
function _implicit_em_storage(st::ImplicitEMPICStepper{D,T}) where {D,T}
    g, w = st.g, st.work
    ns = length(st.owners)
    for a in (st._species_scalars, st._population_counts, st.source_bindings, w.mobile)
        axes(a, 1) == Base.OneTo(ns) ||
            throw(DimensionMismatch("implicit EM binding/cache row count differs"))
    end
    all(w.mobile) || throw(ArgumentError("implicit EM every population is mobile"))
    for f in (st.fields, w.old, w.midpoint, w.endpoint), a in _implicit_em_field_arrays(f)
        size(a) == g.n || throw(DimensionMismatch("implicit EM field dimensions differ"))
    end
    for a in (w.density, w.curlE..., w.curlB..., w.divergence..., g.cbuf, g.tbuf, g.abuf, g.sbuf)
        size(a) == g.n || throw(DimensionMismatch("implicit EM grid work dimensions differ"))
    end
    for d = 1:D, a in (g.midx[d], g.kfull[d], g.ik[d], g.kvec[d])
        axes(a, 1) == Base.OneTo(g.n[d]) ||
            throw(DimensionMismatch("implicit EM Fourier metadata dimensions differ"))
    end
    for stages in (w.old_species, w.midpoint_species, w.endpoint_species)
        length(stages) == ns ||
            throw(DimensionMismatch("implicit EM private species count differs"))
        foreach(_check_particle_storage_axes, stages)
    end
    total = foldl((n, ps) -> Base.checked_add(n, nparticles(ps)), w.old_species; init = 0)
    n = _implicit_em_dimension(g, total)
    length(w.z) == n &&
        all(a -> length(a) == total, (w.Ep..., w.Bp...)) &&
        all(
            a -> length(a) == n,
            (
                w.residual,
                w.probe,
                w.probe_residual,
                w.trial,
                w.trial_residual,
                w.rhs,
                w.zero,
                w.action,
            ),
        ) || throw(DimensionMismatch("implicit EM nonlinear/gather dimensions differ"))
    _validate_gmres_workspace_layout(w.solver, n)
    all(a -> length(a) == n, _gmres_workspace_vectors(w.solver)) ||
        throw(DimensionMismatch("implicit EM GMRES dimensions differ"))
    arrays = _implicit_em_owned_arrays(st)
    metadata = (g.midx..., g.kfull..., g.ik..., g.kvec...)
    for i in eachindex(arrays)
        for j = 1:(i-1)
            _particle_arrays_mightalias(arrays[i], arrays[j]) &&
                throw(ArgumentError("implicit EM public/private/source roles must be independent"))
        end
        for a in metadata
            _particle_arrays_mightalias(arrays[i], a) &&
                throw(ArgumentError("implicit EM writers must not alias Fourier metadata"))
        end
    end
    refs = (st.initialized, st.time, st.step, st.last_dt)
    for i in eachindex(refs), j = 1:(i-1)
        refs[i] === refs[j] &&
            throw(ArgumentError("implicit EM clock references must be independent"))
    end
    return nothing
end
function _implicit_em_particles(st::ImplicitEMPICStepper{D,T}, species) where {D,T}
    for ps in species
        ps isa _ElectrostaticHostParticle{D,T} ||
            throw(ArgumentError("implicit EM requires matching ordinary Vector particles"))
        _check_particle_storage_axes(ps)
        _validated_boris_factor(ps)
        _validate_particle_positions(ps)
        _require_finite_particle_velocities(ps, "implicit EM")
        all(w -> isfinite(w) && _stored_nonnegative(w), ps.weight) ||
            throw(ArgumentError("implicit EM weights must be finite and nonnegative"))
        for d = 1:D
            all(x -> zero(T) <= x < st.g.L[d], ps.x[d]) ||
                throw(ArgumentError("implicit EM positions must lie in the periodic domain"))
        end
    end
    charges = _pic_collection_charges(st.g, zero(T), species)
    iszero(charges.net) ||
        throw(ArgumentError("implicit EM requires exact physical total charge zero"))
    return charges
end
function _implicit_em_binding(
    st::ImplicitEMPICStepper{D,T},
    species;
    values = true,
    counts = values,
    sources = values,
) where {D,T}
    _implicit_em_mode(T)
    all(isodd, st.g.n) ||
        throw(ArgumentError("implicit EM first family requires odd spatial grid lengths"))
    _require_finite_positive_real("implicit EM c", st.c, T)
    c = st.controls
    _implicit_pic_controls(
        T,
        c.rtol,
        c.atol,
        c.max_newton,
        c.max_linear,
        c.memory,
        c.max_backtrack,
        c.max_jv_shrink,
        "implicit EM",
    )
    _implicit_em_constraint_controls(T, st.constraint_tolerances...)
    isequal(
        _implicit_em_configuration(st.g, st.c, st._species_scalars, c, st.constraint_tolerances),
        st._configuration,
    ) || throw(ArgumentError("implicit EM fixed configuration binding changed"))
    _implicit_em_storage(st)
    species isa Vector{<:ParticleSet{D,T}} && length(species) == length(st.owners) ||
        throw(ArgumentError("implicit EM requires its bound ordered ordinary population Vector"))
    for (s, ps) in enumerate(species)
        ps === st.owners[s] && ps isa _ElectrostaticHostParticle{D,T} ||
            throw(ArgumentError("implicit EM population identity/order/storage changed"))
        _require_particle_ordinary_stores(ps, "implicit EM")
        _check_particle_storage_axes(ps)
        q, m = st._species_scalars[s]
        _validated_charge_mass(q, m, T)
        values &&
            !isequal((ps.q, ps.m), st._species_scalars[s]) &&
            throw(ArgumentError("implicit EM original charge/mass changed"))
        counts &&
            nparticles(ps) != st._population_counts[s] &&
            throw(ArgumentError("implicit EM count changed outside initialization/restore"))
        if sources
            binding = st.source_bindings[s]
            isequal(ps.weight, binding.weight) &&
                isequal(ps.id, binding.id) &&
                isequal(ps.tag, binding.tag) ||
                throw(ArgumentError("implicit EM closed weight/provenance binding changed"))
        end
    end
    values && _implicit_em_particles(st, species)
    return nothing
end
function _implicit_em_bindings(species::Vector{_ElectrostaticHostParticle{D,T}}) where {D,T}
    return _ImplicitEMBinding{T}[
        _ImplicitEMBinding(copy(ps.weight), copy(ps.id), copy(ps.tag)) for ps in species
    ]
end
function _implicit_em_phase(
    st::ImplicitEMPICStepper{D,T},
    initialized,
    time,
    step,
    h,
    report,
) where {D,T}
    initialized isa Bool &&
        time isa T &&
        step isa Int &&
        h isa T &&
        isfinite(time) &&
        time >= 0 &&
        step >= 0 &&
        isfinite(h) &&
        _stored_nonnegative(h) ||
        throw(ArgumentError("implicit EM clock/history types or values are invalid"))
    if step == 0
        _particle_stored_zero(time) && _particle_stored_zero(h) && report === nothing ||
            throw(ArgumentError("unstepped implicit EM requires zero clock/history and no report"))
    else
        initialized && time > 0 && _stored_positive(h) && h <= time || throw(
            ArgumentError(
                "advanced implicit EM requires positive initialized phase with last dt <= time",
            ),
        )
        _require_resolved_dt("implicit EM history", h)
        _implicit_em_validate_report(st, report, time, step, h)
    end
    return nothing
end
function _implicit_em_boundary(st::ImplicitEMPICStepper)
    st.initialized[] || throw(ArgumentError("initialize implicit EM before ready use"))
    _implicit_em_phase(st, true, st.time[], st.step[], st.last_dt[], st.last_report)
    all(a -> all(isfinite, a), _implicit_em_field_arrays(st.fields)) ||
        throw(ArgumentError("implicit EM fields must be finite"))
    if st.step[] == 0
        all(a -> all(_particle_stored_zero, a), (st.fields.J..., st.fields.Jraw...)) &&
            isequal(st.fields.rho_n, st.fields.rho_np1) || throw(
            ArgumentError("unstepped implicit EM requires zero current and equal charge snapshots"),
        )
    end
    return nothing
end
function _acquire_implicit_em!(st::ImplicitEMPICStepper)
    st.busy && throw(ArgumentError("implicit EM permits one nonreentrant caller"))
    st.busy = true
    return nothing
end
function _implicit_em_failure(st, reason, cause = nothing)
    w = st.work
    return ImplicitEMConvergenceError(
        reason,
        (w.velocity_residual, w.electric_residual, w.magnetic_residual),
        w.tolerance,
        w.newton_iterations,
        w.linear_iterations,
        w.residual_evaluations,
        w.backtracks,
        cause,
    )
end
_implicit_pic_residual!(out, st::ImplicitEMPICStepper, z) = _implicit_em_residual!(out, st, z)
_implicit_pic_failure(st::ImplicitEMPICStepper, reason, cause = nothing) =
    _implicit_em_failure(st, reason, cause)
_implicit_pic_trial(::ImplicitEMPICStepper, reason) = _ImplicitEMTrialError(reason)
_implicit_pic_trial_error(::ImplicitEMPICStepper, err) = err isa _ImplicitEMTrialError
function _implicit_pic_scale(st::ImplicitEMPICStepper, i, nv)
    return i <= nv ? st.work.velocity_scale :
           i <= nv + 3length(st.fields.rho_n) ? st.work.electric_scale : st.work.magnetic_scale
end
function _implicit_pic_set_metrics!(st::ImplicitEMPICStepper, metrics)
    st.work.velocity_residual, st.work.electric_residual, st.work.magnetic_residual = metrics
    return max(metrics...)
end
Base.eltype(::_ImplicitHybridJacobian{<:ImplicitEMPICStepper{D,T}}) where {D,T} = T
function _implicit_em_mean(a::T, b::T) where {T}
    value = _implicit_pic_mean_value(a, b)
    isfinite(value) || throw(_ImplicitEMTrialError(:affine_range))
    return value
end
function _implicit_em_charge!(rho, st::ImplicitEMPICStepper{D,T}, species) where {D,T}
    return _empic_charge_sum!(
        rho,
        st.work.density,
        species,
        st.g,
        st.shape,
        zero(T),
        zero(Rational{BigInt}),
    )
end
function _implicit_em_residual!(
    out::Vector{T},
    st::ImplicitEMPICStepper{D,T},
    z::Vector{T},
) where {D,T}
    w, g = st.work, st.g
    w.residual_evaluations = Base.checked_add(w.residual_evaluations, 1)
    w.velocity_residual = w.electric_residual = w.magnetic_residual = T(Inf)
    all(isfinite, z) || throw(_ImplicitEMTrialError(:iterate_range))
    i = 0
    for (old, mid, endpoint) in zip(w.old_species, w.midpoint_species, w.endpoint_species)
        for c = 1:3, p in eachindex(old.weight)
            i += 1
            endpoint.v[c][p] = z[i]
            mid.v[c][p] = _implicit_em_mean(old.v[c][p], z[i])
        end
        for d = 1:D, p in eachindex(old.weight)
            half = _coordinate_eval(Val(true), old.x[d][p], mid.v[d][p], w.h / 2)
            full = _coordinate_eval(Val(true), old.x[d][p], mid.v[d][p], w.h)
            isfinite(half) && isfinite(full) || throw(_ImplicitEMTrialError(:lift_range))
            mid.x[d][p] = _periodic_coordinate(half, zero(T), g.L[d], g.L[d])
            endpoint.x[d][p] = _periodic_coordinate(full, zero(T), g.L[d], g.L[d])
        end
    end
    for (old, mid, endpoint) in
        ((w.old.E, w.midpoint.E, w.endpoint.E), (w.old.B, w.midpoint.B, w.endpoint.B))
        for c = 1:3, I in eachindex(old[c])
            i += 1
            endpoint[c][I] = z[i]
            mid[c][I] = _implicit_em_mean(old[c][I], z[i])
        end
    end
    _implicit_em_charge!(w.endpoint.rho_np1, st, w.endpoint_species)
    data = (work = w.current_work, density = w.density)
    _empic_current_sum!(
        w.endpoint.Jraw,
        data,
        w.midpoint_species,
        w.midpoint_species,
        w.mobile,
        g,
        st.shape,
        Val(:stored_positions),
    )
    _empic_correct_current!(w.endpoint.J, w.endpoint.Jraw, w.old.rho_n, w.endpoint.rho_np1, g, w.h)
    all(a -> all(isfinite, a), w.endpoint.J) || throw(_ImplicitEMTrialError(:current_range))
    curl!(w.curlE, w.midpoint.E, g)
    curl!(w.curlB, w.midpoint.B, g)
    all(a -> all(isfinite, a), (w.curlE..., w.curlB...)) ||
        throw(_ImplicitEMTrialError(:curl_range))
    Q = _exact_rational
    i, offset = 0, 0
    vmax = zero(T)
    for (old, mid) in zip(w.old_species, w.midpoint_species)
        n = nparticles(old)
        indices = (offset+1):(offset+n)
        Ep = ntuple(c -> view(w.Ep[c], indices), 3)
        Bp = ntuple(c -> view(w.Bp[c], indices), 3)
        gather_vector!(Ep, w.midpoint.E, mid, g, st.shape)
        gather_vector!(Bp, w.midpoint.B, mid, g, st.shape)
        all(a -> all(isfinite, a), (Ep..., Bp...)) || throw(_ImplicitEMTrialError(:gather_range))
        factor = _validated_boris_factor(old)
        for p in eachindex(old.weight)
            kicked = _kick_eval(
                Val(true),
                old.v[1][p],
                old.v[2][p],
                old.v[3][p],
                Ep[1][p],
                Ep[2][p],
                Ep[3][p],
                Bp[1][p],
                Bp[2][p],
                Bp[3][p],
                factor,
                w.h,
            )
            all(isfinite, kicked) || throw(_ImplicitEMTrialError(:boris_range))
            for c = 1:3
                index = i + (c - 1) * n + p
                value = _from_exact(T, (Q(z[index]) - Q(kicked[c])) / Q(w.velocity_scale))
                isfinite(value) || throw(_ImplicitEMTrialError(:residual_range))
                out[index] = value
                vmax = max(vmax, abs(value))
            end
        end
        i += 3n
        offset += n
    end
    emax, bmax = zero(T), zero(T)
    for c = 1:3, I in eachindex(w.old.E[c])
        i += 1
        value = _from_exact(
            T,
            (
                Q(z[i]) - Q(w.old.E[c][I]) -
                Q(w.h) * (Q(st.c)^2 * Q(w.curlB[c][I]) - Q(w.endpoint.J[c][I]))
            ) / Q(w.electric_scale),
        )
        isfinite(value) || throw(_ImplicitEMTrialError(:residual_range))
        out[i] = value
        emax = max(emax, abs(value))
    end
    for c = 1:3, I in eachindex(w.old.B[c])
        i += 1
        value = _from_exact(
            T,
            (Q(z[i]) - Q(w.old.B[c][I]) + Q(w.h) * Q(w.curlE[c][I])) / Q(w.magnetic_scale),
        )
        isfinite(value) || throw(_ImplicitEMTrialError(:residual_range))
        out[i] = value
        bmax = max(bmax, abs(value))
    end
    return vmax, emax, bmax
end
function _implicit_em_stage!(st::ImplicitEMPICStepper{D,T}, species, h::T) where {D,T}
    _implicit_em_reserve!(st, species)
    w = st.work
    for fields in (w.old, w.midpoint, w.endpoint)
        foreach(copyto!, _implicit_em_field_arrays(fields), _implicit_em_field_arrays(st.fields))
    end
    w.h = h
    w.velocity_scale = max(
        one(T),
        maximum((maximum(abs, a; init = zero(T)) for ps in species for a in ps.v); init = zero(T)),
    )
    w.electric_scale = max(one(T), maximum(a -> maximum(abs, a), w.old.E))
    w.magnetic_scale = max(one(T), maximum(a -> maximum(abs, a), w.old.B))
    w.initial_residual = w.tolerance = zero(T)
    w.velocity_residual = w.electric_residual = w.magnetic_residual = T(Inf)
    w.newton_iterations = w.linear_iterations = w.residual_evaluations = w.backtracks = 0
    i = 0
    for ps in species, c = 1:3, value in ps.v[c]
        i += 1
        w.z[i] = value
    end
    for field in (w.old.E, w.old.B), c = 1:3, value in field[c]
        i += 1
        w.z[i] = value
    end
    return nothing
end
function _implicit_em_replacements(species::Vector{_ElectrostaticHostParticle{D,T}}) where {D,T}
    return (
        _ElectrostaticHostParticle{D,T}[
            _prepare_particle_replacement(ps, collect(eachindex(ps.weight))) for ps in species
        ],
        Int[nparticles(ps) for ps in species],
    )
end
function _implicit_em_commit!(st, prepared, counts, bindings, fields, time, step, h, report)
    for s in eachindex(st.owners)
        _publish_particle_replacement!(st.owners[s], prepared[s])
        st.owners[s].q, st.owners[s].m = st._species_scalars[s]
    end
    foreach(copyto!, _implicit_em_field_arrays(st.fields), _implicit_em_field_arrays(fields))
    copyto!(st._population_counts, counts)
    bindings === st.source_bindings || copyto!(st.source_bindings, bindings)
    st.last_report = report
    st.time[] = time
    st.step[] = step
    st.last_dt[] = h
    st.initialized[] = true
    return st
end

"""
    ImplicitEMPICStepper(g, CIC(), species; c=5, nonlinear_rtol=sqrt(eps(T)),
        nonlinear_atol=64eps(T), max_newton=20, max_linear=200, memory=20,
        max_backtrack=12, max_jv_shrink=8, constraint_rtol=256eps(T), constraint_atol=0)

Distinct nonrelativistic synchronous CPU Float32/64 Vlasov–Maxwell collection on
odd periodic grids. All ordinary Vector populations are mobile, with exact zero
physical total charge and frozen weights/count/provenance between init/restore.
Seed `fields.E/B`, then initialize full-time particles. No exact conservation,
nonlinear stability, native-bit compatibility or speed guarantee is implied.
See docs/IMPLICIT_EM.md. One instance/population/grid permits one caller.
"""
function ImplicitEMPICStepper(
    g::FourierGrid{D,T},
    shape::CIC,
    species;
    c = 5,
    nonlinear_rtol = sqrt(eps(T)),
    nonlinear_atol = 64eps(T),
    max_newton = 20,
    max_linear = 200,
    memory = 20,
    max_backtrack = 12,
    max_jv_shrink = 8,
    constraint_rtol = 256eps(T),
    constraint_atol = zero(T),
) where {D,T}
    _implicit_em_mode(T)
    1 <= D <= 3 && all(isodd, g.n) ||
        throw(ArgumentError("implicit EM requires D1–3 and odd spatial grid lengths"))
    species isa Vector{<:ParticleSet{D,T}} &&
        all(ps -> ps isa _ElectrostaticHostParticle{D,T}, species) || throw(
        ArgumentError(
            "implicit EM requires an ordinary typed Vector of Vector-backed matching populations",
        ),
    )
    speed = _require_finite_positive_real("implicit EM c", c, T)
    controls = _implicit_pic_controls(
        T,
        nonlinear_rtol,
        nonlinear_atol,
        max_newton,
        max_linear,
        memory,
        max_backtrack,
        max_jv_shrink,
        "implicit EM",
    )
    constraints = _implicit_em_constraint_controls(T, constraint_rtol, constraint_atol)
    owners = _ElectrostaticHostParticle{D,T}[ps for ps in species]
    scalars = Tuple{T,T}[(ps.q, ps.m) for ps in owners]
    bindings = _implicit_em_bindings(owners)
    config = _implicit_em_configuration(g, speed, scalars, controls, constraints)
    st = ImplicitEMPICStepper(
        g,
        shape,
        speed,
        _ImplicitEMFields(g),
        owners,
        scalars,
        Int[nparticles(ps) for ps in owners],
        bindings,
        controls,
        constraints,
        config,
        _implicit_em_workspace(g, owners, controls),
        Ref(false),
        Ref(zero(T)),
        Ref(0),
        Ref(zero(T)),
        nothing,
        false,
    )
    _implicit_em_reserve!(st, owners)
    _implicit_em_binding(st, species)
    return st
end
function init!(st::ImplicitEMPICStepper{D,T}, species::AbstractVector) where {D,T}
    _acquire_implicit_em!(st)
    try
        _implicit_em_binding(st, species; counts = false, sources = false)
        all(a -> all(isfinite, a), (st.fields.E..., st.fields.B...)) ||
            throw(ArgumentError("implicit EM initial E/B must be finite"))
        _implicit_em_stage!(st, species, zero(T))
        w = st.work
        foreach(
            copyto!,
            _implicit_em_field_arrays(w.endpoint),
            _implicit_em_field_arrays(st.fields),
        )
        _implicit_em_charge!(w.endpoint.rho_n, st, w.endpoint_species)
        copyto!(w.endpoint.rho_np1, w.endpoint.rho_n)
        _empic_gauss!(w.endpoint.E, w.endpoint.rho_n, st.g)
        foreach(a -> fill!(a, zero(T)), (w.endpoint.J..., w.endpoint.Jraw...))
        _implicit_em_initial_constraints!(st, w.endpoint, w.endpoint_species)
        prepared, counts = _implicit_em_replacements(w.endpoint_species)
        bindings = _implicit_em_bindings(w.endpoint_species)
        return _implicit_em_commit!(
            st,
            prepared,
            counts,
            bindings,
            w.endpoint,
            zero(T),
            0,
            zero(T),
            nothing,
        )
    finally
        st.busy = false
    end
end
function step!(st::ImplicitEMPICStepper{D,T}, species::AbstractVector, dt::Real) where {D,T}
    _acquire_implicit_em!(st)
    try
        _implicit_em_binding(st, species)
        _implicit_em_boundary(st)
        h = _validated_step_dt(T, dt, 1; min_NB = 1, name = "implicit EM step!")
        iszero(h) &&
            !iszero(dt) &&
            throw(ArgumentError("implicit EM nonzero dt must survive conversion"))
        _implicit_em_stage!(st, species, h)
        _implicit_em_ready_constraints!(st)
        iszero(h) && return st
        time, count = _validated_step_clock(st.time[], st.step[], h, 1)
        time > st.time[] ||
            throw(ArgumentError("implicit EM clock must advance at stored precision"))
        _implicit_pic_solve!(st, T)
        # The shared driver has freshly rebuilt all endpoint data/gathers.
        copyto!(st.work.endpoint.rho_n, st.work.old.rho_n)
        report = try
            _implicit_em_report(st, time, count)
        catch err
            err isa _ImplicitEMTrialError || rethrow()
            throw(
                _implicit_em_failure(
                    st,
                    :endpoint_domain,
                    CapturedException(err, catch_backtrace()),
                ),
            )
        end
        prepared, counts = _implicit_em_replacements(st.work.endpoint_species)
        return _implicit_em_commit!(
            st,
            prepared,
            counts,
            st.source_bindings,
            st.work.endpoint,
            time,
            count,
            h,
            report,
        )
    finally
        st.busy = false
    end
end
