# Full-time periodic Yee/cloud midpoint system. The bounded nonlinear algorithm
# is shared; this file owns this family's physics, storage and publication.
struct _CompatibleYeeTrialError <: Exception
    reason::Symbol
end
Base.showerror(io::IO, e::_CompatibleYeeTrialError) =
    print(io, "compatible Yee trial outside admitted ", e.reason, " domain")
struct CompatibleYeeSolveReport{T}
    step::Int
    time::T
    dt::T
    initial_residual::T
    tolerance::T
    scales::NTuple{3,T}
    residuals::NTuple{3,T}
    exact_residuals::NTuple{3,Rational{BigInt}}
    newton_iterations::Int
    linear_iterations::Int
    residual_evaluations::Int
    backtracks::Int
    charge_scale::Rational{BigInt}
    continuity::Rational{BigInt}
    gauss::Rational{BigInt}
    kinematic_bound::Rational{BigInt}
    current_rounding_bound::Rational{BigInt}
    path_error::Rational{BigInt}
    delta_kinetic::Rational{BigInt}
    delta_field::Rational{BigInt}
    particle_work::Rational{BigInt}
    exact_current_work::Rational{BigInt}
    current_work::Rational{BigInt}
    particle_residual_work::Rational{BigInt}
    field_residual_work::Rational{BigInt}
    curl_adjoint_work::Rational{BigInt}
    defects::NTuple{4,Rational{BigInt}}
    energy_change::Rational{BigInt}
    momentum_change::NTuple{3,Rational{BigInt}}
end
mutable struct _CompatibleYeeWork{T}
    old_fields::YeeFields{1,T}
    endpoint::YeeFields{1,T}
    current::NTuple{3,Vector{T}}
    charge0::Vector{Rational{BigInt}}
    charge1::Vector{Rational{BigInt}}
    old_species::Vector{_ElectrostaticHostParticle{1,T}}
    endpoint_species::Vector{_ElectrostaticHostParticle{1,T}}
    mean_E::NTuple{3,Vector{Rational{BigInt}}}
    mean_B::NTuple{3,Vector{Rational{BigInt}}}
    exact_current::NTuple{3,Vector{Rational{BigInt}}}
    average_velocity::NTuple{3,Vector{Rational{BigInt}}}
    Ep::NTuple{3,Vector{Rational{BigInt}}}
    Bp::NTuple{3,Vector{Rational{BigInt}}}
    force_residual::NTuple{3,Vector{Rational{BigInt}}}
    exact_electric_residual::NTuple{3,Vector{Rational{BigInt}}}
    exact_magnetic_residual::NTuple{3,Vector{Rational{BigInt}}}
    lift::Vector{Rational{BigInt}}
    path_error::Vector{Rational{BigInt}}
    kinematic_bound::Vector{Rational{BigInt}}
    factors::Vector{Rational{BigInt}}
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
    tau::Rational{BigInt}
    velocity_scale::T
    electric_scale::T
    magnetic_scale::T
    charge_scale::Rational{BigInt}
    exact_residuals::NTuple{3,Rational{BigInt}}
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
mutable struct ImplicitYeeEMPICStepper{T,S,K}
    const g::YeeGrid{1,T}
    const shape::S
    const c::T
    const fields::YeeFields{1,T}
    const J::NTuple{3,Vector{T}}
    const charge_n::Vector{Rational{BigInt}}
    const charge_np1::Vector{Rational{BigInt}}
    const owners::Vector{_ElectrostaticHostParticle{1,T}}
    const _species_scalars::Vector{Tuple{T,T}}
    const _population_counts::Vector{Int}
    const source_bindings::Vector{_ImplicitEMBinding{T}}
    const controls::_ImplicitHybridControls{T}
    const constraint_tolerances::Tuple{T,T}
    const _configuration::K
    const work::_CompatibleYeeWork{T}
    const initialized::Base.RefValue{Bool}
    const time::Base.RefValue{T}
    const step::Base.RefValue{Int}
    const last_dt::Base.RefValue{T}
    last_report::Union{Nothing,CompatibleYeeSolveReport{T}}
    busy::Bool
end
_pyee_fields(f::YeeFields) = (f.E..., f.B...)
function _pyee_geometry(g::YeeGrid{1,T}, shape::YeeTriangularCloud) where {T}
    _implicit_em_mode(T)
    g.periodic[1] || throw(ArgumentError("compatible implicit Yee requires a periodic grid"))
    L, h, a = _exact_rational(g.L[1]), _exact_rational(g.dx[1]), _exact_rational(shape.half_width)
    L == g.n[1] * h && 0 < a <= L / 2 || throw(
        ArgumentError("compatible implicit Yee requires exact stored tiling and width <= L/2"),
    )
    return nothing
end
function _pyee_configuration(g, shape, c, scalars, controls, constraints)
    return (
        T = typeof(c),
        ncell = g.n,
        L = g.L,
        dx = g.dx,
        periodic = g.periodic,
        width_type = typeof(shape.half_width),
        half_width = shape.half_width,
        c = c,
        species_scalars = Tuple(scalars),
        mobility = :all_mobile,
        background = :none,
        geometry = :periodic_three_image_exact_tiling_v1,
        current = :orbit_adjoint_original_velocity_v1,
        force = :original_Q_lorentz_midpoint_v1,
        path = :straight_less_than_one_period_v1,
        phase = :synchronous_integer,
        norm = :fixed_three_block_infinity_v1,
        solver = :newton_forward_jv_gmres_identity_armijo_v1,
        nonlinear_rtol = controls.rtol,
        nonlinear_atol = controls.atol,
        max_newton = controls.max_newton,
        max_linear = controls.max_linear,
        memory = controls.memory,
        max_backtrack = controls.max_backtrack,
        max_jv_shrink = controls.max_jv_shrink,
        constraint_tolerances = constraints,
    )
end
_pyee_dimension(n, particles) =
    Base.checked_add(Base.checked_mul(3, particles), Base.checked_mul(5, n))
function _pyee_work(g::YeeGrid{1,T}, species, controls) where {T}
    n = g.n[1]
    total = foldl((n, ps) -> Base.checked_add(n, nparticles(ps)), species; init = 0)
    dimension = _pyee_dimension(n, total)
    qvec(k) = zeros(Rational{BigInt}, k)
    qfield() = ntuple(_ -> qvec(n), 3)
    qparticle() = ntuple(_ -> qvec(total), 3)
    vector() = zeros(T, dimension)
    return _CompatibleYeeWork(
        YeeFields(g),
        YeeFields(g),
        ntuple(_ -> zeros(T, n), 3),
        qvec(n),
        qvec(n),
        _ElectrostaticHostParticle{1,T}[],
        _ElectrostaticHostParticle{1,T}[],
        qfield(),
        qfield(),
        qfield(),
        qparticle(),
        qparticle(),
        qparticle(),
        qparticle(),
        qfield(),
        qfield(),
        qvec(total),
        qvec(total),
        qvec(n),
        qvec(length(species)),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        vector(),
        _implicit_hybrid_solver(T, dimension, controls.memory),
        zero(T),
        zero(Rational{BigInt}),
        one(T),
        one(T),
        one(T),
        one(Rational{BigInt}),
        ntuple(_ -> zero(Rational{BigInt}), 3),
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
function _pyee_reserve!(st::ImplicitYeeEMPICStepper{T}, species) where {T}
    w = st.work
    _pic_prepare_species_stages!(w.old_species, species)
    _pic_prepare_species_stages!(w.endpoint_species, species)
    total = foldl((n, ps) -> Base.checked_add(n, nparticles(ps)), species; init = 0)
    for a in (w.average_velocity..., w.Ep..., w.Bp..., w.force_residual..., w.lift, w.path_error)
        resize!(a, total)
    end
    dimension = _pyee_dimension(st.g.n[1], total)
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
        resize!(a, dimension)
    end
    length(w.solver.x) == dimension ||
        (w.solver = _implicit_hybrid_solver(T, dimension, st.controls.memory))
    return nothing
end
function _pyee_owned_arrays(st)
    w = st.work
    return (
        _pyee_fields(st.fields)...,
        st.J...,
        st.charge_n,
        st.charge_np1,
        st.owners,
        st._species_scalars,
        st._population_counts,
        st.source_bindings,
        _implicit_em_binding_arrays(st.source_bindings)...,
        _implicit_hybrid_particle_arrays(st.owners)...,
        _pyee_fields(w.old_fields)...,
        _pyee_fields(w.endpoint)...,
        w.current...,
        w.charge0,
        w.charge1,
        w.old_species,
        w.endpoint_species,
        _implicit_hybrid_particle_arrays(w.old_species)...,
        _implicit_hybrid_particle_arrays(w.endpoint_species)...,
        w.mean_E...,
        w.mean_B...,
        w.exact_current...,
        w.average_velocity...,
        w.Ep...,
        w.Bp...,
        w.force_residual...,
        w.exact_electric_residual...,
        w.exact_magnetic_residual...,
        w.lift,
        w.path_error,
        w.kinematic_bound,
        w.factors,
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
function _pyee_storage(st)
    n, w, ns = st.g.n[1], st.work, length(st.owners)
    all(
        a -> axes(a, 1) == Base.OneTo(ns),
        (st._species_scalars, st._population_counts, st.source_bindings, w.factors),
    ) || throw(DimensionMismatch("compatible Yee binding/cache row count differs"))
    for f in (st.fields, w.old_fields, w.endpoint)
        _yee_check_shape(f.E, st.g, true)
        _yee_check_shape(f.B, st.g, false)
    end
    all(
        a -> axes(a, 1) == Base.OneTo(n),
        (
            st.J...,
            st.charge_n,
            st.charge_np1,
            w.current...,
            w.charge0,
            w.charge1,
            w.mean_E...,
            w.mean_B...,
            w.exact_current...,
            w.exact_electric_residual...,
            w.exact_magnetic_residual...,
            w.kinematic_bound,
        ),
    ) || throw(DimensionMismatch("compatible Yee grid work differs"))
    for stages in (w.old_species, w.endpoint_species)
        length(stages) == ns ||
            throw(DimensionMismatch("compatible Yee private species count differs"))
        foreach(_check_particle_storage_axes, stages)
        foreach(
            ps -> _require_particle_ordinary_stores(ps, "compatible Yee private particles"),
            stages,
        )
    end
    total = foldl((n, ps) -> Base.checked_add(n, nparticles(ps)), w.old_species; init = 0)
    dimension = _pyee_dimension(n, total)
    all(
        a -> length(a) == total,
        (w.average_velocity..., w.Ep..., w.Bp..., w.force_residual..., w.lift, w.path_error),
    ) && all(
        a -> length(a) == dimension,
        (
            w.z,
            w.residual,
            w.probe,
            w.probe_residual,
            w.trial,
            w.trial_residual,
            w.rhs,
            w.zero,
            w.action,
        ),
    ) || throw(DimensionMismatch("compatible Yee nonlinear work differs"))
    _validate_gmres_workspace_layout(w.solver, dimension)
    all(a -> length(a) == dimension, _gmres_workspace_vectors(w.solver)) ||
        throw(DimensionMismatch("compatible Yee GMRES work differs"))
    arrays = _pyee_owned_arrays(st)
    for i in eachindex(arrays), j = 1:(i-1)
        _particle_arrays_mightalias(arrays[i], arrays[j]) &&
            throw(ArgumentError("compatible Yee public/private/source roles must be independent"))
    end
    refs = (st.initialized, st.time, st.step, st.last_dt)
    for i in eachindex(refs), j = 1:(i-1)
        refs[i] === refs[j] &&
            throw(ArgumentError("compatible Yee clock references must be independent"))
    end
    return nothing
end
function _pyee_particles(st, species)
    total = zero(Rational{BigInt})
    for ps in species
        _check_particle_storage_axes(ps)
        _require_particle_ordinary_stores(ps, "compatible implicit Yee")
        _validated_boris_factor(ps)
        _require_finite_particle_velocities(ps, "compatible implicit Yee")
        for p in eachindex(ps.weight)
            isfinite(ps.x[1][p]) && 0 <= _exact_rational(ps.x[1][p]) < _exact_rational(st.g.L[1]) ||
                throw(ArgumentError("compatible Yee positions must be finite in [0,L)"))
            isfinite(ps.weight[p]) && _exact_rational(ps.weight[p]) >= 0 ||
                throw(ArgumentError("compatible Yee weights must be finite nonnegative"))
            total += _exact_rational(ps.q) * _exact_rational(ps.weight[p])
        end
    end
    iszero(total) ||
        throw(ArgumentError("compatible periodic Yee requires exact physical neutrality"))
    return nothing
end
function _pyee_binding(
    st::ImplicitYeeEMPICStepper{T},
    species;
    values = true,
    counts = true,
    sources = true,
) where {T}
    _pyee_geometry(st.g, st.shape)
    _require_finite_positive_real("compatible Yee c", st.c, T)
    isequal(
        st._configuration,
        _pyee_configuration(
            st.g,
            st.shape,
            st.c,
            st._species_scalars,
            st.controls,
            st.constraint_tolerances,
        ),
    ) || throw(ArgumentError("compatible Yee fixed configuration changed"))
    species isa AbstractVector && length(species) == length(st.owners) ||
        throw(ArgumentError("compatible Yee requires its ordered species"))
    _pyee_storage(st)
    for s in eachindex(st.owners)
        ps = species[s]
        ps === st.owners[s] ||
            throw(ArgumentError("compatible Yee population identity/order changed"))
        ps isa _ElectrostaticHostParticle{1,T} ||
            throw(ArgumentError("compatible Yee requires ordinary Vector populations"))
        _check_particle_storage_axes(ps)
        _require_particle_ordinary_stores(ps, "compatible Yee destination")
        values &&
            !isequal((ps.q, ps.m), st._species_scalars[s]) &&
            throw(ArgumentError("compatible Yee original q/m changed"))
        counts &&
            nparticles(ps) != st._population_counts[s] &&
            throw(ArgumentError("compatible Yee count changed outside init/restore"))
        if sources
            b = st.source_bindings[s]
            isequal(ps.weight, b.weight) && isequal(ps.id, b.id) && isequal(ps.tag, b.tag) ||
                throw(ArgumentError("compatible Yee frozen weights/provenance changed"))
        end
    end
    values && _pyee_particles(st, species)
    return nothing
end
function _pyee_field_values(fields, g)
    _yee_validate_fields(fields, g)
    first_value = _exact_rational(first(fields.B[1]))
    all(v -> _exact_rational(v) == first_value, fields.B[1]) ||
        throw(ArgumentError("periodic one-dimensional Bx must be exactly constant"))
    return nothing
end
function _pyee_phase(st::ImplicitYeeEMPICStepper{T}, initialized, time, step, h, report) where {T}
    initialized isa Bool &&
        time isa T &&
        step isa Int &&
        h isa T &&
        isfinite(time) &&
        time >= 0 &&
        step >= 0 &&
        isfinite(h) &&
        _stored_nonnegative(h) || throw(ArgumentError("compatible Yee clock/history is invalid"))
    if step == 0
        _particle_stored_zero(time) && _particle_stored_zero(h) && report === nothing ||
            throw(ArgumentError("unstepped compatible Yee requires zero history and no report"))
    else
        initialized && time > 0 && _stored_positive(h) && h <= time ||
            throw(ArgumentError("advanced compatible Yee requires 0 < last dt <= time"))
        _require_resolved_dt("compatible Yee history", h)
        _pyee_validate_report(st, report, time, step, h)
    end
    return nothing
end
function _pyee_boundary(st)
    st.initialized[] || throw(ArgumentError("initialize compatible Yee before ready use"))
    _pyee_phase(st, true, st.time[], st.step[], st.last_dt[], st.last_report)
    _pyee_field_values(st.fields, st.g)
    all(a -> all(isfinite, a), st.J) ||
        throw(ArgumentError("compatible Yee current must be finite"))
    foreach(_pyee_finite, (st.charge_n, st.charge_np1))
    if st.step[] == 0
        all(a -> all(_particle_stored_zero, a), st.J) && isequal(st.charge_n, st.charge_np1) ||
            throw(ArgumentError("unstepped compatible Yee requires zero current/equal charges"))
    else
        actual = _pyee_constraints(
            st,
            st.fields,
            st.J,
            st.charge_n,
            st.charge_np1,
            st.last_dt[];
            scale = st.last_report.charge_scale,
        )
        actual == (st.last_report.continuity, st.last_report.gauss) ||
            throw(ArgumentError("compatible Yee interval history differs from its report"))
    end
    return nothing
end
function _acquire_pyee!(st)
    st.busy && throw(ArgumentError("compatible Yee instance is already in use"))
    st.busy = true
    return nothing
end
function _pyee_failure(st, reason, cause = nothing)
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
_implicit_pic_residual!(out, st::ImplicitYeeEMPICStepper, z) = _pyee_residual!(out, st, z)
_implicit_pic_failure(st::ImplicitYeeEMPICStepper, reason, cause = nothing) =
    _pyee_failure(st, reason, cause)
_implicit_pic_trial(::ImplicitYeeEMPICStepper, reason) = _CompatibleYeeTrialError(reason)
_implicit_pic_trial_error(::ImplicitYeeEMPICStepper, err) = err isa _CompatibleYeeTrialError
_implicit_pic_scale(st::ImplicitYeeEMPICStepper, i, nv) =
    i <= nv ? st.work.velocity_scale :
    i <= nv + 3st.g.n[1] ? st.work.electric_scale : st.work.magnetic_scale
function _implicit_pic_set_metrics!(st::ImplicitYeeEMPICStepper, metrics)
    st.work.velocity_residual, st.work.electric_residual, st.work.magnetic_residual = metrics
    return max(metrics...)
end
Base.eltype(::_ImplicitHybridJacobian{<:ImplicitYeeEMPICStepper{T}}) where {T} = T
function _pyee_residual!(out::Vector{T}, st::ImplicitYeeEMPICStepper{T}, z::Vector{T}) where {T}
    w, Q, g = st.work, _exact_rational, st.g
    w.residual_evaluations = Base.checked_add(w.residual_evaluations, 1)
    w.velocity_residual = w.electric_residual = w.magnetic_residual = T(Inf)
    all(isfinite, z) || throw(_CompatibleYeeTrialError(:iterate_range))
    L, tau = Q(g.L[1]), w.tau
    i, offset = 0, 0
    for (old, endpoint) in zip(w.old_species, w.endpoint_species)
        for c = 1:3, p in eachindex(old.weight)
            i += 1
            endpoint.v[c][p] = z[i]
            w.average_velocity[c][offset+p] = (Q(old.v[c][p]) + Q(z[i])) / 2
        end
        for p in eachindex(old.weight)
            index = offset + p
            x = Q(old.x[1][p])
            drift = tau * w.average_velocity[1][index]
            abs(drift) < L || throw(_CompatibleYeeTrialError(:transport))
            ideal = _coordinate_eval(Val(true), x, w.average_velocity[1][index], tau)
            wrapped = _periodic_coordinate(ideal, zero(L), L, L)
            winding = Int((ideal - wrapped) / L)
            rounded = _from_exact(T, wrapped)
            rounded == g.L[1] && (winding += 1)
            endpoint.x[1][p] = _periodic_coordinate(rounded, zero(T), g.L[1], g.L[1])
            lifted = Q(endpoint.x[1][p]) + winding * L
            abs(lifted - x) < L || throw(_CompatibleYeeTrialError(:represented_transport))
            w.lift[index] = lifted
            w.path_error[index] = lifted - ideal
        end
        offset += nparticles(old)
    end
    for c = 1:3, index in eachindex(w.endpoint.E[c])
        i += 1
        w.endpoint.E[c][index] = z[i]
        w.mean_E[c][index] = (Q(w.old_fields.E[c][index]) + Q(z[i])) / 2
    end
    copyto!(w.endpoint.B[1], w.old_fields.B[1])
    for index in eachindex(w.endpoint.B[1])
        w.mean_B[1][index] = Q(w.old_fields.B[1][index])
        w.exact_magnetic_residual[1][index] = zero(Rational{BigInt})
    end
    for c = 2:3, index in eachindex(w.endpoint.B[c])
        i += 1
        w.endpoint.B[c][index] = z[i]
        w.mean_B[c][index] = (Q(w.old_fields.B[c][index]) + Q(z[i])) / 2
    end
    _pyee_charge!(w.charge1, w.endpoint_species, g, st.shape)
    mean_ex = sum(Q, w.endpoint.E[1]) / g.n[1]
    _pyee_gauss_ex!(st, w.charge1, mean_ex, _CompatibleYeeTrialError(:field_range))
    for index in eachindex(w.mean_E[1])
        w.mean_E[1][index] = (Q(w.old_fields.E[1][index]) + Q(w.endpoint.E[1][index])) / 2
    end
    _pyee_couple!(w, w.old_species, g, st.shape)
    for c = 1:3, index in eachindex(w.current[c])
        value = _from_exact(T, w.exact_current[c][index])
        isfinite(value) || throw(_CompatibleYeeTrialError(:current_range))
        w.current[c][index] = value
    end
    i, offset = 0, 0
    exact_v = zero(Rational{BigInt})
    for (s, old) in enumerate(w.old_species)
        n = nparticles(old)
        for p in eachindex(old.weight)
            index = offset + p
            velocity = ntuple(c -> w.average_velocity[c][index], 3)
            electric = ntuple(c -> w.Ep[c][index], 3)
            magnetic = ntuple(c -> w.Bp[c][index], 3)
            _, parts = _euler_kick_parts(
                velocity...,
                electric...,
                magnetic...,
                tau * w.factors[s],
                Val(:predict),
            )
            for c = 1:3
                target = i + (c - 1) * n + p
                exact = Q(z[target]) - Q(old.v[c][p]) - parts[c+3]
                w.force_residual[c][index] = exact
                normalized = exact / Q(w.velocity_scale)
                value = _from_exact(T, normalized)
                isfinite(value) || throw(_CompatibleYeeTrialError(:residual_range))
                out[target] = value
                exact_v = max(exact_v, abs(normalized))
            end
        end
        i += 3n
        offset += n
    end
    exact_e = exact_b = zero(Rational{BigInt})
    for c = 1:3, index in eachindex(w.endpoint.E[c])
        i += 1
        curl = _yee_curl_exact(w.mean_B, CartesianIndex(index), c, g, false)
        exact =
            Q(z[i]) - Q(w.old_fields.E[c][index]) -
            tau * (Q(st.c)^2 * curl - Q(w.current[c][index]))
        w.exact_electric_residual[c][index] =
            Q(w.endpoint.E[c][index]) - Q(w.old_fields.E[c][index]) -
            tau * (Q(st.c)^2 * curl - Q(w.current[c][index]))
        normalized = exact / Q(w.electric_scale)
        value = _from_exact(T, normalized)
        isfinite(value) || throw(_CompatibleYeeTrialError(:residual_range))
        out[i] = value
        exact_e = max(exact_e, abs(normalized))
    end
    for c = 2:3, index in eachindex(w.endpoint.B[c])
        i += 1
        curl = _yee_curl_exact(w.mean_E, CartesianIndex(index), c, g, true)
        exact = Q(w.endpoint.B[c][index]) - Q(w.old_fields.B[c][index]) + tau * curl
        w.exact_magnetic_residual[c][index] = exact
        normalized = exact / Q(w.magnetic_scale)
        value = _from_exact(T, normalized)
        isfinite(value) || throw(_CompatibleYeeTrialError(:residual_range))
        out[i] = value
        exact_b = max(exact_b, abs(normalized))
    end
    w.exact_residuals = (exact_v, exact_e, exact_b)
    return map(x -> _from_exact(T, x), w.exact_residuals)
end
function _pyee_charge_scale(st, fields, species, charge)
    Q = _exact_rational
    magnitude = sum(
        (abs(Q(ps.q) * Q(weight)) for ps in species for weight in ps.weight);
        init = zero(Rational{BigInt}),
    )
    incidence = maximum(
        i -> abs(Q(fields.E[1][i]) - Q(fields.E[1][i == 1 ? st.g.n[1] : i - 1])),
        eachindex(fields.E[1]),
    )
    return max(one(Rational{BigInt}), magnitude, maximum(abs, charge), incidence)
end
function _pyee_stage!(st::ImplicitYeeEMPICStepper{T}, species, h::T) where {T}
    _pyee_reserve!(st, species)
    w, Q = st.work, _exact_rational
    foreach(copyto!, _pyee_fields(w.old_fields), _pyee_fields(st.fields))
    foreach(copyto!, _pyee_fields(w.endpoint), _pyee_fields(st.fields))
    w.h, w.tau = h, Q(h)
    for s in eachindex(species)
        w.factors[s] = _boris_exact_factor(_validated_boris_factor(species[s]))
    end
    w.velocity_scale = max(
        one(T),
        maximum((maximum(abs, a; init = zero(T)) for ps in species for a in ps.v); init = zero(T)),
    )
    w.electric_scale = max(one(T), maximum(a -> maximum(abs, a), st.fields.E))
    w.magnetic_scale = max(one(T), maximum(a -> maximum(abs, a), st.fields.B))
    w.initial_residual = w.tolerance = zero(T)
    w.velocity_residual = w.electric_residual = w.magnetic_residual = T(Inf)
    w.exact_residuals = ntuple(_ -> zero(Rational{BigInt}), 3)
    w.newton_iterations = w.linear_iterations = w.residual_evaluations = w.backtracks = 0
    _pyee_charge!(w.charge0, w.old_species, st.g, st.shape)
    w.charge_scale = _pyee_charge_scale(st, st.fields, species, w.charge0)
    i = 0
    for ps in species, c = 1:3, value in ps.v[c]
        i += 1
        w.z[i] = value
    end
    for component in (st.fields.E..., st.fields.B[2], st.fields.B[3]), value in component
        i += 1
        w.z[i] = value
    end
    return nothing
end
function _pyee_gauss_ex!(
    st::ImplicitYeeEMPICStepper{T},
    charge::Vector{Rational{BigInt}},
    mean::Rational{BigInt},
    failure,
) where {T}
    w, n = st.work, st.g.n[1]
    prefix = w.exact_electric_residual[1]
    prefix[1] = zero(Rational{BigInt})
    for i = 2:n
        prefix[i] = prefix[i-1] + charge[i]
    end
    prefix_mean = sum(prefix) / n
    for i = 1:n
        value = _from_exact(T, mean + prefix[i] - prefix_mean)
        isfinite(value) || throw(failure)
        w.endpoint.E[1][i] = value
    end
    return nothing
end
function _pyee_commit!(
    st::ImplicitYeeEMPICStepper{T},
    prepared::Vector{_ElectrostaticHostParticle{1,T}},
    counts::Vector{Int},
    bindings::Vector{_ImplicitEMBinding{T}},
    fields::YeeFields{1,T},
    current::NTuple{3,Vector{T}},
    charge0::Vector{Rational{BigInt}},
    charge1::Vector{Rational{BigInt}},
    time::T,
    count::Int,
    h::T,
    report::Union{Nothing,CompatibleYeeSolveReport{T}},
) where {T}
    for s in eachindex(st.owners)
        _publish_particle_replacement!(st.owners[s], prepared[s])
        st.owners[s].q, st.owners[s].m = st._species_scalars[s]
    end
    foreach(copyto!, _pyee_fields(st.fields), _pyee_fields(fields))
    foreach(copyto!, st.J, current)
    copyto!(st.charge_n, charge0)
    copyto!(st.charge_np1, charge1)
    copyto!(st._population_counts, counts)
    bindings === st.source_bindings || copyto!(st.source_bindings, bindings)
    st.last_report = report
    st.time[] = time
    st.step[] = count
    st.last_dt[] = h
    st.initialized[] = true
    return st
end

"""
    ImplicitYeeEMPICStepper(g::YeeGrid{1}, species; shape::YeeTriangularCloud, c=1, kwargs...)

Nonrelativistic synchronous periodic CPU1D3V cloud midpoint PIC, with directly
paired orbit current/electric gather and complete source/work residuals. Requires
exact stored tiling, width<=L/2, exact neutrality and less-than-one-period trial
transport. Constant Bx may be nonzero. No exact floating-point conservation,
nonlinear stability or speed claim is implied. See docs/IMPLICIT_YEE_EM.md.
"""
function ImplicitYeeEMPICStepper(
    g::YeeGrid{1,T},
    species;
    shape::YeeTriangularCloud,
    c = 1,
    nonlinear_rtol = sqrt(eps(T)),
    nonlinear_atol = 64eps(T),
    max_newton = 20,
    max_linear = 200,
    memory = 20,
    max_backtrack = 12,
    max_jv_shrink = 8,
    constraint_rtol = 256eps(T),
    constraint_atol = zero(T),
) where {T}
    _pyee_geometry(g, shape)
    species isa Vector{<:ParticleSet{1,T}} &&
        all(ps -> ps isa _ElectrostaticHostParticle{1,T}, species) || throw(
        ArgumentError(
            "compatible Yee requires a typed ordinary Vector of Vector-backed populations",
        ),
    )
    speed = _require_finite_positive_real("compatible Yee c", c, T)
    controls = _implicit_pic_controls(
        T,
        nonlinear_rtol,
        nonlinear_atol,
        max_newton,
        max_linear,
        memory,
        max_backtrack,
        max_jv_shrink,
        "compatible implicit Yee",
    )
    constraints = _implicit_em_constraint_controls(T, constraint_rtol, constraint_atol)
    owners = _ElectrostaticHostParticle{1,T}[ps for ps in species]
    scalars = Tuple{T,T}[(ps.q, ps.m) for ps in owners]
    bindings = _implicit_em_bindings(owners)
    n = g.n[1]
    config = _pyee_configuration(g, shape, speed, scalars, controls, constraints)
    st = ImplicitYeeEMPICStepper(
        g,
        shape,
        speed,
        YeeFields(g),
        ntuple(_ -> zeros(T, n), 3),
        zeros(Rational{BigInt}, n),
        zeros(Rational{BigInt}, n),
        owners,
        scalars,
        Int[nparticles(ps) for ps in owners],
        bindings,
        controls,
        constraints,
        config,
        _pyee_work(g, owners, controls),
        Ref(false),
        Ref(zero(T)),
        Ref(0),
        Ref(zero(T)),
        nothing,
        false,
    )
    _pyee_reserve!(st, owners)
    _pyee_binding(st, species)
    return st
end
function init!(st::ImplicitYeeEMPICStepper{T}, species::AbstractVector) where {T}
    _acquire_pyee!(st)
    try
        _pyee_binding(st, species; counts = false, sources = false)
        _pyee_field_values(st.fields, st.g)
        _pyee_stage!(st, species, zero(T))
        w = st.work
        _pyee_gauss_ex!(
            st,
            w.charge0,
            sum(_exact_rational, w.old_fields.E[1]) / st.g.n[1],
            ArgumentError("compatible Yee initial Ex does not round finite"),
        )
        copyto!(w.charge1, w.charge0)
        foreach(a -> fill!(a, zero(T)), w.current)
        _pyee_constraints(
            st,
            w.endpoint,
            w.current,
            w.charge0,
            w.charge1,
            zero(T);
            continuity = false,
        )
        prepared, counts = _implicit_em_replacements(w.endpoint_species)
        bindings = _implicit_em_bindings(w.endpoint_species)
        return _pyee_commit!(
            st,
            prepared,
            counts,
            bindings,
            w.endpoint,
            w.current,
            w.charge0,
            w.charge1,
            zero(T),
            0,
            zero(T),
            nothing,
        )
    finally
        st.busy = false
    end
end
function step!(st::ImplicitYeeEMPICStepper{T}, species::AbstractVector, dt::Real) where {T}
    _acquire_pyee!(st)
    try
        _pyee_binding(st, species)
        _pyee_boundary(st)
        h = _validated_step_dt(T, dt, 1; min_NB = 1, name = "compatible implicit Yee step!")
        iszero(h) &&
            !iszero(dt) &&
            throw(ArgumentError("compatible Yee nonzero dt must survive conversion"))
        _pyee_stage!(st, species, h)
        _pyee_ready!(st)
        iszero(h) && return st
        time, count = _validated_step_clock(st.time[], st.step[], h, 1)
        time > st.time[] || throw(ArgumentError("compatible Yee clock must advance finitely"))
        _implicit_pic_solve!(st, T)
        report = try
            _pyee_report(st, time, count)
        catch err
            err isa _CompatibleYeeTrialError || rethrow()
            throw(_pyee_failure(st, :endpoint_domain, CapturedException(err, catch_backtrace())))
        end
        prepared, counts = _implicit_em_replacements(st.work.endpoint_species)
        return _pyee_commit!(
            st,
            prepared,
            counts,
            st.source_bindings,
            st.work.endpoint,
            st.work.current,
            st.work.charge0,
            st.work.charge1,
            time,
            count,
            h,
            report,
        )
    finally
        st.busy = false
    end
end
