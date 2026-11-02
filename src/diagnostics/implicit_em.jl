# Actual represented constraints and complete signed work, not a conservation
# assertion inferred from the nonlinear tolerance.
function _implicit_em_finite_scalar(value::Rational{BigInt})
    denominator(value) > 0 || throw(ArgumentError("implicit EM exact report scalar must be finite"))
    return nothing
end
function _implicit_em_finite_scalar(value::Real)
    isfinite(value) || throw(ArgumentError("implicit EM report scalar must be finite"))
    return nothing
end
_implicit_em_finite_scalar(values::Tuple) = foreach(_implicit_em_finite_scalar, values)
function _implicit_em_finite_scalar(report::ImplicitEMSolveReport)
    for i = 1:fieldcount(typeof(report))
        _implicit_em_finite_scalar(getfield(report, i))
    end
    return nothing
end
function _implicit_em_validate_report(st::ImplicitEMPICStepper{D,T}, r, time, step, h) where {D,T}
    r isa ImplicitEMSolveReport{T} && r.step == step && isequal(r.time, time) && isequal(r.dt, h) ||
        throw(ArgumentError("implicit EM report must match actual completed clock/history"))
    _implicit_em_finite_scalar(r)
    c = st.controls
    Q = _exact_rational
    r.initial_residual >= 0 &&
        r.tolerance >= 0 &&
        r.tolerance == _from_exact(T, Q(c.atol) + Q(c.rtol) * Q(r.initial_residual)) &&
        all(
            x -> 0 <= x <= r.tolerance,
            (r.velocity_residual, r.electric_residual, r.magnetic_residual),
        ) &&
        all(>=(one(T)), r.scales) ||
        throw(ArgumentError("implicit EM report residual/scale policy differs"))
    0 <= r.newton_iterations <= c.max_newton &&
        0 <= r.linear_iterations <= BigInt(r.newton_iterations) * c.max_linear &&
        0 <= r.backtracks <= BigInt(r.newton_iterations) * c.max_backtrack &&
        r.residual_evaluations >= 2 ||
        throw(ArgumentError("implicit EM report counters differ from fixed solver policy"))
    # Exact physical gates below compare before conversion. Their report values
    # are rounded once; monotonic rounding permits that same rounded boundary.
    tolerance = _from_exact(T, Q(st.constraint_tolerances[1]) + Q(st.constraint_tolerances[2]))
    all(x -> 0 <= x <= tolerance, (r.continuity, r.gauss, r.magnetic_divergence)) &&
        all(>(0), r.constraint_scales) ||
        throw(ArgumentError("implicit EM report constraint metrics are invalid"))
    r.defects == (
        r.delta_kinetic - r.particle_work,
        r.particle_work - r.raw_work,
        r.raw_work - r.current_work,
        r.delta_field + r.current_work,
    ) && r.energy_change == r.delta_kinetic + r.delta_field == sum(r.defects) ||
        throw(ArgumentError("implicit EM complete signed work ledger is inconsistent"))
    return nothing
end
function _implicit_em_set_constraint_scales!(st::ImplicitEMPICStepper, fields, species)
    w, Q = st.work, _exact_rational
    divergence!(w.divergence[1], fields.E, st.g)
    all(isfinite, w.divergence[1]) || throw(_ImplicitEMTrialError(:constraint_range))
    charges = _pic_collection_charges(st.g, zero(st.c), species)
    volume = prod(Q, st.g.dx)
    w.charge_scale =
        max(one(Rational{BigInt}), charges.scale / volume, maximum(x -> abs(Q(x)), w.divergence[1]))
    magnetic = maximum(a -> maximum(x -> abs(Q(x)), a), fields.B)
    w.divergence_scale = max(one(Rational{BigInt}), magnetic / minimum(Q, st.g.dx))
    return nothing
end
function _implicit_em_constraint_values(
    st::ImplicitEMPICStepper{D,T},
    fields,
    rho0,
    h;
    continuity::Bool = true,
) where {D,T}
    w, Q = st.work, _exact_rational
    all(a -> all(isfinite, a), _implicit_em_field_arrays(fields)) ||
        throw(_ImplicitEMTrialError(:field_range))
    divergence!(w.divergence[1], fields.E, st.g)
    divergence!(w.divergence[2], fields.B, st.g)
    continuity && divergence!(w.divergence[3], fields.J, st.g)
    all(a -> all(isfinite, a), continuity ? w.divergence : (w.divergence[1], w.divergence[2])) ||
        throw(_ImplicitEMTrialError(:constraint_range))
    C = G = M = zero(Rational{BigInt})
    for I in eachindex(rho0)
        continuity &&
            (C = max(C, abs(Q(fields.rho_np1[I]) - Q(rho0[I]) + Q(h) * Q(w.divergence[3][I]))))
        G = max(G, abs(Q(w.divergence[1][I]) - Q(fields.rho_np1[I])))
        M = max(M, abs(Q(w.divergence[2][I])))
    end
    exact = (C / w.charge_scale, G / w.charge_scale, M / w.divergence_scale)
    tolerance = Q(st.constraint_tolerances[1]) + Q(st.constraint_tolerances[2])
    all(x -> x <= tolerance, exact) || throw(_ImplicitEMTrialError(:constraint_tolerance))
    values = map(x -> _from_exact(T, x), exact)
    all(isfinite, values) || throw(_ImplicitEMTrialError(:constraint_range))
    return values
end
function _implicit_em_initial_constraints!(st::ImplicitEMPICStepper, fields, species)
    _implicit_em_set_constraint_scales!(st, fields, species)
    return _implicit_em_constraint_values(st, fields, fields.rho_n, zero(st.c); continuity = false)
end
function _implicit_em_ready_constraints!(st::ImplicitEMPICStepper{D,T}) where {D,T}
    w, Q = st.work, _exact_rational
    _implicit_em_set_constraint_scales!(st, w.old, w.old_species)
    # Cached rho_np1 is the completed endpoint. Re-deposit only into private
    # old.rho_n, so the next interval uses actual particle positions as source.
    _implicit_em_charge!(w.old.rho_n, st, w.old_species)
    mismatch = zero(Rational{BigInt})
    for I in eachindex(w.old.rho_n)
        mismatch = max(mismatch, abs(Q(w.old.rho_n[I]) - Q(w.old.rho_np1[I])))
    end
    mismatch / w.charge_scale <= Q(st.constraint_tolerances[1]) + Q(st.constraint_tolerances[2]) ||
        throw(ArgumentError("implicit EM cached charge differs from current particles"))
    # Use the newly deposited boundary source when checking Gauss.
    copyto!(w.midpoint.rho_np1, w.old.rho_n)
    foreach(copyto!, w.midpoint.E, w.old.E)
    foreach(copyto!, w.midpoint.B, w.old.B)
    return _implicit_em_constraint_values(st, w.midpoint, w.old.rho_n, zero(T); continuity = false)
end
function _implicit_em_dot_exact(a, b)
    Q = _exact_rational
    value = zero(Rational{BigInt})
    for c = 1:3, I in eachindex(a[c])
        value += Q(a[c][I]) * Q(b[c][I])
    end
    return value
end
function _implicit_em_report(st::ImplicitEMPICStepper{D,T}, time::T, step::Int) where {D,T}
    w, Q = st.work, _exact_rational
    constraints = _implicit_em_constraint_values(st, w.endpoint, w.old.rho_n, w.h)
    DK = _kinetic_energy_rational(w.endpoint_species) - _kinetic_energy_rational(w.old_species)
    DU =
        _field_energy_rational(w.endpoint.E, w.endpoint.B, Q(st.c)^2, st.g.dx) -
        _field_energy_rational(w.old.E, w.old.B, Q(st.c)^2, st.g.dx)
    Wp = zero(Rational{BigInt})
    offset = 0
    for (old, endpoint) in zip(w.old_species, w.endpoint_species)
        for p in eachindex(old.weight), c = 1:3
            Wp +=
                Q(w.h) * Q(old.q) * Q(old.weight[p]) * (Q(old.v[c][p]) + Q(endpoint.v[c][p])) / 2 *
                Q(w.Ep[c][offset+p])
        end
        offset += nparticles(old)
    end
    volume = prod(Q, st.g.dx)
    Wraw = Q(w.h) * volume * _implicit_em_dot_exact(w.midpoint.E, w.endpoint.Jraw)
    W = Q(w.h) * volume * _implicit_em_dot_exact(w.midpoint.E, w.endpoint.J)
    p0, p1 = _momentum_rational(w.old_species), _momentum_rational(w.endpoint_species)
    r = ImplicitEMSolveReport(
        step,
        time,
        w.h,
        w.initial_residual,
        w.tolerance,
        w.velocity_residual,
        w.electric_residual,
        w.magnetic_residual,
        (w.velocity_scale, w.electric_scale, w.magnetic_scale),
        w.newton_iterations,
        w.linear_iterations,
        w.residual_evaluations,
        w.backtracks,
        constraints...,
        (w.charge_scale, w.divergence_scale),
        DK,
        DU,
        Wp,
        Wraw,
        W,
        (DK - Wp, Wp - Wraw, Wraw - W, DU + W),
        DK + DU,
        ntuple(c -> p1[c] - p0[c], 3),
    )
    _implicit_em_validate_report(st, r, time, step, w.h)
    return r
end

"""
    implicit_em_diagnostics(st, species=st.owners)

Synchronous physical K/U/H and momentum, actual represented constraints and the
last detached solve/work report. The signed energy defect uses complete exact
endpoint reservoirs; continuity-correction work is distinct from solver error.
Rounded diagnostic overflow remains visible. This is not an energy theorem.
"""
function implicit_em_diagnostics(st::ImplicitEMPICStepper{D,T}, species = st.owners) where {D,T}
    _acquire_implicit_em!(st)
    try
        _implicit_em_binding(st, species)
        _implicit_em_boundary(st)
        _implicit_em_stage!(st, species, zero(T))
        constraints = _implicit_em_ready_constraints!(st)
        K = _kinetic_energy_rational(species)
        U = _field_energy_rational(st.fields.E, st.fields.B, _exact_rational(st.c)^2, st.g.dx)
        return (
            kinetic = _from_exact(T, K),
            field = _from_exact(T, U),
            total = _from_exact(T, K + U),
            momentum = map(p -> _from_exact(T, p), _momentum_rational(species)),
            particle_charge = _pic_collection_charges(st.g, zero(T), species).particle,
            gauss = constraints[2],
            magnetic_divergence = constraints[3],
            velocity_time_offset = zero(T),
            magnetic_time_offset = zero(T),
            current_time_offset = st.step[] == 0 ? zero(T) : -st.last_dt[] / 2,
            last_report = deepcopy(st.last_report),
        )
    finally
        st.busy = false
    end
end
