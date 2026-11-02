# All exact report/source values are known numeric records, not opaque metadata.
function _pyee_finite(x::Rational{BigInt})
    denominator(x) > 0 || throw(ArgumentError("compatible Yee exact scalar must be finite"))
    return nothing
end
function _pyee_finite(x::Real)
    isfinite(x) || throw(ArgumentError("compatible Yee scalar must be finite"))
    return nothing
end
_pyee_finite(x::Tuple) = foreach(_pyee_finite, x)
_pyee_finite(x::Vector{Rational{BigInt}}) = foreach(_pyee_finite, x)
function _pyee_finite(report::CompatibleYeeSolveReport)
    for i = 1:fieldcount(typeof(report))
        _pyee_finite(getfield(report, i))
    end
    return nothing
end
function _pyee_validate_report(st::ImplicitYeeEMPICStepper{T}, r, time, step, h) where {T}
    r isa CompatibleYeeSolveReport{T} &&
        r.step == step &&
        isequal(r.time, time) &&
        isequal(r.dt, h) ||
        throw(ArgumentError("compatible Yee report must match completed clock/history"))
    _pyee_finite(r)
    Q, c = _exact_rational, st.controls
    r.initial_residual >= 0 &&
        r.tolerance >= 0 &&
        r.tolerance == _from_exact(T, Q(c.atol) + Q(c.rtol) * Q(r.initial_residual)) &&
        all(>=(one(T)), r.scales) &&
        all(x -> 0 <= x <= Q(r.tolerance), r.exact_residuals) &&
        r.residuals == map(x -> _from_exact(T, x), r.exact_residuals) ||
        throw(ArgumentError("compatible Yee report residual/scale policy differs"))
    0 <= r.newton_iterations <= c.max_newton &&
        0 <= r.linear_iterations <= BigInt(r.newton_iterations) * c.max_linear &&
        0 <= r.backtracks <= BigInt(r.newton_iterations) * c.max_backtrack &&
        r.residual_evaluations >= 2 ||
        throw(ArgumentError("compatible Yee report counters differ from fixed policy"))
    tolerance = Q(st.constraint_tolerances[1]) + Q(st.constraint_tolerances[2])
    r.charge_scale >= 1 &&
        0 <= r.continuity <= tolerance &&
        0 <= r.gauss <= tolerance &&
        r.kinematic_bound >= 0 &&
        r.current_rounding_bound >= 0 &&
        r.continuity <= r.kinematic_bound + r.current_rounding_bound &&
        0 <= r.path_error < Q(st.g.L[1]) ||
        throw(ArgumentError("compatible Yee report constraint/bound values differ"))
    r.defects == (
            r.delta_kinetic - r.particle_work,
            r.particle_work - r.exact_current_work,
            r.exact_current_work - r.current_work,
            r.delta_field + r.current_work,
        ) &&
        r.defects[1] == r.particle_residual_work &&
        iszero(r.defects[2]) &&
        iszero(r.curl_adjoint_work) &&
        r.defects[4] == r.field_residual_work + r.curl_adjoint_work &&
        r.energy_change == r.delta_kinetic + r.delta_field == sum(r.defects) ||
        throw(ArgumentError("compatible Yee complete signed work identities differ"))
    return nothing
end
function _pyee_constraints(
    st,
    fields,
    current,
    charge0,
    charge1,
    h;
    continuity::Bool = true,
    scale = st.work.charge_scale,
)
    _pyee_field_values(fields, st.g)
    foreach(_pyee_finite, (charge0, charge1))
    iszero(sum(charge0)) && iszero(sum(charge1)) ||
        throw(ArgumentError("compatible Yee source snapshots must be exactly neutral"))
    all(a -> all(isfinite, a), current) ||
        throw(ArgumentError("compatible Yee current must be finite"))
    Q, dx = _exact_rational, _exact_rational(st.g.dx[1])
    C = G = zero(Rational{BigInt})
    for i in eachindex(charge1)
        incidence = dx * _yee_difference_exact(fields.E[1], CartesianIndex(i), 1, st.g, false)
        G = max(G, abs(incidence - charge1[i]))
        if continuity
            flux = dx * _yee_difference_exact(current[1], CartesianIndex(i), 1, st.g, false)
            C = max(C, abs(charge1[i] - charge0[i] + Q(h) * flux))
        end
    end
    values = (C / scale, G / scale)
    tolerance = Q(st.constraint_tolerances[1]) + Q(st.constraint_tolerances[2])
    all(x -> x <= tolerance, values) || throw(_CompatibleYeeTrialError(:constraint_tolerance))
    return values
end
function _pyee_ready!(st)
    isequal(st.work.charge0, st.charge_np1) ||
        throw(ArgumentError("compatible Yee cached endpoint charge differs from particles"))
    return _pyee_constraints(
        st,
        st.work.old_fields,
        st.J,
        st.work.charge0,
        st.work.charge0,
        zero(st.c);
        continuity = false,
    )
end
function _pyee_report(st::ImplicitYeeEMPICStepper{T}, time::T, step::Int) where {T}
    w, Q = st.work, _exact_rational
    all(x -> x <= Q(w.tolerance), w.exact_residuals) ||
        throw(_CompatibleYeeTrialError(:exact_final_residual))
    continuity, gauss = _pyee_constraints(st, w.endpoint, w.current, w.charge0, w.charge1, w.h)
    DK = _kinetic_energy_rational(w.endpoint_species) - _kinetic_energy_rational(w.old_species)
    DU =
        _field_energy_rational(w.endpoint.E, w.endpoint.B, Q(st.c)^2, st.g.dx) -
        _field_energy_rational(w.old_fields.E, w.old_fields.B, Q(st.c)^2, st.g.dx)
    Wp = particle_residual = zero(Rational{BigInt})
    offset = 0
    for old in w.old_species
        for p in eachindex(old.weight), c = 1:3
            index = offset + p
            U = w.average_velocity[c][index]
            Wp += w.tau * Q(old.q) * Q(old.weight[p]) * U * w.Ep[c][index]
            particle_residual += Q(old.m) * Q(old.weight[p]) * U * w.force_residual[c][index]
        end
        offset += nparticles(old)
    end
    dx = Q(st.g.dx[1])
    Wstar = w.tau * dx * _implicit_em_dot_exact(w.mean_E, w.exact_current)
    W = w.tau * dx * _implicit_em_dot_exact(w.mean_E, w.current)
    field_residual =
        dx * (
            _implicit_em_dot_exact(w.mean_E, w.exact_electric_residual) +
            Q(st.c)^2 * _implicit_em_dot_exact(w.mean_B, w.exact_magnetic_residual)
        )
    curl_work = zero(Rational{BigInt})
    for c = 1:3, i in eachindex(w.mean_E[c])
        I = CartesianIndex(i)
        curl_work +=
            w.tau *
            dx *
            Q(st.c)^2 *
            (
                w.mean_E[c][i] * _yee_curl_exact(w.mean_B, I, c, st.g, false) -
                w.mean_B[c][i] * _yee_curl_exact(w.mean_E, I, c, st.g, true)
            )
    end
    current_bound = zero(Rational{BigInt})
    for i in eachindex(w.charge1)
        previous = i == 1 ? st.g.n[1] : i - 1
        bound =
            w.tau * (
                abs(Q(w.current[1][i]) - w.exact_current[1][i]) +
                abs(Q(w.current[1][previous]) - w.exact_current[1][previous])
            )
        defect =
            w.charge1[i] - w.charge0[i] + w.tau * (Q(w.current[1][i]) - Q(w.current[1][previous]))
        abs(defect) <= w.kinematic_bound[i] + bound ||
            error("compatible Yee exact continuity bound identity failed")
        current_bound = max(current_bound, bound)
    end
    p0, p1 = _momentum_rational(w.old_species), _momentum_rational(w.endpoint_species)
    report = CompatibleYeeSolveReport(
        step,
        time,
        w.h,
        w.initial_residual,
        w.tolerance,
        (w.velocity_scale, w.electric_scale, w.magnetic_scale),
        (w.velocity_residual, w.electric_residual, w.magnetic_residual),
        w.exact_residuals,
        w.newton_iterations,
        w.linear_iterations,
        w.residual_evaluations,
        w.backtracks,
        w.charge_scale,
        continuity,
        gauss,
        maximum(w.kinematic_bound) / w.charge_scale,
        current_bound / w.charge_scale,
        maximum(abs, w.path_error; init = zero(Rational{BigInt})),
        DK,
        DU,
        Wp,
        Wstar,
        W,
        particle_residual,
        field_residual,
        curl_work,
        (DK - Wp, Wp - Wstar, Wstar - W, DU + W),
        DK + DU,
        ntuple(c -> p1[c] - p0[c], 3),
    )
    _pyee_validate_report(st, report, time, step, w.h)
    return report
end

"""
    implicit_yee_em_diagnostics(st, species=st.owners)

Synchronous physical K/U/H, particle momentum, exact extensive source/Gauss
metrics and detached last solve/work report. True displayed overflow is Inf;
complete signed work uses original-Q endpoints, not Inf-minus-Inf.
"""
function implicit_yee_em_diagnostics(st::ImplicitYeeEMPICStepper{T}, species = st.owners) where {T}
    _acquire_pyee!(st)
    try
        _pyee_binding(st, species)
        _pyee_boundary(st)
        _pyee_stage!(st, species, zero(T))
        _, gauss = _pyee_ready!(st)
        K = _kinetic_energy_rational(species)
        U = _field_energy_rational(st.fields.E, st.fields.B, _exact_rational(st.c)^2, st.g.dx)
        return (
            kinetic = _from_exact(T, K),
            field = _from_exact(T, U),
            total = _from_exact(T, K + U),
            momentum = map(x -> _from_exact(T, x), _momentum_rational(species)),
            charge = sum(st.charge_np1),
            gauss = gauss,
            velocity_time_offset = zero(T),
            magnetic_time_offset = zero(T),
            current_time_offset = st.step[] == 0 ? zero(T) : -st.last_dt[] / 2,
            last_report = deepcopy(st.last_report),
        )
    finally
        st.busy = false
    end
end
