# Exactly two synchronous PIC families share the existing bounded nonlinear
# algorithm. Physics residuals, scales, trial-domain errors and reports stay in
# their owning models; this is not a public solver callback interface.
_implicit_pic_residual!(out, st::ImplicitHybridStepper, z) = _implicit_hybrid_residual!(out, st, z)
_implicit_pic_failure(st::ImplicitHybridStepper, reason, cause = nothing) =
    _implicit_hybrid_failure(st, reason, cause)
_implicit_pic_trial(::ImplicitHybridStepper, reason) = _ImplicitHybridTrialError(reason)
_implicit_pic_trial_error(::ImplicitHybridStepper, err) = err isa _ImplicitHybridTrialError
_implicit_pic_scale(st::ImplicitHybridStepper, i, nv) =
    i <= nv ? st.work.velocity_scale : st.work.magnetic_scale
_implicit_pic_set_metrics!(st::ImplicitHybridStepper, metrics) =
    _implicit_hybrid_set_metrics!(st, metrics)

function _implicit_pic_jv!(out::Vector{T}, A::_ImplicitHybridJacobian, p::Vector{T}) where {T}
    st, w = A.state, A.state.work
    all(isfinite, p) || throw(_implicit_pic_failure(st, :nonfinite_jv_direction))
    all(iszero, p) && return fill!(out, zero(T))
    pnorm = BLAS.nrm2(p)
    isfinite(pnorm) && pnorm > zero(T) || throw(_implicit_pic_failure(st, :jv_scale))
    nv = 3length(w.Ep[1])
    zscale = one(T)
    for i in eachindex(w.z)
        scale = _implicit_pic_scale(st, i, nv)
        zscale = max(zscale, abs(w.z[i] / scale))
    end
    Q = _exact_rational
    epsilon = _from_exact(T, Q(sqrt(eps(T))) * Q(zscale) / Q(pnorm))
    isfinite(epsilon) && epsilon > zero(T) || throw(_implicit_pic_failure(st, :jv_scale))
    cause = nothing
    for attempt = 0:st.controls.max_jv_shrink
        changed = false
        try
            for i in eachindex(w.z)
                scale = _implicit_pic_scale(st, i, nv)
                value = _from_exact(T, Q(w.z[i]) + Q(epsilon) * Q(scale) * Q(p[i]))
                isfinite(value) || throw(_implicit_pic_trial(st, :probe_range))
                changed |= value != w.z[i]
                w.probe[i] = value
            end
            changed || throw(_implicit_pic_failure(st, :unresolved_jv))
            _implicit_pic_residual!(w.probe_residual, st, w.probe)
            for i in eachindex(out)
                value = _from_exact(T, (Q(w.probe_residual[i]) - Q(w.residual[i])) / Q(epsilon))
                isfinite(value) || throw(_implicit_pic_failure(st, :jv_action_range))
                out[i] = value
            end
            return out
        catch err
            _implicit_pic_trial_error(st, err) || rethrow()
            cause = CapturedException(err, catch_backtrace())
        end
        attempt == st.controls.max_jv_shrink && break
        epsilon /= 2
        epsilon > zero(T) || throw(_implicit_pic_failure(st, :unresolved_jv, cause))
    end
    throw(_implicit_pic_failure(st, :jv_domain, cause))
end

function _implicit_pic_solve!(
    st::Union{ImplicitHybridStepper,ImplicitEMPICStepper},
    ::Type{T},
) where {T}
    w, controls = st.work, st.controls
    initial = try
        _implicit_pic_residual!(w.residual, st, w.z)
    catch err
        _implicit_pic_trial_error(st, err) || rethrow()
        throw(_implicit_pic_failure(st, :initial_domain, CapturedException(err, catch_backtrace())))
    end
    rho = _implicit_pic_set_metrics!(st, initial)
    w.initial_residual = rho
    Q = _exact_rational
    w.tolerance = _from_exact(T, Q(controls.atol) + Q(controls.rtol) * Q(rho))
    isfinite(w.tolerance) && (iszero(rho) || w.tolerance > zero(T)) ||
        throw(_implicit_pic_failure(st, :tolerance_range))
    A = _ImplicitHybridJacobian(st)
    for iteration = 0:controls.max_newton
        rho <= w.tolerance && break
        iteration == controls.max_newton && throw(_implicit_pic_failure(st, :nonlinear_limit))
        for i in eachindex(w.rhs)
            w.rhs[i] = -w.residual[i]
        end
        rhsnorm = BLAS.nrm2(w.rhs)
        target = rhsnorm / 2
        estimated_target = target / 2
        isfinite(rhsnorm) && isfinite(target) && estimated_target > zero(T) ||
            throw(_implicit_pic_failure(st, :linear_tolerance_range))
        fill!(w.zero, zero(T))
        # Check the defined action at the explicit zero correction before
        # touching restarted history. No previous correction is retained.
        mul!(w.action, A, w.zero)
        for i in eachindex(w.action)
            w.action[i] = _from_exact(T, Q(w.action[i]) + Q(w.residual[i]))
        end
        zero_residual = BLAS.nrm2(w.action)
        isfinite(zero_residual) || throw(_implicit_pic_failure(st, :linear_initial_residual))
        delta = if zero_residual <= target
            w.zero
        else
            Krylov.gmres!(
                w.solver,
                A,
                w.rhs,
                w.zero;
                restart = true,
                reorthogonalization = true,
                atol = estimated_target,
                rtol = zero(T),
                itmax = controls.max_linear,
                callback = _ -> _implicit_hybrid_linear_tick!(st),
            )
            Krylov.solution(w.solver)
        end
        all(isfinite, delta) || throw(_implicit_pic_failure(st, :nonfinite_linear_solution))
        mul!(w.action, A, delta)
        for i in eachindex(w.action)
            w.action[i] = _from_exact(T, Q(w.action[i]) + Q(w.residual[i]))
        end
        linear_residual = BLAS.nrm2(w.action)
        isfinite(linear_residual) && linear_residual <= target ||
            throw(_implicit_pic_failure(st, :linear_residual))
        accepted = false
        cause = nothing
        nv = 3length(w.Ep[1])
        alpha = one(T)
        for backtrack = 0:controls.max_backtrack
            changed = false
            try
                for i in eachindex(w.z)
                    scale = _implicit_pic_scale(st, i, nv)
                    value = _from_exact(T, Q(w.z[i]) + Q(alpha) * Q(scale) * Q(delta[i]))
                    isfinite(value) || throw(_implicit_pic_trial(st, :line_range))
                    changed |= value != w.z[i]
                    w.trial[i] = value
                end
                changed || throw(_implicit_pic_failure(st, :unresolved_line_search))
                metrics = _implicit_pic_residual!(w.trial_residual, st, w.trial)
                trial_rho = max(metrics...)
                armijo = _from_exact(T, (1 - Q(alpha) / 4) * Q(rho))
                if trial_rho <= w.tolerance || trial_rho <= armijo
                    copyto!(w.z, w.trial)
                    copyto!(w.residual, w.trial_residual)
                    rho = _implicit_pic_set_metrics!(st, metrics)
                    w.newton_iterations = Base.checked_add(w.newton_iterations, 1)
                    accepted = true
                    break
                end
            catch err
                _implicit_pic_trial_error(st, err) || rethrow()
                cause = CapturedException(err, catch_backtrace())
            end
            w.backtracks = Base.checked_add(w.backtracks, 1)
            backtrack == controls.max_backtrack && break
            alpha /= 2
            alpha > zero(T) || throw(_implicit_pic_failure(st, :unresolved_line_search, cause))
        end
        accepted || throw(_implicit_pic_failure(st, :backtrack_limit, cause))
    end
    # Reconstruct every coupled quantity from immutable old snapshots and z.
    # GMRES's action/probe and the accepted trial cache are not acceptance data.
    final = try
        _implicit_pic_residual!(w.residual, st, w.z)
    catch err
        _implicit_pic_trial_error(st, err) || rethrow()
        throw(_implicit_pic_failure(st, :final_domain, CapturedException(err, catch_backtrace())))
    end
    _implicit_pic_set_metrics!(st, final) <= w.tolerance ||
        throw(_implicit_pic_failure(st, :final_residual))
    return nothing
end


function _implicit_pic_controls(
    ::Type{T},
    rtol,
    atol,
    max_newton,
    max_linear,
    memory,
    max_backtrack,
    max_jv_shrink,
    context,
) where {T}
    r = _require_finite_nonnegative_real("nonlinear_rtol", rtol, T)
    a = _require_finite_nonnegative_real("nonlinear_atol", atol, T)
    r < one(T) && (r > zero(T) || a > zero(T)) ||
        throw(ArgumentError("$context requires rtol < 1 and a positive tolerance"))
    (iszero(rtol) || !iszero(r)) && (iszero(atol) || !iszero(a)) ||
        throw(ArgumentError("$context positive tolerances must survive conversion"))
    return _ImplicitHybridControls(
        r,
        a,
        _implicit_hybrid_nonnegative_limit("max_newton", max_newton),
        _require_positive_intlike("max_linear", max_linear),
        _require_positive_intlike("memory", memory),
        _implicit_hybrid_nonnegative_limit("max_backtrack", max_backtrack),
        _implicit_hybrid_nonnegative_limit("max_jv_shrink", max_jv_shrink),
    )
end

_implicit_pic_mean_value(a::T, b::T) where {T} =
    _from_exact(T, (_exact_rational(a) + _exact_rational(b)) / 2)
