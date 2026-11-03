# Specific constrained mixed-field projection; not a shock time integrator.

"A mixed Hall projection failed a true residual, constraint or work gate before publication."
struct MixedHallProjectionError{R} <: Exception
    reason::String
    report::R
end

function Base.showerror(io::IO, err::MixedHallProjectionError)
    print(io, "mixed Hall projection rejected: ", err.reason, "; ", err.report)
end

"""
    MixedHallProjectionWorkspace(op)

Private O(cells + constraint rows) field, face and cold MINRES-QLP storage.
No multiplier, history or workspace is physical checkpoint state. One caller;
public projection outputs must be independent of all operator/workspace roles.
"""
struct MixedHallProjectionWorkspace{D,T,S}
    n::NTuple{D,Int}
    L::NTuple{D,T}
    solver::S
    rhs::Vector{T}
    image::Vector{T}
    volume::Array{T,D}
    face::Array{T,D}
    divergence::Array{T,D}
    current::Array{T,D}
    masked::NTuple{3,Array{T,D}}
    beta1::NTuple{3,Array{T,D}}
    beta::NTuple{3,Array{T,D}}
    candidate::NTuple{3,Array{T,D}}
    adjoint::NTuple{3,Array{T,D}}
end

function _mixed_constraint_count(op)
    cells = _mixed_count(op.n)
    faces = _mixed_count(Base.tail(op.n))
    return Base.Checked.checked_add(cells, Base.Checked.checked_mul(2, faces))
end

function MixedHallProjectionWorkspace(op::MixedShockOperators{D,T}) where {D,T}
    _mixed_storage(op)
    rows = _mixed_constraint_count(op)::Int
    A() = zeros(T, op.n)
    V() = (A(), A(), A())
    solver = Krylov.MinresQlpWorkspace(rows, rows, Vector{T})
    return MixedHallProjectionWorkspace(
        op.n,
        op.L,
        solver,
        zeros(T, rows)::Vector{T},
        zeros(T, rows)::Vector{T},
        A(),
        A(),
        A(),
        A(),
        V(),
        V(),
        V(),
        V(),
        V(),
    )
end

_mixed_projection_fields(w) = (
    w.volume,
    w.face,
    w.divergence,
    w.current,
    w.masked...,
    w.beta1...,
    w.beta...,
    w.candidate...,
    w.adjoint...,
)
_mixed_solver_full(s) = (s.wₖ₋₁, s.wₖ, s.M⁻¹vₖ₋₁, s.M⁻¹vₖ, s.x, s.p)
_mixed_solver_empty(s) = (s.Δx, s.vₖ, s.npc_dir)
_mixed_solver_stats(s) = (s.stats.residuals, s.stats.Aresiduals, s.stats.qvals, s.stats.Acond)
_mixed_projection_arrays(w) = (
    _mixed_projection_fields(w)...,
    w.rhs,
    w.image,
    _mixed_solver_full(w.solver)...,
    _mixed_solver_empty(w.solver)...,
    _mixed_solver_stats(w.solver)...,
)

function _mixed_projection_storage(w::MixedHallProjectionWorkspace{D,T}, op) where {D,T}
    w.n == op.n && w.L == op.L ||
        throw(ArgumentError("mixed projection workspace geometry changed"))
    rows = _mixed_constraint_count(op)::Int
    s = w.solver
    s isa Krylov.MinresQlpWorkspace{T,T,Vector{T}} ||
        throw(ArgumentError("mixed projection requires the concrete CPU MINRES-QLP workspace"))
    s.m == rows && s.n == rows || throw(DimensionMismatch("mixed projection solver shape changed"))
    !s.warm_start || throw(ArgumentError("mixed projection requires a cold solver workspace"))
    for a in (w.rhs, w.image, _mixed_solver_full(s)...)
        length(a) == rows && axes(a, 1) == Base.OneTo(rows) ||
            throw(DimensionMismatch("mixed projection solver vector shape changed"))
    end
    for a in (_mixed_solver_empty(s)..., _mixed_solver_stats(s)...)
        isempty(a) || throw(
            ArgumentError("mixed projection optional/history solver buffers must remain empty"),
        )
    end
    for a in _mixed_projection_fields(w)
        _mixed_array(a, op, "mixed projection private field")
    end
    return nothing
end

function _mixed_projection_preflight(out, inputs, op, w)
    _mixed_projection_storage(w, op)
    _mixed_preflight(op, out, inputs; extra = _mixed_projection_arrays(w))
    return nothing
end

function _mixed_projection_controls(
    ::Type{T},
    rtol,
    atol,
    work_rtol,
    work_atol,
    state_rtol,
    state_atol,
    maxiter,
) where {T}
    function pair(name, relative, absolute)
        r = _require_finite_nonnegative_real(name * "_rtol", relative, T)
        a = _require_finite_nonnegative_real(name * "_atol", absolute, T)
        r < one(T) || throw(ArgumentError("$name relative tolerance must be below one"))
        r > zero(T) || a > zero(T) || throw(ArgumentError("$name needs a positive tolerance"))
        return (_exact_rational(r), _exact_rational(a))
    end
    r, a = pair("projection", rtol, atol)
    wr, wa = pair("projection work", work_rtol, work_atol)
    sr, sa = pair("projection state", state_rtol, state_atol)
    return (
        rtol = r,
        atol = a,
        work_rtol = wr,
        work_atol = wa,
        state_rtol = sr,
        state_atol = sa,
        maxiter = _require_positive_intlike("maxiter", maxiter),
    )
end

function _mixed_mask!(out, input, op::MixedShockOperators{D,T}) where {D,T}
    for c = 1:3, I in CartesianIndices(out[c])
        normal_face = c == 1 && (I[1] == 1 || I[1] == op.n[1])
        out[c][I] = normal_face ? zero(T) : input[c][I]
    end
    return out
end

function _mixed_mask_qinverse!(vector, op::MixedShockOperators{D,T}) where {D,T}
    for c = 1:3, I in CartesianIndices(vector[c])
        face = I[1] == 1 || I[1] == op.n[1]
        vector[c][I] = c == 1 && face ? zero(T) : face ? T(2) * vector[c][I] : vector[c][I]
    end
    foreach(a -> _mixed_finite(a, "mixed projected adjoint"), vector)
    return vector
end

function _mixed_normal_current!(out, vector, op)
    _mixed_derivative!(out, vector[3], op, 2)
    _mixed_derivative!(op.temp, vector[2], op, 3)
    out .-= op.temp
    _mixed_finite(out, "mixed normal current")
    return out
end

# Gr = (dx*div, dx*left/right Jx). Normal B DOFs are masked structurally.
function _mixed_constraints!(rows, vector, op::MixedShockOperators{D,T}, w) where {D,T}
    _mixed_divergence!(w.divergence, vector, op)
    _mixed_normal_current!(w.current, vector, op)
    cells = length(w.divergence)
    for i = 1:cells
        rows[i] = op.dx[1] * w.divergence[i]
    end
    index = cells
    for side in (1, op.n[1]), I in CartesianIndices(Base.tail(op.n))
        index += 1
        rows[index] = op.dx[1] * w.current[CartesianIndex((side, Tuple(I)...))]
    end
    _mixed_finite(rows, "mixed scaled constraints")
    return rows
end

# Ordinary Euclidean Gr^T; Q inverse and Rn follow separately.
function _mixed_constraint_transpose!(out, rows, op::MixedShockOperators{D,T}, w) where {D,T}
    cells = length(w.volume)
    for i = 1:cells
        w.volume[i] = rows[i]
    end
    fill!(w.face, zero(T))
    index = cells
    for side in (1, op.n[1]), I in CartesianIndices(Base.tail(op.n))
        index += 1
        w.face[CartesianIndex((side, Tuple(I)...))] = rows[index]
    end
    for c = 1:3
        _mixed_transpose!(out[c], w.volume, op, c)
    end
    _mixed_transpose!(op.temp, w.face, op, 3)
    out[2] .-= op.temp
    _mixed_transpose!(op.temp, w.face, op, 2)
    out[3] .+= op.temp
    for a in out
        a .*= op.dx[1]
        _mixed_finite(a, "mixed constraint transpose")
    end
    return out
end

struct _MixedHallNormal{T,O,W}
    op::O
    work::W
end
Base.size(A::_MixedHallNormal) = (length(A.work.rhs), length(A.work.rhs))
Base.size(A::_MixedHallNormal, d::Integer) = d <= 2 ? size(A)[d] : 1
Base.eltype(::_MixedHallNormal{T}) where {T} = T

function LinearAlgebra.mul!(out::Vector{T}, A::_MixedHallNormal{T}, input::Vector{T}) where {T}
    _mixed_finite(input, "mixed multiplier matvec input")
    _mixed_constraint_transpose!(A.work.adjoint, input, A.op, A.work)
    _mixed_mask_qinverse!(A.work.adjoint, A.op)
    _mixed_constraints!(out, A.work.adjoint, A.op, A.work)
    return out
end

function _mixed_source_scales(op)
    dx = _exact_rational(op.dx[1])
    transverse = zero(Rational{BigInt})
    for work in op.transverse, k in work.ky
        transverse += abs(_exact_rational(k))
    end
    return (divergence = dx * (2 / dx + transverse), current = dx * transverse)
end

_mixed_maxabs(a) = maximum(abs, a)
_mixed_vector_maxabs(vector) = maximum(_mixed_maxabs, vector)

function _mixed_direct_constraints(vector, op, w, relative, absolute, source)
    _mixed_divergence!(w.divergence, vector, op)
    _mixed_normal_current!(w.current, vector, op)
    div = _exact_rational(_mixed_maxabs(w.divergence))
    current = zero(Rational{BigInt})
    normal_zero = true
    for side in (1, op.n[1]), I in CartesianIndices(Base.tail(op.n))
        K = CartesianIndex((side, Tuple(I)...))
        current = max(current, abs(_exact_rational(w.current[K])))
        normal_zero &= _particle_stored_zero(vector[1][K])
    end
    dx = _exact_rational(op.dx[1])
    scales = _mixed_source_scales(op)
    div_tolerance = absolute + relative * scales.divergence * source
    current_tolerance = absolute + relative * scales.current * source
    return (
        normal_zero = normal_zero,
        divergence = div,
        normal_current = current,
        scaled_divergence = dx * div,
        scaled_normal_current = dx * current,
        divergence_tolerance = div_tolerance,
        normal_current_tolerance = current_tolerance,
        passed = normal_zero && dx * div <= div_tolerance && dx * current <= current_tolerance,
    )
end

# Exact products of finite stored values; no rounded tolerance/subtraction accepts a solve.
function _mixed_true_residual(rhs, image, controls)
    residual = zero(Rational{BigInt})
    magnitude = zero(Rational{BigInt})
    for i in eachindex(rhs, image)
        b, value = _exact_rational(rhs[i]), _exact_rational(image[i])
        residual = max(residual, abs(b - value))
        magnitude = max(magnitude, abs(b))
    end
    tolerance = controls.atol + controls.rtol * magnitude
    return (residual = residual, tolerance = tolerance, passed = residual <= tolerance)
end

function _mixed_work_gate(value, norm1, norm2, relative, absolute)
    magnitude = abs(value)
    return magnitude <= absolute || (magnitude - absolute)^2 <= relative^2 * norm1 * norm2
end

function _mixed_actual_difference(B, candidate, input, op)
    dot, norm = zero(Rational{BigInt}), zero(Rational{BigInt})
    for c = 1:3, I in CartesianIndices(input[c])
        q = I[1] == 1 || I[1] == op.n[1] ? BigInt(1) // BigInt(2) : one(Rational{BigInt})
        difference = _exact_rational(candidate[c][I]) - _exact_rational(input[c][I])
        dot += q * _exact_rational(B[c][I]) * difference
        norm += q * difference^2
    end
    return dot, norm
end

function _mixed_subtraction_defect(input, candidate, beta)
    defect = zero(Rational{BigInt})
    for c = 1:3, i in eachindex(input[c])
        defect = max(
            defect,
            abs(
                _exact_rational(input[c][i]) - _exact_rational(beta[c][i]) -
                _exact_rational(candidate[c][i]),
            ),
        )
    end
    return defect
end

function _mixed_project_candidate!(input, op::MixedShockOperators{D,T}, w, controls) where {D,T}
    _mixed_mask!(w.masked, input, op)
    _mixed_constraints!(w.rhs, w.masked, op, w)
    operator = _MixedHallNormal{T,typeof(op),typeof(w)}(op, w)
    # No x0, preconditioner, shift, linesearch or history. The validated cold
    # workspace keeps the six-vector path; solver status is never sufficient.
    Krylov.minres_qlp!(
        w.solver,
        operator,
        w.rhs;
        M = I,
        λ = zero(T),
        linesearch = false,
        history = false,
        Artol = zero(T),
        rtol = T(controls.rtol),
        atol = T(controls.atol),
        itmax = controls.maxiter,
    )
    _mixed_finite(w.solver.x, "mixed solved multiplier")
    _mixed_constraint_transpose!(w.beta1, w.solver.x, op, w)
    _mixed_mask_qinverse!(w.beta1, op)
    for c = 1:3, I in CartesianIndices(input[c])
        normal_face = c == 1 && (I[1] == 1 || I[1] == op.n[1])
        w.beta[c][I] = normal_face ? input[c][I] : w.beta1[c][I]
        w.candidate[c][I] = w.masked[c][I] - w.beta1[c][I]
    end
    foreach(a -> _mixed_finite(a, "mixed projected candidate"), (w.beta..., w.candidate...))
    mul!(w.image, operator, w.solver.x)
    equation = _mixed_true_residual(w.rhs, w.image, controls)
    source = _exact_rational(max(_mixed_vector_maxabs(w.masked), _mixed_vector_maxabs(w.beta1)))
    constraints =
        _mixed_direct_constraints(w.candidate, op, w, controls.rtol, controls.atol, source)
    return (
        equation = equation,
        constraints = constraints,
        iterations = w.solver.stats.niter,
        solver_solved = w.solver.stats.solved,
        solver_status = w.solver.stats.status,
        subtraction_defect = _mixed_subtraction_defect(input, w.candidate, w.beta),
    )
end

"""
    project_mixed_field!(out, B, op, work; rtol=100eps(T), atol=0,
        work_rtol=100eps(T), work_atol=0, maxiter=200)

Project onto mixed divergence-free fields with both face Bx and Jx zero. Uses
normal-face DOF elimination and a consistent singular Q-orthogonal multiplier
solve, retaining DC/Nyquist nulls. True equation, fresh direct constraints,
orthogonality and energy gates must all pass before any output copies. Returns
exact normalized-Q residual/work/energy ledgers and a borrowed private correction.
This is a restrictive perpendicular boundary space, not a Hall shock solver.
"""
function project_mixed_field!(
    out::NTuple{3,Array{T,D}},
    B::NTuple{3,Array{T,D}},
    op::MixedShockOperators{D,T},
    work::MixedHallProjectionWorkspace{D,T};
    rtol = 100eps(T),
    atol = zero(T),
    work_rtol = 100eps(T),
    work_atol = zero(T),
    maxiter = 200,
) where {D,T}
    _mixed_projection_preflight(out, B, op, work)
    controls = _mixed_projection_controls(T, rtol, atol, work_rtol, work_atol, rtol, atol, maxiter)
    solve = _mixed_project_candidate!(B, op, work, controls)
    before = _mixed_qdot(B, B, op)
    after = _mixed_qdot(work.candidate, work.candidate, op)
    correction = _mixed_qdot(work.beta, work.beta, op)
    dot = _mixed_qdot(work.candidate, work.beta, op)
    orthogonal = _mixed_work_gate(dot, after, correction, controls.work_rtol, controls.work_atol)
    energy = after / 2 <= before / 2 + controls.work_atol + controls.work_rtol * before / 2
    report = (
        solve = solve,
        orthogonality_work = dot,
        energy_before = before / 2,
        energy_after = after / 2,
        correction_energy = correction / 2,
        energy_identity_defect = (before - after - correction) / 2,
        work_passed = orthogonal,
        energy_passed = energy,
    )
    solve.solver_solved &&
        solve.equation.passed &&
        solve.constraints.passed &&
        orthogonal &&
        energy ||
        throw(MixedHallProjectionError("field equation/constraint/work/energy gate", report))
    published_report = merge(report, (correction = work.beta,))
    for c = 1:3
        copyto!(out[c], work.candidate[c])
    end
    return published_report
end

"""
    project_mixed_rate!(out, rate, B, op, work; rtol=100eps(T), atol=0,
        state_rtol=100eps(T), state_atol=0, work_rtol=100eps(T), work_atol=0,
        maxiter=200)

Project a magnetic rate into the same perpendicular constraint space. B must
already have stored-zero normal faces and pass its separate direct constraints.
Gate the actual represented change in dot_Q(B,rate); exact power preservation
only follows for B truly in the algebraic space and an exact projection. All
reported finite-rounding/state/solve defects remain explicit.
"""
function project_mixed_rate!(
    out::NTuple{3,Array{T,D}},
    rate::NTuple{3,Array{T,D}},
    B::NTuple{3,Array{T,D}},
    op::MixedShockOperators{D,T},
    work::MixedHallProjectionWorkspace{D,T};
    rtol = 100eps(T),
    atol = zero(T),
    state_rtol = 100eps(T),
    state_atol = zero(T),
    work_rtol = 100eps(T),
    work_atol = zero(T),
    maxiter = 200,
) where {D,T}
    _mixed_projection_preflight(out, (rate..., B...), op, work)
    controls = _mixed_projection_controls(
        T,
        rtol,
        atol,
        work_rtol,
        work_atol,
        state_rtol,
        state_atol,
        maxiter,
    )
    source = _exact_rational(_mixed_vector_maxabs(B))
    state = _mixed_direct_constraints(B, op, work, controls.state_rtol, controls.state_atol, source)
    state.passed || throw(MixedHallProjectionError("input B state constraint gate", state))
    solve = _mixed_project_candidate!(rate, op, work, controls)
    power, actual_correction_norm = _mixed_actual_difference(B, work.candidate, rate, op)
    state_norm = _mixed_qdot(B, B, op)
    power_ok = _mixed_work_gate(
        power,
        state_norm,
        actual_correction_norm,
        controls.work_rtol,
        controls.work_atol,
    )
    report = (
        solve = solve,
        state = state,
        power_change = power,
        state_norm_squared = state_norm,
        actual_correction_norm_squared = actual_correction_norm,
        work_passed = power_ok,
    )
    solve.solver_solved && solve.equation.passed && solve.constraints.passed && power_ok ||
        throw(MixedHallProjectionError("rate equation/constraint/power gate", report))
    published_report = merge(report, (correction = work.beta,))
    for c = 1:3
        copyto!(out[c], work.candidate[c])
    end
    return published_report
end
