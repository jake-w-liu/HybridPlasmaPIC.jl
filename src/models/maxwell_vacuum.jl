"""
    MaxwellVacuum(g::YeeGrid; c=1)

Source-free Maxwell evolution on a periodic/conducting Yee mesh, with normalized
permittivity ε₀=1 and permeability μ₀=1/c². All three E/B components in `fields`
are synchronized at `time[]`; `step_maxwell!` uses magnetic half-kick, electric
drift, magnetic half-kick. Constant steps have the usual Yee dispersion relation;
changing the step preserves second-order centering. Neither the constant-step
energy invariant nor stability for arbitrary sub-CFL step sequences extends to
variable steps. CPU Float32/Float64 only. In 1D, `step_maxwell!` also accepts
prescribed integrated current arrays. No particle coupling, absorbing layer or
restart adapter is supplied by this field model.

The model owns reusable staging arrays. Concurrent mutation/stepping of one
instance is unsupported. Physical E/B arrays retain their identities on commit.
"""
struct MaxwellVacuum{D,T}
    g::YeeGrid{D,T}
    fields::YeeFields{D,T}
    c::T
    time::Base.RefValue{T}
    step::Base.RefValue{Int}
    _candidate::YeeFields{D,T}
    _cfl::T
end

function MaxwellVacuum(g::YeeGrid{D,T}; c = 1) where {D,T}
    cT = _require_finite_positive_real("c", c, T)
    inverse_squared_limit = _exact_rational(cT)^2 * sum(d -> inv(_exact_rational(g.dx[d])^2), 1:D)
    cfl = _sqrt_from_exact(T, inv(inverse_squared_limit), Int128(0), RoundDown, false)
    return MaxwellVacuum{D,T}(g, YeeFields(g), cT, Ref(zero(T)), Ref(0), YeeFields(g), cfl)
end

"""
    maxwell_cfl(model)

The constant-step Yee stability bound `1/(c*sqrt(sum(1/dx²)))`, rounded downward
once when the model is constructed. A step must be positive and smaller than this returned
bound; equality is rejected because the exact endpoint can have
a repeated amplification eigenvalue with secular growth. This bound does not
guarantee stability under arbitrary variable-step sequences and is not an
accuracy recommendation.
"""
maxwell_cfl(model::MaxwellVacuum) = model._cfl

# Field finiteness is not scanned here: every physical entry is an operand of
# a native stage whose finiteness reduction fails for a nonfinite input, and
# `_yee_advance!` then validates its inputs before any exact completion.
function _yee_validate_storage(model::MaxwellVacuum)
    _yee_check_shape(model.fields.E, model.g, true)
    _yee_check_shape(model.fields.B, model.g, false)
    _yee_validate_pec(model.fields, model.g)
    _yee_check_shape(model._candidate.E, model.g, true)
    _yee_check_shape(model._candidate.B, model.g, false)
    arrays = (model.fields.E..., model.fields.B..., model._candidate.E..., model._candidate.B...)
    _yee_check_independent_arrays(arrays)
    return nothing
end

# ------------------------------------------------------------ stage arithmetic
#
# A stage is `new = old + κ_j*(src_k[+] - src_k[-]) - κ_k*(src_j[+] - src_j[-])`
# with per-axis multipliers κ_d = -dt/(2dx_d) (magnetic half-kick, electric
# source) or dt*c²/dx_d (electric stage), computed once per stage. Ordinary
# arithmetic needs every multiplier intermediate finite and normal. Otherwise
# Float32 stages evaluate the same expression in Float64, which represents
# every binary32 product here without range loss, and Float64 stages keep
# κ_d = m_d*2^s with normal mantissas m_d >= 2^52 (no difference product can
# underflow) and apply 2^s through two exact power-of-two factors. Entries whose
# native result is nonfinite are evaluated exactly from the same inputs and
# rounded once; a nonfinite exact result rejects the step.

struct _YeeStageOp{W,scaled}
    κa::W
    κb::W
    σ1::W
    σ2::W
end

_YeeStageOp{W,scaled}(κa, κb) where {W,scaled} = _YeeStageOp{W,scaled}(κa, κb, one(W), one(W))

@inline _yee_stage_difference(κ, ::Nothing, ::Type) = nothing
@inline _yee_stage_difference(κ, pair::Tuple, ::Type{W}) where {W} = κ * (W(pair[1]) - W(pair[2]))

@inline function _yee_stage_value(op::_YeeStageOp{W,scaled}, old, a, b) where {W,scaled}
    curl = _yee_curl_combination(
        _yee_stage_difference(op.κa, a, W),
        _yee_stage_difference(op.κb, b, W),
    )
    curl === nothing && return W(old)
    return W(old) + (scaled ? (curl * op.σ1) * op.σ2 : curl)
end

@inline (op::_YeeStageOp)(old::T, a, b) where {T<:AbstractFloat} =
    T(_yee_stage_value(op, old, a, b))

# Prescribed impulse: `(old + curl) - impulse`, rounded once. A zero impulse
# entry keeps the vacuum value bitwise, including a signed zero.
@inline function (op::_YeeStageOp{W})(entry::Tuple{T,T}, a, b) where {W,T}
    old, impulse = entry
    value = _yee_stage_value(op, old, a, b)
    return T(ifelse(iszero(impulse), value, value - W(impulse)))
end

# Original entries paired with the impulse entries of the same lattice.
struct _YeeSourcedEntries{A}
    old::A
    impulse::A
end
@inline _yee_load(a::_YeeSourcedEntries, k) = @inbounds (a.old[k], a.impulse[k])

# No representable multiplier: every entry takes the exact expression.
struct _YeeExactStageOp{T} end
@inline (::_YeeExactStageOp{T})(old, a, b) where {T} = T(NaN)

@inline _yee_normal(x) = isfinite(x) & (abs(x) >= floatmin(x))

@inline function _yee_multiplier(dt::W, speed::W, dx::W, electric::Bool) where {W}
    ratio = dt / dx
    if electric
        κ = -(ratio * W(0.5))
        return κ, _yee_normal(ratio) & _yee_normal(κ)
    end
    courant = ratio * speed
    κ = courant * speed
    return κ, _yee_normal(ratio) & _yee_normal(courant) & _yee_normal(κ)
end

# Multiplier as mantissa*2^exponent; mantissa magnitudes lie in [1/8, 2).
@inline function _yee_multiplier_parts(dt::T, speed::T, dx::T, electric::Bool) where {T}
    ft, et = frexp(dt)
    fx, ex = frexp(dx)
    ratio = ft / fx
    electric && return -(ratio * T(0.5)), et - ex
    fc, ec = frexp(speed)
    return (ratio * fc) * fc, et - ex + 2ec
end

# Float64 scaled multipliers, or `nothing` if no two-factor split exists.
function _yee_scaled_multipliers(g::YeeGrid{D,T}, electric::Bool, dt::T, speed::T) where {D,T}
    parts = ntuple(d -> _yee_multiplier_parts(dt, speed, g.dx[d], electric), D)
    shift = minimum(last, parts) - 55
    maximum(last, parts) - shift <= 1020 || return nothing
    first_shift = clamp(shift, -1022, 1023)
    second_shift = shift - first_shift
    -1022 <= second_shift <= 1023 || return nothing
    mantissas = map(part -> ldexp(first(part), last(part) - shift), parts)
    return mantissas, ldexp(one(T), first_shift), ldexp(one(T), second_shift)
end

function _yee_stage_exact(old, source, p, I, c, g, electric, dt::T, speed::T) where {T}
    coefficient =
        electric ? -_exact_rational(dt) / 2 : _exact_rational(dt) * _exact_rational(speed)^2
    exact = _exact_rational(old[I]) + coefficient * _yee_curl_exact(source, I, c, g, electric)
    p === nothing || (exact -= _exact_rational(p[I]))
    value = _from_exact(T, exact)
    isfinite(value) || throw(ArgumentError("Maxwell stage produces a nonfinite field"))
    return value
end

@inline function _yee_stage_component!(
    make_op::F,
    destination,
    original,
    source,
    impulse,
    g::YeeGrid{D,T},
    direction::Val{electric},
    ::Val{c},
    careful::Bool,
    dt::T,
    speed::T,
) where {F,D,T,electric,c}
    ja, jb = mod1(c + 1, 3), mod1(c + 2, 3)
    op = make_op(ja, jb)
    p = impulse === nothing ? nothing : impulse[c]
    old = p === nothing ? original[c] : _YeeSourcedEntries(original[c], p)
    sa, sb = source[jb], source[ja]
    if careful
        exact = I -> _yee_stage_exact(original[c], source, p, I, c, g, electric, dt, speed)
        return _yee_rows!(op, exact, destination[c], old, sa, sb, g, Val(ja), Val(jb), direction)
    end
    return _yee_rows!(op, nothing, destination[c], old, sa, sb, g, Val(ja), Val(jb), direction)
end

@inline function _yee_stage_components!(
    make_op::F,
    destination,
    original,
    source,
    impulse,
    g,
    direction,
    careful,
    dt,
    speed,
) where {F}
    arguments = (destination, original, source, impulse, g, direction)
    ok = _yee_stage_component!(make_op, arguments..., Val(1), careful, dt, speed)
    ok &= _yee_stage_component!(make_op, arguments..., Val(2), careful, dt, speed)
    ok &= _yee_stage_component!(make_op, arguments..., Val(3), careful, dt, speed)
    return ok
end

# One stage into `destination`. The native pass (`careful=false`) returns
# whether every entry is finite. The careful pass reads each original entry
# before writing it (valid in place) and completes nonfinite entries exactly.
function _yee_stage!(
    destination,
    original,
    source,
    g::YeeGrid{D,T},
    direction::Val{electric},
    dt::T,
    speed::T,
    careful::Bool,
    impulse = nothing,
) where {D,T,electric}
    native = ntuple(d -> _yee_multiplier(dt, speed, g.dx[d], electric), D)
    if all(last, native)
        κ = map(first, native)
        make =
            (ja, jb) -> _YeeStageOp{T,false}(ja <= D ? κ[ja] : zero(T), jb <= D ? κ[jb] : zero(T))
        return _yee_stage_components!(
            make,
            destination,
            original,
            source,
            impulse,
            g,
            direction,
            careful,
            dt,
            speed,
        )
    end
    if T === Float32
        wide =
            ntuple(d -> _yee_multiplier(Float64(dt), Float64(speed), Float64(g.dx[d]), electric), D)
        if all(last, wide)
            ω = map(first, wide)
            make =
                (ja, jb) -> _YeeStageOp{Float64,false}(ja <= D ? ω[ja] : 0.0, jb <= D ? ω[jb] : 0.0)
            return _yee_stage_components!(
                make,
                destination,
                original,
                source,
                impulse,
                g,
                direction,
                careful,
                dt,
                speed,
            )
        end
    else
        scaled = _yee_scaled_multipliers(g, electric, dt, speed)
        if scaled !== nothing
            m, σ1, σ2 = scaled
            make =
                (ja, jb) -> _YeeStageOp{T,true}(
                    ja <= D ? m[ja] : zero(T),
                    jb <= D ? m[jb] : zero(T),
                    σ1,
                    σ2,
                )
            return _yee_stage_components!(
                make,
                destination,
                original,
                source,
                impulse,
                g,
                direction,
                careful,
                dt,
                speed,
            )
        end
    end
    make = (ja, jb) -> _YeeExactStageOp{T}()
    return _yee_stage_components!(
        make,
        destination,
        original,
        source,
        impulse,
        g,
        direction,
        careful,
        dt,
        speed,
    )
end

# Complete stage with distinct destination and original arrays: native pass,
# then the careful pass only if some native entry was nonfinite. In-place
# stages have no retained originals and are evaluated carefully. Nonfinite
# inputs are rejected before exact completion.
function _yee_advance!(
    destination,
    original,
    source,
    g::YeeGrid{D,T},
    electric::Bool,
    dt::T,
    speed::T,
    impulse = nothing,
) where {D,T}
    direction = electric ? Val(true) : Val(false)
    in_place = any(c -> Base.mightalias(destination[c], original[c]), 1:3)
    !in_place &&
        _yee_stage!(destination, original, source, g, direction, dt, speed, false, impulse) &&
        return destination
    for a in (original..., source...)
        _yee_all_finite(a) || throw(ArgumentError("Yee fields must be finite"))
    end
    _yee_stage!(destination, original, source, g, direction, dt, speed, true, impulse)
    return destination
end

"""
    step_maxwell!(model, dt)

Advance synchronized vacuum E/B by one positive, strictly sub-CFL step. Shape,
alias, PEC, clock and step-size checks happen before staging writes; nonfinite
fields are detected by the stage finiteness reductions and rejected before any
exact completion. A nonfinite input or exact stage result throws
`ArgumentError`; physical fields, time and step count remain unchanged on every
rejected step. Internal staging buffers may change after a rejected step and
are overwritten on reuse.

Each stage is the native stencil `old + κ_j*Δ_j - κ_k*Δ_k` with per-axis
multipliers computed once per stage (`κ_d = -dt/(2dx_d)` or `dt*c²/dx_d`).
Normwise contract: each stage agrees with its exact expression, from that
stage's inputs, to a few `eps(T)` times `|old| + Σ|κ_d*Δ_d|` plus a few
subnormal quanta. Multipliers outside the normal range are applied in Float64
(Float32 fields) or as exact powers of two (Float64 fields), and entries whose
native value overflows are evaluated exactly and rounded once; finite exact
stages therefore never become nonfinite. This is not a bound to an infinitely
precise trajectory over multiple stages or time steps, and componentwise
relative accuracy of a nearly cancelling entry is not certified.
No Gauss projection is applied: the compatible curls propagate initial bulk
electric and magnetic divergence, including a nonzero initial residual.
"""
function step_maxwell!(model::MaxwellVacuum{D,T}, dt::Real) where {D,T}
    return _yee_step_maxwell!(model, dt, nothing)
end

"""
    step_maxwell!(model, dt, impulse)

Advance a 1D Yee field with prescribed integrated current `impulse = integral(J dt)`
in normalized epsilon0=1 units. Supply a tuple of three ordinary CPU vectors with
the electric lattice shapes and precision. The electric stage is
`Enew = (Eold + dt*c^2*curl(Bhalf)) - impulse`, evaluated natively (Float64 for
out-of-range Float32 multipliers) and rounded once; overflowing entries are
completed exactly with the impulse included. The normwise stage contract of the
vacuum step applies with `|impulse|` added to the scale, as do its CFL, clock and
physical-state rollback contracts.

Sources must be finite, independent of one another and of all physical/staging
arrays, and zero on constrained tangential PEC nodes; they are validated before
any staging write and never modified. A zero impulse entry reproduces the vacuum
value bitwise. This is prescribed forcing, not particle coupling or an automatic
charge/Gauss projection. See [`yee_electric_charge!`](@ref),
[`yee_source_work`](@ref) and [`yee_modified_energy`](@ref) for the associated
diagnostics.
"""
function step_maxwell!(model::MaxwellVacuum{1}, dt::Real, impulse)
    impulse === nothing && throw(ArgumentError("impulse must be a tuple of three current arrays"))
    return _yee_step_maxwell!(model, dt, impulse)
end

function step_maxwell!(model::MaxwellVacuum, dt::Real, impulse)
    throw(ArgumentError("prescribed Yee currents support one spatial dimension"))
end

function _yee_step_maxwell!(model::MaxwellVacuum{D,T}, dt::Real, impulse) where {D,T}
    dtT = _require_finite_positive_real("dt", dt, T)
    dtT < maxwell_cfl(model) ||
        throw(ArgumentError("dt must be strictly below the Maxwell CFL limit"))
    next_time, next_step = _validated_step_clock(model.time[], model.step[], dtT)
    next_time > model.time[] || throw(ArgumentError("Maxwell time must advance finitely"))
    _yee_validate_storage(model)
    if impulse !== nothing
        _yee_validate_impulse(impulse, model.g)
        _yee_check_independent_arrays((
            model.fields.E...,
            model.fields.B...,
            model._candidate.E...,
            model._candidate.B...,
            impulse...,
        ))
    end
    candidate, original = model._candidate, model.fields
    g, speed = model.g, model.c
    _yee_advance!(candidate.B, original.B, original.E, g, true, dtT, speed)
    _yee_advance!(candidate.E, original.E, candidate.B, g, false, dtT, speed, impulse)
    if !_yee_stage!(candidate.B, candidate.B, candidate.E, g, Val(true), dtT, speed, false)
        # The in-place native half-kick overwrote its inputs. Rebuild the
        # deterministic first half-kick from the unchanged physical fields and
        # evaluate the final half-kick carefully.
        _yee_advance!(candidate.B, original.B, original.E, g, true, dtT, speed)
        _yee_stage!(candidate.B, candidate.B, candidate.E, g, Val(true), dtT, speed, true)
    end
    # Homogeneous PEC data stay exactly zero: every constrained entry is
    # `0 + κ*(0 - 0)` (or its exact replay) because each difference it uses is
    # of constrained samples or is absent, and constrained impulses are zero;
    # no clamp or re-validation is needed.
    for c = 1:3
        copyto!(original.E[c], candidate.E[c])
        copyto!(original.B[c], candidate.B[c])
    end
    model.time[] = next_time
    model.step[] = next_step
    return model
end
