# Directed bounds used only by the BigFloat recommended_dt calculation.
# Each operation encloses the exact result at the task's working precision.
struct _TimestepBounds
    lo::BigFloat
    hi::BigFloat
    function _TimestepBounds(lo::BigFloat, hi::BigFloat)
        _timestep_check_exponent(lo)
        _timestep_check_exponent(hi)
        (!isnan(lo) && !isnan(hi) && lo <= hi) || _timestep_bound_error()
        return new(lo, hi)
    end
end

@noinline function _timestep_bound_error()
    throw(
        ArgumentError(
            "cannot certify a finite positive recommended timestep in the current BigFloat exponent range",
        ),
    )
end

# mpfr_get_exp permits retained values outside the active range; arithmetic
# and comparisons do not. Admit their metadata before either can be called.
@inline _timestep_check_exponent(x) = nothing
function _timestep_check_exponent(x::BigFloat)
    if isfinite(x) && !iszero(x)
        mpfr_exponent = exponent(x) + 1
        Base.MPFR.get_emin() <= mpfr_exponent <= Base.MPFR.get_emax() || _timestep_bound_error()
    end
    return nothing
end

function _TimestepBounds(x::Real)
    _timestep_check_exponent(x)
    # Base's directed Rational conversion is enclosing for positive values.
    # Negate a positive enclosure for negative rationals, including typemin.
    if x isa Rational && x < 0
        return -_TimestepBounds(-big(x))
    end
    lo = setrounding(BigFloat, RoundDown) do
        BigFloat(x)
    end
    hi = setrounding(BigFloat, RoundUp) do
        BigFloat(x)
    end
    return _TimestepBounds(lo, hi)
end

Base.convert(::Type{_TimestepBounds}, x::Real) = _TimestepBounds(x)
Base.one(::_TimestepBounds) = _TimestepBounds(1)
Base.zero(::Type{_TimestepBounds}) = _TimestepBounds(0)

function Base.:+(x::_TimestepBounds, y::_TimestepBounds)
    lo = setrounding(BigFloat, RoundDown) do
        x.lo + y.lo
    end
    hi = setrounding(BigFloat, RoundUp) do
        x.hi + y.hi
    end
    return _TimestepBounds(lo, hi)
end

function Base.:-(x::_TimestepBounds, y::_TimestepBounds)
    lo = setrounding(BigFloat, RoundDown) do
        x.lo - y.hi
    end
    hi = setrounding(BigFloat, RoundUp) do
        x.hi - y.lo
    end
    return _TimestepBounds(lo, hi)
end
Base.:-(x::_TimestepBounds) = _TimestepBounds(-x.hi, -x.lo)

function Base.:*(x::_TimestepBounds, y::_TimestepBounds)
    lo = setrounding(BigFloat, RoundDown) do
        min(x.lo * y.lo, x.lo * y.hi, x.hi * y.lo, x.hi * y.hi)
    end
    hi = setrounding(BigFloat, RoundUp) do
        max(x.lo * y.lo, x.lo * y.hi, x.hi * y.lo, x.hi * y.hi)
    end
    return _TimestepBounds(lo, hi)
end

function Base.:/(x::_TimestepBounds, y::_TimestepBounds)
    # The formula divides only by positive lengths, constants and rate sums.
    y.lo > 0 || _timestep_bound_error()
    lo = setrounding(BigFloat, RoundDown) do
        min(x.lo / y.lo, x.lo / y.hi, x.hi / y.lo, x.hi / y.hi)
    end
    hi = setrounding(BigFloat, RoundUp) do
        max(x.lo / y.lo, x.lo / y.hi, x.hi / y.lo, x.hi / y.hi)
    end
    return _TimestepBounds(lo, hi)
end

for operation in (:+, :-, :*, :/)
    @eval begin
        Base.$operation(x::_TimestepBounds, y::Real) = $operation(x, _TimestepBounds(y))
        Base.$operation(x::Real, y::_TimestepBounds) = $operation(_TimestepBounds(x), y)
    end
end

Base.min(x::_TimestepBounds, y::_TimestepBounds) = _TimestepBounds(min(x.lo, y.lo), min(x.hi, y.hi))
Base.max(x::_TimestepBounds, y::_TimestepBounds) = _TimestepBounds(max(x.lo, y.lo), max(x.hi, y.hi))
Base.minmax(x::_TimestepBounds, y::_TimestepBounds) = min(x, y), max(x, y)
for operation in (:min, :max, :minmax)
    @eval begin
        Base.$operation(x::_TimestepBounds, y::Real) = $operation(x, _TimestepBounds(y))
        Base.$operation(x::Real, y::_TimestepBounds) = $operation(_TimestepBounds(x), y)
    end
end

function Base.abs2(x::_TimestepBounds)
    nearest = x.lo <= 0 <= x.hi ? BigFloat(0) : min(abs(x.lo), abs(x.hi))
    farthest = max(abs(x.lo), abs(x.hi))
    lo = setrounding(BigFloat, RoundDown) do
        nearest * nearest
    end
    hi = setrounding(BigFloat, RoundUp) do
        farthest * farthest
    end
    return _TimestepBounds(lo, hi)
end

function Base.hypot(x::_TimestepBounds, y::_TimestepBounds)
    # Both arguments are nonnegative in the two Hall-wave identities.
    (x.lo >= 0 && y.lo >= 0) || _timestep_bound_error()
    lo = setrounding(BigFloat, RoundDown) do
        hypot(x.lo, y.lo)
    end
    hi = setrounding(BigFloat, RoundUp) do
        hypot(x.hi, y.hi)
    end
    return _TimestepBounds(lo, hi)
end

function Base.log(x::_TimestepBounds)
    x.lo >= 0 || _timestep_bound_error()
    lo = setrounding(BigFloat, RoundDown) do
        log(x.lo)
    end
    hi = setrounding(BigFloat, RoundUp) do
        log(x.hi)
    end
    return _TimestepBounds(lo, hi)
end

function Base.log1p(x::_TimestepBounds)
    x.lo >= -1 || _timestep_bound_error()
    lo = setrounding(BigFloat, RoundDown) do
        log1p(x.lo)
    end
    hi = setrounding(BigFloat, RoundUp) do
        log1p(x.hi)
    end
    return _TimestepBounds(lo, hi)
end

function Base.exp(x::_TimestepBounds)
    lo = setrounding(BigFloat, RoundDown) do
        exp(x.lo)
    end
    hi = setrounding(BigFloat, RoundUp) do
        exp(x.hi)
    end
    return _TimestepBounds(lo, hi)
end

# Preserve the native log dispatch and wave branch; the bounded core promotes
# these constants/subcycles before taking logs and chooses a valid identity.
@inline _timestep_log(::Type, x) = log(x)
@inline _timestep_log(::Type{_TimestepBounds}, x) = log(_TimestepBounds(x))
@inline _timestep_nonnegative(x::Real) = x >= 0
@inline _timestep_nonnegative(x::_TimestepBounds) = x.lo >= 0

function _recommended_dt_bigfloat(g, NB, integrator, safety, d_i, Ω_ci, η, ηH, C, Creal)
    caller_precision = precision(BigFloat)
    # The bounds prove safety independently of this width allowance. The extra
    # bits retain useful resolution after subtracting C-long-scale logarithms.
    working_precision = max(256, Base.checked_add(caller_precision, 96))
    for x in (g.dx..., safety, d_i, Ω_ci, η, ηH)
        x isa BigFloat && (working_precision = max(working_precision, precision(x)))
    end
    return setprecision(BigFloat, working_precision) do
        setrounding(BigFloat, RoundNearest) do
            bounds = _recommended_dt_core(
                g,
                NB,
                integrator,
                safety,
                d_i,
                Ω_ci,
                η,
                ηH,
                C,
                Creal,
                _TimestepBounds,
            )
            result = BigFloat(bounds.lo, RoundToZero; precision = caller_precision)
            (isfinite(result) && result > 0) || _timestep_bound_error()
            return result
        end
    end
end
