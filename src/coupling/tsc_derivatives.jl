# Dimensionless derivative of the existing TSC basis. The caller obtains base
# and values from _stencil1d(TSC(), s); no second value-stencil owner is added.
@inline function _tsc_derivative_pair(s::T, base::Int) where {T}
    delta = s - T(base + 1)
    return (delta = delta, left = delta - T(0.5), right = delta + T(0.5))
end
