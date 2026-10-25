# A single-exponential fit of measured complex Fourier coefficients. Theory and
# solver state deliberately do not enter this owner.

function _mode_line_fit(values, dt)
    n = length(values)
    center = sum(values) / n
    xx = xy = 0.0
    for (i, value) in enumerate(values)
        x = (i - 1) / (n - 1) - 0.5
        xx += x * x
        xy += x * (value - center)
    end
    slope = xy / xx
    residual2 = 0.0
    for (i, value) in enumerate(values)
        x = (i - 1) / (n - 1) - 0.5
        residual2 += (value - (center + slope * x))^2
    end
    duration = (n - 1) * dt
    rate = slope / duration
    stderr = sqrt(residual2 / ((n - 2) * xx)) / duration
    all(isfinite, (rate, stderr)) ||
        throw(ArgumentError("mode fit rate or standard error is outside Float64 range"))
    return (; rate, center, rms = sqrt(residual2 / n), stderr)
end

"""
    _fit_mode_history(series, dt)

Fit `A exp((gamma - im*omega)*t)` to a uniformly sampled complex history without
modifying it. Return signed `omega`, `gamma`, phase/log-magnitude residual RMS,
conditional OLS slope standard errors and absolute differences between slopes
fitted to the first and second halves. Those differences measure window
sensitivity; neither diagnostic is a statistical confidence bound for correlated
simulation noise.

At least eight finite, nonzero, Float64-representable samples are required.
`dt` and the window duration must be finite and positive in Float64. Adjacent
wrapped phase changes must be less than pi/2; sampling refinement is still
required to rule out an aliased mode. Multimode histories may fit poorly: callers
must gate residuals and window sensitivity instead of accepting rates alone.
"""
function _fit_mode_history(series::AbstractVector{<:Complex}, dt::Real)
    n = length(series)
    n >= 8 || throw(ArgumentError("mode fit requires at least eight samples"))
    interval = Float64(dt)
    duration = (n - 1) * interval
    isfinite(interval) && interval > 0 && isfinite(duration) && duration > 0 ||
        throw(ArgumentError("mode fit interval and duration must be finite and positive"))
    phase = Vector{Float64}(undef, n)
    logamp = similar(phase)
    previous = 0.0
    phase_compensation = 0.0
    for (i, value) in enumerate(series)
        re, im = Float64(real(value)), Float64(imag(value))
        isfinite(re) && isfinite(im) ||
            throw(ArgumentError("mode fit sample $i must be finite and representable in Float64"))
        scale = max(abs(re), abs(im))
        scale > 0 || throw(ArgumentError("mode fit sample $i has zero amplitude"))
        # Keep the magnitude finite even when hypot(re,im) would overflow.
        logamp[i] = log(scale) + log(hypot(re / scale, im / scale))
        angle = atan(im, re)
        if i == 1
            phase[i] = angle
        else
            difference = angle - previous
            delta = atan(sin(difference), cos(difference))
            abs(delta) < pi / 2 ||
                throw(ArgumentError("mode fit phase change at sample $i requires finer sampling"))
            # Compensated accumulation prevents a long history from inheriting
            # the rounding error of every preceding phase addition.
            adjusted = delta - phase_compensation
            phase[i] = phase[i-1] + adjusted
            phase_compensation = (phase[i] - phase[i-1]) - adjusted
        end
        previous = angle
    end
    pf, af = _mode_line_fit(phase, interval), _mode_line_fit(logamp, interval)
    middle = n ÷ 2
    p1 = _mode_line_fit(@view(phase[1:middle]), interval)
    p2 = _mode_line_fit(@view(phase[(middle+1):n]), interval)
    a1 = _mode_line_fit(@view(logamp[1:middle]), interval)
    a2 = _mode_line_fit(@view(logamp[(middle+1):n]), interval)
    return (
        omega = -pf.rate,
        gamma = af.rate,
        phase_rms = pf.rms,
        logamp_rms = af.rms,
        omega_stderr = pf.stderr,
        gamma_stderr = af.stderr,
        omega_window_difference = abs(p1.rate - p2.rate),
        gamma_window_difference = abs(a1.rate - a2.rate),
        periods = abs(pf.rate) * duration / (2pi),
        samples = n,
        duration = duration,
    )
end
