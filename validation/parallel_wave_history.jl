# Production cold parallel-wave histories. The independent 2x2 reference is
# used for initialization and comparison, never to construct the measurements.

function _parallel_wave_reference(k, eta, branch; coupling = 1.0)
    branch in (:whistler, :ioncyclotron) || throw(ArgumentError("unknown parallel wave branch"))
    operator = ComplexF64[im*(coupling-1) im*k*coupling; im*k im*k^2-eta*k^2]
    roots = eigvals(operator)
    index = branch == :whistler ? argmax(imag.(roots)) : argmin(imag.(roots))
    rate = roots[index]
    velocity_ratio = im * k * coupling / (rate - im * (coupling - 1))
    return (; rate, velocity_ratio, coupling)
end

_advance_parallel_wave!(st::HallMHDState, ::Nothing, dt, subcycles) = step_hall_mhd!(st, dt)
_advance_parallel_wave!(st::HybridStepper, ps::ParticleSet, dt, subcycles) =
    step!(st, ps, dt; NB = subcycles)

function _sample_parallel_wave!(state, particles, mode, dt, steps, subcycles)
    times = Vector{Float64}(undef, steps + 1)
    samples = Vector{ComplexF64}(undef, steps + 1)
    density_error = 0.0
    for j = 0:steps
        j > 0 && _advance_parallel_wave!(state, particles, dt, subcycles)
        times[j+1] = state.time[]
        samples[j+1] =
            (
                mode_amplitude(state.fields.B[2], state.g, (mode,)) +
                im * mode_amplitude(state.fields.B[3], state.g, (mode,))
            ) / state.g.n[1]
        density_error = max(density_error, maximum(x -> abs(x - 1), state.fields.n))
    end
    return (; times, samples, density_error)
end

"""Run an initialized cold circular mode through the production periodic CPU solver.

The return value contains measured times and complex amplitudes, plus a separate
independent continuum reference. Initial particle velocities are at integer time;
the production stepper owns its first leapfrog priming. No fitted/reference
frequency is used during evolution or measurement.
"""
function _parallel_wave_history(;
    solver::Symbol,
    branch::Symbol,
    mode::Int,
    n::Int = 16,
    dt::Float64 = 0.01,
    duration::Float64 = 40.0,
    nppc::Int = 16,
    eta::Float64 = 0.0,
    amplitude::Float64 = 1e-4,
    subcycles::Int = 8,
)
    solver in (:hall, :hybrid) || throw(ArgumentError("solver must be :hall or :hybrid"))
    n >= 8 && 0 < mode < n / 2 || throw(ArgumentError("mode must lie below the grid Nyquist mode"))
    nppc > 0 && subcycles > 0 ||
        throw(ArgumentError("particle count and subcycles must be positive"))
    isfinite(dt) && dt > 0 && isfinite(duration) && duration > 0 ||
        throw(ArgumentError("wave interval and duration must be finite and positive"))
    isfinite(eta) && eta >= 0 || throw(ArgumentError("resistivity must be finite and nonnegative"))
    solver == :hybrid &&
        eta != 0 &&
        throw(ArgumentError("this hybrid comparison requires zero resistivity"))
    isfinite(amplitude) && 0 < amplitude <= 1e-2 ||
        throw(ArgumentError("wave amplitude must be positive and at most 0.01"))
    count = duration / dt
    isfinite(count) && 8 <= count < typemax(Int) || throw(
        ArgumentError("wave duration must contain a representable number of at least eight steps"),
    )
    steps = round(Int, count)
    abs(steps * dt - duration) <= 16eps(duration) ||
        throw(ArgumentError("wave duration must be an integer multiple of dt"))
    state, particles, reference, discrete_reference =
        _parallel_wave_state(solver, branch, mode, n, nppc, eta, amplitude)
    history = _sample_parallel_wave!(state, particles, mode, dt, steps, subcycles)
    return merge(history, (; reference, discrete_reference, dt, steps))
end

function _parallel_wave_state(solver, branch, mode, n, nppc, eta, amplitude)
    grid = FourierGrid((n,), (2pi,))
    k = 2pi * mode / grid.L[1]
    reference = _parallel_wave_reference(k, eta, branch)
    if solver == :hall
        discrete_reference = reference
        state = HallMHDState(grid, HallMHDModel(IsothermalElectrons(0.0); η = eta))
        fill!(state.fields.n, 1)
        for j = 1:n
            wave = amplitude * cis(k * (j - 1) * grid.dx[1])
            velocity = reference.velocity_ratio * wave
            state.fields.B[1][j] = 1
            state.fields.B[2][j], state.fields.B[3][j] = real(wave), imag(wave)
            state.fields.ui[2][j], state.fields.ui[3][j] = real(velocity), imag(velocity)
        end
        particles = nothing
    else
        np = Base.checked_mul(n, nppc)
        particles = ParticleSet{1,Float64}(np)
        load_lattice_1d!(particles, 0.0, grid.L[1])
        set_density_weight!(particles, 1.0, grid)
        # Analytic CIC interpolation of the grid Fourier basis. Deposition is
        # its adjoint; this norm fixes the cold discrete eigenvector without
        # using any measured history or fitted frequency.
        gains = map(particles.x[1]) do x
            coordinate = x / grid.dx[1]
            left = floor(Int, coordinate)
            fraction = coordinate - left
            cis(k * left * grid.dx[1]) * ((1 - fraction) + fraction * cis(k * grid.dx[1]))
        end
        coupling = sum(abs2, gains) / np
        discrete_reference = _parallel_wave_reference(k, eta, branch; coupling)
        for component in particles.v
            fill!(component, 0)
        end
        for j = 1:np
            velocity = amplitude * discrete_reference.velocity_ratio * gains[j] / coupling
            particles.v[2][j], particles.v[3][j] = real(velocity), imag(velocity)
        end
        state = HybridStepper(grid, HybridModel(IsothermalElectrons(0.0)), CIC(), np)
        for j = 1:n
            wave = amplitude * cis(k * (j - 1) * grid.dx[1])
            state.fields.B[1][j] = 1
            state.fields.B[2][j], state.fields.B[3][j] = real(wave), imag(wave)
        end
        init!(state, particles)
    end
    return state, particles, reference, discrete_reference
end
