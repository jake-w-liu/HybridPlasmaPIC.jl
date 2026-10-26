# Validation-only initialization and measurement of oblique production Hall waves.
using LinearAlgebra
if !isdefined(@__MODULE__, :HybridDispersionOracle)
    include(joinpath(@__DIR__, "..", "test", "oracles", "hybrid_dispersion_oracle.jl"))
end

const _OBLIQUE_COMPONENTS = ("density", "ux", "uy", "uz", "bx", "by", "bz")

function _oblique_spectrum(k, angle, pressure, direction)
    background = (cosd(angle), sind(angle), 0.0)
    matrix = HybridDispersionOracle.hybrid_dispersion_matrix(
        (k, 0.0, 0.0),
        background;
        cs = sqrt(pressure),
    )
    decomposition = eigen(matrix)
    scale = max(1.0, maximum(abs, decomposition.values))
    all(z -> abs(imag(z)) <= 1e-10scale, decomposition.values) ||
        throw(ArgumentError("ideal Hall spectrum contains nonreal frequencies"))
    indices = findall(z -> direction * real(z) > 1e-9scale, decomposition.values)
    length(indices) == (iszero(pressure) ? 2 : 3) ||
        throw(ArgumentError("Hall spectrum lacks the required separated propagating roots"))
    roots = real.(decomposition.values[indices])
    gap =
        minimum(abs(roots[i] - roots[j]) for i in eachindex(roots), j in eachindex(roots) if i != j)
    gap > 1e-6scale || throw(ArgumentError("Hall propagating roots are degenerate or too close"))
    vectors = Matrix{ComplexF64}(undef, 7, length(indices))
    for (column, index) in enumerate(indices)
        v = decomposition.vectors[:, index]
        density = k * v[1] / roots[column]
        abs(v[7] - pressure * density) <= 1e-10max(1.0, abs(v[7]), abs(pressure * density)) ||
            throw(ArgumentError("Hall eigenvector violates the scalar pressure closure"))
        abs(v[4]) <= 1e-10norm(v) || throw(ArgumentError("Hall eigenvector is not divergence free"))
        physical = vcat(density, v[1:6])
        physical ./= norm(physical)
        all(isfinite, physical) || throw(ArgumentError("Hall eigenvector is not finite"))
        vectors[:, column] = physical
    end
    return (; roots, vectors, gap)
end

function _oblique_assignment(previous, current)
    count = size(previous, 2)
    # There are exactly two cold or three warm propagating branches per sign.
    permutations =
        count == 2 ? ((1, 2), (2, 1)) :
        ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
    overlap = abs.(previous' * current)
    scores = [sum(overlap[j, p[j]] for j = 1:count) for p in permutations]
    order = sortperm(scores; rev = true)
    gap = scores[order[1]] - scores[order[2]]
    assignment = permutations[order[1]]
    selected = [overlap[j, assignment[j]] for j = 1:count]
    minimum(selected) >= 0.8 && gap > 1e-3 ||
        throw(ArgumentError("Hall branch continuation is ambiguous; refine the reference path"))
    return (; assignment, overlap = selected, gap)
end

"""Independent, explicitly tracked ideal Hall eigenmode; no production state is read."""
function _oblique_hall_reference(k, angle, Ti, Te, branch, direction)
    all(isfinite, (k, angle, Ti, Te)) && k > 0 && 0 <= angle < 90 && Ti >= 0 && Te >= 0 || throw(
        ArgumentError(
            "Hall reference requires finite k>0, 0<=angle<90 and nonnegative temperatures",
        ),
    )
    direction in (-1, 1) || throw(ArgumentError("frequency direction must be -1 or 1"))
    pressure = Float64(Ti + Te)
    isfinite(pressure) || throw(ArgumentError("total temperature must be finite"))
    names = iszero(pressure) ? (:whistler, :ioncyclotron) : (:whistler, :ioncyclotron, :acoustic)
    index = findfirst(==(branch), names)
    index === nothing && throw(ArgumentError("branch is not propagating for this temperature"))
    seeds = [(sqrt(k^4 + 4k^2) + k^2) / 2, (sqrt(k^4 + 4k^2) - k^2) / 2]
    iszero(pressure) || push!(seeds, k * sqrt(pressure))
    seeds .*= direction
    all(isfinite, seeds) || throw(ArgumentError("Hall seed frequency is outside Float64 range"))
    scale = max(1.0, maximum(abs, seeds))
    minimum(abs(seeds[i] - seeds[j]) for i in eachindex(seeds), j in eachindex(seeds) if i != j) >
    1e-6scale || throw(ArgumentError("parallel Hall seed branches are degenerate or too close"))
    seed_spectrum = _oblique_spectrum(k, 0.0, pressure, direction)
    assignment = [argmin(abs.(seed_spectrum.roots .- seed)) for seed in seeds]
    length(unique(assignment)) == length(seeds) ||
        throw(ArgumentError("parallel Hall seeds do not identify distinct roots"))
    maximum(abs.(seed_spectrum.roots[assignment] .- seeds)) <= 1e-10scale ||
        throw(ArgumentError("parallel Hall roots disagree with closed-form seeds"))
    roots = seed_spectrum.roots[assignment]
    vectors = seed_spectrum.vectors[:, assignment]
    trace = NamedTuple[]
    for j in eachindex(names)
        # No overlap assignment is made at the closed-form seed node.
        push!(
            trace,
            (;
                angle = 0.0,
                branch = names[j],
                omega = roots[j],
                overlap = 1.0,
                assignment_gap = NaN,
                spectral_gap = seed_spectrum.gap,
            ),
        )
    end
    nodes = ceil(Int, 10angle)
    for node = 1:nodes
        theta = min(node / 10, Float64(angle))
        spectrum = _oblique_spectrum(k, theta, pressure, direction)
        matched = _oblique_assignment(vectors, spectrum.vectors)
        next_vectors = spectrum.vectors[:, collect(matched.assignment)]
        roots = spectrum.roots[collect(matched.assignment)]
        for j in eachindex(names)
            column = @view next_vectors[:, j]
            phase = dot(@view(vectors[:, j]), column)
            column .*= conj(phase) / abs(phase)
            push!(
                trace,
                (;
                    angle = theta,
                    branch = names[j],
                    omega = roots[j],
                    overlap = matched.overlap[j],
                    assignment_gap = matched.gap,
                    spectral_gap = spectrum.gap,
                ),
            )
        end
        vectors = next_vectors
    end
    vector = vectors[:, index]
    component = argmax(abs.(vector))
    vector .*= conj(vector[component]) / abs(vector[component])
    return (;
        omega = roots[index],
        gamma = 0.0,
        rate = -im * roots[index],
        vector,
        component,
        trace,
        background = (cosd(angle), sind(angle), 0.0),
    )
end

function _oblique_hall_state(n, mode, angle, Ti, Te, branch, direction, amplitude)
    reference = _oblique_hall_reference(Float64(mode), angle, Ti, Te, branch, direction)
    grid = FourierGrid((n,), (2pi,))
    state = HallMHDState(grid, HallMHDModel(IsothermalElectrons(Te); Ti))
    fields = (state.fields.n, state.fields.ui..., state.fields.B...)
    background = (1.0, 0.0, 0.0, 0.0, reference.background...)
    for j = 1:n
        phase = cis(mode * (j - 1) * grid.dx[1])
        for c = 1:7
            fields[c][j] = background[c] + amplitude * real(reference.vector[c] * phase)
        end
    end
    hall_mhd_ohms_law!(state)
    return state, reference
end

function _oblique_mode_vector(state, mode)
    fields = (state.fields.n, state.fields.ui..., state.fields.B...)
    return ComplexF64[2mode_amplitude(f, state.g, (mode,)) / state.g.n[1] for f in fields]
end

function _sample_oblique_hall!(state, reference, mode, dt, steps, io)
    times = Vector{Float64}(undef, steps + 1)
    samples = Matrix{ComplexF64}(undef, 7, steps + 1)
    mass0 = sum(state.fields.n)
    mass_error = longitudinal_error = 0.0
    density_min = minimum(state.fields.n)
    density_max = maximum(state.fields.n)
    floor_count = 0
    if io !== nothing
        header = vcat(
            "time",
            ["$(name)_$(part)" for name in _OBLIQUE_COMPONENTS for part in ("real", "imag")],
        )
        println(io, join(header, ','))
    end
    for j = 0:steps
        j > 0 && step_hall_mhd!(state, dt)
        times[j+1] = state.time[]
        samples[:, j+1] = _oblique_mode_vector(state, mode)
        mass_error = max(mass_error, abs(sum(state.fields.n) / mass0 - 1))
        longitudinal_error = max(
            longitudinal_error,
            maximum(x -> abs(x - reference.background[1]), state.fields.B[1]),
        )
        density_min = min(density_min, minimum(state.fields.n))
        density_max = max(density_max, maximum(state.fields.n))
        floor_count = max(floor_count, state.fields.floor_count[])
        if io !== nothing
            print(io, times[j+1])
            for z in @view samples[:, j+1]
                print(io, ',', real(z), ',', imag(z))
            end
            println(io)
        end
    end
    return (; times, samples, mass_error, longitudinal_error, density_min, density_max, floor_count)
end

"""Evolve the actual periodic Hall-MHD solver and retain measured Fourier vectors.

Only the initial state uses the independent eigenmode. An optional output path
streams successful raw samples and is closed even when evolution throws.
"""
function _oblique_hall_history(;
    branch::Symbol,
    angle::Float64,
    Ti::Float64,
    Te::Float64,
    direction::Int = 1,
    mode::Int = 1,
    n::Int = 16,
    dt::Float64 = 0.01,
    duration::Float64 = 80.0,
    amplitude::Float64 = 1e-5,
    output::Union{Nothing,AbstractString} = nothing,
)
    n >= 8 && 0 < mode < n / 2 || throw(ArgumentError("mode must lie below the grid Nyquist mode"))
    isfinite(dt) && dt > 0 && isfinite(duration) && duration > 0 ||
        throw(ArgumentError("history interval and duration must be finite and positive"))
    isfinite(amplitude) && 0 < amplitude <= 1e-3 ||
        throw(ArgumentError("oblique perturbation amplitude must be positive and at most 0.001"))
    count = duration / dt
    isfinite(count) && 8 <= count < typemax(Int) ||
        throw(ArgumentError("history requires a representable count of at least eight steps"))
    steps = round(Int, count)
    abs(steps * dt - duration) <= 16eps(duration) ||
        throw(ArgumentError("history duration must be an integer multiple of dt"))
    output === nothing || !ispath(output) || throw(ArgumentError("raw history path already exists"))
    state, reference = _oblique_hall_state(n, mode, angle, Ti, Te, branch, direction, amplitude)
    history = if output === nothing
        _sample_oblique_hall!(state, reference, mode, dt, steps, nothing)
    else
        open(io -> _sample_oblique_hall!(state, reference, mode, dt, steps, io), output, "w")
    end
    return merge(history, (; reference, dt, steps))
end

function _oblique_history_errors(history)
    initial = history.samples[:, 1]
    scale = norm(initial)
    scale > 0 || throw(ArgumentError("measured initial oblique mode is zero"))
    error = leakage = 0.0
    for j in eachindex(history.times)
        measured = @view history.samples[:, j]
        expected = initial * exp(history.reference.rate * history.times[j])
        error = max(error, norm(measured - expected) / scale)
        coefficient = dot(history.reference.vector, measured)
        leakage =
            max(leakage, norm(measured - coefficient * history.reference.vector) / norm(measured))
    end
    return (; history_error = error, leakage)
end
