module ModeHistoryTests
using Test, LinearAlgebra
include(joinpath(@__DIR__, "..", "validation", "mode_history.jl"))

@testset "Signed frequency and rate from independent exponential histories" begin
    for T in (Float32, Float64),
        n in (9, 16, 257),
        omega in (-7.0, 0.0, 3.5),
        gamma in (-0.2, 0.0, 0.15),
        angle in (-2.0, 0.7)

        dt = 0.05
        for scale in (T(1e-20), one(T), T(1e20))
            data = Complex{T}[
                scale * exp(gamma * j * dt) * cis(angle - omega * j * dt) for j = 0:(n-1)
            ]
            before = copy(data)
            fit = _fit_mode_history(data, dt)
            tolerance = 512eps(T)
            @test isapprox(fit.omega, omega; atol = tolerance, rtol = tolerance)
            @test isapprox(fit.gamma, gamma; atol = tolerance, rtol = tolerance)
            @test fit.phase_rms < tolerance
            @test fit.logamp_rms < tolerance
            @test fit.omega_window_difference < tolerance
            @test fit.gamma_window_difference < tolerance
            @test data == before
            @test fit.samples == n
        end
    end
end

@testset "Mode fit agrees with independent dense least squares" begin
    n, dt = 129, 0.125
    t = collect(0:(n-1)) .* dt
    # Fixed nonlinear contamination makes the two OLS residuals nonzero.
    phases = @. 0.3 - 1.7t + 0.02sin(0.4t)
    logs = @. -0.8 + 0.03t + 0.01cos(0.7t)
    signal = exp.(logs) .* cis.(phases)
    fit = _fit_mode_history(signal, dt)
    design = hcat(ones(n), t)
    for (values, rate, rms, stderr) in (
        (phases, -fit.omega, fit.phase_rms, fit.omega_stderr),
        (logs, fit.gamma, fit.logamp_rms, fit.gamma_stderr),
    )
        beta = design \ values
        residual = values - design * beta
        covariance = (sum(abs2, residual) / (n - 2)) * inv(design' * design)
        @test isapprox(rate, beta[2]; rtol = 1e-12)
        @test isapprox(rms, norm(residual) / sqrt(n); rtol = 1e-10)
        @test isapprox(stderr, sqrt(covariance[2, 2]); rtol = 1e-10)
    end
    part = @view signal[17:112]
    a = _fit_mode_history(part, dt)
    b = _fit_mode_history(copy(part), dt)
    @test a == b
    reversefit = _fit_mode_history(@view(signal[end:-1:1]), dt)
    @test isapprox(reversefit.omega, -fit.omega; rtol = 1e-12)
    @test isapprox(reversefit.gamma, -fit.gamma; rtol = 1e-12)
    conjugatefit = _fit_mode_history(conj.(signal), dt)
    @test conjugatefit.omega == -fit.omega
    @test conjugatefit.gamma == fit.gamma
end

@testset "Mode fit reports window sensitivity and spectral contamination" begin
    dt = 0.025
    t = collect(0:511) .* dt
    chirp = cis.(@. -2t - 0.2t^2)
    changing_rate = @. exp(-0.05t - 0.01t^2) * cis(-2t)
    two_modes = @. cis(-2t) + 0.3cis(-3.1t)
    @test _fit_mode_history(chirp, dt).omega_window_difference > 1
    @test _fit_mode_history(chirp, dt).phase_rms > 0.1
    @test _fit_mode_history(changing_rate, dt).gamma_window_difference > 0.1
    @test _fit_mode_history(two_modes, dt).phase_rms > 0.1
    @test _fit_mode_history(two_modes, dt).logamp_rms > 0.1
end

@testset "Long history phase accumulation" begin
    for T in (Float32, Float64), reversed in (false, true)
        n = 10000
        signal = Complex{T}[exp(-0.0001j) * cis(-0.1j) for j = 0:(n-1)]
        data = reversed ? view(signal, n:-1:1) : signal
        fit = _fit_mode_history(data, T(0.1))
        @test isapprox(fit.omega, reversed ? -1 : 1; rtol = 128eps(T))
        @test isapprox(fit.gamma, reversed ? 0.001 : -0.001; atol = 128eps(T))
    end
end

@testset "Mode fit representability and failure contract" begin
    for amplitude in (1e-280, 1e280)
        signal = [amplitude * cis(-0.1j) for j = 0:63]
        fit = _fit_mode_history(signal, 0.1)
        @test isapprox(fit.omega, 1; rtol = 1e-12)
        @test abs(fit.gamma) < 1e-12
    end
    for value in
        (complex(floatmax(Float64), floatmax(Float64)), complex(nextfloat(0.0), nextfloat(0.0)))
        fit = _fit_mode_history(fill(value, 16), 0.25)
        @test fit.omega == 0
        @test abs(fit.gamma) < 1e-12
        @test fit.periods == 0
    end
    valid = ComplexF64[cis(-0.1j) for j = 0:15]
    for n in (0, 1, 7)
        @test_throws ArgumentError _fit_mode_history(valid[1:n], 0.1)
    end
    for dt in (0.0, -0.1, NaN, Inf, floatmax(Float64))
        @test_throws ArgumentError _fit_mode_history(valid, dt)
    end
    for bad in (0.0 + 0im, Inf + 0im, NaN + 0im)
        signal = copy(valid)
        signal[8] = bad
        @test_throws ArgumentError _fit_mode_history(signal, 0.1)
        @test isequal(signal[8], bad)
    end
    @test_throws ArgumentError _fit_mode_history(repeat(ComplexF64[1, im, -1, -im], 4), 1)
    @test_throws ArgumentError _fit_mode_history([cis(-2j) for j = 0:15], 1)
end
end # module
