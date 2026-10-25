module ValidationFrequencyTests
using Test
include(joinpath(@__DIR__, "..", "validation", "common.jl"))

@testset "Complex frequency bin, sign and offset invariants" begin
    for T in (Float32, Float64), n in (5, 6, 7, 8, 17, 32, 63, 64)
        dt = T(0.25)
        for bin = 1:(n-1), phase in (0.0, 0.37), dc in (0.0 + 0.0im, 10.0 + 20.0im)
            signal = Complex{T}[dc + 0.25cispi(2bin * j / n + phase) for j = 0:(n-1)]
            original = copy(signal)
            # Sampling aliases bin n-k to -k; the physical frequency is unsigned.
            expected = 2pi * min(bin, n - bin) / (n * dt)
            @test isapprox(_peak_frequency(signal, dt), expected; rtol = 128eps(T))
            @test signal == original
        end
    end
end

@testset "Frequency evidence and interval failures" begin
    for n in (0, 1, 4, 5, 7, 64), value in (0.0 + 0.0im, 1.0 + 2.0im)
        @test isnan(_peak_frequency(fill(value, n), 0.25))
    end
    for value in (NaN + 0im, Inf + 0im, 0.0 - Inf * im)
        signal = ComplexF64[cispi(j / 4) for j = 0:31]
        signal[7] = value
        @test isnan(_peak_frequency(signal, 0.25))
    end
    for dt in (0.0, -0.25, Inf, -Inf, NaN)
        @test_throws ArgumentError _peak_frequency(ones(ComplexF64, 16), dt)
    end
    result = _result(
        id = "flat",
        category = "em_pic",
        reference_kind = "analytic",
        reference = "synthetic constant",
        metric = "frequency",
        measured = NaN,
        expected = 1.0,
        error_kind = "relative",
        error = NaN,
        tolerance = 0.03,
    )
    @test result.status == "fail"
end

@testset "Frequency amplitude range and independent histories" begin
    for T in (Float32, Float64), amplitude in (nextfloat(zero(T)), T(0.125), floatmax(T))
        # Exact quarter-cycle samples avoid contaminating subnormal tests with
        # trigonometric rounding. Neither the period nor amplitude is fitted.
        signal = repeat(Complex{T}[amplitude, im*amplitude, -amplitude, -im*amplitude], 16)
        @test isapprox(_peak_frequency(signal, T(0.25)), 2pi; rtol = 8eps(T))
        @test isapprox(_peak_frequency(conj.(signal), T(0.25)), 2pi; rtol = 8eps(T))
    end
    for T in (Float32, Float64)
        # A large constant component must not erase a separately represented
        # small quadrature oscillation when the history is centered/scaled.
        signal = repeat(
            Complex{T}[
                complex(floatmax(T), y) for
                y in (nextfloat(zero(T)), zero(T), -nextfloat(zero(T)), zero(T))
            ],
            16,
        )
        @test isapprox(_peak_frequency(signal, T(0.25)), 2pi; rtol = 8eps(T))
    end
    n = 256
    dt = 0.2
    for bin in (3, 17, 52), phase in (0.0, 0.4)
        # A real standing wave exercises the two conjugate peaks used by case 20.
        standing = ComplexF64[cospi(2bin * j / n + phase) for j = 0:(n-1)]
        @test isapprox(_peak_frequency(standing, dt), 2pi * bin / (n * dt); rtol = 1e-12)
    end
    for bin in (3.2, 7.35, 17.2), direction in (-1, 1)
        signal = ComplexF64[cispi(direction * 2bin * j / n) for j = 0:(n-1)]
        # Rectangular-window parabolic interpolation is biased off-bin. Its
        # absolute frequency error must remain below half a Fourier bin.
        @test abs(_peak_frequency(signal, dt) - 2pi * bin / (n * dt)) < pi / (n * dt)
    end
    strong = ComplexF64[3cispi(-14j / n) + cispi(38j / n) for j = 0:(n-1)]
    @test isapprox(_peak_frequency(strong, dt), 2pi * 7 / (n * dt); rtol = 1e-12)
    @test isapprox(_peak_frequency(@view(strong[end:-1:1]), dt), 2pi * 7 / (n * dt); rtol = 1e-12)
end
end # module
