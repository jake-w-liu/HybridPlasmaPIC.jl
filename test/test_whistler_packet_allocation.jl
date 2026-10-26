module TestWhistlerPacketAllocation
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC

bits(x::Float64) = reinterpret(UInt64, x)
bits(x::ComplexF64) = (bits(real(x)), bits(imag(x)))
function same_packet(actual, initial, evolved)
    energy = sum(abs2, evolved)
    @test bits(actual.energy0) == bits(initial.energy0)
    @test bits(actual.energyN) == bits(energy)
    @test bits(actual.energy_ratio) == bits(energy / initial.energy0)
    @test bits.(actual.bhat) == bits.(evolved)
    @test bits.(actual.ω) == bits.(initial.ω)
end

function check_history(; steps = 5, method = :cn, dt = 0.05, kwargs...)
    initial = run_whistler(; method, dt, nsteps = 0, kwargs...)
    evolved = copy(initial.bhat)
    factor = method === :cn ? cn_multiplier : H.euler_multiplier
    for step = 0:steps
        if step > 0
            # Recompute in the original step-major order. Powers of a factor
            # or deferred Float64 stores do not preserve this rounding history.
            for m in eachindex(evolved)
                evolved[m] *= factor(initial.ω[m], dt)
            end
        end
        actual = run_whistler(; method, dt, nsteps = step, kwargs...)
        same_packet(actual, initial, evolved)
    end
end

@testset "Whistler packet ordered histories" begin
    for method in (:cn, :euler),
        T in (Float16, Float32, Float64),
        dt in (T(0), -T(0), nextfloat(T(0)), T(0.05), -T(0.3), floatmax(T)),
        (n, band) in ((1, nothing), (4, nothing), (7, 1), (8, 0))

        check_history(; method, dt, n, band, seed = 4)
    end
    cases = (
        (floatmax(Float64), floatmin(Float64), floatmax(Float64)),
        (floatmin(Float64), floatmax(Float64), floatmin(Float64)),
        (floatmin(Float64), 1.0, sqrt(floatmin(Float64))),
        (floatmax(Float64), 1.0, floatmax(Float64)),
        (0.0, 1.0, floatmin(Float64)),
        (big(1) // big(3), big(7) // big(11), Float64(2π)),
    )
    for (B0, n0, L) in cases, method in (:cn, :euler), sign in (-1, 1)
        check_history(; n = 4, method, B0 = sign * B0, n0, L, dt = floatmax(Float64))
    end
    for I in (Int8, UInt8, Int64, UInt64, Int128, UInt128, BigInt), step in (0, 1, 5)
        @test isequal(
            run_whistler(; n = I(7), nsteps = I(step), band = I(1)),
            run_whistler(; n = 7, nsteps = step, band = 1),
        )
    end
    for count in (false, true)
        @test isequal(
            run_whistler(; n = 4, nsteps = count),
            run_whistler(; n = 4, nsteps = Int(count)),
        )
    end
end

struct ObservedStep <: Real
    value::Float64
    calls::Vector{Float64}
end
Base.isfinite(dt::ObservedStep) = isfinite(dt.value)
function observed_factor(f, omega, dt)
    push!(dt.calls, omega)
    return f(omega, dt.value + length(dt.calls) * eps(dt.value))
end
H.cn_multiplier(omega::Float64, dt::ObservedStep) = observed_factor(cn_multiplier, omega, dt)
H.euler_multiplier(omega::Float64, dt::ObservedStep) =
    observed_factor(H.euler_multiplier, omega, dt)

@testset "Generic timestep evaluation and settings compatibility" begin
    for method in (:cn, :euler)
        for dt in (0, 1, 1 // 3, big(1), big(1) // big(3), π)
            check_history(; method, dt, n = 4, steps = 3)
        end
        for p in (24, 128), mode in (RoundNearest, RoundToZero, RoundUp, RoundDown, RoundFromZero)
            dt = BigFloat(1 // 3; precision = 113)
            stored = (precision(dt), Base.decompose(dt))
            setprecision(BigFloat, p) do
                setrounding(BigFloat, mode) do
                    check_history(; method, dt, n = 4, steps = 3)
                    @test precision(BigFloat) == p
                    @test rounding(BigFloat) === mode
                    @test (precision(dt), Base.decompose(dt)) == stored
                end
            end
        end
        dt = ObservedStep(0.2, Float64[])
        initial = run_whistler(; method, dt, n = 4, nsteps = 0)
        @test isempty(dt.calls)
        expected = copy(initial.bhat)
        factor = method === :cn ? cn_multiplier : H.euler_multiplier
        for _ = 1:3, m in eachindex(expected)
            expected[m] *= factor(initial.ω[m], dt)
        end
        calls = copy(dt.calls)
        empty!(dt.calls)
        actual = run_whistler(; method, dt, n = 4, nsteps = 3)
        @test dt.calls == calls == repeat(initial.ω, 3)
        same_packet(actual, initial, expected)
    end
    for method in (:cn, :euler, :invalid), steps in (0, 1)
        @test_throws ArgumentError run_whistler(;
            method,
            nsteps = steps,
            n = 4,
            B0 = floatmax(Float64),
            L = floatmin(Float64),
        )
    end
    @test_throws ArgumentError run_whistler(; method = :invalid, nsteps = 0)
    @test_throws ArgumentError run_whistler(; dt = Inf, nsteps = 0)
end

packet_bytes(n, nsteps, dt, B0, n0, L) = @allocated run_whistler(; n, nsteps, dt, B0, n0, L)

@testset "Hardware packet allocation is independent of step count" begin
    for (B0, n0, L, dt) in (
        (1.0, 1.0, Float64(2π), 1e-6),
        (1, 1, Float64(2π), 1e-6),
        (1.0, 1.0, Float64(2π), floatmax(Float64)),
        (floatmin(Float64), 1.0, sqrt(floatmin(Float64)), floatmax(Float64)),
        (big(1) // big(3), big(7) // big(11), Float64(2π), floatmax(Float64)),
    )
        for _ = 1:2
            packet_bytes(4, 1, dt, B0, n0, L)
            packet_bytes(4, 100, dt, B0, n0, L)
        end
        # Coverage instrumentation (CI) can add a few hundred fixed bytes per call;
        # the regression guarded here allocates per step and mode (megabytes).
        extra = packet_bytes(4, 100, dt, B0, n0, L) - packet_bytes(4, 1, dt, B0, n0, L)
        @test abs(extra) <= 4096
    end
end

@testset "BigInt step iteration cost is independent of mode count" begin
    for steps in (0, 1, 100)
        count = BigInt(steps)
        extras = Int[]
        for n in (8, 128)
            for _ = 1:2
                packet_bytes(n, count, 1e-6, 1.0, 1.0, Float64(2π))
                packet_bytes(n, steps, 1e-6, 1.0, 1.0, Float64(2π))
            end
            big_bytes = packet_bytes(n, count, 1e-6, 1.0, 1.0, Float64(2π))
            fixed_bytes = packet_bytes(n, steps, 1e-6, 1.0, 1.0, Float64(2π))
            push!(extras, big_bytes - fixed_bytes)
        end
        # Same tolerance: a per-mode BigInt iteration cost would differ by megabytes.
        @test abs(extras[1] - extras[2]) <= 4096
    end
end
end
