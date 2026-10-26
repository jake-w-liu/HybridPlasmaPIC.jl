module TestWhistlerBigIntCount
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC

@testset "Whistler exact fixed-counter admission" begin
    for count in BigInt.((0, 1, 2, 100, typemax(Int) - 1, typemax(Int)))
        actual = H._whistler_fixed_count(count)
        if count > 0
            @test actual isa Int
            @test BigInt(actual) == count
        else
            @test actual === count
        end
    end
    for count in (BigInt(-1), -(big(1) << 128), BigInt(typemax(Int)) + 1, big(1) << 128)
        # Test the admission boundary without executing astronomical histories.
        @test H._whistler_fixed_count(count) === count
    end
    for I in (Int8, UInt8, Int, UInt, Int128, UInt128, Bool), count in (0, 1)
        value = I(count)
        @test H._whistler_fixed_count(value) === value
    end
end

@testset "BigInt packet ordered histories" begin
    for method in (:cn, :euler),
        T in (Float16, Float32, Float64),
        dt in (T(0), -T(0), nextfloat(T(0)), T(0.05), -T(0.3), floatmax(T)),
        (n, band) in ((1, nothing), (7, 1), (8, 0))

        initial = run_whistler(; method, dt, n, band, nsteps = BigInt(0), seed = 9)
        expected = copy(initial.bhat)
        factor = method === :cn ? cn_multiplier : H.euler_multiplier
        for steps = 0:12
            if steps > 0
                for m in eachindex(expected)
                    expected[m] *= factor(initial.ω[m], dt)
                end
            end
            actual = run_whistler(; method, dt, n, band, nsteps = BigInt(steps), seed = 9)
            @test reinterpret(UInt64, actual.bhat) == reinterpret(UInt64, expected)
            @test reinterpret(UInt64, actual.ω) == reinterpret(UInt64, initial.ω)
            energy = sum(abs2, expected)
            @test reinterpret(UInt64, actual.energy0) == reinterpret(UInt64, initial.energy0)
            @test reinterpret(UInt64, actual.energyN) == reinterpret(UInt64, energy)
            @test reinterpret(UInt64, actual.energy_ratio) ==
                  reinterpret(UInt64, energy / initial.energy0)
        end
    end
end

struct ObservedStep <: Real
    calls::Vector{Float64}
end
Base.isfinite(::ObservedStep) = true
function factor(f, omega, dt)
    push!(dt.calls, omega)
    return f(omega, 0.2 + length(dt.calls) * eps(0.2))
end
H.cn_multiplier(omega::Float64, dt::ObservedStep) = factor(cn_multiplier, omega, dt)
H.euler_multiplier(omega::Float64, dt::ObservedStep) = factor(H.euler_multiplier, omega, dt)

@testset "BigInt counts preserve generic timestep evaluation" begin
    for method in (:cn, :euler)
        dt = ObservedStep(Float64[])
        initial = run_whistler(; method, dt, n = 7, nsteps = BigInt(0))
        @test isempty(dt.calls)
        expected = copy(initial.bhat)
        multiplier = method === :cn ? cn_multiplier : H.euler_multiplier
        for _ = 1:5, m in eachindex(expected)
            expected[m] *= multiplier(initial.ω[m], dt)
        end
        expected_calls = copy(dt.calls)
        empty!(dt.calls)
        actual = run_whistler(; method, dt, n = 7, nsteps = BigInt(5))
        @test dt.calls == expected_calls == repeat(initial.ω, 5)
        @test reinterpret(UInt64, actual.bhat) == reinterpret(UInt64, expected)
    end
end

packet_bytes(n, nsteps, dt, method) = @allocated run_whistler(; n, nsteps, dt, method)

@testset "Bounded BigInt packet count allocation" begin
    for method in (:cn, :euler),
        T in (Float16, Float32, Float64),
        n in (8, 128),
        steps in (1, 2, 100)

        count, dt = BigInt(steps), T(0.05)
        for _ = 1:2
            packet_bytes(n, count, dt, method)
            packet_bytes(n, steps, dt, method)
        end
        for _ = 1:3
            big_bytes = packet_bytes(n, count, dt, method)
            fixed_bytes = packet_bytes(n, steps, dt, method)
            # Fixed per-call bookkeeping is allowed; per-step BigInt iteration
            # is not. Retain every declared sample and the old 4096-byte gate.
            @test abs(big_bytes - fixed_bytes) <= 512
        end
    end
end
end
