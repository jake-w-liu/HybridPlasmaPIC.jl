module TestLinearCNRange
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
dyadic(e) = e >= 0 ? (big(1) << e) // big(1) : big(1) // (big(1) << -e)
function reference(a, b)
    h = Q(a) * Q(b) / 2
    return ((1 - h^2) / (1 + h^2), -2h / (1 + h^2))
end

@testset "Bounded exact dyadic comparisons and conversion" begin
    rng = MersenneTwister(5701)
    for index = 1:2000
        a, b, c = BigInt.(rand(rng, UInt128, 3))
        ea, eb, ec = Int128.(rand(rng, -2000:2000, 3))
        @test sign(H._compare_scaled(a, ea, b, eb)) == sign(cmp(a * dyadic(ea), b * dyadic(eb)))
        @test sign(H._compare_scaled_sum(a, ea, b, eb, c, ec)) ==
              sign(cmp(a * dyadic(ea), b * dyadic(eb) + c * dyadic(ec)))
        if index <= 200
            ratio = (a + 1) // (b + 1)
            for T in (Float16, Float32, Float64), sign in (-1, 1)
                @test isequal(
                    H._from_exact(T, sign * ratio, ea),
                    H._from_exact(T, sign * ratio * dyadic(ea)),
                )
            end
        end
    end
    huge = Int128(typemax(Int64)) * 4
    for base in (-huge, huge)
        @test H._compare_scaled(big(3), base, big(6), base - 1) == 0
        @test H._compare_scaled_sum(big(3), base, big(1), base + 1, big(1), base) == 0
        @test H._compare_scaled_sum(big(3), base, big(1), base + 1, big(1), -huge - 100) > 0
        @test H._compare_scaled_sum(big(1), base, big(1), base, big(1), -huge - 100) < 0
    end
    @test H._compare_scaled(big(0), huge, big(0), -huge) == 0
    @test H._compare_scaled_sum(big(0), huge, big(0), huge, big(1), -huge) < 0
end

@testset "Finite linear CN components and compatibility" begin
    rng = MersenneTwister(7215)
    for T in (Float16, Float32, Float64)
        values = T[
            0,
            nextfloat(zero(T)),
            floatmin(T),
            sqrt(floatmin(T)),
            0.125,
            prevfloat(one(T)),
            1,
            nextfloat(one(T)),
            prevfloat(T(2)),
            2,
            nextfloat(T(2)),
            4,
            sqrt(floatmax(T)),
            floatmax(T),
        ]
        for a in values, b in values, sign in (-1, 1)
            expected = T.(reference(sign * a, b))
            actual = cn_multiplier(sign * a, b)
            @test typeof(actual) === Complex{T}
            @test isfinite(actual)
            @test real(actual) ≈ expected[1] rtol = 8eps(T) atol = 0
            @test imag(actual) ≈ expected[2] rtol = 8eps(T) atol = 0
        end
        for _ = 1:1000
            a, b = T.(exp2.(rand(rng, 2) .* 20 .- 10))
            rand(rng, Bool) && (a = -a)
            actual = cn_multiplier(a, b)
            re, imaginary = T.(reference(a, b))
            @test real(actual) ≈ re rtol = 8eps(T) atol = 0
            @test imag(actual) ≈ imaginary rtol = 8eps(T) atol = 0
        end
        for a in (zero(T), -zero(T)), b in (one(T), -one(T))
            @test isequal(cn_multiplier(a, b), complex(one(T), -a * b))
            @test isequal(cn_multiplier(b, a), complex(one(T), -a * b))
        end
        for a in (T(Inf), T(-Inf), T(NaN)), b in (zero(T), one(T))
            @test isequal(cn_multiplier(a, b), (1 - im * a * b / 2) / (1 + im * a * b / 2))
        end
    end
    for I in (Int8, UInt8, Int32, UInt32, Int64, UInt64, Int128, UInt128)
        for (a, b) in ((typemax(I), I(2)), (I(2), typemax(I)))
            expected = complex(Float64.(reference(a, b))...)
            actual = cn_multiplier(a, b)
            @test actual isa ComplexF64
            @test real(actual) ≈ real(expected) rtol = 8eps(Float64) atol = 0
            @test imag(actual) ≈ imag(expected) rtol = 8eps(Float64) atol = 0
        end
    end
    for (a, b) in (
        (-2, UInt128(3)),
        (UInt128(3), -2),
        (3 // 5, 0.25f0),
        (0.25f0, 3 // 5),
        (big(3), big(2)),
        (typemax(UInt128), -2),
    )
        actual = cn_multiplier(a, b)
        T = typeof(float(zero(promote_type(typeof(a), typeof(b)))))
        @test typeof(actual) === Complex{T}
        expected = complex(map(x -> H._from_exact(T, x), reference(a, b))...)
        @test real(actual) ≈ real(expected) rtol = 8eps(T) atol = 0
        @test imag(actual) ≈ imag(expected) rtol = 8eps(T) atol = 0
    end
    for I in (Int8, UInt8, Int64, UInt64, Int128, UInt128, BigInt)
        a, b = I(1) // I(2), I(3) // I(4)
        @test cn_multiplier(a, b) == complex(reference(a, b)...)
        if I !== BigInt
            for (a, b) in ((typemax(I) // I(1), I(2) // I(1)), (I(1) // typemax(I), I(2) // I(1)))
                @test cn_multiplier(a, b) == complex(reference(a, b)...)
            end
        end
    end
    @test cn_multiplier(3 // 2, 1 // 8) isa Complex{Rational{Int}}
    @test cn_multiplier(typemax(Int) // 1, 2 // 1) isa Complex{Rational{BigInt}}
end

# The unchanged ordinary arithmetic remains byte-for-byte compatible, while
# corrected cancellation/range regions above use an independent component oracle.
scalar_bytes(a, b) = @allocated cn_multiplier(a, b)
legacy_scalar(a, b) = (1 - im * a * b / 2) / (1 + im * a * b / 2)
legacy_scalar_bytes(a, b) = @allocated legacy_scalar(a, b)
packet_bytes(n, B0, n0, L) = @allocated run_whistler(; n, B0, n0, L, nsteps = 1, dt = 1e-6)
frequency_bytes(n, L, c) = @allocated H._whistler_omega(n, L, c)
@testset "Ordinary CN compatibility and resources" begin
    for T in (Float16, Float32, Float64), (x, y) in ((0.3, 0.05), (40, 1.3), (-0.7, 0.1), (3, 2))
        a, b = T(x), T(y)
        expected = (1 - im * a * b / 2) / (1 + im * a * b / 2)
        @test isequal(cn_multiplier(a, b), expected)
        scalar_bytes(a, b)
        @test scalar_bytes(a, b) == 0
    end
    for (a, b) in ((3.0, 2), (2, 3.0), (3.0f0, 2), (2, 3.0f0))
        scalar_bytes(a, b)
        @test scalar_bytes(a, b) == 0
    end
    for T in (Float16, Float32, Float64),
        (a, b) in ((one(T), T(2)), (nextfloat(one(T)), prevfloat(T(2))))

        scalar_bytes(a, b)
        @test scalar_bytes(a, b) == 0
    end
    for p in (24, 53, 128, 256)
        setprecision(BigFloat, p) do
            for (x, y) in ((3, 0.125), (3, 2), (1.5, 0.75), (0.5, 0.125))
                a, b = BigFloat(x), BigFloat(y)
                scalar_bytes(a, b)
                legacy_scalar_bytes(a, b)
                @test scalar_bytes(a, b) <= legacy_scalar_bytes(a, b)
                expected = BigFloat.(reference(a, b))
                @test isequal(cn_multiplier(a, b), complex(expected...))
            end
            for a in (BigFloat(0), BigFloat(-0.0)), b in (BigFloat(1), BigFloat(-1))
                @test isequal(cn_multiplier(a, b), complex(one(BigFloat), -a * b))
            end
            for a in (BigFloat(Inf), BigFloat(-Inf), BigFloat(NaN)), b in (BigFloat(0), BigFloat(1))
                @test isequal(cn_multiplier(a, b), legacy_scalar(a, b))
            end
        end
    end
    for n in (8, 128), L in (2π, 8, 8.0)
        packet_bytes(n, 1.0, 1.0, Float64(L))
        ceiling = packet_bytes(n, 1.0, 1.0, Float64(L))
        for (B0, n0) in ((1, 1), (1, 1.0), (1.0, 1), (UInt128(1), Int8(1)))
            packet_bytes(n, B0, n0, L)
            @test packet_bytes(n, B0, n0, L) == ceiling
            @test run_whistler(; n, B0, n0, L, nsteps = 1) ==
                  run_whistler(; n, B0 = 1.0, n0 = 1.0, L = Float64(L), nsteps = 1)
        end
        frequency_bytes(n, L, 1.0)
        frequency_bytes(n, L, 1)
        @test frequency_bytes(n, L, 1) == frequency_bytes(n, L, 1.0)
        @test H._whistler_omega(n, L, 1) == H._whistler_omega(n, L, 1.0, 1.0)
    end
end

@testset "Complete whistler frequencies and trajectories" begin
    cases = (
        (1.0, 1.0, Float64(2π)),
        (Float16(1), Float16(3), Float64(2π)),
        (Float32(1), Float32(3), Float64(2π)),
        (floatmax(Float64), floatmin(Float64), floatmax(Float64)),
        (floatmin(Float64), floatmax(Float64), floatmin(Float64)),
        (floatmin(Float64), 1.0, sqrt(floatmin(Float64))),
        (floatmax(Float64), 1.0, floatmax(Float64)),
        (0.0, 1.0, floatmin(Float64)),
        (big(1) // big(3), big(7) // big(11), 2π),
        (typemax(Int64), typemax(UInt64), typemax(Int32)),
        (typemax(UInt128), 1, UInt128(1) << 64),
        (1, typemax(UInt128), 2π),
    )
    for (B0, n0, L) in cases, n in (1, 4, 7, 8), sign in (-1, 1)
        expected = [
            Float64(Q(sign * B0) / Q(n0) * (Q(Float64(2π)) * (m <= n ÷ 2 ? m : m - n) / Q(L))^2) for m = 0:n-1
        ]
        packet = run_whistler(; n, L, B0 = sign * B0, n0, dt = floatmax(Float64), nsteps = 3)
        @test packet.ω ≈ expected rtol = 16eps(Float64) atol = 0
        @test all(isfinite, packet.bhat)
        @test packet.energy_ratio ≈ 1 rtol = 128eps(Float64)
        initial = run_whistler(; n, L, B0 = sign * B0, n0, nsteps = 0)
        evolved = initial.bhat .* cn_multiplier.(packet.ω, floatmax(Float64)) .^ 3
        @test packet.bhat ≈ evolved rtol = 32eps(Float64)
    end
    for I in (Int8, UInt8, Int64, UInt64, Int128, UInt128, BigInt), n in (7, 8), band in (0, 1, 2)
        ordinary = run_whistler(; n, band, nsteps = 3)
        typed = run_whistler(; n = I(n), band = I(band), nsteps = 3)
        @test typed == ordinary
    end
    rng = MersenneTwister(4189)
    for T in (Float16, Float32, Float64), _ = 1:40
        limit = min(600, exponent(floatmax(T)) - 2)
        B0, n0, L = [ldexp(T(1) + rand(rng, T), rand(rng, -limit:limit)) for _ = 1:3]
        rand(rng, Bool) && (B0 = -B0)
        expected = [Float64(Q(B0) / Q(n0) * (Q(Float64(2π)) * m / Q(L))^2) for m in (0, 1, 2, -1)]
        if all(isfinite, expected)
            @test H._whistler_omega(4, L, B0, n0) ≈ expected rtol = 16eps(Float64) atol = 0
        else
            @test_throws ArgumentError H._whistler_omega(4, L, B0, n0)
        end
    end
    @test H._whistler_omega(8, 2π, 1.0)[5] == 16.0
    @test_throws ArgumentError run_whistler(; n = 4, B0 = floatmax(Float64), L = floatmin(Float64))
    @test_throws ArgumentError compare_integrators_whistler(;
        n = 4,
        B0 = floatmax(Float64),
        L = floatmin(Float64),
    )
    comparison = compare_integrators_whistler(; n = 4, nsteps = 1, dt_stiff = floatmax(Float64))
    @test comparison.cn_ratio_stiff ≈ 1 rtol = 64eps(Float64)
end

@testset "BigFloat CN settings and exponent limits in an isolated process" begin
    project = Base.active_project()
    @test project !== nothing
    fixture = joinpath(@__DIR__, "linear_cn_bigfloat_fixture.jl")
    package = dirname(@__DIR__)
    command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(project)) $fixture $package`
    log = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = log, stderr = log))
    transcript = String(take!(log))
    success(process) || @error "BigFloat CN fixture failed" transcript
    @test success(process)
    @test occursin("BigFloat CN settings restored", transcript)
end
end
