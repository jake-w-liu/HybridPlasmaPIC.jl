module LocalCoulombCompletionCancellationTests
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
Q(x) = Rational{BigInt}(x)
function independent(va, vb, ma, mb, speed, sample)
    T = typeof(ma)
    a, b = Q(ma), Q(mb)
    g = Q.(va) .- Q.(vb)
    radius = Q(hypot(va[1] - vb[1], va[2] - vb[2]))
    mag = Q(speed)
    if sample.isotropic
        cosine, sine = Q(sample.cosine), Q(sample.sine)
    else
        d = Q(sample.tangent) * big(2)^sample.tangent_power
        cosine, sine = (1 - d * d) / (1 + d * d), 2d / (1 + d * d)
    end
    cp, sp = Q(sample.cosφ), Q(sample.sinφ)
    rotated =
        iszero(radius) ? (mag * cp, mag * sp, zero(mag)) :
        (
            (g[1] * g[3] * cp - g[2] * mag * sp) / radius,
            (g[2] * g[3] * cp + g[1] * mag * sp) / radius,
            -radius * cp,
        )
    after = g .* cosine .+ rotated .* sine
    center = (a .* Q.(va) .+ b .* Q.(vb)) ./ (a + b)
    map(x -> H._from_exact(T, x), center .+ b / (a + b) .* after),
    map(x -> H._from_exact(T, x), center .- a / (a + b) .* after)
end
counts = Dict{String,Int}()
failures = Any[]
@testset "Structured complete-angle cancellation and certificate boundaries" begin
    for T in (Float32, Float64),
        scale_power in (-100, 0, 100),
        k in (-100, -52, -30, -12, 0),
        mass_case = 1:3,
        center_case = 0:2,
        axis = 0:2,
        angle_case = 1:5

        unit = ldexp(one(T), scale_power)
        small = ldexp(one(T), k)
        g = circshift([unit, unit * small, unit], axis)
        center = center_case == 0 ? zero(T) : (center_case == 1 ? unit : -unit)
        va = (center, center / 2, center / 4)
        vb = Tuple(va .- g)
        rel = va .- vb
        speed = hypot(rel...)
        ma, mb =
            mass_case == 1 ? (one(T), T(3)) :
            mass_case == 2 ? (T(3), one(T)) : (ldexp(one(T), -100), ldexp(one(T), 100))
        tangent =
            angle_case == 1 ? prevfloat(one(T)) :
            angle_case == 2 ? one(T) :
            angle_case == 3 ? nextfloat(one(T)) :
            angle_case == 4 ? ldexp(one(T), -100) : ldexp(one(T), 100)
        sp, cp = sincos(small)
        sample = (;
            tangent,
            tangent_power = 0,
            isotropic = false,
            cosine = zero(T),
            sine = one(T),
            cosφ = cp,
            sinφ = sp,
        )
        adapter =
            T === Float32 ? H._coulomb_widened_update(va, vb, ma, mb, speed, sample) :
            H._coulomb_compensated_update(va, vb, ma, mb, speed, sample)
        key = string(T) * (adapter === nothing ? "_exact" : "_hardware")
        counts[key] = get(counts, key, 0) + 1
        actual = H._coulomb_complete_update(va, vb, ma, mb, speed, sample)
        expected = independent(va, vb, ma, mb, speed, sample)
        for s = 1:2, c = 1:3
            pass =
                isequal(actual[s][c], expected[s][c]) ||
                isapprox(actual[s][c], expected[s][c]; rtol = 8eps(T), atol = nextfloat(zero(T)))
            @test pass
            if !pass
                push!(
                    failures,
                    (;
                        T = string(T),
                        scale_power,
                        k,
                        mass_case,
                        center_case,
                        axis,
                        angle_case,
                        s,
                        c,
                        actual = Float64(actual[s][c]),
                        expected = Float64(expected[s][c]),
                        hardware = adapter !== nothing,
                    ),
                )
            end
        end
    end
    for T in (Float32, Float64)
        @test get(counts, string(T) * "_exact", 0) > 0
        @test get(counts, string(T) * "_hardware", 0) > 0
    end
end

struct FixedAngle <: AbstractRNG
    phase::Float32
end
Random.rand(r::FixedAngle, ::Type{Float32}) = r.phase
Random.randn(::FixedAngle, ::Type{Float32}) = 1.0f0
@testset "Widening must retain cancellation inside the rotated increment" begin
    for k in (-30, -40, -50), sign in (-1.0f0, 1.0f0), (ma, mb) in ((1.0f0, 3.0f0), (3.0f0, 1.0f0))
        tiny = sign * ldexp(1.0f0, k)
        va = (0.0f0, 0.0f0, 0.0f0)
        vb = (-1.0f0, -tiny, -1.0f0)
        speed = hypot((va .- vb)...)
        phase = tiny / Float32(2pi)
        rng = FixedAngle(phase)
        # gcdt and denominator use identical stored products, so the sampled tangent is exactly one.
        base = speed * speed * speed
        sample = H._coulomb_sample_increment(
            1.0f0,
            tiny,
            1.0f0,
            speed,
            base,
            1.0f0,
            1.0f0,
            0.001f0,
            Float32(2pi),
            rng,
            (speed, speed, speed),
        )
        expected = independent(va, vb, ma, mb, speed, sample)
        actual = H._coulomb_complete_update(va, vb, ma, mb, speed, sample)
        @test sample.complete
        @test sample.tangent == 1.0f0
        @test !iszero(expected[1][1])
        @test isapprox(actual[1][1], expected[1][1]; rtol = 8eps(Float32), atol = nextfloat(0.0f0))
    end
end

end
