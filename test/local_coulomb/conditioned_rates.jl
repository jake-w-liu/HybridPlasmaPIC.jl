module LocalCoulombConditionedRateTests
using HybridPlasmaPIC, Test, Random
struct DirectedDraw{T} <: AbstractRNG
    normal::T
    phase::T
end
Random.rand(::DirectedDraw, r::UnitRange{Int}) = last(r)
Random.rand(r::DirectedDraw, T::Union{Type{Float32},Type{Float64}}) = T(r.phase)
Random.randn(r::DirectedDraw, T::Union{Type{Float32},Type{Float64}}) = T(r.normal)
rows = Any[]
function profile(T, masses, direction, power, sign, exchange, phase, vanishing)
    exponent = T === Float32 ? 100 : 900
    scale = ldexp(one(T), exponent)
    tiny = ldexp(one(T), -exponent)
    g = map(
        x -> scale * x,
        direction == 1 ? (zero(T), zero(T), one(T)) :
        direction == 2 ? (T(3), zero(T), T(4)) :
        direction == 3 ? (one(T), T(2), T(2)) :
        direction == 4 ? (T(-3), T(4), T(-12)) : (tiny, zero(T), one(T)),
    )
    ma, mb = T.(masses)
    q = vanishing ? one(T) : ldexp(one(T), div(exponent, 2))
    dt = vanishing ? one(T) : scale
    mu = ma * mb / (ma + mb)
    logarithm = T(4 * T(pi) * mu^2)
    draw = sign * ldexp(one(T), power)
    grid = FourierGrid((4,), (T(4),))
    a, b = ParticleSet{1,T}(1; q, m = ma), ParticleSet{1,T}(1; q, m = mb)
    for ps in (a, b)
        ps.x[1][1] = T(0.25)
        ps.weight[1] = one(T)
    end
    for c = 1:3
        a.v[c][1] = exchange ? (2mb / (ma + mb)) * g[c] : zero(T)
        b.v[c][1] = a.v[c][1] - g[c]
    end
    original = (ntuple(c -> a.v[c][1], 3), ntuple(c -> b.v[c][1], 3))
    physical = setprecision(BigFloat, 4096) do
        aa, bb = map(x -> BigFloat.(x), original)
        relative = aa .- bb
        speed = sqrt(sum(abs2, relative))
        reduced = BigFloat(ma) * BigFloat(mb) / (BigFloat(ma) + BigFloat(mb))
        variance =
            BigFloat(q)^4 * BigFloat(logarithm) * BigFloat(dt) /
            (8 * BigFloat(pi) * reduced^2 * speed^3)
        tangent = sqrt(variance) * BigFloat(draw)
        cosine_change = -2tangent^2 / (1 + tangent^2)
        sine = 2tangent / (1 + tangent^2)
        phi = BigFloat(T(2 * T(pi) * T(phase)))
        sp, cp = sincos(phi)
        radius = hypot(relative[1], relative[2])
        rotation =
            iszero(radius) ? (speed * cp, speed * sp, BigFloat(0)) :
            (
                (relative[1] * relative[3] * cp - relative[2] * speed * sp) / radius,
                (relative[2] * relative[3] * cp + relative[1] * speed * sp) / radius,
                -radius * cp,
            )
        change = relative .* cosine_change .+ rotation .* sine
        fraction = BigFloat(mb) / (BigFloat(ma) + BigFloat(mb))
        (T.(aa .+ fraction .* change), T.(bb .- (1 - fraction) .* change))
    end
    H = HybridPlasmaPIC
    relative = original[1] .- original[2]
    speed = hypot(relative...)
    channel = H._coulomb_channel_rate(a, b, grid, dt, one(T), one(T), logarithm)
    rate = H._CoulombLocalPairRate(channel, 1, 1, false)
    fraction, power = H._coulomb_local_scaled_coefficient(rate)
    sample = H._coulomb_sample_increment(
        relative...,
        speed,
        ldexp(fraction, power),
        one(T),
        one(T),
        one(T),
        2 * T(pi),
        DirectedDraw(draw, T(phase)),
        rate,
    )
    errors_in_inputs = setprecision(BigFloat, 512) do
        aa, bb = map(x -> BigFloat.(x), original)
        exact_speed = sqrt(sum(abs2, aa .- bb))
        reduced = BigFloat(ma) * BigFloat(mb) / (BigFloat(ma) + BigFloat(mb))
        variance =
            BigFloat(q)^4 * BigFloat(logarithm) * BigFloat(dt) /
            (8 * BigFloat(pi) * reduced^2 * exact_speed^3)
        exact_square = variance * BigFloat(draw)^2
        sampled_square = ldexp(BigFloat(sample.tangent), sample.tangent_power)^2
        exact_sin, exact_cos = sincos(BigFloat(T(2 * T(pi) * T(phase))))
        rate_error = abs((sampled_square - exact_square) / exact_square)
        speed_error = abs((BigFloat(speed) - exact_speed) / exact_speed)
        phase_error =
            max(abs(BigFloat(sample.sinφ) - exact_sin), abs(BigFloat(sample.cosφ) - exact_cos))
        @test rate_error <= 4096eps(T)
        @test speed_error <= 4eps(T)
        @test phase_error <= 2eps(T)
        Float64.((rate_error, speed_error, phase_error))
    end
    expected = setprecision(BigFloat, 4096) do
        # Independent center-of-mass reconstruction of the actual sampled law.
        # This separates integration error from conditioning of the physical rate.
        aa, bb = map(x -> BigFloat.(x), original)
        g = aa .- bb
        tangent = ldexp(BigFloat(sample.tangent), sample.tangent_power)
        cosine = (1 - tangent^2) / (1 + tangent^2)
        sine = 2tangent / (1 + tangent^2)
        radius = BigFloat(hypot(relative[1], relative[2]))
        magnitude = BigFloat(speed)
        cp, sp = BigFloat(sample.cosφ), BigFloat(sample.sinφ)
        rotation =
            iszero(radius) ? (magnitude * cp, magnitude * sp, BigFloat(0)) :
            (
                (g[1] * g[3] * cp - g[2] * magnitude * sp) / radius,
                (g[2] * g[3] * cp + g[1] * magnitude * sp) / radius,
                -radius * cp,
            )
        after = g .* cosine .+ rotation .* sine
        total = BigFloat(ma) + BigFloat(mb)
        center = (BigFloat(ma) .* aa .+ BigFloat(mb) .* bb) ./ total
        (T.(center .+ BigFloat(mb) / total .* after), T.(center .- BigFloat(ma) / total .* after))
    end
    collide_coulomb!(
        (a, b),
        grid,
        dt;
        coulomb_log = T[0 logarithm; logarithm 0],
        epsilon0 = 1,
        weight_to_number = 1,
        u_floor = 1,
        rng = DirectedDraw(draw, T(phase)),
    )
    actual = (ntuple(c -> a.v[c][1], 3), ntuple(c -> b.v[c][1], 3))
    errors = Float64[]
    for s = 1:2, c = 1:3
        got, want = actual[s][c], expected[s][c]
        @test isfinite(got)
        @test isapprox(got, want; rtol = 4096eps(T), atol = nextfloat(zero(T)))
        # The physical rate is rounded before sampling. A cancellation residual
        # is conditioned by both original velocity and the complete physical kick.
        before = BigFloat(original[s][c])
        physical_value = BigFloat(physical[s][c])
        scale = max(abs(before), abs(physical_value - before))
        @test abs(BigFloat(got) - physical_value) <=
              4096eps(T) * scale + BigFloat(nextfloat(zero(T)))
        push!(
            errors,
            iszero(want) ? (iszero(got) ? 0.0 : Inf) :
            Float64(abs((BigFloat(got) - BigFloat(want)) / BigFloat(want))),
        )
    end
    push!(
        rows,
        Dict(
            "type" => string(T),
            "masses" => collect(masses),
            "direction" => direction,
            "power" => power,
            "sign" => sign,
            "exchange" => exchange,
            "phase" => phase,
            "vanishing" => vanishing,
            "max_relative_error" => maximum(errors),
            "actual" => collect(Iterators.flatten(actual)),
            "expected" => collect(Iterators.flatten(expected)),
            "physical" => collect(Iterators.flatten(physical)),
            "sample_tangent" => Float64(sample.tangent),
            "sample_power" => sample.tangent_power,
            "input_roundoff" => collect(errors_in_inputs),
        ),
    )
end
@testset "Sampled update and conditioned physical-rate attribution" begin
    for T in (Float32, Float64),
        masses in ((1, 1), (1, 3), (3, 1)),
        direction = 1:5,
        power in (T === Float32 ? (-80, -12, 0, 12, 80) : (-600, -20, 0, 20, 600)),
        sign in (-1, 1),
        exchange in (false, true),
        phase in (0.0, 0.125, 0.3125)

        profile(T, masses, direction, power, sign, exchange, phase, false)
    end
    for T in (Float32, Float64), direction = 1:5, sign in (-1, 1), phase in (0.0, 0.125, 0.3125)
        profile(T, (1, 1), direction, 0, sign, false, phase, true)
    end
end

end
