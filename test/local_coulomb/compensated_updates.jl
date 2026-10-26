module LocalCoulombCompensatedUpdateTests
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
rational(x) = Rational{BigInt}(x)
function oracle(va, vb, ma, mb, speed, sample)
    ga, gb = map(rational, va), map(rational, vb)
    g = ga .- gb
    if sample.isotropic
        cosine = rational(sample.cosine)
        sine = rational(sample.sine)
    else
        tangent = rational(sample.tangent)
        k = sample.tangent_power
        tangent = k >= 0 ? tangent * (big(1) << k) : tangent / (big(1) << -k)
        cosine = (1 - tangent^2) / (1 + tangent^2)
        sine = 2tangent / (1 + tangent^2)
    end
    radius = rational(hypot(va[1] - vb[1], va[2] - vb[2]))
    mag = rational(speed)
    cp, sp = rational(sample.cosφ), rational(sample.sinφ)
    rotation =
        iszero(radius) ? (mag * cp, mag * sp, zero(mag)) :
        (
            (g[1] * g[3] * cp - g[2] * mag * sp) / radius,
            (g[2] * g[3] * cp + g[1] * mag * sp) / radius,
            -radius * cp,
        )
    after = g .* cosine .+ rotation .* sine
    total = rational(ma) + rational(mb)
    center = (rational(ma) .* ga .+ rational(mb) .* gb) ./ total
    (
        map(x -> H._from_exact(Float64, x), center .+ rational(mb) / total .* after),
        map(x -> H._from_exact(Float64, x), center .- rational(ma) / total .* after),
    )
end
rng = MersenneTwister(7613)
randomvalue() = ldexp(rand(rng, Float64) + 0.5f0, rand(rng, -101:101)) * rand(rng, (-1.0f0, 1.0f0))
rows = Any[]
counts = Dict("widened" => 0, "exact" => 0)
@testset "Independent compensated Float64 implementation boundary cells" begin
    for i = 1:10000
        va = ntuple(_ -> randomvalue(), 3)
        vb = ntuple(_ -> randomvalue(), 3)
        ma, mb = abs(randomvalue()), abs(randomvalue())
        relative = va .- vb
        speed = hypot(relative...)
        if !isfinite(speed) || iszero(speed)
            continue
        end
        tangent = randomvalue()
        phi = Float64(2pi) * rand(rng, Float64)
        sp, cp = sincos(phi)
        isotropic = iszero(i % 3)
        cosine = 2rand(rng, Float64) - 1.0f0
        sine = sqrt(max(0.0f0, 1.0f0 - cosine^2))
        power = iszero(i % 5) ? rand(rng, -1000:0) : 0
        sample = (; tangent, tangent_power = power, isotropic, cosine, sine, cosφ = cp, sinφ = sp)
        candidate = H._coulomb_compensated_update(va, vb, ma, mb, speed, sample)
        counts[candidate === nothing ? "exact" : "widened"] += 1
        actual =
            candidate === nothing ? H._coulomb_complete_update(va, vb, ma, mb, speed, sample) :
            candidate
        expected = oracle(va, vb, ma, mb, speed, sample)
        for s = 1:2, c = 1:3
            @test isequal(actual[s][c], expected[s][c]) || isapprox(
                actual[s][c],
                expected[s][c];
                rtol = 8eps(Float64),
                atol = nextfloat(0.0),
            )
        end
        if !all(
            s -> all(
                c ->
                    isequal(actual[s][c], expected[s][c]) || isapprox(
                        actual[s][c],
                        expected[s][c];
                        rtol = 8eps(Float64),
                        atol = nextfloat(0.0),
                    ),
                1:3,
            ),
            1:2,
        )
            push!(
                rows,
                Dict(
                    "index" => i,
                    "va" => collect(va),
                    "vb" => collect(vb),
                    "ma" => ma,
                    "mb" => mb,
                    "tangent" => tangent,
                    "power" => power,
                    "isotropic" => isotropic,
                    "cosine" => cosine,
                    "sine" => sine,
                    "phase" => phi,
                    "widened" => candidate !== nothing,
                    "actual" => collect(Iterators.flatten(actual)),
                    "expected" => collect(Iterators.flatten(expected)),
                ),
            )
        end
    end
    @test counts["widened"] > 100
    @test counts["exact"] > 100
end

end
