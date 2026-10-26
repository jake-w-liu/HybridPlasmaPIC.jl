module LocalCoulombWidenedUpdateTests
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
        map(x -> H._from_exact(Float32, x), center .+ rational(mb) / total .* after),
        map(x -> H._from_exact(Float32, x), center .- rational(ma) / total .* after),
    )
end
rng = MersenneTwister(7613)
randomvalue() = ldexp(rand(rng, Float32) + 0.5f0, rand(rng, -140:120)) * rand(rng, (-1.0f0, 1.0f0))
rows = Any[]
counts = Dict("widened" => 0, "exact" => 0)
@testset "Independent Float32 widened-update boundary cells" begin
    for i = 1:2400
        va = ntuple(_ -> randomvalue(), 3)
        vb = ntuple(_ -> randomvalue(), 3)
        ma, mb = abs(randomvalue()), abs(randomvalue())
        relative = va .- vb
        speed = hypot(relative...)
        if !isfinite(speed) || iszero(speed)
            continue
        end
        tangent = randomvalue()
        phi = Float32(2pi) * rand(rng, Float32)
        sp, cp = sincos(phi)
        isotropic = iszero(i % 3)
        cosine = 2rand(rng, Float32) - 1.0f0
        sine = sqrt(max(0.0f0, 1.0f0 - cosine^2))
        power = iszero(i % 5) ? rand(rng, -1000:0) : 0
        sample = (; tangent, tangent_power = power, isotropic, cosine, sine, cosφ = cp, sinφ = sp)
        candidate = H._coulomb_widened_update(va, vb, ma, mb, speed, sample)
        counts[candidate === nothing ? "exact" : "widened"] += 1
        actual = H._coulomb_complete_update(va, vb, ma, mb, speed, sample)
        expected = oracle(va, vb, ma, mb, speed, sample)
        for s = 1:2, c = 1:3
            @test isequal(actual[s][c], expected[s][c]) || isapprox(
                actual[s][c],
                expected[s][c];
                rtol = 8eps(Float32),
                atol = nextfloat(0.0f0),
            )
        end
        if !all(
            s -> all(
                c ->
                    isequal(actual[s][c], expected[s][c]) || isapprox(
                        actual[s][c],
                        expected[s][c];
                        rtol = 8eps(Float32),
                        atol = nextfloat(0.0f0),
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
@testset "Float32 completion on both sides of the overflow midpoint" begin
    for k in (-20, -18, -16, -14, -12),
        sign in (-1.0f0, 1.0f0),
        coordinate = 1:3,
        side in (-1, 0, 1)

        delta = ldexp(1.0f0, k)
        ma, mb = 1.0f0 + delta, 1.0f0 - delta
        base = ldexp(1.0f0, 104) * (1.0f0 + delta)
        # Four Float32 ulps survive subtraction from the larger common center.
        magnitude = side < 0 ? prevfloat(base, 4) : side > 0 ? nextfloat(base, 4) : base
        center = ldexp(1.0f0, 107)
        va = Base.setindex((center, center, center), sign * floatmax(Float32), coordinate)
        relative_axis = coordinate == 3 ? 1 : 3
        vb = Base.setindex(va, va[relative_axis] - magnitude, relative_axis)
        speed = hypot((va .- vb)...)
        cp, sp = coordinate == 1 ? (sign, 0.0f0) : coordinate == 2 ? (0.0f0, sign) : (-sign, 0.0f0)
        sample = (;
            tangent = 1.0f0,
            tangent_power = 0,
            isotropic = false,
            cosine = 0.0f0,
            sine = 1.0f0,
            cosφ = cp,
            sinφ = sp,
        )
        @test speed == magnitude
        expected = oracle(va, vb, ma, mb, speed, sample)
        @test isfinite(expected[1][coordinate]) == (side <= 0)
        @test H._coulomb_widened_update(va, vb, ma, mb, speed, sample) === nothing
        actual = H._coulomb_complete_update(va, vb, ma, mb, speed, sample)
        for s = 1:2, c = 1:3
            @test isequal(actual[s][c], expected[s][c])
        end
    end
end

struct OverflowMidpointRNG <: AbstractRNG end
Random.rand(::OverflowMidpointRNG, ::Type{Float32}) = 0.0f0
Random.randn(::OverflowMidpointRNG, ::Type{Float32}) = sqrt(2.0f0)
@testset "Public collision finite completion and genuine-overflow rollback" begin
    for k in (-20, -18, -16), side in (-1, 0, 1)
        delta = ldexp(1.0f0, k)
        ma, mb = 1.0f0 + delta, 1.0f0 - delta
        base = ldexp(1.0f0, 104) * (1.0f0 + delta)
        magnitude = side < 0 ? prevfloat(base, 4) : side > 0 ? nextfloat(base, 4) : base
        center = ldexp(1.0f0, 107)
        va = (floatmax(Float32), center, center)
        vb = (floatmax(Float32), center, center - magnitude)
        grid = FourierGrid((1,), (1.0f0,))
        species = (
            ParticleSet{1,Float32}(1; m = ma, q = ldexp(1.0f0, 50)),
            ParticleSet{1,Float32}(1; m = mb, q = ldexp(1.0f0, 50)),
        )
        for (ps, v) in zip(species, (va, vb))
            fill!(ps.x[1], 0.25f0)
            for c = 1:3
                ps.v[c][1] = v[c]
            end
        end
        before = map(ps -> map(copy, ps.v), species)
        work = CoulombCollisionWorkspace(species, grid)
        # These binary scales give variance approximately1/2, keeping the
        # sqrt(2) Gaussian draw near tangent1 without overflowing physical inputs.
        call =
            () -> collide_coulomb!(
                species,
                grid,
                1.0f0;
                coulomb_log = Float32(pi),
                epsilon0 = 1.0f0,
                weight_to_number = ldexp(1.0f0, 112),
                workspace = work,
                rng = OverflowMidpointRNG(),
            )
        if side <= 0
            @test call() === species
            @test species[1].v[1][1] == floatmax(Float32)
            @test all(ps -> all(v -> all(isfinite, v), ps.v), species)
            @test any(s -> !isequal(species[s].v, before[s]), 1:2)
        else
            @test_throws ArgumentError call()
            @test all(s -> isequal(species[s].v, before[s]), 1:2)
        end
    end
end

end
