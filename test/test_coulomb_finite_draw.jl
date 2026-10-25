module TestCoulombFiniteDraw
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC
mutable struct FiniteDrawRNG{T} <: AbstractRNG
    draw::T
    uniforms::Int
    events::Vector{Symbol}
end
FiniteDrawRNG(draw::T) where {T} = FiniteDrawRNG(draw, 0, Symbol[])
Random.shuffle!(rng::FiniteDrawRNG, values::AbstractVector) = (push!(rng.events, :shuffle); values)
function uniform(rng::FiniteDrawRNG, ::Type{T}) where {T}
    push!(rng.events, :uniform)
    rng.uniforms += 1
    return rng.uniforms == 1 ? T(0.25) : zero(T)
end
for T in (Float32, Float64)
    @eval Random.randn(rng::FiniteDrawRNG, ::Type{$T}) = (push!(rng.events, :normal); $T(rng.draw))
    @eval Random.rand(rng::FiniteDrawRNG, ::Type{$T}) = uniform(rng, $T)
end
# Stored inputs from independent unequal-weight variance-one witnesses. Keeping
# bits makes the original boundary reproducible without a stochastic search.
const CASES = (
    (; speed = 0x3f2cb2ad, dt = 0x3f97c112, coefficient = 0x3e5d879f, weight = 0x3f2bc048),
    (; speed = 0x3f29c711, dt = 0x3faf8330, coefficient = 0x3e4af509, weight = 0x3f5d0931),
    (; speed = 0x3f895432, dt = 0x3f25b1fe, coefficient = 0x3fed4722, weight = 0x3f716df8),
    (;
        speed = 0x3fe16266893b05d4,
        dt = 0x3fe70dd70d1932a6,
        coefficient = 0x3fc7f4b7aaa1bbb8,
        weight = 0x3fe5d237457a8356,
    ),
    (;
        speed = 0x3ff19a23a456d71d,
        dt = 0x3ff2aca7c03652ff,
        coefficient = 0x3fed1cf184bffcad,
        weight = 0x3fe30a42cd5fdce2,
    ),
    (;
        speed = 0x3ff5e813d73459a2,
        dt = 0x3fe91d6c8542dcb8,
        coefficient = 0x40046744d2241a48,
        weight = 0x3fe1e9b0544aab64,
    ),
)
scalar(bits::UInt32) = reinterpret(Float32, bits)
scalar(bits::UInt64) = reinterpret(Float64, bits)
function oracle(speed, coefficient, dt, average, draw)
    T = typeof(speed)
    setprecision(BigFloat, 4096) do
        variance = BigFloat(coefficient) * BigFloat(dt) / (BigFloat(average) * BigFloat(speed)^3)
        tangent = sqrt(variance) * BigFloat(draw)
        denominator = 1 + tangent^2
        sine, cosine = sincos((2T(pi)) * T(0.25))
        kick = BigFloat(speed) * tangent / denominator
        first = (
            BigFloat(speed) * tangent^2 / denominator,
            -kick * BigFloat(sine),
            -kick * BigFloat(cosine),
        )
        second = (BigFloat(speed) / denominator, -first[2], -first[3])
        return map(x -> T.(x), (first, second))
    end
end
@testset "Finite Gaussian draw survives rounded variance overshoot" begin
    for bits in CASES
        speed, dt, coefficient, weight =
            map(scalar, (bits.speed, bits.dt, bits.coefficient, bits.weight))
        T = typeof(speed)
        average = Float64((1 + Rational{BigInt}(weight)) / 2)
        exact =
            Rational{BigInt}(coefficient) * Rational{BigInt}(dt) /
            (Rational{BigInt}(average) * Rational{BigInt}(speed)^3)
        approximate = ((coefficient * dt) * T(one(T) / average)) / ((speed * speed) * speed)
        @test exact <= 1
        @test sqrt(approximate) > one(T)
        @test isinf(sqrt(approximate) * floatmax(T))
        @test H._coulomb_normalized_mean(T[1, weight], one(T)) == average
        for coefficient in (coefficient, prevfloat(coefficient)),
            draw in (
                zero(T),
                floatmax(T),
                -floatmax(T),
                prevfloat(floatmax(T)),
                -prevfloat(floatmax(T)),
                floatmax(T) / 4,
                -floatmax(T) / 4,
            )

            ps = ParticleSet{1,T}(2)
            fill!(ps.x[1], T(0.25))
            ps.v[1][2] = speed
            ps.weight[2] = weight
            before =
                (x = deepcopy(ps.x), weight = copy(ps.weight), id = copy(ps.id), tag = copy(ps.tag))
            rng = FiniteDrawRNG(draw)
            scratch = map(copy, ps.v)
            collide_coulomb!(ps, coefficient, dt; work = scratch, index_work = zeros(Int, 2), rng)
            expected = oracle(speed, coefficient, dt, average, draw)
            for p = 1:2, c = 1:3
                @test isfinite(ps.v[c][p])
                @test isapprox(
                    ps.v[c][p],
                    expected[p][c];
                    rtol = 64eps(T),
                    atol = 4nextfloat(zero(T)),
                )
            end
            @test isequal(ps.v, scratch)
            @test ps.x == before.x &&
                  ps.weight == before.weight &&
                  ps.id == before.id &&
                  ps.tag == before.tag
            @test rng.events == [:shuffle, :normal, :uniform, :uniform, :uniform]
        end
    end
end
@testset "Nonfinite Gaussian draws retain transactional rejection" begin
    for bits in CASES
        speed, dt, coefficient, weight =
            map(scalar, (bits.speed, bits.dt, bits.coefficient, bits.weight))
        T = typeof(speed)
        for draw in (T(NaN), T(Inf), T(-Inf))
            ps = ParticleSet{1,T}(2)
            fill!(ps.x[1], T(0.25))
            ps.v[1][2] = speed
            ps.weight[2] = weight
            before = map(copy, ps.v)
            rng = FiniteDrawRNG(draw)
            @test_throws ArgumentError collide_coulomb!(
                ps,
                coefficient,
                dt;
                index_work = zeros(Int, 2),
                rng,
            )
            @test isequal(before, ps.v)
            @test rng.events == [:shuffle, :normal]
        end
    end
end
end
