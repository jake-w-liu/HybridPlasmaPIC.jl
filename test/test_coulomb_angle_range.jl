module TestCoulombAngleRange
using HybridPlasmaPIC, Random, Test

mutable struct AngleRNG{T} <: AbstractRNG
    draw::T
    phase::T
    uniforms::Int
    events::Vector{Symbol}
end
AngleRNG(draw::T, phase::T) where {T} = AngleRNG(draw, phase, 0, Symbol[])
function Random.shuffle!(rng::AngleRNG, values::AbstractVector)
    push!(rng.events, :shuffle)
    return values
end
function normal(rng::AngleRNG, ::Type{T}) where {T}
    push!(rng.events, :normal)
    return T(rng.draw)
end
function uniform(rng::AngleRNG, ::Type{T}) where {T}
    push!(rng.events, :uniform)
    rng.uniforms += 1
    return rng.uniforms == 1 ? T(rng.phase) : zero(T)
end
for T in (Float16, Float32, Float64, BigFloat)
    @eval Random.randn(rng::AngleRNG, ::Type{$T}) = normal(rng, $T)
    @eval Random.rand(rng::AngleRNG, ::Type{$T}) = uniform(rng, $T)
end

function cases(::Type{T}) where {T}
    exponents = if T === Float16
        (8, 8, -20, -20, -10, 12, 10)
    elseif T === Float32
        (40, 30, -110, -100, -60, 100, 70)
    else
        (350, 100, -600, -800, -500, 900, 620)
    end
    speed_exp, rate_exp, draw_exp, tiny_rate, tiny_speed, large_rate, large_speed = exponents
    speed = T === Float16 ? T(100) : T === Float32 ? T(1e12) : T(1e100)
    rate = T === Float16 ? T(0.001) : T === Float32 ? T(1e-12) : T(1e-100)
    return (
        (name = :variance, speed, rate, dt = one(T), draw = one(T), exchange = false),
        (
            name = :tangent,
            speed = ldexp(one(T), speed_exp),
            rate = ldexp(one(T), rate_exp),
            dt = one(T),
            draw = ldexp(one(T), draw_exp),
            exchange = false,
        ),
        (
            name = :rate_underflow,
            speed = ldexp(one(T), tiny_speed),
            rate = ldexp(one(T), tiny_rate),
            dt = ldexp(one(T), tiny_rate),
            draw = one(T),
            exchange = false,
        ),
        (
            name = :rate_overflow,
            speed = ldexp(one(T), large_speed),
            rate = ldexp(one(T), large_rate),
            dt = ldexp(one(T), large_rate),
            draw = one(T),
            exchange = false,
        ),
        (
            name = :small_draw,
            speed = one(T),
            rate = T(0.25),
            dt = one(T),
            draw = eps(one(T)),
            exchange = false,
        ),
        (
            name = :near_backscatter,
            speed = one(T),
            rate = T(0.25),
            dt = one(T),
            draw = inv(eps(one(T))),
            exchange = true,
        ),
    )
end

# Independent complete physical formula, preserving the original finite rate
# factors. Geometry and sin/cos are the native sampled quantities; the rational
# Cayley rotation and final velocity addition are evaluated at high precision.
function reference(ps, row, phase, floor)
    T = eltype(ps.weight)
    va, vb = ntuple(c -> ps.v[c][1], 3), ntuple(c -> ps.v[c][2], 3)
    g = va .- vb
    magnitude = hypot(g...)
    radius = hypot(g[1], g[2])
    scale = maximum(ps.weight)
    normalized = ps.weight ./ scale
    average = sum(Rational{BigInt}, normalized) / 2
    maxweight = maximum(normalized)
    sine_phase, cosine_phase = sincos((2T(pi)) * phase)
    return setprecision(BigFloat, 2048) do
        a, b = BigFloat.(va), BigFloat.(vb)
        difference = a .- b
        speed = BigFloat(max(magnitude, floor))
        variance =
            BigFloat(row.rate) * BigFloat(row.dt) * BigFloat(maxweight) /
            (BigFloat(average) * speed^3)
        @assert 0 < variance <= 1
        tangent = sqrt(variance) * BigFloat(row.draw)
        squared = tangent^2
        longitudinal = -2squared / (1 + squared)
        transverse = 2tangent / (1 + squared)
        cp, sp, m, r =
            BigFloat(cosine_phase), BigFloat(sine_phase), BigFloat(magnitude), BigFloat(radius)
        orthogonal = if iszero(radius)
            (m * cp, m * sp, zero(m))
        else
            (
                (difference[1] * difference[3] * cp - difference[2] * m * sp) / r,
                (difference[2] * difference[3] * cp + difference[1] * m * sp) / r,
                -r * cp,
            )
        end
        change = difference .* longitudinal .+ orthogonal .* transverse
        a .+ change ./ 2, b .- change ./ 2
    end
end

@testset "Complete Coulomb rate, tangent and velocity expression" begin
    for T in (Float16, Float32, Float64, BigFloat),
        row in cases(T),
        D in (1, 3),
        axis in (1, 2, 3),
        oblique in (false, true),
        direction in (-1, 1),
        unequal in (false, true)

        ps = ParticleSet{D,T}(2)
        foreach(x -> fill!(x, T(0.25)), ps.x)
        for c = 1:3
            component = mod1(c + axis - 1, 3)
            value = direction * row.speed * (oblique ? T(c) / T(3) : T(c == 1))
            ps.v[component][1] = row.exchange ? value : zero(T)
            ps.v[component][2] = row.exchange ? zero(T) : value
        end
        unequal && (ps.weight[2] = T(3))
        phase = isodd(axis) ? zero(T) : T(0.125)
        floor = row.speed / T(16)
        original = deepcopy(ps)
        expected = reference(ps, row, phase, floor)
        relative_tolerance = BigFloat(32eps(one(T)))
        rng = AngleRNG(row.draw, phase)
        @test collide_coulomb!(
            ps,
            row.rate,
            row.dt;
            rng,
            u_floor = floor,
            index_work = Vector{Int}(undef, 2),
        ) === ps
        @test rng.events == (
            unequal ? [:shuffle, :normal, :uniform, :uniform, :uniform] :
            [:shuffle, :normal, :uniform]
        )
        @test isequal(ps.x, original.x) && ps.weight == original.weight && ps.id == original.id
        for p = 1:2, c = 1:3
            actual, target = ps.v[c][p], expected[p][c]
            rounded = T(target)
            tolerance = max(relative_tolerance * abs(target), BigFloat(eps(rounded)) / 2)
            @test isfinite(actual)
            @test abs(BigFloat(actual) - target) <= tolerance
        end
    end
end

@testset "Coulomb subnormal weight quotient retains the original factors" begin
    for T in (Float16, Float32, Float64), units in (1, 2, 3, 7)
        ps = ParticleSet{1,T}(6)
        ps.v[1][2] = one(T)
        weight = T(units) * nextfloat(zero(T))
        ps.weight[1:2] .= weight
        gc = floatmax(T) / T(64)
        # Include all six native normalized weights in the exact population mean.
        # The first pair's quotient is subnormal before its large rate restores
        # a normal variance; compute the independent reference from the factors.
        average = (4 + 2Rational{BigInt}(weight)) / 6
        expected = setprecision(BigFloat, 4096) do
            variance = BigFloat(gc) * BigFloat(weight) / BigFloat(average)
            delta = sqrt(variance)
            x = delta^2 / (1 + delta^2)
            z = -delta / (1 + delta^2)
            ((x, zero(x), z), (1 - x, zero(x), -z))
        end
        collide_coulomb!(
            ps,
            gc,
            one(T);
            rng = AngleRNG(one(T), zero(T)),
            index_work = Vector{Int}(undef, 6),
        )
        for p = 1:2, c = 1:3
            target = expected[p][c]
            bound = max(BigFloat(32eps(T)) * abs(target), BigFloat(eps(T(target))) / 2)
            @test abs(BigFloat(ps.v[c][p]) - target) <= bound
        end
    end
end

@testset "Coulomb regime boundary retains the original rate factors" begin
    for T in (Float16, Float32, Float64, BigFloat)
        powers = T === Float16 ? (-10, 6) : T === Float32 ? (-80, 50) : (-600, 400)
        for power in powers, side in (:below, :equal, :above)
            # Binary powers make the exact variance one. Its immediate rate
            # neighbors lie strictly below/above one even when the product or
            # speed cube is outside the hardware exponent range.
            rate_power = fld(3power, 2)
            gc = ldexp(one(T), rate_power)
            gc = side === :below ? prevfloat(gc) : side === :above ? nextfloat(gc) : gc
            dt = ldexp(one(T), 3power - rate_power)
            ps = ParticleSet{1,T}(2)
            fill!(ps.x[1], T(0.25))
            ps.v[1][2] = ldexp(one(T), power)
            rng = AngleRNG(zero(T), T(0.25))
            collide_coulomb!(
                ps,
                gc,
                dt;
                rng,
                u_floor = ps.v[1][2] / 2,
                index_work = Vector{Int}(undef, 2),
            )
            expected =
                side === :above ? [:shuffle, :uniform, :uniform] : [:shuffle, :normal, :uniform]
            @test rng.events == expected
            @test all(v -> all(isfinite, v), ps.v)
        end
    end
end
end
