module LocalCoulombCompleteAngleTests
using HybridPlasmaPIC, Test, Random
struct PrescribedNormal{T} <: AbstractRNG
    draw::T
end
Random.rand(::PrescribedNormal, r::UnitRange{Int}) = last(r)
Random.rand(::PrescribedNormal, T::Union{Type{Float32},Type{Float64}}) = zero(T)
Random.randn(r::PrescribedNormal, T::Union{Type{Float32},Type{Float64}}) = T(r.draw)
rows = Any[]
function profile(T, ma, mb, power, direction, exchange, vanishing_variance)
    gexp = T === Float32 ? 100 : 900
    speed = ldexp(one(T), gexp)
    charge = vanishing_variance ? one(T) : ldexp(one(T), div(gexp, 2))
    dt = vanishing_variance ? one(T) : speed
    draw = direction * ldexp(one(T), power)
    mass_a, mass_b = T(ma), T(mb)
    fraction = mass_b / (mass_a + mass_b)
    mu = mass_a * mass_b / (mass_a + mass_b)
    logarithm = T(4 * T(pi) * mu^2)
    grid = FourierGrid((4,), (T(4),))
    a = ParticleSet{1,T}(1; q = charge, m = mass_a)
    b = ParticleSet{1,T}(1; q = charge, m = mass_b)
    for ps in (a, b)
        ps.x[1][1] = T(0.25)
        ps.weight[1] = one(T)
        for v in ps.v
            fill!(v, zero(T))
        end
    end
    a.v[3][1] = exchange ? 2fraction * speed : zero(T)
    b.v[3][1] = a.v[3][1] - speed
    expected = setprecision(BigFloat, 512) do
        mass_fraction = BigFloat(mass_b) / (BigFloat(mass_a) + BigFloat(mass_b))
        reduced = BigFloat(mass_a) * BigFloat(mass_b) / (BigFloat(mass_a) + BigFloat(mass_b))
        variance =
            BigFloat(charge)^4 * BigFloat(logarithm) * BigFloat(dt) /
            (8 * BigFloat(T(pi)) * reduced^2 * BigFloat(speed)^3)
        q = variance * BigFloat(draw)^2
        tangent = copysign(sqrt(q), BigFloat(draw))
        transverse = 2mass_fraction * BigFloat(speed) * tangent / (1 + q)
        # These algebraic forms avoid subtracting nearly equal high-precision values.
        longitudinal = 2mass_fraction * BigFloat(speed) * (exchange ? one(q) : -q) / (1 + q)
        (T(transverse), T(longitudinal))
    end
    collide_coulomb!(
        (a, b),
        grid,
        dt;
        coulomb_log = T[0 logarithm; logarithm 0],
        epsilon0 = 1,
        weight_to_number = 1,
        u_floor = 1,
        rng = PrescribedNormal(draw),
    )
    push!(
        rows,
        Dict(
            "type" => string(T),
            "ma" => ma,
            "mb" => mb,
            "draw_power" => power,
            "direction" => direction,
            "exchange" => exchange,
            "vanishing_variance" => vanishing_variance,
            "actual_x" => Float64(a.v[1][1]),
            "expected_x" => Float64(expected[1]),
            "actual_z" => Float64(a.v[3][1]),
            "expected_z" => Float64(expected[2]),
        ),
    )
    @test all(ps -> all(v -> all(isfinite, v), ps.v), (a, b))
    @test a.v[1][1] ≈ expected[1] rtol = 4096eps(T) atol = zero(T)
    @test a.v[3][1] ≈ expected[2] rtol = 4096eps(T) atol = zero(T)
end
@testset "Complete angle and physical-update matrix" begin
    for T in (Float32, Float64),
        (ma, mb) in ((1, 1), (1, 3), (3, 1)),
        power in
        (T === Float32 ? (-80, -24, -12, 0, 12, 24, 80) : (-600, -40, -20, 0, 20, 40, 600)),
        direction in (-1, 1),
        exchange in (false, true)

        profile(T, ma, mb, power, direction, exchange, false)
    end
    for T in (Float32, Float64), direction in (-1, 1)
        profile(T, 1, 1, 0, direction, false, true)
    end
end

end
