module LocalCoulombComponentRangeTests
using HybridPlasmaPIC, Test, Random
struct ComponentDraw{T} <: AbstractRNG
    normal::T
end
Random.rand(::ComponentDraw, r::UnitRange{Int}) = last(r)
Random.rand(::ComponentDraw, T::Union{Type{Float32},Type{Float64}}) = zero(T)
Random.randn(r::ComponentDraw, T::Union{Type{Float32},Type{Float64}}) = T(r.normal)
rows = Any[]
@testset "Finite component kicks after unit-direction underflow" begin
    for T in (Float32, Float64),
        masses in ((1, 1), (1, 3), (3, 1)),
        small_axis in (1, 2),
        direction in (-1, 1),
        normal_sign in (-1, 1)

        large_power, small_power, draw_power = T === Float32 ? (100, 0, -30) : (900, 200, -250)
        large = ldexp(one(T), large_power)
        small = direction * ldexp(one(T), small_power)
        draw = normal_sign * ldexp(one(T), draw_power)
        g = small_axis == 1 ? (small, large, zero(T)) : (large, small, zero(T))
        ma, mb = T.(masses)
        q = ldexp(one(T), div(large_power, 2))
        dt = large
        mu = ma * mb / (ma + mb)
        logarithm = T(4 * T(pi) * mu^2)
        grid = FourierGrid((4,), (T(4),))
        a, b = ParticleSet{1,T}(1; q, m = ma), ParticleSet{1,T}(1; q, m = mb)
        for ps in (a, b)
            ps.x[1][1] = T(0.25)
            ps.weight[1] = one(T)
        end
        for c = 1:3
            a.v[c][1] = zero(T)
            b.v[c][1] = -g[c]
        end
        expected = setprecision(BigFloat, 4096) do
            speed = sqrt(sum(x -> BigFloat(x)^2, g))
            reduced = BigFloat(ma) * BigFloat(mb) / (BigFloat(ma) + BigFloat(mb))
            variance =
                BigFloat(q)^4 * BigFloat(logarithm) * BigFloat(dt) /
                (8 * BigFloat(T(pi)) * reduced^2 * speed^3)
            angle_square = variance * BigFloat(draw)^2
            T(
                -2BigFloat(mb) / (BigFloat(ma) + BigFloat(mb)) * BigFloat(small) * angle_square / (1 + angle_square),
            )
        end
        collide_coulomb!(
            (a, b),
            grid,
            dt;
            coulomb_log = T[0 logarithm; logarithm 0],
            epsilon0 = 1,
            weight_to_number = 1,
            u_floor = 1,
            rng = ComponentDraw(draw),
        )
        got = a.v[small_axis][1]
        push!(
            rows,
            Dict(
                "type" => string(T),
                "masses" => collect(masses),
                "small_axis" => small_axis,
                "direction" => direction,
                "normal_sign" => normal_sign,
                "actual" => Float64(got),
                "expected" => Float64(expected),
            ),
        )
        @test isfinite(got)
        @test !iszero(expected)
        @test isapprox(got, expected; rtol = 4096eps(T), atol = zero(T))
    end
end

end
