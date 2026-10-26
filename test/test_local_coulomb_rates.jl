module LocalCoulombPhysicalRateTests
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

mutable struct RateRNG{T} <: AbstractRNG
    draw::T
    events::Vector{Symbol}
end
Random.rand(rng::RateRNG, range::UnitRange{Int}) = (push!(rng.events, :shuffle); last(range))
for T in (Float32, Float64)
    @eval Random.randn(rng::RateRNG, ::Type{$T}) = (push!(rng.events, :normal); $T(rng.draw))
    @eval Random.rand(rng::RateRNG, ::Type{$T}) = (push!(rng.events, :uniform); $T(0.25))
end

function stored_rational(value::BigFloat)
    significand, exponent, sign = Base.decompose(value)
    numerator = BigInt(sign) * significand
    return exponent >= 0 ? (numerator << exponent) // BigInt(1) :
           numerator // (BigInt(1) << -exponent)
end

@testset "Physical pi comparison preserves directed high-precision endpoints" begin
    original = (precision(BigFloat), rounding(BigFloat))
    for bits in (24, 53, 64, 128, 256, 512, 1024)
        low = stored_rational(BigFloat(pi, Base.MPFR.MPFRRoundDown; precision = bits))
        high = stored_rational(BigFloat(pi, Base.MPFR.MPFRRoundUp; precision = bits))
        @test !H._coulomb_above_pi(low)
        @test H._coulomb_above_pi(high)
    end
    @test (precision(BigFloat), rounding(BigFloat)) == original
end

@testset "Public physical variance-one boundary retains original factors" begin
    for T in (Float32, Float64), side = -2:2
        dt = 2T(pi)
        for _ = 1:abs(side)
            dt = side < 0 ? prevfloat(dt) : nextfloat(dt)
        end
        ps = ParticleSet{1,T}(2)
        fill!(ps.x[1], T(0.25))
        ps.v[1][2] = one(T)
        grid = FourierGrid((4,), (T(4),))
        variance = setprecision(BigFloat, 512) do
            BigFloat(dt) / (2big(pi))
        end
        rng = RateRNG(zero(T), Symbol[])
        before = map(copy, ps.v)
        collide_coulomb!(
            (ps,),
            grid,
            dt;
            coulomb_log = one(T),
            epsilon0 = one(T),
            weight_to_number = one(T),
            rng,
        )
        @test rng.events[2] === (variance > 1 ? :uniform : :normal)
        @test isequal(before, ps.v) == (variance <= 1)
    end
end

function independent_variance(channel, count, reuse, triplet, weight, speed)
    setprecision(BigFloat, 1024) do
        mass_a, mass_b = BigFloat(channel.ma), BigFloat(channel.mb)
        reduced = mass_a * mass_b / (mass_a + mass_b)
        number_density =
            BigFloat(channel.number_scale) * BigFloat(weight) * BigFloat(count) / BigFloat(reuse) /
            prod(BigFloat.(channel.widths))
        triplet && (number_density /= 2)
        return (BigFloat(channel.qa) * BigFloat(channel.qb))^2 *
               BigFloat(channel.logarithm) *
               BigFloat(channel.dt) *
               number_density /
               (8big(pi) * BigFloat(channel.epsilon)^2 * reduced^2 * BigFloat(speed)^3)
    end
end

@testset "Sharp local branch retains masses, integer counts and physical pi" begin
    rejected = Dict(:rounded_mass => 0, :rounded_count => 0, :native_pi => 0)
    for T in (Float32, Float64),
        masses in ((1, 1), (1, 3), (1, 10), (7, 13), (0.125, 3.25)),
        count in (1, 3, 2^24 + 1, 2^53 + 1, 2^60 + 3),
        reuse in (1, 7),
        triplet in (false, true)

        grid = FourierGrid((4,), (T(4),))
        a = ParticleSet{1,T}(1; m = T(masses[1]))
        b = ParticleSet{1,T}(1; m = T(masses[2]))
        unit = H._coulomb_channel_rate(a, b, grid, one(T), one(T), one(T), one(T))
        boundary = T(inv(independent_variance(unit, count, reuse, triplet, one(T), one(T))))
        @test isfinite(boundary) && boundary > 0
        for side = -2:2
            dt = boundary
            for _ = 1:abs(side)
                dt = side < 0 ? prevfloat(dt) : nextfloat(dt)
            end
            channel = H._coulomb_channel_rate(a, b, grid, dt, one(T), one(T), one(T))
            rate = H._CoulombLocalPairRate(channel, count, reuse, triplet)
            exact = independent_variance(channel, count, reuse, triplet, one(T), one(T))
            isotropic = exact > 1
            @test H._coulomb_large_angle_exact(rate, one(T), one(T), one(T)) == isotropic
            fraction, power = H._coulomb_local_scaled_coefficient(rate)
            rng = RateRNG(zero(T), Symbol[])
            sample = H._coulomb_sample_increment(
                -one(T),
                zero(T),
                zero(T),
                one(T),
                ldexp(fraction, power),
                one(T),
                one(T),
                T(0.001),
                2T(pi),
                rng,
                rate,
            )
            @test sample.isotropic == isotropic
            @test rng.events == (isotropic ? [:uniform, :uniform] : [:normal, :uniform])
            if !isotropic
                @test begin
                    extreme = H._coulomb_sample_increment(
                        -one(T),
                        zero(T),
                        zero(T),
                        one(T),
                        ldexp(fraction, power),
                        one(T),
                        one(T),
                        T(0.001),
                        2T(pi),
                        RateRNG(floatmax(T), Symbol[]),
                        rate,
                    )
                    first, second = H._coulomb_complete_update(
                        (zero(T), zero(T), zero(T)),
                        (one(T), zero(T), zero(T)),
                        a.m,
                        b.m,
                        one(T),
                        extreme,
                    )
                    all(isfinite, (first..., second...))
                end
            end
            setprecision(BigFloat, 1024) do
                reduced = BigFloat(a.m) * BigFloat(b.m) / (BigFloat(a.m) + BigFloat(b.m))
                rounded = a.m * b.m / (a.m + b.m)
                rejected[:rounded_mass] +=
                    (exact * (reduced / BigFloat(rounded))^2 > 1) != isotropic
                rejected[:rounded_count] +=
                    (
                        exact * BigFloat(T(count) / T(reuse)) /
                        (BigFloat(count) / BigFloat(reuse)) > 1
                    ) != isotropic
                rejected[:native_pi] += (exact * big(pi) / BigFloat(T(pi)) > 1) != isotropic
            end
        end
    end
    for patch in keys(rejected)
        @test rejected[patch] > 0
    end
end

@noinline function decision_allocation(rate::R, wmax::T, wbar::S, speed::T) where {R,T,S}
    return @allocated H._coulomb_large_angle_exact(rate, wmax, wbar, speed)
end

@testset "Certified binary64 local decision retains original factors" begin
    original = (precision(BigFloat), rounding(BigFloat))
    for D = 1:3, masses in ((1, 1), (1, 3), (7, 13)), count in (1, 2^24 + 1, 2^53 + 1)
        T = Float32
        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        a = ParticleSet{D,T}(1; m = T(masses[1]))
        b = ParticleSet{D,T}(1; m = T(masses[2]))
        unit = H._coulomb_channel_rate(a, b, grid, one(T), one(T), one(T), one(T))
        boundary = T(inv(independent_variance(unit, count, 7, true, one(T), one(T))))
        for dt in (prevfloat(boundary), boundary, nextfloat(boundary))
            channel = H._coulomb_channel_rate(a, b, grid, dt, one(T), one(T), one(T))
            rate = H._CoulombLocalPairRate(channel, count, 7, true)
            exact = independent_variance(channel, count, 7, true, one(T), one(T))
            wide = H._coulomb_local_variance64(rate, one(T), one(T), one(T))
            setprecision(BigFloat, 1024) do
                u = BigFloat(eps(Float64)) / 2
                @test abs(BigFloat(wide) / exact - 1) <= (1 - u)^(-(32 + D)) - 1
            end
            @test abs(wide - 1) > H._coulomb_variance_roundoff(rate, Float64)
            @test H._coulomb_large_angle_exact(rate, one(T), one(T), one(T)) == (exact > 1)
            decision_allocation(rate, one(T), one(T), one(T))
            @test decision_allocation(rate, one(T), one(T), one(T)) == 0
        end
    end
    # Consecutive rational convergents of mathematical pi put both sides of
    # this physical rate inside the binary64 certificate's unresolved interval.
    sides = Set{Bool}()
    for D = 1:3, (numerator, denominator) in ((165707065, 52746197), (245850922, 78256779))
        T = Float32
        grid = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        a = ParticleSet{D,T}(1)
        channel = H._coulomb_channel_rate(a, a, grid, one(T), one(T), one(T), one(T))
        rate = H._CoulombLocalPairRate(channel, 2numerator, denominator, false)
        exact = independent_variance(channel, rate.count, rate.reuse, false, one(T), one(T))
        wide = H._coulomb_local_variance64(rate, one(T), one(T), one(T))
        @test abs(wide - 1) <= H._coulomb_variance_roundoff(rate, Float64)
        @test H._coulomb_large_angle_exact(rate, one(T), one(T), one(T)) == (exact > 1)
        @test H._coulomb_large_angle_exact(rate, one(T), 1.0, one(T)) == (exact > 1)
        @test H._coulomb_large_angle_exact(rate, one(T), BigFloat(1), one(T)) == (exact > 1)
        push!(sides, exact > 1)
    end
    @test sides == Set((false, true))
    @test (precision(BigFloat), rounding(BigFloat)) == original
end

end
