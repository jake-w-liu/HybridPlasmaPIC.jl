module LocalCoulombBoundaryTests

using HybridPlasmaPIC, Test, Random

struct ZeroAngleRNG <: AbstractRNG end
Random.rand(::ZeroAngleRNG, range::UnitRange{Int}) = last(range)
Random.rand(::ZeroAngleRNG, T::Union{Type{Float32},Type{Float64}}) = T(0.125)
Random.randn(::ZeroAngleRNG, T::Union{Type{Float32},Type{Float64}}) = zero(T)
@testset "Zero scattering angle preserves unequal velocity scales" begin
    for T in (Float32, Float64),
        mass in (1, 4),
        scale in (sqrt(floatmin(T)), T(4) * nextfloat(zero(T)))

        grid = FourierGrid((4,), (T(4),))
        a, b = ParticleSet{1,T}(1), ParticleSet{1,T}(1; m = mass)
        a.x[1][1] = b.x[1][1] = T(0.25)
        a.v[3][1], b.v[3][1] = scale, one(T)
        before = map(ps -> map(copy, ps.v), (a, b))
        collide_coulomb!(
            (a, b),
            grid,
            T(0.001);
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = ZeroAngleRNG(),
        )
        @info "zero angle identity" T mass scale actual = a.v[3][1]
        @test isequal((a.v, b.v), before)
    end
end

const H = HybridPlasmaPIC
struct FixedRNG <: AbstractRNG end
Random.rand(::FixedRNG, range::UnitRange{Int}) = last(range)
Random.rand(::FixedRNG, T::Union{Type{Float32},Type{Float64}}) = T(0.125)
Random.randn(::FixedRNG, T::Union{Type{Float32},Type{Float64}}) = one(T)

@testset "Noncubic physical bins and range-safe Coulomb prefactor" begin
    for T in (Float32, Float64), D = 1:3, mass in (floatmin(T), one(T), floatmax(T) / 8)
        lengths = ntuple(d -> T(2d + 1), D)
        grid = FourierGrid(ntuple(_ -> 4, D), lengths)
        a, b = ParticleSet{D,T}(1; m = mass, q = sqrt(mass)),
        ParticleSet{D,T}(1; m = mass, q = -sqrt(mass))
        for ps in (a, b), x in ps.x
            fill!(x, T(0.1))
        end
        a.v[3][1] = one(T)
        dt = T(0.001)
        collide_coulomb!(
            (a, b),
            grid,
            dt;
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = FixedRNG(),
        )
        mu = BigFloat(mass) / 2
        variance =
            (BigFloat(a.q) * BigFloat(b.q))^2 * BigFloat(dt) /
            (8big(π) * mu^2 * prod(BigFloat.(grid.dx)))
        delta = sqrt(variance)
        cosine = (1 - delta^2) / (1 + delta^2)
        sine = 2delta / (1 + delta^2)
        relative = (sine / sqrt(big(2)), sine / sqrt(big(2)), cosine)
        for c = 1:3
            expected_a = (c == 3 ? big(0.5) : big(0)) + relative[c] / 2
            expected_b = (c == 3 ? big(0.5) : big(0)) - relative[c] / 2
            @test a.v[c][1] ≈ T(expected_a) atol = 1024eps(T) rtol = 1024eps(T)
            @test b.v[c][1] ≈ T(expected_b) atol = 1024eps(T) rtol = 1024eps(T)
        end
    end
end

@testset "Exact cold mass-update arithmetic" begin
    for T in (Float32, Float64),
        (ma, mb) in (
            (one(T), T(4)),
            (nextfloat(zero(T)), floatmax(T) / 4),
            (floatmax(T) / 4, nextfloat(zero(T))),
        )

        fa, fb = H._coulomb_mass_fractions(ma, mb)
        for (va, vb, g, gp) in (
            (zero(T), zero(T), one(T), -one(T)),
            (nextfloat(zero(T)), one(T), one(T), prevfloat(one(T))),
            (T(0.25), T(0.5), -floatmax(T) / 2, floatmax(T) / 2),
        )
            increment = gp - g
            @test isfinite(increment)
            a, b = H._coulomb_mass_update(va, vb, ma, mb, fa, fb, increment)
            qa, qb = Rational{BigInt}(ma), Rational{BigInt}(mb)
            delta = Rational{BigInt}(gp) - Rational{BigInt}(g)
            @test Rational{BigInt}(increment) == delta
            ea = T(Rational{BigInt}(va) + qb * delta / (qa + qb))
            eb = T(Rational{BigInt}(vb) - qa * delta / (qa + qb))
            @test a ≈ ea rtol = 8eps(T) atol = nextfloat(zero(T))
            @test b ≈ eb rtol = 8eps(T) atol = nextfloat(zero(T))
        end
    end
end

@testset "Strided particle storage and inactive channels" begin
    for T in (Float32, Float64)
        grid = FourierGrid((4,), (T(4),))
        source = ParticleSet{1,T}(12)
        fill!(source.x[1], T(0.25))
        rng = MersenneTwister(32)
        foreach(v -> randn!(rng, v), source.v)
        select = 1:2:11
        ps = ParticleSet{1,T}(
            map(x -> view(x, select), source.x),
            map(v -> view(v, select), source.v),
            view(source.weight, select),
            view(source.id, select),
            view(source.tag, select),
            1,
            1,
        )
        before = map(copy, source.v)
        reference = ParticleSet{1,T}(
            map(copy, ps.x),
            map(copy, ps.v),
            copy(ps.weight),
            copy(ps.id),
            copy(ps.tag),
            1,
            1,
        )
        for particle in (ps, reference)
            collide_coulomb!(
                (particle,),
                grid,
                T(0.01);
                coulomb_log = 1,
                epsilon0 = 1,
                weight_to_number = 1,
                rng = MersenneTwister(22),
            )
        end
        @test ps.v == reference.v
        @test all(c -> source.v[c][2:2:12] == before[c][2:2:12], 1:3)
        # Zero-weight slots cannot change the physical density or pair graph.
        padded = ParticleSet{1,T}(8)
        fill!(padded.x[1], T(0.25))
        for c = 1:3
            padded.v[c][1:6] .= reference.v[c]
        end
        padded.weight[7:8] .= 0
        for particle in (padded, reference)
            collide_coulomb!(
                (particle,),
                grid,
                T(0.01);
                coulomb_log = 1,
                epsilon0 = 1,
                weight_to_number = 1,
                rng = MersenneTwister(55),
            )
        end
        @test all(c -> padded.v[c][1:6] == reference.v[c], 1:3)
        @test all(c -> all(iszero, padded.v[c][7:8]), 1:3)
        for channel in (:neutral, :zero_log)
            particle = deepcopy(reference)
            channel == :neutral && (particle.q = 0)
            before = map(copy, particle.v)
            rng = MersenneTwister(37)
            copyrng = copy(rng)
            collide_coulomb!(
                (particle,),
                grid,
                T(1);
                coulomb_log = channel == :zero_log ? 0 : 1,
                epsilon0 = 1,
                weight_to_number = 1,
                rng,
            )
            @test particle.v == before
            @test rand(rng, UInt64) == rand(copyrng, UInt64)
        end
        @test_throws UndefKeywordError collide_coulomb!(
            (reference,),
            grid,
            T(0.01);
            coulomb_log = 1,
        )
    end
end
end

module LocalCoulombAxesWitnessTests

using HybridPlasmaPIC, Test, Random
struct ShiftedSpecies{P} <: AbstractVector{P}
    data::Vector{P}
end
Base.size(v::ShiftedSpecies) = size(v.data)
Base.axes(v::ShiftedSpecies) = (0:length(v.data)-1,)
Base.getindex(v::ShiftedSpecies, i::Int) = v.data[i+1]
Base.iterate(v::ShiftedSpecies, state = 1) = iterate(v.data, state)
@testset "Reject non-one-based species collections before RNG" begin
    ps = ParticleSet{1,Float64}(3)
    fill!(ps.x[1], 0.25)
    species = ShiftedSpecies([ps])
    grid = FourierGrid((4,), (4.0,))
    rng = MersenneTwister(1)
    unchanged = copy(rng)
    @test_throws ArgumentError collide_coulomb!(
        species,
        grid,
        0.01;
        coulomb_log = 1,
        epsilon0 = 1,
        weight_to_number = 1,
        rng,
    )
    @test rand(rng, UInt64) == rand(unchanged, UInt64)
end

end

module LocalCoulombSmallAngleWitnessTests

using Test, Random, HybridPlasmaPIC
struct SmallAngleRNG <: AbstractRNG
    normal::Float64
end
Random.rand(::SmallAngleRNG, r::UnitRange{Int}) = last(r)
Random.rand(::SmallAngleRNG, T::Union{Type{Float32},Type{Float64}}) = T(0.125)
Random.randn(r::SmallAngleRNG, T::Union{Type{Float32},Type{Float64}}) = T(r.normal)
@testset "Oblique small-angle physical increment accuracy" begin
    for T in (Float32, Float64), normal in (0.0, 1e-15, 1e-8, 1e-4, 1.0), mass in (1, 4)
        grid = FourierGrid((4,), (T(4),))
        a, b = ParticleSet{1,T}(1), ParticleSet{1,T}(1; m = mass)
        a.x[1][1] = b.x[1][1] = T(0.25)
        for c = 1:3
            b.v[c][1] = T(c)
        end
        dt = T(0.001)
        collide_coulomb!(
            (a, b),
            grid,
            dt;
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = SmallAngleRNG(normal),
        )
        ma, mb = big(1), BigFloat(mass)
        mu = ma * mb / (ma + mb)
        gx, gy, gz = -big(1), -big(2), -big(3)
        g = sqrt(gx^2 + gy^2 + gz^2)
        perp = sqrt(gx^2 + gy^2)
        delta = BigFloat(T(normal)) * sqrt(BigFloat(dt) / (8big(π) * mu^2 * g^3))
        cm1 = -2delta^2 / (1 + delta^2)
        sine = 2delta / (1 + delta^2)
        az = sqrt(big(0.5))
        change = (
            gx * cm1 + gx / perp * gz * sine * az - gy / perp * g * sine * az,
            gy * cm1 + gy / perp * gz * sine * az + gx / perp * g * sine * az,
            gz * cm1 - perp * sine * az,
        )
        expected = mb / (ma + mb) .* change
        tolerance = 128eps(T) * maximum(abs, expected)
        for c = 1:3
            @test abs(BigFloat(a.v[c][1]) - expected[c]) <= tolerance
        end
    end
end

end

module LocalCoulombExtremeIncrementWitnessTests

using Test, Random, HybridPlasmaPIC
mutable struct ReversalRNG <: AbstractRNG
    draw::Int
end
Random.rand(::ReversalRNG, r::UnitRange{Int}) = last(r)
function Random.rand(r::ReversalRNG, T::Union{Type{Float32},Type{Float64}})
    r.draw += 1
    return r.draw == 1 ? zero(T) : T(0.5)
end
Random.randn(::ReversalRNG, T::Union{Type{Float32},Type{Float64}}) =
    error("strong-scattering sample expected")
@testset "Unrepresentable relative increment with finite physical updates" begin
    for T in (Float32, Float64), mass in (1, 4)
        M = floatmax(T)
        g = FourierGrid((1,), (one(T),))
        a, b = ParticleSet{1,T}(1; q = M / 4), ParticleSet{1,T}(1; q = M / 4, m = mass)
        a.x[1][1] = b.x[1][1] = T(0.25)
        a.v[1][1] = -T(3 / 8) * M
        b.v[1][1] = T(3 / 8) * M
        va, vb = BigFloat(a.v[1][1]), BigFloat(b.v[1][1])
        collide_coulomb!(
            (a, b),
            g,
            one(T);
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = ReversalRNG(0),
        )
        change = -2 * (va - vb)
        @test a.v[1][1] ≈ T(va + BigFloat(mass) / (1 + mass) * change) rtol = 8eps(T)
        @test b.v[1][1] ≈ T(vb - change / (1 + mass)) rtol = 8eps(T)
    end
end
@testset "Rejected macro update does not impose a fictitious range error" begin
    for T in (Float32, Float64)
        M = floatmax(T)
        g = FourierGrid((1,), (one(T),))
        a, b = ParticleSet{1,T}(1; q = M / 4), ParticleSet{1,T}(1; q = M / 4, m = M / 4)
        a.x[1][1] = b.x[1][1] = T(0.25)
        a.v[1][1] = -T(3 / 8) * M
        b.v[1][1] = T(3 / 8) * M
        a.weight[1] = T(100)
        va = copy(a.v[1])
        vb = copy(b.v[1])
        @test collide_coulomb!(
            (a, b),
            g,
            one(T);
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = ReversalRNG(0),
        ) == (a, b)
        @test a.v[1] == va
        @test all(isfinite, b.v[1])
        # The same out-of-range proposal must fail atomically when it is accepted.
        a.weight[1] = one(T)
        before = map(ps -> map(copy, ps.v), (a, b))
        @test_throws ArgumentError collide_coulomb!(
            (a, b),
            g,
            one(T);
            coulomb_log = T[0 1; 1 0],
            epsilon0 = 1,
            weight_to_number = 1,
            rng = ReversalRNG(0),
        )
        @test isequal((a.v, b.v), before)
    end
end

end
