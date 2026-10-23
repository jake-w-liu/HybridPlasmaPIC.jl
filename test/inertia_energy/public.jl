module InertiaEnergyPublicTests
include("support.jl")
include("legacy.jl")
old_public(ps, B, n, c, g, d) = @timed original_budget(ps, B, n, c, g; de2 = d)
new_public(ps, B, n, c, g, d) = @timed energy_budget(ps, B, n, c, g; de2 = d)
@testset "Public inertia recovers finite spectral-current reservoirs" begin
    for T in (Float32, Float64), sign in (-1, 1)
        e = T === Float32 ? 80 : 600
        ce = T === Float32 ? 100 : 800
        amplitude = ldexp(one(T), sign * e)
        de2 = ldexp(one(T), -sign * ce)
        g = FourierGrid((4,), (T(2pi),))
        B = (zeros(T, 4), zeros(T, 4), amplitude .* T[1, 0, -1, 0])
        n = ones(T, 4)
        J = ntuple(_ -> zeros(T, 4), 3)
        H.curl!(J, B, g)
        expected = exact(J, de2, g.dx)
        before = map(copy, (B..., n))
        for ps in (ParticleSet{1,T}(0), ParticleSet{1,T}[])
            result = energy_budget(ps, B, n, IsothermalElectrons(zero(T)), g; de2 = de2)
            @test rounded(result.electron_inertia, expected)
            @test isfinite(result.electron_inertia) && result.electron_inertia > 0
            @test all(isequal.((B..., n), before))
        end
    end
    for T in (Float32, Float64)
        side = ldexp(one(T), T === Float32 ? 50 : 400)
        g = FourierGrid((4, 4, 4), (side, side, side))
        B = ntuple(_ -> zeros(T, 4, 4, 4), 3)
        result = energy_budget(
            ParticleSet{3,T}(0),
            B,
            ones(T, 4, 4, 4),
            IsothermalElectrons(zero(T)),
            g;
            de2 = T(0.25),
        )
        @test all(iszero, result)
    end
end
@testset "Public budget uses accepted components and retains allocation and validation" begin
    rng = MersenneTwister(678302)
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(d -> T(2pi), D))
        B = ntuple(_ -> randn(rng, T, g.n), 3)
        n = ones(T, g.n)
        ps = ParticleSet{D,T}(32)
        foreach(a -> fill!(a, zero(T)), ps.x)
        foreach(a -> copyto!(a, randn(rng, T, length(a))), ps.v)
        c = IsothermalElectrons(T(0.5))
        before = deepcopy((ps.x, ps.v, ps.weight, B, n))
        for population in (ps, [ps], view([ps], :)), de2 in (zero(T), T(0.25), T(2))
            old_public(population, B, n, c, g, de2)
            new_public(population, B, n, c, g, de2)
            old = old_public(population, B, n, c, g, de2)
            new = new_public(population, B, n, c, g, de2)
            if T === Float32
                @test new.value.kinetic === old.value.kinetic
                @test new.value.magnetic === old.value.magnetic
                @test new.value.electron_internal === old.value.electron_internal
                J = ntuple(_ -> zeros(T, g.n), 3)
                H.curl!(J, B, g)
                @test rounded(new.value.electron_inertia, exact(J, de2, g.dx))
                v = new.value
                @test v.total === v.kinetic + v.magnetic + v.electron_internal + v.electron_inertia
            else
                @test isequal(new.value, old.value)
            end
            @test new.bytes == old.bytes
            @test isequal((ps.x, ps.v, ps.weight, B, n), before)
        end
        for bad in (T(-1), T(Inf), T(NaN))
            @test_throws ArgumentError energy_budget(ps, B, n, c, g; de2 = bad)
            @test isequal((ps.x, ps.v, ps.weight, B, n), before)
        end
        @test_throws ArgumentError energy_budget(
            ps,
            B,
            n,
            EvolvingElectrons(T(2)),
            g;
            de2 = T(0.25),
            pe = zeros(T, g.n),
        )
        @test isequal((ps.x, ps.v, ps.weight, B, n), before)
    end
end
end
