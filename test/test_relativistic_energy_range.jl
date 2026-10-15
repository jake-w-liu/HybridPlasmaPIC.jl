# Complete relativistic energy, ordinary-order compatibility and range certificate.
using HybridPlasmaPIC, Test, Random

@testset "Relativistic radical cancellation and midpoint termination" begin
    Q = Rational{BigInt}
    for T in (Float16, Float32, Float64)
        U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
        for lower in (one(T), nextfloat(one(T)), prevfloat(floatmax(T)))
            midpoint = (Q(lower) + Q(nextfloat(lower))) / 2
            result = HybridPlasmaPIC._finish_relativistic_energy(T, midpoint, Q[2, 8], Q[1, -1//2])
            @test result === (iseven(reinterpret(U, lower)) ? lower : nextfloat(lower))
        end
        @test HybridPlasmaPIC._finish_relativistic_energy(T, zero(Q), Q[2, 8], Q[1, -1//2]) ===
              zero(T)
    end
    for bits in (32, 128, 256), mode in (RoundNearest, RoundDown, RoundUp)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                @test iszero(
                    HybridPlasmaPIC._finish_relativistic_energy(
                        BigFloat,
                        zero(Q),
                        Q[2, 8],
                        Q[1, -1//2],
                    ),
                )
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
            end
        end
    end
end

function _relativistic_energy_ordinary_reference(ps::ParticleSet{D,T}, c::Real) where {D,T}
    cT = HybridPlasmaPIC._require_finite_positive_real("c", c, T)
    c2 = cT * cT
    s = zero(T)
    vx, vy, vz = ps.v
    @inbounds for p in eachindex(ps.weight)
        v2 = vx[p]^2 + vy[p]^2 + vz[p]^2
        v2 < c2 || throw(
            ArgumentError(
                "particle $p has |v| ≥ c (|v|² = $v2, c² = $c2): not a valid relativistic state",
            ),
        )
        # (γ−1) = β²/(√(1−β²)(1+√(1−β²))): algebraically exact, and avoids the
        # catastrophic cancellation of computing 1/√(1−β²) − 1 at small β.
        β2 = v2 / c2
        r = sqrt(one(T) - β2)
        s += ps.weight[p] * β2 / (r * (one(T) + r))
    end
    return ps.m * c2 * s
end

function _relativistic_energy_reference(ps, c)
    setprecision(BigFloat, 8192) do
        cw = BigFloat(c)
        sum(eachindex(ps.weight); init = BigFloat(0)) do p
            v2 = sum(BigFloat(ps.v[d][p])^2 for d = 1:3)
            root = sqrt(1 - v2 / cw^2)
            BigFloat(ps.m) * BigFloat(ps.weight[p]) * v2 / (root * (1 + root))
        end
    end
end

@testset "Public complete relativistic energy range" begin
    for T in (Float16, Float32, Float64), D = 1:3
        large = T === Float16 ? T(1e3) : T === Float32 ? T(1e30) : T(1e200)
        for c in (inv(large), large)
            ps = ParticleSet{D,T}(1; m = inv(c))
            foreach(x -> fill!(x, zero(T)), ps.x)
            ps.v[1][1] = c / 2
            before = map(copy, (ps.x..., ps.v..., ps.weight))
            @test kinetic_energy_relativistic(ps, c) === T(_relativistic_energy_reference(ps, c))
            @test all(isequal(a, b) for (a, b) in zip((ps.x..., ps.v..., ps.weight), before))
            @test kinetic_energy_relativistic(ParticleSet{D,T}(0), c) === zero(T)
        end
        ps = ParticleSet{D,T}(1; m = large)
        ps.v[1][1] = T(0.5)
        ps.weight[1] = nextfloat(zero(T))
        @test kinetic_energy_relativistic(ps, one(T)) ===
              T(_relativistic_energy_reference(ps, one(T)))
        ps = ParticleSet{D,T}(1; m = floatmax(T))
        ps.weight[1] = floatmax(T)
        ps.v[1][1] = nextfloat(zero(T))
        @test kinetic_energy_relativistic(ps, one(T)) ===
              T(_relativistic_energy_reference(ps, one(T)))
        ps = ParticleSet{D,T}(16; m = inv(large))
        ps.weight .= floatmax(T) / 2
        ps.v[1] .= T(0.5)
        @test kinetic_energy_relativistic(ps, one(T)) ===
              T(_relativistic_energy_reference(ps, one(T)))
        ps = ParticleSet{D,T}(3; m = floatmax(T))
        ps.weight .= T[1, -0.5, 2]
        ps.v[1] .= T[2.5, 3.75, 3]
        ps.v[2] .= T[2.5, 2.5, 0]
        ps.v[3] .= T[0, 1.25, 0]
        @test kinetic_energy_relativistic(ps, T(5)) === zero(T)
    end
end

@testset "Public ordinary results and nonfinite compatibility" begin
    rng = Xoshiro(98413)
    for T in (Float16, Float32, Float64), D = 1:3, count in (0, 1, 13, 256)
        ps = ParticleSet{D,T}(count; m = T(3))
        foreach(x -> fill!(x, zero(T)), ps.x)
        for v in ps.v
            # Keep component squares normal; the original subnormal fixtures
            # are exercised separately with an independent complete-expression oracle.
            v .= T.(0.125 .+ 0.25rand(rng, count))
        end
        ps.weight .= T.(rand(rng, count) .+ 0.5)
        @test kinetic_energy_relativistic(ps, T(5)) ===
              _relativistic_energy_ordinary_reference(ps, T(5))
    end
    for T in (Float16, Float32, Float64),
        m in T[-Inf, NaN, Inf, 0, -1, 1],
        weight in T[-Inf, NaN, Inf, 0, -1, 1],
        speed in T[0, 0.5, 0.9]

        ps = ParticleSet{1,T}(2; m)
        ps.weight .= T[weight, -1]
        ps.v[1] .= T[speed, 0.5]
        @test isequal(
            kinetic_energy_relativistic(ps, one(T)),
            _relativistic_energy_ordinary_reference(ps, one(T)),
        )
    end
    for T in (Float16, Float32, Float64)
        ps = ParticleSet{1,T}(1)
        ps.v[1][1] = prevfloat(sqrt(nextfloat(zero(T))))
        ps.weight[1] = T(Inf)
        @test kinetic_energy_relativistic(ps, one(T)) === T(Inf)
    end
end

@testset "Public exact speed boundary" begin
    invalid = (
        (Float16, (0.6143, 0.621, -0.4868)),
        (Float32, (-0.86362565, -0.49116498, 0.113612115)),
        (Float64, (-0.474852586251129, 0.8763033506698198, -0.0812862776577999)),
    )
    valid = (
        (Float16, (-0.747, 0.5933, 0.2998)),
        (Float32, (-0.037062123, 0.16837542, -0.98502594)),
        (Float64, (0.9659915805226796, -0.21859548531915907, 0.13811690757244974)),
    )
    for (T, v) in invalid, factor in T[0, 1, Inf, NaN]
        ps = ParticleSet{1,T}(1; m = factor)
        ps.weight[1] = factor
        for d = 1:3
            ps.v[d][1] = v[d]
        end
        @test_throws ArgumentError kinetic_energy_relativistic(ps, one(T))
    end
    for (T, v) in valid
        ps = ParticleSet{1,T}(1)
        for d = 1:3
            ps.v[d][1] = v[d]
        end
        @test kinetic_energy_relativistic(ps, one(T)) ===
              T(_relativistic_energy_reference(ps, one(T)))
    end
    for T in (Float16, Float32, Float64),
        velocity in T[NaN, Inf, -Inf, 1, 2],
        factor in T[0, 1, Inf, NaN]

        ps = ParticleSet{1,T}(1; m = factor)
        ps.weight[1] = factor
        ps.v[1][1] = velocity
        @test_throws ArgumentError kinetic_energy_relativistic(ps, one(T))
    end
end

@testset "Public BigFloat compatibility and retained precision" begin
    for bits in (32, 128, 256), mode in (RoundNearest, RoundDown, RoundUp)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                ps = ParticleSet{1,BigFloat}(2; m = BigFloat(3))
                ps.v[1] .= BigFloat[1, 2]
                ps.v[2] .= BigFloat[0.25, -0.5]
                ps.weight .= BigFloat[2, -0.5]
                @test isequal(
                    kinetic_energy_relativistic(ps, BigFloat(10)),
                    _relativistic_energy_ordinary_reference(ps, BigFloat(10)),
                )
                boundary = ParticleSet{1,BigFloat}(1)
                boundary.v[1][1] = prevfloat(BigFloat(1))
                expected = BigFloat(
                    _relativistic_energy_reference(boundary, BigFloat(1)),
                    mode;
                    precision = bits,
                )
                @test kinetic_energy_relativistic(boundary, BigFloat(1)) == expected
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
            end
        end
    end
    ps, c = setprecision(BigFloat, 256) do
        ps = ParticleSet{1,BigFloat}(1)
        ps.v[1][1] = 1 - BigFloat(2)^(-200)
        ps, BigFloat(1)
    end
    before = copy(ps.v[1])
    for mode in (RoundNearest, RoundDown, RoundUp)
        setprecision(BigFloat, 32) do
            setrounding(BigFloat, mode) do
                @test kinetic_energy_relativistic(ps, c) ==
                      BigFloat(_relativistic_energy_reference(ps, c), mode; precision = 32)
                @test ps.v[1] == before
            end
        end
    end
end

@testset "Public diagnostic views and invalid light speed" begin
    for T in (Float16, Float32, Float64), D = 1:3
        source = ParticleSet{D,T}(5; m = T(2))
        foreach(x -> fill!(x, zero(T)), source.x)
        source.v[1] .= T[0.125, 0.25, 0.375, 0.5, 0.625]
        source.weight .= T[1, -2, 3, -4, 5]
        ps = ParticleSet{D,T}(
            map(a -> view(a, 1:2:5), source.x),
            map(a -> view(a, 1:2:5), source.v),
            view(source.weight, 1:2:5),
            view(source.id, 1:2:5),
            view(source.tag, 1:2:5),
            source.q,
            source.m,
        )
        arrays = (source.x..., source.v..., source.weight, source.id, source.tag)
        before = map(copy, arrays)
        @test kinetic_energy_relativistic(ps, T(4)) ===
              _relativistic_energy_ordinary_reference(ps, T(4))
        large = T === Float16 ? T(1e3) : T === Float32 ? T(1e30) : T(1e200)
        @test kinetic_energy_relativistic(ps, large) ===
              T(_relativistic_energy_reference(ps, large))
        for c in T[0, -1, Inf, NaN]
            @test_throws ArgumentError kinetic_energy_relativistic(ps, c)
        end
        @test all(isequal(a, b) for (a, b) in zip(arrays, before))
    end
end

@testset "Range certificate retains checked energy semantics" begin
    rng = Xoshiro(781956)
    for T in (Float32, Float64), D in (1, 3), count in (0, 1, 13, 256)
        # The independently evaluated precision/exponent bounds give these
        # limits for these population sizes (see the source certificate bounds).
        k = T === Float32 ? 11 : 120
        lower, upper = ldexp(one(T), -k), ldexp(one(T), k)
        for mass in (lower, -lower, upper, -upper, zero(T))
            ps = ParticleSet{D,T}(count; m = mass)
            foreach(x -> fill!(x, zero(T)), ps.x)
            for v in ps.v
                v .= ldexp.(ones(T, count), rand(rng, -k:k-1, count))
            end
            ps.weight .=
                ldexp.(ones(T, count), rand(rng, -k:k, count)) .* rand(rng, T[-1, 0, 1], count)
            @test HybridPlasmaPIC._relativistic_range_certificate(ps, upper)
            checked = HybridPlasmaPIC._kinetic_energy_relativistic(ps, upper, Val(true))
            @test isequal(
                HybridPlasmaPIC._kinetic_energy_relativistic(ps, upper, Val(false)),
                checked,
            )
            @test isequal(kinetic_energy_relativistic(ps, upper), checked)
        end
        ps = ParticleSet{D,T}(max(1, count))
        foreach(x -> fill!(x, zero(T)), ps.x)
        for value in (lower, -lower, upper, -upper, zero(T))
            for field in (ps.v..., ps.weight)
                old = field[1]
                field[1] = value
                @test HybridPlasmaPIC._relativistic_range_certificate(ps, upper)
                field[1] = old
            end
        end
        for value in (
            prevfloat(lower),
            -prevfloat(lower),
            nextfloat(upper),
            -nextfloat(upper),
            T(NaN),
            T(Inf),
        )
            for field in (ps.v..., ps.weight)
                old = field[1]
                field[1] = value
                @test !HybridPlasmaPIC._relativistic_range_certificate(ps, upper)
                field[1] = old
            end
            ps.m = value
            @test !HybridPlasmaPIC._relativistic_range_certificate(ps, upper)
            ps.m = one(T)
        end
        ps.v[1] .= lower
        @test HybridPlasmaPIC._relativistic_range_certificate(ps, 2lower)
        @test kinetic_energy_relativistic(ps, 2lower) ===
              HybridPlasmaPIC._kinetic_energy_relativistic(ps, 2lower, Val(true))
        ps.v[1] .= 2lower
        @test_throws ArgumentError kinetic_energy_relativistic(ps, 2lower)
        ps.v[1] .= prevfloat(2lower)
        @test kinetic_energy_relativistic(ps, 2lower) ===
              HybridPlasmaPIC._kinetic_energy_relativistic_exact(ps, 2lower)
    end
    for T in (Float16, BigFloat)
        @test !HybridPlasmaPIC._relativistic_range_certificate(ParticleSet{1,T}(1), one(T))
    end
end

@testset "One subnormal component selects the complete population energy" begin
    rng = Xoshiro(781957)
    for T in (Float16, Float32, Float64), D = 1:3, count in (3, 257)
        ps = ParticleSet{D,T}(count; m = T(3))
        foreach(x -> fill!(x, zero(T)), ps.x)
        for v in ps.v
            v .= T.(0.125 .+ 0.25rand(rng, count))
        end
        ps.v[1][1] = nextfloat(zero(T))
        ps.weight .= T.(rand(rng, count) .+ 0.5)
        before = map(copy, (ps.v..., ps.weight))
        @test kinetic_energy_relativistic(ps, T(5)) === T(_relativistic_energy_reference(ps, T(5)))
        @test all(isequal(a, b) for (a, b) in zip((ps.v..., ps.weight), before))
    end
end

@testset "Public independent complete finite-range oracle" begin
    rng = Xoshiro(781955)
    for T in (Float16, Float32, Float64),
        D = 1:3,
        count in (1, 3, 11),
        direction in (-1, 1),
        trial = 1:3

        scale = T === Float16 ? 14 : T === Float32 ? 100 : 700
        c = ldexp(one(T), direction * scale)
        ps = ParticleSet{D,T}(count; m = inv(c))
        foreach(x -> fill!(x, zero(T)), ps.x)
        for v in ps.v
            v .= T.(rand(rng, count) .- 0.5) .* c
        end
        ps.weight .= T.(3rand(rng, count) .- 1.5)
        arrays = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
        before = map(copy, arrays)
        @test kinetic_energy_relativistic(ps, c) === T(_relativistic_energy_reference(ps, c))
        @test all(isequal(a, b) for (a, b) in zip(arrays, before))
    end
end
