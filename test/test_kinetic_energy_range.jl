using HybridPlasmaPIC, Test, Random

module KineticRangeOracle
using HybridPlasmaPIC
Q = Rational{BigInt}
function exact(species)
    result = Q(0)
    for p in species, i in eachindex(p.weight)
        result += Q(p.m) * Q(p.weight[i]) * sum(Q(v[i])^2 for v in p.v) / 2
    end
    return result
end
function ordinary(p::ParticleSet{D,T}) where {D,T}
    total = zero(T)
    for i in eachindex(p.weight)
        total += p.weight[i] * (p.v[1][i]^2 + p.v[2][i]^2 + p.v[3][i]^2)
    end
    return T(0.5) * p.m * total
end
function rounded(x::T, value) where {T}
    isfinite(x) || isinf(x) || return false
    largest = Q(floatmax(T))
    edge = largest + (largest - Q(prevfloat(floatmax(T)))) / 2
    if isinf(x)
        return signbit(x) ? value <= -edge : value >= edge
    end
    value = signbit(x) ? -value : value
    x = abs(x)
    if iszero(x)
        return 0 <= value <= Q(nextfloat(zero(T))) / 2
    end
    low = (Q(prevfloat(x)) + Q(x)) / 2
    high = isfinite(nextfloat(x)) ? (Q(x) + Q(nextfloat(x))) / 2 : edge
    U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
    even = iseven(reinterpret(U, x))
    return low < value < high || even && (low == value || high == value)
end
end

@testset "Kinetic complete finite range, populations and views" begin
    oracle = KineticRangeOracle
    for T in (Float16, Float32, Float64), D = 1:3
        large = T === Float16 ? T(1000) : T === Float32 ? T(1e25) : T(1e200)
        tiny = T === Float16 ? T(1e-4) : T === Float32 ? T(1e-30) : T(1e-300)
        for (mass, weight, speed) in (
            (tiny, one(T), large),
            (floatmax(T), floatmax(T), nextfloat(zero(T))),
            (nextfloat(zero(T)), floatmax(T), one(T)),
            (zero(T), one(T), large),
            (-tiny, one(T), large),
        )
            ps = ParticleSet{D,T}(2; m = mass)
            ps.weight .= (weight, -weight / T(2))
            ps.v[1] .= speed
            before = map(copy, (ps.v..., ps.weight))
            ref = oracle.exact((ps,))
            @test oracle.rounded(kinetic_energy(ps), ref)
            @test oracle.rounded(kinetic_energy([ps]), ref)
            viewed = ParticleSet{D,T}(
                map(a -> view(a, :), ps.x),
                map(a -> view(a, :), ps.v),
                view(ps.weight, :),
                view(ps.id, :),
                view(ps.tag, :),
                ps.q,
                ps.m,
            )
            @test isequal(kinetic_energy(viewed), kinetic_energy(ps))
            @test all(isequal(a, b) for (a, b) in zip((ps.v..., ps.weight), before))
        end
        a, b = [ParticleSet{D,T}(1; m = floatmax(T)) for _ = 1:2]
        a.v[1][1] = b.v[1][1] = T(2)
        b.weight[1] = -one(T)
        c = ParticleSet{D,T}(1; m = T(2))
        c.v[1][1] = one(T)
        @test kinetic_energy([a, b]) === zero(T)
        @test kinetic_energy([a, b, c]) === one(T)
        @test kinetic_energy(view([a, b, c], :)) === one(T)
        @test oracle.exact([a, b, c]) == 1
        @test isinf(kinetic_energy(a)) && isinf(kinetic_energy(b))

        quarter = (floatmax(T) - prevfloat(floatmax(T))) / T(4)
        for tail in (prevfloat(quarter), quarter, nextfloat(quarter))
            ps = ParticleSet{D,T}(3; m = T(2))
            ps.weight .= (floatmax(T), quarter, tail)
            ps.v[1] .= one(T)
            ref = oracle.exact((ps,))
            @test oracle.rounded(kinetic_energy(ps), ref)
            parts = [ParticleSet{D,T}(1; m = T(2)) for _ = 1:3]
            for (p, w) in zip(parts, ps.weight)
                p.weight[1] = w
                p.v[1][1] = one(T)
            end
            @test oracle.rounded(kinetic_energy(parts), ref)
        end
    end
end

function kinetic_range_alloc(ps)
    kinetic_energy(ps)
    return @allocated kinetic_energy(ps)
end

@testset "Kinetic ordinary arithmetic and nonfinite compatibility" begin
    oracle = KineticRangeOracle
    rng = MersenneTwister(78192)
    for T in (Float16, Float32, Float64), D = 1:3
        # Small dyadic factors are inside the positive certificate even at
        # Float16 precision; the original reduction order remains observable.
        ps = ParticleSet{D,T}(4; m = one(T))
        for v in ps.v
            v .= T.(rand(rng, (0.5, 1, 1.5), 4))
        end
        @test isequal(kinetic_energy(ps), oracle.ordinary(ps))
        @test kinetic_range_alloc(ps) == 0
        @test kinetic_range_alloc([ps]) == 0
        @test kinetic_energy(ParticleSet{D,T}[]) === zero(T)
        for m in (-floatmax(T), -one(T), -zero(T), zero(T), floatmax(T), T(NaN), T(Inf))
            empty = ParticleSet{D,T}(0; m)
            @test isequal(kinetic_energy(empty), oracle.ordinary(empty))
        end
        for m in T[-Inf, Inf, NaN, -1, 0, 1],
            w in T[-Inf, Inf, NaN, -1, 0, 1],
            speed in T[0, 0.5, Inf, NaN]

            p = ParticleSet{D,T}(1; m)
            p.weight[1] = w
            p.v[1][1] = speed
            if !all(isfinite, (m, w, speed))
                @test isequal(kinetic_energy(p), oracle.ordinary(p))
                @test isequal(kinetic_energy([p]), zero(T) + oracle.ordinary(p))
            end
        end
    end
end

@testset "Kinetic exact subnormal rounding boundaries" begin
    oracle = KineticRangeOracle
    for T in (Float16, Float32, Float64), sign in (-one(T), one(T)), count = 1:4
        ps = ParticleSet{1,T}(count)
        ps.weight .= sign * nextfloat(zero(T))
        ps.v[1] .= one(T)
        exact = oracle.exact((ps,))
        @test oracle.rounded(kinetic_energy(ps), exact)
        @test isequal(kinetic_energy(ps), T(exact))
        @test oracle.rounded(kinetic_energy([ps]), exact)
    end
end

@testset "Kinetic range certificates against exact quadratic sums" begin
    oracle = KineticRangeOracle
    Q = Rational{BigInt}
    rng = MersenneTwister(398512)
    for T in (Float16, Float32, Float64), count in (1, 2, 4, 16, 128), mixed in (false, true)
        bounds = HybridPlasmaPIC._kinetic_energy_bounds(T, Int128(count))
        lower, upper = mixed ? bounds[3:4] : bounds[1:2]
        isfinite(lower) || continue
        for trial = 1:8
            ps = ParticleSet{1,T}(count)
            values = T[lower, upper, one(T), prevfloat(upper)]
            filter!(x -> lower <= x <= upper, values)
            ps.m = rand(rng, values)
            for v in ps.v
                v .= rand(rng, values, count)
            end
            ps.weight .= rand(rng, values, count)
            mixed && (ps.weight[1:2:end] .*= -one(T))
            raw, positive_ok, mixed_ok, plus, minus =
                HybridPlasmaPIC._kinetic_energy_probe(ps, bounds, Val(true))
            @test mixed ? mixed_ok : positive_ok && !(plus && minus)
            exact = oracle.exact((ps,))
            @test iszero(exact) || abs(exact) >= Q(floatmin(T))
            @test abs(exact) < Q(floatmax(T))
            @test isequal(kinetic_energy(ps), oracle.ordinary(ps))
            # Independently replay every original intermediate, checking the
            # range contract rather than reproducing the exponent inequalities.
            total = zero(T)
            for i in eachindex(ps.weight)
                x2, y2, z2 = (v[i]^2 for v in ps.v)
                xy = x2 + y2
                norm = xy + z2
                term = ps.weight[i] * norm
                next = total + term
                @test all(
                    x -> isfinite(x) && (iszero(x) || abs(x) >= floatmin(T)),
                    (x2, y2, z2, xy, norm, term, next),
                )
                total = next
            end
            @test isfinite(raw) && (iszero(raw) || abs(raw) >= floatmin(T))
            @test isfinite(T(0.5) * ps.m)
        end
    end
end
