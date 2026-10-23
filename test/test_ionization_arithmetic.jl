module IonizationArithmeticTests
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

mutable struct MidpointRNG <: AbstractRNG
    uniform::Int
    normal::Int
end
MidpointRNG() = MidpointRNG(0, 0)
function uniform(rng, T)
    rng.uniform += 1
    return T(0.5)
end
function normal(rng, T)
    rng.normal += 1
    return zero(T)
end
for T in (Float16, Float32, Float64, BigFloat)
    @eval Random.rand(rng::MidpointRNG, ::Type{$T}) = uniform(rng, $T)
    @eval Random.randn(rng::MidpointRNG, ::Type{$T}) = normal(rng, $T)
end

# This oracle reads stored significands, independently of production conversion.
function exact(x::AbstractFloat)
    significand, power, sign = Base.decompose(x)
    n = BigInt(sign) * significand
    return power >= 0 ? (n << power) // BigInt(1) : n // (BigInt(1) << -power)
end
energy(v, mass) = exact(mass) * sum(x -> exact(x)^2, v) / 2
state(p) = (map(copy, (p.x..., p.v..., p.weight, p.id, p.tag)), p.q, p.m)

function oracle(v::NTuple{3,T}, mass, loss) where {T}
    e = energy(v, mass)
    # 8192 bits resolve the deliberately wide Float64 exponent cases. The
    # square-root primitive is tested separately against exact rounding cells.
    return setprecision(BigFloat, 8192) do
        setrounding(BigFloat, RoundNearest) do
            ntuple(3) do c
                squared = exact(v[c])^2 * (e - exact(loss)) / e
                copysign(T(sqrt(BigFloat(squared))), v[c])
            end
        end
    end
end

function exercise(D, v::NTuple{3,T}, mass, loss; rate = one(T), dt = T(16)) where {T}
    e = ParticleSet{D,T}(1; q = -one(T), m = mass)
    i = ParticleSet{D,T}(1; q = one(T), m = T(100))
    for d = 1:D
        e.x[d][1] = T(d) / T(8)
        i.x[d][1] = T(d) / T(4)
    end
    for c = 1:3
        e.v[c][1] = v[c]
    end
    e.weight[1] = T(2)
    e.id[1], i.id[1] = UInt64(7), UInt64(11)
    e.tag[1], i.tag[1] = 3, 4
    before = state(e), state(i)
    ec, ic = Ref(UInt64(20)), Ref(UInt64(30))
    rng = MidpointRNG()
    eligible = energy(v, mass) > exact(loss)
    # Evaluate probability once at high precision before invoking production.
    depth = setprecision(BigFloat, 8192) do
        BigFloat(exact(rate) * exact(dt)) * sqrt(BigFloat(sum(x -> exact(x)^2, v)))
    end
    selected = eligible && depth > log(BigFloat(2))
    expected = selected ? oracle(v, mass, loss) : v
    count =
        ionize_mcc!(e, i, dt; nσ_iz = rate, E_iz = loss, e_nextid = ec, i_nextid = ic, rng = rng)
    @test count == Int(selected)
    @test (nparticles(e), nparticles(i)) == (1 + count, 1 + count)
    @test (rng.uniform, rng.normal) ==
          (Int(eligible && !iszero(rate) && !iszero(dt)), 6 * Int(selected))
    @test (ec[], ic[]) == (UInt64(20 + Int(selected)), UInt64(30 + Int(selected)))
    # Keep a failed birth-count assertion from causing unrelated bounds errors.
    count == Int(selected) || return nothing
    if !selected
        @test isequal((state(e), state(i)), before)
    else
        actual = ntuple(c -> e.v[c][1], 3)
        @test all(isfinite, actual)
        for c = 1:3
            @test signbit(actual[c]) == signbit(v[c])
            @test abs(exact(actual[c]) - exact(expected[c])) <=
                  8 * exact(eps(T)) * abs(exact(expected[c]))
        end
        @test e.id == UInt64[7, 20]
        @test i.id == UInt64[11, 30]
        @test e.weight == T[2, 2]
        @test i.weight == T[1, 2]
        @test all(e.x[d][2] == i.x[d][2] == e.x[d][1] for d = 1:D)
        @test all(iszero(e.v[c][2]) && iszero(i.v[c][2]) for c = 1:3)
        @test all(isequal(i.v[c][1], before[2][1][D+c][1]) for c = 1:3)
        @test e.tag[1] == 3 && i.tag[1] == 4
    end
    return nothing
end

@testset "Ionization uses complete stored energy and optical depth" begin
    for T in (Float32, Float64), D = 1:3
        maximum = floatmax(T)
        huge = (maximum, maximum, maximum)
        exercise(D, huge, floatmin(T) / T(4), maximum)
        exercise(D, huge, floatmin(T), maximum)
        exercise(D, huge, floatmin(T), maximum; rate = floatmin(T) / T(128), dt = one(T))
        for speed in (prevfloat(one(T)), one(T), nextfloat(one(T)), T(2))
            exercise(D, (speed, zero(T), -zero(T)), one(T), T(0.5))
        end
        tiny = ldexp(one(T), T === Float32 ? -80 : -600)
        exercise(D, (one(T), tiny, zero(T)), one(T), T(0.5))
        k = T === Float32 ? 100 : 1000
        wide = (ldexp(one(T), k), ldexp(one(T), -k), zero(T))
        exercise(D, wide, ldexp(one(T), -k), ldexp(one(T), k - 1))
        for v in ((T(2), T(3), T(4)), huge), loss in (zero(T), T(0.5))
            exercise(D, v, one(T), loss)
        end
        exercise(D, (T(2), zero(T), zero(T)), one(T), T(0.5); rate = zero(T))
        exercise(D, (T(2), zero(T), zero(T)), one(T), T(0.5); dt = zero(T))
    end

    @testset "Sign, permutation and exponent variants" begin
        rng = MersenneTwister(0x19e)
        for T in (Float32, Float64), D = 1:3, trial = 1:24
            power = rand(rng, T === Float32 ? (-40:40) : (-400:400))
            components = (ldexp(one(T), power), ldexp(T(0.5), power), ldexp(T(0.25), power))
            order = randperm(rng, 3)
            v = ntuple(c -> (rand(rng, Bool) ? one(T) : -one(T)) * components[order[c]], 3)
            mass = ldexp(one(T), -power)
            loss = ldexp(T(0.5), power)
            exercise(D, v, mass, loss; dt = ldexp(T(16), -power))
        end
    end

    @testset "Other native types and caller contexts" begin
        for D = 1:3
            exercise(D, (Float16(1), nextfloat(Float16(0)), Float16(0)), Float16(1), Float16(0.5))
            exercise(D, (Float16(2), Float16(3), Float16(4)), Float16(1), Float16(0.5))
        end
        stored_velocity = setprecision(BigFloat, 512) do
            (BigFloat(1) + ldexp(BigFloat(1), -400), BigFloat(0), BigFloat(0))
        end
        for bits in (96, 256), mode in (RoundNearest, RoundUp, RoundDown)
            setprecision(BigFloat, bits) do
                setrounding(BigFloat, mode) do
                    exercise(
                        1,
                        (nextfloat(BigFloat(1)), BigFloat(0), BigFloat(0)),
                        BigFloat(1),
                        BigFloat(0.5),
                    )
                    exercise(1, stored_velocity, BigFloat(1), BigFloat(0.5))
                    @test precision(BigFloat) == bits
                    @test rounding(BigFloat) == mode
                end
            end
        end
    end

    @testset "Concurrent threshold decisions" begin
        # Build all oracle results serially; production IEEE paths must not
        # change the process-wide BigFloat precision or rounding context.
        for T in (Float32, Float64)
            v = (nextfloat(one(T)), zero(T), zero(T))
            expected = oracle(v, one(T), T(0.5))
            outputs = Vector{NTuple{3,T}}(undef, 32)
            counts = Vector{Int}(undef, 32)
            Threads.@threads for n = 1:32
                e = ParticleSet{1,T}(1; q = -one(T), m = one(T))
                i = ParticleSet{1,T}(0; q = one(T), m = T(100))
                e.x[1][1] = zero(T)
                for c = 1:3
                    e.v[c][1] = v[c]
                end
                counts[n] =
                    ionize_mcc!(e, i, T(16); nσ_iz = one(T), E_iz = T(0.5), rng = MidpointRNG())
                outputs[n] = ntuple(c -> e.v[c][1], 3)
            end
            @test all(==(1), counts)
            @test all(==(expected), outputs)
        end
    end
end
end
