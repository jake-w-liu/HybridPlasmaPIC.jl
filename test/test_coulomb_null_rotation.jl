module TestCoulombNullRotation
using HybridPlasmaPIC, Random, Test

mutable struct ScriptedRNG <: AbstractRNG
    normal::Float64
    uniforms::Vector{Float64}
    draws::Int
    events::Vector{Symbol}
end
ScriptedRNG(normal = 0.0; uniforms = [0.125]) = ScriptedRNG(normal, uniforms, 0, Symbol[])
function Random.shuffle!(rng::ScriptedRNG, values::AbstractVector)
    push!(rng.events, :shuffle)
    return values
end
function scripted_normal(rng::ScriptedRNG, ::Type{T}) where {T<:AbstractFloat}
    push!(rng.events, :normal)
    return T(rng.normal)
end
function scripted_uniform(rng::ScriptedRNG, ::Type{T}) where {T<:AbstractFloat}
    push!(rng.events, :uniform)
    rng.draws += 1
    return T(rng.uniforms[mod1(rng.draws, length(rng.uniforms))])
end

# Exact type methods avoid intersections with Random's hardware-float methods.
for T in (Float16, Float32, Float64, BigFloat)
    @eval Random.randn(rng::ScriptedRNG, ::Type{$T}) = scripted_normal(rng, $T)
    @eval Random.rand(rng::ScriptedRNG, ::Type{$T}) = scripted_uniform(rng, $T)
end

function pair_indices(n)
    pairs = isodd(n) ? [(1, 2), (2, 3), (3, 1)] : Tuple{Int,Int}[]
    first = isodd(n) ? 4 : 1
    append!(pairs, [(i, i + 1) for i = first:2:n-1])
    return pairs
end

@testset "Exact Coulomb identity preserves original particle bits and draws" begin
    for T in (Float16, Float32, Float64, BigFloat),
        D in (1, 2, 3),
        n in (2, 3, 4, 5),
        oblique in (false, true),
        axis = 1:3,
        negative in (false, true),
        weighting in (:equal, :unequal, :zero_partner)

        ps = ParticleSet{D,T}(n)
        for x in ps.x
            fill!(x, T(0.25))
        end
        tiny = T === BigFloat ? ldexp(one(T), -3precision(T)) : nextfloat(zero(T))
        sign = negative ? -one(T) : one(T)
        for p = 1:n, c = 1:3
            component = mod1(c + axis - 1, 3)
            ps.v[component][p] = (oblique || c == 1) ? sign * T(c * (p - 1)) : -zero(T)
        end
        ps.v[axis][1] = sign * tiny
        weighting === :unequal && (ps.weight .= T.(1:n))
        weighting === :zero_partner && (ps.weight[1] = zero(T))
        before = deepcopy(ps)
        rng = ScriptedRNG(negative ? -0.0 : 0.0)
        work = map(similar, ps.v)
        indices = Vector{Int}(undef, n)
        @test collide_coulomb!(ps, T(0.001), T(0.01); rng, work, index_work = indices) === ps
        @test isequal(ps.v, before.v)
        @test isequal(work, before.v)
        @test isequal(ps.x, before.x)
        @test ps.weight == before.weight && ps.id == before.id && ps.tag == before.tag
        expected = Symbol[:shuffle]
        for (i, j) in pair_indices(n)
            append!(expected, [:normal, :uniform])
            ps.weight[i] != ps.weight[j] && append!(expected, [:uniform, :uniform])
        end
        @test rng.events == expected
    end
end

@testset "Null angle retains invalid draw failure and rollback" begin
    for T in (Float32, Float64),
        (normal, uniforms) in ((NaN, [0.125]), (0.0, [NaN]), (0.0, [0.125, NaN]))

        ps = ParticleSet{1,T}(4)
        fill!(ps.x[1], T(0.25))
        ps.v[1] .= T[0, 1, 2, 3]
        before = deepcopy(ps)
        rng = ScriptedRNG(normal; uniforms)
        @test_throws ArgumentError collide_coulomb!(
            ps,
            T(0.001),
            T(0.01);
            rng,
            index_work = Vector{Int}(undef, 4),
        )
        @test isequal(ps.v, before.v)
        @test isequal(ps.x, before.x)
    end
end

@testset "Zero sine at pi is a swap, and small nonzero angles still rotate" begin
    for T in (Float32, Float64)
        ps = ParticleSet{1,T}(2)
        fill!(ps.x[1], T(0.25))
        ps.v[1] .= T[0, 1]
        # Large variance selects isotropic sampling; a uniform zero gives cosθ=-1.
        rng = ScriptedRNG(; uniforms = [0.0])
        collide_coulomb!(ps, T(10), one(T); rng, index_work = Vector{Int}(undef, 2))
        @test ps.v[1] == T[1, 0]
        @test all(iszero, ps.v[2]) && all(iszero, ps.v[3])
        @test rng.events == [:shuffle, :uniform, :uniform]

        for normal in (Float64(eps(T)), -Float64(eps(T)))
            ps.v[1] .= T[0, 1]
            fill!(ps.v[2], zero(T))
            fill!(ps.v[3], zero(T))
            rng = ScriptedRNG(normal; uniforms = [0.0])
            collide_coulomb!(ps, T(0.25), one(T); rng, index_work = Vector{Int}(undef, 2))
            # With g=(-1,0,0), φ=0, the y kick vanishes and z kick has normal's sign.
            @test !iszero(ps.v[3][1]) && signbit(ps.v[3][1]) != signbit(T(normal))
            @test ps.v[3][2] == -ps.v[3][1]
            @test all(iszero, ps.v[2])
        end
    end
end
end
