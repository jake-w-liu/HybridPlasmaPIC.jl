module IonizationGrowthAliasTests
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

# Parent metadata is available even when a custom vector has not opted into
# Base.dataids. Distinct mutable wrappers can therefore hide shared backing
# from Base.mightalias while still exposing it to the particle growth guard.
mutable struct ParentGrowthVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(a::ParentGrowthVector) = size(a.data)
Base.getindex(a::ParentGrowthVector, i::Int) = a.data[i]
Base.setindex!(a::ParentGrowthVector, value, i::Int) = (a.data[i] = value; a)
Base.parent(a::ParentGrowthVector) = a.data
Base.resize!(a::ParentGrowthVector, n::Integer) = (resize!(a.data, n); a)
Base.append!(a::ParentGrowthVector, b::AbstractVector) = (append!(a.data, b); a)

struct SelectBirthRNG <: AbstractRNG end
Random.rand(::SelectBirthRNG, ::Type{T}) where {T<:AbstractFloat} = zero(T)
Random.randn(::SelectBirthRNG, ::Type{Float32}) = 0.0f0
Random.randn(::SelectBirthRNG, ::Type{Float64}) = 0.0

components(p) = (p.x..., p.v..., p.weight, p.id, p.tag)
state(p) = (map(copy, components(p)), p.q, p.m)
function rebuild(p::ParticleSet{D,T}, a) where {D,T}
    return ParticleSet{D,T}(
        ntuple(d -> a[d], D),
        ntuple(c -> a[D+c], 3),
        a[D+4],
        a[D+5],
        a[D+6],
        p.q,
        p.m,
    )
end
function populations(D, T)
    e = ParticleSet{D,T}(1; q = -one(T), m = one(T))
    i = ParticleSet{D,T}(1; q = one(T), m = T(100))
    for x in e.x
        fill!(x, T(0.25))
    end
    for x in i.x
        fill!(x, T(0.75))
    end
    e.v[1][1] = T(2)
    return rebuild(e, map(ParentGrowthVector, components(e))),
    rebuild(i, map(ParentGrowthVector, components(i)))
end
function birth!(e, i, dt, ec, ic)
    T = eltype(e)
    return ionize_mcc!(
        e,
        i,
        dt;
        nσ_iz = one(T),
        E_iz = T(0.5),
        rng = SelectBirthRNG(),
        e_nextid = ec,
        i_nextid = ic,
    )
end

@testset "Ionization preserves independent resizable parent wrappers" begin
    for T in (Float32, Float64), D = 1:3
        @test rand(SelectBirthRNG(), T) == zero(T)
        @test randn(SelectBirthRNG(), T) == zero(T)
        e, i = populations(D, T)
        ec, ic = Ref(UInt64(2)), Ref(UInt64(2))
        @test birth!(e, i, T(0.1), ec, ic) == 1
        @test (nparticles(e), nparticles(i)) == (2, 2)
        @test e.id == i.id == UInt64[1, 2]
        @test (ec[], ic[]) == (UInt64(3), UInt64(3))
        @test all(axes(a, 1) == Base.OneTo(2) for p in (e, i) for a in components(p))
    end
end

@testset "Ionization rejects every cross-species shared backing before mutation" begin
    for T in (Float32, Float64), D = 1:3
        # Float fields can share with every other float field. IDs and tags have
        # distinct element types, so each shares only with its matching field.
        pairs = [(a, b) for a = 1:(D+4) for b = 1:(D+4)]
        append!(pairs, [(D + 5, D + 5), (D + 6, D + 6)])
        for (a, b) in pairs
            e, i = populations(D, T)
            shared = ParentGrowthVector(parent(components(e)[a]))
            i = rebuild(i, Base.setindex(components(i), shared, b))
            @test components(e)[a] !== components(i)[b]
            @test parent(components(e)[a]) === parent(components(i)[b])
            @test !H._particle_sets_mightalias(e, i)
            before = state(e), state(i)
            ec, ic = Ref(UInt64(2)), Ref(UInt64(2))
            @test birth!(e, i, zero(T), ec, ic) == 0
            @test isequal((state(e), state(i)), before)
            @test (ec[], ic[]) == (UInt64(2), UInt64(2))
            err = try
                birth!(e, i, T(0.1), ec, ic)
                nothing
            catch value
                value
            end
            @test err isa ArgumentError
            @test err isa ArgumentError &&
                  err.msg == "ionize_mcc!: electron and ion particle storage must not alias"
            @test isequal((state(e), state(i)), before)
            @test (ec[], ic[]) == (UInt64(2), UInt64(2))
        end
    end
end
end
