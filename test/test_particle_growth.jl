module ParticleGrowthTests
using HybridPlasmaPIC, Test, Random

struct SelectedGrowthRNG <: AbstractRNG end
Random.rand(::SelectedGrowthRNG, ::Type{T}) where {T<:AbstractFloat} = zero(T)
Random.randn(::SelectedGrowthRNG, ::Type{Float32}) = 0.0f0
Random.randn(::SelectedGrowthRNG, ::Type{Float64}) = 0.0

struct GrowthVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(a::GrowthVector) = size(a.data)
Base.getindex(a::GrowthVector, i::Int) = a.data[i]
Base.setindex!(a::GrowthVector, v, i::Int) = (a.data[i] = v)
Base.parent(a::GrowthVector) = a.data
Base.resize!(a::GrowthVector, n::Integer) = (resize!(a.data, n); a)
Base.append!(a::GrowthVector, b::AbstractVector) = (append!(a.data, b); a)

struct LengthOnlyVector{T} <: AbstractVector{T}
    n::Int
end
Base.size(a::LengthOnlyVector) = (a.n,)
Base.getindex(::LengthOnlyVector, ::Int) = error("capacity rejection must not read elements")

function species(D, T, N; q = one(T), m = one(T))
    ps = ParticleSet{D,T}(N; q, m)
    for d = 1:D
        ps.x[d] .= T(d) / T(4)
    end
    ps.v[1] .= T(2)
    return ps
end
components(ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
state(ps) = (map(copy, components(ps)), ps.q, ps.m)
function rebuild(ps::ParticleSet{D,T}, arrays) where {D,T}
    return ParticleSet{D,T}(
        ntuple(d -> arrays[d], D),
        ntuple(c -> arrays[D+c], 3),
        arrays[D+4],
        arrays[D+5],
        arrays[D+6],
        ps.q,
        ps.m,
    )
end
function view_component(ps::ParticleSet{D,T}, index) where {D,T}
    a = components(ps)
    arrays = ntuple(k -> k == index ? view(copy(a[k]), :) : a[k], length(a))
    # The convenience constructor requires uniform position/velocity storage
    # types. Explicit abstract component fields let each late slot be tested.
    x = ntuple(d -> arrays[d], D)
    v = ntuple(c -> arrays[D+c], 3)
    weight, id, tag = arrays[D+4], arrays[D+5], arrays[D+6]
    @test all(axes(a, 1) == Base.OneTo(nparticles(ps)) for a in arrays)
    return ParticleSet{
        D,
        T,
        AbstractVector{T},
        AbstractVector{T},
        typeof(weight),
        typeof(id),
        typeof(tag),
    }(
        x,
        v,
        weight,
        id,
        tag,
        ps.q,
        ps.m,
    )
end

function rejected_unchanged(dest, src)
    before = state(dest), state(src)
    @test_throws ArgumentError append_particles!(dest, src)
    @test isequal((state(dest), state(src)), before)
end

@testset "GROWTH-001 append validates the complete storage batch" begin
    for T in (Float32, Float64), D = 1:3
        for k = 1:(D+6)
            dest = view_component(species(D, T, 2), k)
            rejected_unchanged(dest, species(D, T, 1))
        end
        # Every floating-point destination pair, including empty shared storage.
        for N in (0, 2), i = 1:(D+4), j = (i+1):(D+4)
            base = species(D, T, N)
            a = components(base)
            dest = rebuild(base, ntuple(k -> k == j ? a[i] : a[k], length(a)))
            rejected_unchanged(dest, species(D, T, 1))
        end
        # Every same-element-type component pairing across distinct populations.
        for i = 1:(D+4), j = 1:(D+4)
            dest, base = species(D, T, 2), species(D, T, 2)
            a = components(base)
            src = rebuild(base, ntuple(k -> k == j ? components(dest)[i] : a[k], length(a)))
            rejected_unchanged(dest, src)
        end
        for k = 1:(D+6), side in (:source, :destination)
            dest, src = species(D, T, 2), species(D, T, 1)
            a = components(side == :source ? src : dest)[k]
            resize!(a, length(a) + 1)
            fill!(a, zero(eltype(a)))
            rejected_unchanged(dest, src)
        end
        dest = species(D, T, 2)
        huge = ParticleSet{D,T}(
            ntuple(_ -> LengthOnlyVector{T}(typemax(Int)), D),
            ntuple(_ -> LengthOnlyVector{T}(typemax(Int)), 3),
            LengthOnlyVector{T}(typemax(Int)),
            LengthOnlyVector{UInt64}(typemax(Int)),
            LengthOnlyVector{UInt32}(typemax(Int)),
            one(T),
            one(T),
        )
        before = state(dest)
        @test_throws ArgumentError append_particles!(dest, huge)
        @test isequal(state(dest), before)
    end
end

@testset "Growth sees backing arrays and mutable view indices" begin
    for T in (Float32, Float64), N in (0, 2)
        dest = species(1, T, N)
        base = species(1, T, 1)
        src = ParticleSet{1,T}(
            (view(dest.x[1], 1:0),),
            (T[], T[], T[]),
            T[],
            UInt64[],
            UInt32[],
            one(T),
            one(T),
        )
        # Empty-source calls are unchanged no-ops, even with a shared empty view.
        before = state(dest), state(src)
        @test append_particles!(dest, src) === dest
        @test isequal((state(dest), state(src)), before)
    end
    for T in (Float32, Float64)
        dest = species(1, T, 2)
        base = species(1, T, 2)
        indexed = view(T[1, 2, 3, 4], dest.id)
        src = ParticleSet{1,T}((indexed,), base.v, base.weight, base.id, base.tag, base.q, base.m)
        rejected_unchanged(dest, src)
        # Wrappers expose a shared resizable parent: reject before growth.
        raw = T[1, 2, 3, 4]
        dest = rebuild(
            species(1, T, 2),
            (
                GrowthVector(raw),
                GrowthVector(T[0, 0, 0, 0]),
                GrowthVector(T[0, 0, 0, 0]),
                GrowthVector(T[0, 0, 0, 0]),
                GrowthVector(T[1, 1, 1, 1]),
                UInt64[1, 2, 3, 4],
                UInt32[0, 0, 0, 0],
            ),
        )
        base = species(1, T, 2)
        src = ParticleSet{1,T}(
            (view(raw, 1:2),),
            base.v,
            base.weight,
            base.id,
            base.tag,
            one(T),
            one(T),
        )
        rejected_unchanged(dest, src)
    end
end

@testset "Growth preserves independent, self and read-only-source appends" begin
    for T in (Float32, Float64), D = 1:3, wrapped in (false, true)
        dest, src = species(D, T, 2), species(D, T, 1)
        if wrapped
            dest = rebuild(dest, map(a -> GrowthVector(a), components(dest)))
        end
        before, source_before = state(dest), state(src)
        @test append_particles!(dest, src) === dest
        @test map(collect, components(dest)) == map(vcat, before[1], source_before[1])
        @test isequal(state(src), source_before)
        before = state(dest)
        @test append_particles!(dest, dest) === dest
        @test map(collect, components(dest)) == map(a -> vcat(a, a), before[1])

        base = species(D, T, 2)
        a = components(base)
        src = rebuild(base, ntuple(k -> k == D + 4 ? a[1] : a[k], length(a)))
        dest = species(D, T, 0)
        before = state(src)
        @test append_particles!(dest, src) === dest
        @test isequal(state(dest), before)
        @test isequal(state(src), before)
        dest = view_component(species(D, T, 2), D + 6)
        before = state(dest)
        @test append_particles!(dest, species(D, T, 0)) === dest
        @test isequal(state(dest), before)

        src = rebuild(base, map(a -> view(copy(a), :), components(base)))
        dest = species(D, T, 0)
        before = state(src)
        @test append_particles!(dest, src) === dest
        @test isequal(state(dest), before)
        @test isequal(state(src), before)
    end
end

@testset "Ionization validates both destinations before cooling or growth" begin
    for T in (Float32, Float64), D = 1:3, side in (:electron, :ion), kind in (:view, :alias)
        electrons, ions = species(D, T, 1; q = -one(T)), species(D, T, 0; m = T(100))
        function invalid(ps)
            kind == :view && return view_component(ps, D + 6)
            a = components(ps)
            return rebuild(ps, ntuple(k -> k == D + 4 ? a[1] : a[k], length(a)))
        end
        if side == :electron
            electrons = invalid(electrons)
        else
            ions = invalid(ions)
        end
        e_nextid, i_nextid = Ref(UInt64(2)), Ref(UInt64(1))
        before = state(electrons), state(ions)
        @test_throws ArgumentError ionize_mcc!(
            electrons,
            ions,
            one(T);
            nσ_iz = T(100),
            E_iz = T(0.5),
            rng = SelectedGrowthRNG(),
            e_nextid,
            i_nextid,
        )
        @test isequal((state(electrons), state(ions)), before)
        @test (e_nextid[], i_nextid[]) == (UInt64(2), UInt64(1))
        @test ionize_mcc!(
            electrons,
            ions,
            zero(T);
            nσ_iz = T(100),
            E_iz = T(0.5),
            rng = SelectedGrowthRNG(),
        ) == 0
        @test isequal((state(electrons), state(ions)), before)
    end
    for T in (Float32, Float64), D = 1:3
        electrons, ions = species(D, T, 1; q = -one(T)), species(D, T, 0; m = T(100))
        e_nextid, i_nextid = Ref(UInt64(2)), Ref(UInt64(1))
        @test ionize_mcc!(
            electrons,
            ions,
            one(T);
            nσ_iz = T(100),
            E_iz = T(0.5),
            rng = SelectedGrowthRNG(),
            e_nextid,
            i_nextid,
        ) == 1
        @test (nparticles(electrons), nparticles(ions)) == (2, 1)
        @test electrons.v[1][1] ≈ sqrt(T(3))
        @test electrons.id == UInt64[1, 2]
        @test ions.id == UInt64[1]
        @test (e_nextid[], i_nextid[]) == (UInt64(3), UInt64(2))
    end
end

@testset "Injection storage rejection preserves particle and bookkeeping state" begin
    for T in (Float32, Float64), k = 1:7
        ps = view_component(species(1, T, 2), k)
        before = state(ps)
        acc, nextid = Ref(0.25), Ref(UInt64(3))
        @test_throws ArgumentError inject_face_1d!(
            ps,
            SelectedGrowthRNG(),
            zero(T),
            1,
            one(T),
            one(T),
            zero(T),
            (zero(T), zero(T)),
            zero(T),
            one(T),
            one(T),
            acc,
            nextid,
        )
        @test isequal(state(ps), before)
        @test (acc[], nextid[]) == (0.25, UInt64(3))
    end
end
end
