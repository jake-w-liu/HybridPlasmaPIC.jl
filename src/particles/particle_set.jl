# particles.jl — structure-of-arrays particle storage, loading, the Boris mover,
# and particle boundaries. Positions carry exactly D spatial coordinates; the
# velocity always has all three components (D-dimensional space, 3-D velocity).
#
# Units are normalized: the mover takes q/m directly, so the caller supplies
# whatever normalization (the hybrid model uses Ω_ci-normalized units).

"""
    ParticleSet{D,T}(N; q=1, m=1)

Structure-of-arrays storage for `N` particles of one species in `D` spatial
dimensions. Fields: `x` (D position vectors), `v` (3 velocity vectors),
`weight`, `id`, `tag`, and species charge `q` / mass `m`.
"""
mutable struct ParticleSet{
    D,
    T,
    X<:AbstractVector{T},
    V<:AbstractVector{T},
    W<:AbstractVector{T},
    I<:AbstractVector{UInt64},
    G<:AbstractVector{UInt32},
}
    x::NTuple{D,X}
    v::NTuple{3,V}
    weight::W
    id::I
    tag::G      # provenance/flags; widen when shock diagnostics need more
    q::T
    m::T
end

function _particle_length(N::Integer)
    N >= 0 || throw(ArgumentError("N must be ≥ 0"))
    N <= typemax(Int) || throw(ArgumentError("N must fit in Int"))
    return Int(N)
end

function _particle_ids(N::Int)
    ids = Vector{UInt64}(undef, N)
    @inbounds for i in eachindex(ids)
        ids[i] = UInt64(i)
    end
    return ids
end

function _check_particle_vector_axes(name::Symbol, a::AbstractVector, N::Int)
    axes(a, 1) == Base.OneTo(N) || throw(
        ArgumentError("particle array $(name) must use one-based axes 1:$N, got $(axes(a, 1))"),
    )
    return nothing
end

function ParticleSet{D,T}(N::Integer; q = one(T), m = one(T)) where {D,T}
    D >= 1 || throw(ArgumentError("D must be ≥ 1"))
    Np = _particle_length(N)
    x = ntuple(_ -> Vector{T}(undef, Np), Val(D))
    v = ntuple(_ -> zeros(T, Np), Val(3))
    return ParticleSet{D,T,Vector{T},Vector{T},Vector{T},Vector{UInt64},Vector{UInt32}}(
        x,
        v,
        ones(T, Np),
        _particle_ids(Np),
        zeros(UInt32, Np),
        T(q),
        T(m),
    )
end

function ParticleSet{D,T}(
    x::NTuple{D,X},
    v::NTuple{3,V},
    weight::W,
    id::I,
    tag::G,
    q,
    m,
) where {
    D,
    T,
    X<:AbstractVector{T},
    V<:AbstractVector{T},
    W<:AbstractVector{T},
    I<:AbstractVector{UInt64},
    G<:AbstractVector{UInt32},
}
    D >= 1 || throw(ArgumentError("D must be ≥ 1"))
    N = length(weight)
    all(length(x[d]) == N for d = 1:D) ||
        throw(DimensionMismatch("all position arrays must have length $N"))
    all(length(v[c]) == N for c = 1:3) ||
        throw(DimensionMismatch("all velocity arrays must have length $N"))
    length(id) == N || throw(DimensionMismatch("id length $(length(id)) must equal $N"))
    length(tag) == N || throw(DimensionMismatch("tag length $(length(tag)) must equal $N"))
    for d = 1:D
        _check_particle_vector_axes(Symbol(:x, d), x[d], N)
    end
    for c = 1:3
        _check_particle_vector_axes(Symbol(:v, c), v[c], N)
    end
    _check_particle_vector_axes(:weight, weight, N)
    _check_particle_vector_axes(:id, id, N)
    _check_particle_vector_axes(:tag, tag, N)
    return ParticleSet{D,T,X,V,W,I,G}(x, v, weight, id, tag, T(q), T(m))
end

nparticles(ps::ParticleSet) = length(ps.weight)
Base.eltype(::ParticleSet{D,T}) where {D,T} = T

# ---------------------------------------------------------------- loading
