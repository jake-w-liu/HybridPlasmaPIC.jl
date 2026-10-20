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
    _check_spatial_dimension(D)
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
    _check_spatial_dimension(D)
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

function ParticleSet{D,T}(
    x::Tuple,
    v::Tuple,
    weight::AbstractVector,
    id::AbstractVector,
    tag::AbstractVector,
    q,
    m,
) where {D,T}
    _check_spatial_dimension(D)
    length(x) == D || throw(DimensionMismatch("position tuple length $(length(x)) must equal D=$D"))
    length(v) == 3 || throw(DimensionMismatch("velocity tuple length $(length(v)) must equal 3"))
    for d = 1:D
        x[d] isa AbstractVector{T} || throw(ArgumentError("x[$d] must be an AbstractVector{$T}"))
    end
    for c = 1:3
        v[c] isa AbstractVector{T} || throw(ArgumentError("v[$c] must be an AbstractVector{$T}"))
    end
    weight isa AbstractVector{T} || throw(ArgumentError("weight must be an AbstractVector{$T}"))
    id isa AbstractVector{UInt64} || throw(ArgumentError("id must be an AbstractVector{UInt64}"))
    tag isa AbstractVector{UInt32} || throw(ArgumentError("tag must be an AbstractVector{UInt32}"))
    throw(ArgumentError("ParticleSet custom arrays must have compatible vector storage"))
end

nparticles(ps::ParticleSet) = length(ps.weight)
Base.eltype(::ParticleSet{D,T}) where {D,T} = T

# ---------------------------------------------------------------- loading

# Shared storage checks also serve collision scratch, collections and checkpoint restore.
function _particle_backing_mightalias(a::AbstractArray, b::AbstractArray)
    (a === b || Base.mightalias(a, b)) && return true
    pa, pb = parent(a), parent(b)
    pa === a || (_particle_backing_mightalias(pa, b) && return true)
    pb === b || (_particle_backing_mightalias(a, pb) && return true)
    # A view's mutable indices can grow along with another particle array even
    # when its data parent is independent. Ranges contain no resizable storage.
    if a isa SubArray
        for index in parentindices(a)
            index isa AbstractArray &&
                !(index isa AbstractRange) &&
                _particle_backing_mightalias(index, b) &&
                return true
        end
    end
    if b isa SubArray
        for index in parentindices(b)
            index isa AbstractArray &&
                !(index isa AbstractRange) &&
                _particle_backing_mightalias(a, index) &&
                return true
        end
    end
    return false
end

function _particle_arrays_mightalias(a::AbstractArray, b::AbstractArray)
    (a === b || Base.mightalias(a, b)) && return true
    # Empty arrays have no overlapping elements, but wrappers can still share
    # storage whose later resize would change another particle component.
    return (isempty(a) || isempty(b)) && _particle_backing_mightalias(a, b)
end

function _particle_components_mightalias(ps::ParticleSet, alias_test = _particle_arrays_mightalias)
    for (d, x) in enumerate(ps.x)
        for k = 1:(d-1)
            alias_test(x, ps.x[k]) && return true
        end
        for v in ps.v
            alias_test(x, v) && return true
        end
        alias_test(x, ps.weight) && return true
        alias_test(x, ps.id) && return true
        alias_test(x, ps.tag) && return true
    end
    for (c, v) in enumerate(ps.v)
        for k = 1:(c-1)
            alias_test(v, ps.v[k]) && return true
        end
        alias_test(v, ps.weight) && return true
        alias_test(v, ps.id) && return true
        alias_test(v, ps.tag) && return true
    end
    alias_test(ps.weight, ps.id) && return true
    alias_test(ps.weight, ps.tag) && return true
    alias_test(ps.id, ps.tag) && return true
    return false
end

function _mightalias_particle_storage(
    candidate,
    ps::ParticleSet,
    alias_test = _particle_arrays_mightalias,
)
    for d = 1:length(ps.x)
        alias_test(candidate, ps.x[d]) && return true
    end
    for c = 1:3
        alias_test(candidate, ps.v[c]) && return true
    end
    alias_test(candidate, ps.weight) && return true
    alias_test(candidate, ps.id) && return true
    alias_test(candidate, ps.tag) && return true
    return false
end

function _particle_sets_mightalias(
    a::ParticleSet,
    b::ParticleSet,
    alias_test = _particle_arrays_mightalias,
)
    for d = 1:length(a.x)
        _mightalias_particle_storage(a.x[d], b, alias_test) && return true
    end
    for c = 1:3
        _mightalias_particle_storage(a.v[c], b, alias_test) && return true
    end
    _mightalias_particle_storage(a.weight, b, alias_test) && return true
    _mightalias_particle_storage(a.id, b, alias_test) && return true
    _mightalias_particle_storage(a.tag, b, alias_test) && return true
    return false
end

function _require_particle_resize_support(a::AbstractVector, N::Int, context::AbstractString)
    if length(a) != N && !applicable(resize!, a, N)
        throw(
            ArgumentError(
                "$context requires $N particles, but destination $(typeof(a)) cannot resize from $(length(a))",
            ),
        )
    end
    return nothing
end

function _require_particle_resize_support(ps::ParticleSet, N::Int, context::AbstractString)
    # Check homogeneous position/velocity tuples separately from id/tag storage.
    # This keeps the applicable(resize!, ...) calls concretely typed.
    for x in ps.x
        _require_particle_resize_support(x, N, context)
    end
    for v in ps.v
        _require_particle_resize_support(v, N, context)
    end
    _require_particle_resize_support(ps.weight, N, context)
    _require_particle_resize_support(ps.id, N, context)
    _require_particle_resize_support(ps.tag, N, context)
    return nothing
end

function _check_particle_storage_axes(ps::ParticleSet)
    N = nparticles(ps)
    for (d, x) in enumerate(ps.x)
        axes(x, 1) == Base.OneTo(N) || _check_particle_vector_axes(Symbol(:x, d), x, N)
    end
    for (c, v) in enumerate(ps.v)
        axes(v, 1) == Base.OneTo(N) || _check_particle_vector_axes(Symbol(:v, c), v, N)
    end
    _check_particle_vector_axes(:weight, ps.weight, N)
    _check_particle_vector_axes(:id, ps.id, N)
    _check_particle_vector_axes(:tag, ps.tag, N)
    return N
end

function _validate_particle_growth(ps::ParticleSet, count::Int, context::AbstractString)
    N = _check_particle_storage_axes(ps)
    0 <= count <= typemax(Int) - N ||
        throw(ArgumentError("$context particle count exceeds Int capacity"))
    # Growth affects backing storage and mutable view indices, even where element
    # ranges do not overlap. Ordinary fixed-size writes keep the less strict predicate.
    _particle_components_mightalias(ps, _particle_backing_mightalias) &&
        throw(ArgumentError("$context particle components must not share backing storage"))
    _require_particle_resize_support(ps, N + count, context)
    return nothing
end
