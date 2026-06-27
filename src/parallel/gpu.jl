# gpu.jl (§5.3 parallel/) — backend-resident particle storage helpers.
#
# The production CUDA/Metal particle kernels are added by package extensions.
# This file owns the backend-neutral contract: particle arrays can be copied to
# a loaded backend's array type and copied back to host memory without changing
# particle semantics. CPU routines remain the reference implementation.

const _PARTICLE_ARRAY_FIELDS = (:x, :v, :weight, :id, :tag)

_particle_arrays(ps::ParticleSet) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)

"""
    BackendMemoryStatus

Memory telemetry for a compute backend. `total_bytes`, `free_bytes`,
`used_bytes`, `cached_bytes`, and `reserved_bytes` are `nothing` when the loaded
backend does not expose that counter. `pool_supported` reports whether the
backend exposes explicit pool/cache telemetry or reclamation hooks.
"""
Base.@kwdef struct BackendMemoryStatus
    backend::Symbol
    device_available::Bool
    pool_supported::Bool
    total_bytes::Union{Nothing,Int} = nothing
    free_bytes::Union{Nothing,Int} = nothing
    used_bytes::Union{Nothing,Int} = nothing
    cached_bytes::Union{Nothing,Int} = nothing
    reserved_bytes::Union{Nothing,Int} = nothing
    note::String = ""
end

function _nonnegative_bytes(value, field::Symbol)
    value === nothing && return nothing
    value isa Integer || throw(ArgumentError("$(field) must be an integer byte count or nothing"))
    value >= 0 || throw(ArgumentError("$(field) must be nonnegative"))
    value <= typemax(Int) || throw(ArgumentError("$(field) exceeds Int range"))
    return Int(value)
end

function BackendMemoryStatus(
    backend::Symbol,
    device_available::Bool,
    pool_supported::Bool;
    total_bytes = nothing,
    free_bytes = nothing,
    used_bytes = nothing,
    cached_bytes = nothing,
    reserved_bytes = nothing,
    note::AbstractString = "",
)
    return BackendMemoryStatus(
        backend,
        device_available,
        pool_supported,
        _nonnegative_bytes(total_bytes, :total_bytes),
        _nonnegative_bytes(free_bytes, :free_bytes),
        _nonnegative_bytes(used_bytes, :used_bytes),
        _nonnegative_bytes(cached_bytes, :cached_bytes),
        _nonnegative_bytes(reserved_bytes, :reserved_bytes),
        String(note),
    )
end

function memory_pressure(status::BackendMemoryStatus)
    total = status.total_bytes
    total === nothing && return nothing
    total == 0 && return nothing
    used = status.used_bytes
    if used === nothing && status.free_bytes !== nothing
        used = max(total - status.free_bytes, 0)
    end
    used === nothing && return nothing
    return clamp(used / total, 0.0, 1.0)
end

function backend_memory_status(::Val{:cpu})
    total = _nonnegative_bytes(Sys.total_memory(), :total_bytes)
    free = _nonnegative_bytes(Sys.free_memory(), :free_bytes)
    used = max(total - free, 0)
    return BackendMemoryStatus(
        :cpu,
        true,
        false;
        total_bytes = total,
        free_bytes = free,
        used_bytes = used,
        note = "CPU backend has no GPU memory pool",
    )
end

function backend_memory_status(::Val{name}) where {name}
    require_extension(Val(name))
    error(
        "HybridPlasmaPIC extension $(extension_name(Val(name))) loaded but did not provide backend_memory_status",
    )
end

function backend_memory_status(::Val{:mixed})
    return BackendMemoryStatus(
        :mixed,
        true,
        false;
        note = "particle arrays are split across multiple backends",
    )
end

function backend_memory_status(::Val{:unknown})
    return BackendMemoryStatus(
        :unknown,
        false,
        false;
        note = "particle array backend could not be identified",
    )
end

backend_memory_status(ps::ParticleSet) = backend_memory_status(Val(particle_storage_backend(ps)))

reclaim_backend_memory!(::Val{:cpu}) = false

function reclaim_backend_memory!(::Val{name}) where {name}
    require_extension(Val(name))
    error(
        "HybridPlasmaPIC extension $(extension_name(Val(name))) loaded but did not provide reclaim_backend_memory!",
    )
end

particle_array_backend(A::Array) = :cpu

function particle_array_backend(A::AbstractArray)
    for name in (:cuda, :metal)
        if extension_loaded(Val(name))
            AT = extension_device_array_type(Val(name))
            A isa AT && return name
        end
    end
    return :unknown
end

"""
    particle_storage_backend(ps) -> Symbol

Return `:cpu`, `:cuda`, `:metal`, `:mixed`, or `:unknown` for the storage backing
all particle arrays in `ps`. `:mixed` means arrays are backed by more than one
known backend; `:unknown` means at least one non-`Array` backend cannot be
identified from the loaded extensions.
"""
function particle_storage_backend(ps::ParticleSet)
    backends = map(particle_array_backend, _particle_arrays(ps))
    any(==(:unknown), backends) && return :unknown
    first_backend = first(backends)
    all(==(first_backend), backends) && return first_backend
    return :mixed
end

function _copy_particle_array_to_backend(::Val{:cpu}, A::AbstractArray)
    return Array(A)
end

function _copy_particle_array_to_backend(::Val{name}, A::AbstractArray) where {name}
    AT = extension_device_array_type(Val(name))
    return AT(A)
end

validate_particle_backend_eltype(::Val{name}, ::Type{T}) where {name,T} = nothing

"""
    copy_particles_to_backend(Val(backend), ps; disallow_scalar_indexing=true)

Return a new [`ParticleSet`](@ref) whose arrays live on `backend`. `backend` may
be `:cpu` or a loaded optional GPU extension key such as `:cuda` or `:metal`.
GPU backends require importing the matching weak dependency first.
"""
function copy_particles_to_backend(
    ::Val{name},
    ps::ParticleSet{D,T};
    disallow_scalar_indexing::Bool = true,
) where {name,D,T}
    if name !== :cpu
        prepare_gpu_backend!(Val(name); disallow_scalar_indexing)
        validate_particle_backend_eltype(Val(name), T)
    end
    x = ntuple(d -> _copy_particle_array_to_backend(Val(name), ps.x[d]), D)
    v = ntuple(c -> _copy_particle_array_to_backend(Val(name), ps.v[c]), 3)
    weight = _copy_particle_array_to_backend(Val(name), ps.weight)
    id = _copy_particle_array_to_backend(Val(name), ps.id)
    tag = _copy_particle_array_to_backend(Val(name), ps.tag)
    return ParticleSet{D,T}(x, v, weight, id, tag, ps.q, ps.m)
end

"""
    copy_particles_to_host(ps) -> ParticleSet

Return a CPU-backed copy of `ps`. This is the required staging operation before
calling CPU-only diagnostics or serialization on backend-resident particles.
"""
copy_particles_to_host(ps::ParticleSet) = copy_particles_to_backend(Val(:cpu), ps)

"""
    prepare_gpu_backend!(Val(backend); disallow_scalar_indexing=true)

Require and initialize a loaded GPU extension. When supported by the backend
extension, scalar indexing is disabled so accidental CPU fallbacks fail fast.
"""
function prepare_gpu_backend!(::Val{name}; disallow_scalar_indexing::Bool = true) where {name}
    name === :cpu && return nothing
    ext = require_extension(Val(name))
    disallow_scalar_indexing && disallow_scalar_indexing!(Val(name))
    return ext
end
