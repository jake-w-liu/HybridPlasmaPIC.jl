module HybridPlasmaPICCUDAExt

import CUDA
import HybridPlasmaPIC

HybridPlasmaPIC.extension_dependency_module(::Val{:cuda}) = CUDA
HybridPlasmaPIC.extension_device_array_type(::Val{:cuda}) = CUDA.CuArray

_cuda_functional() =
    try
        CUDA.functional()
    catch err
        false
    end

function _cuda_optional_bytes(f)
    try
        return f()
    catch err
        return nothing
    end
end

function HybridPlasmaPIC.backend_memory_status(::Val{:cuda})
    functional = _cuda_functional()
    if !functional
        return HybridPlasmaPIC.BackendMemoryStatus(
            :cuda,
            false,
            isdefined(CUDA, :pool_status) || isdefined(CUDA, :reclaim);
            note = "CUDA package is loaded, but no functional CUDA device is available",
        )
    end

    free = _cuda_optional_bytes(CUDA.free_memory)
    total = _cuda_optional_bytes(CUDA.total_memory)
    used = isdefined(CUDA, :used_memory) ? _cuda_optional_bytes(CUDA.used_memory) : nothing
    cached = isdefined(CUDA, :cached_memory) ? _cuda_optional_bytes(CUDA.cached_memory) : nothing
    reserved = used === nothing || cached === nothing ? nothing : used + cached
    return HybridPlasmaPIC.BackendMemoryStatus(
        :cuda,
        true,
        isdefined(CUDA, :pool_status) || isdefined(CUDA, :reclaim);
        total_bytes = total,
        free_bytes = free,
        used_bytes = used,
        cached_bytes = cached,
        reserved_bytes = reserved,
        note = "CUDA memory telemetry",
    )
end

function HybridPlasmaPIC.reclaim_backend_memory!(::Val{:cuda})
    functional = _cuda_functional()
    functional || return false
    if isdefined(CUDA, :reclaim)
        CUDA.reclaim()
        return true
    end
    return false
end

function HybridPlasmaPIC.disallow_scalar_indexing!(::Val{:cuda})
    CUDA.allowscalar(false)
    return nothing
end

end # module
