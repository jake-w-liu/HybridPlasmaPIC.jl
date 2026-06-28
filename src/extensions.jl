# extensions.jl — explicit optional package-extension surface.

const _SUPPORTED_EXTENSIONS = (:cuda, :metal, :io, :pencilfft)

"Supported optional extension keys."
supported_extensions() = _SUPPORTED_EXTENSIONS

extension_name(::Val{:cuda}) = :HybridPlasmaPICCUDAExt
extension_name(::Val{:metal}) = :HybridPlasmaPICMetalExt
extension_name(::Val{:io}) = :HybridPlasmaPICIOExt
extension_name(::Val{:pencilfft}) = :HybridPlasmaPICPencilFFTSExt

function extension_name(::Val{name}) where {name}
    throw(
        ArgumentError(
            "unsupported HybridPlasmaPIC extension $(name); expected one of $(_SUPPORTED_EXTENSIONS)",
        ),
    )
end

extension_dependency_name(::Val{:cuda}) = :CUDA
extension_dependency_name(::Val{:metal}) = :Metal
extension_dependency_name(::Val{:io}) = :HDF5
extension_dependency_name(::Val{:pencilfft}) = :PencilFFTs

function extension_dependency_name(::Val{name}) where {name}
    throw(
        ArgumentError(
            "unsupported HybridPlasmaPIC extension $(name); expected one of $(_SUPPORTED_EXTENSIONS)",
        ),
    )
end

"""
    extension_loaded(Val(name)) -> Bool

Return whether Julia has loaded HybridPlasmaPIC's package extension for `name`.
This uses Julia's package-extension loader and does not import optional
dependencies by itself.
"""
function extension_loaded(::Val{name}) where {name}
    return Base.get_extension(@__MODULE__, extension_name(Val(name))) !== nothing
end

"Tuple of currently loaded HybridPlasmaPIC extension keys."
loaded_extensions() = Tuple(name for name in _SUPPORTED_EXTENSIONS if extension_loaded(Val(name)))

"""
    require_extension(Val(name)) -> Module

Return the loaded extension module for `name`, or throw a clear error naming the
dependency that must be loaded first.
"""
function require_extension(::Val{name}) where {name}
    extname = extension_name(Val(name))
    ext = Base.get_extension(@__MODULE__, extname)
    if ext === nothing
        dep = extension_dependency_name(Val(name))
        error(
            "HybridPlasmaPIC extension $(extname) is not loaded; load/import $(dep) with HybridPlasmaPIC first",
        )
    end
    return ext
end

function _extension_missing(name::Symbol, feature::Symbol)
    require_extension(Val(name))
    error(
        "HybridPlasmaPIC extension $(extension_name(Val(name))) loaded but did not provide $(feature)",
    )
end

extension_dependency_module(::Val{name}) where {name} =
    _extension_missing(name, :extension_dependency_module)
extension_device_array_type(::Val{name}) where {name} =
    _extension_missing(name, :extension_device_array_type)
disallow_scalar_indexing!(::Val{name}) where {name} =
    _extension_missing(name, :disallow_scalar_indexing!)

write_field_hdf5(args...; kwargs...) = _extension_missing(:io, :write_field_hdf5)
read_field_hdf5(args...; kwargs...) = _extension_missing(:io, :read_field_hdf5)

distributed_fft_plan(args...; kwargs...) = _extension_missing(:pencilfft, :distributed_fft_plan)
distributed_fft_input(args...; kwargs...) = _extension_missing(:pencilfft, :distributed_fft_input)
distributed_fft_output(args...; kwargs...) = _extension_missing(:pencilfft, :distributed_fft_output)
distributed_fft_forward!(args...; kwargs...) =
    _extension_missing(:pencilfft, :distributed_fft_forward!)
distributed_fft_inverse!(args...; kwargs...) =
    _extension_missing(:pencilfft, :distributed_fft_inverse!)
distributed_fft_roundtrip_error(args...; kwargs...) =
    _extension_missing(:pencilfft, :distributed_fft_roundtrip_error)
