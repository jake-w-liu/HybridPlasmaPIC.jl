module HybridPlasmaPICIOExt

import HDF5
import HybridPlasmaPIC

HybridPlasmaPIC.extension_dependency_module(::Val{:io}) = HDF5

function _field_write_buffer(A::AbstractArray)
    # HDF5 writes a flat pointer and checks only the first stride; views can
    # have gaps between columns even when that first stride is one.
    if A isa StridedArray && stride(A, 1) == 1
        expected_stride = 1
        contiguous = true
        for d = 1:ndims(A)
            if size(A, d) > 1 && stride(A, d) != expected_stride
                contiguous = false
                break
            end
            expected_stride *= size(A, d)
        end
        contiguous && return A
    end
    return Array(A)
end

function HybridPlasmaPIC.write_field_hdf5(
    path::AbstractString,
    dataset::AbstractString,
    A::AbstractArray,
)
    data = _field_write_buffer(A)
    HDF5.h5open(path, "cw") do fid
        fid[dataset] = data
    end
    return path
end

function HybridPlasmaPIC.read_field_hdf5(path::AbstractString, dataset::AbstractString)
    HDF5.h5open(path, "r") do fid
        return read(fid[dataset])
    end
end

end # module
