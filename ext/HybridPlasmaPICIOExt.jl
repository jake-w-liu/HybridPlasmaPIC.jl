module HybridPlasmaPICIOExt

import HDF5
import HybridPlasmaPIC

HybridPlasmaPIC.extension_dependency_module(::Val{:io}) = HDF5

function HybridPlasmaPIC.write_field_hdf5(
    path::AbstractString,
    dataset::AbstractString,
    A::AbstractArray,
)
    HDF5.h5open(path, "w") do fid
        fid[dataset] = A
    end
    return path
end

function HybridPlasmaPIC.read_field_hdf5(path::AbstractString, dataset::AbstractString)
    HDF5.h5open(path, "r") do fid
        return read(fid[dataset])
    end
end

end # module
