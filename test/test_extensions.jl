using HybridPlasmaPIC
using Test

const HDF5_AVAILABLE = try
    @eval import HDF5
    true
catch err
    @warn "HDF5 unavailable; skipping I/O extension tests" exception = err
    false
end

@testset "explicit package extensions" begin
    @test supported_extensions() == (:cuda, :metal, :io, :pencilfft)
    @test extension_name(Val(:cuda)) == :HybridPlasmaPICCUDAExt
    @test extension_name(Val(:metal)) == :HybridPlasmaPICMetalExt
    @test extension_name(Val(:io)) == :HybridPlasmaPICIOExt
    @test extension_name(Val(:pencilfft)) == :HybridPlasmaPICPencilFFTSExt
    @test extension_dependency_name(Val(:cuda)) == :CUDA
    @test extension_dependency_name(Val(:metal)) == :Metal
    @test extension_dependency_name(Val(:io)) == :HDF5
    @test extension_dependency_name(Val(:pencilfft)) == :PencilFFTs
    @test_throws ArgumentError extension_name(Val(:rocm))
    @test_throws ArgumentError extension_dependency_name(Val(:rocm))

    @test !extension_loaded(Val(:cuda))
    @test !extension_loaded(Val(:metal))

    @test !(:mpi in supported_extensions())

    if HDF5_AVAILABLE
        import HDF5
        @test extension_loaded(Val(:io))
        @test extension_dependency_module(Val(:io)) === HDF5

        mktempdir() do dir
            path = joinpath(dir, "field.h5")
            A = reshape(collect(Float64, 1:12), 3, 4)
            @test write_field_hdf5(path, "density", A) == path
            @test read_field_hdf5(path, "density") == A
            B = reshape(collect(Float64, 13:24), 3, 4)
            @test write_field_hdf5(path, "pressure", B) == path
            @test read_field_hdf5(path, "density") == A
            @test read_field_hdf5(path, "pressure") == B
        end

        @testset "HDF5 field buffers preserve logical array order" begin
            base = reshape(collect(1.0:30.0), 6, 5)
            io_ext = Base.get_extension(HybridPlasmaPIC, :HybridPlasmaPICIOExt)
            @test io_ext._field_write_buffer(base) === base
            contiguous_view = @view base[:, 2:4]
            @test io_ext._field_write_buffer(contiguous_view) === contiguous_view
            inputs = (
                reshape(1.0:12.0, 3, 4),
                @view(base[2:4, 2:4]),
                @view(base[1:2:6, 2:4]),
                @view(base[:, 2:4]),
                transpose(base),
                adjoint(complex.(base, -base)),
                zeros(Float64, 0, 3),
                reshape([42.0]),
                @view(base[:, 1:0]),
            )
            mktempdir() do dir
                path = joinpath(dir, "layouts.h5")
                for (i, A) in enumerate(inputs)
                    dataset = "fields/value_$i"
                    @test write_field_hdf5(path, dataset, A) == path
                    B = read_field_hdf5(path, dataset)
                    @test size(B) == size(A)
                    @test eltype(B) === eltype(A)
                    # HDF5 scalar dataspaces read as scalars, including 0D inputs.
                    @test B == (ndims(A) == 0 ? A[] : A)
                end
                @test_throws Exception read_field_hdf5(path, "missing")
                @test read_field_hdf5(path, "fields/value_2") == inputs[2]
            end
        end
    else
        @test_skip "HDF5 not available in this environment"
    end
end
