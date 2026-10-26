using HybridPlasmaPIC, PencilArrays, PencilFFTs, MPI, Test

Base.@nospecializeinfer function _reject_metadata_alias(
    @nospecialize(operation),
    @nospecialize(arguments),
    origin,
    g,
    fields,
)
    metadata = (g.kfull..., g.kvec..., g.ik..., g.midx...)
    metadata_before = map(copy, metadata)
    fields_before = map(field -> copy(parent(field)), fields)
    workspaces = (g.physical, g.spectral, g.spectral_work..., g.staging, g.staging_work...)
    for (index, workspace) in enumerate(workspaces)
        fill!(parent(workspace), eltype(workspace)(19 + index))
    end
    workspace_before = map(workspace -> copy(parent(workspace)), workspaces)
    rank = MPI.Comm_rank(g.fft.comm)
    caught = try
        operation(arguments...)
        nothing
    catch error
        error
    end
    @test caught isa (rank == origin ? ArgumentError : ErrorException)
    @test rank != origin || (
        caught isa ArgumentError &&
        occursin("must not alias private Fourier workspace storage", sprint(showerror, caught))
    )
    @test all(i -> isequal(metadata[i], metadata_before[i]), eachindex(metadata))
    @test all(i -> isequal(parent(fields[i]), fields_before[i]), eachindex(fields))
    @test all(i -> isequal(parent(workspaces[i]), workspace_before[i]), eachindex(workspaces))
    @test MPI.Allreduce(1, +, g.fft.comm) == MPI.Comm_size(g.fft.comm)
    return nothing
end

# Both spans share initialized aligned storage, including Int32/Float64 where
# the physical span extends beyond the shorter logical metadata span.
function _with_metadata_span(test, values::Vector{S}, ::Type{T}, shape) where {S,T}
    field_bytes = Base.checked_mul(sizeof(T), prod(shape))
    metadata_bytes = Base.checked_mul(sizeof(S), length(values))
    alignment = max(Base.datatype_alignment(T), Base.datatype_alignment(S))
    capacity = Base.checked_add(max(field_bytes, metadata_bytes), alignment - 1)
    backing = zeros(UInt8, capacity)
    GC.@preserve backing values begin
        remainder = Int(UInt(pointer(backing)) % UInt(alignment))
        offset = mod(alignment - remainder, alignment)
        @test offset + field_bytes <= length(backing)
        @test offset + metadata_bytes <= length(backing)
        start = pointer(backing, offset + 1)
        @test iszero(UInt(start) % UInt(alignment))
        metadata = unsafe_wrap(Vector{S}, Ptr{S}(start), length(values); own = false)
        copyto!(metadata, values)
        @test isequal(metadata, values)
        storage = unsafe_wrap(Array, Ptr{T}(start), shape; own = false)
        GC.@preserve metadata test(metadata, storage)
    end
    return nothing
end

function _with_metadata_alias_fixture(test, g, f, family, axis)
    values = getproperty(g, family)[axis]
    T = eltype(f)
    shape = size(parent(f))
    if eltype(values) === T
        GC.@preserve values test(g, values, reshape(values, shape))
        return nothing
    end
    _with_metadata_span(values, T, shape) do metadata, storage
        replacement = Base.setindex(getproperty(g, family), metadata, axis)
        arguments = map(fieldnames(typeof(g))) do name
            name === family ? replacement : getfield(g, name)
        end
        grid = typeof(g)(arguments...)
        @test getproperty(grid, family)[axis] === metadata
        test(grid, metadata, storage)
    end
    return nothing
end

@testset "Metadata spans cover both integer widths without an architecture skip" begin
    for S in (Int32, Int64, ComplexF32, ComplexF64), T in (Float32, Float64)
        original = S.(1:8)
        _with_metadata_span(original, T, (8, 1, 1)) do metadata, storage
            @test size(storage) == (8, 1, 1)
            @test eltype(storage) === T
            @test Base.mightalias(storage, metadata)
            @test Base.dataids(storage) == Base.dataids(metadata)
            bytes = reinterpret(UInt8, vec(storage))
            extra_start = min(sizeof(S) * length(metadata), length(bytes)) + 1
            @test all(iszero, view(bytes, extra_start:length(bytes)))
            @test isequal(metadata, original)
        end
    end
end

@testset "Distributed operands reject every metadata family collectively" begin
    ensure_mpi_initialized!()
    comm = MPI.COMM_WORLD
    rank, ranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    dims = MPI.Dims_create(ranks, (Cint(0), Cint(0)))
    p, q = Int(dims[1]), Int(dims[2])
    for T in (Float32, Float64), family in (:kfull, :kvec, :ik, :midx)
        configurations = [(1, (8, p, q)), (2, (p, 8p, q))]
        p == q && push!(configurations, (3, (p, p, 8q)))
        for (axis, n) in configurations
            g = distributed_fourier_grid(n, (T(2), T(3), T(5)); comm)
            f = distributed_field(g)
            view = PencilArrays.global_view(f)
            for I in CartesianIndices(view)
                view[I] = sin(T(2pi) * T(I[axis] - 1) / T(n[axis]))
            end
            A = (copy(f), copy(f), copy(f))
            O = (distributed_field(g), distributed_field(g), distributed_field(g))
            foreach(field -> fill!(parent(field), T(7)), O)
            fields = (f, A..., O...)
            reference = copy(parent(deriv(f, g, axis)))
            _with_metadata_alias_fixture(g, f, family, axis) do grid, values, storage
                # Replacing only a numerically identical metadata vector keeps
                # valid geometry and transform behavior before alias rejection.
                @test parent(deriv!(O[1], f, grid, axis)) == reference
                alias = PencilArrays.PencilArray(PencilArrays.pencil(f), storage)
                @test typeof(alias) === typeof(f)
                @test axes(alias) == axes(f)
                @test PencilArrays.pencil(alias) === PencilArrays.pencil(f)
                @test Base.mightalias(alias, values)
                @test Base.dataids(alias) == Base.dataids(values)
                for origin = 0:ranks-1
                    selected(real, invalid) = rank == origin ? invalid : real
                    reject(operation, valid, invalid) = _reject_metadata_alias(
                        operation,
                        selected(valid, invalid),
                        origin,
                        grid,
                        fields,
                    )
                    # Every mutating output surface rejects the shared geometry.
                    reject(deriv!, (O[1], f, grid, axis), (alias, f, grid, axis))
                    reject(gradient!, (O, f, grid), ((O[1], alias, O[3]), f, grid))
                    reject(divergence!, (O[1], A, grid), (alias, A, grid))
                    reject(curl!, (O, A, grid), ((O[1], O[2], alias), A, grid))
                    reject(laplacian!, (O[1], f, grid), (alias, f, grid))
                    reject(project_divfree!, (A, grid), ((A[1], alias, A[3]), grid))
                    reject(exp_filter!, (f, grid), (alias, grid))
                    reject(dealias_two_thirds!, (f, grid), (alias, grid))
                    # Metadata is also inadmissible as a physical input.
                    reject(deriv!, (O[1], f, grid, axis), (O[1], alias, grid, axis))
                    reject(deriv, (f, grid, axis), (alias, grid, axis))
                    reject(gradient!, (O, f, grid), (O, alias, grid))
                    reject(divergence!, (O[1], A, grid), (O[1], (alias, A[2], A[3]), grid))
                    reject(curl!, (O, A, grid), (O, (A[1], alias, A[3]), grid))
                    reject(laplacian!, (O[1], f, grid), (O[1], alias, grid))
                    # The same workspace remains usable after every rejection path.
                    @test parent(deriv!(O[1], f, grid, axis)) == reference
                end
            end
        end
    end
end
