using HybridPlasmaPIC, MPI, Test
import PencilArrays
import PencilFFTs

function _distributed_derivative_failure(f, origin, comm)
    rank = MPI.Comm_rank(comm)
    caught = try
        f()
        nothing
    catch err
        err
    end
    @test caught isa (rank == origin ? ArgumentError : ErrorException)
    @test MPI.Allreduce(1, +, comm) == MPI.Comm_size(comm)
end

# Independent finite DFT: no FFT library or shared multiplier/metadata helper.
function _dense_derivative_reference(values, lengths, axis)
    n = size(values)
    points = collect(CartesianIndices(values))
    coefficients = ComplexF64[
        sum(
            ComplexF64(values[X]) * cis(-2pi * sum((K[d] - 1) * (X[d] - 1) / n[d] for d = 1:3))
            for X in points
        ) for K in points
    ]
    result = zeros(Float64, n)
    for X in points
        value = 0.0im
        for (i, K) in enumerate(points)
            mode = K[axis] - 1
            mode = mode <= n[axis] ÷ 2 ? mode : mode - n[axis]
            if iseven(n[axis]) && K[axis] == n[axis] ÷ 2 + 1
                mode = 0
            end
            wave_number = 2pi * mode / Float64(lengths[axis])
            value +=
                coefficients[i] *
                im *
                wave_number *
                cis(2pi * sum((K[d] - 1) * (X[d] - 1) / n[d] for d = 1:3))
        end
        result[X] = real(value) / length(values)
    end
    return result
end

@testset "Rank-owned periodic derivatives" begin
    ensure_mpi_initialized!()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    ranks = MPI.Comm_size(comm)
    for T in (Float32, Float64), n in ((8, 6, 4), (5, 7, 9), (5, 4, 3))
        L = (T(2), T(3), T(5))
        g = distributed_fourier_grid(n, L; comm)
        field = distributed_field(g)
        out = distributed_field(g)
        @test eltype(field) === T
        @test all(iszero, parent(field))
        @test PencilArrays.pencil(field) === PencilArrays.pencil(g.physical)
        @test g.n == n && g.L == L
        @test MPI.Allreduce(length(field), +, comm) == prod(n)
        global_field = PencilArrays.global_view(field)
        global_out = PencilArrays.global_view(out)
        for axis = 1:3
            for I in CartesianIndices(global_field)
                global_field[I] = sum(T(d) * sin(T(2pi) * T(I[d] - 1) / T(n[d])) for d = 1:3)
            end
            original = copy(parent(field))
            @test deriv!(out, field, g, axis) === out
            bound = T(4096) * eps(T) * T(axis) * T(2pi) / L[axis]
            for I in CartesianIndices(global_out)
                expected = T(axis) * T(2pi) / L[axis] * cos(T(2pi) * T(I[axis] - 1) / T(n[axis]))
                @test isapprox(global_out[I], expected; atol = bound, rtol = T(4096) * eps(T))
            end
            @test isequal(parent(field), original)
            @test parent(deriv(field, g, axis)) == parent(out)
            @test deriv!(field, field, g, axis) === field
            @test parent(field) == parent(out)
            if n == (5, 4, 3)
                for I in CartesianIndices(global_field)
                    global_field[I] = sin(T(0.17) * T(I[1])) + cos(T(0.31) * T(I[2] * I[3]))
                end
                deriv!(out, field, g, axis)
                full_field = PencilArrays.gather(field)
                full_out = PencilArrays.gather(out)
                if rank == 0
                    dense = _dense_derivative_reference(full_field, L, axis)
                    @test isapprox(full_out, dense; atol = 4096eps(T), rtol = 4096eps(T))
                    serial = deriv(full_field, FourierGrid(n, L), axis)
                    @test isapprox(full_out, serial; atol = 4096eps(T), rtol = 4096eps(T))
                end
            end
        end
        for axis = 1:3
            fill!(parent(field), T(3))
            deriv!(out, field, g, axis)
            @test maximum(abs, parent(out)) <= T(4096) * eps(T)
            if iseven(n[axis])
                for I in CartesianIndices(global_field)
                    global_field[I] = iseven(I[axis] - 1) ? one(T) : -one(T)
                end
                deriv!(out, field, g, axis)
                @test maximum(abs, parent(out)) <= T(4096) * eps(T)
            end
        end
    end
    for T in (Float32, Float64)
        n = (8, 8, 8)
        L = (T(2), T(3), T(5))
        g = distributed_fourier_grid(n, L; comm)
        field = distributed_field(g)
        out = distributed_field(g)
        fill!(parent(field), floatmax(T) / T(2))
        for axis = 1:3
            deriv!(out, field, g, axis)
            @test all(iszero, parent(out))
        end
        for origin = 0:ranks-1, invalid in (:axis, :type, :nan, :workspace, :pencil)
            fill!(parent(field), one(T))
            fill!(parent(out), T(7))
            original = copy(parent(out))
            other = distributed_fourier_grid(n, L; comm)
            bad_field =
                invalid == :pencil ? distributed_field(other) :
                invalid == :type ? similar(field, Int) : field
            invalid == :nan && rank == origin && (parent(field)[1] = T(NaN))
            _distributed_derivative_failure(origin, comm) do
                input = rank == origin ? bad_field : field
                output = rank == origin && invalid == :workspace ? g.staging : out
                axis = rank == origin && invalid == :axis ? 0 : 1
                deriv!(output, input, g, axis)
            end
            @test isequal(parent(out), original)
            fill!(parent(field), one(T))
            deriv!(out, field, g, 1)
            @test all(iszero, parent(out))
        end
        fill!(parent(field), one(T))
        rank == 0 && (parent(field)[1] = nextfloat(zero(T)))
        rank == ranks - 1 && (parent(field)[end] = floatmax(T) / T(2))
        fill!(parent(out), T(9))
        _distributed_derivative_failure(0, comm) do
            deriv!(out, field, g, 1)
        end
        @test all(==(T(9)), parent(out))
        if ranks > 1
            _distributed_derivative_failure(ranks - 1, comm) do
                distributed_fourier_grid(n, rank == ranks - 1 ? (T(4), L[2], L[3]) : L; comm)
            end
            fill!(parent(field), one(T))
            _distributed_derivative_failure(ranks - 1, comm) do
                deriv!(out, field, g, rank == ranks - 1 ? 2 : 1)
            end
        end
    end
end

@testset "Distributed derivative boundary and transaction contract" begin
    comm = MPI.COMM_WORLD
    rank, ranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    n = (8, 8, 8)
    for origin = 0:ranks-1
        for bad in ((8, 8), (8, 0, 8), [8, 8, 8], (8, 8, 1.5))
            _distributed_derivative_failure(origin, comm) do
                distributed_fourier_grid(rank == origin ? bad : n, (2.0, 3.0, 5.0); comm)
            end
        end
        for bad in (
            (2.0, 3.0),
            [2.0, 3.0, 5.0],
            (2, 3, 5),
            (Float32(2), 3.0, 5.0),
            (NaN, 3.0, 5.0),
            (Inf, 3.0, 5.0),
            (0.0, 3.0, 5.0),
            (-2.0, 3.0, 5.0),
        )
            _distributed_derivative_failure(origin, comm) do
                distributed_fourier_grid(n, rank == origin ? bad : (2.0, 3.0, 5.0); comm)
            end
        end
    end
    @test_throws OverflowError distributed_fourier_grid((typemax(Int), 2, 2), (2.0, 3.0, 5.0); comm)
    @test MPI.Allreduce(1, +, comm) == ranks
    dims = MPI.Dims_create(ranks, (Cint(0), Cint(0)))
    thin = (Int(dims[1]), Int(maximum(dims)), Int(dims[2]))
    for T in (Float32, Float64)
        for lengths in ((nextfloat(zero(T)), T(3), T(5)),)
            @test_throws ArgumentError distributed_fourier_grid(n, lengths; comm)
            @test MPI.Allreduce(1, +, comm) == ranks
        end
        g = distributed_fourier_grid(thin, (T(2), T(3), T(5)); comm)
        field = distributed_field(g)
        fill!(parent(field), one(T))
        for axis = 1:3
            @test all(iszero, parent(deriv(field, g, axis)))
        end
        tiny_length = T === Float32 ? T(1e-30) : T(1e-250)
        g = distributed_fourier_grid(n, (tiny_length, tiny_length, tiny_length); comm)
        field = distributed_field(g)
        fill!(parent(field), floatmax(T) / T(2))
        for axis = 1:3
            @test all(iszero, parent(deriv(field, g, axis)))
        end
        g = distributed_fourier_grid(n, (T(2), T(3), T(5)); comm)
        field, out = distributed_field(g), distributed_field(g)
        for origin = 0:ranks-1, bad_axis in (1.0, nothing, "1", big(typemax(Int)) + 1)
            fill!(parent(field), one(T))
            fill!(parent(out), T(11))
            _distributed_derivative_failure(origin, comm) do
                deriv!(out, field, g, rank == origin ? bad_axis : 1)
            end
            @test all(==(T(11)), parent(out))
        end
        for origin = 0:ranks-1
            fill!(parent(out), T(13))
            _distributed_derivative_failure(origin, comm) do
                deriv!(out, rank == origin ? g.staging : field, g, 1)
            end
            @test all(==(T(13)), parent(out))
        end
        # Different wrapper objects may share valid public physical storage.
        global_field = PencilArrays.global_view(field)
        for I in CartesianIndices(global_field)
            global_field[I] = sin(T(2pi) * T(I[1] - 1) / T(n[1]))
        end
        reference = deriv(field, g, 1)
        alias = PencilArrays.PencilArray(PencilArrays.pencil(field), parent(field))
        @test parent(alias) === parent(field) && Base.mightalias(alias, field)
        @test deriv!(alias, field, g, 1) === alias
        @test parent(alias) == parent(reference)
        if ranks > 1
            _distributed_derivative_failure(ranks - 1, comm) do
                distributed_fourier_grid(
                    n,
                    rank == ranks - 1 ? (T === Float32 ? (2.0, 3.0, 5.0) : (2.0f0, 3.0f0, 5.0f0)) :
                    (T(2), T(3), T(5));
                    comm,
                )
            end
        end
    end
end
