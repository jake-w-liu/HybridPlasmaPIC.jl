using HybridPlasmaPIC, Test

function moment_storage_population(::Type{T}, ::Val{D}, count) where {T,D}
    ps = ParticleSet{D,T}(count)
    for d = 1:D
        ps.x[d] .= T(0.2d)
    end
    ps.weight .= T(0.4)
    for c = 1:3
        ps.v[c] .= T(c + 1)
    end
    return ps
end

moment_storage_arrays(ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
function moment_storage_reject(operation, arrays)
    before = map(copy, arrays)
    @test_throws Union{ArgumentError,DimensionMismatch} operation()
    @test all(isequal(a, b) for (a, b) in zip(arrays, before))
end

@testset "Moment scratch cannot overwrite particle inputs" begin
    for T in (Float32, Float64), D = 1:3, shape in (NGP(), CIC(), TSC()), count in (0, 2)
        g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> one(T), D))
        for component in (:weight, :position, :velocity),
            operation in (:momentum, :current, :single, :multi, :pressure)

            ps = moment_storage_population(T, Val(D), count)
            f = HybridFields{D,T}(g.n)
            fill!(f.n, T(7))
            foreach(a -> fill!(a, T(8)), f.ui)
            ntmp = fill(T(9), g.n)
            mtmp = ntuple(_ -> fill(T(10), g.n), 3)
            P = ntuple(_ -> fill(T(11), g.n), 6)
            work = component === :weight ? ps.weight : component === :position ? ps.x[D] : ps.v[3]
            run = if operation === :momentum
                () -> momentum!(f.ui, ps, g, shape; work)
            elseif operation === :current
                () -> current!(f.ui, ps, g, shape; work)
            elseif operation === :single
                () -> compute_moments!(f, ps, g, shape, T(0.1); work)
            elseif operation === :multi
                () -> compute_moments_multi!(f, [ps], g, shape, T(0.1); ntmp, mtmp, works = [work])
            else
                () -> pressure_tensor!(P, ps, g, shape; work, nbuf = ntmp, mom = mtmp)
            end
            moment_storage_reject(
                run,
                (f.n, f.ui..., ntmp, mtmp..., P..., moment_storage_arrays(ps)...),
            )
        end
    end
end

@testset "Moment operations preflight every grid and population" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> one(T), D))
        count = prod(g.n)
        for failure in (
            :output_particle,
            :output_pair,
            :output_work,
            :late_shape,
            :late_work,
            :cross_work,
            :scratch_output,
            :scratch_pair,
            :matrix_work,
            :late_position,
            :late_velocity_length,
        )
            species = [moment_storage_population(T, Val(D), count) for _ = 1:2]
            f = HybridFields{D,T}(g.n)
            fill!(f.n, T(7))
            foreach(a -> fill!(a, T(8)), f.ui)
            ntmp = fill(T(9), g.n)
            mtmp = ntuple(_ -> fill(T(10), g.n), 3)
            works = Any[fill(T(11), count), fill(T(12), count)]
            if failure === :output_particle
                f.ui = (f.ui[1], f.ui[2], reshape(species[2].weight, g.n))
            elseif failure === :output_pair
                f.ui = (f.ui[1], f.ui[2], f.n)
            elseif failure === :output_work
                works[2] = vec(f.ui[3])
            elseif failure === :late_shape
                f.ui = (f.ui[1], f.ui[2], zeros(T, ntuple(_ -> 2, D)))
            elseif failure === :late_work
                works[2] = zeros(T, count + 1)
            elseif failure === :cross_work
                works[1] = species[2].v[3]
            elseif failure === :scratch_output
                ntmp = f.n
            elseif failure === :scratch_pair
                mtmp = (mtmp[1], mtmp[2], ntmp)
            elseif failure === :matrix_work
                works[2] = zeros(T, count, 1)
            elseif failure === :late_position
                species[2].x[D][end] = T(NaN)
            else
                resize!(species[2].v[3], count - 1)
            end
            arrays = (
                f.n,
                f.ui...,
                ntmp,
                mtmp...,
                works...,
                moment_storage_arrays(species[1])...,
                moment_storage_arrays(species[2])...,
            )
            moment_storage_reject(arrays) do
                compute_moments_multi!(f, species, g, CIC(), T(0.1); ntmp, mtmp, works)
            end
        end
        for operation in (:momentum, :single), failure in (:late_shape, :output_pair, :position)
            ps = moment_storage_population(T, Val(D), count)
            f = HybridFields{D,T}(g.n)
            fill!(f.n, T(7))
            foreach(a -> fill!(a, T(8)), f.ui)
            work = fill(T(9), count)
            if failure === :late_shape
                f.ui = (f.ui[1], f.ui[2], zeros(T, ntuple(_ -> 2, D)))
            elseif failure === :output_pair
                f.ui = (f.ui[1], f.ui[2], f.ui[1])
            else
                ps.x[D][end] = T(Inf)
            end
            moment_storage_reject((f.n, f.ui..., work, moment_storage_arrays(ps)...)) do
                operation === :momentum ? momentum!(f.ui, ps, g, CIC(); work) :
                compute_moments!(f, ps, g, CIC(), T(0.1); work)
            end
        end
    end
end

@testset "Scalar deposition rejects overlapping output" begin
    for T in (Float32, Float64), D = 1:3, source in (:weight, :position, :values)
        g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> one(T), D))
        ps = moment_storage_population(T, Val(D), prod(g.n))
        values = fill(T(2), prod(g.n))
        input = source === :weight ? ps.weight : source === :position ? ps.x[D] : values
        out = reshape(input, g.n)
        moment_storage_reject((out, values, moment_storage_arrays(ps)...)) do
            deposit_scalar!(out, ps, values, g, CIC())
        end
        if source !== :values
            moment_storage_reject((out, moment_storage_arrays(ps)...)) do
                density!(out, ps, g, CIC())
            end
        end
    end
end

@testset "Independent views and sequentially reused scratch preserve moments" begin
    for T in (Float32, Float64), D = 1:3, shape in (NGP(), CIC(), TSC()), count in (0, 3^D)
        g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> one(T), D))
        species = [moment_storage_population(T, Val(D), count) for _ = 1:2]
        for (i, ps) in enumerate(species)
            ps.q = T(i)
            ps.v[1] .*= T(i)
        end
        # Disjoint data views in the same parent remain valid scratch/input.
        backing = zeros(T, 2count + 1)
        copyto!(view(backing, 1:count), species[1].x[1])
        p = species[1]
        viewed = ParticleSet{D,T}(
            ntuple(d -> d == 1 ? view(backing, 1:count) : view(p.x[d], 1:count), D),
            map(a -> view(a, :), p.v),
            view(p.weight, :),
            view(p.id, :),
            view(p.tag, :),
            p.q,
            p.m,
        )
        # Empty views sharing a parent are deliberately rejected by the storage
        # contract. Use independent empty work for the zero-particle case.
        work = count == 0 ? T[] : view(backing, (count+2):(2count+1))
        populations = ParticleSet{D,T}[viewed, species[2]]
        f = HybridFields{D,T}(g.n)
        ntmp = similar(f.n)
        mtmp = map(similar, f.ui)
        before = map(ps -> map(copy, moment_storage_arrays(ps)), populations)
        compute_moments_multi!(f, populations, g, shape, T(0.01); ntmp, mtmp, works = (work, work))
        expected_charge = sum(ps.q * sum(ps.weight) for ps in populations)
        volume = prod(g.dx)
        tolerance = T(256) * eps(T)
        @test sum(f.n) * volume ≈ expected_charge rtol = tolerance atol = tolerance
        for c = 1:3
            expected = sum(ps.q * sum(ps.weight .* ps.v[c]) for ps in populations)
            @test sum(f.ui[c] .* max.(f.n, T(0.01))) * volume ≈ expected rtol = tolerance atol =
                tolerance
        end
        @test all(
            all(isequal(a, b) for (a, b) in zip(moment_storage_arrays(ps), old)) for
            (ps, old) in zip(populations, before)
        )
        for ps in populations
            momentum!(mtmp, ps, g, shape; work)
            for c = 1:3
                @test sum(mtmp[c]) * volume ≈ sum(ps.weight .* ps.v[c]) rtol = tolerance atol =
                    tolerance
            end
            current!(mtmp, ps, g, shape; work)
            for c = 1:3
                @test sum(mtmp[c]) * volume ≈ ps.q * sum(ps.weight .* ps.v[c]) rtol = tolerance atol =
                    tolerance
            end
        end
    end
end
