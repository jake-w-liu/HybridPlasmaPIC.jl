# Deposition/gather benchmarks DEP-001..005 plus deposit/gather adjoint
# consistency. Oracles: partition of unity, exact weight sum, linear-field
# reproduction, translation symmetry, and the N^{-1/2} noise law.

using HybridPlasmaPIC, Test, LinearAlgebra, Random, Statistics

struct OffsetTestVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    first_index::Int
end

Base.size(v::OffsetTestVector) = size(v.data)
Base.axes(v::OffsetTestVector) = (v.first_index:(v.first_index+length(v.data)-1),)
Base.IndexStyle(::Type{<:OffsetTestVector}) = IndexLinear()
Base.getindex(v::OffsetTestVector, i::Int) = v.data[i-v.first_index+1]
Base.setindex!(v::OffsetTestVector, x, i::Int) = (v.data[i-v.first_index+1] = x)

@testset "DEP-001 partition of unity" begin
    T = Float64
    for shape in (NGP(), CIC(), TSC())
        for s in range(-2.3, 5.7; length = 41)
            st = HybridPlasmaPIC._stencil1d(shape, T(s))
            @test sum(st[2]) ≈ 1
        end
    end
    # global: deposit unit weights → Σ_g field = N
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 12, D)
        g = FourierGrid(n, ntuple(_ -> 1.0, D))
        N = 200
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(D), ntuple(_ -> 0.0, D), ntuple(_ -> 1.0, D))
        field = Array{T,D}(undef, n)
        deposit_scalar!(field, ps, ps.weight, g, shape)
        @test sum(field) ≈ N
    end
end

@testset "DEP-002 global particle number" begin
    T = Float64
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 10, D)
        g = FourierGrid(n, ntuple(_ -> 2.0, D))
        N = 300
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(7), ntuple(_ -> 0.0, D), ntuple(_ -> 2.0, D))
        ps.weight .= 0.5 .+ rand(MersenneTwister(8), N)        # nonuniform weights
        nfield = Array{T,D}(undef, n)
        density!(nfield, ps, g, shape)
        @test sum(nfield) * prod(g.dx) ≈ sum(ps.weight)
    end
end

@testset "DEP-003 linear field gather (CIC exact)" begin
    T = Float64
    for D = 1:3
        n = ntuple(_ -> 16, D)
        g = FourierGrid(n, ntuple(_ -> 1.0, D))
        a = 0.5
        b = ntuple(_ -> 1.3, D)
        field = Array{T,D}(undef, n)
        for I in CartesianIndices(field)
            t = Tuple(I)
            field[I] = a + sum(b[d] * (t[d] - 1) * g.dx[d] for d = 1:D)
        end
        M = 50
        ps = ParticleSet{D,T}(M)
        rng = MersenneTwister(3)
        for d = 1:D, p = 1:M
            ps.x[d][p] = (3 + 10 * rand(rng)) * g.dx[d]        # interior: no seam wrap
        end
        out = zeros(T, M)
        gather_scalar!(out, field, ps, g, CIC())
        exact = [a + sum(b[d] * ps.x[d][p] for d = 1:D) for p = 1:M]
        @test maximum(abs.(out .- exact)) < 1e-12
    end
end

@testset "DEP-004 translation invariance" begin
    T = Float64
    ncell = 32
    L = 1.0
    g = FourierGrid((ncell,), (L,))
    for shape in (NGP(), CIC(), TSC())
        N = 500
        ps = ParticleSet{1,T}(N)
        load_uniform!(ps, MersenneTwister(5), (0.0,), (L,))
        n1 = zeros(T, ncell)
        density!(n1, ps, g, shape)
        for p = 1:N                                           # shift by exactly one cell
            ps.x[1][p] = mod(ps.x[1][p] + g.dx[1], L)
        end
        n2 = zeros(T, ncell)
        density!(n2, ps, g, shape)
        @test maximum(abs.(n2 .- circshift(n1, 1))) < 1e-12
    end
end

@testset "DEP-005 noise scaling N^(-1/2)" begin
    T = Float64
    ncell = 64
    L = 1.0
    g = FourierGrid((ncell,), (L,))
    nppcs = [16, 64, 256, 1024]
    seeds = 1:8
    rms = Float64[]
    for nppc in nppcs
        N = nppc * ncell
        acc = 0.0
        for sd in seeds
            ps = ParticleSet{1,T}(N)
            load_uniform!(ps, MersenneTwister(sd * 1000 + nppc), (0.0,), (L,))
            nf = zeros(T, ncell)
            density!(nf, ps, g, CIC())
            acc += std(nf) / mean(nf)
        end
        push!(rms, acc / length(seeds))
    end
    slope = (log(rms[end]) - log(rms[1])) / (log(nppcs[end]) - log(nppcs[1]))
    @test -0.6 < slope < -0.4
end

@testset "deposit/gather adjoint consistency" begin
    T = Float64
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 8, D)
        g = FourierGrid(n, ntuple(_ -> 1.0, D))
        N = 100
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(11), ntuple(_ -> 0.0, D), ntuple(_ -> 1.0, D))
        vals = randn(MersenneTwister(12), N)
        field = randn(MersenneTwister(13), n...)
        dep = Array{T,D}(undef, n)
        deposit_scalar!(dep, ps, vals, g, shape)
        gat = zeros(T, N)
        gather_scalar!(gat, field, ps, g, shape)
        @test dot(vec(dep), vec(field)) ≈ dot(vals, gat)       # ⟨deposit·field⟩ = ⟨vals·gather⟩
    end
end

@testset "deposit/gather defensive shape checks" begin
    T = Float64
    g = FourierGrid((4, 4), (1.0, 1.0))
    ps = ParticleSet{2,T}(1)
    ps.x[1][1] = 0.2
    ps.x[2][1] = 0.3
    vals = ones(T, 1)

    @test_throws DimensionMismatch deposit_scalar!(zeros(T, 3, 4), ps, vals, g, CIC())
    @test_throws DimensionMismatch gather_scalar!(zeros(T, 1), zeros(T, 3, 4), ps, g, CIC())
end

@testset "deposit/gather reject non-finite particle positions" begin
    T = Float64
    g = FourierGrid((4,), (1.0,))
    vals = T[1.0]
    field = ones(T, g.n)

    for shape in (NGP(), CIC(), TSC())
        ps = ParticleSet{1,T}(1)
        ps.x[1][1] = T(NaN)
        @test_throws ArgumentError deposit_scalar!(zeros(T, g.n), ps, vals, g, shape)
        @test_throws ArgumentError gather_scalar!(zeros(T, 1), field, ps, g, shape)

        ps.x[1][1] = T(Inf)
        @test_throws ArgumentError density!(zeros(T, g.n), ps, g, shape)
        @test_throws ArgumentError momentum!(ntuple(_ -> zeros(T, g.n), 3), ps, g, shape)
        @test_throws ArgumentError current!(ntuple(_ -> zeros(T, g.n), 3), ps, g, shape)
        @test_throws ArgumentError pressure_tensor!(ntuple(_ -> zeros(T, g.n), 6), ps, g, shape)
    end
end

@testset "deposit/gather AbstractVector indices are decoupled from particle indices" begin
    T = Float64
    g = FourierGrid((8,), (1.0,))
    ps = ParticleSet{1,T}(3)
    ps.x[1] .= T[0.10, 0.35, 0.90]

    vals_data = T[1.0, 2.0, 3.0]
    vals_offset = OffsetTestVector(vals_data, -2)
    dep_plain = zeros(T, g.n)
    dep_offset = zeros(T, g.n)
    deposit_scalar!(dep_plain, ps, vals_data, g, CIC())
    deposit_scalar!(dep_offset, ps, vals_offset, g, CIC())
    @test dep_offset ≈ dep_plain

    field = randn(MersenneTwister(42), g.n...)
    out_plain = zeros(T, nparticles(ps))
    out_data = fill(T(NaN), nparticles(ps))
    out_offset = OffsetTestVector(out_data, 5)
    gather_scalar!(out_plain, field, ps, g, CIC())
    gather_scalar!(out_offset, field, ps, g, CIC())
    @test out_data ≈ out_plain
end

@testset "deposit/gather allocate nothing in steady state" begin
    T = Float64
    g = FourierGrid((16, 16), (1.0, 1.0))
    ps = ParticleSet{2,T}(100)
    load_uniform!(ps, MersenneTwister(21), (0.0, 0.0), (1.0, 1.0))
    vals = copy(ps.weight)
    field = randn(MersenneTwister(22), g.n...)
    dep = zeros(T, g.n)
    out = zeros(T, nparticles(ps))
    deposit_scalar!(dep, ps, vals, g, CIC())
    gather_scalar!(out, field, ps, g, CIC())
    @test (@allocated deposit_scalar!(dep, ps, vals, g, CIC())) == 0
    @test (@allocated gather_scalar!(out, field, ps, g, CIC())) == 0
end
