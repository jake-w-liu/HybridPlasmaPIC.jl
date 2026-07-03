# Threaded deposition (§21.4). Oracle: deposit_scalar_threaded! reproduces the
# serial deposit_scalar! for every D∈{1,2,3} and every shape, to within an
# rtol tolerating only floating-point reduction-order differences; conservation
# (Σ out = Σ vals) is preserved; density_threaded! matches density!.

using HybridPlasmaPIC, Test, LinearAlgebra, Random

@info "test_threaded: Threads.nthreads() = $(Threads.nthreads())"

struct ThreadOffsetVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    first_index::Int
end

Base.size(v::ThreadOffsetVector) = size(v.data)
Base.axes(v::ThreadOffsetVector) = (v.first_index:(v.first_index+length(v.data)-1),)
Base.IndexStyle(::Type{<:ThreadOffsetVector}) = IndexLinear()
Base.getindex(v::ThreadOffsetVector, i::Int) = v.data[i-v.first_index+1]
Base.setindex!(v::ThreadOffsetVector, x, i::Int) = (v.data[i-v.first_index+1] = x)

@testset "THR-001 threaded == serial deposit_scalar!" begin
    T = Float64
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 12, D)
        g = FourierGrid(n, ntuple(_ -> 1.7, D))
        N = 5000
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(100 + D), ntuple(_ -> 0.0, D), ntuple(d -> g.L[d], D))
        vals = randn(MersenneTwister(200 + D), N)

        serial = Array{T,D}(undef, n)
        deposit_scalar!(serial, ps, vals, g, shape)
        thr = Array{T,D}(undef, n)
        deposit_scalar_threaded!(thr, ps, vals, g, shape)

        # Match to within reduction-order roundoff.
        @test isapprox(thr, serial; rtol = 1e-10, atol = 1e-12)
        @test maximum(abs.(thr .- serial)) < 1e-10 * (maximum(abs.(serial)) + 1)
        # Conservation: total deposited weight unchanged (partition of unity).
        @test sum(thr) ≈ sum(vals)
        @test sum(thr) ≈ sum(serial)
    end
end

@testset "THR-001b threaded deposit accepts offset-indexed values" begin
    T = Float64
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 8, D)
        g = FourierGrid(n, ntuple(_ -> 1.0, D))
        N = 127
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(300 + D), ntuple(_ -> 0.0, D), ntuple(d -> g.L[d], D))
        vals = randn(MersenneTwister(400 + D), N)
        vals_offset = ThreadOffsetVector(vals, -5)

        serial = zeros(T, n)
        threaded = zeros(T, n)
        deposit_scalar!(serial, ps, vals, g, shape)
        deposit_scalar_threaded!(threaded, ps, vals_offset, g, shape)

        @test threaded ≈ serial
        @test sum(threaded) ≈ sum(vals)
    end
end

@testset "THR-002 nonuniform weights + interior particles" begin
    T = Float64
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 9, D)
        g = FourierGrid(n, ntuple(_ -> 0.3, D))
        N = 3333   # not divisible by typical nthreads → exercises ranged partition
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(11 * D), ntuple(_ -> 0.0, D), ntuple(d -> g.L[d], D))
        vals = 0.5 .+ rand(MersenneTwister(17 * D), N)

        serial = zeros(T, n)
        deposit_scalar!(serial, ps, vals, g, shape)
        thr = zeros(T, n)
        deposit_scalar_threaded!(thr, ps, vals, g, shape)
        @test isapprox(thr, serial; rtol = 1e-10, atol = 1e-12)
        @test sum(thr) ≈ sum(vals)
    end
end

@testset "THR-003 density_threaded! == density!" begin
    T = Float64
    for D = 1:3, shape in (NGP(), CIC(), TSC())
        n = ntuple(_ -> 10, D)
        g = FourierGrid(n, ntuple(_ -> 2.0, D))
        N = 2000
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(900 + D), ntuple(_ -> 0.0, D), ntuple(d -> g.L[d], D))
        ps.weight .= 0.5 .+ rand(MersenneTwister(901 + D), N)

        nser = Array{T,D}(undef, n)
        density!(nser, ps, g, shape)
        nthr = Array{T,D}(undef, n)
        density_threaded!(nthr, ps, g, shape)
        @test isapprox(nthr, nser; rtol = 1e-10, atol = 1e-12)
        @test sum(nthr) * prod(g.dx) ≈ sum(ps.weight)
    end
end

@testset "THR-004 edge cases (N=0, single particle)" begin
    T = Float64
    D = 2
    n = ntuple(_ -> 8, D)
    g = FourierGrid(n, ntuple(_ -> 1.0, D))
    for shape in (NGP(), CIC(), TSC())
        # Zero particles → all zeros, conserves trivially.
        ps0 = ParticleSet{D,T}(0)
        z = ones(T, n)
        deposit_scalar_threaded!(z, ps0, T[], g, shape)
        @test all(==(zero(T)), z)

        # One particle.
        ps1 = ParticleSet{D,T}(1)
        ps1.x[1][1] = 3.2
        ps1.x[2][1] = 5.1
        v1 = T[2.5]
        s = zeros(T, n)
        deposit_scalar!(s, ps1, v1, g, shape)
        t = zeros(T, n)
        deposit_scalar_threaded!(t, ps1, v1, g, shape)
        @test isapprox(t, s; rtol = 1e-12, atol = 1e-14)
        @test sum(t) ≈ 2.5
    end
end

# DimensionMismatch guard mirrors the serial routine.
@testset "THR-005 length guard" begin
    T = Float64
    g = FourierGrid((6,), (1.0,))
    ps = ParticleSet{1,T}(10)
    @test_throws DimensionMismatch deposit_scalar_threaded!(zeros(T, 6), ps, zeros(T, 9), g, CIC())
    @test_throws DimensionMismatch deposit_scalar_threaded!(zeros(T, 5), ps, zeros(T, 10), g, CIC())
    @test_throws DimensionMismatch density_threaded!(zeros(T, 5), ps, g, CIC())
end

@testset "THR-005b non-finite particle positions are rejected" begin
    T = Float64
    g = FourierGrid((6,), (1.0,))
    ps = ParticleSet{1,T}(1)
    ps.x[1][1] = T(Inf)
    @test_throws ArgumentError deposit_scalar_threaded!(zeros(T, g.n), ps, T[1.0], g, CIC())
    @test_throws ArgumentError density_threaded!(zeros(T, g.n), ps, g, CIC())
end

# When multiple threads are active, the threaded path (not the serial fallback)
# must still match serial — this asserts the partial-accumulator reduction.
if Threads.nthreads() > 1
    @testset "THR-006 multi-thread path matches serial" begin
        T = Float64
        D = 3
        shape = TSC()
        n = ntuple(_ -> 11, D)
        g = FourierGrid(n, ntuple(_ -> 1.0, D))
        N = 20000
        ps = ParticleSet{D,T}(N)
        load_uniform!(ps, MersenneTwister(42), ntuple(_ -> 0.0, D), ntuple(d -> g.L[d], D))
        vals = randn(MersenneTwister(43), N)
        s = zeros(T, n)
        deposit_scalar!(s, ps, vals, g, shape)
        t = zeros(T, n)
        deposit_scalar_threaded!(t, ps, vals, g, shape)
        @test isapprox(t, s; rtol = 1e-10, atol = 1e-12)
        @test sum(t) ≈ sum(vals)
    end
end
