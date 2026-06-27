# Tests for the explicit binomial (1,2,1)/4 smoothing operator.
# Oracle: a pure cosine mode cos(k x_g) on a periodic grid is an exact
# eigenfunction of the stencil with eigenvalue cos²(k dx/2) per pass per axis.
# The oracle (smoothing_transfer) is closed-form and independent of the
# in-place stencil implementation under test.

using HybridPlasmaPIC, Test, LinearAlgebra, Random

@testset "smoothing_transfer analytic form" begin
    for T in (Float64, Float32)
        dx = T(0.1)
        for k in T.((0.0, 0.5, 1.3, 3.7))
            # matches cos²(k dx/2)
            @test smoothing_transfer(k, dx; passes = 1) ≈ cos(k * dx / 2)^2
            # passes compose multiplicatively
            for p in (0, 1, 2, 5)
                @test smoothing_transfer(k, dx; passes = p) ≈
                      smoothing_transfer(k, dx; passes = 1)^p
            end
        end
        # k=0 → exactly 1 for any number of passes (mean/total preserved)
        @test smoothing_transfer(zero(T), dx; passes = 3) == one(T)
        # Nyquist k dx = π → fully removed
        @test smoothing_transfer(T(π) / dx, dx; passes = 1) ≈ zero(T) atol = 1e-12
        # passes = 0 is the identity (transfer ≡ 1)
        @test smoothing_transfer(T(2.0), dx; passes = 0) == one(T)
        @test_throws ArgumentError smoothing_transfer(T(1.0), dx; passes = -1)
    end
end

@testset "1-D single mode multiplied by transfer to ~1e-12" begin
    T = Float64
    n = (32,)
    L = (T(2π),)
    g = FourierGrid(n, L)
    dx = g.dx[1]
    xg = [(i - 1) * dx for i = 1:n[1]]   # grid points
    for mode in (1, 2, 3, 7)               # integer modes resolved on the grid
        k = T(mode)                        # since L = 2π, k = integer mode number
        for passes in (1, 2, 4)
            f = cos.(k .* xg)
            f0 = copy(f)
            binomial_smooth!(f, g; passes = passes)
            factor = smoothing_transfer(k, dx; passes = passes)
            @test maximum(abs.(f .- factor .* f0)) < 1e-12
        end
    end
end

@testset "k=0 uniform field unchanged & total conserved" begin
    for T in (Float64, Float32)
        for n in ((24,), (16, 20), (8, 10, 12))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            c = T(2.5)
            f = fill(c, n...)
            f0 = copy(f)
            binomial_smooth!(f, g; passes = 3)
            @test f == f0                              # uniform field exactly unchanged
            # total (k=0 content) conserved for an arbitrary field too
            Random.seed!(7)
            h = randn(T, n...)
            s_before = sum(h)
            binomial_smooth!(h, g; passes = 2)
            @test isapprox(sum(h), s_before; rtol = (T == Float64 ? 1e-12 : 1e-5))
        end
    end
end

@testset "multi-D separable transfer (product over axes)" begin
    T = Float64
    n = (16, 24)
    L = ntuple(_ -> T(2π), 2)
    g = FourierGrid(n, L)
    kx, ky = T(3), T(5)
    f = Array{T,2}(undef, n)
    for I in CartesianIndices(f)
        i, j = Tuple(I)
        f[I] = cos(kx * (i - 1) * g.dx[1]) * cos(ky * (j - 1) * g.dx[2])
    end
    f0 = copy(f)
    passes = 2
    binomial_smooth!(f, g; passes = passes)
    factor =
        smoothing_transfer(kx, g.dx[1]; passes = passes) *
        smoothing_transfer(ky, g.dx[2]; passes = passes)
    @test maximum(abs.(f .- factor .* f0)) < 1e-12
end

@testset "passes = 0 is identity; argument checking" begin
    T = Float64
    n = (16, 16)
    g = FourierGrid(n, ntuple(_ -> T(2π), 2))
    f = randn(T, n...)
    f0 = copy(f)
    @test binomial_smooth!(f, g; passes = 0) === f
    @test f == f0
    @test_throws ArgumentError binomial_smooth!(f, g; passes = -1)
    bad = randn(T, 8, 8)
    @test_throws DimensionMismatch binomial_smooth!(bad, g; passes = 1)
end

@testset "filter never amplifies (transfer in [0,1])" begin
    T = Float64
    n = (40,)
    g = FourierGrid(n, (T(2π),))
    dx = g.dx[1]
    Random.seed!(11)
    f = randn(T, n...)
    e_before = norm(f)
    binomial_smooth!(f, g; passes = 1)
    @test norm(f) <= e_before + 1e-12
    for m = 0:(n[1]÷2)
        @test 0 <= smoothing_transfer(T(m), dx; passes = 1) <= 1
    end
end
