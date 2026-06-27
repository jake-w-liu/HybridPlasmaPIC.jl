# Operator benchmarks OP-001/002/003/005/006 (periodic Fourier).
# OP-004 (mixed FD/Fourier) is deferred to the mixed-operator phase.
# Oracles are closed-form and never call the production operator under test.

using HybridPlasmaPIC, Test, LinearAlgebra, Random, Statistics

relL2(a, b) = norm(a .- b) / norm(b)

# Build f = sin(kx x) cos(ky y) cos(kz z) and exact ∂x f on a 2π^D grid (integer modes).
function sincos_field(::Type{T}, n::NTuple{D,Int}, modes::NTuple{D,Int}) where {T,D}
    L = ntuple(_ -> T(2π), D)
    g = FourierGrid(n, L)
    f = Array{T,D}(undef, n)
    dfx = similar(f)
    @inbounds for I in CartesianIndices(f)
        t = Tuple(I)
        x = (t[1] - 1) * g.dx[1]
        val = sin(modes[1] * x)
        dval = modes[1] * cos(modes[1] * x)
        if D >= 2
            y = (t[2] - 1) * g.dx[2]
            val *= cos(modes[2] * y)
            dval *= cos(modes[2] * y)
        end
        if D >= 3
            z = (t[3] - 1) * g.dx[3]
            val *= cos(modes[3] * z)
            dval *= cos(modes[3] * z)
        end
        f[I] = val
        dfx[I] = dval
    end
    return g, f, dfx
end

@testset "OP-001 periodic scalar derivative" begin
    for (T, tol) in ((Float64, 1e-11), (Float32, 1e-5))
        for (n, modes) in (((16,), (3,)), ((16, 16), (3, 2)), ((12, 16, 8), (2, 3, 2)))
            g, f, dfx = sincos_field(T, n, modes)
            out = similar(f)
            deriv!(out, f, g, 1)
            e = relL2(out, dfx)
            @test e < tol
        end
    end
end

@testset "OP-002 div(curl)=0" begin
    Random.seed!(1)
    for T in (Float64, Float32)
        tol = T == Float64 ? 1e-10 : 1e-3
        for n in ((16,), (16, 16), (8, 12, 16))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            A = ntuple(_ -> randn(T, n...), 3)
            B = ntuple(_ -> similar(first(A)), 3)
            curl!(B, A, g)
            divB = similar(B[1])
            divergence!(divB, B, g)
            kmax = maximum(maximum(abs, g.kvec[d]) for d = 1:D)
            resid = norm(divB) / (kmax * norm(B[1]) + eps(T))
            @test resid < tol
        end
    end
end

@testset "OP-003 divergence-free projection" begin
    Random.seed!(2)
    for T in (Float64, Float32)
        tol = T == Float64 ? 1e-10 : 1e-3
        for n in ((16,), (16, 16), (8, 12, 16))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            # transverse part (already divergence-free) + uniform background
            A = ntuple(_ -> randn(T, n...), 3)
            Btrans = ntuple(_ -> similar(first(A)), 3)
            curl!(Btrans, A, g)
            bg = (T(0.7), T(-0.3), T(1.1))
            target = ntuple(c -> Btrans[c] .+ bg[c], 3)
            # add a longitudinal (gradient) part that projection must remove
            φ = randn(T, n...)
            gradφ = ntuple(_ -> similar(φ), D)
            gradient!(gradφ, φ, g)
            B = ntuple(c -> copy(target[c]), 3)
            for c = 1:D
                B[c] .+= gradφ[c]
            end
            project_divfree!(B, g)
            # longitudinal removed → recover transverse + background
            err = norm(B[1] .- target[1]) + norm(B[2] .- target[2]) + norm(B[3] .- target[3])
            scale = norm(target[1]) + norm(target[2]) + norm(target[3])
            @test err / scale < tol
            # result is divergence-free
            divB = similar(B[1])
            divergence!(divB, B, g)
            kmax = maximum(maximum(abs, g.kvec[d]) for d = 1:D)
            @test norm(divB) / (kmax * scale + eps(T)) < tol
            # k=0 mean preserved on every component
            for c = 1:3
                @test isapprox(mean(B[c]), bg[c]; atol = (T == Float64 ? 1e-10 : 1e-3))
            end
        end
    end
end

@testset "OP-005 discrete integration by parts" begin
    Random.seed!(3)
    for T in (Float64, Float32)
        tol = T == Float64 ? 1e-10 : 1e-4
        for n in ((16,), (16, 16), (8, 12, 16))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            f = randn(T, n...)
            h = randn(T, n...)
            for j = 1:D
                lhs = dot(f, deriv(h, g, j))
                rhs = -dot(deriv(f, g, j), h)
                @test abs(lhs - rhs) / (abs(lhs) + abs(rhs) + eps(T)) < tol
            end
        end
    end
end

@testset "OP-006 Nyquist handling" begin
    for T in (Float64, Float32)
        n = (16,)
        g = FourierGrid(n, (T(2π),))
        # pure Nyquist mode cos(N/2 · 2π x / L) = (-1)^i : unrepresentable derivative → 0
        f = T[(-1.0)^(i - 1) for i = 1:n[1]]
        out = deriv!(similar(f), f, g, 1)
        @test maximum(abs, out) < (T == Float64 ? 1e-10 : 1e-4)
        # real input → real output (deriv! returns a real array by construction)
        @test eltype(out) == T
        # a resolved mode is still differentiated correctly
        g2, fr, dfx = sincos_field(T, n, (3,))
        @test relL2(deriv(fr, g2, 1), dfx) < (T == Float64 ? 1e-11 : 1e-5)
    end
end

@testset "operators allocate nothing in steady state" begin
    T = Float64
    n = (16, 16, 16)
    g = FourierGrid(n, ntuple(_ -> T(2π), 3))
    f = randn(T, n...)
    out = similar(f)
    grad = ntuple(_ -> similar(f), 3)
    A = ntuple(_ -> randn(T, n...), 3)
    B = ntuple(_ -> similar(f), 3)
    divB = similar(f)
    deriv!(out, f, g, 1)
    gradient!(grad, f, g)
    curl!(B, A, g)
    divergence!(divB, B, g)
    project_divfree!(B, g)   # warm up
    @test (@allocated deriv!(out, f, g, 1)) == 0
    @test (@allocated gradient!(grad, f, g)) == 0
    @test (@allocated curl!(B, A, g)) == 0
    @test (@allocated divergence!(divB, B, g)) == 0
    @test (@allocated project_divfree!(B, g)) == 0
end

@testset "laplacian! includes the Nyquist mode (even N)" begin
    # Regression: ∇² is an even derivative, so the Nyquist mode must get −k_nyq²·f,
    # not 0. The first-derivative kvec zeroes Nyquist; laplacian! must not reuse it.
    g = FourierGrid((8,), (2π,))
    f = Float64[(-1)^(i - 1) for i = 1:8]          # pure Nyquist mode
    out = similar(f)
    laplacian!(out, f, g)
    knyq = 2π * (8 ÷ 2) / 2π                        # = 4
    @test out ≈ (-knyq^2) .* f rtol = 1e-10        # was identically 0 before the fix
end
