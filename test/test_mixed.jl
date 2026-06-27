# OP-004 / OP-005(SBP) — mixed FD/Fourier shock operator. SBP summation-by-parts
# identity, x-grid convergence of the non-periodic SBP derivative, and the 2-D
# mixed (SBP-x, Fourier-y) manufactured solution.

using HybridPlasmaPIC, Test, LinearAlgebra, Random

@testset "OP-005 SBP summation by parts" begin
    T = Float64
    for n in (17, 33, 64)
        s = SBP1D(n, 1.0)
        Random.seed!(n)
        u = randn(n)
        v = randn(n)
        Du = sbp_deriv(u, s)
        Dv = sbp_deriv(v, s)
        lhs = sum(s.H .* u .* Dv) + sum(s.H .* Du .* v)
        rhs = u[n] * v[n] - u[1] * v[1]               # uᵀ B v
        @test abs(lhs - rhs) < 1e-12
    end
end

@testset "OP-004 SBP derivative convergence" begin
    L = 1.0
    f(x) = sin(2x) + 0.5x^2
    df(x) = 2cos(2x) + x                              # non-periodic on [0,L]
    errs = Float64[]
    ns = [33, 65, 129, 257]
    for n in ns
        s = SBP1D(n, L)
        x = range(0, L; length = n)
        fa = f.(x)
        D = sbp_deriv(collect(fa), s)
        e = sqrt(sum(s.H .* (D .- df.(x)) .^ 2))      # H-norm error
        push!(errs, e)
    end
    # diagonal-norm SBP-(2,1): 2nd-order interior, 1st-order one-sided boundary →
    # H-norm derivative error converges at the (well-known) rate 3/2
    rate = log(errs[1] / errs[end]) / log((ns[end] - 1) / (ns[1] - 1))
    @test 1.3 < rate < 1.7
end

@testset "OP-004 mixed SBP-x / Fourier-y manufactured solution" begin
    T = Float64
    Lx = 1.0
    Ly = 2π
    ny = 16
    ky = 2.0
    g(x) = sin(2x) + 0.5x^2
    dg(x) = 2cos(2x) + x
    # ∂y is spectrally exact; check it once, then refine x for ∂x convergence
    errs_x = Float64[]
    ns = [33, 65, 129]
    for nx in ns
        s = SBP1D(nx, Lx)
        x = collect(range(0, Lx; length = nx))
        y = [(j - 1) * Ly / ny for j = 1:ny]
        F = [g(x[i]) * cos(ky * y[j]) for i = 1:nx, j = 1:ny]
        # ∂x via SBP
        Fx = similar(F)
        sbp_deriv_x!(Fx, F, s)
        Fx_exact = [dg(x[i]) * cos(ky * y[j]) for i = 1:nx, j = 1:ny]
        ex = sqrt(sum(s.H[i] * (Fx[i, j] - Fx_exact[i, j])^2 for i = 1:nx, j = 1:ny))
        push!(errs_x, ex)
        # ∂y via Fourier (spectrally exact, band-limited)
        Fy = similar(F)
        fourier_deriv_y!(Fy, F, Ly)
        Fy_exact = [-ky * g(x[i]) * sin(ky * y[j]) for i = 1:nx, j = 1:ny]
        @test maximum(abs, Fy .- Fy_exact) < 1e-11
    end
    rate = log(errs_x[1] / errs_x[end]) / log((ns[end] - 1) / (ns[1] - 1))
    @test 1.3 < rate < 1.7                            # SBP-(2,1) H-norm rate 3/2
end
