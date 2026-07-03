# KDV-001 — KdV soliton: shape, propagation speed, conserved quantities, and
# dealiasing sensitivity, against the analytic sech² soliton.

using HybridPlasmaPIC, Test

struct OffsetKDVVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    firstindex::Int
end
Base.size(v::OffsetKDVVector) = size(v.data)
Base.axes(v::OffsetKDVVector) = (v.firstindex:(v.firstindex+length(v.data)-1),)
Base.IndexStyle(::Type{<:OffsetKDVVector}) = IndexLinear()
Base.getindex(v::OffsetKDVVector, i::Int) = v.data[i-v.firstindex+1]

@testset "KDV input validation and vector axes" begin
    u0 = [1.0, 0.5, -0.25, 0.125, -0.0625, 0.0, 0.25, -0.5]
    @test_throws ArgumentError kdv_solve(u0, 1.0, 1.0, 1.0, 1.0, 0.01, 1.5)
    @test_throws ArgumentError kdv_solve(u0, Inf, 1.0, 1.0, 1.0, 0.01, 1)
    @test_throws ArgumentError kdv_solve(u0, 1.0, 1.0, 1.0, 1.0, NaN, 1)
    @test_throws ArgumentError kdv_solve(u0, 1.0, 1.0, 1.0, 1.0, Inf, 1)
    @test_throws ArgumentError kdv_solve([1.0, NaN, 0.0], 1.0, 1.0, 1.0, 1.0, 0.01, 1)
    @test_throws ArgumentError kdv_solve([1.0, Inf, 0.0], 1.0, 1.0, 1.0, 1.0, 0.01, 0)

    uoff = OffsetKDVVector(u0, -3)
    @test kdv_solve(uoff, 1.0, 1.0, 1.0, 1.0, 0.01, 0) == u0
    @test kdv_solve(uoff, 1.0, 1.0, 1.0, 1.0, 0.01, 2) ≈ kdv_solve(u0, 1.0, 1.0, 1.0, 1.0, 0.01, 2)
end

@testset "KDV-001 soliton" begin
    Ld = 40.0
    n = 512
    x = collect(range(0, Ld; length = n + 1)[1:n])
    c0 = 0.0
    α = 6.0
    β = 1.0
    A = 1.0
    x0 = 10.0
    V = c0 + α * A / 3
    u0 = kdv_soliton(A, c0, α, β, x, 0.0, x0, Ld)
    Tend = 4.0
    dt = 0.004
    ns = round(Int, Tend / dt)

    uf = kdv_solve(u0, Ld, c0, α, β, dt, ns)
    ua = kdv_soliton(A, c0, α, β, x, Tend, x0, Ld)

    shape_err = sqrt(sum((uf .- ua) .^ 2) / sum(ua .^ 2))
    @test shape_err < 1e-3                              # shape preserved
    @test abs(maximum(uf) - A) / A < 1e-2              # amplitude preserved
    dx = Ld / n
    @test abs(x[argmax(uf)] - mod(x0 + V * Tend, Ld)) < 2dx   # propagation speed

    # conserved quantities
    @test abs(sum(uf) - sum(u0)) / abs(sum(u0)) < 1e-8        # ∫u dx (mass, exact at k=0)
    @test abs(sum(uf .^ 2) - sum(u0 .^ 2)) / sum(u0 .^ 2) < 1e-3  # ∫u² dx (momentum)

    # dealiasing sensitivity: removing the 2/3 filter degrades the shape
    uf_nd = kdv_solve(u0, Ld, c0, α, β, dt, ns; dealias = false)
    shape_err_nd = sqrt(sum((uf_nd .- ua) .^ 2) / sum(ua .^ 2))
    @test shape_err <= shape_err_nd
end

@testset "KDV-001 long-time phase error" begin
    Ld = 40.0
    n = 512
    x = collect(range(0, Ld; length = n + 1)[1:n])
    c0 = 0.0
    α = 6.0
    β = 1.0
    A = 1.0
    x0 = 10.0
    V = c0 + α * A / 3
    u0 = kdv_soliton(A, c0, α, β, x, 0.0, x0, Ld)
    Tend = 15.0
    dt = 0.004
    ns = round(Int, Tend / dt)   # soliton travels 30 (wraps the box once)
    uf = kdv_solve(u0, Ld, c0, α, β, dt, ns)
    ua = kdv_soliton(A, c0, α, β, x, Tend, x0, Ld)
    @test sqrt(sum((uf .- ua) .^ 2) / sum(ua .^ 2)) < 3e-3   # shape held over long time
    @test abs(maximum(uf) - A) / A < 1e-2                    # amplitude held
    dx = Ld / n
    @test abs(x[argmax(uf)] - mod(x0 + V * Tend, Ld)) < 2dx  # phase (position) error bounded
end
