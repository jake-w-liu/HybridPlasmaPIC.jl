module TestCoulombWeightMean
using HybridPlasmaPIC, Random, Test
const H = HybridPlasmaPIC

# Count equal native inputs before constructing the independent exact mean.
# This keeps large, sparse reference populations inexpensive without reusing
# the implementation's accumulation order or floating-point accumulator.
function exact_mean(weights)
    scale = maximum(weights)
    counts = Dict{eltype(weights),Int}()
    for weight in weights
        normalized = weight / scale
        counts[normalized] = get(counts, normalized, 0) + 1
    end
    return sum(Rational{BigInt}(value) * count for (value, count) in counts) / length(weights)
end

function population_weights(::Type{T}, n, pattern) where {T}
    w = ones(T, n)
    if pattern === :sparse
        fill!(w, zero(T))
        w[1] = one(T)
    elseif pattern === :mixed
        for i in eachindex(w)
            w[i] = T(mod(i, 7)) / T(8)
        end
        w[1] = one(T)
    elseif pattern === :collective
        fill!(w, eps(one(T)) / T(4))
        w[1] = one(T)
    elseif pattern === :subnormal
        fill!(w, nextfloat(zero(T)))
        w[1] = one(T)
    elseif pattern === :large
        fill!(w, floatmax(T) / T(8))
        for i = 2:2:n
            w[i] /= T(8)
        end
    end
    return w
end

@testset "Compensated population means preserve count and small contributions" begin
    for T in (Float16, Float32, Float64, BigFloat)
        counts =
            T === BigFloat ? (2, 3, 127, 128, 129, 513) :
            (2, 3, 2047, 2048, 2049, 4095, 4096, 4097, 65503, 65504, 65505, 65535, 65536, 65537)
        # The minimum BigFloat subnormal has a huge default binary exponent;
        # exact-rational materialization of it is outside this hardware repair.
        patterns =
            T === BigFloat ? (:uniform, :sparse, :mixed, :collective) :
            (:uniform, :sparse, :mixed, :collective, :subnormal, :large)
        for n in counts, pattern in patterns
            weights = population_weights(T, n, pattern)
            target = exact_mean(weights)
            for order in (:forward, :reverse, :shuffled)
                input =
                    order === :forward ? weights :
                    order === :reverse ? reverse(weights) : shuffle(MersenneTwister(471), weights)
                mean = H._coulomb_normalized_mean(input, maximum(input))
                A = promote_type(T, Float64)
                @test mean isa A
                @test isfinite(mean) && mean > zero(mean)
                @test abs(Rational{BigInt}(mean) - target) <=
                      4Rational{BigInt}(eps(one(A))) * target
                pattern === :uniform && @test mean == one(A)
            end
        end
    end
end

struct SparseWeights{T} <: AbstractVector{T}
    n::Int
end
Base.size(w::SparseWeights) = (w.n,)
Base.IndexStyle(::Type{<:SparseWeights}) = IndexLinear()
function Base.getindex(w::SparseWeights{T}, i::Int) where {T}
    @boundscheck checkbounds(w, i)
    return i == 1 ? one(T) : zero(T)
end

@testset "Population mean remains wider than sparse Float16 storage" begin
    n = 2^25 + 1
    weights = SparseWeights{Float16}(n)
    mean = H._coulomb_normalized_mean(weights, Float16(1))
    @test mean === 1.0 / n
    @test mean > 0 && iszero(Float16(mean))
end

mutable struct MeanRNG{T} <: AbstractRNG
    draw::T
    events::Vector{Symbol}
end
Random.shuffle!(rng::MeanRNG, values::AbstractVector) = (push!(rng.events, :shuffle); values)
for T in (Float16, Float32, Float64, BigFloat)
    @eval Random.randn(rng::MeanRNG, ::Type{$T}) = (push!(rng.events, :normal); $T(rng.draw))
    @eval Random.rand(rng::MeanRNG, ::Type{$T}) = (push!(rng.events, :uniform); zero($T))
end

@testset "Wide sparse mean composes with scaled sampler and complete update" begin
    for T in (Float16, Float32, Float64), power in (25, 30, 40, 50, 60), strong in (false, true)
        count = (Int(1) << power) + 1
        mean = 2.0 / count
        rate = strong ? one(T) : ldexp(one(T), T === Float16 ? -24 : T === Float32 ? -80 : -500)
        variance = Rational{BigInt}(rate)^2 / Rational{BigInt}(mean)
        sample = H._coulomb_sample_increment(
            -one(T),
            zero(T),
            zero(T),
            one(T),
            rate * rate,
            one(T),
            mean,
            T(0.001),
            2T(pi),
            MeanRNG(one(T), Symbol[]),
            (rate, rate),
        )
        @test sample.isotropic == (variance > 1)
        a, b = H._coulomb_complete_update(
            (zero(T), zero(T), zero(T)),
            (one(T), zero(T), zero(T)),
            one(T),
            one(T),
            one(T),
            sample,
        )
        expected = setprecision(BigFloat, 4096) do
            if variance > 1
                ((BigFloat(1), BigFloat(0), BigFloat(0)), (BigFloat(0), BigFloat(0), BigFloat(0)))
            else
                delta = sqrt(BigFloat(variance))
                x = delta^2 / (1 + delta^2)
                z = -delta / (1 + delta^2)
                ((x, zero(x), z), (1 - x, zero(x), -z))
            end
        end
        for (actual, target) in zip((a..., b...), (expected[1]..., expected[2]...))
            @test isfinite(actual)
            @test abs(BigFloat(actual) - target) <=
                  max(BigFloat(32eps(T)) * abs(target), BigFloat(eps(T(target))) / 2)
        end
    end
end

@testset "Public collision uses the whole-population mean" begin
    for T in (Float16, Float32, Float64),
        n in (4094, 4096, 4098, 65534, 65536, 65538),
        pattern in (:uniform, :sparse),
        draw in (zero(T), T(0.5))

        ps = ParticleSet{1,T}(n)
        ps.v[1][2] = one(T)
        if pattern === :sparse
            fill!(ps.weight, zero(T))
            ps.weight[1:2] .= one(T)
        end
        gc = pattern === :uniform ? T(0.5) : T(0.001)
        dt = gc
        variance = Rational{BigInt}(gc) * Rational{BigInt}(dt) / exact_mean(ps.weight)
        @test 0 < variance < 1
        expected = setprecision(BigFloat, 256) do
            tangent = sqrt(BigFloat(variance)) * BigFloat(draw)
            longitudinal = tangent^2 / (1 + tangent^2)
            transverse = -tangent / (1 + tangent^2)
            (
                (longitudinal, zero(longitudinal), transverse),
                (1 - longitudinal, zero(longitudinal), -transverse),
            )
        end
        rng = MeanRNG(draw, Symbol[])
        before = map(copy, ps.v)
        @test collide_coulomb!(ps, gc, dt; rng, index_work = zeros(Int, n)) === ps
        @test rng.events == [:shuffle, :normal, :uniform]
        if iszero(draw)
            @test isequal(ps.v, before)
        else
            for p = 1:2, c = 1:3
                target = expected[p][c]
                bound = max(BigFloat(32eps(T)) * abs(target), BigFloat(eps(T(target))) / 2)
                @test abs(BigFloat(ps.v[c][p]) - target) <= bound
            end
            @test all(c -> all(iszero, @view ps.v[c][3:end]), 1:3)
        end
    end
end

@testset "Normalization failures leave particles, RNG and caller scratch untouched" begin
    for T in (Float16, Float32, Float64)
        ps = ParticleSet{1,T}(4)
        ps.v[1][2] = one(T)
        ps.weight .= T[floatmax(T), nextfloat(zero(T)), 0, 1]
        before = deepcopy(ps)
        work = ntuple(_ -> fill(T(-7), 4), 3)
        indices = fill(-3, 4)
        rng = MeanRNG(one(T), Symbol[])
        @test_throws ArgumentError collide_coulomb!(
            ps,
            one(T),
            one(T);
            rng,
            work,
            index_work = indices,
        )
        @test isequal(ps.v, before.v) && isequal(ps.weight, before.weight)
        @test isempty(rng.events)
        @test all(v -> all(==(T(-7)), v), work) && all(==(-3), indices)
        fill!(ps.weight, zero(T))
        @test collide_coulomb!(ps, one(T), one(T); rng, work, index_work = indices) === ps
        @test isempty(rng.events) && isequal(ps.v, before.v)
    end
end

original_precision = precision(BigFloat)
@testset "Coarse BigFloat population mean follow-up" begin
    for prec in (4, 8, 16), multiple in (1, 2, 4)
        setprecision(BigFloat, prec) do
            n = multiple * (1 << prec)
            ps = ParticleSet{1,BigFloat}(n)
            ps.v[1][2] = 1
            rng = MeanRNG(zero(BigFloat), Symbol[])
            before = map(copy, ps.v)
            mean = HybridPlasmaPIC._coulomb_normalized_mean(ps.weight, BigFloat(1))
            collide_coulomb!(ps, BigFloat(1), BigFloat(1); rng, index_work = zeros(Int, n))
            @test mean == 1
            @test rng.events == [:shuffle, :normal, :uniform]
            @test isequal(ps.v, before)
        end
    end
    @test precision(BigFloat) == original_precision
end

initial = (precision(BigFloat), rounding(BigFloat))
@testset "Local accumulator precision, independent means and rounding settings" begin
    for prec in (4, 8, 16, 24, 53, 64, 128, 256),
        mode in (RoundNearest, RoundToZero, RoundDown, RoundUp, RoundFromZero)

        setprecision(BigFloat, prec) do
            setrounding(BigFloat, mode) do
                for n in (3, 16, 1024), pattern in (:mixed, :tiny, :uniform)
                    w =
                        pattern === :mixed ? [BigFloat(mod(i, 7)) / 8 for i = 1:n] :
                        pattern === :tiny ? fill(eps(BigFloat) / 4, n) : ones(BigFloat, n)
                    w[1] = 1
                    scale = maximum(w)
                    counts = Dict{BigFloat,Int}()
                    for v in w
                        value = v / scale
                        counts[value] = get(counts, value, 0) + 1
                    end
                    target = sum(Rational{BigInt}(v) * n for (v, n) in counts) / n
                    for input in (w, reverse(w))
                        mean = HybridPlasmaPIC._coulomb_normalized_mean(input, scale)
                        @test isfinite(mean) && mean > 0 && precision(mean) == prec
                        @test abs(Rational{BigInt}(mean) - target) <=
                              4Rational{BigInt}(eps(BigFloat)) * target
                        @test (precision(BigFloat), rounding(BigFloat)) == (prec, mode)
                    end
                end
            end
        end
    end
    @test (precision(BigFloat), rounding(BigFloat)) == initial
end

@testset "Exact integer divisor survives compact MPFR count range" begin
    original = (
        ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ()),
        ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ()),
    )
    try
        @test ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), -20) == 0
        @test ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), 4) == 0
        for p in (4, 8, 64, 256)
            setprecision(BigFloat, p) do
                w = zeros(BigFloat, 128)
                w[1:2] .= 1
                @test !isfinite(BigFloat(length(w)))
                @test HybridPlasmaPIC._coulomb_normalized_mean(w, BigFloat(1)) ==
                      ldexp(BigFloat(1), -6)
            end
        end
    finally
        ccall((:mpfr_set_emin, Base.MPFR.libmpfr), Cint, (Clong,), original[1])
        ccall((:mpfr_set_emax, Base.MPFR.libmpfr), Cint, (Clong,), original[2])
    end
    @test (
        ccall((:mpfr_get_emin, Base.MPFR.libmpfr), Clong, ()),
        ccall((:mpfr_get_emax, Base.MPFR.libmpfr), Clong, ()),
    ) == original
end

end
