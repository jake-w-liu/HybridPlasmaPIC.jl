module EMPICNyquistRangeTests
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# Independent complete alternating-mean expression and nearest-even intervals.
function exact_filter(values)
    n = length(values)
    (n == 0 || isodd(n)) && return Q.(values)
    coefficient =
        sum((isodd(i) ? Q(values[i]) : -Q(values[i]) for i in eachindex(values)); init = Q(0)) / n
    return [Q(values[i]) - (isodd(i) ? coefficient : -coefficient) for i in eachindex(values)]
end
function rounded_cell(value::T, exact::Q) where {T}
    isnan(value) && return false
    maximum = Q(floatmax(T))
    threshold = maximum + (maximum - Q(prevfloat(floatmax(T)))) / 2
    isinf(value) && return signbit(value) ? exact <= -threshold : exact >= threshold
    before, after = prevfloat(value), nextfloat(value)
    lower = isfinite(before) ? (Q(before) + Q(value)) / 2 : -threshold
    upper = isfinite(after) ? (Q(after) + Q(value)) / 2 : threshold
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, value))
    return (exact > lower || (even && exact == lower)) &&
           (exact < upper || (even && exact == upper))
end
# The previous ordered arithmetic is a compatibility reference for ordinary and
# nonfinite inputs, not the exact-expression oracle for the selected cold path.
function ordered_filter!(values)
    n = length(values)
    isodd(n) && return values
    T = eltype(values)
    coefficient, sign = zero(T), one(T)
    for i in eachindex(values)
        coefficient += sign * values[i]
        sign = -sign
    end
    coefficient /= n
    sign = one(T)
    for i in eachindex(values)
        values[i] -= sign * coefficient
        sign = -sign
    end
    return values
end

@testset "Complete Nyquist subtraction retains finite range and ties" begin
    rng = MersenneTwister(80421)
    for T in (Float32, Float64)
        m, u = floatmax(T), nextfloat(zero(T))
        cases = [
            T[m, -m],
            # A tiny prefix can be lost by later normal cancellation; retaining
            # its range flag is needed for the final subnormal midpoint.
            T[u, 1, u, -1],
            T[3u, 1, 3u, -1],
            T[m, -m, -m, m],
            T[m, m, m, -m],
            T[m, m, m, m],
            T[(m/4)*3, -m/4, (m/4)*3, -m/4],
        ]
        for n in (2, 4, 6, 8, 16, 64)
            push!(cases, T[isodd(i) ? m / 2 : -m / 2 for i = 1:n])
            push!(cases, T[isodd(i) ? u : zero(T) for i = 1:n])
            push!(cases, T[isodd(i) ? 3u : zero(T) for i = 1:n])
            for _ = 1:8
                push!(cases, T.(rand(rng, -100:100, n)) .* u)
            end
        end
        for input in cases, polarity in (one(T), -one(T))
            values = polarity .* input
            expected = exact_filter(values)
            @test H._zero_nyquist_mode!(values) === values
            @test all(rounded_cell(value, exact) for (value, exact) in zip(values, expected))
        end
        # Preserve the represented nonzero DC while removing a large Nyquist mode.
        values = T[(m/4)*3, -m/4, (m/4)*3, -m/4]
        H._zero_nyquist_mode!(values)
        @test all(isfinite, values)
        @test all(x -> x > zero(T), values)
        for bits in (2, 53, 127, 1024)
            previous = precision(BigFloat)
            setprecision(BigFloat, bits) do
                values = T[m, -m, -m, m]
                expected = exact_filter(values)
                H._zero_nyquist_mode!(values)
                @test all(rounded_cell(value, exact) for (value, exact) in zip(values, expected))
                @test precision(BigFloat) == bits
            end
            @test precision(BigFloat) == previous
        end
    end
end

@testset "Nyquist ordinary, odd, empty and nonfinite compatibility" begin
    rng = MersenneTwister(62714)
    for T in (Float32, Float64), n in (0, 1, 2, 3, 4, 5, 16, 33, 128)
        values = randn(rng, T, n)
        expected = ordered_filter!(copy(values))
        @test H._zero_nyquist_mode!(values) === values
        @test isequal(values, expected)
        for bad in (T(NaN), T(Inf), -T(Inf))
            n == 0 && continue
            values = randn(rng, T, n)
            values[end] = bad
            expected = ordered_filter!(copy(values))
            @test H._zero_nyquist_mode!(values) === values
            @test isequal(values, expected)
        end
    end
    for T in (Float16, BigFloat), values in (T[], T[1, 2, 3], T[1, 3, 2, 7], T[Inf, 1])
        expected = ordered_filter!(copy(values))
        @test H._zero_nyquist_mode!(values) === values
        @test isequal(values, expected)
    end
    for T in (ComplexF32, ComplexF64),
        input in (T[], T[1+im, 2-im, 3im], T[1+2im, 3-im, 2im, -2+im], T[Inf+im, 2-im])

        values = copy(input)
        expected = ordered_filter!(copy(input))
        @test H._zero_nyquist_mode!(values) === values
        @test isequal(values, expected)
    end
end

@testset "Legacy public steps remove finite extreme transverse Nyquist current" begin
    for T in (Float32, Float64), shape in (NGP(), CIC(), TSC()), mobile in (false, true)
        n = 16
        g = FourierGrid((n,), (T(n),))
        es = EMPIC1D(g, n; shape, mobile, n0 = one(T))
        e = ParticleSet{1,T}(n; q = -one(T))
        ions = mobile ? ParticleSet{1,T}(n; q = one(T), m = es.mi) : nothing
        for ps in (e, ions)
            ps === nothing && continue
            ps.x[1] .= T.(0:n-1) .+ T(0.25)
            fill!(ps.weight, mobile ? T(0.5) : one(T))
            foreach(v -> fill!(v, zero(T)), ps.v)
            ps.v[2] .= T[isodd(i) ? floatmax(T) / 2 : -floatmax(T) / 2 for i = 1:n]
            ps === ions && (ps.v[2] .*= -one(T))
        end
        init_empic!(es, e, ions)
        storage =
            Tuple(a for ps in (e, ions) if ps !== nothing for a in (ps.x..., ps.v..., ps.weight))
        before = map(copy, storage)
        for (dt, ns) in ((T(0.01), 1), (T(0.02), 2))
            es.n_sub = ns
            @test step_empic!(es, e, ions, dt) === es
            @test all(iszero, es.Jy)
            @test all(iszero, es.Ey)
            @test all(iszero, es.Bz)
            @test all(isequal(a, b) for (a, b) in zip(storage, before))
        end
    end
end
end
