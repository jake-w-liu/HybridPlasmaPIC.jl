# Phase-4 forcing operators: Landau-Lifshitz radiation reaction (synchrotron cooling) and
# the external antenna/RF field source (divergence-free Faraday injection). Uniform applied
# fields + gravity + a uniform RF drive are already provided by push_uniform!.

using HybridPlasmaPIC, Test, Random
using HybridPlasmaPIC: FourierGrid, ParticleSet, divergence!, curl!

include("test_radiation_inputs.jl")
include("test_radiation_arithmetic.jl")

@testset "RR-001 radiation reaction: synchrotron cooling rate vs analytic" begin
    T = Float64
    γof(vx) = 1 / sqrt(1 - vx^2)
    K = 0.001
    B0 = 1.0
    dt = 0.05
    ps = ParticleSet{1,T}(1)
    ps.v[1][1] = 0.995                                  # ultrarelativistic, ⊥ B
    γ0 = γof(ps.v[1][1])
    inv = Float64[1/γ0]
    for _ = 1:60
        apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, B0), dt; K = K)
        push!(inv, 1 / γof(ps.v[1][1]))
    end
    @test γ0 > 9                                        # started ultrarelativistic
    @test all(diff(inv) .>= 0)                          # 1/γ increases ⇒ monotone cooling
    # early-time slope d(1/γ)/dt ≈ K·B²·v0⁴ (Landau-Lifshitz synchrotron)
    slope = (inv[21] - inv[1]) / (20 * dt)
    @test isapprox(slope, K * B0^2 * 0.995^4; rtol = 0.02)
end

@testset "RR-002 radiation reaction: no-op, guards, determinism" begin
    T = Float64
    ps = ParticleSet{1,T}(10)
    for i = 1:10
        ps.v[1][i] = 0.5
    end
    snap = copy(ps.v[1])
    apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.1; K = 0.0)   # K=0
    @test ps.v[1] == snap
    apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.0; K = 1.0)   # dt=0
    @test ps.v[1] == snap
    @test_throws ArgumentError apply_radiation_reaction!(
        ps,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = -1.0,
    )
    @test_throws ArgumentError apply_radiation_reaction!(
        ps,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = 1.0,
        c = 0.0,
    )
    a = ParticleSet{1,T}(5)
    b = ParticleSet{1,T}(5)
    for i = 1:5
        a.v[1][i] = 0.8
        b.v[1][i] = 0.8
    end
    apply_radiation_reaction!(a, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.05; K = 0.01)
    apply_radiation_reaction!(b, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.05; K = 0.01)
    @test a.v[1] == b.v[1]
    # exponential integrator is unconditionally stable: a large K·dt cools, never reverses/heats
    ps3 = ParticleSet{1,T}(1)
    ps3.v[1][1] = 0.99
    γ0 = 1 / sqrt(1 - 0.99^2)
    apply_radiation_reaction!(ps3, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 0.5; K = 1.0)   # overshoot regime
    @test ps3.v[1][1] > 0 && 1 / sqrt(1 - ps3.v[1][1]^2) < γ0
    @test_throws ArgumentError apply_radiation_reaction!(
        ps3,
        (NaN, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = 0.01,
    )

    bad = ParticleSet{1,T}(2)
    bad.v[1] .= (0.8, 1.0)
    bad0 = map(copy, bad.v)
    @test_throws ArgumentError apply_radiation_reaction!(
        bad,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = 0.01,
    )
    @test all(bad.v[c] == bad0[c] for c = 1:3)

    bad.v[1][2] = NaN
    bad0 = map(copy, bad.v)
    @test_throws ArgumentError apply_radiation_reaction!(
        bad,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        0.1;
        K = 0.01,
    )
    @test all(isequal(bad.v[c], bad0[c]) for c = 1:3)
end

@testset "ANT-001 antenna injects −dt·∇×E_ant, divergence-free" begin
    T = Float64
    N = (16, 12)
    L = (2π, 2π)
    g = FourierGrid(N, L)
    xs = [(i - 1) * L[1] / N[1] for i = 1:N[1]]
    ys = [(j - 1) * L[2] / N[2] for j = 1:N[2]]
    # a smooth structured antenna E-field (all components), exercises every curl term
    Eant = (
        [sin(2xs[i]) * cos(ys[j]) for i = 1:N[1], j = 1:N[2]],
        [cos(xs[i]) * sin(3ys[j]) for i = 1:N[1], j = 1:N[2]],
        [sin(xs[i] + 2ys[j]) for i = 1:N[1], j = 1:N[2]],
    )
    B = ntuple(_ -> zeros(T, N), 3)
    dt = 0.1
    apply_antenna!(B, Eant, dt, g)
    # matches −dt·∇×E_ant
    ce = ntuple(_ -> zeros(T, N), 3)
    curl!(ce, Eant, g)
    for c = 1:3
        @test maximum(abs.(B[c] .- (-dt .* ce[c]))) < 1e-12
    end
    # divergence-free: ∇·B ≈ 0 (a curl has zero divergence)
    dv = zeros(T, N)
    divergence!(dv, B, g)
    @test maximum(abs.(dv)) < 1e-12
end

@testset "ANT-002 antenna: no-op (dt=0) & validation" begin
    T = Float64
    N = (8, 8)
    g = FourierGrid(N, (2π, 2π))
    Eant = ntuple(_ -> [sin(i + j) for i = 1:N[1], j = 1:N[2]] .* 1.0, 3)
    B = ntuple(_ -> ones(T, N), 3)
    snap = map(copy, B)
    apply_antenna!(B, Eant, 0.0, g)                     # dt=0 ⇒ untouched
    @test all(B[c] == snap[c] for c = 1:3)
    @test_throws ArgumentError apply_antenna!(B, Eant, -0.1, g)
    Bbad = ntuple(_ -> zeros(T, (N[1] + 1, N[2])), 3)   # B not matching the grid ⇒ guarded
    @test_throws DimensionMismatch apply_antenna!(Bbad, Eant, 0.1, g)
end

@testset "ANT-003 antenna preserves old field and grid storage" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2π),))
        @test g.kvec[1][2] == one(T)
        wave, slope = T[0, 1, 0, -1], T[1, 0, -1, 0]
        Eant = (zeros(T, 4), copy(wave), copy(wave))
        curlE = (zeros(T, 4), -slope, slope)
        B = ntuple(c -> fill(T(c), 4), 3)
        work = ntuple(_ -> fill(T(7), 4), 3)
        dt = T(0.125)
        expected = ntuple(c -> B[c] .- dt .* curlE[c], 3)
        @test apply_antenna!(B, Eant, dt, g; work) === B
        @test all(isapprox(B[c], expected[c]; rtol = 8eps(T), atol = 8eps(T)) for c = 1:3)
        # Input sharing is safe: the curl finishes before the field update.
        shared = map(copy, Eant)
        expected = ntuple(c -> shared[c] .- dt .* curlE[c], 3)
        apply_antenna!(shared, shared, dt, g; work)
        @test all(isapprox(shared[c], expected[c]; rtol = 8eps(T), atol = 8eps(T)) for c = 1:3)
        for (badB, badwork) in (
            (B, B),
            (B, (work[1], B[3], work[3])),
            ((B[1], B[1], B[3]), work),
            ((g.kfull[1], B[2], B[3]), work),
            (B, (work[1], g.kvec[1], work[3])),
        )
            storage = (
                badB...,
                badwork...,
                Eant...,
                g.cbuf,
                g.tbuf,
                g.abuf,
                g.sbuf,
                g.midx...,
                g.kfull...,
                g.ik...,
                g.kvec...,
            )
            before = copy.(storage)
            @test_throws ArgumentError apply_antenna!(badB, Eant, dt, g; work = badwork)
            @test isequal(storage, before)
        end
    end
end

@testset "ANT-004 finite stored-curl range and staged publication" begin
    H = HybridPlasmaPIC
    Q = Rational{BigInt}
    oldmode = get_zero_subnormals()
    try
        @test set_zero_subnormals(false)
        for T in (Float32, Float64)
            U = T === Float32 ? UInt32 : UInt64
            @test rounding(T) === RoundNearest
            g = FourierGrid((4,), (T(2π),))
            @test g.kvec[1][2] == one(T)
            amplitude = ldexp(T(3), exponent(floatmax(T)) - 2)
            Eant = (zeros(T, 4), amplitude .* T[0, 1, 0, -1], zeros(T, 4))
            electric_before = map(copy, Eant)
            samplecurl = ntuple(_ -> zeros(T, 4), 3)
            curl!(samplecurl, Eant, g)
            maximum_exact = Q(floatmax(T))
            observation = (
                type = T,
                julia = VERSION,
                spectral_operators = Base.pkgversion(H.SpectralOperators),
                flush = get_zero_subnormals(),
                electric_words = collect(reinterpret(U, Eant[2])),
                curl_words = collect(reinterpret(U, samplecurl[3])),
                native_product_words = collect(reinterpret(U, T(3) .* samplecurl[3])),
            )
            @info "Antenna finite-curl fixture" observation
            # Qualify actual FFT storage; the oracle below concerns this stored curl.
            @test all(a -> all(isfinite, a), samplecurl)
            @test all(a -> all(x -> abs(Q(x)) <= maximum_exact / 2, a), samplecurl)
            @test any(x -> abs(Q(x)) > maximum_exact / 3, samplecurl[3])
            @test any(x -> !isfinite(T(3) * x), samplecurl[3])
            expected = map(a -> -T(2) .* a, samplecurl)
            @test all(a -> all(isfinite, a), expected)
            for supplied_work in (false, true)
                B = map(copy, samplecurl)
                work = ntuple(_ -> fill(T(7), 4), 3)
                returned =
                    supplied_work ? apply_antenna!(B, Eant, T(3), g; work) :
                    apply_antenna!(B, Eant, T(3), g)
                @test returned === B
                # B=c and dt=3 give -2c, an exact finite binary scaling.
                @test B == expected
                @test isequal(Eant, electric_before)
            end

            # The late active component overflows genuinely; earlier candidates stay private.
            B = map(a -> -a, samplecurl)
            before = map(a -> collect(reinterpret(U, a)), B)
            work = ntuple(_ -> fill(T(7), 4), 3)
            @test_throws ArgumentError apply_antenna!(B, Eant, T(3), g; work)
            @test all(collect(reinterpret(U, B[c])) == before[c] for c = 1:3)
            @test isequal(Eant, electric_before)
            B = ntuple(_ -> zeros(T, 4), 3)
            B[3][end] = T(NaN)
            before = map(a -> collect(reinterpret(U, a)), B)
            @test_throws ArgumentError apply_antenna!(B, Eant, T(0.125), g; work)
            @test all(collect(reinterpret(U, B[c])) == before[c] for c = 1:3)

            # Finite E can have a nonfinite materialized curl despite a small dt.
            short = FourierGrid((4,), (T(2π) / T(4),))
            @test short.kvec[1][2] == T(4)
            curl!(work, Eant, short)
            observation = (type = T, curl_words = collect(reinterpret(U, work[3])))
            @info "Antenna materialized-curl boundary" observation
            @test any(x -> !isfinite(x), work[3])
            B = ntuple(_ -> fill(one(T), 4), 3)
            before = map(a -> collect(reinterpret(U, a)), B)
            @test_throws ArgumentError apply_antenna!(B, Eant, T(0.0625), short; work)
            @test all(collect(reinterpret(U, B[c])) == before[c] for c = 1:3)
            @test isequal(Eant, electric_before)

            # Ordinary arithmetic keeps the historical native two-operation bits.
            ordinary_E = (zeros(T, 4), T[0, 1, 0, -1], T[0, -2, 0, 2])
            curl!(work, ordinary_E, g)
            B = ntuple(c -> fill(T(c), 4), 3)
            dt = T(0.125)
            native = ntuple(c -> T[B[c][i] - dt * work[c][i] for i in eachindex(B[c])], 3)
            apply_antenna!(B, ordinary_E, dt, g; work)
            @test isequal(B, native)
            badE = map(copy, ordinary_E)
            badE[3][end] = T(Inf)
            before = map(a -> collect(reinterpret(U, a)), B)
            @test_throws ArgumentError apply_antenna!(B, badE, dt, g; work)
            @test all(collect(reinterpret(U, B[c])) == before[c] for c = 1:3)
            for dtzero in (zero(T), -zero(T))
                @test apply_antenna!(B, badE, dtzero, g; work) === B
                @test all(collect(reinterpret(U, B[c])) == before[c] for c = 1:3)
            end

            # Unsupported arithmetic rejects before any FFT or workspace write.
            storage = (B..., ordinary_E..., work..., g.cbuf, g.tbuf, g.abuf, g.sbuf)
            stored_words = map(a -> collect(reinterpret(U, a)), storage)
            @test set_zero_subnormals(true)
            @test get_zero_subnormals()
            observation = (type = T, requested_flush = true, active_flush = get_zero_subnormals())
            @info "Antenna context admission" observation
            for dtzero in (zero(T), -zero(T))
                @test apply_antenna!(B, ordinary_E, dtzero, g; work) === B
            end
            @test_throws ArgumentError apply_antenna!(B, ordinary_E, dt, g; work)
            @test all(collect(reinterpret(U, a)) == b for (a, b) in zip(storage, stored_words))
            @test set_zero_subnormals(false)
        end

        # U=2^971. Stored |c|=M/5+U/20, so 5|c|=M+U/4 rounds to M.
        # Native M+3U/8 also rounds M, but exact M+5U/8 exceeds M+U/2.
        T = Float64
        b, c, dt = ldexp(T(3), 968), -reinterpret(T, 0x7fc9999999999999), T(5)
        product = dt * c
        native = b - product
        exact = Q(b) - Q(dt) * Q(c)
        threshold = Q(floatmax(T)) + Q(BigInt(1) << 971) / 2
        observation = (
            b_word = reinterpret(UInt64, b),
            curl_word = reinterpret(UInt64, c),
            native_product = product,
            native_output = native,
            exact_overflow = exact > threshold,
        )
        @info "Antenna finite-native overflow boundary" observation
        @test product === -floatmax(T)
        @test native === floatmax(T)
        @test exact > threshold
        @test_throws ArgumentError H._antenna_update_value(b, c, dt)
    finally
        set_zero_subnormals(oldmode)
    end
    @test get_zero_subnormals() == oldmode
end
