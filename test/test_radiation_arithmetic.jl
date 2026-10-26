using HybridPlasmaPIC, Test, Random

# Independent direct momentum expression at 8192 bits, evaluated serially (MPFR's shared
# constant caches are unsafe under concurrent use, NUM-001). It damps p = γv and recovers v;
# the production velocity factor and its fallbacks are not used.
function radiation_momentum_oracle(v, E, B, dt, K, c)
    return setprecision(BigFloat, 8192) do
        setrounding(BigFloat, RoundNearest) do
            u, e, b = map(BigFloat, v), map(BigFloat, E), map(BigFloat, B)
            cwide = BigFloat(c)
            γ = inv(sqrt(1 - sum(abs2, u) / cwide^2))
            w = (
                e[1] + u[2] * b[3] - u[3] * b[2],
                e[2] + u[3] * b[1] - u[1] * b[3],
                e[3] + u[1] * b[2] - u[2] * b[1],
            )
            χ2 = sum(abs2, w) - sum(u .* e)^2 / cwide^2
            momentum = map(x -> γ * x * exp(-BigFloat(K) * BigFloat(dt) * γ * χ2), u)
            map(x -> x / sqrt(1 + sum(abs2, momentum) / cwide^2), momentum)
        end
    end
end

# The update before this repair (momentum round trip in particle precision), for showing
# that the witnesses below exercise RR-002/RR-003 and for BigFloat identity.
function legacy_radiation_velocity(v::NTuple{3,T}, E, B, dt, K, c) where {T}
    ux, uy, uz = v
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    c2 = c * c
    γ = one(T) / sqrt(one(T) - (ux * ux + uy * uy + uz * uz) / c2)
    wx = Ex + (uy * Bz - uz * By)
    wy = Ey + (uz * Bx - ux * Bz)
    wz = Ez + (ux * By - uy * Bx)
    vE = ux * Ex + uy * Ey + uz * Ez
    χ2 = (wx * wx + wy * wy + wz * wz) - (vE * vE) / c2
    χ2 > 0 || return v
    s = γ * exp(-K * γ * χ2 * dt)
    px, py, pz = s * ux, s * uy, s * uz
    γn = sqrt(one(T) + (px * px + py * py + pz * pz) / c2)
    return (px / γn, py / γn, pz / γn)
end

# Largest componentwise error in ulps of the stored speed |v| (normwise contract).
function radiation_error_ulps(actual, expected, v::NTuple{3,T}) where {T}
    speed = setprecision(() -> sqrt(sum(x -> BigFloat(x)^2, v)), BigFloat, 256)
    ulp = BigFloat(eps(T(speed)))
    return maximum(abs(BigFloat(a) - b) for (a, b) in zip(actual, expected)) / ulp
end

# The radiation preflight accepts a stored velocity when it is exactly subluminal and its
# rounded |v/c|² stays below one.
radiation_accepted(v::NTuple{3,T}, c) where {T} =
    sum(x -> BigFloat(x)^2, v) < BigFloat(c)^2 && sum(x -> (x / c)^2, v) < one(T)

function radiation_apply(v::NTuple{3,T}, E, B, dt, K, c; D = 1) where {T}
    ps = ParticleSet{D,T}(1)
    foreach(x -> fill!(x, T(0.2)), ps.x)
    for d = 1:3
        ps.v[d][1] = v[d]
    end
    before = deepcopy((ps.x, ps.weight, ps.id, ps.tag, ps.q, ps.m))
    same = apply_radiation_reaction!(ps, E, B, dt; K, c) === ps
    unchanged = (ps.x, ps.weight, ps.id, ps.tag, ps.q, ps.m) == before
    return ntuple(d -> ps.v[d][1], 3), same && unchanged
end

# Normwise accuracy, finiteness, no componentwise speed increase, preserved signs and other
# particle state, for D = 1..3. Returns the D = 1 result.
function check_radiation(v::NTuple{3,T}, E, B, dt, K, c, expected; ulps = 8) where {T}
    results = ntuple(D -> radiation_apply(v, E, B, dt, K, c; D), 3)
    actual = results[1][1]
    @test all(r -> r[2] && r[1] === actual, results)
    @test all(isfinite, actual)
    err = radiation_error_ulps(actual, expected, v)
    @test err <= ulps || (println("radiation error $err ulps: ", (v, E, B, dt, K, c)); false)
    @test all(abs(actual[d]) <= abs(v[d]) for d = 1:3) &&
          all(iszero(actual[d]) || signbit(actual[d]) == signbit(v[d]) for d = 1:3)
    return actual
end
check_radiation(v, E, B, dt, K, c; kwargs...) =
    check_radiation(v, E, B, dt, K, c, radiation_momentum_oracle(v, E, B, dt, K, c); kwargs...)

# Exact symmetries of the update: a cyclic axis rotation, and v, E → -v, -E.
radiation_rotate(t, axis) = ntuple(d -> t[mod1(d - axis + 1, 3)], 3)

normalize3(x) = x ./ sqrt(sum(abs2, x))

# γ and χ² of the stored inputs (256 bits), used only to choose sensitive optical depths.
function radiation_invariants(v, E, B, c)
    return setprecision(BigFloat, 256) do
        u, e, b = map(BigFloat, v), map(BigFloat, E), map(BigFloat, B)
        cwide = BigFloat(c)
        w = (
            e[1] + u[2] * b[3] - u[3] * b[2],
            e[2] + u[3] * b[1] - u[1] * b[3],
            e[3] + u[1] * b[2] - u[2] * b[1],
        )
        inv(sqrt(1 - sum(abs2, u) / cwide^2)), sum(abs2, w) - sum(u .* e)^2 / cwide^2
    end
end

radiation_allocated(ps, E, B, dt, K, c) = @allocated apply_radiation_reaction!(ps, E, B, dt; K, c)

radiation_extremes(::Type{Float16}) = (Float16(2)^-5, Float16(2)^5)
radiation_extremes(::Type{Float32}) = (2.0f0^-40, 2.0f0^40)
radiation_extremes(::Type{Float64}) = (2.0^-300, 2.0^300)

legacy_radiation_wrong(v, E, B, dt, K, c, expected) =
    let legacy = legacy_radiation_velocity(v, E, B, dt, K, c)
        !all(isfinite, legacy) || radiation_error_ulps(legacy, expected, v) > 8
    end

@testset "Radiation extreme-scale matrix matches the momentum reference normwise" begin
    for T in (Float16, Float32, Float64)
        s, L = radiation_extremes(T)
        z, o = zero(T), one(T)
        legacy_wrong = 0
        for c in (s, o, L), Em in (z, s, o, L), Bm in (z, s, o, L)
            (iszero(Em) && iszero(Bm)) && continue
            # depth scales K·dt to the force: α ≈ depth·γ·χ²/scale² (all powers of two).
            for parallel in (true, false), depth in (T(2)^-(precision(T) + 2), o, T(16))
                scale = max(Em, c * Bm / 2)
                dt, K = inv(scale), depth / scale
                v0 = (c / 2, z, z)
                E0 = parallel ? (Em, z, z) : (z, z, -Em)   # along v, or against v×B
                B0 = (z, Bm, z)
                reference = radiation_momentum_oracle(v0, E0, B0, dt, K, c)
                for axis = 1:3, sign in (o, -o)
                    v = map(x -> sign * x, radiation_rotate(v0, axis))
                    E = map(x -> sign * x, radiation_rotate(E0, axis))
                    B = radiation_rotate(B0, axis)
                    expected = map(x -> sign * x, radiation_rotate(reference, axis))
                    check_radiation(v, E, B, dt, K, c, expected)
                    legacy_wrong += legacy_radiation_wrong(v, E, B, dt, K, c, expected)
                end
            end
        end
        # K and dt at independent extremes, mostly saturating to complete or no damping.
        for c in (s, o, L), (Em, Bm) in ((s, z), (o, o), (L, L), (z, L), (L, s))
            for (K, dt) in ((s, s), (s, o), (o, o), (L, o), (L, L), (L, s))
                v0, E0, B0 = (c / 2, z, z), (Em, z, z), (z, Bm, z)
                reference = radiation_momentum_oracle(v0, E0, B0, dt, K, c)
                for axis = 1:3, sign in (o, -o)
                    v = map(x -> sign * x, radiation_rotate(v0, axis))
                    E = map(x -> sign * x, radiation_rotate(E0, axis))
                    B = radiation_rotate(B0, axis)
                    expected = map(x -> sign * x, radiation_rotate(reference, axis))
                    check_radiation(v, E, B, dt, K, c, expected)
                    legacy_wrong += legacy_radiation_wrong(v, E, B, dt, K, c, expected)
                end
            end
        end
        # The matrix contains cases the previous arithmetic got wrong (RR-002).
        @test legacy_wrong > 0
    end
end

@testset "Radiation near the light cone and along E keeps the frozen-rate result" begin
    for T in (Float16, Float32, Float64)
        z, o = zero(T), one(T)
        s, L = radiation_extremes(T)
        for c in (s, o, L)
            # Largest accepted speed on an axis, and a planar velocity nearly on the cone.
            fast = (prevfloat(c), z, z)
            a = T(0.6) * c
            b = T(sqrt(max(0.0, Float64(c)^2 - Float64(a)^2)))
            while !radiation_accepted((a, b, z), c)
                b = prevfloat(b)
            end
            for v in (fast, (a, b, z))
                speed = sqrt(sum(x -> BigFloat(x)^2, v))
                n = map(x -> x / speed, map(BigFloat, v))
                γ, _ = radiation_invariants(v, (z, z, z), (z, z, z), c)
                # Synchrotron geometry (v ⊥ B), and a field along v plus a normal part of
                # relative size 1/γ: |w|² and (v·E/c)² then cancel to χ² ≈ 2/γ².
                fields = (((z, z, z), (z, z, T(4) / c)), ((T(n[1]), T(n[2]), T(inv(γ))), (z, z, z)))
                for (E, B) in fields
                    _, χ2 = radiation_invariants(v, E, B, c)
                    # Optical depths around ln γ, where the factor is most sensitive.
                    for m in (0.5, 1.0, 2.0)
                        depth = BigFloat(m) * max(log(γ), 1) / (γ * χ2)
                        K = T(depth * BigFloat(c)^2)
                        dt = inv(T(c))^2
                        (isfinite(K) && K > 0 && dt > 0) || continue
                        check_radiation(v, E, B, dt, K, c)
                    end
                end
            end
        end
    end
end

# A stored velocity with exact 0 < 1 - |v|² ≤ 2^-k: two components put the remaining gap
# near 2^-(p+1), and a small third component then tunes it with granularity ≈ 2^-2p.
function radiation_near_null(rng, ::Type{T}, k) where {T}
    p = precision(T)
    for _ = 1:200_000
        a = T(0.6 + 0.1rand(rng))
        b = T(sqrt(1 - Float64(a)^2 - 2.0^-(p + 1) * (1 + rand(rng))))
        gap = setprecision(() -> 1 - BigFloat(a)^2 - BigFloat(b)^2, BigFloat, 512)
        gap > 0 || continue
        c = T(sqrt(max(gap - BigFloat(2)^-k, zero(gap))))
        r = setprecision(() -> gap - BigFloat(c)^2, BigFloat, 512)
        v = (a, b, c)
        (0 < r <= BigFloat(2)^-k && radiation_accepted(v, one(T))) && return v
    end
    error("no near-null velocity found for $T, k = $k")
end

@testset "Radiation resolves three-component speeds near the light cone" begin
    rng = Xoshiro(5)
    for T in (Float32, Float64), k in (T === Float32 ? (30, 40) : (60, 90, 110))
        z, o = zero(T), one(T)
        v = radiation_near_null(rng, T, k)
        γ, _ = radiation_invariants(v, (z, z, z), (z, z, z), o)
        @test γ >= BigFloat(2)^(k ÷ 2)
        for (E, B) in (((z, z, z), (z, z, o)), ((v[1], v[2], v[3] + T(0.01)), (z, z, z)))
            _, χ2 = radiation_invariants(v, E, B, o)
            for m in (0.5, 2.0)
                check_radiation(v, E, B, o, T(m * max(log(γ), 1) / (γ * χ2)), o)
            end
        end
    end
end

@testset "Radiation random drag matches the momentum reference normwise" begin
    rng = Xoshiro(20261026)
    for T in (Float16, Float32, Float64)
        o = one(T)
        draws = 0
        while draws < 240
            kind = draws % 3
            if kind == 0       # ordinary speeds, fields and strong drag
                v = ntuple(_ -> T(rand(rng) - 0.5), 3)
                E = ntuple(_ -> T(randn(rng)), 3)
                B = ntuple(_ -> T(randn(rng)), 3)
                K, dt = T(rand(rng)), T(rand(rng))
            elseif kind == 1   # relativistic, 1 - β² down to near eps, sensitive depth
                r = exp(log(T === Float16 ? 1e-3 : T === Float32 ? 1e-6 : 1e-14) * rand(rng))
                n = normalize3(randn(rng, 3))
                v = ntuple(i -> T(sqrt(1 - r) * n[i]), 3)
                E = ntuple(_ -> T(randn(rng)), 3)
                B = ntuple(_ -> T(randn(rng)), 3)
                K, dt = T(exp(log(100) * rand(rng)) * 0.1 * sqrt(r)), o
            else               # near an E×B drift: inexact cancellation times large K·dt
                v = ntuple(_ -> T(rand(rng) - 0.5), 3)
                B = ntuple(_ -> T(randn(rng)), 3)
                drift = (
                    v[2] * B[3] - v[3] * B[2],
                    v[3] * B[1] - v[1] * B[3],
                    v[1] * B[2] - v[2] * B[1],
                )
                δ = exp(log(T === Float16 ? 1e-2 : 1e-6) * rand(rng))
                E = ntuple(i -> T(-drift[i] + δ * randn(rng)), 3)
                K, dt = T(min(rand(rng) / δ^2, floatmax(T) / 4)), o
            end
            radiation_accepted(v, o) || continue
            check_radiation(v, E, B, dt, K, o)
            draws += 1
        end
    end
end

@testset "Radiation weak drag never increases a stored speed" begin
    rng = Xoshiro(534)
    for T in (Float16, Float32, Float64)
        legacy_heats = 0
        for _ = 1:1000
            v = ntuple(_ -> T(0.8rand(rng) - 0.4), 3)
            E = ntuple(_ -> T(rand(rng)), 3)
            B = ntuple(_ -> T(rand(rng)), 3)
            K = eps(T)^2 * T(rand(rng))
            actual, intact = radiation_apply(v, E, B, one(T), K, one(T))
            @test intact
            @test all(abs(actual[d]) <= abs(v[d]) for d = 1:3)
            @test sum(x -> BigFloat(x)^2, actual) <= sum(x -> BigFloat(x)^2, v)
            legacy = legacy_radiation_velocity(v, E, B, one(T), K, one(T))
            legacy_heats += sum(x -> BigFloat(x)^2, legacy) > sum(x -> BigFloat(x)^2, v)
        end
        @test legacy_heats > 0   # RR-003 witness for the previous arithmetic
    end
end

@testset "Radiation exact zero force returns the stored velocity" begin
    for T in (Float16, Float32, Float64)
        z, o = zero(T), one(T)
        s, L = radiation_extremes(T)
        for c in (s, o, L), field in (s, o, L), K in (s, o, L), D = 1:3
            v = (c / 2, -c / 4, z)
            # No fields; v ∥ B; and E = -v×B with exactly representable products.
            for (E, B) in (
                ((z, z, z), (z, z, z)),
                ((z, z, z), (-field, field / 2, z)),
                ((z, z, -c * field / 2), (z, field, z)),
            )
                ps = ParticleSet{D,T}(2)
                for d = 1:3
                    ps.v[d] .= (v[d], z)
                end
                apply_radiation_reaction!(ps, E, B, o; K, c)
                @test all(ps.v[d][1] === v[d] for d = 1:3)
                @test all(ps.v[d][2] === z for d = 1:3)
            end
        end
    end
end

@testset "Radiation inexact cancellation amplified by K·dt uses the complete force" begin
    for (T, field, coefficient) in ((Float32, 1.0f-25, 1.0f37), (Float64, 1e-170, 1e300))
        v = (T(0.3), T(0.3), zero(T))
        force = v[1] * field
        E, B = (-force, force, zero(T)), (zero(T), zero(T), field)
        residual = Rational{BigInt}(v[1]) * Rational{BigInt}(field) - Rational{BigInt}(force)
        @test !iszero(residual)
        @test 2residual^2 * Rational{BigInt}(coefficient)^2 > 10000
        actual = check_radiation(v, E, B, coefficient, coefficient, one(T))
        @test all(iszero, actual)
        # One-ulp perturbation of an exact zero force with a huge coefficient.
        for field in (T(2)^-20, one(T), T(2)^20)
            E = (zero(T), nextfloat(field / 2), zero(T))
            check_radiation(
                (T(0.5), zero(T), zero(T)),
                E,
                (zero(T), zero(T), field),
                one(T),
                floatmax(T) / 2,
                one(T),
            )
        end
    end
end

@testset "Radiation BigFloat particles keep the ambient-precision momentum update" begin
    stored = setprecision(BigFloat, 512) do
        (sqrt(BigFloat(2)) / 4, BigFloat(1) / 3, -BigFloat(1) / 5)
    end
    for bits in (24, 53, 256), mode in (RoundNearest, RoundDown, RoundUp)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                o, z = one(BigFloat), zero(BigFloat)
                for (v, E, B, dt, K) in (
                    (stored, (o, o / 2, z), (z, o, o / 4), o / 10, o / 3),
                    (
                        (o / 2, o / 4, z),
                        (BigFloat("1e-1000"), z, z),
                        (z, BigFloat("1e-1000"), z),
                        BigFloat("1e1000"),
                        BigFloat("1e1000"),
                    ),
                )
                    ps = ParticleSet{1,BigFloat}(1)
                    for d = 1:3
                        ps.v[d][1] = v[d]
                    end
                    apply_radiation_reaction!(ps, E, B, dt; K)
                    actual = ntuple(d -> ps.v[d][1], 3)
                    # Stored velocities keep their own precision; arithmetic is ambient.
                    expected = legacy_radiation_velocity(v, E, B, dt, K, o)
                    @test all(actual[d] == expected[d] for d = 1:3)
                    @test all(precision(x) == bits for x in actual)
                end
                @test precision(BigFloat) == bits
                @test rounding(BigFloat) == mode
            end
        end
    end
end

@testset "Radiation frozen-rate substep converges to synchrotron evolution" begin
    c, B, K, duration, initial = 2.0, 0.7, 0.08, 2.0, 1.6
    expected = setprecision(BigFloat, 512) do
        # y = 1/γ obeys y' = K*c²*B²*(1-y²)². Integrate its reciprocal
        # analytically and invert the strictly increasing primitive by bisection.
        primitive(y) = y / (2 * (1 - y^2)) + atanh(y) / 2
        y0 = sqrt(1 - (BigFloat(initial) / c)^2)
        target = primitive(y0) + BigFloat(K) * BigFloat(c)^2 * BigFloat(B)^2 * duration
        lo, hi = y0, one(y0)
        for _ = 1:400
            mid = (lo + hi) / 2
            if primitive(mid) < target
                lo = mid
            else
                hi = mid
            end
        end
        Float64(BigFloat(c) * sqrt(1 - ((lo + hi) / 2)^2))
    end
    errors = map((4, 8, 16, 32)) do steps
        ps = ParticleSet{1,Float64}(1)
        ps.x[1][1] = 0
        ps.v[1][1] = initial
        for _ = 1:steps
            apply_radiation_reaction!(ps, (0.0, 0.0, 0.0), (0.0, 0.0, B), duration / steps; K, c)
        end
        abs(ps.v[1][1] - expected)
    end
    @test all(>(0), errors)
    @test all(1.8 < errors[i] / errors[i+1] < 2.2 for i = 1:3)
end

@testset "Radiation evaluation is allocation-free and task-independent" begin
    rng = Xoshiro(77)
    for T in (Float16, Float32, Float64), D = 1:3
        s, L = radiation_extremes(T)
        cases = (
            (T.((0.3, 0.2, -0.1)), T.((0.1, -0.2, 0.4)), T(0.05), T(0.01), one(T), one(T)),  # ordinary
            (T.((0.9, -0.6, 0.3)), T.((0.5, 1.1, -0.7)), T(0.5), T(8), one(T), one(T)),     # strong
            ((L, zero(T), zero(T)), (zero(T), zero(T), s), s, s, L, L),                   # out of range
        )
        for (E, B, dt, K, c, speed) in cases
            ps = ParticleSet{D,T}(64)
            for d = 1:3
                ps.v[d] .= T.(speed * 0.55 * (rand(rng, 64) .- 0.5))
            end
            initial = map(copy, ps.v)
            apply_radiation_reaction!(ps, E, B, dt; K, c)
            serial = map(copy, ps.v)
            foreach(d -> copyto!(ps.v[d], initial[d]), 1:3)
            @test radiation_allocated(ps, E, B, dt, K, c) == 0
            @test ps.v == serial
            copies = [deepcopy(ps) for _ = 1:4]
            foreach(q -> foreach(d -> copyto!(q.v[d], initial[d]), 1:3), copies)
            tasks = [Threads.@spawn apply_radiation_reaction!(q, E, B, dt; K, c) for q in copies]
            foreach(wait, tasks)
            @test all(q.v == serial for q in copies)
        end
    end
end

@testset "Radiation retains small amplified forces beside cancelling large fields" begin
    for exponent in (250, 260, 270, 300, 520, 600, 900),
        axis = 1:3,
        sign in (-1.0, 1.0),
        transverse in (false, true)

        small, large = ldexp(1.0, -exponent), ldexp(1.0, exponent)
        velocity = (sign * 0.5, 0.0, 0.0)
        electric = transverse ? (0.0, sign * large / 2, small) : (small, sign * large / 2, 0.0)
        magnetic = (0.0, 0.0, large)
        rotate(t) = ntuple(d -> t[mod1(d + axis - 1, 3)], 3)
        v, E, B = rotate(velocity), rotate(electric), rotate(magnetic)
        expected = radiation_momentum_oracle(v, E, B, large, large, 1.0)
        actual, unchanged = radiation_apply(v, E, B, large, large, 1.0)
        @test unchanged
        @test radiation_error_ulps(actual, expected, v) <= 8
        @test all(abs(actual[d]) <= abs(v[d]) for d = 1:3)
    end
    for exponent in (300, 600, 900), cpower in (-100, 0, 100), coefficient_shift in (-100, 0, 20)
        c = ldexp(1.0, cpower)
        v = (c / 2, 0.0, 0.0)
        E = (ldexp(1.0, -exponent), ldexp(1.0, exponent + cpower - 1), 0.0)
        B = (0.0, 0.0, ldexp(1.0, exponent))
        K, dt = ldexp(1.0, exponent + coefficient_shift), ldexp(1.0, exponent)
        expected = radiation_momentum_oracle(v, E, B, dt, K, c)
        actual, unchanged = radiation_apply(v, E, B, dt, K, c)
        @test unchanged
        @test radiation_error_ulps(actual, expected, v) <= 8
    end
    setprecision(BigFloat, 8192) do
        high, low = HybridPlasmaPIC._RADIATION_LOG2
        exact = log(BigFloat(2))
        @test high == Float64(exact)
        @test low == Float64(exact - BigFloat(high))
    end
end

@testset "Amplified radiation completion is independent across tasks" begin
    initial_precision, initial_rounding = precision(BigFloat), rounding(BigFloat)
    for (T, power) in ((Float16, 12), (Float32, 80), (Float64, 900))
        z, c = zero(T), one(T)
        large, small = ldexp(c, power), ldexp(c, -power)
        v, E, B = (T(0.5), z, z), (small, large / 2, z), (z, z, large)
        expected = radiation_momentum_oracle(v, E, B, large, large, c)
        initial = ParticleSet{3,T}(16)
        for d = 1:3
            fill!(initial.x[d], T(0.2))
            fill!(initial.v[d], v[d])
        end
        serial = deepcopy(initial)
        apply_radiation_reaction!(serial, E, B, large; K = large, c)
        @test radiation_error_ulps(ntuple(d -> serial.v[d][1], 3), expected, v) <= 8
        copies = [deepcopy(initial) for _ = 1:4]
        tasks =
            [Threads.@spawn apply_radiation_reaction!(q, E, B, large; K = large, c) for q in copies]
        @test all(fetch(task) === q for (task, q) in zip(tasks, copies))
        @test all(q.v == serial.v for q in copies)
        @test all(
            (q.x, q.weight, q.id, q.tag, q.q, q.m) ==
            (initial.x, initial.weight, initial.id, initial.tag, initial.q, initial.m) for
            q in copies
        )
    end
    @test precision(BigFloat) == initial_precision
    @test rounding(BigFloat) == initial_rounding
end
