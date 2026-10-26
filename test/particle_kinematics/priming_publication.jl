# Residual certificates contain the exact endpoint results, including uncertain
# denominators and underflow. These references use no production EFT/DD helper.
@testset "Priming residual intervals enclose exact arithmetic" begin
    exact_center(a) = Q(a[1][1]) + Q(a[1][2])
    endpoints(a) = (exact_center(a) - Q(a[2]), exact_center(a) + Q(a[2]))
    values = (
        0.0,
        nextfloat(0.0),
        -nextfloat(0.0),
        floatmin(Float64),
        1.0,
        prevfloat(1.0),
        nextfloat(1.0),
        -1.0,
        0x1p-500,
        0x1p500,
    )
    pairs = [((x, 0.0), (y, 0.0)) for x in values for y in values]
    rng = Xoshiro(0x726573696475616c)
    for _ = 1:128
        sample() = copysign(ldexp(1 + rand(rng), rand(rng, -500:500)), rand(rng, Bool) ? 1.0 : -1.0)
        x, y = sample(), sample()
        push!(pairs, ((x, eps(x) / 4), (y, -eps(y) / 4)))
    end
    accepted = 0
    for (aw, bw) in pairs
        a, b = (aw, abs(aw[1]) * eps(Float64)^2), (bw, abs(bw[1]) * eps(Float64)^2)
        al, ah = endpoints(a)
        bl, bh = endpoints(b)
        for (operation, reference) in (
            (H._rel_prime_ball_add, +),
            (H._rel_prime_ball_sub, -),
            (H._rel_prime_ball_mul, *),
            (H._rel_prime_ball_div, /),
        )
            result = operation(a, b)
            if H._rel_prime_ball_valid(result)
                operation === H._rel_prime_ball_div && @test !(bl <= 0 <= bh)
                targets = [reference(x, y) for x in (al, ah) for y in (bl, bh)]
                lo, hi = endpoints(result)
                @test lo <= minimum(targets) <= maximum(targets) <= hi
                accepted += 1
            end
        end
        if al > 0
            result = H._rel_prime_ball_sqrt(a)
            if H._rel_prime_ball_valid(result)
                lo, hi = endpoints(result)
                @test hi >= 0 && (lo <= 0 || lo^2 <= al) && hi^2 >= ah
                accepted += 1
            end
        end
    end
    @test accepted > length(pairs)
    @test !H._rel_prime_ball_valid(H._rel_prime_ball_div(((1.0, 0.0), 0.0), ((1.0, 0.0), 1.0)))
    @test !H._rel_prime_ball_valid(H._rel_prime_ball_sqrt(((0.0, 0.0), 1.0)))
    for T in (Float16, Float32, Float64)
        half = (Float64(nextfloat(one(T))) - 1.0) / 2
        radius = half / 16
        @test !last(H._rel_prime_round_cell(T, (1.0, half - radius / 2), radius))
        @test last(H._rel_prime_round_cell(T, (1.0, half - radius / 2), 0.0))
    end
end

@testset "Ill-conditioned priming publishes the independently rounded stage" begin
    for T in (Float16, Float32, Float64)
        rng = Xoshiro(0x7072696d65)
        span, fspan = T === Float16 ? (8, 6) : (40, 20)
        for i = 1:120
            planar = isodd(i)
            direction = randn(rng, 3)
            planar && (direction[3] = 0)
            direction /= norm(direction)
            speed = T(rand(rng, (0.125, 0.9, 0.99)))
            c = ldexp(one(T), rand(rng, -span:span))
            v = Tuple(T.(direction) .* speed .* c)
            qspan = T === Float16 ? 2 : 10
            qm =
                copysign(ldexp(one(T), rand(rng, -qspan:qspan)), rand(rng, Bool) ? one(T) : -one(T))
            h = ldexp(one(T), rand(rng, -qspan:qspan))
            sample() =
                copysign(ldexp(one(T), rand(rng, -fspan:fspan)), rand(rng, Bool) ? one(T) : -one(T))
            E = ntuple(j -> planar && j == 3 ? zero(T) : sample(), 3)
            B = ntuple(j -> planar && j != 3 ? zero(T) : sample(), 3)
            reference, _ = rel_prime_reference(v, E, B, qm, h, c)
            result, resolved = H._rel_prime_dd_attempt(ntuple(_ -> zero(T), 3), v, E, B, qm, h, c)
            @test resolved
            @test result == T.(reference)
        end
    end
end

@testset "Rounded-native priming threshold retains local normwise contract" begin
    rows = []
    for T in (Float16, Float32, Float64),
        k in (T === Float16 ? (-4, 0, 4) : T === Float32 ? (-40, 0, 40) : (-400, 0, 400)),
        raw in (
            (T(0.6), T(0.8), zero(T)),
            (one(T), zero(T), zero(T)),
            (nextfloat(one(T)), zero(T), zero(T)),
        ),
        order in ((1, 2, 3), (2, 1, 3), (3, 1, 2)),
        sign in (-one(T), one(T))

        c = ldexp(one(T), k)
        m = ntuple(i -> sign * raw[order[i]] * c, 3)
        E = map(-, m)
        z = ntuple(_ -> zero(T), 3)
        ref, scale = rel_prime_reference(z, E, z, one(T), one(T), c)
        for planar in (false, true)
            planar && !iszero(E[3]) && continue
            actual = ntuple(_ -> T[0], 3)
            old = map(copy, actual)
            if planar
                H._prime_empic1d!(actual, [E[1]], [E[2]], T[0], one(T), one(T), 1, true, c)
                old_prime_empic1d!(old, [E[1]], [E[2]], T[0], one(T), one(T), 1, true, c)
            else
                H._prime_leapfrog!(
                    actual,
                    ntuple(i -> [E[i]], 3),
                    ntuple(_ -> T[0], 3),
                    one(T),
                    one(T),
                    1,
                    true,
                    c,
                )
                old_prime!(
                    old,
                    ntuple(i -> [E[i]], 3),
                    ntuple(_ -> T[0], 3),
                    one(T),
                    one(T),
                    1,
                    true,
                    c,
                )
            end
            a, o = map(only, actual), map(only, old)
            @test rel_ok(a, ref, scale)
            native_beta = sum(abs2, o ./ sqrt(c * c))
            if native_beta <= T(0.5)
                @test isequal(a, o)
            else
                @test a == T.(ref)
            end
            physical = Q.(m)
            physical_beta = sum(abs2, physical) / (Q(c)^2 + sum(abs2, physical))
            push!(
                rows,
                Dict(
                    "type" => string(T),
                    "c_exponent" => k,
                    "planar" => planar,
                    "physical_above_half" => physical_beta > 1 // 2,
                    "native_beta" => Float64(native_beta),
                    "historical_bits" => isequal(a, o),
                    "correctly_rounded" => a == T.(ref),
                ),
            )
        end
    end
    @test any(
        row ->
            row["type"] == "Float64" &&
                row["physical_above_half"] &&
                row["native_beta"] <= 0.5 &&
                !row["correctly_rounded"],
        rows,
    )
end
