using HybridPlasmaPIC, Test, Random
include("sourced_yee_oracle.jl")

const SY_Q = Rational{BigInt}
const sy_rat = SourcedYeeOracle.rat

# Normwise stage contract: 8eps(T)*(|old| + |κΔ| + |impulse|) + one subnormal quantum.
sy_stage_close(value::T, exact, scale) where {T} =
    isfinite(value) && abs(SY_Q(value) - exact) <= 8SY_Q(eps(T)) * scale + SY_Q(nextfloat(zero(T)))

# Reduced scalar contract: (8 + terms)eps(T)*scale + one subnormal quantum. A
# nonfinite result must be the rounded exact value (outside the range).
function sy_scalar_close(value::T, exact, scale, terms) where {T}
    isfinite(value) || return isequal(value, T(exact))
    return abs(SY_Q(value) - exact) <= (8 + terms) * SY_Q(eps(T)) * scale + SY_Q(nextfloat(zero(T)))
end
sy_terms(E, B) = sum(length, E) + 2sum(length, B)

function sy_populate!(m, rng)
    for a in (m.fields.E..., m.fields.B...)
        a .= rand(rng, eltype(a), size(a)) .- eltype(a)(0.5)
    end
    apply_yee_pec!(m.fields, m.g)
    return m
end

sy_snapshot(m) = (map(copy, m.fields.E), map(copy, m.fields.B), m.time[], m.step[])
sy_unchanged(m, old) =
    isequal(m.fields.E, old[1]) &&
    isequal(m.fields.B, old[2]) &&
    isequal(m.time[], old[3]) &&
    m.step[] == old[4]

# Check each stage from its actual stored inputs with the independent rational
# oracle, then require the public step to publish exactly their composition.
function sy_normwise_step!(m, h, current)
    H = HybridPlasmaPIC
    g = m.g
    n, periodic, dx = g.n[1], g.periodic[1], g.dx[1]
    E, B = map(copy, m.fields.E), map(copy, m.fields.B)
    half = map(copy, B)
    H._yee_advance!(half, B, E, g, true, h, m.c)
    next_E = map(copy, E)
    H._yee_advance!(next_E, E, half, g, false, h, m.c, current)
    next_B = map(copy, half)
    H._yee_advance!(next_B, half, next_E, g, true, h, m.c)
    q, κ = sy_rat(h) / 2, sy_rat(h) * sy_rat(m.c)^2
    ce = SourcedYeeOracle.curl(E, n, periodic, true, dx)
    cb = SourcedYeeOracle.curl(half, n, periodic, false, dx)
    cn = SourcedYeeOracle.curl(next_E, n, periodic, true, dx)
    for c = 1:3
        old_B, old_E, impulse = sy_rat.(B[c]), sy_rat.(E[c]), sy_rat.(current[c])
        @test all(sy_stage_close.(half[c], old_B .- q .* ce[c], abs.(old_B) .+ q .* abs.(ce[c])))
        @test all(
            sy_stage_close.(
                next_E[c],
                old_E .+ κ .* cb[c] .- impulse,
                abs.(old_E) .+ κ .* abs.(cb[c]) .+ abs.(impulse),
            ),
        )
        mid = sy_rat.(half[c])
        @test all(sy_stage_close.(next_B[c], mid .- q .* cn[c], abs.(mid) .+ q .* abs.(cn[c])))
    end
    @test step_maxwell!(m, h, current) === m
    @test isequal(m.fields.E, next_E) && isequal(m.fields.B, next_B)
    return half
end

@testset "Prescribed-current Yee fields" begin
    rng = MersenneTwister(260170)
    @testset "Independent stages, charge and work" begin
        for T in (Float32, Float64), periodic in (false, true), n in (2, 3, 8, 17), trial = 1:12
            g = YeeGrid((n,), (T(2),); boundaries = (periodic ? :periodic : :pec,))
            m = sy_populate!(MaxwellVacuum(g; c = T(1.25)), rng)
            current = map(a -> T(0.01) .* randn(rng, T, size(a)), m.fields.E)
            if !periodic
                for c = 2:3
                    current[c][1] = current[c][end] = zero(T)
                end
            end
            old = sy_snapshot(m)
            saved_current = map(copy, current)
            ids = map(objectid, (m.fields.E..., m.fields.B...))
            h = T(0.375) * maxwell_cfl(m)
            exact = SourcedYeeOracle.step(old[1], old[2], current, n, periodic, g.dx[1], h, m.c)
            terms = sy_terms(old[1], old[2])
            h0 = yee_modified_energy(m, h)
            @test sy_scalar_close(
                h0,
                SourcedYeeOracle.energy(old[1], old[2], n, periodic, g.dx[1], h, m.c),
                SourcedYeeOracle.energy_scale(old[1], old[2], n, periodic, g.dx[1], h, m.c),
                terms,
            )
            sy_normwise_step!(m, h, current)
            for c = 1:3
                @test maximum(abs, m.fields.E[c] - exact[1][c]) <=
                      T(64) * eps(T) * max(one(T), maximum(abs, exact[1][c]))
            end
            @test ids == map(objectid, (m.fields.E..., m.fields.B...))
            @test current == saved_current
            @test m.time[] == h && m.step[] == 1
            q = zeros(T, n + Int(!periodic))
            q0 = copy(q)
            qi = copy(q)
            @test yee_electric_charge!(q, m.fields.E, g) === q
            yee_electric_charge!(q0, old[1], g)
            yee_electric_charge!(qi, current, g)
            @test q == T.(SourcedYeeOracle.charge(m.fields.E, periodic))
            @test maximum(abs, q - q0 + qi) <= T(8) * eps(T)
            work = yee_source_work(old[1], m.fields.E, current, g)
            @test sy_scalar_close(
                work,
                SourcedYeeOracle.work(old[1], m.fields.E, current, g.dx[1]),
                SourcedYeeOracle.work_scale(old[1], m.fields.E, current, g.dx[1]),
                sum(length, current),
            )
            h1 = yee_modified_energy(m, h)
            @test sy_scalar_close(
                h1,
                SourcedYeeOracle.energy(m.fields.E, m.fields.B, n, periodic, g.dx[1], h, m.c),
                SourcedYeeOracle.energy_scale(m.fields.E, m.fields.B, n, periodic, g.dx[1], h, m.c),
                terms,
            )
            @test abs(h1 - h0 + work) <= T(64) * eps(T) * max(one(T), h0, h1, abs(work))
        end
    end
    @testset "Vacuum equality and periodic mean current" begin
        for T in (Float32, Float64), p in (false, true), n in (2, 8, 32)
            g = YeeGrid((n,), (one(T),); boundaries = (p ? :periodic : :pec,))
            a = sy_populate!(MaxwellVacuum(g), rng)
            b = deepcopy(a)
            c = deepcopy(a)
            current = map(x -> zero(x), a.fields.E)
            h = maxwell_cfl(a) / T(4)
            step_maxwell!(a, h)
            step_maxwell!(b, h, current)
            @test isequal(a.fields.E, b.fields.E) && isequal(a.fields.B, b.fields.B)
            # Negative zero impulses also reproduce the vacuum values bitwise.
            step_maxwell!(c, h, map(x -> -x, current))
            @test isequal(a.fields.E, c.fields.E) && isequal(a.fields.B, c.fields.B)
            m = MaxwellVacuum(g)
            current[1] .= T(0.125)
            step_maxwell!(m, h, current)
            @test all(==(-T(0.125)), m.fields.E[1])
            charge = zeros(T, n + Int(!p))
            yee_electric_charge!(charge, m.fields.E, g)
            @test charge == (p ? zeros(T, n) : [T(-0.125); zeros(T, n - 1); T(0.125)])
        end
    end
    @testset "Complete-expression range and cancellation" begin
        # Both derivative numerators and old-impulse overflow here, while the
        # complete electric fields are exactly a (exact completion).
        for T in (Float32, Float64), sign in (-1, 1), n in (2, 4, 8)
            g = YeeGrid((n,), (T(n),); boundaries = (:periodic,))
            m = MaxwellVacuum(g)
            a = T(sign) * T(0.75) * floatmax(T)
            m.fields.E[2] .= a
            current = map(x -> zero(x), m.fields.E)
            for i = 1:n
                m.fields.B[3][i] = isodd(i) ? a : -a
                current[2][i] = isodd(i) ? -a : a
            end
            h = T(0.5)
            @test !isfinite(a - current[2][1])
            sy_normwise_step!(m, h, current)
            @test m.fields.E[2] == fill(a, n)
            @test all(isfinite, m.fields.B[3])
        end
        for T in (Float32, Float64), sign in (-1, 1), scale in (T(0.5), T(0.75))
            g = YeeGrid((2,), (T(2),); boundaries = (:periodic,))
            m = MaxwellVacuum(g)
            a = T(sign) * scale * floatmax(T)
            m.fields.E[2] .= a
            m.fields.B[3] .= [-a, a]
            current = map(x -> zero(x), m.fields.E)
            current[2] .= a
            sy_normwise_step!(m, T(0.25), current)
            @test all(isfinite, m.fields.E[2]) && all(isfinite, m.fields.B[3])
        end
        for T in (Float32, Float64), offset in (-8, -2, -1, 0, 1, 2, 8), exp in (-10, 0, 10)
            g = YeeGrid((2,), (T(2),); boundaries = (:periodic,))
            m = MaxwellVacuum(g)
            e = eps(T)
            m.fields.E[2] .= ldexp(one(T) + e, exp)
            m.fields.B[3] .= T[0, ldexp(one(T) - e, exp)]
            current = map(x -> zero(x), m.fields.E)
            current[2] .= ldexp(T(1.25) + T(offset) * e, exp)
            sy_normwise_step!(m, T(0.25), current)
        end
        for T in (Float32, Float64), k in (-4, -1, 0, 1, 4)
            g = YeeGrid((2,), (T(2),); boundaries = (:periodic,))
            m = MaxwellVacuum(g)
            impulse = map(x -> zero(x), m.fields.E)
            impulse[1] .= T(k) * nextfloat(zero(T))
            step_maxwell!(m, T(0.25), impulse)
            @test m.fields.E[1] == -impulse[1]
        end
    end
    @testset "Scalar original-input diagnostic oracles" begin
        for T in (Float32, Float64), n in (2, 4, 8), exponent in (-100, 0, 100), sign in (-1, 1)
            g = YeeGrid((n,), (T(n),); boundaries = (:periodic,))
            m = MaxwellVacuum(g)
            x = ldexp(one(T), exponent)
            m.fields.E[2] .= x
            old = map(copy, m.fields.E)
            new = map(copy, m.fields.E)
            impulse = map(a -> zero(a), old)
            impulse[2] .= T(sign) * x
            impulse[2][1] = -T(sign) * x * T(n - 1)
            new[2][end] = nextfloat(x)
            @test sy_scalar_close(
                yee_source_work(old, new, impulse, g),
                SourcedYeeOracle.work(old, new, impulse, g.dx[1]),
                SourcedYeeOracle.work_scale(old, new, impulse, g.dx[1]),
                sum(length, impulse),
            )
            for i = 1:n
                m.fields.E[2][i] = isodd(i) ? x : -x
            end
            h = prevfloat(maxwell_cfl(m))
            @test sy_scalar_close(
                yee_modified_energy(m, h),
                SourcedYeeOracle.energy(m.fields.E, m.fields.B, n, true, g.dx[1], h, m.c),
                SourcedYeeOracle.energy_scale(m.fields.E, m.fields.B, n, true, g.dx[1], h, m.c),
                sy_terms(m.fields.E, m.fields.B),
            )
        end
        for T in (Float32, Float64)
            g = YeeGrid((2,), (T(2),); boundaries = (:periodic,))
            m = MaxwellVacuum(g)
            old = map(a -> zero(a), m.fields.E)
            new = map(copy, old)
            current = map(copy, old)
            old[1] .= floatmax(T)
            new[1] .= floatmax(T)
            current[1] .= T[1, -1]
            @test yee_source_work(old, new, current, g) == zero(T)
            current[1] .= one(T)
            @test yee_source_work(old, new, current, g) == T(Inf)
            m.fields.E[1] .= floatmax(T)
            @test yee_modified_energy(m, T(0.25)) == T(Inf)
        end
    end
    @testset "Validation and physical rollback" begin
        for T in (Float32, Float64), p in (false, true)
            g = YeeGrid((4,), (one(T),); boundaries = (p ? :periodic : :pec,))
            m = sy_populate!(MaxwellVacuum(g), rng)
            h = maxwell_cfl(m) / T(4)
            valid = map(a -> zero(a), m.fields.E)
            invalid = Any[
                nothing,
                [valid...],
                (valid[1], valid[2]),
                (zeros(T, 3), valid[2], valid[3]),
                map(a -> Float16.(a), valid),
                (view(valid[1], :), valid[2], valid[3]),
                (m.fields.E[1], valid[2], valid[3]),
                (m._candidate.E[1], valid[2], valid[3]),
                (valid[1], valid[2], valid[2]),
            ]
            nonfinite = map(copy, valid)
            nonfinite[1][end] = T(NaN)
            push!(invalid, nonfinite)
            if !p
                for c = 2:3, side in (1, length(valid[c]))
                    wall = map(copy, valid)
                    wall[c][side] = one(T)
                    push!(invalid, wall)
                end
            end
            for current in invalid
                saved = sy_snapshot(m)
                scratch = deepcopy(m._candidate)
                @test_throws Union{ArgumentError,DimensionMismatch} step_maxwell!(m, h, current)
                @test sy_unchanged(m, saved)
                @test isequal(m._candidate.E, scratch.E) && isequal(m._candidate.B, scratch.B)
            end
            for dt in (zero(T), -h, T(NaN), T(Inf), maxwell_cfl(m))
                saved = sy_snapshot(m)
                @test_throws ArgumentError step_maxwell!(m, dt, valid)
                @test sy_unchanged(m, saved)
                @test_throws ArgumentError yee_modified_energy(m, dt)
            end
            m.fields.E[1] .= -floatmax(T)
            valid[1] .= floatmax(T)
            saved = sy_snapshot(m)
            @test_throws ArgumentError step_maxwell!(m, h, valid)
            @test sy_unchanged(m, saved)
            m.fields.E[1] .= zero(T)
            valid[1] .= zero(T)
            @test step_maxwell!(m, h, valid) === m
            m.time[] = floatmax(T)
            saved = sy_snapshot(m)
            @test_throws ArgumentError step_maxwell!(m, h, valid)
            @test sy_unchanged(m, saved)
            m.time[] = zero(T)
            m.step[] = typemax(Int)
            saved = sy_snapshot(m)
            @test_throws ArgumentError step_maxwell!(m, h, valid)
            @test sy_unchanged(m, saved)
            out = fill(T(7), 4 + Int(!p))
            prior = copy(out)
            broken = map(copy, m.fields.E)
            broken[1][1] = T(Inf)
            @test_throws ArgumentError yee_electric_charge!(out, broken, g)
            @test out == prior
            @test_throws DimensionMismatch yee_electric_charge!(zeros(T, 1), m.fields.E, g)
            alias = p ? m.fields.E[1] : m.fields.E[2]
            @test_throws ArgumentError yee_electric_charge!(alias, m.fields.E, g)
            @test_throws ArgumentError yee_source_work(broken, m.fields.E, valid, g)
        end
        m = MaxwellVacuum(YeeGrid((2, 2), (1.0, 1.0)))
        @test_throws ArgumentError step_maxwell!(m, 0.1, map(a -> zero(a), m.fields.E))
    end
end
