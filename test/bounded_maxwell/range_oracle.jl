# Independent Rational stencil and stage oracle. It does not call production
# difference/curl or exact-arithmetic helpers.
const YEE_Q = Rational{BigInt}

# Exact derivative terms of each curl component: curl = first - second.
function reference_terms(fields, g, electric)
    D = length(g.n)
    shapes =
        ntuple(c -> ntuple(d -> g.n[d] + Int(!g.periodic[d] && (electric ? c == d : c != d)), D), 3)
    first = ntuple(c -> Array{YEE_Q}(undef, shapes[c]), 3)
    second = ntuple(c -> Array{YEE_Q}(undef, shapes[c]), 3)
    function derivative(a, I, d)
        d > D && return zero(YEE_Q)
        i = I[d]
        if electric
            j = i == g.n[d] && g.periodic[d] ? 1 : i + 1
            neighbor = CartesianIndex(ntuple(k -> k == d ? j : I[k], D))
            return (YEE_Q(a[neighbor]) - YEE_Q(a[I])) / YEE_Q(g.dx[d])
        else
            !g.periodic[d] && (i == 1 || i == g.n[d] + 1) && return zero(YEE_Q)
            j = i == 1 ? g.n[d] : i - 1
            neighbor = CartesianIndex(ntuple(k -> k == d ? j : I[k], D))
            return (YEE_Q(a[I]) - YEE_Q(a[neighbor])) / YEE_Q(g.dx[d])
        end
    end
    for c = 1:3, I in CartesianIndices(first[c])
        j = c == 3 ? 1 : c + 1
        k = c == 1 ? 3 : c - 1
        first[c][I] = derivative(fields[k], I, j)
        second[c][I] = derivative(fields[j], I, k)
    end
    return first, second
end
reference_curl(fields, g, electric) = map((a, b) -> a .- b, reference_terms(fields, g, electric)...)

function reference_step(m, dt)
    T = eltype(m.fields.E[1])
    E, B = map(copy, m.fields.E), map(copy, m.fields.B)
    half = YEE_Q(dt) / 2
    CE = reference_curl(E, m.g, true)
    Bhalf = ntuple(c -> map((b, k) -> T(YEE_Q(b) - half * k), B[c], CE[c]), 3)
    CB = reference_curl(Bhalf, m.g, false)
    coefficient = YEE_Q(dt) * YEE_Q(m.c)^2
    Enext = ntuple(c -> map((e, k) -> T(YEE_Q(e) + coefficient * k), E[c], CB[c]), 3)
    CE = reference_curl(Enext, m.g, true)
    Bnext = ntuple(c -> map((b, k) -> T(YEE_Q(b) - half * k), Bhalf[c], CE[c]), 3)
    return Enext, Bnext
end

# Normwise stage contract from the stage's actual inputs: with
# κ_d = coefficient/dx_d and Δ_d the sample differences,
# |actual - exact| <= 8eps(T)*(|old| + Σ|κ_d Δ_d|) + nextfloat(0).
# Returns the largest error divided by that allowance (<= 1 passes).
function normwise_stage_ratio(actual, old, source, g, electric, dt, speed)
    T = eltype(old[1])
    coefficient = electric ? -YEE_Q(dt) / 2 : YEE_Q(dt) * YEE_Q(speed)^2
    first, second = reference_terms(source, g, electric)
    worst = 0.0
    for c = 1:3, I in CartesianIndices(old[c])
        exact = YEE_Q(old[c][I]) + coefficient * (first[c][I] - second[c][I])
        scale = abs(YEE_Q(old[c][I])) + abs(coefficient) * (abs(first[c][I]) + abs(second[c][I]))
        allowance = 8YEE_Q(eps(T)) * scale + YEE_Q(nextfloat(zero(T)))
        isfinite(actual[c][I]) || return Inf
        worst = max(worst, Float64(abs(YEE_Q(actual[c][I]) - exact) / allowance))
    end
    return worst
end

# Check each stage of one step from its actual stored inputs, then require
# the public step to publish exactly the composition of those stages.
function normwise_step_check(m, dt)
    H = HybridPlasmaPIC
    E, B = map(copy, m.fields.E), map(copy, m.fields.B)
    half = map(copy, B)
    H._yee_advance!(half, B, E, m.g, true, dt, m.c)
    @test normwise_stage_ratio(half, B, E, m.g, true, dt, m.c) <= 1
    next_E = map(copy, E)
    H._yee_advance!(next_E, E, half, m.g, false, dt, m.c)
    @test normwise_stage_ratio(next_E, E, half, m.g, false, dt, m.c) <= 1
    next_B = map(copy, half)
    H._yee_advance!(next_B, half, next_E, m.g, true, dt, m.c)
    @test normwise_stage_ratio(next_B, half, next_E, m.g, true, dt, m.c) <= 1
    step_maxwell!(m, dt)
    @test isequal(m.fields.E, next_E) && isequal(m.fields.B, next_B)
    @test all(a -> all(isfinite, a), (m.fields.E..., m.fields.B...))
    return next_E, next_B
end

@testset "Complete finite Maxwell stage range" begin
    for T in (Float32, Float64),
        D = 1:3,
        periodic in (false, true),
        profile in (:c_large, :c_small, :derivative_large, :derivative_small, :half_subnormal)

        wide = T === Float32 ? T(1e20) : T(1e200)
        small = inv(wide)
        speed = profile === :c_large ? wide : profile === :c_small ? small : one(T)
        length =
            profile === :derivative_large ? small : profile === :derivative_small ? wide : one(T)
        g = YeeGrid(
            ntuple(_ -> 4, D),
            ntuple(_ -> length, D);
            boundaries = ntuple(_ -> periodic ? :periodic : :pec, D),
        )
        m = MaxwellVacuum(g; c = speed)
        dt = profile === :half_subnormal ? nextfloat(zero(T)) : T(0.1) * maxwell_cfl(m)
        amplitude =
            profile === :c_large ? small :
            profile === :c_small ? wide :
            profile === :derivative_large ? wide :
            profile === :derivative_small ? small : floatmax(T) / 16
        for c = 1:3, I in CartesianIndices(m.fields.B[c])
            m.fields.B[c][I] = T(mod(sum(Tuple(I)) + c, 3) - 1) * amplitude
        end
        # Derivative and half-step cases must also exercise the magnetic update.
        if profile in (:derivative_large, :derivative_small, :half_subnormal)
            for c = 1:3, I in CartesianIndices(m.fields.E[c])
                m.fields.E[c][I] = T(mod(sum(Tuple(I)) + c, 3) - 1) * amplitude
            end
        end
        apply_yee_pec!(m.fields, g)
        expectedE, expectedB = reference_step(m, dt)
        @test all(a -> all(isfinite, a), (expectedE..., expectedB...))
        normwise_step_check(m, dt)
        # Whole-step global-norm comparison with the once-rounded exact stages.
        for (actual, expected) in zip(components(m.fields), (expectedE..., expectedB...))
            scale = maximum(abs, expected)
            @test maximum(abs, actual - expected) <= 32eps(T) * scale + 4nextfloat(zero(T))
        end
        @test m.time[] == dt && m.step[] == 1
    end
end
