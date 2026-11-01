# Manual observation only: exit 0 does not certify relativistic accuracy.
using HybridPlasmaPIC, LinearAlgebra
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# Existing normal-range nonrelativistic failures from test_boris_accuracy.jl.
# This finite-c adaptation is an unverified candidate until this script runs.
const ROWS = (
    (
        0x3e80, 0x39d7, (0x38b0, 0x3d4c, 0xbe81),
        (0x3a90, 0xb8e1, 0xbdc8), (0x3e7a, 0x3eb6, 0x3ffc),
    ),
    (
        0x3c93, 0x3cea, (0x3ba7, 0xbe0d, 0x3fbb),
        (0x3def, 0xbfad, 0x3a42), (0xbac0, 0x39bb, 0x3e47),
    ),
    (
        0x3af5, 0x38c1, (0xbefd, 0x3fa6, 0x3ab6),
        (0xb8e1, 0x3e23, 0x3ebb), (0xbaf9, 0x3ec5, 0xbff8),
    ),
    (
        0x3f8f, 0x3ddb, (0xbd53, 0x3c71, 0xbeea),
        (0xbcad, 0x3f01, 0xbc14), (0x3e7c, 0x3d3e, 0x3fcf),
    ),
)

# Independent Cayley rotation, not the production two-cross-product kernel.
# Original factors and the original-component Lorentz gap enter as rationals.
function reference(bits, ::Type{T}, v, E, B, ratio, duration, c, prime) where {T}
    return setprecision(BigFloat, bits) do
        setrounding(BigFloat, RoundNearest) do
            gap = 1 - sum(x -> Q(x)^2, v) / Q(c)^2
            V, A, C = map(x -> BigFloat.(collect(x)), (v, E, B))
            initial = V / sqrt(BigFloat(gap))
            coefficient = BigFloat(Q(duration) * ratio)
            if prime
                impulse = coefficient .* (A + cross(V, C))
                momentum = initial - impulse
            else
                electric = (coefficient / 2) .* A
                minus = initial + electric
                gamma_minus = sqrt(1 + dot(minus, minus) / BigFloat(c)^2)
                t = (coefficient / 2) .* C / gamma_minus
                t2 = dot(t, t)
                rotated =
                    ((1 - t2) .* minus + 2dot(t, minus) .* t + 2cross(minus, t)) / (1 + t2)
                momentum = rotated + electric
                impulse = electric
            end
            gamma_final = sqrt(1 + dot(momentum, momentum) / BigFloat(c)^2)
            values = Tuple(momentum / gamma_final)
            scale = max(norm(collect(values)), max(norm(initial), norm(impulse)) / gamma_final)
            power = max(
                iszero(scale) ? exponent(floatmin(T)) : exponent(scale),
                exponent(floatmin(T)),
            )
            budget = BigFloat(2)^(power - precision(T) + 3)
            (; values, scale, budget, rounded = T.(values))
        end
    end
end

function comparison(actual, ref)
    return setprecision(BigFloat, precision(ref.budget)) do
        setrounding(BigFloat, RoundNearest) do
            error = map((a, r) -> abs(BigFloat(a) - r), actual, ref.values)
            (; error, violation = !all(x -> isfinite(x) && x <= ref.budget, error))
        end
    end
end

# The general priming owner's historical native trial, for range classification.
# Actual publication below calls the shared production priming entry point.
function native_prime(v, E, B, factor, h, c)
    qm = H._boris_native_factor(factor)
    beta = v ./ c
    beta2 = beta[1]^2 + beta[2]^2 + beta[3]^2
    beta2 <= oftype(c, 0.5) || error("Fixture exceeds native Lorentz-gap conditioning slice")
    gamma = one(c) / sqrt(one(c) - beta2)
    force = (
        E[1] + (v[2] * B[3] - v[3] * B[2]),
        E[2] + (v[3] * B[1] - v[1] * B[3]),
        E[3] + (v[1] * B[2] - v[2] * B[1]),
    )
    a = qm .* force
    initial, impulse = gamma .* v, h .* a
    momentum = initial .- impulse
    gamma_final = sqrt(one(c) + (momentum[1]^2 + momentum[2]^2 + momentum[3]^2) / (c * c))
    out = momentum ./ gamma_final
    admitted = H._rel_prime_in_range(out, initial, a, impulse, gamma_final, h, factor, c * c, true)
    return out, admitted
end

function observe(::Type{T}, index, q, m, v, E, B, dt, c, meaning) where {T}
    qm = q / m
    factor = meaning === :original ? H._BorisSpeciesFactors(q, m, qm) : qm
    ratio = meaning === :original ? Q(q) / Q(m) : Q(qm)
    Ep, Bp = map(x -> [x], E), map(x -> [x], B)
    kick, parts, u, c2, gap = H._boris_rel_parts(v..., E..., B..., qm, dt, c, Val(false))
    kick_range = H._euler_factor_in_range(factor, T) &&
                 H._boris_rel_in_range(kick, parts, u, qm, dt, c2, gap)
    scales = H._boris_batch_scales(map(x -> [x], v), E, B, qm * dt / 2, 1:1)
    kick_batch = H._boris_rel_batch_certificate(qm * dt / 2, factor, dt, c, scales...)
    scalar = H._boris_rel(v..., E..., B..., factor, dt, c)
    native_published = H._boris_rel_native(v..., E..., B..., factor, dt, c)
    prime, prime_range = native_prime(v, E, B, factor, dt / 2, c)
    prime_scales = H._euler_batch_scales(map(x -> [x], v), E, B, (dt / 2, qm), 1:1, T)
    prime_batch = H._rel_prime_batch_certificate(dt / 2, factor, c, prime_scales...)
    primed = map(x -> [x], v)
    H._prime_leapfrog!(primed, Ep, Bp, factor, dt / 2, 1, true, c)
    observed = (
        (:kick, dt, kick_range, kick_batch,
         (raw_native = kick, guarded = scalar, batch_native = native_published)),
        (:prime, dt / 2, prime_range, prime_batch,
         (raw_native = prime, guarded = map(only, primed))),
    )
    violations = unresolved = 0
    for (stage, duration, range_admitted, batch_admitted, actual) in observed
        low = reference(512, T, v, E, B, ratio, duration, c, stage === :prime)
        high = reference(1024, T, v, E, B, ratio, duration, c, stage === :prime)
        comparisons = map(x -> comparison(x, low), actual)
        refined = map(x -> comparison(x, high), actual)
        stable = isequal(low.rounded, high.rounded) && low.budget == high.budget &&
                 map(x -> x.violation, comparisons) == map(x -> x.violation, refined)
        status = !stable ? :UNRESOLVED :
                 any(x -> x.violation, refined) ? :OBSERVED_VIOLATION : :NO_VIOLATION_OBSERVED
        raw_admission_violation = stable && range_admitted && refined.raw_native.violation
        production_violation = stable && (
            refined.guarded.violation ||
            (stage === :kick && batch_admitted && refined.batch_native.violation)
        )
        unresolved += !stable
        violations += production_violation
        println((;
            T, index, meaning, stage, q, m, qm, v, E, B, dt, c, ratio, duration,
            range_admitted, batch_admitted, actual, reference_512 = low,
            reference_1024 = high, errors_512 = comparisons, errors_1024 = refined,
            stable, raw_admission_violation, production_violation, status,
        ))
    end
    return violations, unresolved
end

println("Observation only. Raw kick/scalar/native publication and raw/shared priming are separate outputs. A raw native violation outside its admission is not a production admission failure.")
println((; julia = VERSION, rounding = rounding(Float64), zero_subnormals = get_zero_subnormals()))
let violations = 0, unresolved = 0, cases = 0
    for T in (Float16, Float32, Float64), meaning in (:original, :explicit)
        cv(x) = T(reinterpret(Float16, x))
        for (index, row) in enumerate(ROWS)
            q, m = cv(row[1]), cv(row[2])
            v, E, B = cv.(row[3]) ./ T(4), cv.(row[4]) ./ T(4), cv.(row[5])
            found, undecided = observe(T, index, q, m, v, E, B, T(2), T(128), meaning)
            violations += found
            unresolved += undecided
            cases += 1
        end
        found, undecided = observe(
            T, :ordinary, one(T), one(T), T.((0.5, -0.25, 0.125)),
            T.((0.125, 0.0625, -0.03125)), T.((0.25, -0.125, 0.0625)),
            T(0.125), T(128), meaning,
        )
        violations += found
        unresolved += undecided
        cases += 1
    end
    println((;
        cases, observed_production_stage_violations = violations, unresolved_stages = unresolved,
    ))
    iszero(unresolved) || error(
        "Independent 512/1024-bit references disagree; no accuracy classification accepted for unresolved stages",
    )
end
println("Completed observation; successful job status is not accuracy acceptance. Review each admission flag and each actual result separately.")
