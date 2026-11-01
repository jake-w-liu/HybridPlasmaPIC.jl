module BorisAccuracyTests
using HybridPlasmaPIC, Test, LinearAlgebra
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
uinttype(::Type{Float16}) = UInt16
uinttype(::Type{Float32}) = UInt32
uinttype(::Type{Float64}) = UInt64
function binade(q::Q)
    n, d = numerator(q), denominator(q)
    e = ndigits(n; base = 2) - ndigits(d; base = 2)
    return e - (e >= 0 ? n < d << e : n << -e < d)
end
function exact_reference(v, E, B, qm, dt)
    v, e, b = map(x -> Q.(collect(x)), (v, E, B))
    a = Q(qm) * Q(dt) / 2
    k, t = a .* e, a .* b
    w = v .+ k
    t2 = dot(t, t)
    out = ((1 - t2) .* w .+ 2dot(t, w) .* t .+ 2cross(w, t)) ./ (1 + t2) .+ k
    return Tuple(out), max(dot(v, v), dot(k, k))
end
function euler_reference(v, E, B, q, m, h, sign)
    u, e, b = map(x -> Q.(collect(x)), (v, E, B))
    k = (Q(h) * Q(q) / Q(m)) .* (e .+ cross(u, b))
    return Tuple(u .+ sign .* k)
end
function contract(actual::NTuple{3,T}, exact, scale2) where {T}
    emin = exponent(floatmin(T))
    e = iszero(scale2) ? emin : max(fld(binade(scale2), 2), emin)
    tol = 4Q(2)^(e - precision(T) + 1)
    threshold = Q(floatmax(T)) + Q(2)^(exponent(floatmax(T)) - precision(T))
    return all(1:3) do c
        x = exact[c]
        abs(x) >= threshold ? isinf(actual[c]) && signbit(actual[c]) == (x < 0) :
        isfinite(actual[c]) && abs(Q(actual[c]) - x) <= tol
    end
end
# Independent exact nearest-even reference: search monotone positive IEEE
# encodings instead of using the production significand/divrem converter.
function nearest(::Type{T}, x::Q) where {T}
    iszero(x) && return zero(T)
    U = uinttype(T)
    a = abs(x)
    maximum = Q(floatmax(T))
    threshold = maximum + Q(2)^(exponent(floatmax(T)) - precision(T))
    a >= threshold && return x < 0 ? -T(Inf) : T(Inf)
    low, high = zero(U), reinterpret(U, floatmax(T))
    while low < high
        middle = low + (high - low + one(U)) ÷ U(2)
        Q(reinterpret(T, middle)) <= a ? (low = middle) : (high = middle - one(U))
    end
    if low == reinterpret(U, floatmax(T))
        value = floatmax(T)
    else
        lower, upper = reinterpret(T, low), reinterpret(T, low + one(U))
        dl, du = a - Q(lower), Q(upper) - a
        value = (dl < du || (dl == du && iseven(low))) ? lower : upper
    end
    return x < 0 ? -value : value
end

function stored_rows(::Type{T}) where {T}
    # These are the unchanged four original Float16 physical failures. Wider
    # tests are exact widenings of those same stored values.
    raw = (
        (
            0x3e80,
            0x39d7,
            (0x38b0, 0x3d4c, 0xbe81),
            (0x3a90, 0xb8e1, 0xbdc8),
            (0x3e7a, 0x3eb6, 0x3ffc),
        ),
        (
            0x3c93,
            0x3cea,
            (0x3ba7, 0xbe0d, 0x3fbb),
            (0x3def, 0xbfad, 0x3a42),
            (0xbac0, 0x39bb, 0x3e47),
        ),
        (
            0x3af5,
            0x38c1,
            (0xbefd, 0x3fa6, 0x3ab6),
            (0xb8e1, 0x3e23, 0x3ebb),
            (0xbaf9, 0x3ec5, 0xbff8),
        ),
        (
            0x3f8f,
            0x3ddb,
            (0xbd53, 0x3c71, 0xbeea),
            (0xbcad, 0x3f01, 0xbc14),
            (0x3e7c, 0x3d3e, 0x3fcf),
        ),
    )
    cv(x) = T(reinterpret(Float16, x))
    return [
        (q = cv(r[1]), m = cv(r[2]), v = cv.(r[3]), E = cv.(r[4]), B = cv.(r[5]), dt = T(2)) for
        r in raw
    ]
end
function check_case(::Type{T}, q, m, v, E, B, dt) where {T}
    factor = H._BorisSpeciesFactors(q, m, q / m)
    @test isbitstype(typeof(factor)) && sizeof(factor) == 3sizeof(T)
    native, parts, nodes = H._boris_kick_parts(v..., E..., B..., q / m, dt, Val(true))
    @test isequal(H._boris_kick_parts(v..., E..., B..., q / m, dt), (native, parts))
    @test length(nodes) == 48
    for meaning in (:original, :explicit)
        f = meaning === :original ? factor : q / m
        !isfinite(q / m) && meaning === :explicit && continue
        exact, scale2 = exact_reference(v, E, B, meaning === :original ? Q(q) / Q(m) : Q(q / m), dt)
        actual =
            meaning === :original ? H._kick_eval(Val(true), v..., E..., B..., factor, dt) :
            boris_kick(v..., E..., B..., q / m, dt)
        @test contract(actual, exact, scale2)
        range = H._boris_kick_in_range(native, parts, v..., q / m, dt)
        if range && contract(native, exact, scale2)
            @test isequal(actual, native)
        else
            @test isequal(actual, map(x -> nearest(T, x), exact))
        end
    end
    return nothing
end

function check_public_case(::Type{T}, row) where {T}
    for D = 1:3, route in (:uniform, :gathered), stride in (false, true)
        n = 3
        base = ParticleSet{D,T}(n; q = row.q, m = row.m)
        foreach(x -> fill!(x, zero(T)), base.x)
        ps = if stride
            x = ntuple(_ -> view(zeros(T, 2n + 1), 2:2:2n), D)
            v = ntuple(_ -> view(zeros(T, 2n + 1), 2:2:2n), 3)
            ParticleSet{D,T}(x, v, base.weight, base.id, base.tag, row.q, row.m)
        else
            base
        end
        for c = 1:3
            ps.v[c] .= row.v[c]
        end
        E = map(e -> fill(e, n), row.E)
        B = map(b -> fill(b, n), row.B)
        midpoint = ntuple(_ -> stride ? view(zeros(T, 2n + 1), 2:2:2n) : zeros(T, n), D)
        if route === :uniform
            push_uniform!(ps, row.E, row.B, row.dt)
        else
            push_gathered!(ps, E, B, row.dt; xmid = midpoint)
        end
        exact, scale2 = exact_reference(row.v, row.E, row.B, Q(row.q) / Q(row.m), row.dt)
        @test all(p -> contract(ntuple(c -> ps.v[c][p], 3), exact, scale2), 1:n)
        for p = 1:n, d = 1:D
            @test isequal(ps.x[d][p], nearest(T, Q(row.dt) * Q(ps.v[d][p])))
            if route === :gathered
                @test isequal(midpoint[d][p], nearest(T, Q(row.dt / T(2)) * Q(ps.v[d][p])))
            end
        end
    end
end

# Keep neighboring generic arithmetic's observable variadic dispatch. The
# wrapped exact rationals also permit an independent Cayley physics check.
const generic_binary_adds = Ref(0)
const generic_ternary_adds = Ref(0)
struct GenericBorisNumber <: Number
    value::Q
end
Base.:+(a::GenericBorisNumber, b::GenericBorisNumber) =
    (generic_binary_adds[] += 1; GenericBorisNumber(a.value + b.value))
Base.:+(a::Integer, b::GenericBorisNumber) =
    (generic_binary_adds[] += 1; GenericBorisNumber(Q(a) + b.value))
Base.:+(a::GenericBorisNumber, b::GenericBorisNumber, c::GenericBorisNumber) =
    (generic_ternary_adds[] += 1; GenericBorisNumber(a.value + b.value + c.value))
Base.:-(a::GenericBorisNumber, b::GenericBorisNumber) = GenericBorisNumber(a.value - b.value)
Base.:*(a::GenericBorisNumber, b::GenericBorisNumber) = GenericBorisNumber(a.value * b.value)
Base.:/(a::GenericBorisNumber, b::GenericBorisNumber) = GenericBorisNumber(a.value / b.value)
Base.:/(a::GenericBorisNumber, b::Integer) = GenericBorisNumber(a.value / Q(b))
Base.:/(a::Integer, b::GenericBorisNumber) = GenericBorisNumber(Q(a) / b.value)

@noinline function retention_in_mode(
    out::NTuple{3,T},
    exact::NTuple{3,Q},
    scale2::Q,
    flush::Bool,
) where {T}
    saved = get_zero_subnormals()
    try
        set_zero_subnormals(flush)
        verdict = H._boris_exact_native_ok(out, exact, scale2)
        return verdict, get_zero_subnormals()
    finally
        set_zero_subnormals(saved)
    end
end

@testset "Boris retention compares stored subnormal bits" begin
    exact = (Q(0), Q(0), Q(0))
    scale2 = Q(0)
    saved = get_zero_subnormals()
    for T in (Float16, Float32, Float64), magnitude in (4, 5)
        # At zero scale the contract permits four minimum-subnormal quanta.
        # Construct the trial output from bits before changing the mode.
        out = (reinterpret(T, uinttype(T)(magnitude)), zero(T), zero(T))
        for flush in (false, true) # The normal call also prewarms this typed owner.
            verdict, active = retention_in_mode(out, exact, scale2, flush)
            @test active == flush
            @test get_zero_subnormals() == saved
            @test verdict == (magnitude <= 4)
        end
    end
end

@testset "Generic Boris arithmetic dispatch" begin
    values = (
        Q(1) / 2,
        -Q(1) / 4,
        Q(1) / 8,
        Q(1) / 8,
        Q(1) / 16,
        -Q(1) / 32,
        Q(1) / 4,
        -Q(1) / 8,
        Q(1) / 16,
        Q(2) / 3,
        Q(1) / 2,
    )
    inputs = map(GenericBorisNumber, values)
    generic_binary_adds[] = 0
    generic_ternary_adds[] = 0
    out, parts = H._boris_kick_parts(inputs...)
    default_binary = generic_binary_adds[]
    @test generic_ternary_adds[] == 1 # The historical three-square sum.
    exact, _ = exact_reference(values[1:3], values[4:6], values[7:9], values[10], values[11])
    @test map(x -> x.value, out) == exact
    @test map(x -> x.value, parts) == last(H._boris_kick_parts(values...))
    generic_binary_adds[] = 0
    generic_ternary_adds[] = 0
    traced_out, traced_parts, nodes = H._boris_kick_parts(inputs..., Val(true))
    @test generic_ternary_adds[] == 1
    @test generic_binary_adds[] == default_binary + 1 # Explicit diagnostic prefix only.
    @test map(x -> x.value, traced_out) == exact
    @test map(x -> x.value, traced_parts) == map(x -> x.value, parts)
    @test nodes[15].value == nodes[12].value + nodes[13].value
    @test nodes[16].value == nodes[12].value + nodes[13].value + nodes[14].value
    generic_ternary_adds[] = 0
    @test map(x -> x.value, boris_kick(inputs...)) == exact
    @test generic_ternary_adds[] == 1
    for bits in (64, 128, 256), mode in (RoundNearest, RoundUp, RoundDown, RoundToZero)
        setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                # At this precision the two small squares expose an incorrectly
                # reassociated square sum next to the first large square.
                z, o = zero(BigFloat), one(BigFloat)
                t = (ldexp(o, bits), ldexp(o, bits ÷ 2), ldexp(o, bits ÷ 2))
                args = (o, -o, o, z, z, z, t..., BigFloat(2), o)
                actual, detail = H._boris_kick_parts(args...)
                @test detail[8] == t[1] * t[1] + t[2] * t[2] + t[3] * t[3]
                @test isequal(actual, boris_kick(args...))
            end
        end
    end
end

@testset "Original-factor and explicit-qm Boris accuracy" begin
    for T in (Float16, Float32, Float64)
        for row in stored_rows(T)
            check_case(T, row.q, row.m, row.v, row.E, row.B, row.dt)
            check_public_case(T, row)
            scales = H._boris_batch_scales(
                map(x -> [x], row.v),
                row.E,
                row.B,
                (row.q / row.m) * row.dt / 2,
                1:1,
            )
            joint, _ = H._boris_batch_admission(
                (row.q / row.m) * row.dt / 2,
                H._BorisSpeciesFactors(row.q, row.m, row.q / row.m),
                row.dt,
                scales...,
            )
            if T === Float16
                @test !joint # A range-only bypass cannot admit the retained failures.
            end
        end
        z = (zero(T), zero(T), zero(T))
        D = nextfloat(zero(T))
        for q in (one(T), -one(T), zero(T), -zero(T)), dt in (T(2), -T(2), zero(T), -zero(T))
            check_case(
                T,
                q,
                T(3),
                (T(0.25), T(-0.5), T(1)),
                (T(0.125), T(-0.25), T(0.5)),
                (T(0.125), T(0.25), T(-0.5)),
                dt,
            )
        end
        for B in ((T(2) * sqrt(floatmax(T)), zero(T), zero(T)), (T(1), T(-1), T(1)))
            check_case(T, one(T), one(T), (floatmax(T), floatmax(T), -floatmax(T)), z, B, T(2))
        end
        check_case(
            T,
            one(T),
            one(T),
            z,
            (eps(T) / T(2), -eps(T) / T(2), eps(T) / T(2)),
            (one(T) / floatmin(T), zero(T), zero(T)),
            T(2) * floatmin(T),
        )
        # Generic and explicit nonfinite arithmetic retains the native expression.
        @test isequal(
            boris_kick(T(NaN), z..., z..., z[1], z[2], one(T), T(2)),
            H._boris_kick_native(T(NaN), z..., z..., z[1], z[2], one(T), T(2)),
        )
        for n in (0, 1, 4), route in (:uniform, :gathered), stride in (false, true)
            q = D
            m =
                T === Float16 ? ldexp(one(T), 10) :
                T === Float32 ? ldexp(one(T), 80) : ldexp(one(T), 600)
            dt =
                T === Float16 ? ldexp(one(T), 10) :
                T === Float32 ? ldexp(one(T), 60) : ldexp(one(T), 400)
            for sign in (-one(T), one(T))
                ps = ParticleSet{3,T}(n; q = sign * q, m)
                foreach(x -> fill!(x, zero(T)), ps.x)
                E = (zero(T), zero(T), dt)
                fields = map(e -> fill(e, n), E)
                magnetic = ntuple(_ -> zeros(T, n), 3)
                midpoint = ntuple(_ -> stride ? view(zeros(T, 2n + 1), 2:2:2n) : zeros(T, n), 3)
                if route === :uniform
                    push_uniform!(ps, E, z, dt)
                else
                    push_gathered!(ps, fields, magnetic, dt; xmid = midpoint)
                end
                exact, scale2 = exact_reference(z, E, z, Q(sign * q) / Q(m), dt)
                @test all(p -> contract(ntuple(c -> ps.v[c][p], 3), exact, scale2), 1:n)
                expected = map(x -> nearest(T, x), exact)
                @test all(c -> all(isequal(expected[c]), ps.v[c]), 1:3)
                @test all(c -> all(isequal(nearest(T, Q(dt) * Q(expected[c]))), ps.x[c]), 1:3)
                if route === :gathered
                    @test all(
                        c -> all(isequal(nearest(T, Q(dt / T(2)) * Q(expected[c]))), midpoint[c]),
                        1:3,
                    )
                end
            end
        end
    end
    v = (1 // 2, -1 // 4, 1 // 8)
    E = (1 // 8, 1 // 16, -1 // 32)
    B = (1 // 4, -1 // 8, 1 // 16)
    @test boris_kick(v..., E..., B..., 2 // 3, 1 // 2) ==
          H._boris_kick_native(v..., E..., B..., 2 // 3, 1 // 2)
end

@testset "Model kick owners retain original species factors" begin
    for T in (Float32, Float64)
        q = nextfloat(zero(T))
        m = ldexp(one(T), T === Float32 ? 80 : 600)
        dt = ldexp(one(T), T === Float32 ? 60 : 400)
        z = (zero(T), zero(T), zero(T))
        E = (dt, zero(T), zero(T))
        @test iszero(q / m)
        exact, scale2 = exact_reference(z, E, z, Q(q) / Q(m), dt)
        expected = map(x -> nearest(T, x), exact)
        @test !iszero(expected[1])
        grid = FourierGrid((4,), (one(T),))
        for (route, primed) in (
            (:hybrid, true),
            (:hybrid, false),
            (:electrostatic, true),
            (:electrostatic, false),
            (:empic1d, true),
            (:empic, true),
        )
            ps = ParticleSet{1,T}(1; q, m)
            fill!(ps.x[1], zero(T))
            Ep = map(e -> [e], E)
            Bp = ntuple(_ -> zeros(T, 1), 3)
            start = primed ? z : map(x -> nearest(T, x), euler_reference(z, E, z, q, m, dt / 2, -1))
            reference, reference_scale = exact_reference(start, E, z, Q(q) / Q(m), dt)
            expected = map(x -> nearest(T, x), reference)
            if route === :hybrid
                st = HybridStepper(grid, HybridModel(IsothermalElectrons(zero(T))), CIC(), 1)
                fill!(st.fields.E[1], E[1])
                st.step[] = primed ? 1 : 0
                H._push_hybrid_population!(st, st, ps, dt, dt)
                half = nearest(T, Q(dt / T(2)) * Q(expected[1]))
                @test isequal(st.xmid[1][1], mod(half, grid.L[1]))
            elseif route === :electrostatic
                H._espic_push!(ps, (Ep[1],), grid, dt, primed, primed ? dt : zero(T))
            elseif route === :empic1d
                es = EMPIC1D(grid, 1; n0 = zero(T))
                H._push_velocities!(es, ps, Ep[1], Ep[2], Bp[3], dt)
            else
                es = EMPIC(grid, 1; n0 = zero(T))
                H._push_velocities!(es, ps, Ep, Bp, dt)
            end
            actual = ntuple(c -> ps.v[c][1], 3)
            @test contract(actual, reference, reference_scale)
            @test isequal(actual, expected)
            if route in (:hybrid, :electrostatic)
                moved = nearest(T, Q(dt) * Q(expected[1]))
                @test isequal(ps.x[1][1], mod(moved, grid.L[1]))
            end
        end
    end
end

@testset "Original-factor Euler and electrostatic priming range" begin
    for T in (Float16, Float32, Float64), kind in (:underflow, :overflow)
        z = (zero(T), zero(T), zero(T))
        q = kind === :underflow ? nextfloat(zero(T)) : one(T)
        m =
            kind === :underflow ? ldexp(one(T), T === Float16 ? 10 : T === Float32 ? 80 : 600) :
            nextfloat(zero(T))
        h =
            kind === :underflow ? ldexp(one(T), T === Float16 ? 10 : T === Float32 ? 60 : 400) :
            floatmin(T)
        E = (kind === :underflow ? h : one(T), kind === :underflow ? -h : -one(T), zero(T))
        B = (zero(T), zero(T), one(T))
        u = (zero(T), zero(T), T(0.25))
        factor = H._BorisSpeciesFactors(q, m, q / m)
        hq = h * factor.qm
        @test kind === :underflow ? iszero(factor.qm) : isinf(factor.qm)
        expected = map(x -> nearest(T, x), euler_reference(u, E, B, q, m, h, -1))
        v, Ep, Bp = map(x -> [x], u), map(x -> [x], E), map(x -> [x], B)
        H._prime_leapfrog!(v, Ep, Bp, factor, h, 1)
        @test isequal(map(only, v), expected)
        v = map(x -> [x], u)
        H._prime_empic1d!(v, Ep[1], Ep[2], Bp[3], factor, h, 1, false, one(T))
        @test isequal(map(only, v), expected) # Planar priming must retain stored vz.
        v = map(x -> [x], u)
        H._predict_half_kick!(v, v, Ep, Bp, factor, h, 1)
        predicted = map(x -> nearest(T, x), euler_reference(u, E, B, q, m, h, 1))
        @test isequal(map(only, v), predicted)
        @test isequal(H._es_prime_component(z[1], E[1], h, factor, hq), expected[1])
        # A vanishing Lorentz force keeps a finite unchanged velocity even
        # when the native quotient cannot represent the original factors.
        balanced = (one(T), -one(T), T(0.25))
        force = (one(T), one(T), zero(T))
        actual = H._euler_kick(balanced, force, B, h, factor, hq, Val(:prime))
        @test isequal(actual, balanced)
        zero_factor = H._BorisSpeciesFactors(-zero(T), T(2), -zero(T))
        signed = (-zero(T), zero(T), -zero(T))
        native = ntuple(c -> signed[c] - (one(T) * zero_factor.qm) * force[c], 3)
        @test isequal(
            H._euler_kick(signed, force, z, one(T), zero_factor, -zero(T), Val(:prime)),
            native,
        )
        empty = ntuple(_ -> T[], 3)
        @test H._prime_leapfrog!(empty, empty, empty, factor, h, 0) === empty
        @test H._predict_half_kick!(empty, empty, empty, empty, factor, h, 0) === empty
        huge = H._BorisSpeciesFactors(floatmax(T), nextfloat(zero(T)), T(Inf))
        @test_throws ArgumentError H._euler_kick(
            z,
            (one(T), z[2], z[3]),
            z,
            one(T),
            huge,
            T(Inf),
            Val(:prime),
        )
        @test H._euler_kick(z, (one(T), z[2], z[3]), z, one(T), huge, T(Inf), Val(:predict))[1] ==
              T(Inf)
        @test_throws ArgumentError H._es_prime_component(z[1], one(T), one(T), huge, T(Inf))
    end
end

@testset "Euler priming and predictor preserve aliased sequential reads" begin
    T, n = Float16, 14
    h = T(3) / T(1024)
    for route in (:prime, :planar, :predict)
        backing = fill(T(96), n + 1)
        v = (view(backing, 1:n), view(backing, 2:(n+1)), zeros(T, n))
        E = ntuple(_ -> zeros(T, n), 3)
        B = (zeros(T, n), zeros(T, n), v[1])
        factor = route === :predict ? -one(T) : one(T)
        sign = route === :predict ? 1 : -1
        scales = H._euler_batch_scales(v, E, B, h * factor, 1:n, T)
        @test H._euler_batch_certificate(h, factor, h * factor, scales...)
        @test !H._euler_batch_storage_safe(v, v, E, B)

        # Replay the historical particle/store order. Each scalar result is
        # checked independently against exact stored-input Euler physics before
        # it becomes the next particle's input through the shifted views.
        expected = copy(backing)
        sequential = (view(expected, 1:n), view(expected, 2:(n+1)), zeros(T, n))
        witnessed = false
        for p = 1:n
            u = ntuple(c -> sequential[c][p], 3)
            e, b = (zero(T), zero(T), zero(T)), (zero(T), zero(T), expected[p])
            exact = euler_reference(u, e, b, factor, one(T), h, sign)
            U, A, C = map(x -> Q.(collect(x)), (u, e, b))
            scale2 = max(
                sum(abs2, U),
                Q(h * factor)^2 * sum(abs2, A),
                Q(h * factor)^2 * sum(abs2, U) * sum(abs2, C),
            )
            mode = route === :predict ? Val(:predict) : Val(:prime)
            actual = H._euler_kick(u, e, b, h, factor, h * factor, mode)
            @test contract(actual, exact, scale2)
            native = first(H._euler_kick_parts(u..., e..., b..., h * factor, mode))
            if !all(isfinite, native)
                rounded = map(x -> nearest(T, x), exact)
                @test all(isfinite, rounded)
                @test isequal(actual, rounded)
                witnessed = true
            end
            sequential[1][p], sequential[2][p], sequential[3][p] = actual
        end
        @test witnessed
        if route === :prime
            H._prime_leapfrog!(v, E, B, factor, h, n)
        elseif route === :planar
            H._prime_empic1d!(v, E[1], E[2], B[3], factor, h, n, false, one(T))
        else
            H._predict_half_kick!(v, v, E, B, factor, h, n)
        end
        @test isequal(backing, expected)
        @test isequal(v[3], sequential[3])
    end
    independent = ntuple(_ -> zeros(T, 2), 3)
    fields = ntuple(_ -> zeros(T, 2), 3)
    @test H._euler_batch_storage_safe(independent, independent, fields, fields)
    repeated = view(T[1], [1, 1])
    @test !H._euler_batch_storage_safe((repeated,), (repeated,), (), ())
end

@testset "Original-factor hybrid phase and checkpoint continuation" begin
    grid = FourierGrid((4,), (1.0,))
    ps = ParticleSet{1,Float64}(1; m = 2.0)
    fill!(ps.x[1], 0.25)
    fill!(ps.weight, 0.0)
    st = HybridStepper(grid, HybridModel(IsothermalElectrons(0.0)), CIC(), 1)
    init!(st, ps)
    step!(st, ps, 0.0)
    @test (st.step[], st.last_dt[], st.time[]) == (0, 0.0, 0.0)
    elapsed = 0.0
    for dt in (0.125, 0.0625)
        fill!(st.fields.E[1], 1.0)
        step!(st, ps, dt)
        elapsed += dt
        @test ps.v[1][1] == 0.5 * (elapsed - dt / 2)
        @test ps.x[1][1] == 0.25 + 0.25 * elapsed^2
        @test st.last_dt[] == dt
    end
    mktempdir() do dir
        path = joinpath(dir, "original-factor-phase.jld2")
        save_checkpoint(path, st, ps)
        restored = HybridStepper(grid, st.model, CIC(), 1)
        particles = ParticleSet{1,Float64}(1)
        load_checkpoint!(restored, particles, path)
        @test (particles.q, particles.m, restored.last_dt[]) == (ps.q, ps.m, st.last_dt[])
        for (solver, population) in ((st, ps), (restored, particles))
            fill!(solver.fields.E[1], 1.0)
            step!(solver, population, 0.03125)
        end
        @test particles.x == ps.x && particles.v == ps.v
        @test (restored.step[], restored.time[], restored.last_dt[]) ==
              (st.step[], st.time[], st.last_dt[])
    end
end

@testset "Euler and ES nearest-overflow endpoints" begin
    for T in (Float16, Float32, Float64), kind in (:finite_native, :overflowing_product)
        ulp = Q(2)^(exponent(floatmax(T)) - precision(T) + 1)
        threshold = Q(floatmax(T)) + ulp / 2
        if kind === :finite_native
            power = Q(2)^exponent(floatmax(T))
            # Retained coordinate fixture: 3 divides 2^p+1 for p=11/53;
            # 97 divides 2^24+1. Euler prediction is the same x+hE expression.
            h = T === Float32 ? T(97) / T(64) : T(3) / T(2)
            x, electric = T(power - ulp), (power + ulp / 2) / Q(h)
        else
            # The exact product is the overflow tie; ±one subnormal in u
            # distinguishes it from its finite and overflowing neighbors.
            h = T === Float32 ? T(31) : T(3)
            x, electric = zero(T), threshold / Q(h)
        end
        E0 = T(electric)
        @test Q(E0) == electric
        z = zero(T)
        for direction in (-one(T), one(T)), offset in (-1, 0, 1)
            original = offset < 0 ? prevfloat(x) : offset > 0 ? nextfloat(x) : x
            u, E = direction * original, direction * E0
            exact = Q(u) + Q(h) * Q(E)
            scale2 = max(Q(u)^2, (Q(h) * Q(E))^2)
            native = u + (h * one(T)) * E
            kind === :finite_native && offset == 0 && @test isfinite(native)
            for factor in (one(T), H._BorisSpeciesFactors(one(T), one(T), one(T)))
                predicted =
                    H._euler_kick((u, z, z), (E, z, z), (z, z, z), h, factor, h, Val(:predict))
                @test contract(predicted, (exact, Q(0), Q(0)), scale2)
                if abs(exact) < threshold
                    expected = isfinite(native) ? native : nearest(T, exact)
                    @test isequal(predicted[1], expected)
                    primed =
                        H._euler_kick((u, z, z), (-E, z, z), (z, z, z), h, factor, h, Val(:prime))
                    @test isequal(primed[1], expected)
                    @test isequal(H._es_prime_component(u, -E, h, factor, h), expected)
                    planar = ([u], [z], [T(0.25)])
                    H._prime_empic1d!(planar, [-E], [z], [z], factor, h, 1, false, one(T))
                    @test isequal(only(planar[1]), expected)
                    @test only(planar[3]) == T(0.25)
                else
                    @test_throws ArgumentError H._euler_kick(
                        (u, z, z),
                        (-E, z, z),
                        (z, z, z),
                        h,
                        factor,
                        h,
                        Val(:prime),
                    )
                    @test_throws ArgumentError H._es_prime_component(u, -E, h, factor, h)
                    planar = ([u], [z], [T(0.25)])
                    @test_throws ArgumentError H._prime_empic1d!(
                        planar,
                        [-E],
                        [z],
                        [z],
                        factor,
                        h,
                        1,
                        false,
                        one(T),
                    )
                end
            end
        end
    end
end

@testset "Euler ordinary accuracy admission retains native arithmetic" begin
    for T in (Float16, Float32, Float64)
        u, E, B = (T(0.5), T(-0.25), T(0.125)),
        (T(0.125), T(0.0625), T(-0.03125)),
        (T(0.25), T(-0.125), T(0.0625))
        h = T(0.0625)
        for factor in (one(T), H._BorisSpeciesFactors(one(T), one(T), one(T)))
            native = ntuple(c -> u[c] - h * E[c], 3) # u and B are parallel; cross(u,B)=0 exactly.
            actual = H._euler_kick(u, E, B, h, factor, h, Val(:prime))
            @test isequal(actual, native)
            scales = H._euler_batch_scales(map(x -> [x], u), E, B, h, 1:1, T)
            @test H._euler_batch_certificate(h, factor, h, scales...)
            exact = euler_reference(u, E, B, one(T), one(T), h, -1)
            U, A, C = map(x -> Q.(collect(x)), (u, E, B))
            scale2 = max(dot(U, U), Q(h)^2 * dot(A, A), Q(h)^2 * dot(U, U) * dot(C, C))
            @test contract(actual, exact, scale2)
        end
    end
end

@testset "Repeated particle storage retains sequential public pushes" begin
    # Four small dyadic increments and every stored intermediate are exact in
    # Float16 as well as the wider types. Uniform push kicks the whole batch
    # before drifting; gathered push kicks and drifts each particle in order.
    for T in (Float16, Float32, Float64), gathered in (false, true)
        n, dt, amplitude = 4, T(1 / 16), T(1 / 8)
        positions = view(T[1/4], fill(1, n))
        velocity = ntuple(3) do c
            view(T[c == 1 ? 1 / 2 : 0], fill(1, n))
        end
        ps = ParticleSet{1,T}(
            (positions,),
            velocity,
            ones(T, n),
            UInt64.(1:n),
            zeros(UInt32, n),
            one(T),
            T(2),
        )
        if gathered
            electric = (fill(amplitude, n), zeros(T, n), zeros(T, n))
            magnetic = ntuple(_ -> zeros(T, n), 3)
            push_gathered!(ps, electric, magnetic, dt)
        else
            push_uniform!(ps, (amplitude, zero(T), zero(T)), ntuple(_ -> zero(T), 3), dt)
        end
        increment = dt * amplitude / T(2)
        expected_v = T(1 / 2) + T(n) * increment
        velocity_sum = T(n) * T(1 / 2) + T(n * (n + 1) ÷ 2) * increment
        expected_x = T(1 / 4) + dt * (gathered ? velocity_sum : T(n) * expected_v)
        @test all(==(expected_v), ps.v[1])
        @test all(iszero, ps.v[2]) && all(iszero, ps.v[3])
        @test all(==(expected_x), ps.x[1])
    end
end

@testset "Stored scalar factor and timestep admission across flush modes" begin
    for T in (Float16, Float32, Float64)
        U = uinttype(T)
        sign = one(U) << (8sizeof(T) - 1)
        tiny, negative_tiny = reinterpret(T, one(U)), reinterpret(T, sign | one(U))
        z, negative_zero = reinterpret(T, zero(U)), reinterpret(T, sign)
        unresolved = T[tiny, negative_tiny, floatmin(T), prevfloat(2floatmin(T))]
        invalid_mass = T[z, negative_zero, -one(T), negative_tiny, T(Inf), T(NaN)]
        oldmode = get_zero_subnormals()
        try
            @test set_zero_subnormals(false)
            ps = ParticleSet{1,T}(1; q = tiny, m = tiny)
            fill!(ps.x[1], z)
            electric = (T(0.125), z, z)
            magnetic = (z, z, z)
            gathered = map(x -> T[x], electric), map(x -> T[x], magnetic)
            midpoint = (T[z],)
            for flush in (false, true)
                @test set_zero_subnormals(flush)
                @test get_zero_subnormals() == flush
                ps.q, ps.m = tiny, tiny
                factor = H._validated_boris_factor(ps)
                observation = (
                    type = T,
                    requested_flush = flush,
                    active_flush = get_zero_subnormals(),
                    charge_word = reinterpret(U, factor.q),
                    mass_word = reinterpret(U, factor.m),
                    quotient_word = reinterpret(U, factor.qm),
                    unresolved_words = collect(reinterpret(U, unresolved)),
                )
                @info "Stored scalar admission" observation
                @test reinterpret(U, factor.q) == one(U)
                @test reinterpret(U, factor.m) == one(U)
                @test H._validated_charge_mass(ps.q, ps.m, T) === (tiny, tiny)
                @test H._stored_positive(tiny) && !H._stored_positive(negative_tiny)
                @test !H._stored_nonnegative(negative_tiny)
                for zero_charge in (z, negative_zero)
                    ps.q = zero_charge
                    @test H._validated_native_boris_qm(ps) === zero_charge
                    @test H._particle_stored_zero(zero_charge)
                    @test H._stored_nonnegative(zero_charge)
                end
                ps.q = one(T)
                for mass in invalid_mass
                    ps.m = mass
                    @test_throws ArgumentError H._validated_boris_factor(ps)
                    @test_throws ArgumentError H._validated_native_boris_qm(ps)
                end
                for (q, m) in ((tiny, floatmax(T)), (one(T), floatmax(T)), (floatmax(T), tiny))
                    ps.q, ps.m = q, m
                    @test_throws ArgumentError H._validated_native_boris_qm(ps)
                end
                ps.q, ps.m = one(T), one(T)
                @test H._validated_native_boris_qm(ps) === one(T)
                ps.q, ps.m = tiny, tiny
                for use_gathered in (false, true)
                    foreach(a -> fill!(a, z), (ps.x..., ps.v...))
                    if use_gathered
                        push_gathered!(ps, gathered..., one(T); xmid = midpoint)
                        @test only(midpoint[1]) == T(0.0625)
                    else
                        push_uniform!(ps, electric, magnetic, one(T))
                    end
                    # Independent exact ratio tiny/tiny=1, E=1/8, dt=1.
                    @test only(ps.v[1]) == T(0.125)
                    @test only(ps.x[1]) == T(0.125)
                    @test all(H._particle_stored_zero, ps.v[2])
                    @test all(H._particle_stored_zero, ps.v[3])
                end
                arrays = (ps.x..., ps.v..., midpoint...)
                before = map(a -> collect(reinterpret(U, a)), arrays)
                for dt in unresolved
                    @test_throws ArgumentError H._require_resolved_dt("stored dt", dt)
                    @test_throws ArgumentError push_uniform!(ps, electric, magnetic, dt)
                    @test_throws ArgumentError push_gathered!(ps, gathered..., dt; xmid = midpoint)
                    @test all(collect(reinterpret(U, a)) == b for (a, b) in zip(arrays, before))
                end
                @test_throws ArgumentError H._validated_nonnegative_dt(
                    T, negative_tiny; name = "stored dt",
                )
                @test H._validated_nonnegative_dt(T, tiny; name = "stored dt") === tiny
                for zero_dt in (z, negative_zero)
                    @test H._require_resolved_dt("stored dt", zero_dt) === zero_dt
                    @test push_uniform!(ps, electric, magnetic, zero_dt) === ps
                    @test all(collect(reinterpret(U, a)) == b for (a, b) in zip(arrays, before))
                end
                @test H._require_resolved_dt("stored dt", 2floatmin(T)) === 2floatmin(T)
            end
        finally
            set_zero_subnormals(oldmode)
        end
        @test get_zero_subnormals() == oldmode
    end
    @test H._particle_stored_zero(BigFloat(0))
    @test H._stored_positive(BigFloat(1)) && !H._stored_nonnegative(BigFloat(-1))
end

@testset "Stored collection charge and weight signs across flush modes" begin
    for T in (Float32, Float64)
        U = uinttype(T)
        tiny = reinterpret(T, one(U))
        negative_tiny = reinterpret(T, (one(U) << (8sizeof(T) - 1)) | one(U))
        oldmode = get_zero_subnormals()
        try
            @test set_zero_subnormals(false)
            g = FourierGrid((4,), (T(4),))
            ps = ParticleSet{1,T}(4)
            ps.x[1] .= T.(0:3)
            species = [ps]
            hybrid = HybridStepper(g, HybridModel(IsothermalElectrons(T(0.1))), CIC(), 4)
            electrostatic = ElectrostaticPIC(g, species; neutrality = :neutralize)
            electromagnetic = EMPIC(g, species; neutrality = :neutralize, c = one(T))
            init!(hybrid, species)
            init_espic!(electrostatic, species)
            init_empic!(electromagnetic, species)
            for flush in (false, true)
                @test set_zero_subnormals(flush)
                @test get_zero_subnormals() == flush
                ps.q = tiny
                @test H._validate_hybrid_species(species) == 4
                for charge in (negative_tiny, zero(T), -zero(T))
                    ps.q = charge
                    @test_throws ArgumentError H._validate_hybrid_species(species)
                end
                ps.q = one(T)
                ps.weight[2] = negative_tiny
                observation = (
                    type = T,
                    requested_flush = flush,
                    active_flush = get_zero_subnormals(),
                    positive_charge_word = reinterpret(U, tiny),
                    negative_weight_word = reinterpret(U, ps.weight[2]),
                )
                @info "Stored collection admission" observation
                arrays = (
                    hybrid.fields.B..., hybrid.fields.E..., hybrid.fields.J..., hybrid.fields.n,
                    hybrid.fields.ui..., hybrid.fields.pe,
                    electrostatic.E..., electrostatic.ne,
                    electromagnetic.E..., electromagnetic.B..., electromagnetic.J...,
                    electromagnetic.rho_n, electromagnetic.rho_np1,
                    ps.x..., ps.v..., ps.weight,
                )
                before = map(a -> collect(reinterpret(U, a)), arrays)
                phase = (
                    hybrid.time[], hybrid.step[], hybrid.last_dt[],
                    electrostatic.primed[], electrostatic.last_dt[],
                    electromagnetic.time[], electromagnetic.step[],
                    electromagnetic.last_dt[], electromagnetic.last_sub_dt[],
                )
                @test_throws ArgumentError init!(hybrid, species)
                @test_throws ArgumentError init_espic!(electrostatic, species)
                @test_throws ArgumentError init_empic!(electromagnetic, species)
                @test all(collect(reinterpret(U, a)) == b for (a, b) in zip(arrays, before))
                @test (
                    hybrid.time[], hybrid.step[], hybrid.last_dt[],
                    electrostatic.primed[], electrostatic.last_dt[],
                    electromagnetic.time[], electromagnetic.step[],
                    electromagnetic.last_dt[], electromagnetic.last_sub_dt[],
                ) === phase
                ps.weight[2] = one(T)
            end
        finally
            set_zero_subnormals(oldmode)
        end
        @test get_zero_subnormals() == oldmode
    end
end
end
