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
        for route in (:hybrid, :electrostatic, :empic1d, :empic)
            ps = ParticleSet{1,T}(1; q, m)
            fill!(ps.x[1], zero(T))
            Ep = map(e -> [e], E)
            Bp = ntuple(_ -> zeros(T, 1), 3)
            if route === :hybrid
                st = HybridStepper(grid, HybridModel(IsothermalElectrons(zero(T))), CIC(), 1)
                fill!(st.fields.E[1], E[1])
                st.step[] = 1 # Already staggered: priming is a separate arithmetic owner.
                H._push_hybrid_population!(st, st, ps, dt, dt)
                half = nearest(T, Q(dt / T(2)) * Q(expected[1]))
                @test isequal(st.xmid[1][1], mod(half, grid.L[1]))
            elseif route === :electrostatic
                H._espic_push!(ps, (Ep[1],), grid, dt, true, dt)
            elseif route === :empic1d
                es = EMPIC1D(grid, 1; n0 = zero(T))
                H._push_velocities!(es, ps, Ep[1], Ep[2], Bp[3], dt)
            else
                es = EMPIC(grid, 1; n0 = zero(T))
                H._push_velocities!(es, ps, Ep, Bp, dt)
            end
            actual = ntuple(c -> ps.v[c][1], 3)
            @test contract(actual, exact, scale2)
            @test isequal(actual, expected)
            if route in (:hybrid, :electrostatic)
                moved = nearest(T, Q(dt) * Q(expected[1]))
                @test isequal(ps.x[1][1], mod(moved, grid.L[1]))
            end
        end
    end
end
end
