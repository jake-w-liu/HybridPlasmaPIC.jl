module EMPICChargeConservationTests

using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# Independent nearest-even cells; no conversion through the production owner.
function rounding_cell(value::T, exact::Q) where {T<:Union{Float32,Float64}}
    isnan(value) && return false
    maximum = Q(floatmax(T))
    overflow = maximum + (maximum - Q(prevfloat(floatmax(T)))) / 2
    isinf(value) && return signbit(value) ? exact <= -overflow : exact >= overflow
    before, after = prevfloat(value), nextfloat(value)
    lower = isfinite(before) ? (Q(before) + Q(value)) / 2 : -overflow
    upper = isfinite(after) ? (Q(after) + Q(value)) / 2 : overflow
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, value))
    return (exact > lower || (even && exact == lower)) &&
           (exact < upper || (even && exact == upper))
end

@testset "Exact represented charge sums" begin
    rng = MersenneTwister(372915)
    for T in (Float16, Float32, Float64)
        U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
        u, m = nextfloat(zero(T)), floatmax(T)
        cases = [
            T[],
            T[0, -0.0],
            T[m],
            T[m, m, -m],
            T[m, u, -m],
            T[1, 1, -1, -1, u],
            fill(T(1.25), 257),
            fill(zero(T), 17),
        ]
        # Exercise both sides of the integer-accumulator exponent bound, as
        # well as unrelated repeated runs and arbitrary cancellation order.
        for power in (60, 90, 110), count in (1, 3, 17)
            value = ldexp(one(T), min(power, exponent(floatmax(T))))
            push!(cases, vcat(fill(value, count), T[1], fill(-value, count)))
        end
        for n in (1, 2, 3, 16, 127), trial = 1:20
            values = T[]
            while length(values) < n
                value = reinterpret(T, rand(rng, U))
                isfinite(value) && push!(values, value)
            end
            push!(cases, values)
        end
        for values in cases
            expected = sum(Q, values; init = Q(0))
            before = copy(values)
            @test H._exact_float_sum(values) == expected
            @test H._exact_float_sum(view(values, :)) == expected
            @test isequal(values, before)
        end
        for bad in (T(NaN), T(Inf), -T(Inf))
            @test_throws ArgumentError H._exact_float_sum(T[bad])
            @test_throws ArgumentError H._exact_float_sum(T[1, bad, 1])
        end
        for bits in (2, 53, 127, 1024)
            old = precision(BigFloat)
            setprecision(BigFloat, bits) do
                @test H._exact_float_sum(T[m, u, -m]) == Q(u)
                @test precision(BigFloat) == bits
            end
            @test precision(BigFloat) == old
        end
    end
end

@testset "Complete uniform charge shift" begin
    rng = MersenneTwister(203914)
    for T in (Float32, Float64)
        u, f, m = nextfloat(zero(T)), floatmin(T), floatmax(T)
        values = T[0, -0.0, u, -u, f, -f, prevfloat(f), nextfloat(f), 1, -1, m, -m, prevfloat(m)]
        for _ = 1:200
            value = reinterpret(T, rand(rng, T === Float32 ? UInt32 : UInt64))
            isfinite(value) && push!(values, value)
        end
        shifts = Q[
            0,
            Q(u)/2,
            -Q(u)/2,
            3Q(u)/2,
            Q(f)/3,
            Q(1)/3,
            -Q(1)/3,
            Q(m),
            Q(m)+Q(u),
            2Q(m),
            -2Q(m),
            Q(m)+Q(m)/2,
        ]
        for shift in shifts
            lower, upper = H._empic_charge_shift_bounds(T, shift)
            @test !isfinite(lower) || Q(lower) <= shift
            @test !isfinite(upper) || shift <= Q(upper)
            out = copy(values)
            @test H._empic_complete_charge_shift!(out, shift) === out
            for (value, input) in zip(out, values)
                @test rounding_cell(value, Q(input) - shift)
            end
        end
        out = T[m]
        @test H._empic_complete_charge_shift!(out, 2Q(m)) === out
        @test out == T[-m]
        storage = T[8, 1, 2, 4]
        cells = view(storage, 2:3)
        @test H._empic_complete_charge_shift!(cells, Q(1) / 3) === cells
        @test storage[[1, 4]] == T[8, 4]
        @test all(rounding_cell(cells[i], Q(i) - Q(1) / 3) for i in eachindex(cells))
        for bad in (T(NaN), T(Inf), -T(Inf))
            out = T[1, bad, 2]
            before = copy(out)
            @test_throws ArgumentError H._empic_complete_charge_shift!(out, Q(1))
            @test isequal(out, before)
        end
        out = T[1, 2]
        before = copy(out)
        @test_throws ArgumentError H._empic_complete_charge_shift!(out, BigInt(1) // BigInt(0))
        @test isequal(out, before)
        for bits in (2, 53, 127, 1024)
            old = precision(BigFloat)
            setprecision(BigFloat, bits) do
                out = T[1, -1, u]
                H._empic_complete_charge_shift!(out, Q(1) / 3)
                @test all(rounding_cell(x, Q(y) - Q(1) / 3) for (x, y) in zip(out, T[1, -1, u]))
                @test precision(BigFloat) == bits
            end
            @test precision(BigFloat) == old
        end
        # A rational target need not equal the exact sum of rounded cells.
        out = zeros(T, 2)
        H._empic_complete_charge_shift!(out, -Q(1) / 6)
        @test all(rounding_cell(x, Q(1) / 6) for x in out)
        @test sum(Q, out) != Q(1) / 3
    end
end

function fixture(
    T,
    widths;
    mobile = true,
    ew = T[1],
    iw = T[],
    shape = NGP(),
    n0 = zero(T),
    legacy = false,
    cells = 2,
)
    D = length(widths)
    g = FourierGrid(ntuple(_ -> cells, D), map(x -> T(cells) * x, widths))
    e = ParticleSet{D,T}(length(ew); q = -one(T), m = one(T))
    ions = ParticleSet{D,T}(length(iw); q = one(T), m = T(1836))
    for (ps, weights) in ((e, ew), (ions, iw))
        foreach(x -> fill!(x, zero(T)), ps.x)
        copyto!(ps.weight, weights)
    end
    es =
        legacy ? EMPIC1D(g, length(ew); mobile, shape, n0) : EMPIC(g, length(ew); mobile, shape, n0)
    return (; es, e, ions, g)
end

function exact_target(s)
    volume = prod(Q(dx) for dx in s.g.dx)
    return s.es.mobile ?
           (
        Q(s.e.q) * sum(Q, s.e.weight; init = Q(0)) +
        Q(s.ions.q) * sum(Q, s.ions.weight; init = Q(0))
    ) / volume : length(s.es.rho_n) * Q(s.es.n0) - sum(Q, s.e.weight; init = Q(0)) / volume
end

function check_capture(s)
    es, e, ions = s.es, s.e, s.es.mobile ? s.ions : nothing
    raw, scratch = similar(es.rho_n), similar(es.ne)
    state = H._empic_raw_charge_density!(raw, scratch, es, e, ions)
    arrays = Tuple(a for ps in (s.e, s.ions) for a in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
    before = map(copy, arrays)
    @test H._charge_density!(es.rho_n, es.ne, es, e, ions) === es.rho_n
    @test isequal(es.ne, scratch)
    @test all(isequal(a, b) for (a, b) in zip(arrays, before))
    @test state.eligible && all(isfinite, raw)
    target = exact_target(s)
    shift = (sum(Q, raw) - target) / length(raw)
    @test all(rounding_cell(es.rho_n[i], Q(raw[i]) - shift) for i in eachindex(raw))
    if all(isfinite, es.rho_n)
        errors = Q.(es.rho_n) .- (Q.(raw) .- shift)
        @test sum(Q, es.rho_n) - target == sum(errors)
        # Uniform correction cancels from every exact cell difference. The
        # remaining non-DC differences are precisely final rounding errors.
        @test all(
            Q(es.rho_n[i]) - Q(es.rho_n[1]) == Q(raw[i]) - Q(raw[1]) + errors[i] - errors[1] for
            i in eachindex(raw)
        )
    end
    return raw, shift
end

@testset "Actual signed charge and complete stored geometry" begin
    rng = MersenneTwister(63720)
    for T in (Float32, Float64),
        D = 1:3,
        legacy in (false, true),
        mobile in (false, true),
        shape in (NGP(), CIC(), TSC())

        legacy && D != 1 && continue
        s = fixture(
            T,
            ntuple(_ -> T(0.7), D);
            ew = randn(rng, T, 13),
            iw = randn(rng, T, 9),
            mobile,
            shape,
            legacy,
            n0 = T(0.3),
            cells = 4,
        )
        for ps in (s.e, s.ions), d = 1:D
            ps.x[d] .= rand(rng, T, length(ps.weight)) .* s.g.L[d]
        end
        check_capture(s)
        first_target = exact_target(s)
        s.e.weight[1] += one(T)
        @test exact_target(s) != first_target
        check_capture(s)
    end
    for T in (Float32, Float64), legacy in (false, true)
        for (ew, iw, mobile, n0, target) in (
            (T[3, -2, 0], T[2], true, T(0), Q(1)),
            (T[], T[1, -1], true, T(0), Q(0)),
            (T[], T[], false, T(0.25), Q(1) / 2),
            (T[3, -2], T[], false, T(0.25), -Q(1) / 2),
        )
            s = fixture(T, (one(T),); ew, iw, mobile, n0, legacy)
            check_capture(s)
            @test sum(Q, s.es.rho_n) == target
        end
        # The complete raw replay and final reconciliation are distinct even
        # when the large species cancel exactly before cell conversion.
        changed = false
        for surviving in (one(T), T(1.3))
            s = fixture(
                T,
                (one(T),);
                ew = fill(floatmax(T) / 4, 32),
                iw = vcat(fill(floatmax(T) / 4, 32), surviving),
                legacy,
                shape = TSC(),
                cells = 5,
            )
            fill!(s.e.x[1], T(0.1))
            fill!(s.ions.x[1], T(0.1))
            raw, shift = check_capture(s)
            changed |= raw != s.es.rho_n
            @test iszero(shift) ||
                  any(iszero(raw[i]) && !iszero(s.es.rho_n[i]) for i in eachindex(raw))
        end
        @test changed
        s = fixture(T, (one(T),); ew = T[], iw = fill(floatmax(T) / 2, 4), legacy, cells = 4)
        s.ions.x[1] .= T[0, 1, 2, 3]
        check_capture(s)
        @test exact_target(s) > Q(floatmax(T))
        @test all(==(floatmax(T) / 2), s.es.rho_n)
    end
    for T in (Float32, Float64)
        recovered = T === Float32 ? (T(1e-22), T(1e-22), T(1e30)) : (T(1e-160), T(1e-160), T(1e200))
        subnormal =
            T === Float32 ? (ldexp(one(T), -70), ldexp(T(3), -80)) :
            (ldexp(one(T), -500), ldexp(T(3), -575))
        huge, tiny = T === Float32 ? (T(1e25), T(1e-25)) : (T(1e200), T(1e-200))
        cases = (
            (recovered, one(T)),
            (subnormal, floatmin(T)),
            ((huge, huge), T === Float32 ? T(1e30) : T(1e300)),
            ((tiny, tiny), T === Float32 ? T(1e-30) : T(1e-300)),
        )
        for (widths, weight) in cases
            s = fixture(T, widths; ew = T[weight])
            @test H._moment_volume(s.g).lost
            raw, _ = check_capture(s)
            @test maximum(abs, vec(s.es.rho_n)[2:end]) <= eps(abs(raw[1]))
        end
    end
end

end
