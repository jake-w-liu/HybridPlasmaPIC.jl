module CoordinateRangeEndpointTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function endpoint_inputs(::Type{T}, kind) where {T}
    ulp = Q(2)^(exponent(floatmax(T)) - precision(T) + 1)
    threshold = Q(floatmax(T)) + ulp / 2
    if kind === :finite_native
        power = Q(2)^exponent(floatmax(T))
        # 3 divides 2^p+1 for p=11/53; 97 divides 2^24+1.
        dt = T === Float32 ? T(97) / T(64) : T(3) / T(2)
        x, exact_v = T(power - ulp), (power + ulp / 2) / Q(dt)
    else
        # 3 divides 2^12-1 and 2^54-1; 31 divides 2^25-1.
        dt = T === Float32 ? T(31) : T(3)
        x, exact_v = zero(T), threshold / Q(dt)
    end
    v = T(exact_v)
    @test Q(v) == exact_v
    return x, v, dt, threshold
end

function check_range(actual::T, x, v, dt, threshold) where {T}
    exact = Q(x) + Q(dt) * Q(v)
    if abs(exact) >= threshold
        @test actual == (exact < 0 ? -T(Inf) : T(Inf))
    else
        @test isfinite(actual)
        native = H._particle_coordinate_native(x, v, dt)
        if isfinite(native)
            @test isequal(actual, native)
        else
            # These cold endpoint rows lie inside the floatmax rounding cell.
            ulp = Q(2)^(exponent(floatmax(T)) - precision(T) + 1)
            @test Q(floatmax(T)) - ulp / 2 < abs(exact) < threshold
            @test actual == (exact < 0 ? -floatmax(T) : floatmax(T))
        end
    end
end

function population(::Type{T}, n, strided, x, v) where {T}
    vector(a) = strided ? view(fill(a, 2n), 1:2:2n) : fill(a, n)
    return ParticleSet{1,T}(
        (vector(x),),
        (vector(v), vector(zero(T)), vector(zero(T))),
        ones(T, n),
        UInt64.(1:n),
        zeros(UInt32, n),
        zero(T),
        one(T),
    )
end

@testset "Coordinate nearest-overflow classification and finite native bits" begin
    for T in (Float16, Float32, Float64), kind in (:finite_native, :overflowing_product)
        x, v, dt, threshold = endpoint_inputs(T, kind)
        fields = ntuple(_ -> zero(T), 3)
        for direction in (-1, 1), offset in (-1, 0, 1)
            x0 = offset < 0 ? prevfloat(x) : offset > 0 ? nextfloat(x) : x
            xx, vv = T(direction) * x0, T(direction) * v
            check_range(H._particle_coordinate(xx, vv, dt), xx, vv, dt, threshold)
            for n in (1, 17, 64), strided in (false, true)
                # Every index, including the last, must keep the guarded path
                # when one boundary row is mixed with ordinary particles.
                for index in unique((1, cld(n, 2), n))
                    ps = population(T, n, strided, one(T), one(T))
                    ps.x[1][index], ps.v[1][index] = xx, vv
                    push_uniform!(ps, fields, fields, dt)
                    check_range(ps.x[1][index], xx, vv, dt, threshold)
                    @test ps.v[1][index] == vv

                    ps = population(T, n, strided, one(T), one(T))
                    ps.x[1][index], ps.v[1][index] = xx, vv
                    electric = ntuple(_ -> zeros(T, n), 3)
                    magnetic = ntuple(_ -> zeros(T, n), 3)
                    midpoint = (zeros(T, n),)
                    push_gathered!(ps, electric, magnetic, 2dt; xmid = midpoint)
                    check_range(midpoint[1][index], xx, vv, dt, threshold)
                    @test ps.v[1][index] == vv
                end
            end
        end
    end
end

@testset "Coordinate signed zero and nonfinite input propagation" begin
    for T in (Float16, Float32, Float64)
        for x in (zero(T), -zero(T)), v in (zero(T), -zero(T)), dt in (one(T), -one(T))
            @test isequal(H._particle_coordinate(x, v, dt), x + dt * v)
        end
        for operands in (
            (T(Inf), one(T), one(T)),
            (T(Inf), -T(Inf), one(T)),
            (one(T), T(Inf), zero(T)),
            (T(NaN), one(T), one(T)),
        )
            @test isequal(
                H._particle_coordinate(operands...),
                H._particle_coordinate_native(operands...),
            )
        end
    end
end

@testset "Shifted particle views cannot retain a stale coordinate bound" begin
    T = Float16
    for gathered in (false, true)
        # Earlier position writes become later velocity reads: 600 -> 2400 ->
        # 7280 -> 21840, exceeding the initial bound 9600. The last exact drift
        # is the overflow tie minus one minsubnormal and rounds to floatmax.
        backing = T[600, 600, 80, 0, -nextfloat(zero(T))]
        ps = ParticleSet{1,T}(
            (view(backing, 2:5),),
            (view(backing, 1:4), view(zeros(T, 4), 1:4), view(zeros(T, 4), 1:4)),
            ones(T, 4),
            UInt64.(1:4),
            zeros(UInt32, 4),
            zero(T),
            one(T),
        )
        if gathered
            E, B = ntuple(_ -> zeros(T, 4), 3), ntuple(_ -> zeros(T, 4), 3)
            push_gathered!(ps, E, B, T(3); xmid = (zeros(T, 4),))
        else
            fields = ntuple(_ -> zero(T), 3)
            push_uniform!(ps, fields, fields, T(3))
        end
        @test backing[2:4] == T[2400, 7280, 21840]
        @test backing[5] == floatmax(T)
    end
end

# Independent endpoint oracle: exact overflow threshold, native bits for a
# valid finite native result, and the known floatmax cell for rescued rows.
function endpoint_expected(::Type{T}, x, v, dt, threshold) where {T}
    exact = Q(x) + Q(dt) * Q(v)
    abs(exact) >= threshold && return exact < 0 ? -T(Inf) : T(Inf)
    native = x + dt * v
    isfinite(native) && return native
    ulp = Q(2)^(exponent(floatmax(T)) - precision(T) + 1)
    @test Q(floatmax(T)) - ulp / 2 < abs(exact) < threshold
    return exact < 0 ? -floatmax(T) : floatmax(T)
end

@testset "Coordinate endpoint classification reaches model drift owners" begin
    for T in (Float16, Float32, Float64), kind in (:finite_native, :overflowing_product)
        x, v, dt, threshold = endpoint_inputs(T, kind)
        for direction in (-1, 1), offset in (-1, 0, 1), D = 1:3, strided in (false, true)
            x0 = offset < 0 ? prevfloat(x) : offset > 0 ? nextfloat(x) : x
            xx, vv = T(direction) * x0, T(direction) * v
            vector(a) = strided ? view(fill(a, 34), 1:2:34) : fill(a, 17)
            function make_population()
                ps = ParticleSet{D,T}(
                    ntuple(_ -> vector(one(T)), D),
                    ntuple(_ -> vector(zero(T)), 3),
                    ones(T, 17),
                    UInt64.(1:17),
                    zeros(UInt32, 17),
                    zero(T),
                    one(T),
                )
                ps.x[D][17], ps.v[D][17] = xx, vv
                return ps
            end
            # The periodic midpoint owner records a native finite half-step;
            # the unwrapped full position retains the exact overflow decision.
            ps = make_population()
            midpoint = ntuple(_ -> zeros(T, 17), D)
            L = ntuple(_ -> T(8), D)
            D == 1 ? H._drift_record_mid!(ps, midpoint[1], dt, L[1]) :
            H._drift_record_mid!(ps, midpoint, dt, L)
            check_range(ps.x[D][17], xx, vv, dt, threshold)
            half = endpoint_expected(T, xx, vv, dt / 2, threshold)
            @test isfinite(half)
            @test isequal(midpoint[D][17], mod(half, L[D]))
            @test all(==(one(T)), ps.x[D][1:16])

            # The hybrid kernel publishes an unwrapped midpoint before its
            # existing boundary owner runs. Exercise that exact common owner.
            ps = make_population()
            work = (
                Ep = ntuple(_ -> zeros(T, 17), 3),
                Bp = ntuple(_ -> zeros(T, 17), 3),
                xmid = ntuple(_ -> zeros(T, 17), D),
            )
            @test !H._coordinate_batch_certified(ps, H._coordinate_batch_max(ps.x, T), dt, abs(vv))
            H._hybrid_kick_drift!(ps, work, zero(T), dt, dt, Val(true))
            check_range(ps.x[D][17], xx, vv, dt, threshold)
            @test isequal(work.xmid[D][17], half)
            @test ps.v[D][17] == vv

            # FFTW-backed electrostatics supports Float32/64. Its existing
            # periodic owner rejects true overflow and wraps finite rescue.
            if T !== Float16
                ps = make_population()
                grid = FourierGrid(ntuple(_ -> 4, D), L)
                electric = ntuple(_ -> zeros(T, 17), 3)
                expected = endpoint_expected(T, xx, vv, dt, threshold)
                if isfinite(expected)
                    H._espic_push!(ps, electric, grid, dt, true, dt)
                    @test isequal(ps.x[D][17], mod(expected, L[D]))
                    @test ps.v[D][17] == vv
                else
                    @test_throws ArgumentError H._espic_push!(ps, electric, grid, dt, true, dt)
                end
            end
        end
    end
end

@testset "Zero-factor drift keeps exact range and native IEEE bits" begin
    for T in (Float16, Float32, Float64)
        M = floatmax(T)
        coordinates = (
            zero(T),
            -zero(T),
            one(T),
            -one(T),
            T(0.75) * M,
            -T(0.75) * M,
            prevfloat(M),
            -prevfloat(M),
            M,
            -M,
        )
        factors = (-M, -one(T), -nextfloat(zero(T)), nextfloat(zero(T)), one(T), M)
        for x in coordinates,
            z in (zero(T), -zero(T)),
            factor in factors,
            (v, dt) in ((z, factor), (factor, z))

            native = H._particle_coordinate_native(x, v, dt)
            actual = H._particle_coordinate(x, v, dt)
            @test isfinite(actual)
            @test Q(actual) == Q(x)
            @test bitstring(actual) == bitstring(native)
        end
        for x in (M, -M, T(0.75) * M, -T(0.75) * M),
            n in (1, 17),
            (v, dt) in ((zero(T), M), (M, zero(T)))

            fields = ntuple(_ -> zero(T), 3)
            ps = population(T, n, false, x, v)
            push_uniform!(ps, fields, fields, dt)
            @test all(y -> bitstring(y) == bitstring(x), ps.x[1])
            @test all(y -> bitstring(y) == bitstring(v), ps.v[1])
            ps = population(T, n, false, x, v)
            gathered = ntuple(_ -> zeros(T, n), 3)
            midpoint = (zeros(T, n),)
            push_gathered!(ps, gathered, gathered, dt; xmid = midpoint)
            @test all(y -> bitstring(y) == bitstring(x), ps.x[1])
            @test all(y -> bitstring(y) == bitstring(x), midpoint[1])
            @test all(y -> bitstring(y) == bitstring(v), ps.v[1])
        end
        for x in (zero(T), M, T(Inf), -T(Inf), T(NaN)),
            z in (zero(T), -zero(T)),
            factor in (T(Inf), -T(Inf), T(NaN)),
            (v, dt) in ((z, factor), (factor, z))

            @test bitstring(H._particle_coordinate(x, v, dt)) ==
                  bitstring(H._particle_coordinate_native(x, v, dt))
        end
    end
end
end
