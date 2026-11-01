module TestRK4OriginalDuration
using HybridPlasmaPIC, Test
include("yee_charge_oracle.jl")
const H = HybridPlasmaPIC
const O = YeeChargeOracle
const Q = Rational{BigInt}

# Complete original operands; the rounding-cell checker is independent of the
# production exact conversion and does not form a rounded RK coefficient.
function affine_reference(initial, h, rates, weights, divisor)
    total = sum(weight * Q(rate) for (rate, weight) in zip(rates, weights))
    return Q(initial) + Q(h) * total / divisor
end

@testset "RK4 original duration survives coefficient range loss" begin
    for T in (Float32, Float64), quanta in (1, 3)
        h = T(quanta) * nextfloat(zero(T))
        rate = ldexp(one(T), T === Float32 ? 100 : 1000)
        @test isfinite(6 * rate)
        for divisor in (1, 2)
            measured = H._RK4Affine(h, divisor)(zero(T), rate)
            @test O.rounded(measured, affine_reference(zero(T), h, (rate,), (1,), divisor))
            @test measured > zero(T)
        end
        measured = H._RK4Affine(h, 6)(zero(T), rate, rate, rate, rate)
        @test O.rounded(measured, Q(h) * Q(rate))
        @test measured > zero(T)
    end
end

@testset "RK4 complete affine products and signed weighted sums" begin
    for T in (Float32, Float64)
        largest, quantum = floatmax(T), nextfloat(zero(T))
        h, rate = T(0.1), T(0.2)
        initial = -(h * rate)
        @test affine_reference(initial, h, (rate,), (1,), 1) != 0
        for (x, duration, k, divisor) in (
            (-largest, T(2), largest, 1),
            (quantum, T(0.5), quantum, 1),
            (initial, h, rate, 1),
        )
            measured = H._RK4Affine(duration, divisor)(x, k)
            @test O.rounded(measured, affine_reference(x, duration, (k,), (1,), divisor))
        end
        large = ldexp(one(T), T === Float32 ? 100 : 1000)
        for (x, duration, rates) in (
            (zero(T), one(T), (largest, largest, -largest, -largest)),
            (zero(T), one(T), (largest, largest / 2, -largest / 2, -prevfloat(largest))),
            (-largest, T(12), (largest, zero(T), zero(T), zero(T))),
            (zero(T), largest, (quantum, quantum, quantum, quantum)),
            (zero(T), one(T), (large, one(T), -large / 2, zero(T))),
        )
            measured = H._RK4Affine(duration, 6)(x, rates...)
            @test O.rounded(measured, affine_reference(x, duration, rates, (1, 2, 2, 1), 6))
        end
        @test H._RK4Affine(T(2), 1)(largest, largest) === T(Inf)
        @test H._RK4Affine(T(12), 6)(largest, largest, largest, largest, largest) === T(Inf)
        # Nonfinite inputs retain the original IEEE expression; no rational
        # conversion attempts to reinterpret an invalid physical rate.
        @test isnan(H._RK4Affine(zero(T), 2)(zero(T), T(Inf)))
        @test isnan(H._RK4Affine(one(T), 6)(zero(T), T(Inf), -T(Inf), zero(T), zero(T)))
    end
end

@testset "Ordinary RK4 scalar and fused-array expressions retain their bits" begin
    for T in (Float32, Float64)
        h = T(0.125)
        initial = T[zero(T), -zero(T), 1, -2, 0.125, nextfloat(zero(T))]
        rates = (
            T[zero(T), -zero(T), 0.2, -0.2, 0.3, zero(T)],
            T[zero(T), -zero(T), 0.3, -0.3, 0.4, zero(T)],
            T[zero(T), -zero(T), 0.4, -0.4, 0.5, zero(T)],
            T[zero(T), -zero(T), 0.5, -0.5, 0.6, zero(T)],
        )
        for divisor in (1, 2), k in rates
            operation = H._RK4Affine(h, divisor)
            native = @. initial + (h / divisor) * k
            @test isequal(operation.(initial, k), native)
            @test all(operation(x, r) === y for (x, r, y) in zip(initial, k, native))
        end
        k1, k2, k3, k4 = rates
        operation = H._RK4Affine(h, 6)
        native = @. initial + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        @test isequal(operation.(initial, k1, k2, k3, k4), native)
        @test all(
            operation(initial[i], k1[i], k2[i], k3[i], k4[i]) === native[i] for
            i in eachindex(initial)
        )
    end
end

@testset "Public Ohmic RK4 heating completes overflowing finite rate sum" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2pi),))
        ps = ParticleSet{1,T}(4)
        for i = 1:4
            ps.x[1][i] = T(i - 1) * g.dx[1]
        end
        set_density_weight!(ps, one(T), g)
        st = HybridStepper(g, HybridModel(TransportElectrons(T(2)); η = one(T)), NGP(), 4)
        amplitude = ldexp(one(T), T === Float32 ? 63 : 511)
        st.fields.B[2] .= amplitude .* T[1, 0, -1, 0]
        st.fields.B[3] .= amplitude .* T[0, 1, 0, -1]
        fill!(st.fields.pe, zero(T))
        init!(st, ps)
        @test all(a -> all(iszero, a), st.pressure_work.Eforce)
        initialB, initialp = map(copy, st.fields.B), copy(st.fields.pe)
        h = T(2) * floatmin(T)
        @test step!(st, ps, h) === st
        @test st.time[] === h && st.step[] == 1
        @test all(a -> all(isfinite, a), st.pressure_work.kp)
        k1, k2, k3, k4 = st.pressure_work.kp
        native_sum = @. k1 + 2 * k2 + 2 * k3 + k4
        @test all(!isfinite, native_sum)
        for i in eachindex(initialp)
            reference = affine_reference(
                initialp[i], h, (k1[i], k2[i], k3[i], k4[i]), (1, 2, 2, 1), 6,
            )
            @test O.rounded(st.fields.pe[i], reference)
        end
        # Independent physical Beltrami rate: J=-k B, gamma-1=eta=1.
        physical = Q(h) * Q(g.kvec[1][2])^2 * Q(amplitude)^2
        @test all(p -> abs(Q(p) - physical) <= 64 * Q(eps(T)) * abs(physical), st.fields.pe)
        @test all(>(zero(T)), st.fields.pe)
        @test isequal(st.fields.B, initialB)
    end
end

end
