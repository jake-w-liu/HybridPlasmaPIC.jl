module OriginalIdealFaradayRangeTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function stage(::Type{T}, D, closure; n = 4, η = zero(T)) where {T}
    ek, eu = T === Float32 ? (60, -100) : (300, -600)
    g = FourierGrid(ntuple(_ -> n, D), ntuple(_ -> ldexp(T(2pi), -ek), D))
    st = HybridStepper(g, HybridModel(closure; η), NGP(), 0)
    fill!(st.fn, one(T))
    fill!(st.fields.n, one(T))
    fill!(st.fields.ninv, one(T))
    return st, ldexp(one(T), eu)
end

# Independent exact four-point centered collocation row: Df[i]=k(f[i+1]-f[i-1])/2.
function D4(field, g, axis)
    out = zeros(Q, size(field))
    for I in CartesianIndices(field)
        plus = CartesianIndex(ntuple(d -> d == axis ? mod1(I[d] + 1, 4) : I[d], ndims(field)))
        minus = CartesianIndex(ntuple(d -> d == axis ? mod1(I[d] - 1, 4) : I[d], ndims(field)))
        out[I] = Q(g.kvec[axis][2]) * (Q(field[plus]) - Q(field[minus])) / 2
    end
    return out
end

@testset "Original product survives pointwise storage before Hybrid Faraday" begin
    for T in (Float32, Float64),
        D = 1:3,
        closure in (IsothermalElectrons(zero(T)), EvolvingElectrons(T(2)), TransportElectrons(T(2)))

        st, a = stage(T, D, closure)
        fill!(st.fields.B[2], a)
        for I in CartesianIndices(st.fn)
            st.fui[1][I] = a * T((0, 1, 0, -1)[I[1]])
        end
        @test iszero(a * a)
        expected = -Q(a) .* D4(st.fui[1], st.g, 1)
        if st.pressure_work === nothing
            H._bfield_rhs!(st.k1, st.fields.B, st)
        else
            H._bpressure_rhs!(st.k1, st.pressure_work.kp[1], st.fields.B, st.fields.pe, st, 1, 1)
            @test all(iszero, st.pressure_work.kp[1])
        end
        @test all(all(iszero, e) for e in st.Escr)
        @test isequal(st.k1[2], T.(expected))
        @test any(!iszero, st.k1[2])
        @test all(iszero, st.k1[1]) && all(iszero, st.k1[3])
    end
end

@testset "Product-range source cancellation, constants and centered nulls" begin
    for T in (Float32, Float64)
        st, a = stage(T, 1, IsothermalElectrons(zero(T)))
        fill!(st.fields.B[1], one(T))
        fill!(st.fields.B[2], a)
        fill!(st.fui[2], floatmax(T) / 4)
        st.fui[1] .= a .* T[0, 1, 0, -1]
        H._bfield_rhs!(st.k1, st.fields.B, st)
        @test all(==(floatmax(T) / 4), st.Escr[3])
        @test isequal(st.k1[2], T.(-Q(a) .* D4(st.fui[1], st.g, 1)))

        # The aliased collocation derivative is not a continuum product rule:
        # sampled sin(x)cos(x) vanishes on this four-point grid.
        aliased, a = stage(T, 1, IsothermalElectrons(zero(T)))
        aliased.fui[1] .= a .* T[0, 1, 0, -1]
        aliased.fields.B[2] .= a .* T[2, 1, 0, 1]
        u, b = aliased.fui[1], aliased.fields.B[2]
        collocation = -D4(Q.(u) .* Q.(b), aliased.g, 1)
        product_rule = -Q.(b) .* D4(u, aliased.g, 1) - Q.(u) .* D4(b, aliased.g, 1)
        @test product_rule != collocation
        H._bfield_rhs!(aliased.k1, aliased.fields.B, aliased)
        @test isequal(aliased.k1[2], T.(collocation))

        for (n, pattern) in ((1, [1]), (2, [1, -1]), (4, [1, 1, 1, 1]), (4, [1, -1, 1, -1]))
            null, a = stage(T, 1, IsothermalElectrons(zero(T)); n)
            fill!(null.fields.B[2], a)
            null.fui[1] .= a .* T.(pattern)
            H._bfield_rhs!(null.k1, null.fields.B, null)
            @test all(all(iszero, rate) for rate in null.k1)
        end
        # Separate overflowing products cancel in the constitutive expression.
        canceled, _ = stage(T, 1, IsothermalElectrons(zero(T)))
        large = ldexp(one(T), exponent(floatmax(T)) ÷ 2 + 1)
        for c in (2, 3)
            fill!(canceled.fields.B[c], large)
            fill!(canceled.fui[c], large)
        end
        H._bfield_rhs!(canceled.k1, canceled.fields.B, canceled)
        @test all(all(iszero, rate) for rate in canceled.k1)
        fill!(canceled.fui[3], zero(T))
        # An unrepresentable required physical E is rejected even on a null curl axis.
        @test_throws ArgumentError H._bfield_rhs!(canceled.k1, canceled.fields.B, canceled)

        # Equal sampled Ex,Ey have equal x/y derivatives on this diagonal mode;
        # the signed modal curl cancels before inverse transformation.
        matched, a = stage(T, 2, IsothermalElectrons(zero(T)))
        fill!(matched.fields.B[3], a)
        for I in CartesianIndices(matched.fn)
            u = a * T((0, 1, 0, -1)[mod1(I[1] + I[2] - 1, 4)])
            matched.fui[1][I], matched.fui[2][I] = u, -u
        end
        H._bfield_rhs!(matched.k1, matched.fields.B, matched)
        @test all(all(iszero, rate) for rate in matched.k1)
    end
end

@testset "Ohmic product source and ordinary native neighbor" begin
    for T in (Float32, Float64)
        η = ldexp(one(T), T === Float32 ? -100 : -600)
        st, _ = stage(T, 1, IsothermalElectrons(zero(T)); η)
        amplitude = ldexp(one(T), T === Float32 ? -140 : -1000)
        st.fields.B[2] .= amplitude .* T[0, 1, 0, -1]
        H._bfield_rhs!(st.k1, st.fields.B, st)
        expected = Q(η) .* D4(st.fields.J[3], st.g, 1)
        @test all(iszero, st.Escr[3]) && any(!iszero, expected)
        @test isequal(st.k1[2], T.(expected))

        overflow, a = stage(T, 1, IsothermalElectrons(zero(T)))
        fill!(overflow.fields.B[1], one(T))
        fill!(overflow.fields.B[2], a)
        overflow.fui[1] .= a .* T[0, 1, 0, -1]
        overflow.fui[2] .= (floatmax(T) / 4) .* T[0, 1, 0, -1]
        @test_throws ArgumentError H._bfield_rhs!(overflow.k1, overflow.fields.B, overflow)

        underflow = HybridStepper(
            FourierGrid((4,), (T(2pi),)),
            HybridModel(IsothermalElectrons(zero(T))),
            NGP(),
            0,
        )
        fill!(underflow.fields.ninv, one(T))
        fill!(underflow.fields.B[2], a)
        underflow.fui[1] .= a .* T[0, 1, 0, -1]
        @test_throws ArgumentError H._bfield_rhs!(underflow.k1, underflow.fields.B, underflow)

        ordinary = HybridStepper(
            FourierGrid((4,), (T(2pi),)),
            HybridModel(IsothermalElectrons(zero(T))),
            NGP(),
            0,
        )
        fill!(ordinary.fields.ninv, one(T))
        fill!(ordinary.fields.B[2], T(0.25))
        ordinary.fui[1] .= T[0, 0.125, 0, -0.125]
        field = ntuple(_ -> zeros(T, 4), 3)
        @test !H._ohm_Efield!(
            field,
            ordinary.fui,
            ordinary.fields.B,
            ordinary.fields.J,
            ordinary.fields.lapJ,
            ordinary.fields.gradp,
            ordinary.fields.ninv,
            zero(T),
            zero(T),
            ordinary.g,
        )
        expected_rate = ntuple(_ -> zeros(T, 4), 3)
        faraday_rhs!(expected_rate, field, ordinary.g)
        H._bfield_rhs!(ordinary.k1, ordinary.fields.B, ordinary)
        @test isequal(ordinary.k1, expected_rate)
        args = (
            T(NaN),
            zero(T),
            zero(T),
            one(T),
            zero(T),
            zero(T),
            one(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            nothing,
        )
        value, risk = H._ohm_scalar_field_value(args..., Val(true))
        @test isnan(value) && !risk
    end
end
end
