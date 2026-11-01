module ScalarPressureRangeTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function centered_four(g, axis = 1)
    k = Q(g.kvec[axis][2])
    return [mod(j - i, 4) == 1 ? k / 2 : mod(j - i, 4) == 3 ? -k / 2 : zero(k) for i = 1:4, j = 1:4]
end

function fixture(::Type{T}, kind) where {T}
    large = T === Float32 ? 80 : 600
    ek =
        kind in (:ohmic_cancel, :overflow) ? large ÷ 2 :
        kind === :derivative_rescue ? (T === Float32 ? -100 : -600) : 0
    g = FourierGrid((4,), (ldexp(T(2pi), -ek),))
    @test g.kvec[1][2] == ldexp(one(T), ek)
    p, u, b, η = if kind in (:ohmic_cancel, :overflow)
        (
            ldexp.(T[1, 2, 4, 8], large),
            ldexp.(T[-4, 34, -16, 46], large - ek),
            ldexp.(T[1, 0, -1, 0], large - ek),
            T(kind === :overflow ? 54 : 27),
        )
    elseif kind in (:half_quantum, :half_quantum_mixed)
        e = exponent(floatmin(T)) - precision(T) + 1
        j = fld(e, 2)
        (
            fill(ldexp(one(T), j), 4),
            ldexp.(T[0, 1, 0, -1], e - 2 - j),
            ldexp.(T[0, 1, 0, -1], j),
            ldexp(one(T), e - 1 - 2j),
        )
    elseif kind === :derivative_rescue
        (ldexp.(T[1, 2, 4, 8], ek), fill(ldexp(one(T), -ek), 4), zeros(T, 4), zero(T))
    elseif kind === :ideal_cancel
        (ldexp.(T[1, 2, 1, 4], large), ldexp.(T[2, 2, 2, 0], large), zeros(T, 4), zero(T))
    elseif kind === :all_underflow
        (
            fill(ldexp(one(T), T === Float32 ? -100 : -800), 4),
            ldexp.(T[0, 1, 0, -1], T === Float32 ? -60 : -400),
            zeros(T, 4),
            zero(T),
        )
    elseif kind in (:ordinary, :ideal_ordinary)
        (
            T[1, 2, 4, 8] ./ T(16),
            T[-4, 34, -16, 46] ./ T(256),
            T[1, 0, -1, 0] ./ T(8),
            kind === :ideal_ordinary ? zero(T) : T(0.125),
        )
    else
        error("unknown pressure fixture: $kind")
    end
    closure =
        kind in (:ideal_cancel, :ideal_ordinary, :derivative_rescue) ? EvolvingElectrons(T(2)) :
        TransportElectrons(T(2))
    st = HybridStepper(g, HybridModel(closure; η), CIC(), 0)
    fill!(st.fields.n, one(T))
    fill!(st.fn, one(T))
    st.fields.pe .= p
    st.fields.ui[1] .= u
    st.fui[1] .= u
    st.fields.B[2] .= b
    kind === :half_quantum_mixed && (st.fields.B[3] .= T[1, 0, -1, 0])
    ohms_law!(st.fields, st.model, g)
    return st
end

function exact_pressure(st)
    p, u = Q.(st.fields.pe), Q.(st.fui[1])
    d = centered_four(st.g)
    square = sum(component -> Q.(component) .^ 2, st.fields.J)
    return -u .* (d * p) .- 2 * p .* (d * u) .+ Q(st.model.η) .* square
end

function native_pressure!(st)
    f, w = st.fields, st.pressure_work
    electron_velocity!(w.ue, st.fui, f.J, st.fn)
    divergence!(w.divue, w.ue, st.g)
    rate = similar(f.pe)
    H._electron_pressure_rhs!(rate, f.pe, w.ue, f.gradp, w.divue, eltype(rate)(2), Val(1))
    iszero(st.model.η) ||
        H._add_ohmic_pressure_heating!(rate, f.J, eltype(rate)(2), eltype(rate)(st.model.η))
    return rate
end

function faraday_observation(st)
    f, g = st.fields, st.g
    T = eltype(f.pe)
    η = T(st.model.η)
    gradp = map(copy, f.gradp)
    J, lapJ = map(copy, f.J), map(copy, f.lapJ)
    gradient!(gradp, f.pe, g)
    field = map(similar, f.B)
    field_risk = H._ohm_Efield!(
        field, st.fui, f.B, J, lapJ, gradp, f.ninv, zero(T), zero(T), g,
    )
    ohmic_risk = H._ohm_ohmic_source_risk(J, η)
    # Capture before any other owner can reuse gradient scratch.
    source = H._ohm_original_source(st.fui, f.B, J, gradp, f.ninv, η, g)
    H._add_ohmic_field!(field, field, J, η)
    # These existing fixtures are D1, so this is curl!'s actual one-axis gain.
    gain = maximum(abs, g.kvec[1])
    scale = H.SpectralOperators._fft_input_scale(field, 3, g, gain)
    native, dyadic, without_ex = ntuple(_ -> map(similar, f.B), 3)
    H.faraday_rhs!(native, field, g)
    H._ohm_original_faraday_rhs!(dyadic, source, field, g)
    field_without_ex = map(copy, field)
    fill!(field_without_ex[1], zero(T))
    scale_without_ex = H.SpectralOperators._fft_input_scale(field_without_ex, 3, g, gain)
    H.faraday_rhs!(without_ex, field_without_ex, g)
    return (
        canonical_replay = (
            field_risk = field_risk,
            ohmic_risk = ohmic_risk,
            nearest_gradual = H._boris_accuracy_mode(),
            scale = scale,
            scale_without_ex = scale_without_ex,
        ),
        input = (
            pressure = copy(f.pe),
            ui = map(copy, st.fui),
            B = map(copy, f.B),
            J = J,
            lapJ = lapJ,
            gradp = gradp,
            ninv = copy(f.ninv),
            η = η,
            kvec = copy(g.kvec[1]),
            ik = copy(g.ik[1]),
        ),
        physical_E = field,
        native = native,
        dyadic = dyadic,
        without_ex = without_ex,
    )
end

@testset "Complete ideal and Ohmic pressure before range classification" begin
    for T in (Float32, Float64),
        kind in
        (:ohmic_cancel, :half_quantum, :half_quantum_mixed, :ideal_cancel, :derivative_rescue)

        st = fixture(T, kind)
        exact = exact_pressure(st)
        legacy = native_pressure!(st)
        if kind in (:half_quantum, :half_quantum_mixed)
            @test exact[1] == 0 && exact[3] == Q(nextfloat(zero(T)))
            @test iszero(legacy[3])
            kind === :half_quantum ? (@test all(iszero, legacy)) :
            (@test maximum(abs, legacy) >= floatmin(T))
        elseif kind === :derivative_rescue
            @test all(iszero, st.fields.gradp[1]) && all(iszero, legacy)
            @test maximum(abs, exact) > Q(floatmin(T))
        else
            @test all(iszero, exact)
            @test any(!isfinite, legacy)
        end
        observation = faraday_observation(st)
        H._bpressure_rhs!(st.k1, st.pressure_work.kp[1], st.fields.B, st.fields.pe, st, 1, 1)
        @test st.pressure_work.kp[1] == T.(exact)
        @test all(a -> all(isfinite, a), (st.fields.E..., st.Escr..., st.fields.J..., st.k1...))
        # Independent original-product Faraday reference. Longitudinal pressure
        # force is curl-null; Bx=Jx=0 for these four-point fixtures.
        ez = -Q.(st.fui[1]) .* Q.(st.fields.B[2]) .+ Q(st.model.η) .* Q.(st.fields.J[3])
        ey = Q.(st.fui[1]) .* Q.(st.fields.B[3]) .+ Q(st.model.η) .* Q.(st.fields.J[2])
        print("PRESSURE_FARADAY_OBSERVATION ")
        show(IOContext(stdout, :limit => false), (
            julia = VERSION,
            type = T,
            kind = kind,
            observation = observation,
            stage = st.k1,
            oracle_y = T.(centered_four(st.g) * ez),
            oracle_z = T.(-centered_four(st.g) * ey),
        ))
        println()
        @test st.k1[2] == T.(centered_four(st.g) * ez)
        @test st.k1[3] == T.(-centered_four(st.g) * ey)
        @test all(iszero, st.k1[1])
    end
end

@testset "Pressure completion retains admitted bits and true failures" begin
    for T in (Float32, Float64), ordinary in (:ordinary, :ideal_ordinary)
        st = fixture(T, ordinary)
        rate = native_pressure!(st)
        before = copy(rate)
        f, w = st.fields, st.pressure_work
        @test !H._complete_scalar_pressure_if_needed!(
            rate,
            f.pe,
            w.ue,
            f.gradp,
            w.divue,
            f.J,
            T(2),
            T(st.model.η),
            true,
            st.g,
        )
        @test isequal(rate, before)
        fill!(f.pe, one(T))
        foreach(a -> fill!(a, zero(T)), w.ue)
        fill!(w.divue, zero(T))
        foreach(a -> fill!(a, zero(T)), f.gradp)
        foreach(a -> fill!(a, zero(T)), f.J)
        fill!(rate, -zero(T))
        @test !H._complete_scalar_pressure_if_needed!(
            rate,
            f.pe,
            w.ue,
            f.gradp,
            w.divue,
            f.J,
            T(2),
            zero(T),
            false,
            st.g,
        )
        @test all(signbit, rate)
        for c = 2:3
            fill!(w.ue[c], T(NaN))
        end
        foreach(a -> fill!(a, T(NaN)), f.J)
        @test H._complete_scalar_pressure_if_needed!(
            rate,
            f.pe,
            w.ue,
            f.gradp,
            w.divue,
            f.J,
            T(2),
            zero(T),
            false,
            st.g;
            derivative_risk = true,
        )
        @test all(iszero, rate)
        @test_throws ArgumentError H._complete_scalar_pressure_if_needed!(
            rate,
            f.pe,
            w.ue,
            f.gradp,
            w.divue,
            f.J,
            T(2),
            one(T),
            true,
            st.g;
            derivative_risk = true,
        )
        for kind in (:overflow, :all_underflow)
            rejected = fixture(T, kind)
            exact = exact_pressure(rejected)
            @test any(!iszero, exact)
            kind === :overflow ? (@test any(x -> abs(x) > Q(floatmax(T)), exact)) :
            (@test all(iszero, T.(exact)))
            @test_throws ArgumentError H._bpressure_rhs!(
                rejected.k1,
                rejected.pressure_work.kp[1],
                rejected.fields.B,
                rejected.fields.pe,
                rejected,
                1,
                1,
            )
        end
        # A nonfinite original operand must be rejected before exact extraction.
        w.ue[1][1] = T(Inf)
        rate[1] = T(NaN)
        @test_throws ArgumentError H._complete_scalar_pressure_if_needed!(
            rate,
            f.pe,
            w.ue,
            f.gradp,
            w.divue,
            f.J,
            T(2),
            zero(T),
            true,
            st.g,
        )
    end
end

@testset "Pressure derivative risk is checked in its consumed axis" begin
    for T in (Float32, Float64), axis = 1:2, direction in (-one(T), one(T))
        a = T === Float32 ? 100 : 600
        tiny = ldexp(one(T), -a)
        lengths = ntuple(d -> d == axis ? ldexp(T(2pi), a) : T(2pi), 2)
        g = FourierGrid((4, 4), lengths)
        @test g.kvec[axis][2] == tiny
        st = HybridStepper(g, HybridModel(EvolvingElectrons(T(2))), CIC(), 0)
        f, w = st.fields, st.pressure_work
        fill!(f.n, one(T))
        fill!(st.fn, one(T))
        fill!(f.pe, one(T))
        if axis == 1
            f.pe[:, 1] .= tiny .* T[2, 3, 2, 1]
        else
            f.pe[1, :] .= tiny .* T[2, 3, 2, 1]
        end
        fill!(f.ui[axis], direction / tiny)
        fill!(st.fui[axis], direction / tiny)
        original_pressure = copy(f.pe)
        ohms_law!(f, st.model, g)
        electron_velocity!(w.ue, st.fui, f.J, st.fn)
        divergence!(w.divue, w.ue, g)
        native = similar(f.pe)
        active = H._transport_pressure_rhs_with_activity!(
            native,
            f.pe,
            w.ue,
            f.gradp,
            w.divue,
            T(2),
            f.J,
            zero(T),
            false,
            Val(2),
        )
        @test all(iszero, f.gradp[axis])
        @test maximum(abs, f.gradp[3-axis]) >= T(0.25)
        @test !active && all(iszero, native)
        @test H._transport_pressure_derivative_range_risk(
            f.pe,
            w.ue,
            f.gradp,
            g,
            H._transport_derivative_gains(g),
        )
        derivative = centered_four(g, axis)
        original = Q.(f.pe)
        reference =
            -Q(direction) .*
            (axis == 1 ? derivative * original : original * transpose(derivative)) ./ Q(tiny)
        @test maximum(abs, reference) == Q(tiny) > Q(floatmin(T))
        H._bpressure_rhs!(st.k1, w.kp[1], f.B, f.pe, st, 1, 1)
        @test w.kp[1] == T.(reference)
        @test all(component -> all(iszero, component), st.k1)
        @test isequal(f.pe, original_pressure)
        @test all(==(direction / tiny), st.fui[axis])
    end
end

@testset "Complete Faraday consumes a tiny axis-active field beneath a null background" begin
    for T in (Float32, Float64), axis = 1:2
        a = T === Float32 ? 60 : 300
        tiny = ldexp(one(T), -a)
        lengths = ntuple(d -> d == axis ? ldexp(T(2pi), -a) : T(2pi), 2)
        g = FourierGrid((4, 4), lengths)
        @test g.kvec[axis][2] == inv(tiny)
        ideal, current, out = ntuple(_ -> ntuple(_ -> zeros(T, g.n), 3), 3)
        # This component's curl consumes only `axis`. The normal transverse
        # background is constant in that axis and contributes exactly zero.
        field = ideal[3-axis]
        fill!(field, one(T))
        axis == 1 ? (field[:, 1] .= tiny .* T[0, 1, 0, -1]) :
        (field[1, :] .= tiny .* T[0, 1, 0, -1])
        original = map(copy, ideal)
        derivative = centered_four(g, axis)
        exact = axis == 1 ? -derivative * Q.(field) : Q.(field) * transpose(derivative)
        @test maximum(abs, exact) == 1
        H._complete_transport_faraday_rhs!(out, ideal, current, zero(T), zero(T), g)
        @test isequal(out[3], T.(exact))
        @test all(iszero, out[1]) && all(iszero, out[2])
        @test isequal(ideal, original)
        @test all(component -> all(iszero, component), current)
    end
end
end
