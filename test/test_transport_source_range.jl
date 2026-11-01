using HybridPlasmaPIC, Test, Serialization

tsr_rk4(z) = 1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24

# Four-point real centered differentiation, independent of package FFTs. The
# even Nyquist mode has zero odd derivative and a separate positive projector.
function tsr_D4(g)
    k = Rational{BigInt}(g.kvec[1][2])
    return [mod(j - i, 4) == 1 ? k / 2 : mod(j - i, 4) == 3 ? -k / 2 : zero(k) for i = 1:4, j = 1:4]
end

@testset "Scaled public damping and independently paired stage energy" begin
    setprecision(BigFloat, 8192) do
        n = 8
        phase = [cis(2BigFloat(pi) * (i - 1) * (m - 1) / n) for i = 1:n, m = 1:n]
        for T in (Float32, Float64)
            scale = T === Float32 ? T(1e12) : T(1e100)
            amplitude = T === Float32 ? T(1e-12) : T(1e-80)
            etaH = T === Float32 ? T(1e38) : T(1e300)
            pressure = T === Float32 ? T(1e-38) : T(1e-270)
            dt = T === Float32 ? T(1e8) : T(1e98)
            gamma = T(5 / 3)
            g = FourierGrid((n,), (T(2pi) * scale,))
            ps = ParticleSet{1,T}(n)
            for i = 1:n
                ps.x[1][i] = T(i - 1) * g.dx[1]
            end
            set_density_weight!(ps, one(T), g)
            st = HybridStepper(
                g,
                HybridModel(TransportElectrons(gamma); η = zero(T), ηH = etaH),
                NGP(),
                n,
            )
            fill!(st.fields.pe, pressure)
            for i = 1:n
                st.fields.B[2][i] = amplitude * sin(T(2pi * (i - 1) / n))
            end
            init!(st, ps)
            foreach(copyto!, st.fui, st.fields.ui)
            copyto!(st.fn, st.fields.n)
            initial = copy(st.fields.B[2])
            current_hat =
                [sum(j -> BigFloat(st.fields.J[3][j]) * conj(phase[j, m]), 1:n) for m = 1:n]
            field_reference = [
                real(
                    sum(
                        m ->
                            BigFloat(etaH) *
                            BigFloat(g.kfull[1][m])^2 *
                            current_hat[m] *
                            phase[i, m],
                        1:n,
                    ),
                ) / n for i = 1:n
            ]
            @test all(iszero, st.fields.lapJ[3])
            @test all(iszero, st.pressure_work.Eforce[3])
            @test any(!iszero, st.fields.E[3])
            @test maximum(abs, BigFloat.(st.fields.E[3]) .- field_reference) /
                  maximum(abs, field_reference) < 2048eps(T)
            HybridPlasmaPIC._bpressure_rhs!(
                st.k1,
                st.pressure_work.kp[1],
                st.fields.B,
                st.fields.pe,
                st,
                1,
                1,
            )
            loss = -sum(sum(BigFloat.(st.fields.B[c]) .* BigFloat.(st.k1[c])) for c = 1:3)
            heating = sum(BigFloat.(st.pressure_work.kp[1])) / (BigFloat(gamma) - 1)
            @test loss > 0 && heating > 0
            @test abs(loss - heating) / max(loss, heating) < 4096eps(T)
            # Modal RK4 evolution is an independent time-composition oracle;
            # use the actual stored full/centered metadata and initial samples.
            hat = [sum(j -> BigFloat(initial[j]) * conj(phase[j, m]), 1:n) for m = 1:n]
            damping = [
                BigFloat(dt) *
                BigFloat(etaH) *
                BigFloat(g.kvec[1][m])^2 *
                BigFloat(g.kfull[1][m])^2 for m = 1:n
            ]
            @test all(a -> abs(tsr_rk4(-a)) <= 1, damping)
            reference =
                [real(sum(m -> tsr_rk4(-damping[m]) * hat[m] * phase[i, m], 1:n)) / n for i = 1:n]
            step!(st, ps, dt)
            @test maximum(abs, BigFloat.(st.fields.B[2]) .- reference) / BigFloat(amplitude) <
                  4096eps(T)
            @test all(all(isfinite, a) for a in (st.fields.B..., st.fields.E..., st.fields.pe))
            @test maximum(st.fields.pe) > pressure
        end
    end
end

function tsr_pressure_reference(g, p, u, J, gamma, etaH)
    Q = Rational{BigInt}
    derivative = tsr_D4(g)
    k = Q(g.kfull[1][3])
    projector = [isodd(i - j) ? -one(k) / 4 : one(k) / 4 for i = 1:4, j = 1:4]
    heat = zeros(Q, 4)
    for current in J
        jq = Q.(current)
        heat .+= (derivative * jq) .^ 2 .+ (k .* (projector * jq)) .^ 2
    end
    return -Q.(u[1]) .* (derivative * Q.(p)) .- Q(gamma) .* Q.(p) .* (derivative * Q.(u[1])) .+
           Q(gamma - one(gamma)) * Q(etaH) .* heat
end


@testset "Complete mixed-source public field and exact scalar context" begin
    for T in (Float32, Float64)
        maximum_value = floatmax(T)
        g = FourierGrid((4,), (T(4pi),))
        st = HybridStepper(
            g,
            HybridModel(TransportElectrons(T(5 / 3)); η = maximum_value, ηH = maximum_value),
            NGP(),
            0,
        )
        f = st.fields
        fill!(f.n, one(T))
        fill!(f.pe, one(T))
        fill!(f.B[1], one(T))
        f.B[2] .= T[0, 3, 0, -3]
        f.ui[2] .= -maximum_value .* T[1, 0, -1, 0]
        expected_current = T[1.5, 0, -1.5, 0]
        expected_laplacian = T[-0.375, 0, 0.375, 0]
        @test any(!isfinite, maximum_value .* expected_current)
        reference = setprecision(BigFloat, 8192) do
            T.(
                -BigFloat(maximum_value) .* BigFloat.(T[1, 0, -1, 0]) .+
                BigFloat(maximum_value) .* BigFloat.(expected_current) .-
                BigFloat(maximum_value) .* BigFloat.(expected_laplacian)
            )
        end
        @test ohms_law!(f, st.model, g) === f.E
        @test f.J[3] == expected_current
        @test f.lapJ[3] == expected_laplacian
        @test isequal(f.E[3], reference)
        @test all(all(isfinite, a) for a in (f.E..., f.J..., f.lapJ...))
        # This cold composition must use the canonical exact owner, independent
        # of the caller's MPFR precision/rounding and without mutating either.
        expected = reference[1]
        for bits in (11, 53, 4096), mode in (RoundNearest, RoundDown, RoundUp)
            setprecision(BigFloat, bits) do
                setrounding(BigFloat, mode) do
                    actual = HybridPlasmaPIC._complete_transport_field_value(
                        -maximum_value,
                        maximum_value,
                        T(1.5),
                        maximum_value,
                        Rational{BigInt}(T(-0.375)),
                    )
                    @test isequal(actual, expected)
                    @test precision(BigFloat) == bits && rounding(BigFloat) === mode
                end
            end
        end
        for value in (zero(T), nextfloat(zero(T)), floatmin(T), T(1.3), floatmax(T)),
            power in (-3000, 0, 3000)

            actual = HybridPlasmaPIC._exact_binary_scale(value, power)
            expected = setprecision(BigFloat, 8192) do
                Rational{BigInt}(BigFloat(value) * BigFloat(2)^power)
            end
            @test actual == expected
        end
    end
end

@testset "Public hidden induction with constant E and canceled complete pressure" begin
    setprecision(BigFloat, 8192) do
        for T in (Float32, Float64), variant in (:hidden_normal, :below_minimum)
            ek = T === Float32 ? 24 : 200
            ej = variant === :hidden_normal ? -10 : (T === Float32 ? -101 : -674)
            eh = T === Float32 ? -106 : -810
            ep =
                variant === :hidden_normal ? (T === Float32 ? -40 : -100) :
                (T === Float32 ? -140 : -1000)
            eu = eh + 2ej + ek - ep
            en = T === Float32 ? 60 : 500
            g = FourierGrid((4,), (ldexp(T(2pi), -ek),))
            etaH, dt, gamma = ldexp(T(27), eh), T(0.5), T(2)
            @test g.kvec[1][2] == ldexp(one(T), ek)
            ps = ParticleSet{1,T}(4)
            for i = 1:4
                ps.x[1][i] = T(i - 1) * g.dx[1]
                ps.v[1][i] = ldexp(T((-4, 34, -16, 46)[i]), eu)
                ps.v[2][i] = ldexp(one(T), 2ek + 10)
            end
            set_density_weight!(ps, ldexp(one(T), en), g)
            st = HybridStepper(
                g,
                HybridModel(TransportElectrons(gamma); η = zero(T), ηH = etaH),
                NGP(),
                4,
            )
            st.fields.pe .= ldexp.(T[1, 2, 4, 8], ep)
            fill!(st.fields.B[1], ldexp(one(T), -2ek - 10))
            st.fields.B[2] .= ldexp.(T[0, 1, 0, -1], ej - ek)
            init!(st, ps)
            foreach(copyto!, st.fui, st.fields.ui)
            copyto!(st.fn, st.fields.n)
            # A caller of the raw stage owns its midpoint moments. Supply exact
            # dyadic data rather than assume deposition is an exact round trip.
            st.fui[1] .= ldexp.(T[-4, 34, -16, 46], eu)
            HybridPlasmaPIC._bpressure_rhs!(
                st.k1,
                st.pressure_work.kp[1],
                st.fields.B,
                st.fields.pe,
                st,
                1,
                1,
            )
            pressure = tsr_pressure_reference(
                g,
                st.fields.pe,
                st.pressure_work.ue,
                st.fields.J,
                gamma,
                etaH,
            )
            @test all(iszero, pressure)
            @test all(iszero, st.pressure_work.kp[1])
            @test all(isone, st.Escr[3]) && all(isone, st.fields.E[3])
            @test all(iszero, tsr_D4(g) * Rational{BigInt}.(st.Escr[3]))
            @test any(!iszero, st.k1[2])
            @test maximum(abs, st.fields.lapJ[3]) > zero(T)
            initial = copy(st.fields.B[2])
            damping =
                BigFloat(dt) * BigFloat(etaH) * BigFloat(g.kvec[1][2])^2 * BigFloat(g.kfull[1][2])^2
            @test abs(tsr_rk4(-4damping)) < 1
            reference = BigFloat.(initial) .* tsr_rk4(-damping)
            if variant === :below_minimum
                # The first stage cancels exactly, but damping changes J in the
                # next stage while compression remains fixed. Its complete
                # pressure rate is then nonzero and below T's physical range.
                trial = initial .+ (dt / 2) .* st.k1[2]
                Q = Rational{BigInt}
                trial_current = (zeros(Q, 4), zeros(Q, 4), tsr_D4(g) * Q.(trial))
                next_pressure = tsr_pressure_reference(
                    g,
                    st.fields.pe,
                    st.pressure_work.ue,
                    trial_current,
                    gamma,
                    etaH,
                )
                @test any(!iszero, next_pressure)
                @test all(iszero, T.(next_pressure))
                tr_rejected(() -> step!(st, ps, dt), st, ps, "complete transport pressure rate")
            else
                step!(st, ps, dt)
                @test maximum(abs, BigFloat.(st.fields.B[2]) .- reference) /
                      ldexp(BigFloat(1), ej - ek) < 4096eps(T)
            end
            # The cold replay's derivative scratch must still describe the
            # carried pressure/current, and ions retain the independently
            # evaluated ideal force rather than a resistive-field subtraction.
            ideal = deepcopy(st.fields)
            ideal_model = HybridModel(st.model.closure; η = zero(T), ηH = zero(T))
            ohms_law!(ideal, ideal_model, g)
            @test isequal(st.pressure_work.Eforce, ideal.E)
            @test isequal(st.fields.J, ideal.J)
            @test isequal(st.fields.gradp, ideal.gradp)
            @test all(all(isfinite, a) for a in (st.fields.E..., st.fields.lapJ...))
        end
    end
end

function tsr_rate_fixture(::Type{T}, collection) where {T}
    scale = T === Float32 ? T(1e19) : T(1e150)
    amplitude = T === Float32 ? T(1e-22) : T(1e-160)
    etaH = T === Float32 ? T(1e38) : T(1e308)
    dt = T === Float32 ? T(1e36) : T(1e290)
    g = FourierGrid((8,), (T(2pi) * scale,))
    ps = ParticleSet{1,T}(8)
    for i = 1:8
        ps.x[1][i] = T(i - 1) * g.dx[1]
    end
    set_density_weight!(ps, one(T), g)
    populations = collection ? [ps, ParticleSet{1,T}(0; q = T(2), m = T(4))] : ps
    st = HybridStepper(
        g,
        HybridModel(TransportElectrons(T(5 / 3)); η = zero(T), ηH = etaH),
        NGP(),
        8,
    )
    h = sqrt(T(0.5))
    st.fields.B[2] .= amplitude .* T[0, h, 1, h, 0, -h, -1, -h]
    init!(st, populations)
    return st, populations, dt
end

@testset "Required physical rate range, rollback and restart" begin
    for T in (Float32, Float64), collection in (false, true)
        st, particles, dt = tsr_rate_fixture(T, collection)
        @test any(any(!iszero, a) for a in st.fields.E)
        @test all(all(isfinite, a) for a in (st.fields.E..., st.fields.J..., st.fields.lapJ...))
        mktempdir() do directory
            checkpoint = joinpath(directory, "physical-rate-domain.ser")
            save_checkpoint(checkpoint, st, particles)
            saved = deserialize(checkpoint)
            tr_rejected(
                () -> step!(st, particles, dt),
                st,
                particles,
                "complete transport pressure rate",
            )
            @test isequal(HybridPlasmaPIC._checkpoint_state(st, particles), saved)
            @test load_checkpoint!(st, particles, checkpoint) === st
            @test isequal(HybridPlasmaPIC._checkpoint_state(st, particles), saved)
            # The stored physical state is valid; this particular derivative
            # domain remains unsupported without dt-scaled RK storage.
            tr_rejected(
                () -> step!(st, particles, dt),
                st,
                particles,
                "complete transport pressure rate",
            )
            foreach(a -> fill!(a, zero(T)), st.fields.B)
            init!(st, particles)
            expected, other = deepcopy((st, particles))
            step!(st, particles, dt)
            step!(expected, other, dt)
            @test all(isequal(a, b) for (a, b) in zip(tr_state(st), tr_state(expected)))
            @test tr_clock(st) == tr_clock(expected)
            @test tr_same_particles(particles, tr_particles(other))
        end
    end
    # Exact cancellation is admissible. Classification of magnetic components
    # is joint and independent of the separate pressure-heating admission.
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2pi),))
        current = (zeros(T, 4), T[1, 0, -1, 0], zeros(T, 4))
        ideal = map(a -> -a, current) # etaH=1 and L(J)=-J for this mode.
        output = ntuple(_ -> fill(T(NaN), 4), 3)
        HybridPlasmaPIC._complete_transport_faraday_rhs!(output, ideal, current, zero(T), one(T), g)
        @test all(all(iszero, a) for a in output)
        long = FourierGrid((4,), (ldexp(T(2pi), T === Float32 ? 30 : 400),))
        tiny = ldexp(one(T), T === Float32 ? -100 : -1000)
        current = (zeros(T, 4), tiny .* T[1, 0, -1, 0], zeros(T, 4))
        ideal = ntuple(_ -> zeros(T, 4), 3)
        @test_throws ArgumentError HybridPlasmaPIC._complete_transport_faraday_rhs!(
            output,
            ideal,
            current,
            zero(T),
            one(T),
            long,
        )
    end
end

@testset "Inactive current still admits complete pressure and ideal source checks" begin
    for T in (Float32, Float64), kind in (:pressure, :ideal_induction), collection in (false, true)
        ek = kind === :pressure ? 0 : T === Float32 ? -60 : -400
        eu = kind === :pressure ? (T === Float32 ? -60 : -400) : (T === Float32 ? -40 : -400)
        ep = T === Float32 ? -100 : -800
        eb = T === Float32 ? -60 : -400
        g = FourierGrid((4,), (ldexp(T(2pi), -ek),))
        ps = ParticleSet{1,T}(4)
        for i = 1:4
            ps.x[1][i] = T(i - 1) * g.dx[1]
            ps.v[1][i] = ldexp(T((0, 1, 0, -1)[i]), eu)
        end
        set_density_weight!(ps, one(T), g)
        particles = collection ? [ps, ParticleSet{1,T}(0; q = T(2), m = T(4))] : ps
        st = HybridStepper(
            g,
            HybridModel(TransportElectrons(T(2)); η = zero(T), ηH = one(T)),
            NGP(),
            4,
        )
        kind === :pressure ? fill!(st.fields.pe, ldexp(one(T), ep)) :
        fill!(st.fields.B[2], ldexp(one(T), eb))
        init!(st, particles)
        @test all(all(iszero, a) for a in st.fields.J)
        derivative = tsr_D4(g)
        reference =
            kind === :pressure ?
            -Rational{BigInt}(2) .* Rational{BigInt}.(st.fields.pe) .*
            (derivative * Rational{BigInt}.(st.fields.ui[1])) :
            derivative * Rational{BigInt}.(st.fields.E[3])
        @test any(!iszero, reference) && all(iszero, T.(reference))
        pattern =
            kind === :pressure ? "complete transport pressure rate" :
            "complete transport magnetic rate"
        tr_rejected(() -> step!(st, particles, one(T)), st, particles, pattern)
        # A retry with the exact-zero source remains valid after rejection and
        # retains the same public arrays and transaction phase semantics.
        for population in tr_populations(particles), component in population.v
            fill!(component, zero(T))
        end
        init!(st, particles)
        expected, other = deepcopy((st, particles))
        step!(st, particles, one(T))
        step!(expected, other, one(T))
        @test all(isequal(a, b) for (a, b) in zip(tr_state(st), tr_state(expected)))
        @test tr_clock(st) == tr_clock(expected)
        @test tr_same_particles(particles, tr_particles(other))
    end
end

@testset "Source activity respects curl axes, constants and centered nulls" begin
    for T in (Float32, Float64), D = 1:3, n in (1, 2, 3, 4, 5, 6)
        g = FourierGrid(ntuple(_ -> n, D), ntuple(_ -> T(2pi), D))
        ideal = ntuple(_ -> zeros(T, g.n), 3)
        current = ntuple(_ -> zeros(T, g.n), 3)
        for I in CartesianIndices(g.n)
            ideal[1][I] = T(I[1])
        end
        @test !HybridPlasmaPIC._transport_faraday_sources_active(ideal, current, zero(T), one(T), g)
        fill!(ideal[1], zero(T))
        for I in CartesianIndices(g.n)
            ideal[3][I] = T(I[1])
        end
        @test HybridPlasmaPIC._transport_faraday_sources_active(
            ideal,
            current,
            zero(T),
            one(T),
            g,
        ) == (n > 2)
        if D >= 2
            fill!(ideal[3], zero(T))
            for I in CartesianIndices(g.n)
                ideal[1][I] = T(I[2])
            end
            @test HybridPlasmaPIC._transport_faraday_sources_active(
                ideal,
                current,
                zero(T),
                one(T),
                g,
            ) == (n > 2)
        end
        p = ones(T, g.n)
        velocity = ntuple(_ -> zeros(T, g.n), 3)
        gradient = ntuple(_ -> zeros(T, g.n), D)
        divergence = zeros(T, g.n)
        out = fill(T(NaN), g.n)
        @test !HybridPlasmaPIC._transport_pressure_rhs_with_activity!(
            out,
            p,
            velocity,
            gradient,
            divergence,
            T(2),
            current,
            zero(T),
            false,
            Val(D),
        )
        @test all(iszero, out)
    end
end
