using HybridPlasmaPIC, Test

function thr_expected(initial::T, ηH, lap) where {T}
    setprecision(BigFloat, 8192) do
        T(BigFloat(initial) - BigFloat(ηH) * BigFloat(lap))
    end
end

function thr_fixture(::Type{T}, ::Val{D}, kind) where {T,D}
    short = kind === :lapJ_overflow
    length1 = short ? T(2pi) * (T === Float32 ? T(1e-12) : T(1e-100)) : T(2pi)
    n = ntuple(d -> d == 1 ? 8 : 2, D)
    g = FourierGrid(n, ntuple(d -> d == 1 ? length1 : T(2pi), D))
    ηH =
        short ? (T === Float32 ? T(1e-26) : T(1e-300)) :
        (kind === :sum_overflow ? floatmax(T) / T(2) : floatmax(T))
    γ = short ? T(5 / 3) : T(1.1)
    st = HybridStepper(g, HybridModel(TransportElectrons(γ); η = zero(T), ηH), CIC(), 0)
    f = st.fields
    for a in (f.n, f.ninv, st.fn, f.pe)
        fill!(a, one(T))
    end
    for I in CartesianIndices(n)
        θ = T(2pi * (I[1] - 1) / n[1])
        if short
            f.B[2][I] = (T === Float32 ? T(1e4) : T(1e80)) * sin(θ)
        else
            f.B[2][I] = T(3) + T(1.2) * sin(θ)
            scale =
                kind === :sum_overflow ? -T(0.75) * floatmax(T) :
                (kind === :product_cancel ? T(0.6) * floatmax(T) : zero(T))
            f.ui[1][I] = scale * (cos(θ) / f.B[2][I])
            st.fui[1][I] = f.ui[1][I]
        end
    end
    return st
end

# Direct wide DFT on one nonconstant axis, independent of package FFT execution.
function thr_lap_reference(values, full_modes, phase)
    n = length(values)
    hat = [sum(j -> BigFloat(values[j]) * conj(phase[j, m]), 1:n) for m = 1:n]
    return [real(sum(m -> -BigFloat(full_modes[m])^2 * hat[m] * phase[i, m], 1:n)) / n for i = 1:n]
end
thr_line(a) = view(a, :, ntuple(_ -> 1, ndims(a) - 1)...)

@testset "Hyperresistive scalar range and ordinary arithmetic" begin
    for T in (Float32, Float64)
        mx, tiny = floatmax(T), nextfloat(zero(T))
        initials = T[-mx, -1, -0.0, 0, tiny, floatmin(T), 1, nextfloat(one(T)), mx]
        coefficients = T[tiny, floatmin(T), eps(T), 1, nextfloat(one(T)), 2, mx]
        laplacians = T[-mx, -1, -tiny, 0, tiny, prevfloat(one(T)), 1, mx]
        for initial in initials, ηH in coefficients, lap in laplacians
            ordinary = initial - ηH * lap
            if isfinite(ordinary)
                actual = HybridPlasmaPIC._hyperresistive_field_value(initial, ηH, lap)
                @test actual isa T && isequal(actual, ordinary)
            else
                expected = thr_expected(initial, ηH, lap)
                if isfinite(expected)
                    actual = HybridPlasmaPIC._hyperresistive_field_value(initial, ηH, lap)
                    @test actual isa T && isequal(actual, expected)
                else
                    @test_throws ArgumentError HybridPlasmaPIC._hyperresistive_field_value(
                        initial,
                        ηH,
                        lap,
                    )
                end
            end
        end
        # An unconditional fused/exact expression would alter a finite old result.
        initial, ηH, lap = one(T), nextfloat(one(T)), prevfloat(one(T))
        @test !isequal(initial - ηH * lap, fma(-ηH, lap, initial))
        @test isequal(
            HybridPlasmaPIC._hyperresistive_field_value(initial, ηH, lap),
            initial - ηH * lap,
        )
        for bad in (T(Inf), -T(Inf), T(NaN))
            @test_throws ArgumentError HybridPlasmaPIC._hyperresistive_field_value(
                bad,
                one(T),
                one(T),
            )
            @test_throws ArgumentError HybridPlasmaPIC._hyperresistive_field_value(
                one(T),
                bad,
                one(T),
            )
            @test_throws ArgumentError HybridPlasmaPIC._hyperresistive_field_value(
                one(T),
                one(T),
                bad,
            )
        end
        # The existing exact owner must neither depend on nor change MPFR context.
        expected = thr_expected(mx, T(2), T(0.75) * mx)
        for bits in (11, 53, 4096), mode in (RoundNearest, RoundDown, RoundUp)
            setprecision(BigFloat, bits) do
                setrounding(BigFloat, mode) do
                    actual = HybridPlasmaPIC._hyperresistive_field_value(mx, T(2), T(0.75) * mx)
                    @test isequal(actual, expected)
                    @test precision(BigFloat) == bits && rounding(BigFloat) === mode
                end
            end
        end
    end
    # Exact threshold neighbor where minimum Julia's software fma returned -Inf.
    initial = reinterpret(Float64, 0x7c8ffffffffffff9)
    @test HybridPlasmaPIC._hyperresistive_field_value(
        initial,
        nextfloat(1.0),
        prevfloat(floatmax(Float64)),
    ) === -floatmax(Float64)
end

@testset "Public hyperresistive workspace and field range" begin
    setprecision(BigFloat, 8192) do
        n = 8
        phase = [cis(2BigFloat(pi) * (i - 1) * (m - 1) / n) for i = 1:n, m = 1:n]
        for T in (Float32, Float64),
            D = 1:3,
            kind in (:lapJ_overflow, :product_cancel, :product_overflow, :sum_overflow)

            st = thr_fixture(T, Val(D), kind)
            f, g = st.fields, st.g
            ideal = HybridModel(st.model.closure; η = zero(T), ηH = zero(T))
            ohms_law!(f, ideal, g)
            force = map(copy, f.E)
            reference = thr_lap_reference(thr_line(f.J[3]), g.kfull[1], phase)
            before = map(copy, tr_fields(f))
            if kind === :lapJ_overflow
                @test maximum(abs, reference) > BigFloat(floatmax(T))
                @test maximum(abs, BigFloat(st.model.ηH) .* reference) < BigFloat(floatmax(T))
                error = try
                    ohms_law!(f, st.model, g)
                    nothing
                catch caught
                    caught
                end
                @test error isa ArgumentError && occursin("Laplacian", sprint(showerror, error))
                @test any(!isequal(a, saved) for (a, saved) in zip(tr_fields(f), before))
                @test_throws ArgumentError HybridPlasmaPIC._bpressure_rhs!(
                    st.k1,
                    st.pressure_work.kp[1],
                    f.B,
                    f.pe,
                    st,
                    1,
                    1,
                )
            elseif kind === :product_cancel
                @test ohms_law!(f, st.model, g) === f.E
                @test all(all(isfinite, a) for a in (f.E..., f.lapJ...))
                @test any(!isfinite(T(st.model.ηH) * x) for x in f.lapJ[3])
                for c = 1:3, I in eachindex(f.E[c])
                    ordinary = force[c][I] - T(st.model.ηH) * f.lapJ[c][I]
                    expected =
                        isfinite(ordinary) ? ordinary :
                        thr_expected(force[c][I], T(st.model.ηH), f.lapJ[c][I])
                    @test isequal(f.E[c][I], expected)
                end
                actual = thr_line(f.lapJ[3])
                @test maximum(abs, BigFloat.(actual) .- reference) <=
                      256eps(T) * maximum(abs, reference)
                HybridPlasmaPIC._bpressure_rhs!(st.k1, st.pressure_work.kp[1], f.B, f.pe, st, 1, 1)
                @test all(all(isfinite, a) for a in (st.k1..., st.pressure_work.kp[1]))
            else
                error = try
                    ohms_law!(f, st.model, g)
                    nothing
                catch caught
                    caught
                end
                @test error isa ArgumentError &&
                      occursin("electric field", sprint(showerror, error))
                @test all(all(isfinite, a) for a in f.lapJ)
                @test_throws ArgumentError HybridPlasmaPIC._bpressure_rhs!(
                    st.k1,
                    st.pressure_work.kp[1],
                    f.B,
                    f.pe,
                    st,
                    1,
                    1,
                )
            end
        end
    end
    # Inactive hyperresistivity retains its complete prior no-op, even for scratch.
    for T in (Float32, Float64), D = 1:3
        st = thr_fixture(T, Val(D), :product_cancel)
        for a in st.fields.J
            fill!(a, T(NaN))
        end
        for a in st.fields.lapJ
            fill!(a, T(Inf))
        end
        before = map(copy, (st.fields.E..., st.fields.J..., st.fields.lapJ...))
        @test HybridPlasmaPIC._add_hyperresistive_field!(
            st.fields.E,
            st.fields.J,
            st.fields.lapJ,
            zero(T),
            st.g,
        ) === st.fields.E
        @test all(
            isequal(a, b) for
            (a, b) in zip((st.fields.E..., st.fields.J..., st.fields.lapJ...), before)
        )
    end
end

@testset "Hyperresistive range rejection restores coupled public state" begin
    for T in (Float32, Float64), collection in (false, true)
        st = thr_fixture(T, Val(1), :lapJ_overflow)
        ps = ParticleSet{1,T}(16)
        load_lattice!(ps, (zero(T),), st.g.L, (16,))
        set_density_weight!(ps, one(T), st.g)
        populations = collection ? [ps, ParticleSet{1,T}(0; q = T(2), m = T(4))] : ps
        identities = tr_state(st)
        tr_rejected(() -> init!(st, populations), st, populations, "Laplacian")
        @test all(a === b for (a, b) in zip(identities, tr_state(st)))
        badB = map(copy, st.fields.B)
        foreach(a -> fill!(a, zero(T)), st.fields.B)
        init!(st, populations)
        mktempdir() do directory
            checkpoint = joinpath(directory, "valid-before-range.ser")
            save_checkpoint(checkpoint, st, populations)
            foreach(copyto!, st.fields.B, badB)
            dt = T === Float32 ? T(1e-33) : T(1e-285)
            identities = tr_state(st)
            tr_rejected(() -> step!(st, populations, dt), st, populations, "Laplacian")
            @test all(a === b for (a, b) in zip(identities, tr_state(st)))
            # Existing persistence restores physical state; derived scratch was
            # restored by the failed step's existing transaction, not the payload.
            load_checkpoint!(st, populations, checkpoint)
            control = thr_fixture(T, Val(1), :lapJ_overflow)
            foreach(a -> fill!(a, zero(T)), control.fields.B)
            other = ParticleSet{1,T}(16)
            load_lattice!(other, (zero(T),), control.g.L, (16,))
            set_density_weight!(other, one(T), control.g)
            other_populations =
                collection ? [other, ParticleSet{1,T}(0; q = T(2), m = T(4))] : other
            init!(control, other_populations)
            step!(st, populations, dt)
            step!(control, other_populations, dt)
            @test all(isequal(a, b) for (a, b) in zip(tr_state(st), tr_state(control)))
            @test tr_clock(st) == tr_clock(control)
            @test tr_same_particles(populations, tr_particles(other_populations))
        end
    end
end
