module ElectronConductionTests
using HybridPlasmaPIC, Test, LinearAlgebra, Serialization
import JET
include("transport_hyperresistivity/reference.jl")
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function fixture(::Type{T}; κ = T(0.07), η = 0.0, ηH = 0.0, collection = false,
    density = one(T), constant_pressure = false, legacy = false) where {T}
    g = FourierGrid((8,), (T(2pi),))
    ps = ParticleSet{1,T}(16)
    load_lattice_1d!(ps, zero(T), g.L[1])
    set_density_weight!(ps, density, g)
    closure = legacy ? TransportElectrons(T(2)) : ConductingElectrons(T(2); κ)
    st = HybridStepper(g, HybridModel(closure; η, ηH), NGP(), 0)
    for i = 1:8
        x = T(i - 1) * g.dx[1]
        st.fields.pe[i] = constant_pressure ? T(0.2) : T(0.2) + T(0.02) * cos(x)
        st.fields.B[3][i] = constant_pressure ? zero(T) : T(0.03) * sin(x)
    end
    particles = collection ? [ps, ParticleSet{1,T}(0; q = T(2), m = T(4))] : ps
    init!(st, particles)
    return st, particles
end

function public_state(st, particles)
    f = st.fields
    values = map(fieldnames(typeof(f))) do name
        value = getfield(f, name)
        value isa Ref ? value[] : value
    end
    populations = particles isa ParticleSet ? (particles,) : particles
    return deepcopy((values, st.pressure_work.Eforce,
        [(p.x, p.v, p.weight, p.id, p.tag, p.q, p.m) for p in populations],
        st.time[], st.step[], st.last_dt[]))
end

function stage!(st)
    H._bpressure_rhs!(st.k1, st.pressure_work.kp[1], st.fields.B, st.fields.pe, st, 1, 1)
    return st.pressure_work.kp[1]
end

# Exact four-point Fourier matrices, independent of production FFTs. Full L
# retains Nyquist; the centered first derivative annihilates it.
function matrices(g)
    k = Q(g.kfull[1][2])
    D = [mod(j-i, 4) == 1 ? k/2 : mod(j-i, 4) == 3 ? -k/2 : zero(k) for i = 1:4, j = 1:4]
    L = [mod(j-i, 4) == 0 ? -3k^2/2 : mod(j-i, 4) == 2 ? -k^2/2 : k^2 for i = 1:4, j = 1:4]
    return D, L
end

@testset "Conducting closure and original quotient admission" begin
    for T in (Float32, Float64)
        c = ConductingElectrons(T(2); κ = T(0.1))
        @test closure_gamma(c) == T(2)
        @test fieldnames(TransportElectrons{T}) === (:γ,)
        for bad in (-one(T), T(Inf), T(NaN))
            @test_throws ArgumentError ConductingElectrons(T(2); κ = bad)
        end
        @test_throws ArgumentError ConductingElectrons(one(T); κ = one(T))
        @test_throws ArgumentError HybridModel(c; de2 = 0.1)
        g = FourierGrid((4,), (T(2pi),))
        @test_throws ArgumentError CAMCLStepper(g, HybridModel(c), NGP(), 0)
        @test_throws ArgumentError HallMHDModel(c)
        @test_throws ArgumentError electron_pressure!(ones(T, 4), ones(T, 4), c)
        quotient, p, n = zeros(Q, 4), T[2, 6, 10, 14], T[1, 3, 5, 7]
        @test H._conduction_quotient!(quotient, p, n, g) === quotient
        @test quotient == fill(Q(2), 4)
        n .= T[3, 5, 7, 9]
        H._conduction_quotient!(quotient, p, n, g)
        @test quotient == Q.(p) ./ Q.(n)
        @test_throws ArgumentError H._conduction_quotient!(p, p, n, g)
        @test_throws DimensionMismatch H._conduction_quotient!(zeros(Q, 3), p, n, g)
        for (pressure, density) in ((floatmax(T), T(0.5)), (nextfloat(zero(T)), T(2)))
            H._conduction_quotient!(quotient, fill(pressure, 4), fill(density, 4), g)
            @test quotient == fill(Q(pressure) / Q(density), 4)
        end
        for pe in (zeros(T, 4), fill(-zero(T), 4))
            H._conduction_quotient!(quotient, pe, n, g)
            @test all(iszero, quotient)
        end
        for (pressure, density) in ((-one(T), one(T)), (T(Inf), one(T)),
            (one(T), zero(T)), (one(T), T(NaN)))
            @test_throws ArgumentError H._conduction_quotient!(quotient, fill(pressure, 4), fill(density, 4), g)
        end
        u = Q(1) / (BigInt(1) << precision(T))
        for value in (Q(1) / 3, -Q(7) / 11, Q(floatmax(T)) / Q(nextfloat(zero(T))),
            Q(nextfloat(zero(T))) / Q(floatmax(T)), zero(Q))
            prepared = H._conduction_dyadic(T, value)
            @test ispow2(denominator(prepared))
            @test abs(value - prepared) <= u^2 * abs(value)
            @test iszero(prepared) == iszero(value)
        end
    end
    g = FourierGrid((4,), (Float32(2pi),))
    @test_throws ArgumentError HybridStepper(g, HybridModel(ConductingElectrons(2.0; κ = 1e-100)), NGP(), 0)
    @test_throws ArgumentError HybridStepper(g, HybridModel(ConductingElectrons(2.0; κ = 1e100)), NGP(), 0)
end

@testset "Independent variable-density conduction, full Nyquist and signed stages" begin
    for T in (Float32, Float64), dims in ((5,), (4,), (4, 5)), ηH in (0.0, 0.02)
        D = length(dims)
        g = FourierGrid(dims, ntuple(_ -> T(2pi), D))
        κ, η, γ = T(0.07), T(0.03), T(2)
        st = HybridStepper(g, HybridModel(ConductingElectrons(γ; κ); η, ηH), NGP(), 0)
        f = st.fields
        for I in CartesianIndices(dims)
            phase = sum(T(I[d] - 1) * g.dx[d] for d = 1:D)
            f.n[I] = T(1.1) + T(0.03) * sin(phase)
            st.fn[I], f.ninv[I] = f.n[I], inv(f.n[I])
            f.pe[I] = T(0.7) + T(0.04) * cos(phase)
            for c = 1:3
                f.B[c][I] = T(0.1 / c) * sin(phase + T(c))
                st.fui[c][I] = T(0.02 / c) * cos(phase + T(c / 3))
            end
        end
        legacy = th_reference(f.B, f.pe, f.n, st.fui, g.L, γ, η, T(ηH))
        # Remove exact axis-null backgrounds before converting the independent
        # direct-DFT oracle operands. No production preparation/FFT supplies it.
        quotient = Q.(f.pe) ./ Q.(f.n)
        residuals = [map(CartesianIndices(dims)) do I
            reference = CartesianIndex(ntuple(c -> c == d ? 1 : I[c], D))
            quotient[I] - quotient[reference]
        end for d = 1:D]
        lap = sum(th_apply(th_matrix(dims[d], g.L[d], :lap), Float64.(residuals[d]), d) for d = 1:D)
        conductive = Float64(γ - 1) * Float64(κ) * lap
        actual = stage!(st)
        @test st.pressure_work.quotient == quotient
        @test actual ≈ legacy.outp + conductive rtol = 2048eps(T) atol = 2048eps(T)
        computed = zeros(Q, dims)
        H._add_conduction_pressure!(computed, quotient, one(Q), g)
        scale = sum(maximum(abs2, Float64.(g.kfull[d])) * norm(Float64.(residuals[d])) for d = 1:D)
        @test norm(Float64.(computed) - lap) <= (2048eps(T) + 2.0^(-2precision(T))) * scale
        @test iszero(sum(computed))
        for c = 1:3
            @test st.Escr[c] ≈ legacy.E[c] rtol = 2048eps(T) atol = 2048eps(T)
            @test st.k1[c] ≈ legacy.outB[c] rtol = 2048eps(T) atol = 2048eps(T)
        end
        @test abs(sum(conductive)) <= 4096eps(T) * max(norm(conductive), 1)
    end
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2pi),))
        st = HybridStepper(g, HybridModel(ConductingElectrons(T(2); κ = T(0.1))), NGP(), 0)
        st.fn .= T[3, 6, 9, 12]
        st.fields.pe .= T[1, 2, 3, 4]
        st.fields.ninv .= inv.(st.fn)
        rate = stage!(st)
        @test st.pressure_work.quotient == fill(Q(1) / 3, 4)
        @test all(iszero, rate) # L(pe) is nonzero, but original p/n is constant.
        st.fn .= one(T)
        st.fields.ninv .= one(T)
        st.fields.pe .= T[2, 1, 2, 1]
        _, L = matrices(g)
        rate = stage!(st)
        expected = Q(T(0.1)) .* (L * Q.(st.fields.pe))
        @test rate ≈ T.(expected) rtol = 16eps(T)
        @test any(!iszero, rate)
        @test sum(rate) == zero(T)
        @test dot(st.pressure_work.quotient, Q.(rate)) < 0
    end
end

@testset "Original quotient variation, range rescue and exact null projection" begin
    for T in (Float32, Float64), ηH in (0.0, 0.02)
        g = FourierGrid((4,), (T(2pi),))
        M = ldexp(one(T), precision(T) - 1)
        st = HybridStepper(g, HybridModel(ConductingElectrons(T(2); κ = one(T)); ηH), NGP(), 0)
        st.fields.pe .= T[M, M + 1, M, M + 1]
        st.fn .= T[M + 1, M + 2, M + 1, M + 2]
        st.fields.ninv .= inv.(st.fn)
        quotient = Q.(st.fields.pe) ./ Q.(st.fn)
        @test all(==(T(quotient[1])), T.(quotient))
        delta = Q(1) / ((Q(M) + 1) * (Q(M) + 2))
        @test quotient[2] - quotient[1] == delta
        _, L = matrices(g)
        expected = L * quotient
        @test expected == 2 * Q(g.kfull[1][2])^2 * delta .* Q[1, -1, 1, -1]
        actual = stage!(st)
        @test all(x -> abs(x) >= floatmin(T), actual)
        @test actual ≈ T.(expected) rtol = 16eps(T)
        for (p, n, κ) in (
            (fill(floatmax(T), 4), T[0.5, 1, 0.5, 1], ldexp(one(T), -exponent(floatmax(T)))),
            (T[nextfloat(zero(T)), 0, nextfloat(zero(T)), 0], fill(T(2), 4), ldexp(one(T), exponent(floatmax(T)))),
        )
            st.model = HybridModel(ConductingElectrons(T(2); κ); ηH)
            st.fields.pe .= p
            st.fn .= n
            st.fields.ninv .= inv.(n)
            @test all(isfinite, st.fields.ninv)
            expected = Q(κ) .* (L * (Q.(p) ./ Q.(n)))
            @test all(x -> floatmin(T) < abs(x) < Q(floatmax(T)), expected)
            @test stage!(st) ≈ T.(expected) rtol = 16eps(T)
        end
        # The stiff first axis is exactly null; an odd active second axis
        # exercises the required rational projection after the last FFT.
        stiff = T === Float32 ? 30 : 300
        g2 = FourierGrid((4, 5), (ldexp(T(2pi), -stiff), T(2pi)))
        q = [Q(2) + Q(j - 1) / 7 for i = 1:4, j = 1:5]
        projected = zeros(Q, 4, 5)
        H._add_conduction_pressure!(projected, q, one(Q), g2)
        expected2 = th_apply(th_matrix(5, g2.L[2], :lap), Float64.(q .- q[:, 1:1]), 2)
        scale = maximum(abs2, Float64.(g2.kfull[2])) * norm(Float64.(q .- q[:, 1:1]))
        @test norm(Float64.(projected) - expected2) <= 2048eps(T) * scale
        @test all(iszero, sum(projected; dims = 2))
        @test all(projected[1, :] == projected[i, :] for i = 2:4)
    end
end

@testset "Complete original conduction before source range classification" begin
    for T in (Float32, Float64)
        ek, ep = T === Float32 ? (40, 60) : (300, 500)
        g = FourierGrid((4,), (ldexp(T(2pi), -ek),))
        @test g.kfull[1][2] == ldexp(one(T), ek)
        p = ldexp(T(3), ep)
        st = HybridStepper(g, HybridModel(ConductingElectrons(T(2); κ = p)), NGP(), 0)
        st.fields.pe .= p
        st.fn .= ldexp.(T[1, 1.5, 3, 1.5], ep)
        st.fields.ninv .= inv.(st.fn)
        st.fui[1] .= ldexp.(T[0, -1, 0, 1], ek - 1)
        D, L = matrices(g)
        temperature = Q[3, 2, 1, 2]
        compression = -2 * Q(p) .* (D * Q.(st.fui[1]))
        conduction = Q(p) .* (L * temperature)
        @test all(iszero, compression + conduction)
        @test any(x -> abs(x) > Q(floatmax(T)), conduction)
        @test all(iszero, stage!(st))
        # (gamma-1)*kappa underflows separately; full scaled-L source is normal.
        γ, κ = nextfloat(one(T)), nextfloat(zero(T))
        rescued = HybridStepper(g, HybridModel(ConductingElectrons(γ; κ)), NGP(), 0)
        rescued.fn .= one(T)
        rescued.fields.ninv .= one(T)
        rescued.fields.pe .= T[2, 1, 2, 1]
        @test iszero((γ - one(T)) * κ)
        expected = Q(γ - one(T)) * Q(κ) .* (L * Q.(rescued.fields.pe))
        actual = stage!(rescued)
        @test all(isfinite, actual) && any(!iszero, actual)
        @test actual ≈ T.(expected) rtol = 16eps(T)
        # A required physical Laplacian would overflow, while the complete
        # source fits; it remains normalized until its coefficient is consumed.
        exponent_p = T === Float32 ? 80 : 500
        exponent_k = T === Float32 ? 30 : 300
        small = ldexp(one(T), -exponent_p)
        wide = FourierGrid((4,), (ldexp(T(2pi), -exponent_k),))
        scaled = HybridStepper(wide, HybridModel(ConductingElectrons(T(2); κ = small)), NGP(), 0)
        scaled.fn .= one(T)
        scaled.fields.ninv .= one(T)
        scaled.fields.pe .= ldexp.(T[4, 2, 4, 2], exponent_p)
        _, Lwide = matrices(wide)
        physical_lap = Lwide * Q.(scaled.fields.pe)
        @test any(x -> abs(x) > Q(floatmax(T)), physical_lap)
        expected = Q(small) .* physical_lap
        @test stage!(scaled) ≈ T.(expected) rtol = 16eps(T)
        scaled.model = HybridModel(ConductingElectrons(T(2); κ = one(T)))
        @test_throws ArgumentError stage!(scaled)
        scaled.model = HybridModel(ConductingElectrons(T(2); κ = small))
        scaled.fields.pe .*= T(0.5)
        @test stage!(scaled) ≈ T.(expected ./ 2) rtol = 16eps(T)
    end
end

@testset "Frozen conduction RK4, positivity rejection and timestep cap" begin
    function decay_error(h)
        T = Float64
        g = FourierGrid((8,), (T(2pi),))
        st = HybridStepper(g, HybridModel(ConductingElectrons(T(2); κ = 0.2)), NGP(), 0)
        st.fn .= one(T)
        st.fields.ninv .= one(T)
        initial = [0.2cos((i-1) * g.dx[1]) for i = 1:8]
        st.fields.pe .= 1 .+ initial
        for cycle = 1:round(Int, 0.4 / h)
            H._rk4_Bp!(st, h, cycle)
        end
        expected = 1 .+ exp(-0.2 * 0.4) .* initial
        return norm(st.fields.pe - expected)
    end
    errors = decay_error.((0.2, 0.1, 0.05))
    @test errors[1] > errors[2] > errors[3] > 0
    @test 12 < errors[1] / errors[2] < 20
    @test 12 < errors[2] / errors[3] < 20
    for T in (Float32, Float64)
        st, ps = fixture(T; κ = T(1), constant_pressure = true)
        st.fn .= one(T)
        kmax2 = sum(maximum(abs2, axis) for axis in st.g.kfull)
        boundary = T(H._integrator_real_axis_limit(:rk4) / kmax2)
        @test H._validate_conduction_substep(st, prevfloat(boundary)) === nothing
        before = public_state(st, ps)
        @test_throws ArgumentError step!(st, ps, nextfloat(boundary))
        @test isequal(public_state(st, ps), before)
        suggested = recommended_dt(st; NB = 3, safety = 0.8)
        cap = Q(3) * Q(0.8) * Q(H._integrator_real_axis_limit(:rk4)) * Q(T(st.model.nfloor)) /
            (Q(T(st.model.closure.κ)) * sum(maximum(abs, Q.(axis))^2 for axis in st.g.kfull))
        @test suggested isa T && Q(suggested) <= cap < Q(nextfloat(suggested))
        # safety=1 must not recommend exact equality at the strict guard;
        # NB>1 also exercises the actual grid-precision substep division.
        strict = HybridStepper(FourierGrid((4,), (T(2pi),)),
            HybridModel(ConductingElectrons(T(2); κ = T(1024)); nfloor = 1.0), NGP(), 0)
        strict.fn .= one(T)
        for NB in (1, 3)
            suggestion = recommended_dt(strict; NB, safety = 1.0)
            @test suggestion > zero(T)
            @test H._validate_conduction_substep(strict, suggestion / NB) === nothing
            strictcap = Q(NB) * Q(H._integrator_real_axis_limit(:rk4)) /
                (Q(T(1024)) * sum(maximum(abs, Q.(axis))^2 for axis in strict.g.kfull))
            @test Q(suggestion) <= strictcap
            if Q(T(strictcap)) == strictcap
                @test suggestion < T(strictcap)
            end
        end
        # A stored BigFloat safety factor can exceed the caller's MPFR
        # precision. Its exact value, not a new low-precision rationalization,
        # supplies the conservative cap.
        wide_safety = setprecision(BigFloat, 256) do
            BigFloat(3) / 4 + ldexp(BigFloat(1), -100)
        end
        low = setprecision(BigFloat, 2) do
            recommended_dt(strict; safety = wide_safety)
        end
        original_safety = Q(3) / 4 + Q(1) / (BigInt(1) << 100)
        exact_cap = original_safety * Q(H._integrator_real_axis_limit(:rk4)) /
            (Q(T(1024)) * sum(maximum(abs, Q.(axis))^2 for axis in strict.g.kfull))
        @test low > zero(T) && Q(low) <= exact_cap
        # No positive full step is admitted below the existing dt resolution.
        strict.model = HybridModel(ConductingElectrons(T(2); κ = floatmax(T)); nfloor = 1.0)
        @test recommended_dt(strict) == zero(T)
        # Full spectral L has a negative off-diagonal; finite stage pressure is
        # checked, not clipped. Public failure restores fields/force/particles.
        g = FourierGrid((4,), (T(2pi),))
        ions = ParticleSet{1,T}(8)
        load_lattice_1d!(ions, zero(T), g.L[1])
        set_density_weight!(ions, one(T), g)
        negative = HybridStepper(g, HybridModel(ConductingElectrons(T(2); κ = T(0.1))), NGP(), 0)
        negative.fields.pe .= T[0, 0, 1, 0]
        init!(negative, ions)
        before = public_state(negative, ions)
        @test_throws ArgumentError step!(negative, ions, T(0.001))
        @test isequal(public_state(negative, ions), before)
        for (pressure, density) in ((floatmax(T), T(0.5)), (nextfloat(zero(T)), T(2)))
            constant, pops = fixture(T; κ = T(0.001), density, constant_pressure = true)
            constant.fields.pe .= pressure
            init!(constant, pops)
            step!(constant, pops, T(0.0001))
            @test constant.fields.pe == fill(pressure, 8)
            @test constant.time[] == T(0.0001) && constant.step[] == 1
        end
    end
end

@testset "Zero conductivity legacy reduction and complete restart policy" begin
    for T in (Float32, Float64), collection in (false, true), (η, ηH) in ((0.0, 0.0), (0.02, 0.01))
        old, ions = fixture(T; legacy = true, η, ηH, collection)
        inactive, tracers = fixture(T; κ = zero(T), η, ηH, collection)
        @test old.pressure_work.quotient === nothing
        @test inactive.pressure_work.quotient isa Array{Q}
        @test isequal(public_state(old, ions), public_state(inactive, tracers))
        for dt in (T(0.001), T(0.002))
            step!(old, ions, dt; NB = 2)
            step!(inactive, tracers, dt; NB = 2)
            @test isequal(public_state(old, ions), public_state(inactive, tracers))
        end
        st, particles = fixture(T; collection, η, ηH)
        step!(st, particles, T(0.001); NB = 2)
        mktempdir() do dir
            path = joinpath(dir, "state.bin")
            archive = joinpath(dir, "archive.bin")
            save_checkpoint(path, st, particles)
            archive_run(archive, st, particles; rng_seed = 17)
            saved = deserialize(path)
            @test saved.transport_policy_version === 1 && saved.conduction_policy_version === 2
            @test isequal(saved, load_archive(archive).state)
            restored, pop = fixture(T; κ = zero(T), collection, η, ηH)
            theta = restored.pressure_work.quotient
            load_checkpoint!(restored, pop, path)
            @test restored.pressure_work.quotient === theta
            @test restored.model.closure.κ == st.model.closure.κ
            @test isequal(H._checkpoint_state(restored, pop), saved)
            step!(restored, pop, T(0.001); NB = 2)
            step!(st, particles, T(0.001); NB = 2)
            @test isequal(public_state(restored, pop), public_state(st, particles))
            drop_policy = (; (k => v for (k, v) in pairs(saved) if k != :conduction_policy_version)...)
            bad = reinterpret(ConductingElectrons{T}, T[T(2), -one(T)])[1]
            malformed = (drop_policy, merge(saved, (; conduction_policy_version = 1)),
                merge(saved, (; conduction_policy_version = 2.0)),
                merge(saved, (; conduction_policy_version = true)),
                merge(saved, (; conduction_policy_version = BigInt(2))),
                merge(saved, (; model = HybridModel(bad, η, ηH, 1e-6, 0.0))),
                merge(saved, (; Eforce = (saved.Eforce[1], saved.Eforce[1], saved.Eforce[3]))))
            for payload in malformed
                before = public_state(restored, pop)
                @test_throws ArgumentError H._restore_checkpoint_state!(restored, pop, payload)
                @test isequal(public_state(restored, pop), before)
            end
        end
    end
end

@testset "Conduction scratch and numerical context preflight" begin
    for T in (Float32, Float64)
        st, ps = fixture(T)
        theta = st.pressure_work.quotient
        # The concrete Q workspace excludes safe assignment of a physical T
        # array. Exercise its structural guard without manufacturing an alias.
        @test_throws ArgumentError H._validate_conduction_quotient_storage(st, ps, st.fields.pe)
        st.pressure_work.quotient = zeros(Q, 3)
        before = public_state(st, ps)
        @test_throws DimensionMismatch step!(st, ps, T(0.001))
        @test isequal(public_state(st, ps), before)
        st.pressure_work.quotient = theta
        for (mode, flush) in ((RoundDown, false), (RoundNearest, true))
            saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
            before = public_state(st, ps)
            try
                Base.Rounding.setrounding_raw(Float64, Base.Rounding.to_fenv(mode))
                Base.set_zero_subnormals(flush)
                @test_throws ArgumentError step!(st, ps, T(0.001))
            finally
                Base.Rounding.setrounding_raw(Float64, saved[1])
                Base.set_zero_subnormals(saved[2])
            end
            @test isequal(public_state(st, ps), before)
        end
        for P in (typeof(ps), Vector{typeof(ps)})
            for (operation, signature) in ((init!, Tuple{typeof(st),P}),
                (step!, Tuple{typeof(st),P,T}), (H._checkpoint_state, Tuple{typeof(st),P}))
                report = JET.report_call(operation, signature; sourceinfo = :full, stacktrace_types_limit = 0)
                if !isempty(JET.get_reports(report))
                    show(IOContext(stdout, :limit => false), MIME"text/plain"(), report)
                end
                @test isempty(JET.get_reports(report))
            end
        end
    end
end
end
