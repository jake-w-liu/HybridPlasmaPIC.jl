module FieldAlignedConductionTests
using HybridPlasmaPIC, Test, LinearAlgebra, Serialization, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

# Independent four-node real Fourier rows. D kills Nyquist; N is the positive
# real k_N projector. No production derivative, preparation or FFT enters it.
function row4(a, I::CartesianIndex{D}, d, g, nyquist = false) where {D}
    if nyquist
        return abs(Q(g.kfull[d][3])) * sum(
            j ->
                (iseven(I[d] - j) ? 1 : -1) *
                Q(a[CartesianIndex(ntuple(c -> c == d ? j : I[c], D))]),
            1:4,
        ) / 4
    end
    next = CartesianIndex(ntuple(c -> c == d ? mod1(I[c] + 1, 4) : I[c], D))
    previous = CartesianIndex(ntuple(c -> c == d ? mod1(I[c] - 1, 4) : I[c], D))
    return Q(g.kfull[d][2]) * (Q(a[next]) - Q(a[previous])) / 2
end

function reference(theta, B, parallel, perpendicular, g::FourierGrid{D}) where {D}
    grad = ntuple(d -> map(I -> row4(theta, I, d, g), CartesianIndices(g.n)), D)
    difference = Q(parallel) - Q(perpendicular)
    norm2 = map(I -> sum(c -> Q(B[c][I])^2, 1:3), CartesianIndices(g.n))
    A(i, j, I) =
        (i == j ? Q(perpendicular) : zero(Q)) + difference * Q(B[i][I]) * Q(B[j][I]) / norm2[I]
    out = zeros(Q, g.n)
    dissipation = zero(Q)
    for d = 1:D
        flux = map(I -> sum(j -> A(d, j, I) * grad[j][I], 1:D), CartesianIndices(g.n))
        z = map(I -> row4(theta, I, d, g, true), CartesianIndices(g.n))
        heat = map(I -> A(d, d, I) * z[I], CartesianIndices(g.n))
        for I in CartesianIndices(g.n)
            out[I] += row4(flux, I, d, g) - row4(heat, I, d, g, true)
            dissipation += grad[d][I] * flux[I] + z[I] * heat[I]
        end
    end
    return out, dissipation
end

function fixture(
    ::Type{T},
    ::Val{D};
    closure = FieldAlignedConductingElectrons(
        T(2);
        κ_parallel = T(1 / 16),
        κ_perpendicular = T(1 / 64),
    ),
    ηH = 0.0,
    shape = CIC(),
) where {T,D}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
    ps = ParticleSet{D,T}(2 * prod(g.n))
    sine, cosine = T[0, 1, 0, -1], T[1, 0, -1, 0]
    for (p, I) in enumerate(CartesianIndices(g.n)), copy = 0:1
        j = p + copy * prod(g.n)
        for d = 1:D
            ps.x[d][j] = T(I[d] - 1) + T(1 / 4)
        end
        ps.weight[j] = T(1 / 2) + sine[I[1]] / 32
        for c = 1:3
            ps.v[c][j] = (copy == 0 ? T(1 / 32) : -T(1 / 64)) * sine[I[mod1(c, D)]]
        end
    end
    st = HybridStepper(g, HybridModel(closure; η = 0.01, ηH), shape, nparticles(ps))
    for I in CartesianIndices(g.n)
        st.fields.pe[I] = one(T) + sum(cosine[I[d]] / T(32d) for d = 1:D)
        st.fields.B[1][I] = T(1 / 4) + (D > 1 ? sine[I[2]] / 64 : zero(T))
        st.fields.B[2][I] = T(1 / 2) + sine[I[1]] / 32
        st.fields.B[3][I] = T(1) + (D == 3 ? cosine[I[3]] / 64 : sine[I[1]] / 64)
    end
    return st, ps
end

function snapshot(st, ps)
    return deepcopy((
        fields = H._hybrid_field_state(st.fields),
        floor = st.fields.floor_count[],
        force = st.pressure_work.Eforce,
        particles = H._checkpoint_particle_state(ps),
        clock = (st.time[], st.step[], st.last_dt[]),
    ))
end

@testset "Exact tensor target and two-preparation actual channels" begin
    for T in (Float32, Float64), D in (1, 2)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        closure =
            FieldAlignedConductingElectrons(T(2); κ_parallel = T(3 / 8), κ_perpendicular = T(1 / 8))
        p = map(I -> T(3 + I[1] + (D == 2 ? I[2] : 0)), CartesianIndices(g.n))
        n = map(I -> T(2I[1] + 1), CartesianIndices(g.n))
        theta = Q.(p) ./ Q.(n)
        B = ntuple(c -> map(I -> T(c) + T(I[1] - 1) / T(8c), CartesianIndices(g.n)), 3)
        target, dissipation = reference(theta, B, closure.κ_parallel, closure.κ_perpendicular, g)
        @test sum(target) == 0
        @test sum(theta .* target) == -dissipation
        @test dissipation >= 0
        columns = map(1:prod(g.n)) do j
            basis = zeros(Q, g.n)
            basis[j] = 1
            vec(first(reference(basis, B, closure.κ_parallel, closure.κ_perpendicular, g)))
        end
        @test hcat(columns...) == transpose(hcat(columns...))
        actual = zeros(Q, g.n)
        work = H._field_aligned_storage(g, closure)
        measured =
            H._field_aligned_conduction!(actual, theta, one(Q), closure, B, work, g; budget = true)
        # Qualified small-grid FFT/preparation acceptance, not a global FFTW
        # bound or an assertion that the computed map is exactly self-adjoint.
        scale = max(norm(Float64.(target)), 1.0)
        @test norm(Float64.(actual - target)) <= 8192eps(T) * scale
        @test sum(actual) == 0
        @test measured >= 0
        @test abs(Float64(sum(theta .* actual) + measured)) <=
              8192eps(T) * scale * max(norm(Float64.(theta)), 1)
        reversed = zeros(Q, g.n)
        H._field_aligned_conduction!(reversed, theta, one(Q), closure, map(a -> -a, B), work, g)
        @test actual == reversed
        if D == 1
            # Out-of-plane B contributes to the projector denominator.
            transverse = (zeros(T, g.n), ones(T, g.n), fill(T(2), g.n))
            rate = zeros(Q, g.n)
            H._field_aligned_conduction!(rate, theta, one(Q), closure, transverse, work, g)
            isotropic = zeros(Q, g.n)
            H._add_conduction_pressure!(isotropic, theta, Q(closure.κ_perpendicular), g)
            @test norm(Float64.(rate - isotropic)) <= 8192eps(T) * scale
            for scale_exponent in (exponent(floatmax(T)) - 3, exponent(floatmin(T)) + 3)
                huge_or_tiny = ntuple(c -> fill(ldexp(T(c), scale_exponent), g.n), 3)
                rate .= 0
                H._field_aligned_conduction!(rate, theta, one(Q), closure, huge_or_tiny, work, g)
                expected = first(
                    reference(theta, huge_or_tiny, closure.κ_parallel, closure.κ_perpendicular, g),
                )
                @test norm(Float64.(rate - expected)) <= 8192eps(T) * scale
            end
        end
    end
end

@testset "Charged joint lifecycle and instantaneous fresh-density budget" begin
    for (T, D, shape) in ((Float32, 1, NGP()), (Float64, 2, CIC()), (Float64, 3, TSC()))
        st, ps = fixture(T, Val(D); shape, ηH = 0.01)
        init!(st, ps)
        before = snapshot(st, ps)
        scratch =
            (st.pressure_work.quotient, H._field_aligned_arrays(st.pressure_work.conduction)...)
        step!(st, ps, T(1 / 1024); NB = 2)
        @test !isequal(before.fields.pe, st.fields.pe)
        @test any(c -> !isequal(before.fields.B[c], st.fields.B[c]), 1:3)
        @test !isequal(before.particles.v, ps.v)
        state = snapshot(st, ps)
        d = electron_conduction_diagnostics(st, ps)
        @test isequal(snapshot(st, ps), state)
        @test d.scope === :instantaneous && d.time == st.time[] && d.step == st.step[]
        @test d.internal_energy_rate == 0 && d.channel_dissipation >= 0
        @test d.adjoint_defect ==
              d.variance_rate +
              Q(H._pressure_heat_gamma(T(st.model.closure.γ))) * d.channel_dissipation
        # Plain restore omits derived n; the diagnostic must redeposit it.
        saved_n = copy(st.fields.n)
        st.fields.n .= T(99)
        @test electron_conduction_diagnostics(st, ps) == d
        copyto!(st.fields.n, saved_n)
        mktempdir() do dir
            raw, archive = joinpath(dir, "state.bin"), joinpath(dir, "state.hpr")
            save_checkpoint(raw, st, ps)
            archive_run(archive, st, ps; rng_seed = 17)
            copies = [fixture(T, Val(D); shape, ηH = 0.01) for _ = 1:2]
            load_checkpoint!(copies[1][1], copies[1][2], raw)
            load_archive!(copies[2][1], copies[2][2], archive)
            step!(st, ps, T(1 / 2048))
            for (restored, particles) in copies
                step!(restored, particles, T(1 / 2048))
                @test isequal(snapshot(restored, particles), snapshot(st, ps))
            end
        end
        @test all(
            i ->
                scratch[i] === (
                    st.pressure_work.quotient,
                    H._field_aligned_arrays(st.pressure_work.conduction)...,
                )[i],
            eachindex(scratch),
        )
        @test isfinite(
            energy_budget(
                ps,
                st.fields.B,
                st.fields.n,
                st.model.closure,
                st.g;
                pe = st.fields.pe,
            ).total,
        )
    end
end

@testset "Stage direction, exact reductions, rejection and zero-to-active restore" begin
    T = Float64
    for ηH in (0.0, 0.01), κ in (0.0, 0.0625)
        closure = FieldAlignedConductingElectrons(T(2); κ_parallel = κ, κ_perpendicular = κ)
        left, ps = fixture(T, Val(1); closure, ηH)
        legacy = iszero(κ) ? TransportElectrons(T(2)) : ConductingElectrons(T(2); κ)
        right, qs = fixture(T, Val(1); closure = legacy, ηH)
        init!(left, ps)
        init!(right, qs)
        step!(left, ps, 1 / 1024)
        step!(right, qs, 1 / 1024)
        @test isequal(snapshot(left, ps), snapshot(right, qs))
        @test right.pressure_work.conduction === nothing
        @test left.pressure_work.conduction !== nothing
    end
    st, ps = fixture(T, Val(1))
    init!(st, ps)
    st.fn .= st.fields.n
    foreach(c -> copyto!(st.fui[c], st.fields.ui[c]), 1:3)
    st.fields.ninv .= inv.(st.fn)
    orientations =
        ((fill(1.0, 4), fill(2.0, 4), fill(2.0, 4)), (fill(2.0, 4), fill(1.0, 4), fill(2.0, 4)))
    rates = map(orientations) do B
        H._bpressure_rhs!(st.k1, st.pressure_work.kp[1], B, st.fields.pe, st, 1, 1)
        copy(st.pressure_work.kp[1])
    end
    @test rates[1] != rates[2]
    original = snapshot(st, ps)
    for bad in (-1.0, Inf, NaN)
        @test_throws ArgumentError FieldAlignedConductingElectrons(
            2.0;
            κ_parallel = 1.0,
            κ_perpendicular = bad,
        )
        @test_throws ArgumentError FieldAlignedConductingElectrons(
            2.0;
            κ_parallel = bad,
            κ_perpendicular = 1.0,
        )
    end
    @test_throws ArgumentError HybridStepper(
        FourierGrid((4,), (4.0f0,)),
        HybridModel(
            FieldAlignedConductingElectrons(2.0; κ_parallel = 1.0, κ_perpendicular = 1e-100),
        ),
        NGP(),
        0,
    )
    for component in st.fields.B
        component[1] = 0
    end
    null_state = snapshot(st, ps)
    @test_throws ArgumentError step!(st, ps, 1 / 1024)
    @test isequal(snapshot(st, ps), null_state)
    foreach(c -> copyto!(st.fields.B[c], original.fields.B[c]), 1:3)
    @test_throws ArgumentError H._validate_evolving_pressure(fill(-1.0, 4), st.g)
    @test_throws ArgumentError step!(st, ps, 1e6)
    @test isequal(snapshot(st, ps), original)
    for (mode, flush) in ((RoundDown, false), (RoundNearest, true))
        saved = Base.Rounding.rounding_raw(Float64), get_zero_subnormals()
        try
            Base.Rounding.setrounding_raw(Float64, Base.Rounding.to_fenv(mode))
            Base.set_zero_subnormals(flush)
            @test_throws ArgumentError step!(st, ps, 1 / 1024)
        finally
            Base.Rounding.setrounding_raw(Float64, saved[1])
            Base.set_zero_subnormals(saved[2])
        end
        @test isequal(snapshot(st, ps), original)
    end
    old_work = st.pressure_work.conduction
    st.pressure_work.conduction =
        H._FieldAlignedConductionWorkspace{1}(st.pressure_work.quotient, old_work.gradient)
    @test_throws ArgumentError step!(st, ps, 1 / 1024)
    @test isequal(snapshot(st, ps), original)
    st.pressure_work.conduction = old_work
    step!(st, ps, 1 / 1024) # valid retry after the failed active call
    mktempdir() do dir
        path = joinpath(dir, "active.bin")
        save_checkpoint(path, st, ps)
        dst, qs =
            fixture(T, Val(1); closure = FieldAlignedConductingElectrons(2.0; κ_parallel = 0.0))
        identities =
            (dst.pressure_work.quotient, H._field_aligned_arrays(dst.pressure_work.conduction)...)
        load_checkpoint!(dst, qs, path)
        @test dst.model.closure.κ_parallel == st.model.closure.κ_parallel
        @test all(
            i ->
                identities[i] === (
                    dst.pressure_work.quotient,
                    H._field_aligned_arrays(dst.pressure_work.conduction)...,
                )[i],
            eachindex(identities),
        )
        state = deserialize(path)
        before = snapshot(dst, qs)
        for malformed in (
            merge(state, (conduction_policy_version = true,)),
            merge(state, (conduction_direction = :isotropic,)),
            merge(state, (B = ntuple(_ -> zeros(T, 4), 3),)),
        )
            serialize(path, malformed)
            @test_throws ArgumentError load_checkpoint!(dst, qs, path)
            @test isequal(snapshot(dst, qs), before)
        end
        step!(st, ps, 1 / 2048)
        step!(dst, qs, 1 / 2048)
        @test isequal(snapshot(st, ps), snapshot(dst, qs))
    end
    isotropic, qs = fixture(
        T,
        Val(1);
        closure = FieldAlignedConductingElectrons(2.0; κ_parallel = 0.1, κ_perpendicular = 0.1),
    )
    foreach(a -> fill!(a, 0), isotropic.fields.B)
    @test init!(isotropic, qs) === isotropic
    # Full spectral diffusion is not positivity preserving. The new anisotropic
    # operator must reject a negative RK trial and restore the complete state.
    g = FourierGrid((4,), (T(2pi),))
    ions = ParticleSet{1,T}(8)
    load_lattice_1d!(ions, zero(T), g.L[1])
    set_density_weight!(ions, one(T), g)
    negative = HybridStepper(
        g,
        HybridModel(
            FieldAlignedConductingElectrons(T(2); κ_parallel = T(0.1), κ_perpendicular = T(0.025)),
        ),
        NGP(),
        0,
    )
    negative.fields.B[1] .= one(T)
    negative.fields.pe .= T[0, 0, 1, 0]
    init!(negative, ions)
    before = snapshot(negative, ions)
    @test_throws ArgumentError step!(negative, ions, T(0.001))
    @test isequal(snapshot(negative, ions), before)
end

include("collisional_fixtures.jl")
@testset "Ion-ion composition, RNG rollback and archive continuation" begin
    for logs in (0.0, 0.05)
        closure =
            FieldAlignedConductingElectrons(2.0; κ_parallel = 0.0625, κ_perpendicular = 0.015625)
        st, ps = _composition_transport_fixture(Float64, Xoshiro; closure, logs)
        reference_st, reference_ps = _composition_transport_fixture(Float64, Xoshiro; closure, logs)
        init!(st, ps)
        init!(reference_st, reference_ps)
        step!(st, ps, st.dt)
        _composition_component!(reference_st, reference_ps)
        @test isequal(
            _composition_physical(st, ps),
            _composition_physical(reference_st, reference_ps),
        )
        before = _composition_snapshot(st, ps)
        @test electron_conduction_diagnostics(st, ps).internal_energy_rate == 0
        @test isequal(_composition_snapshot(st, ps), before)
        mktempdir() do dir
            path = joinpath(dir, "collision.hpr")
            archive_run(path, st, ps; rng_seed = 712)
            restored, qs = _composition_transport_fixture(Float64, Xoshiro; closure, logs)
            load_archive!(restored, qs, path)
            step!(st, ps, st.dt)
            step!(restored, qs, restored.dt)
            @test isequal(_composition_snapshot(st, ps), _composition_snapshot(restored, qs))
        end
        before = _composition_snapshot(st, ps)
        saved_model = st.inner.model
        # A changed fixed model rejects before advancing fields or RNG.
        st.inner.model = HybridModel(FieldAlignedConductingElectrons(2.0; κ_parallel = 1e20))
        @test_throws ArgumentError step!(st, ps, st.dt)
        st.inner.model = saved_model
        @test isequal(_composition_snapshot(st, ps), before)
        step!(st, ps, st.dt)
    end
    closure = FieldAlignedConductingElectrons(2.0; κ_parallel = 1024.0, κ_perpendicular = 256.0)
    st, ps = _composition_transport_fixture(Float64, Xoshiro; closure)
    init!(st, ps)
    before = _composition_snapshot(st, ps)
    @test_throws ArgumentError step!(st, ps, st.dt)
    @test isequal(_composition_snapshot(st, ps), before)
    step!(st, ps, st.dt; NB = 16)
    @test st.source_count == 1
end
end
