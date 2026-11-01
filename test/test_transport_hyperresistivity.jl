using HybridPlasmaPIC, Test, LinearAlgebra, Serialization
include("transport_hyperresistivity/reference.jl")

@testset "Hyperresistive model contract" begin
    for T in (Float32, Float64)
        g = FourierGrid((8,), (T(2pi),))
        closure = TransportElectrons(T(5 / 3))
        @test HybridStepper(g, HybridModel(closure; η = 0.02, ηH = 0.03), CIC(), 0) isa
              HybridStepper
        for bad in (-0.01, Inf, NaN)
            @test_throws ArgumentError HybridModel(closure; ηH = bad)
        end
        @test_throws ArgumentError HybridModel(closure; ηH = 0.03, de2 = 0.01)
        @test_throws ArgumentError HybridModel(EvolvingElectrons(); ηH = 0.03)
        for constructor in (CAMCLStepper,)
            @test_throws ArgumentError constructor(g, HybridModel(closure; ηH = 0.03), CIC(), 0)
        end
        @test_throws ArgumentError HallMHDModel(closure; ηH = 0.03)
        st, ps = tr_fixture(T, Val(1), closure; ηH = 0.03)
        saved = map(copy, tr_state(st))
        clock = tr_clock(st)
        particles = tr_particles(ps)
        @test step!(st, ps, zero(T)) === st
        @test all(isequal(a, b) for (a, b) in zip(saved, tr_state(st)))
        @test tr_clock(st) == clock && tr_same_particles(ps, particles)
    end
    g = FourierGrid((8,), (Float32(2pi),))
    for bad in (1e-50, 1e300)
        @test_throws ArgumentError HybridStepper(
            g,
            HybridModel(TransportElectrons(); ηH = bad),
            CIC(),
            0,
        )
    end
    @test HybridStepper(g, HybridModel(TransportElectrons(); ηH = 1e-40), CIC(), 0) isa
          HybridStepper
end

@testset "Independent full Fourier hyperresistive stages and SBP" begin
    for T in (Float32, Float64), D = 1:3, even in (false, true), η in (0.0, 0.03)
        n = ntuple(d -> even ? 6 : 5, D)
        L = ntuple(d -> T(2pi), D)
        g = FourierGrid(n, L)
        st = HybridStepper(g, HybridModel(TransportElectrons(T(5 / 3)); η, ηH = 0.025), CIC(), 0)
        f = st.fields
        for I in CartesianIndices(n)
            phase = sum(T(I[d] - 1) * g.dx[d] for d = 1:D)
            f.pe[I] = T(0.7) + T(0.02) * cos(phase)
            f.n[I] = T(1.1) + T(0.03) * sin(phase)
            f.ninv[I] = inv(f.n[I])
            st.fn[I] = f.n[I]
            for c = 1:3
                f.B[c][I] = T(0.3 / c) * sin(phase + T(c))
                st.fui[c][I] = T(0.02 / c) * cos(phase + T(c / 3))
            end
            if even && D >= 2
                f.B[3][I] += T(0.6) * (-one(T))^(I[1] - 1) * sin(T(I[2] - 1) * g.dx[2])
            end
        end
        ref = th_reference(f.B, f.pe, f.n, st.fui, L, T(5 / 3), η, 0.025)
        HybridPlasmaPIC._bpressure_rhs!(st.k1, st.pressure_work.kp[1], f.B, f.pe, st, 1, 1)
        for (actual, expected) in
            ((f.J, ref.J), (f.lapJ, ref.LJ), (st.Escr, ref.E), (st.k1, ref.outB))
            for c = 1:3
                @test maximum(abs, actual[c] .- expected[c]) <=
                      512eps(T) * max(1.0, maximum(abs, expected[c]))
            end
        end
        @test maximum(abs, st.pressure_work.kp[1] .- ref.outp) <=
              512eps(T) * max(1.0, maximum(abs, ref.outp))
        # Isolate the dissipative exchange from ideal flow/pressure work.
        heat = zeros(T, n)
        HybridPlasmaPIC._add_hyperresistive_pressure_heating!(
            heat,
            f.J,
            f.gradp,
            st.pressure_work.divue,
            T(5 / 3),
            T(0.025),
            g,
        )
        loss = -T(0.025) * sum(dot(vec(f.J[c]), vec(f.lapJ[c])) for c = 1:3)
        @test minimum(heat) >= 0
        @test isapprox(sum(heat) / (T(5 / 3) - 1), loss; rtol = 512eps(T), atol = 512eps(T))
        if even && D >= 2
            # Reject omission of the Nyquist channel and reversal of the field sign.
            missing = 0.025sum(ref.nyquist)
            @test missing > 1
            @test abs(loss - 0.025sum(ref.ordinary)) > missing / 2
            @test abs(-loss - sum(heat) / (T(5 / 3) - 1)) > loss
        end
        # A second, distinct stage rejects stale-current/heat implementations.
        nextB = map(a -> T(0.8) .* a, f.B)
        nextp = T(0.9) .* f.pe
        second = th_reference(nextB, nextp, f.n, st.fui, L, T(5 / 3), η, 0.025)
        HybridPlasmaPIC._bpressure_rhs!(st.k2, st.pressure_work.kp[2], nextB, nextp, st, 1, 2)
        @test maximum(abs, st.pressure_work.kp[2] .- second.outp) <=
              512eps(T) * max(1.0, maximum(abs, second.outp))
        for c = 1:3
            @test maximum(abs, st.k2[c] .- second.outB[c]) <=
                  512eps(T) * max(1.0, maximum(abs, second.outB[c]))
        end
        # Direct/public carried field uses the same split and same full Laplacian.
        copyto!.(f.ui, st.fui)
        ohms_law!(f, st.model, g)
        for c = 1:3
            @test maximum(abs, f.E[c] .- ref.E[c]) <= 512eps(T) * max(1.0, maximum(abs, ref.E[c]))
        end
    end
end

@testset "Hyperresistive analytic field/heat and energy refinement" begin
    for T in (Float32, Float64), D = 1:3, NB in (1, 2), η in (0.0, 0.01)
        r = tr_beltrami(T, Val(D), η, T(0.04), 8; NB, ηH = 0.2 / D^2)
        tolerance = T === Float32 ? 128eps(T) : 1e-10
        @test r.field < tolerance
        @test r.pressure < tolerance
        @test r.energy < tolerance * r.total
        @test r.speed < 64eps(T)
        @test r.force < 128eps(T)
        @test r.split
        @test r.heated > 0.003
    end
    results = [tr_beltrami(Float64, Val(2), 0.05, 0.4 / k, k; ηH = 0.1) for k in (2, 4, 8)]
    for name in (:field, :pressure, :energy)
        errors = [getproperty(r, name) for r in results]
        @test all(12 .< errors[1:2] ./ errors[2:3] .< 20)
    end
    for multispecies in (false, true)
        errors =
            [tr_energy_error(k; multispecies, η = 0.02, ηH = 0.03, NB = 128) for k in (8, 16, 32)]
        @test all(3 .< errors[1:2] ./ errors[2:3] .< 5)
    end
end

@testset "Hyperresistive rollback and restart" begin
    for T in (Float32, Float64), D = 1:3, collection in (false, true)
        source, ps =
            tr_fixture(T, Val(D), TransportElectrons(T(5 / 3)); η = 0.02, ηH = 0.03, collection)
        step!(source, ps, T(0.002); NB = 2)
        # Failed direct density, failed step model and late pressure rejection.
        saved = source.model
        source.model = HybridModel(saved.closure, saved.η, -0.1, saved.nfloor, saved.de2)
        tr_rejected(() -> step!(source, ps, T(0.002)), source, ps, "ηH")
        source.model = saved
        target, restored =
            tr_fixture(T, Val(D), TransportElectrons(T(5 / 3)); η = 0.8, ηH = 0.9, collection)
        mktempdir() do dir
            path = joinpath(dir, "hyper.ser")
            save_checkpoint(path, source, ps)
            raw = deserialize(path)
            @test raw.transport_policy_version === 1
            @test raw.model.ηH === 0.03
            load_checkpoint!(target, restored, path)
            @test target.model.ηH === 0.03 && target.model.η === 0.02
            @test isequal(HybridPlasmaPIC._checkpoint_state(target, restored), raw)
            for (dt, NB) in ((T(0.003), 1), (T(0.001), 3))
                step!(source, ps, dt; NB)
                step!(target, restored, dt; NB)
                @test all(isequal(a, b) for (a, b) in zip(tr_state(source), tr_state(target)))
                @test tr_clock(source) == tr_clock(target)
                @test tr_same_particles(ps, tr_particles(restored))
            end
            archive = joinpath(dir, "hyper.arc")
            archive_run(archive, source, ps; rng_seed = 711)
            archived = load_archive(archive).state
            HybridPlasmaPIC._restore_checkpoint_state!(target, restored, archived)
            @test isequal(HybridPlasmaPIC._checkpoint_state(target, restored), archived)
            step!(source, ps, T(0.001))
            step!(target, restored, T(0.001))
            @test isequal(
                transport_io_public_state(source, ps),
                transport_io_public_state(target, restored),
            )
        end
    end
    # Positive sources cannot conceal loss of pressure admissibility.
    for T in (Float32, Float64)
        st, ps = tr_fixture(T, Val(1), TransportElectrons(T(3)); η = 0.02, ηH = 0.03)
        @. ps.v[1] = sin(ps.x[1])
        init!(st, ps)
        tr_rejected(() -> step!(st, ps, T(0.8); NB = 2), st, ps, "subcycle")
    end
end

@testset "Final carried fields after hyperresistive scratch reuse" begin
    for T in (Float32, Float64), D = 1:3, collection in (false, true)
        st, ps =
            tr_fixture(T, Val(D), TransportElectrons(T(5 / 3)); η = 0.02, ηH = 0.03, collection)
        step!(st, ps, T(0.001); NB = 2)
        f = st.fields
        reference = th_reference(f.B, f.pe, f.n, f.ui, st.g.L, T(5 / 3), 0.02, 0.03)
        for (actual, expected) in (
            (st.pressure_work.Eforce, reference.force),
            (f.E, reference.E),
            (f.J, reference.J),
            (f.lapJ, reference.LJ),
        )
            for c = 1:3
                @test maximum(abs, actual[c] .- expected[c]) <=
                      512eps(T) * max(1.0, maximum(abs, expected[c]))
            end
        end
        for d = 1:D
            expected = th_apply(th_matrix(st.g.n[d], st.g.L[d], :first), f.pe, d)
            @test maximum(abs, f.gradp[d] .- expected) <=
                  512eps(T) * max(1.0, maximum(abs, expected))
        end
    end
end
