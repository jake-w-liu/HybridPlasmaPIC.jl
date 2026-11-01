@testset "Hyperresistive singleton and mixed-axis Nyquist modes" begin
    for T in (Float32, Float64), D = 1:3, axis = 1:D
        # In D1 both singleton and pure-Nyquist current vanish. In higher D,
        # keep a separate resolved direction so every axis has a nonzero channel.
        resolved = D == 1 ? 1 : mod1(axis + 1, D)
        for singleton in (false, true)
            n = ntuple(d -> d == axis ? (singleton ? 1 : 4) : (d == resolved ? 5 : 1), D)
            g = FourierGrid(n, ntuple(_ -> T(2pi), D))
            st = HybridStepper(
                g,
                HybridModel(TransportElectrons(T(5 / 3)); η = 0.02, ηH = 0.03),
                CIC(),
                0,
            )
            f = st.fields
            fill!(f.n, one(T))
            fill!(f.ninv, one(T))
            fill!(st.fn, one(T))
            fill!(f.pe, one(T))
            component = resolved == 3 ? 2 : 3
            for I in CartesianIndices(n)
                phase = T(I[resolved] - 1) * g.dx[resolved]
                f.B[component][I] =
                    T(0.5) * (-one(T))^(I[axis] - 1) * (D == 1 ? one(T) : sin(phase))
            end
            reference = th_reference(f.B, f.pe, f.n, st.fui, g.L, T(5 / 3), 0.02, 0.03)
            HybridPlasmaPIC._bpressure_rhs!(st.k1, st.pressure_work.kp[1], f.B, f.pe, st, 1, 1)
            for c = 1:3
                @test maximum(abs, st.k1[c] .- reference.outB[c]) <=
                      512eps(T) * max(1.0, maximum(abs, reference.outB[c]))
            end
            @test maximum(abs, st.pressure_work.kp[1] .- reference.outp) <=
                  512eps(T) * max(1.0, maximum(abs, reference.outp))
            if D == 1
                @test all(iszero, st.pressure_work.kp[1])
            elseif !singleton
                @test sum(reference.nyquist) > 1
            end
        end
    end
end
