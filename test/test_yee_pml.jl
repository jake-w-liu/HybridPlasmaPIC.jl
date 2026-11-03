using HybridPlasmaPIC, Test

# YeePML must construct, initialize and advance a seeded pulse through its
# public operations in every supported dimension and precision.
function yee_pml_pulse(::Type{T}, D) where {T}
    n = ntuple(d -> 24 + 2d, D)
    g = YeeGrid(n, ntuple(d -> T(d + 1), D); boundaries = ntuple(_ -> :pec, D))
    vacuum = MaxwellVacuum(g; c = T(1))
    for c = 1:3, a in (vacuum.fields.E[c], vacuum.fields.B[c])
        for I in CartesianIndices(a)
            a[I] = T(sin(0.3 * sum(Tuple(I)) + c) * exp(-0.01 * sum(abs2, Tuple(I) .- 12)))
        end
    end
    apply_yee_pec!(vacuum.fields, g)
    return vacuum
end

@testset "YeePML constructs and advances" begin
    for T in (Float32, Float64), D = 1:3
        vacuum = yee_pml_pulse(T, D)
        h = T(0.4) * maxwell_cfl(vacuum)
        pml = YeePML(
            vacuum;
            dt = h,
            cells = ntuple(_ -> (4, 4), D),
            peak_rate = ntuple(_ -> (T(2), T(2)), D),
            profile = :quadratic,
        )
        init!(pml)
        for _ = 1:12
            step_maxwell!(pml, h)
        end
        @test all(a -> all(isfinite, a), (vacuum.fields.E..., vacuum.fields.B...))
        @test pml_diagnostics(pml) !== nothing
    end
end
