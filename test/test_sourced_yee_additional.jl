using HybridPlasmaPIC, Test
@testset "Sourced final magnetic failure is transactional" begin
    for T in (Float32, Float64), n in (2, 4, 8), sign in (-1, 1)
        g = YeeGrid((n,), (T(n),); boundaries = (:periodic,))
        m = MaxwellVacuum(g)
        a = T(sign) * T(0.9) * floatmax(T)
        m.fields.B[3] .= a
        j = map(x -> zero(x), m.fields.E)
        for i in eachindex(j[2])
            j[2][i] = isodd(i) ? a : -a
        end
        old = map(copy, (m.fields.E..., m.fields.B...))
        source = map(copy, j)
        @test_throws ArgumentError step_maxwell!(m, T(0.5), j)
        @test all(isequal.(old, (m.fields.E..., m.fields.B...)))
        @test j == source && m.time[] == zero(T) && m.step[] == 0
        # Demonstrate that the valid full E stage had already completed; this
        # distinguishes final-B failure from prevalidation/initial-B rejection.
        @test m._candidate.E[2] == -j[2]
        @test all(isfinite, m._candidate.E[2])
    end
end

# Independent closed form of the forced discrete modal map. This checks both
# transverse components and both precisions without a continuum truncation error.
@testset "Sourced discrete forced modes in both precisions" begin
    for T in (Float32, Float64),
        p in (false, true),
        component in (2, 3),
        mode in (1, 2, 3),
        sgn in (-1, 1)

        n = 16
        steps = 64
        g = YeeGrid((n,), (one(T),); boundaries = (p ? :periodic : :pec,))
        m = MaxwellVacuum(g; c = T(1.25))
        h = maxwell_cfl(m) / T(4)
        jamp = T(sgn) * T(0.375)
        k = (p ? 2pi : pi) * mode
        kappa = 2sin(k * Float64(g.dx[1]) / 2) / Float64(g.dx[1])
        j = map(x -> zero(x), m.fields.E)
        x = only(yee_coordinates(g, :electric, component))
        j[component] .= h * jamp .* T.(sin.(k .* Float64.(x)))
        if !p
            j[component][1] = j[component][end] = zero(T)
        end
        for _ = 1:steps
            step_maxwell!(m, h, j)
        end
        theta = 2asin(Float64(m.c) * kappa * Float64(h) / 2)
        e = -Float64(h * jamp) * sin(steps * theta) / sin(theta)
        b = Float64(jamp) * (1 - cos(steps * theta)) / (Float64(m.c)^2 * kappa)
        magnetic = component == 2 ? 3 : 2
        orientation = component == 2 ? 1 : -1
        exact_e = e .* sin.(k .* Float64.(x))
        xb = only(yee_coordinates(g, :magnetic, magnetic))
        exact_b = orientation * b .* cos.(k .* Float64.(xb))
        @test maximum(abs, Float64.(m.fields.E[component]) - exact_e) <=
              256eps(T) * max(1.0, abs(e))
        @test maximum(abs, Float64.(m.fields.B[magnetic]) - exact_b) <= 256eps(T) * max(1.0, abs(b))
    end
end
