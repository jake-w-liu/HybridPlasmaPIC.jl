using HybridPlasmaPIC, Test

_tuple_maxabs(a, b) = maximum(abs.(a .- b))

@testset "fusion: toroidal curvilinear mesh" begin
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 2, 8, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 8, 2, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 8, 8, 2)

    g = ToroidalGrid(3.0, 1.0, 8, 10, 12)
    @test gridsize(g) == (8, 10, 12)
    @test scale_factors(g, 2, 3) == (1.0, g.r[2], 3.0 + g.r[2] * cos(g.θ[3]))
    @test jacobian(g, 2, 3) ≈ g.r[2] * (3.0 + g.r[2] * cos(g.θ[3]))
    @test _tuple_maxabs(to_cartesian(g, 0.5, π / 2, 0.0), (3.0, 0.0, 0.5)) < 1e-14

    f = [g.r[i]^2 * cos(g.θ[j]) * sin(g.φ[k]) for i = 1:8, j = 1:10, k = 1:12]
    gr, gθ, gφ = metric_gradient(g, f)
    @test size(gr) == gridsize(g)
    @test size(gθ) == gridsize(g)
    @test size(gφ) == gridsize(g)

    div = metric_divergence(g, ones(gridsize(g)), zeros(gridsize(g)), zeros(gridsize(g)))
    @test size(div) == gridsize(g)
    @test all(isfinite, div)

    ffun(r, θ, φ) = r^2 * cos(θ) * sin(φ)
    grad_an(r, θ, φ, R) = (2r * cos(θ) * sin(φ), -r * sin(θ) * sin(φ), (r^2 * cos(θ) * cos(φ)) / R)
    function max_gradient_error(N)
        gg = ToroidalGrid(3.0, 1.0, N, N, N)
        Nr, Nθ, Nφ = gridsize(gg)
        ff = [ffun(gg.r[i], gg.θ[j], gg.φ[k]) for i = 1:Nr, j = 1:Nθ, k = 1:Nφ]
        fr, fθ, fφ = metric_gradient(gg, ff)
        e = 0.0
        for k = 1:Nφ, j = 1:Nθ, i = 2:Nr-1
            R = gg.R0 + gg.r[i] * cos(gg.θ[j])
            ga = grad_an(gg.r[i], gg.θ[j], gg.φ[k], R)
            e = max(e, abs(fr[i, j, k] - ga[1]), abs(fθ[i, j, k] - ga[2]), abs(fφ[i, j, k] - ga[3]))
        end
        return e
    end
    e1 = max_gradient_error(24)
    e2 = max_gradient_error(48)
    @test e1 < 0.05
    @test e2 < e1 / 3

    radial = ToroidalGrid(3.0, 1.0, 20, 24, 16)
    Nr, Nθ, Nφ = gridsize(radial)
    Ar = [sin(radial.φ[k]) for i = 1:Nr, j = 1:Nθ, k = 1:Nφ]
    z = zeros(Nr, Nθ, Nφ)
    d = metric_divergence(radial, Ar, z, z)
    e_exact = 0.0
    for k = 1:Nφ, j = 1:Nθ, i = 2:Nr-1
        r = radial.r[i]
        R = radial.R0 + r * cos(radial.θ[j])
        da = sin(radial.φ[k]) * (radial.R0 + 2r * cos(radial.θ[j])) / (r * R)
        e_exact = max(e_exact, abs(d[i, j, k] - da))
    end
    @test e_exact < 1e-10
end

@testset "fusion: guiding-centre gyrokinetics" begin
    @test _tuple_maxabs(exb_drift((0.0, 1.0, 0.0), (0.0, 0.0, 2.0)), (0.5, 0.0, 0.0)) < 1e-14
    @test_throws ArgumentError exb_drift((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))

    vg = gradb_drift(1.0, 1.0, 1.0, (0.0, 0.0, 2.0), (3.0, 0.0, 0.0))
    @test _tuple_maxabs(vg, (0.0, 6 / 16, 0.0)) < 1e-14

    vc = curvature_drift(1.0, 1.0, 1.0, (0.0, 0.0, 2.0), (0.5, 0.0, 0.0))
    @test _tuple_maxabs(vc, (0.0, 0.25, 0.0)) < 1e-14

    gc = GuidingCentre((0.0, 0.0, 0.0), 0.5, 0.0, 1.0, 1.0)
    push_guiding_centre!(
        gc;
        dt = 0.2,
        E = (0.0, 0.0, 1.0),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = 0.0,
    )
    @test _tuple_maxabs(gc.X, (0.0, 0.0, 0.1)) < 1e-14
    @test gc.vpar ≈ 0.7
    @test gyroaverage(x -> x[1]^2 + x[2]^2 + x[3]^2, (1.0, 2.0, 3.0), 0.5, (0.0, 0.0, 1.0); n = 32) ≈
          14.25
    @test_throws ArgumentError gyroaverage(x -> x[1], (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); n = 2)

    B = (0.3, -1.2, 2.0)
    vtot = drift_velocity(;
        vpar = 0.7,
        vperp = 1.1,
        q = 1.0,
        m = 1.0,
        E = (0.5, 0.2, -0.1),
        B = B,
        gradB = (0.4, 0.0, 0.3),
        κ = (0.1, -0.2, 0.05),
    )
    @test abs(vtot[1] * B[1] + vtot[2] * B[2] + vtot[3] * B[3]) < 1e-12

    gc2 = GuidingCentre((0.0, 0.0, 0.0), 0.5, 0.5, 1.0, 1.0)
    for _ = 1:100
        push_guiding_centre!(
            gc2;
            dt = 0.01,
            E = (0.0, 1.0, 0.0),
            B = (0.0, 0.0, 1.0),
            gradB = (0.0, 0.0, 0.0),
            κ = (0.0, 0.0, 0.0),
            gradpar_B = 0.0,
        )
    end
    @test gc2.X[1] ≈ 1.0 rtol = 1e-12
    @test abs(gc2.X[2]) < 1e-12
    @test gc2.X[3] ≈ 0.5 rtol = 1e-12
    @test gc2.vpar ≈ 0.5
end

@testset "fusion: adaptive mesh refinement" begin
    coarse = AMRGrid([1.0, 2.0, 3.0], 0.25; x0 = -0.5)
    fine = refine(coarse)
    @test ncells(fine) == 2ncells(coarse)
    @test effective_resolution(fine) ≈ 0.125
    @test cell_center(coarse, 1) ≈ -0.375
    @test !any(refine_flags(coarse.u, 1e-12))
    @test_throws ArgumentError refine_flags(coarse.u, -1.0)

    back = AMRGrid(similar(coarse.u), coarse.dx; x0 = coarse.x0)
    restrict!(back, fine)
    @test back.u ≈ coarse.u

    single = refine(AMRGrid([3.0], 1.0))
    @test single.u == [3.0, 3.0]

    dx = 0.2
    x0 = 0.5
    a, b = 3.0, -1.7
    linear = AMRGrid([a + b * (x0 + (i - 0.5) * dx) for i = 1:12], dx; x0 = x0)
    linear_fine = refine(linear)
    for i = 3:2*12-2
        @test linear_fine.u[i] ≈ a + b * cell_center(linear_fine, i) rtol = 1e-13
    end

    mass_coarse = AMRGrid(Float64[exp(-((i - 8.0) / 3)^2) for i = 1:16], 0.1)
    mass_fine = refine(mass_coarse)
    @test sum(mass_fine.u) * mass_fine.dx ≈ sum(mass_coarse.u) * mass_coarse.dx rtol = 1e-13
end
