using HybridPlasmaPIC, Test

_tuple_maxabs(a, b) = maximum(abs.(a .- b))

struct ZeroBasedCube{T} <: AbstractArray{T,3}
    parent::Array{T,3}
end
Base.size(a::ZeroBasedCube) = size(a.parent)
Base.axes(a::ZeroBasedCube) = ntuple(d -> 0:(size(a, d)-1), 3)
Base.getindex(a::ZeroBasedCube, i::Int, j::Int, k::Int) = a.parent[i+1, j+1, k+1]

@testset "fusion: toroidal curvilinear mesh" begin
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 2, 8, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 8, 2, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 8, 8, 2)
    @test_throws ArgumentError ToroidalGrid(Inf, 1.0, 8, 8, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, Inf, 8, 8, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 0.0, 8, 8, 8)
    @test_throws ArgumentError ToroidalGrid(1.0, 1.0, 8, 8, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, big(typemax(Int)) + 1, 8, 8)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 8, 8, 8; T = Int)
    @test_throws ArgumentError ToroidalGrid(3.0, 1.0, 8, 8, 8; T = AbstractFloat)

    g = ToroidalGrid(3.0, 1.0, 8, 10, 12)
    @test gridsize(g) == (8, 10, 12)
    @test scale_factors(g, 2, 3) == (1.0, g.r[2], 3.0 + g.r[2] * cos(g.θ[3]))
    @test jacobian(g, 2, 3) ≈ g.r[2] * (3.0 + g.r[2] * cos(g.θ[3]))
    @test _tuple_maxabs(to_cartesian(g, 0.5, π / 2, 0.0), (3.0, 0.0, 0.5)) < 1e-14
    @test_throws ArgumentError to_cartesian(g, Inf, 0.0, 0.0)
    @test_throws ArgumentError to_cartesian(g, 0.5, NaN, 0.0)
    @test_throws ArgumentError to_cartesian(g, 0.5, 0.0, Inf)

    g32 = ToroidalGrid(3.0, 1.0, 8, 8, 8; T = Float32)
    @test g32 isa ToroidalGrid{Float32}
    @test all(isfinite, (g32.R0, g32.dr, g32.dθ, g32.dφ))

    f = [g.r[i]^2 * cos(g.θ[j]) * sin(g.φ[k]) for i = 1:8, j = 1:10, k = 1:12]
    gr, gθ, gφ = metric_gradient(g, f)
    @test size(gr) == gridsize(g)
    @test size(gθ) == gridsize(g)
    @test size(gφ) == gridsize(g)

    div = metric_divergence(g, ones(gridsize(g)), zeros(gridsize(g)), zeros(gridsize(g)))
    @test size(div) == gridsize(g)
    @test all(isfinite, div)

    large = ToroidalGrid(2.0e200, 1.0e200, 8, 8, 8)
    large_ones = ones(gridsize(large))
    large_zero = zeros(gridsize(large))
    large_divergence = metric_divergence(large, large_ones, large_zero, large_zero)
    large_exact = [
        inv(large.r[i]) + cos(large.θ[j]) / (large.R0 + large.r[i] * cos(large.θ[j])) for
        i = 1:8, j = 1:8, k = 1:8
    ]
    @test all(isfinite, large_divergence)
    @test maximum(abs.((large_divergence .- large_exact) ./ large_exact)) < 1.0e-12

    large_Aφ = [sin(large.φ[k]) for i = 1:8, j = 1:8, k = 1:8]
    large_poloidal = metric_divergence(large, large_zero, large_ones, large_zero)
    large_toroidal = metric_divergence(large, large_zero, large_zero, large_Aφ)
    poloidal_exact = [
        -sin(large.θ[j]) * (sin(large.dθ) / large.dθ) / (large.R0 + large.r[i] * cos(large.θ[j])) for i = 1:8, j = 1:8, k = 1:8
    ]
    toroidal_exact = [
        cos(large.φ[k]) * (sin(large.dφ) / large.dφ) / (large.R0 + large.r[i] * cos(large.θ[j])) for i = 1:8, j = 1:8, k = 1:8
    ]
    @test large_poloidal ≈ poloidal_exact rtol = 5.0e-14 atol = 1.0e-215
    @test large_toroidal ≈ toroidal_exact rtol = 5.0e-14 atol = 1.0e-215

    offset = ZeroBasedCube(ones(gridsize(g)))
    @test_throws ArgumentError metric_gradient(g, offset)
    @test_throws ArgumentError metric_divergence(g, offset, offset, offset)

    alloc_grid = ToroidalGrid(3.0, 1.0, 20, 24, 16)
    alloc_Ar = rand(gridsize(alloc_grid)...)
    alloc_zero = zeros(gridsize(alloc_grid))
    metric_divergence(alloc_grid, alloc_Ar, alloc_zero, alloc_zero)
    divergence_bytes = @allocated metric_divergence(alloc_grid, alloc_Ar, alloc_zero, alloc_zero)
    @test divergence_bytes < 2sizeof(alloc_Ar)

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

    # The metric flux Fr=r*R*Ar vanishes at the coordinate axis.  For
    # Ar=r^2/R it is exactly r^3, so the axis-aware inner stencil must recover
    # dFr/dr=3r^2 without the sign reversal of a regular one-sided stencil.
    for Naxis in (3, 4, 5, 16)
        axis = ToroidalGrid(3.0, 1.0, Naxis, Naxis, Naxis)
        Nr, Nθ, Nφ = gridsize(axis)
        axis_Ar =
            [axis.r[i]^2 / (axis.R0 + axis.r[i] * cos(axis.θ[j])) for i = 1:Nr, j = 1:Nθ, k = 1:Nφ]
        axis_zero = zeros(Nr, Nθ, Nφ)
        axis_div = metric_divergence(axis, axis_Ar, axis_zero, axis_zero)
        for i = 1:Nr-1
            R = axis.R0 + axis.r[i] * cos(axis.θ[2])
            @test axis_div[i, 2, 1] ≈ 3axis.r[i] / R rtol = 1e-13 atol = 1e-14
        end
    end

    function max_divergence_error(N)
        gg = ToroidalGrid(3.0, 1.0, N, N, N)
        Nr, Nθ, Nφ = gridsize(gg)
        Ar = Array{Float64}(undef, Nr, Nθ, Nφ)
        Aθ = similar(Ar)
        Aφ = similar(Ar)
        exact = similar(Ar)
        for k = 1:Nφ, j = 1:Nθ, i = 1:Nr
            r = gg.r[i]
            θ = gg.θ[j]
            φ = gg.φ[k]
            R = gg.R0 + r * cos(θ)
            Ar[i, j, k] = r^2 * sin(θ) * cos(φ)
            Aθ[i, j, k] = r * cos(θ) * sin(φ)
            Aφ[i, j, k] = r * sin(θ) * sin(φ)
            dFr = (3r^2 * R + r^3 * cos(θ)) * sin(θ) * cos(φ)
            dFθ = (-r^2 * sin(θ) * cos(θ) - R * r * sin(θ)) * sin(φ)
            dFφ = r^2 * sin(θ) * cos(φ)
            exact[i, j, k] = (dFr + dFθ + dFφ) / (r * R)
        end
        computed = metric_divergence(gg, Ar, Aθ, Aφ)
        return maximum(abs.(computed .- exact))
    end
    div_e16 = max_divergence_error(16)
    div_e32 = max_divergence_error(32)
    @test div_e32 < div_e16 / 3
end

@testset "fusion: guiding-centre gyrokinetics" begin
    @test _tuple_maxabs(exb_drift((0.0, 1.0, 0.0), (0.0, 0.0, 2.0)), (0.5, 0.0, 0.0)) < 1e-14
    @test_throws ArgumentError exb_drift((1.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    @test_throws ArgumentError exb_drift((1.0, 0.0, 0.0), (0.0, 0.0, Inf))
    @test_throws DimensionMismatch exb_drift((1.0, 0.0), (0.0, 0.0, 1.0))

    vg = gradb_drift(1.0, 1.0, 1.0, (0.0, 0.0, 2.0), (3.0, 0.0, 0.0))
    @test _tuple_maxabs(vg, (0.0, 6 / 16, 0.0)) < 1e-14
    @test_throws ArgumentError gradb_drift(-1.0, 1.0, 1.0, (0.0, 0.0, 2.0), (3.0, 0.0, 0.0))
    @test_throws ArgumentError gradb_drift(1.0, 0.0, 1.0, (0.0, 0.0, 2.0), (3.0, 0.0, 0.0))
    @test_throws ArgumentError gradb_drift(1.0, 1.0, -1.0, (0.0, 0.0, 2.0), (3.0, 0.0, 0.0))

    vc = curvature_drift(1.0, 1.0, 1.0, (0.0, 0.0, 2.0), (0.5, 0.0, 0.0))
    @test _tuple_maxabs(vc, (0.0, 0.25, 0.0)) < 1e-14
    @test_throws ArgumentError curvature_drift(Inf, 1.0, 1.0, (0.0, 0.0, 2.0), (0.5, 0.0, 0.0))

    @testset "scale-safe drift arithmetic" begin
        for scale in (1e-200, 1e200)
            Eextreme = (0.0, scale, 0.0)
            Bextreme = (0.0, 0.0, scale)
            grad_extreme = (scale, 0.0, 0.0)
            κextreme = (scale, 0.0, 0.0)
            @test collect(exb_drift(Eextreme, Bextreme)) ≈ [1.0, 0.0, 0.0] rtol = 8eps()
            @test collect(gradb_drift(1.0, 1.0, 1.0, Bextreme, grad_extreme)) ≈
                  [0.0, inv(2scale), 0.0] rtol = 16eps()
            @test collect(curvature_drift(1.0, 1.0, 1.0, Bextreme, κextreme)) ≈ [0.0, 1.0, 0.0] rtol =
                16eps()
        end

        M = floatmax(Float64)
        cancellation_args = (;
            vpar = 0.5,
            vperp = 1.0,
            q = 1.0,
            m = 2.0,
            E = (0.0, 0.75M, 0.0),
            B = (0.0, 0.0, 1.0),
            gradB = (0.0, -0.75M, 0.0),
            κ = (0.0, M, 0.0),
        )
        cancellation_total = drift_velocity(; cancellation_args...)
        cancellation_oracle = Float64(
            setprecision(2048) do
                BigFloat(cancellation_args.E[2]) - BigFloat(cancellation_args.gradB[2]) -
                BigFloat(cancellation_args.κ[2]) / 2
            end,
        )
        @test all(isfinite, cancellation_total)
        @test cancellation_total == (cancellation_oracle, 0.0, 0.0)

        ratio_numerators = (-7.392147589140323e-207, -5.024715391148077e126)
        ratio_denominators = (5.97155727804291e-38, -1.4093084331804231e281)
        ratio_oracle = Float64(
            (Rational{BigInt}(ratio_numerators[1]) * Rational{BigInt}(ratio_numerators[2])) /
            (Rational{BigInt}(ratio_denominators[1]) * Rational{BigInt}(ratio_denominators[2])),
        )
        @test ratio_oracle == -nextfloat(0.0)
        @test HybridPlasmaPIC._gyro_scaled_ratio(ratio_numerators, ratio_denominators) ==
              ratio_oracle

        Ealloc = (0.5, 0.2, -0.1)
        Balloc = (0.3, -1.2, 2.0)
        grad_alloc = (0.4, 0.0, 0.3)
        κalloc = (0.1, -0.2, 0.05)
        drift_args = (;
            vpar = 0.7,
            vperp = 1.1,
            q = 1.0,
            m = 1.0,
            E = Ealloc,
            B = Balloc,
            gradB = grad_alloc,
            κ = κalloc,
        )
        exb_drift(Ealloc, Balloc)
        gradb_drift(drift_args.vperp, drift_args.q, drift_args.m, Balloc, grad_alloc)
        curvature_drift(drift_args.vpar, drift_args.q, drift_args.m, Balloc, κalloc)
        drift_velocity(; drift_args...)
        @test (@allocated exb_drift(Ealloc, Balloc)) == 0
        @test (@allocated gradb_drift(
            drift_args.vperp,
            drift_args.q,
            drift_args.m,
            Balloc,
            grad_alloc,
        )) == 0
        @test (@allocated curvature_drift(
            drift_args.vpar,
            drift_args.q,
            drift_args.m,
            Balloc,
            κalloc,
        )) == 0
        @test (@allocated drift_velocity(; drift_args...)) == 0
    end

    @test_throws ArgumentError GuidingCentre((0.0, 0.0, 0.0), Inf, 0.0, 1.0, 1.0)
    @test_throws ArgumentError GuidingCentre((0.0, 0.0, 0.0), 0.0, -1.0, 1.0, 1.0)
    @test_throws ArgumentError GuidingCentre((0.0, 0.0, 0.0), 0.0, 0.0, 0.0, 1.0)
    @test_throws ArgumentError GuidingCentre((0.0, 0.0, 0.0), 0.0, 0.0, 1.0, -1.0)
    gc_int = GuidingCentre((0, 0, 0), 0, 0, 1, 1)
    @test gc_int isa GuidingCentre{Float64}

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
    gc_before = (gc.X, gc.vpar, gc.μ, gc.q, gc.m)
    @test_throws ArgumentError push_guiding_centre!(
        gc;
        dt = Inf,
        E = (0.0, 0.0, 0.0),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = 0.0,
    )
    @test (gc.X, gc.vpar, gc.μ, gc.q, gc.m) == gc_before
    @test_throws ArgumentError push_guiding_centre!(
        gc;
        dt = 0.1,
        E = (Inf, 0.0, 0.0),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = 0.0,
    )
    @test (gc.X, gc.vpar, gc.μ, gc.q, gc.m) == gc_before
    @test gyroaverage(
        x -> x[1]^2 + x[2]^2 + x[3]^2,
        (1.0, 2.0, 3.0),
        0.5,
        (0.0, 0.0, 1.0);
        n = 32,
    ) ≈ 14.25
    @test_throws ArgumentError gyroaverage(x -> x[1], (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); n = 2)
    @test_throws ArgumentError gyroaverage(x -> x[1], (0.0, 0.0, 0.0), -1.0, (0.0, 0.0, 1.0))
    @test_throws ArgumentError gyroaverage(x -> x[1], (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, Inf))

    constant_max = _ -> floatmax(Float64)
    @test gyroaverage(constant_max, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); n = 16) ==
          floatmax(Float64)
    smallest = nextfloat(0.0)
    @test gyroaverage(_ -> smallest, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); n = 16) == smallest
    ring_square = x -> x[1]^2 + x[2]^2
    @test gyroaverage(ring_square, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1e-200); n = 16) ≈ 1.0 rtol =
        8eps()
    @test gyroaverage(ring_square, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1e200); n = 16) ≈ 1.0 rtol =
        8eps()
    @test_throws ArgumentError gyroaverage(_ -> Inf, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0))
    @test_throws ArgumentError gyroaverage(_ -> 1 + im, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0))
    constant_one = _ -> 1.0
    gyroaverage(constant_one, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); n = 16)
    @test (@allocated gyroaverage(constant_one, (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0); n = 16)) == 0

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

    tiny_field = 1e-200
    gc3 = GuidingCentre((0.0, 0.0, 0.0), 0.0, 0.0, 1.0, 1.0)
    push_guiding_centre!(
        gc3;
        dt = 0.1,
        E = (tiny_field, 0.0, 0.0),
        B = (0.0, 0.0, tiny_field),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = 0.0,
    )
    @test collect(gc3.X) ≈ [0.0, -0.1, 0.0] rtol = 8eps()
    @test gc3.vpar == 0.0

    gc4 = GuidingCentre((0.0, 0.0, 0.0), 0.0, 1e200, 1e200, 1e200)
    push_guiding_centre!(
        gc4;
        dt = 1.0,
        E = (0.0, 0.0, 1e200),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = 1e200,
    )
    @test gc4.vpar == 0.0
    @test all(isfinite, gc4.X)

    t = nextfloat(0.0)
    M = floatmax(Float64)
    expected_force = Float64(setprecision(4096) do
        BigFloat(2t) * BigFloat(M) - BigFloat(M) * BigFloat(t)
    end)
    @test HybridPlasmaPIC._gyro_scaled_product_difference(2t, M, M, t) == expected_force
    gc5 = GuidingCentre((0.0, 0.0, 0.0), 0.0, M, 2t, M)
    push_guiding_centre!(
        gc5;
        dt = M,
        E = (0.0, 0.0, M),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = t,
    )
    @test gc5.vpar == expected_force
    @test all(isfinite, gc5.X)

    push_args = (;
        dt = 0.01,
        E = (0.0, 1.0, 0.0),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = 0.0,
    )
    push_guiding_centre!(gc3; push_args...)
    @test (@allocated push_guiding_centre!(gc3; push_args...)) == 0
end

@testset "fusion: adaptive mesh refinement" begin
    coarse = AMRGrid([1.0, 2.0, 3.0], 0.25; x0 = -0.5)
    fine = refine(coarse)
    @test ncells(fine) == 2ncells(coarse)
    @test effective_resolution(fine) ≈ 0.125
    @test cell_center(coarse, 1) ≈ -0.375
    @test !any(refine_flags(coarse.u, 1e-12))
    @test_throws ArgumentError refine_flags(coarse.u, -1.0)

    @test !any(refine_flags(fill(floatmax(Float64), 3), 1.0))
    @test refine_flags([1e308, -0.5, -1e308], 0.5) == Bool[false, true, false]
    tiny = nextfloat(0.0)
    crossed_scale = [floatmax(Float64), -tiny, -floatmax(Float64)]
    crossed_oracle = Float64(
        setprecision(4096) do
            BigFloat(crossed_scale[3]) - 2BigFloat(crossed_scale[2]) + BigFloat(crossed_scale[1])
        end,
    )
    @test HybridPlasmaPIC._amr_second_difference(crossed_scale, 2) == crossed_oracle
    @test refine_flags(crossed_scale, 0.0) == Bool[false, true, false]
    integer_extreme = [typemax(Int), typemin(Int), typemax(Int)]
    @test refine_flags(integer_extreme, 10) == Bool[false, true, false]

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

    @testset "extreme-scale transfers stay finite when the result is representable" begin
        for T in (Float32, Float64)
            M = floatmax(T)
            quarter_M = M / T(4)
            extreme_coarse = AMRGrid(T[0, -M, 0, M, 0], one(T))
            extreme_fine = refine(extreme_coarse)
            expected =
                T[quarter_M, -quarter_M, -M, -M, -quarter_M, quarter_M, M, M, quarter_M, -quarter_M]
            @test all(isfinite, extreme_fine.u)
            @test extreme_fine.u == expected

            extreme_back = AMRGrid(similar(extreme_coarse.u), extreme_coarse.dx)
            restrict!(extreme_back, extreme_fine)
            @test extreme_back.u == extreme_coarse.u

            one_coarse = AMRGrid(T[0], one(T))
            same_large = AMRGrid(fill(M, 2), T(0.5))
            restrict!(one_coarse, same_large)
            @test one_coarse.u == T[M]

            smallest = nextfloat(zero(T))
            same_smallest = AMRGrid(fill(smallest, 2), T(0.5))
            restrict!(one_coarse, same_smallest)
            @test one_coarse.u == T[smallest]

            @test (@allocated prolong!(extreme_fine, extreme_coarse)) == 0
            @test (@allocated restrict!(extreme_back, extreme_fine)) == 0
        end
    end

    @testset "geometry validation is finite and transactional" begin
        @test_throws ArgumentError AMRGrid(Float64[], 1.0)
        @test_throws ArgumentError AMRGrid([1.0], 0.0)
        @test_throws ArgumentError AMRGrid([1.0], -1.0)
        @test_throws ArgumentError AMRGrid([1.0], Inf)
        @test_throws ArgumentError AMRGrid([1.0], 1.0; x0 = NaN)
        @test_throws ArgumentError refine_flags([1.0, 2.0, 1.0], Inf)
        @test_throws ArgumentError refine_flags([1.0, 2.0, 1.0], NaN)

        integer_grid = AMRGrid([1, 2, 3], 0.5)
        @test eltype(integer_grid.u) === Float64
        @test integer_grid.dx == 0.5

        bad_dx = AMRGrid(fill(-7.0, 6), 0.2; x0 = coarse.x0)
        bad_dx_before = copy(bad_dx.u)
        @test_throws ArgumentError prolong!(bad_dx, coarse)
        @test bad_dx.u == bad_dx_before

        bad_origin = AMRGrid(fill(-9.0, 6), coarse.dx / 2; x0 = coarse.x0 + 1)
        bad_origin_before = copy(bad_origin.u)
        @test_throws ArgumentError prolong!(bad_origin, coarse)
        @test bad_origin.u == bad_origin_before

        coarse_out = AMRGrid(fill(-11.0, 3), coarse.dx; x0 = coarse.x0)
        coarse_out_before = copy(coarse_out.u)
        @test_throws ArgumentError restrict!(coarse_out, bad_origin)
        @test coarse_out.u == coarse_out_before
    end
end
