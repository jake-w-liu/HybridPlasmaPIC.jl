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
    gc_force = GuidingCentre((0.0, 0.0, 0.0), 0.0, M, 2t, 1.0)
    push_guiding_centre!(
        gc_force;
        dt = 1.0,
        E = (0.0, 0.0, M),
        B = (0.0, 0.0, 1.0),
        gradB = (0.0, 0.0, 0.0),
        κ = (0.0, 0.0, 0.0),
        gradpar_B = t,
    )
    @test gc_force.vpar == expected_force
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

    @testset "magnetic-moment drift avoids perpendicular-speed squaring" begin
        # The stored moment gives v_gradB = μ (b × ∇B)/(q |B|) directly.
        # Evaluate the documented local-field equations independently at high
        # precision, including cases where v_perp^2 is outside the state range.
        function moment_push_oracle(gc; dt, E, B, gradB, κ, gradpar_B)
            return setprecision(2048) do
                eb, bb, gb, kb = map(x -> BigFloat.(x), (E, B, gradB, κ))
                Bmag = sqrt(sum(abs2, bb))
                b = bb ./ Bmag
                cross(a, c) = (
                    a[2] * c[3] - a[3] * c[2],
                    a[3] * c[1] - a[1] * c[3],
                    a[1] * c[2] - a[2] * c[1],
                )
                μ, q, m, vpar = BigFloat.((gc.μ, gc.q, gc.m, gc.vpar))
                velocity =
                    vpar .* b .+ cross(eb, bb) ./ Bmag^2 .+ (μ / q / Bmag) .* cross(b, gb) .+
                    (m * vpar^2 / q / Bmag) .* cross(b, kb)
                X = BigFloat.(gc.X) .+ BigFloat(dt) .* velocity
                v = vpar + BigFloat(dt) * (q * sum(eb .* b) - μ * BigFloat(gradpar_B)) / m
                return X, v
            end
        end

        maximum_error = 0.0
        for T in (Float32, Float64)
            scales = T === Float32 ? (T(1e-30), T(1e30)) : (T(1e-200), T(1e200))
            z = (zero(T), zero(T), zero(T))
            for scale in scales, charge in (-T(2), T(2)), mass in (T(0.5), T(4.25))
                gc_moment = GuidingCentre(T.((0.2, -0.1, 0.3)), T(-0.35), scale, charge, mass)
                args = (;
                    dt = T(0.125),
                    E = z,
                    B = scale .* T.((0.2, -0.6, 0.7)),
                    gradB = T.((0.3, -0.4, 0.5)),
                    κ = z,
                    gradpar_B = zero(T),
                )
                expected_X, expected_v = moment_push_oracle(gc_moment; args...)
                push_guiding_centre!(gc_moment; args...)
                error = maximum(abs.(BigFloat.(gc_moment.X) .- expected_X))
                maximum_error = max(maximum_error, Float64(error))
                @test error < 16eps(T)
                @test gc_moment.vpar == T(expected_v)
                @test gc_moment.μ == scale
                @test gc_moment isa GuidingCentre{T}
            end

            gc_normal = GuidingCentre(T.((0.2, -0.1, 0.3)), T(-0.35), T(0.8), -T(2), T(4.25))
            normal_args = (;
                dt = T(0.125),
                E = T.((0.4, -0.2, 0.3)),
                B = T.((0.2, -0.6, 0.7)),
                gradB = T.((0.3, -0.4, 0.5)),
                κ = T.((0.1, 0.2, -0.3)),
                gradpar_B = T(0.2),
            )
            expected_X, expected_v = moment_push_oracle(gc_normal; normal_args...)
            push_guiding_centre!(gc_normal; normal_args...)
            @test maximum(abs.(BigFloat.(gc_normal.X) .- expected_X)) < 16eps(T)
            @test abs(BigFloat(gc_normal.vpar) - expected_v) < 16eps(T)

            gc_zero = GuidingCentre(z, zero(T), zero(T), one(T), one(T))
            push_guiding_centre!(gc_zero; normal_args..., E = z, κ = z)
            @test gc_zero.X == z
            @test gc_zero.vpar == zero(T)

            gc_overflow =
                GuidingCentre((zero(T), floatmax(T), zero(T)), zero(T), one(T), one(T), one(T))
            before = (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m)
            @test_throws ArgumentError push_guiding_centre!(
                gc_overflow;
                dt = floatmax(T),
                E = z,
                B = (zero(T), zero(T), one(T)),
                gradB = (one(T), zero(T), zero(T)),
                κ = z,
                gradpar_B = zero(T),
            )
            @test (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m) ==
                  before
        end
        @info "guiding-centre magnetic-moment oracle" maximum_error
    end

    @testset "parallel acceleration retains the scaled force difference" begin
        # Exact rational arithmetic supplies the local constant-field Euler
        # oracle without rounding q*E, μ*gradpar_B, or dt/m independently.
        function check_parallel_case(T; q, Epar, μ, gradpar_B, m, dt, vpar = T(0.125))
            z = (zero(T), zero(T), zero(T))
            gc_parallel = GuidingCentre(z, vpar, μ, q, m)
            exact(x) = Rational{BigInt}(x)
            expected_delta =
                exact(dt) * (exact(q) * exact(Epar) - exact(μ) * exact(gradpar_B)) / exact(m)
            expected_v = T(exact(vpar) + expected_delta)
            expected_z = T(exact(dt) * exact(vpar))
            push_guiding_centre!(
                gc_parallel;
                dt,
                E = (zero(T), zero(T), Epar),
                B = (zero(T), zero(T), one(T)),
                gradB = z,
                κ = z,
                gradpar_B,
            )
            @test gc_parallel.vpar == expected_v
            @test gc_parallel.X == (zero(T), zero(T), expected_z)
            @test gc_parallel.μ == μ
        end

        function check_parallel_allocations(T)
            z = (zero(T), zero(T), zero(T))
            gc_parallel = GuidingCentre(z, T(0.125), one(T), one(T), one(T))
            for gradpar_B in (zero(T), T(0.2), T(0.5))
                args = (;
                    dt = T(0.125),
                    E = (zero(T), zero(T), T(0.5)),
                    B = (zero(T), zero(T), one(T)),
                    gradB = z,
                    κ = z,
                    gradpar_B,
                )
                push_guiding_centre!(gc_parallel; args...)
                @test (@allocated push_guiding_centre!(gc_parallel; args...)) == 0
            end
        end

        for T in (Float32, Float64)
            # Squaring these binary scales exceeds either floating-point range;
            # the combined scaled acceleration is exactly representable.
            exponent = T === Float32 ? 80 : 600
            for scale in (ldexp(one(T), -exponent), ldexp(one(T), exponent))
                for sign in (-one(T), one(T))
                    check_parallel_case(
                        T;
                        q = sign * scale,
                        Epar = scale,
                        μ = zero(T),
                        gradpar_B = zero(T),
                        m = scale,
                        dt = inv(scale),
                    )
                    check_parallel_case(
                        T;
                        q = one(T),
                        Epar = zero(T),
                        μ = scale,
                        gradpar_B = sign * scale,
                        m = scale,
                        dt = inv(scale),
                    )
                end
                check_parallel_case(
                    T;
                    q = scale,
                    Epar = scale,
                    μ = scale,
                    gradpar_B = -scale,
                    m = scale,
                    dt = inv(scale),
                )
                check_parallel_case(
                    T;
                    q = scale,
                    Epar = scale,
                    μ = scale,
                    gradpar_B = scale,
                    m = scale,
                    dt = inv(scale),
                )
                check_parallel_case(
                    T;
                    q = scale,
                    Epar = nextfloat(scale),
                    μ = scale,
                    gradpar_B = scale,
                    m = scale,
                    dt = inv(scale),
                    vpar = eps(T) / T(2),
                )
                check_parallel_case(
                    T;
                    q = scale,
                    Epar = scale,
                    μ = zero(T),
                    gradpar_B = zero(T),
                    m = scale,
                    dt = -inv(scale),
                )
                check_parallel_case(
                    T;
                    q = scale,
                    Epar = scale,
                    μ = zero(T),
                    gradpar_B = zero(T),
                    m = scale,
                    dt = zero(T),
                )
            end

            large = ldexp(one(T), exponent)
            small = inv(large)
            check_parallel_case(
                T;
                q = small,
                Epar = one(T),
                μ = zero(T),
                gradpar_B = zero(T),
                m = small,
                dt = large,
            )
            check_parallel_case(
                T;
                q = large,
                Epar = one(T),
                μ = zero(T),
                gradpar_B = zero(T),
                m = large,
                dt = small,
                vpar = small,
            )
            tiny = nextfloat(zero(T))
            maximum = floatmax(T)
            check_parallel_case(
                T;
                q = 2tiny,
                Epar = maximum,
                μ = maximum,
                gradpar_B = tiny,
                m = maximum,
                dt = maximum,
                vpar = zero(T),
            )
            for dt in (T(0.5), one(T), T(1.5)), sign in (-one(T), one(T))
                check_parallel_case(
                    T;
                    q = sign * tiny,
                    Epar = one(T),
                    μ = zero(T),
                    gradpar_B = zero(T),
                    m = one(T),
                    dt,
                    vpar = zero(T),
                )
            end
            check_parallel_case(
                T;
                q = maximum,
                Epar = one(T),
                μ = zero(T),
                gradpar_B = zero(T),
                m = one(T),
                dt = one(T),
                vpar = zero(T),
            )
            for dt in (T(3), -T(3))
                check_parallel_case(
                    T;
                    q = tiny,
                    Epar = T(1.5),
                    μ = tiny,
                    gradpar_B = one(T),
                    m = one(T),
                    dt,
                    vpar = zero(T),
                )
            end
            check_parallel_case(
                T;
                q = maximum,
                Epar = one(T),
                μ = maximum / T(2),
                gradpar_B = one(T),
                m = one(T),
                dt = T(2),
                vpar = zero(T),
            )
            check_parallel_allocations(T)

            z = (zero(T), zero(T), zero(T))
            for (vpar, Epar, m, dt) in
                ((zero(T), one(T), small, maximum), (maximum / T(2), maximum, one(T), one(T)))
                gc_overflow = GuidingCentre(z, vpar, zero(T), one(T), m)
                before =
                    (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m)
                @test_throws ArgumentError push_guiding_centre!(
                    gc_overflow;
                    dt,
                    E = (zero(T), zero(T), Epar),
                    B = (zero(T), zero(T), one(T)),
                    gradB = z,
                    κ = z,
                    gradpar_B = zero(T),
                )
                @test (
                    gc_overflow.X,
                    gc_overflow.vpar,
                    gc_overflow.μ,
                    gc_overflow.q,
                    gc_overflow.m,
                ) == before
            end
        end
    end

    @testset "parallel velocity evaluates the full field expression" begin
        function check_full_parallel_case(
            T;
            E,
            B,
            q,
            m,
            dt,
            vpar = zero(T),
            μ = zero(T),
            gradpar_B = zero(T),
        )
            z = (zero(T), zero(T), zero(T))
            gc_parallel = GuidingCentre(z, vpar, μ, q, m)
            expected = setprecision(4096) do
                eb, bb = BigFloat.(E), BigFloat.(B)
                Epar = sum(eb .* bb) / sqrt(sum(abs2, bb))
                BigFloat(vpar) +
                BigFloat(dt) / BigFloat(m) *
                (BigFloat(q) * Epar - BigFloat(μ) * BigFloat(gradpar_B))
            end
            push_guiding_centre!(gc_parallel; dt, E, B, gradB = z, κ = z, gradpar_B)
            @test abs(BigFloat(gc_parallel.vpar) - expected) <= 16eps(T) * abs(expected)
            @test all(isfinite, gc_parallel.X)
            @test gc_parallel.μ == μ
        end

        for T in (Float32, Float64)
            maximum, tiny = floatmax(T), nextfloat(zero(T))
            z = (zero(T), zero(T), zero(T))
            axial = (zero(T), zero(T), one(T))
            small = ldexp(one(T), T === Float32 ? -80 : -600)
            check_full_parallel_case(
                T;
                E = (maximum, maximum, maximum),
                B = (one(T), one(T), one(T)),
                q = small,
                m = one(T),
                dt = one(T),
            )
            check_full_parallel_case(
                T;
                E = (maximum, tiny, zero(T)),
                B = (zero(T), one(T), zero(T)),
                q = maximum,
                m = tiny,
                dt = one(T),
            )
            check_full_parallel_case(
                T;
                E = (zero(T), maximum, zero(T)),
                B = (maximum, tiny, zero(T)),
                q = maximum,
                m = tiny,
                dt = one(T),
            )
            check_full_parallel_case(
                T;
                E = (maximum, -maximum, tiny),
                B = (one(T), one(T), one(T)),
                q = maximum,
                m = tiny,
                dt = one(T),
            )
            check_full_parallel_case(
                T;
                E = (zero(T), zero(T), T(4)),
                B = axial,
                q = maximum,
                m = one(T),
                dt = T(0.5),
                vpar = -maximum,
            )
            check_full_parallel_case(
                T;
                E = (zero(T), zero(T), T(4)),
                B = axial,
                q = maximum,
                m = one(T),
                dt = T(0.5),
                vpar = -maximum,
                μ = maximum,
                gradpar_B = one(T),
            )
            check_full_parallel_case(
                T;
                E = (zero(T), zero(T), nextfloat(one(T))),
                B = axial,
                q = one(T),
                m = one(T),
                dt = T(0.5),
                vpar = -one(T),
                μ = one(T),
                gradpar_B = -one(T),
            )
            check_full_parallel_case(
                T;
                E = (maximum, maximum, maximum),
                B = (one(T), one(T), one(T)),
                q = maximum,
                m = tiny,
                dt = zero(T),
                vpar = T(0.125),
            )

            gc_overflow = GuidingCentre(z, zero(T), zero(T), maximum, one(T))
            before = (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m)
            @test_throws ArgumentError push_guiding_centre!(
                gc_overflow;
                dt = one(T),
                E = (zero(T), zero(T), T(4)),
                B = axial,
                gradB = z,
                κ = z,
                gradpar_B = zero(T),
            )
            @test (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m) ==
                  before
        end
    end

    @testset "certified stationary force balance stays allocation-free" begin
        function check_stationary_balance(T, q, μ, Epar, gradpar_B, mass, dt, Bz)
            q, μ, Epar, gradpar_B, mass, dt, Bz = T.((q, μ, Epar, gradpar_B, mass, dt, Bz))
            exact(x) = Rational{BigInt}(x)
            @test exact(q) * exact(Epar) * sign(Bz) == exact(μ) * exact(gradpar_B)
            z = (zero(T), zero(T), zero(T))
            gc_balance = GuidingCentre(z, zero(T), μ, q, mass)
            args = (;
                dt,
                E = (zero(T), zero(T), Epar),
                B = (zero(T), zero(T), Bz),
                gradB = z,
                κ = z,
                gradpar_B,
            )
            push_guiding_centre!(gc_balance; args...)
            @test gc_balance.vpar == zero(T)
            @test gc_balance.X == z
            @test (@allocated push_guiding_centre!(gc_balance; args...)) == 0
        end
        for T in (Float16, Float32, Float64)
            check_stationary_balance(T, 1, 1, 1, 1, 1, 0.1, 1)
            check_stationary_balance(T, -2, 4, -0.5, 0.25, 0.5, -0.1, 2)
            check_stationary_balance(T, 0.5, 0.25, -0.5, 1, 2, 0.2, -0.5)

            # A tiny aligned term must not disappear behind an exact large pair.
            tiny = nextfloat(zero(T))
            terms = (
                (T(0.5), zero(T), 1, true),
                (T(0.5), zero(T), exponent(tiny), true),
                (-T(0.5), zero(T), 1, true),
            )
            _, needs_exact = HybridPlasmaPIC._gyro_sum_scaled_parts(terms)
            @test needs_exact
        end

        # Float16 can round a nonzero FMA residual to zero. Such a residual
        # cannot certify a product even when the rounded product is nonzero.
        a, b = Float16(0.006355), Float16(0.006924)
        residual = fma(a, b, -(a * b))
        @test iszero(residual)
        @test Rational{BigInt}(a) * Rational{BigInt}(b) != Rational{BigInt}(a * b)
        @test !HybridPlasmaPIC._gyro_product_is_exact(a, b, residual)
    end

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
