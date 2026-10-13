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
    # Verify exact rounding from neighboring representable values, independently
    # of both the production converter and Julia's ambient BigFloat conversion.
    function nearest_hardware_result(candidate::T, value::Rational{BigInt}) where {T}
        isnan(candidate) && return false
        iszero(value) && return iszero(candidate) && !signbit(candidate)
        signbit(candidate) == (value < 0) || return false
        exact(x) = Rational{BigInt}(x)
        x, magnitude = abs(candidate), abs(value)
        maximum = floatmax(T)
        overflow = exact(maximum) + (exact(maximum) - exact(prevfloat(maximum))) / 2
        isinf(x) && return magnitude >= overflow
        lower = iszero(x) ? zero(value) : (exact(prevfloat(x)) + exact(x)) / 2
        upper = x == maximum ? overflow : (exact(x) + exact(nextfloat(x))) / 2
        U = T === Float16 ? UInt16 : T === Float32 ? UInt32 : UInt64
        even = iseven(reinterpret(U, x))
        return (magnitude > lower || (magnitude == lower && even)) &&
               (magnitude < upper || (magnitude == upper && even))
    end

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
        ratio_exact =
            (Rational{BigInt}(ratio_numerators[1]) * Rational{BigInt}(ratio_numerators[2])) /
            (Rational{BigInt}(ratio_denominators[1]) * Rational{BigInt}(ratio_denominators[2]))
        ratio = HybridPlasmaPIC._gyro_scaled_ratio(ratio_numerators, ratio_denominators)
        @test ratio == -nextfloat(0.0)
        @test nearest_hardware_result(ratio, ratio_exact)

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

    @testset "public drifts retain raw cross components" begin
        function exact_cross_drift(a, b, numerators, denominators)
            exact(x) = Rational{BigInt}(x)
            coefficient =
                prod(exact, numerators; init = big(1) // big(1)) /
                prod(exact, denominators; init = big(1) // big(1))
            return ntuple(3) do i
                j, k = mod1(i + 1, 3), mod1(i + 2, 3)
                coefficient * (exact(a[j]) * exact(b[k]) - exact(a[k]) * exact(b[j]))
            end
        end
        rotate(a, shift) = ntuple(i -> a[mod1(i + shift, 3)], 3)
        for T in (Float16, Float32, Float64)
            z, o, M, tiny = zero(T), one(T), floatmax(T), nextfloat(zero(T))
            for sign in (-o, o), shift = 0:2
                field = rotate((M, sign * tiny, z), shift)
                axis = rotate((z, z, o), shift)
                for (actual, expected) in (
                    (exb_drift(field, axis), exact_cross_drift(field, axis, (), ())),
                    (
                        gradb_drift(o, sign, T(2), axis, field),
                        exact_cross_drift(axis, field, (T(2),), (T(2), sign)),
                    ),
                    (
                        curvature_drift(-o, sign, o, axis, field),
                        exact_cross_drift(axis, field, (o,), (sign,)),
                    ),
                    (
                        gradb_drift(M, sign * T(0.5), M, field, axis),
                        exact_cross_drift(field, axis, (M, M, M), (T(2), sign * T(0.5), M, M, M)),
                    ),
                    (
                        curvature_drift(-M, sign, o, field, axis),
                        exact_cross_drift(field, axis, (M, M), (sign, M, M)),
                    ),
                )
                    @test actual isa NTuple{3,T}
                    @test all(nearest_hardware_result(actual[i], expected[i]) for i = 1:3)
                end
                # The public composition must retain each repaired component.
                total = drift_velocity(;
                    vpar = z,
                    vperp = z,
                    q = sign,
                    m = o,
                    E = field,
                    B = axis,
                    gradB = field,
                    κ = field,
                )
                @test total == exb_drift(field, axis)
            end

            # Cross cancellation occurs before division by the magnetic norm.
            # Include raw products that individually overflow, with finite drifts.
            for (B, field) in
                (((o, T(2), z), (nextfloat(o), T(2), z)), ((M, prevfloat(M), z), (M, M, z)))
                Bs = maximum(abs, B)
                Bn = hypot((B ./ Bs)...)
                cases = (
                    (exb_drift(field, B), exact_cross_drift(field, B, (), (Bs, Bs, Bn, Bn))),
                    (
                        gradb_drift(o, -o, T(2), B, field),
                        exact_cross_drift(B, field, (T(2),), (T(2), -o, Bs, Bs, Bs, Bn, Bn, Bn)),
                    ),
                    (
                        curvature_drift(-o, -o, o, B, field),
                        exact_cross_drift(B, field, (o,), (-o, Bs, Bs, Bn, Bn)),
                    ),
                )
                for (actual, expected) in cases
                    @test all(nearest_hardware_result(actual[i], expected[i]) for i = 1:3)
                end
            end
            axis = (z, z, o)
            huge = (M, M, M)
            @test exb_drift(huge, huge) == (z, z, z)
            @test gradb_drift(z, o, M, (z, z, tiny), huge) == (z, z, z)
            @test curvature_drift(z, -o, M, (z, z, tiny), huge) == (z, z, z)
            @test_throws ArgumentError exb_drift((M, z, z), (z, z, T(0.5)))
            @test_throws ArgumentError gradb_drift(M, o, o, axis, (o, z, z))
            @test_throws ArgumentError curvature_drift(M, o, o, axis, (o, z, z))
            @test gradb_drift(o, 1, 1, axis, (o, z, z)) isa NTuple{3,T}
            @test curvature_drift(o, 1, 1, axis, (o, z, z)) isa NTuple{3,T}
        end
        @test exb_drift((1, 2, 0), (0, 0, 1)) isa NTuple{3,Float64}
        @test exb_drift(Float32.((1, 2, 0)), (0.0, 0.0, 1.0)) isa NTuple{3,Float64}

        # Retained BigFloat inputs can carry more bits than the working precision.
        # Cancellation must use those stored bits and preserve caller rounding.
        retained = setprecision(BigFloat, 256) do
            z, o = zero(BigFloat), one(BigFloat)
            ((o + ldexp(o, -200), o, z), (o, o, z))
        end
        for mode in (RoundNearest, RoundDown, RoundUp)
            setprecision(BigFloat, 64) do
                setrounding(BigFloat, mode) do
                    field, B = retained
                    Bs = maximum(abs, B)
                    Bn = hypot((B ./ Bs)...)
                    reference = setprecision(BigFloat, 2048) do
                        ldexp(one(BigFloat), -200) / (BigFloat(Bs) * BigFloat(Bn))^2
                    end
                    expected = BigFloat(reference)
                    @test exb_drift(field, B)[3] == expected
                    @test curvature_drift(
                        one(BigFloat),
                        -one(BigFloat),
                        one(BigFloat),
                        B,
                        field,
                    )[3] == expected
                    @test precision(BigFloat) == 64
                    @test rounding(BigFloat) == mode
                end
            end
        end
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

    @testset "ring geometry and streaming sample mean" begin
        for n in (8192, 12000, 70000)
            calls, finite_points = Ref(0), Ref(true)
            callback = pt -> begin
                calls[] += 1
                finite_points[] &= all(isfinite, pt)
                one(Float16)
            end
            @test gyroaverage(callback, Float16.((0, 0, 0)), Float16(0), Float16.((0, 0, 1)); n) ===
                  Float16(1)
            @test calls[] == n
            @test finite_points[]
        end

        permutations = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
        for T in (Float16, Float32, Float64)
            z, o, M, tiny = zero(T), one(T), floatmax(T), nextfloat(zero(T))
            center, axis = (z, z, z), (z, z, o)
            for constant in (M, tiny, -M, -tiny), n in (3, 16, 8192)
                @test gyroaverage(_ -> constant, center, z, axis; n) === constant
            end
            for sign in (-o, o)
                points = NTuple{3,T}[]
                @test gyroaverage(p -> (push!(points, p); o), center, T(2), (z, z, sign); n = 4) ===
                      o
                @test points ==
                      [(T(2), z, z), (z, sign * T(2), z), (-T(2), z, z), (z, -sign * T(2), z)]
            end
            # A tiny direction component remains representable after multiplying
            # by the radius. The first point is independently known in this plane.
            points = NTuple{3,T}[]
            gyroaverage(p -> (push!(points, p); z), center, M, (M, tiny, z); n = 4)
            @test first(points) == (-tiny, M, z)
            @test all(p -> all(isfinite, p), points)

            for values in (
                    (M, 3tiny, -M),
                    (T(3), T(3) * eps(T) / T(2), 3tiny),
                    (-T(3), -T(3) * eps(T) / T(2), -3tiny),
                ),
                order in permutations

                index = Ref(0)
                expected = sum(Rational{BigInt}, values) / 3
                actual = gyroaverage(_ -> values[order[(index[]+=1)]], center, z, axis; n = 3)
                @test nearest_hardware_result(actual, expected)
                @test index[] == 3
                @test actual isa T
            end

            # A nonrepresentable point is rejected before its callback; inputs
            # remain intact, and callback errors retain their original identity.
            calls = Ref(0)
            X = (M, z, z)
            @test_throws ArgumentError gyroaverage(_ -> (calls[] += 1; o), X, M, axis; n = 4)
            @test calls[] == 0
            @test X == (M, z, z)
            sentinel = ErrorException("ring callback failed")
            failure = try
                gyroaverage(_ -> begin
                    calls[] += 1
                    calls[] == 2 && throw(sentinel)
                    o
                end, center, z, axis; n = 4)
                nothing
            catch error
                error
            end
            @test failure === sentinel
            @test calls[] == 2
        end

        for bad in (Inf, NaN, 1 + im)
            calls = Ref(0)
            @test_throws ArgumentError gyroaverage(
                _ -> (calls[] += 1; calls[] == 2 ? bad : 1.0),
                (0.0, 0.0, 0.0),
                0.0,
                (0.0, 0.0, 1.0);
                n = 4,
            )
            @test calls[] == 2
        end
        @test_throws ArgumentError gyroaverage(
            _ -> 1.0,
            (0.0, 0.0, 0.0),
            0.0,
            (0.0, 0.0, 1.0);
            n = big(typemax(Int)) + 1,
        )
        values = (Float16(1), Float32(2), 3.0)
        index = Ref(0)
        mixed = gyroaverage(_ -> values[(index[]+=1)], Float16.((0, 0, 0)), 0, (0, 0, 1); n = 3)
        @test mixed === 2.0

        function stored_rational(x::BigFloat)
            setprecision(BigFloat, precision(x)) do
                Rational{BigInt}(x)
            end
        end
        function bigfloat_rounding_cell(value, exact, mode)
            iszero(exact) && return iszero(value)
            (isfinite(value) && !iszero(value)) || return false
            lower, current, upper = stored_rational.((prevfloat(value), value, nextfloat(value)))
            mode === RoundDown && return current <= exact < upper
            mode === RoundUp && return lower < exact <= current
            even = iseven(
                numerator(stored_rational(significand(abs(value))) * big(2)^(precision(value) - 1)),
            )
            left, right = (lower + current) / 2, (current + upper) / 2
            return (left < exact || (left == exact && even)) &&
                   (exact < right || (exact == right && even))
        end
        retained = setprecision(BigFloat, 256) do
            z, o = zero(BigFloat), one(BigFloat)
            δ = ldexp(o, -200)
            (; values = (o + δ, -o, z), X = (-o + δ, z, z), ρ = o + δ, δ)
        end
        for mode in (RoundNearest, RoundDown, RoundUp)
            setprecision(BigFloat, 64) do
                setrounding(BigFloat, mode) do
                    z, o = zero(BigFloat), one(BigFloat)
                    index = Ref(0)
                    average = gyroaverage(
                        _ -> retained.values[(index[]+=1)],
                        (z, z, z),
                        z,
                        (z, z, o);
                        n = 3,
                    )
                    @test bigfloat_rounding_cell(average, stored_rational(retained.δ) / 3, mode)
                    for (X, radius) in ((retained.X, o), ((-o, z, z), retained.ρ))
                        points = NTuple{3,BigFloat}[]
                        gyroaverage(p -> (push!(points, p); z), X, radius, (z, z, o); n = 4)
                        @test first(points)[1] == retained.δ
                    end
                    @test precision(BigFloat) == 64
                    @test rounding(BigFloat) == mode
                end
            end
        end
    end

    @testset "ordinary ring arithmetic retains complete coordinates" begin
        part_value(p) = (Rational{BigInt}(p[1]) + Rational{BigInt}(p[2])) * (big(2) // big(1))^p[3]
        function ordinary_ring_resources(::Type{T}) where {T}
            # Do not capture the enclosing testset's mutable field bindings.
            local callback, X, radius, B
            callback(p) = p[1]^2 + p[2]^2 + p[3]^2 + one(T)
            X, radius, B = T.((0.1, 0.2, 0.3)), T(0.7), T.((0.3, -1.2, 2))
            gyroaverage(callback, X, radius, B; n = 16)
            return @allocated gyroaverage(callback, X, radius, B; n = 16)
        end
        for T in (Float32, Float64)
            z, o = zero(T), one(T)
            @test ordinary_ring_resources(T) == 0
            for scale in (ldexp(o, -10), o, ldexp(o, 10)),
                B in (T.((0.3, -1.2, 2)), T.((2, 0.2, -0.1)))

                X, radius = scale .* T.((0.1, 0.2, 0.3)), scale * T(0.7)
                geometry = HybridPlasmaPIC._gyro_ring_geometry(X, radius, B)
                @test geometry.aligned.usable
                for i = 1:3, j = 1:3
                    original = (geometry.center[i], geometry.cosine[i], geometry.sine[i])[j]
                    cache = geometry.aligned.coordinates[i]
                    @test part_value(cache.terms[j]) * (big(2) // big(1))^cache.power ==
                          part_value(original)
                end
                points = NTuple{3,T}[]
                average = gyroaverage(X, radius, B; n = 16) do point
                    push!(points, point)
                    point[1]^2 + point[2]^2 + point[3]^2 + o
                end
                @test length(points) == 16
                samples = map(p -> p[1]^2 + p[2]^2 + p[3]^2 + o, points)
                @test nearest_hardware_result(average, sum(Rational{BigInt}, samples) / 16)
                reference = geometry.reference
                j, k = mod1(reference + 1, 3), mod1(reference + 2, 3)
                H = Rational{BigInt}(geometry.Bscale) * Rational{BigInt}(geometry.Bnorm)
                P = Rational{BigInt}(geometry.Pscale) * Rational{BigInt}(geometry.Pnorm)
                for phase = 0:15
                    s, c = T.(sincospi(2 * (phase / 16)))
                    for i = 1:3
                        u =
                            i == reference ? P / H :
                            -Rational{BigInt}(B[reference]) * Rational{BigInt}(B[i]) / (H * P)
                        v =
                            i == reference ? big(0) // big(1) :
                            i == j ? Rational{BigInt}(B[k]) / P : -Rational{BigInt}(B[j]) / P
                        exact =
                            Rational{BigInt}(X[i]) +
                            Rational{BigInt}(radius) *
                            (Rational{BigInt}(c) * u + Rational{BigInt}(s) * v)
                        @test abs(Rational{BigInt}(points[phase+1][i]) - exact) <=
                              Rational{BigInt}(eps(T)) * abs(exact)
                    end
                end
            end
            # Binary alignment must retain every low bit. A crossed-scale basis
            # cannot use the cached arithmetic even if its head would round away.
            radius = nextfloat(z)
            geometry = HybridPlasmaPIC._gyro_ring_geometry((o, z, z), radius, (z, z, o))
            @test !geometry.aligned.usable
            @test HybridPlasmaPIC._gyro_ring_point((o, z, z), radius, (z, z, o), geometry, z, o) ==
                  (o, radius, z)
            geometry = HybridPlasmaPIC._gyro_ring_geometry((z, z, z), T(0.5), (z, z, o))
            floor = geometry.aligned.floor
            for phase in (prevfloat(floor), floor, nextfloat(floor)), sign in (-o, o)
                c = sign * phase
                radius = T(0.5)
                for center in (-radius * c, prevfloat(-radius * c), nextfloat(-radius * c))
                    X, B = (center, z, z), (z, z, o)
                    geometry = HybridPlasmaPIC._gyro_ring_geometry(X, radius, B)
                    point = HybridPlasmaPIC._gyro_ring_point(X, radius, B, geometry, c, o)
                    exact =
                        Rational{BigInt}(center) + Rational{BigInt}(radius) * Rational{BigInt}(c)
                    @test nearest_hardware_result(point[1], exact)
                    @test point[2:3] == (radius, z)
                end
            end
            # Straddle the cancellation guard and final-output range boundaries.
            for radius in (floatmin(T), T(2) * floatmin(T), o, floatmax(T) / T(2)),
                center in (-radius, prevfloat(-radius), nextfloat(-radius))

                X, B = (center, z, z), (z, z, o)
                geometry = HybridPlasmaPIC._gyro_ring_geometry(X, radius, B)
                point = HybridPlasmaPIC._gyro_ring_point(X, radius, B, geometry, o, z)
                @test nearest_hardware_result(
                    point[1],
                    Rational{BigInt}(center) + Rational{BigInt}(radius),
                )
            end
            for multiplier in (T(118), T(120), T(122))
                X, B = (-o + multiplier * eps(T), z, z), (z, z, o)
                geometry = HybridPlasmaPIC._gyro_ring_geometry(X, o, B)
                point = HybridPlasmaPIC._gyro_ring_point(X, o, B, geometry, o, z)
                @test nearest_hardware_result(point[1], Rational{BigInt}(X[1]) + 1)
            end
        end
    end

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
            expected_v = exact(vpar) + expected_delta
            expected_z = exact(dt) * exact(vpar)
            push_guiding_centre!(
                gc_parallel;
                dt,
                E = (zero(T), zero(T), Epar),
                B = (zero(T), zero(T), one(T)),
                gradB = z,
                κ = z,
                gradpar_B,
            )
            @test nearest_hardware_result(gc_parallel.vpar, expected_v)
            @test gc_parallel.X[1:2] == (zero(T), zero(T))
            @test nearest_hardware_result(gc_parallel.X[3], expected_z)
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

            # A retained low product can prove cancellation; an underflowed
            # low*low contribution cannot be silently declared exact.
            parts = HybridPlasmaPIC._gyro_ratio_parts((T(0.7), T(0.2)), ())
            represented =
                (Rational{BigInt}(parts[1]) + Rational{BigInt}(parts[2])) *
                (BigInt(2) // BigInt(1))^parts[3]
            @test parts[4]
            @test represented == Rational{BigInt}(T(0.7)) * Rational{BigInt}(T(0.2))
            opposite = (-parts[1], -parts[2], parts[3], parts[4])
            cancelled, uncertain = HybridPlasmaPIC._gyro_sum_parts((parts, opposite))
            @test cancelled[4] && iszero(cancelled[1]) && iszero(cancelled[2])
            @test !uncertain
            low_product = HybridPlasmaPIC._gyro_multiply_parts(
                (T(0.5), tiny, 0, true),
                (T(0.5), tiny, 0, true),
            )
            @test !low_product[4]
            _, uncertain = HybridPlasmaPIC._gyro_sum_parts((
                (T(0.5), tiny, -1, true),
                (-T(0.5), zero(T), -1, true),
                (T(0.5), zero(T), 1, true),
                (-T(0.5), zero(T), 1, true),
            ))
            @test uncertain
        end

        # Float16 can round a nonzero FMA residual to zero. Such a residual
        # cannot certify a product even when the rounded product is nonzero.
        a, b = Float16(0.006355), Float16(0.006924)
        residual = fma(a, b, -(a * b))
        @test iszero(residual)
        @test Rational{BigInt}(a) * Rational{BigInt}(b) != Rational{BigInt}(a * b)
        @test !HybridPlasmaPIC._gyro_product_is_exact(a, b, residual)
    end

    @testset "position retains the complete displacement expression" begin
        function position_oracle(gc; dt, E, B, gradB, κ, gradpar_B)
            return setprecision(8192) do
                eb, bb, gb, kb = map(x -> BigFloat.(x), (E, B, gradB, κ))
                magnitude = sqrt(sum(abs2, bb))
                cross(a, b) = (
                    a[2] * b[3] - a[3] * b[2],
                    a[3] * b[1] - a[1] * b[3],
                    a[1] * b[2] - a[2] * b[1],
                )
                v, μ, q, m = BigFloat.((gc.vpar, gc.μ, gc.q, gc.m))
                velocity =
                    v .* bb ./ magnitude .+
                    (cross(eb, bb) .+ μ / q .* cross(bb, gb) .+ m * v^2 / q .* cross(bb, kb)) ./
                    magnitude^2
                BigFloat.(gc.X) .+ BigFloat(dt) .* velocity
            end
        end

        function check_position_case(
            T;
            X = (zero(T), zero(T), zero(T)),
            vpar = zero(T),
            μ = zero(T),
            q = one(T),
            m = one(T),
            dt = one(T),
            E = (zero(T), zero(T), zero(T)),
            B = (zero(T), zero(T), one(T)),
            gradB = (zero(T), zero(T), zero(T)),
            κ = (zero(T), zero(T), zero(T)),
        )
            gc_position = GuidingCentre(X, vpar, μ, q, m)
            args = (; dt, E, B, gradB, κ, gradpar_B = zero(T))
            expected = map(x -> T(x, RoundNearest), position_oracle(gc_position; args...))
            @test all(isfinite, expected)
            push_guiding_centre!(gc_position; args...)
            @test gc_position.X == expected
            @test gc_position.vpar == vpar
            @test gc_position.μ == μ
        end

        function check_parallel_field_allocations(T)
            z = (zero(T), zero(T), zero(T))
            field = T.((0.3, -1.2, 2))
            for vpar in (zero(T), T(0.25))
                gc_position = GuidingCentre(z, vpar, one(T), one(T), one(T))
                args = (;
                    dt = T(0.01),
                    E = field,
                    B = field,
                    gradB = field,
                    κ = field,
                    gradpar_B = zero(T),
                )
                expected = position_oracle(gc_position; args...)
                push_guiding_centre!(gc_position; args...)
                @test maximum(abs.(BigFloat.(gc_position.X) .- expected)) <= 8eps(T)
                gc_position.X, gc_position.vpar = z, vpar
                @test (@allocated push_guiding_centre!(gc_position; args...)) == 0
            end
            # Electric and μ-gradient drifts balance before multiplication by
            # a generally inexact timestep, for either charge sign.
            field = T.((0.5, 0.2, 0))
            for charge in (T(0.7), -T(0.7))
                gc_position = GuidingCentre(z, zero(T), abs(charge), charge, T(3))
                args = (;
                    dt = T(0.1),
                    E = field,
                    B = (zero(T), zero(T), T(2)),
                    gradB = sign(charge) .* field,
                    κ = z,
                    gradpar_B = zero(T),
                )
                @test all(iszero, position_oracle(gc_position; args...))
                push_guiding_centre!(gc_position; args...)
                @test gc_position.X == z
                @test (@allocated push_guiding_centre!(gc_position; args...)) == 0
            end
            for charge in (T(3), -T(3))
                gc_position = GuidingCentre(z, zero(T), T(2), charge, T(3))
                args = (;
                    dt = T(0.1),
                    E = T.((0.5, 2, 0)),
                    B = (zero(T), zero(T), one(T)),
                    gradB = sign(charge) .* T.((0.75, 3, 0)),
                    κ = z,
                    gradpar_B = zero(T),
                )
                @test all(iszero, position_oracle(gc_position; args...))
                push_guiding_centre!(gc_position; args...)
                @test gc_position.X == z
                @test (@allocated push_guiding_centre!(gc_position; args...)) == 0
            end
        end

        for T in (Float16, Float32, Float64)
            M, tiny, z, o = floatmax(T), nextfloat(zero(T)), zero(T), one(T)
            check_position_case(T; vpar = M, B = (M, tiny, z))
            check_position_case(T; X = (z, z, -M), vpar = M, dt = T(2))
            check_position_case(T; E = (M, tiny, z), X = (z, M, z))
            check_position_case(T; E = (M, z, z), B = (z, z, tiny), dt = tiny)
            check_position_case(T; E = (z, M, z), μ = M, q = T(0.5), gradB = (z, o, z))
            check_position_case(T; E = (z, tiny, z), B = (z, z, T(2)), dt = T(2))
            check_position_case(T; E = (z, M, z), μ = M, gradB = (z, -o, z), dt = T(0.5))
            check_position_case(T; μ = o, q = M, gradB = (M, tiny, z), dt = M)
            check_position_case(T; X = (z, z, -M), vpar = o, q = M, κ = (M, tiny, z), dt = M)

            # Zero time does not require the unadvanced drift to fit in T.
            check_position_case(
                T;
                X = (o, -o, z),
                vpar = M,
                μ = M,
                E = (M, z, z),
                B = (z, z, tiny),
                gradB = (M, z, z),
                κ = (M, z, z),
                dt = z,
            )

            gc_overflow = GuidingCentre((z, z, z), M, z, o, o)
            before = (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m)
            @test_throws ArgumentError push_guiding_centre!(
                gc_overflow;
                dt = T(2),
                E = (z, z, z),
                B = (z, z, o),
                gradB = (z, z, z),
                κ = (z, z, z),
                gradpar_B = z,
            )
            @test (gc_overflow.X, gc_overflow.vpar, gc_overflow.μ, gc_overflow.q, gc_overflow.m) ==
                  before
            check_parallel_field_allocations(T)
        end

        @testset "exact fallbacks preserve retained BigFloat state" begin
            ambient_precision = precision(BigFloat)
            ambient_rounding = rounding(BigFloat)
            δ_exact = big(1) // big(2)^200
            for working_precision in (64, 128),
                mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)

                retained = setprecision(256) do
                    z, o = BigFloat(0), BigFloat(1)
                    δ = BigFloat(δ_exact)
                    high = o + δ
                    Z, axis = (z, z, z), (z, z, o)
                    position_cases = (
                        (GuidingCentre((z, z, -o), high, z, o, o), Z, Z, Z),
                        (GuidingCentre((z, z, -o + δ), o, z, o, o), Z, Z, Z),
                        (GuidingCentre(Z, z, high, o, o), (z, o, z), (z, o, z), Z),
                        (GuidingCentre(Z, z, o, high, o), (z, o, z), (z, o, z), Z),
                        (GuidingCentre(Z, o, z, o, high), (z, o, z), Z, (z, o, z)),
                    )
                    force_cases = (
                        (GuidingCentre(Z, z, o, high, o), o, δ_exact),
                        (GuidingCentre(Z, z, high, o, o), o, -δ_exact),
                        (GuidingCentre(Z, -o + δ, z, o, o), z, δ_exact),
                        (GuidingCentre(Z, -o, z, o, high), z, -δ_exact / (1 + δ_exact)),
                    )
                    (; z, o, Z, axis, high, position_cases, force_cases)
                end
                setprecision(working_precision) do
                    setrounding(BigFloat, mode) do
                        @test HybridPlasmaPIC._exact_rational(retained.high) == 1 + δ_exact
                        for (gc_retained, E, gradB, κ) in retained.position_cases
                            @test all(
                                x -> precision(x) == 256,
                                (
                                    gc_retained.X...,
                                    gc_retained.vpar,
                                    gc_retained.μ,
                                    gc_retained.q,
                                    gc_retained.m,
                                ),
                            )
                            args = (;
                                dt = retained.o,
                                E,
                                B = retained.axis,
                                gradB,
                                κ,
                                gradpar_B = retained.z,
                            )
                            expected = BigFloat.(position_oracle(gc_retained; args...))
                            μ_before = gc_retained.μ
                            push_guiding_centre!(gc_retained; args...)
                            @test gc_retained.X == expected
                            @test all(x -> precision(x) == working_precision, gc_retained.X)
                            @test gc_retained.μ == μ_before
                            @test precision(BigFloat) == working_precision
                            @test rounding(BigFloat) == mode
                        end
                        for (gc_retained, gradpar_B, exact_velocity) in retained.force_cases
                            expected = BigFloat(setprecision(8192) do
                                BigFloat(exact_velocity)
                            end)
                            push_guiding_centre!(
                                gc_retained;
                                dt = retained.o,
                                E = retained.axis,
                                B = retained.axis,
                                gradB = retained.Z,
                                κ = retained.Z,
                                gradpar_B,
                            )
                            @test gc_retained.vpar == expected
                            @test precision(gc_retained.vpar) == working_precision
                            @test precision(BigFloat) == working_precision
                            @test rounding(BigFloat) == mode
                        end
                    end
                end
                @test precision(BigFloat) == ambient_precision
                @test rounding(BigFloat) == ambient_rounding
            end
        end
    end

    @testset "hardware fallback rounding ignores BigFloat settings" begin
        ambient_precision, ambient_rounding = precision(BigFloat), rounding(BigFloat)
        for T in (Float16, Float32, Float64),
            working_precision in (2, precision(T), 256),
            mode in (RoundNearest, RoundDown, RoundUp, RoundToZero)

            z, o = zero(T), one(T)
            tiny, half_ulp = nextfloat(z), eps(o) / 2
            axis, zeros3 = (z, z, o), (z, z, z)
            setprecision(BigFloat, working_precision) do
                setrounding(BigFloat, mode) do
                    for sign in (-1, 1), side in (-1, 0, 1)
                        exact_result =
                            sign * (
                                Rational{BigInt}(o) +
                                Rational{BigInt}(half_ulp) +
                                side * Rational{BigInt}(tiny)
                            )
                        expected = sign * (side > 0 ? nextfloat(o) : o)
                        @test nearest_hardware_result(expected, exact_result)
                        gc_position = GuidingCentre((sign * o, z, z), z, o, o, o)
                        push_guiding_centre!(
                            gc_position;
                            dt = o,
                            E = (z, sign * half_ulp, z),
                            B = axis,
                            gradB = (z, -sign * side * tiny, z),
                            κ = zeros3,
                            gradpar_B = z,
                        )
                        @test gc_position.X == (expected, z, z)
                        @test gc_position.vpar == z
                        gc_force = GuidingCentre(zeros3, sign * o, o, o, o)
                        push_guiding_centre!(
                            gc_force;
                            dt = o,
                            E = (z, z, sign * half_ulp),
                            B = axis,
                            gradB = zeros3,
                            κ = zeros3,
                            gradpar_B = -sign * side * tiny,
                        )
                        @test gc_force.X == (z, z, sign * o)
                        @test gc_force.vpar == expected
                    end
                    # Exercise the older exact product-ratio fallback through
                    # its public drift caller as well as the complete push.
                    drift = exb_drift((z, T(5) * tiny, z), (z, z, T(9)))
                    @test drift == (tiny, z, z)
                    @test nearest_hardware_result(drift[1], 5Rational{BigInt}(tiny) / 9)
                    @test precision(BigFloat) == working_precision
                    @test rounding(BigFloat) == mode
                end
            end
        end
        @test precision(BigFloat) == ambient_precision
        @test rounding(BigFloat) == ambient_rounding
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
