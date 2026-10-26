# Dimension-parametric electromagnetic PIC. These tests cover the 2D/3D paths
# that are intentionally separate from the legacy, more specialized EMPIC1D
# regression suite.

using HybridPlasmaPIC, FFTW, Test, Random

module ContinuityRangeOracle
using HybridPlasmaPIC
const Q = Rational{BigInt}
const CQ = Complex{Q}
exact(x) = Q(x)
complex_exact(x) = complex(exact(real(x)), exact(imag(x)))

# Gaussian-integer roots give an exact, independent DFT on these small grids.
function dft4(values; inverse = false)
    roots = CQ[1, im, -1, -im]
    direction = inverse ? 1 : -1
    return map(CartesianIndices(values)) do I
        total = zero(CQ)
        for K in CartesianIndices(values)
            phase =
                direction *
                sum((K[d] - 1) * (I[d] - 1) * (4 ÷ size(values, d)) for d = 1:ndims(values))
            total += values[K] * roots[mod(phase, 4)+1]
        end
        inverse ? total / length(values) : total
    end
end

# Index-based classification deliberately does not call the production mode helper.
null_index(I, shape) =
    all(d -> I[d] == 1 || (iseven(shape[d]) && I[d] == shape[d] ÷ 2 + 1), 1:length(shape))
function residual(g, rho_n, rho_np1, J, dt; transform = dft4, convert_value = complex_exact)
    delta = map((a, b) -> convert_value(a) - convert_value(b), rho_np1, rho_n)
    density_hat = transform(delta)
    current_hat = map(field -> transform(convert_value.(field)), J[1:length(g.n)])
    residual_hat = similar(density_hat)
    dc = first(CartesianIndices(density_hat))
    for I in CartesianIndices(density_hat)
        I != dc && null_index(I, size(density_hat)) && (density_hat[I] = zero(eltype(density_hat)))
        divergence =
            sum(im * convert_value(g.kvec[d][I[d]]) * current_hat[d][I] for d = 1:length(g.n))
        residual_hat[I] = density_hat[I] + convert_value(dt) * divergence
    end
    density = transform(density_hat; inverse = true)
    residual_values = transform(residual_hat; inverse = true)
    scale = maximum(z -> abs(real(z)), density)
    numerator = maximum(z -> abs(real(z)), residual_values)
    return numerator / (iszero(scale) ? one(scale) : scale), scale
end

function nearest_nonnegative(value::T, target::Q) where {T}
    target >= 0 || return false
    isinf(value) && return target >= exact(floatmax(T)) + exact(eps(floatmax(T))) / 2
    isfinite(value) && value >= 0 || return false
    iszero(value) && return target <= exact(nextfloat(zero(T))) / 2
    current = exact(value)
    below = (exact(prevfloat(value)) + current) / 2
    above =
        value == floatmax(T) ? current + exact(eps(value)) / 2 :
        (current + exact(nextfloat(value))) / 2
    even = iseven(reinterpret(T === Float32 ? UInt32 : UInt64, value))
    return (target > below || (even && target == below)) &&
           (target < above || (even && target == above))
end

function big_dft(values; inverse = false)
    direction = inverse ? 1 : -1
    return map(CartesianIndices(values)) do I
        total = zero(Complex{BigFloat})
        for K in CartesianIndices(values)
            phase = sum(BigFloat((K[d] - 1) * (I[d] - 1)) / size(values, d) for d = 1:ndims(values))
            total += values[K] * complex(cospi(2phase), direction * sinpi(2phase))
        end
        inverse ? total / length(values) : total
    end
end

physical(es::EMPIC) = (
    copy(es.rho_n),
    copy(es.rho_np1),
    map(copy, es.J),
    map(copy, es.Jraw),
    map(copy, es.E),
    map(copy, es.B),
)
physical(es::EMPIC1D) =
    (copy(es.rho_n), copy(es.rho_np1), copy(es.Jx), copy(es.Ex), copy(es.Ey), copy(es.Bz))
end

@testset "EM continuity diagnostic retains complete range-sensitive expression" begin
    H, O = HybridPlasmaPIC, ContinuityRangeOracle
    for T in (Float32, Float64), D = 1:3, axis = 1:D, sign in (-1, 1)
        shape = ntuple(d -> d == axis ? 4 : 2, D)
        M = floatmax(T)
        for kind in (
            :difference,
            :forward,
            :inverse,
            :product,
            :derivative_small,
            :derivative_large,
            :dc,
            :zero_dt,
            :negative_dt,
        )
            desired_k =
                kind == :derivative_small ? (T === Float32 ? T(1e-20) : T(1e-100)) :
                kind in (:derivative_large, :zero_dt) ? T(8) : one(T)
            g = FourierGrid(shape, ntuple(d -> d == axis ? T(2pi) / desired_k : T(2pi), D))
            es = EMPIC(g, 0; n0 = zero(T))
            a, b, j, dt = zero(T), zero(T), zero(T), one(T)
            if kind == :difference || kind == :dc
                a = T(0.75) * M
                b = -a
                dt = M
            elseif kind == :forward
                a = T(0.75) * M
            elseif kind == :inverse
                a = T(0.4) * M / T(2^(D - 1))
            elseif kind == :product || kind == :negative_dt
                a = -T(0.24) * M / T(2^(D - 1))
                j = T(2.75)
                dt = (kind == :negative_dt ? -one(T) : one(T)) * T(0.4) * M
            elseif kind == :derivative_small
                j = T === Float32 ? T(1e-30) : T(1e-300)
                dt = T === Float32 ? T(1e20) : T(1e100)
                a = -(j * dt) * g.kvec[axis][2]
            else
                j = T(0.2) * M
                dt = kind == :zero_dt ? zero(T) : T(1) / 16
                a = -T(0.1) * M
            end
            for I in CartesianIndices(es.rho_n)
                cosine = kind == :dc ? 1 : (1, 0, -1, 0)[I[axis]]
                es.rho_np1[I] = sign * a * cosine
                es.rho_n[I] = sign * b * cosine
                es.J[axis][I] = sign * j * (0, 1, 0, -1)[I[axis]]
            end
            reference, _ = O.residual(g, es.rho_n, es.rho_np1, es.J, dt)
            before = O.physical(es)
            got = charge_conservation_residual(es, dt)
            @test isfinite(got)
            @test O.nearest_nonnegative(got, reference)
            @test O.physical(es) == before
            if D == 1
                legacy = EMPIC1D(g, 0; n0 = zero(T))
                legacy.rho_n .= es.rho_n
                legacy.rho_np1 .= es.rho_np1
                legacy.Jx .= es.J[1]
                old = O.physical(legacy)
                @test charge_conservation_residual(legacy, dt) === got
                @test O.physical(legacy) == old
            end
        end
    end

    # Each coefficient contains widely separated signed bits, not just a
    # collection of differently scaled whole coefficients.
    for T in (Float32, Float64), D = 1:3
        shape = ntuple(d -> d == 1 ? 4 : 2, D)
        g = FourierGrid(shape, ntuple(_ -> T(2pi), D))
        M, tiny = O.exact(floatmax(T)), O.exact(nextfloat(zero(T)))
        values = fill(zero(O.CQ), shape)
        values[1] = complex(M, tiny)
        values[2] = complex(-M + tiny, -tiny)
        values[3] = complex(tiny, M - tiny)
        original = copy(values)
        forward = H._em_dyadic_transform(values, g)
        @test forward == O.dft4(values)
        @test H._em_dyadic_transform(forward, g, true) == values
        @test values == original
        @test all(iszero, H._em_dyadic_transform(fill(zero(O.CQ), shape), g))
    end

    for T in (Float32, Float64), D in (2, 3), kind in (:mixed_nyquist, :solenoidal, :pure_null)
        shape = ntuple(_ -> 4, D)
        g = FourierGrid(shape, ntuple(d -> d == 2 ? T(pi) : T(2pi), D))
        es = EMPIC(g, 0; n0 = zero(T))
        amplitude = T(0.0625) * floatmax(T)
        for I in CartesianIndices(es.rho_n)
            phase = mod(I[1] + I[2] - 2, 4) + 1
            if kind == :mixed_nyquist
                nyquist = isodd(I[1]) ? 1 : -1
                es.rho_np1[I] = -T(2) * amplitude * nyquist * (1, 0, -1, 0)[I[2]]
                es.J[2][I] = amplitude * nyquist * (0, 1, 0, -1)[I[2]]
            elseif kind == :solenoidal
                es.J[1][I] = amplitude * (0, 1, 0, -1)[phase]
                es.J[2][I] = -amplitude / 2 * (0, 1, 0, -1)[phase]
            else
                es.rho_np1[I] = amplitude * (iseven(I[1] + I[2]) ? 1 : -1)
            end
        end
        before = O.physical(es)
        @test charge_conservation_residual(es, one(T)) == zero(T)
        @test O.physical(es) == before
    end

    # Divergence uses only spatial components, even when unused components
    # are nonfinite. Both ordinary and range-sensitive paths keep that rule.
    for T in (Float32, Float64), D in (1, 2), amplitude in (one(T), T(0.75) * floatmax(T))
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
        es = EMPIC(g, 0; n0 = zero(T))
        fill!(es.rho_np1, amplitude)
        fill!(es.J[D+1], T(NaN))
        before = O.physical(es)
        @test charge_conservation_residual(es, one(T)) == one(T)
        @test isequal(O.physical(es), before)
    end

    for T in (Float32, Float64)
        # A genuinely unrepresentable final ratio remains infinity.
        g = FourierGrid((4,), (T(2pi),))
        es = EMPIC(g, 0; n0 = zero(T))
        tiny = nextfloat(zero(T))
        es.rho_np1 .= tiny .* T[1, 0, -1, 0]
        es.J[1] .= T[0, 2, 0, -2]
        @test charge_conservation_residual(es, one(T)) == T(Inf)
        # Null modes remain excluded, whereas a true DC density defect remains.
        fill!(es.J[1], zero(T))
        M = floatmax(T)
        es.rho_np1 .= T(0.75) * M .* T[1, -1, 1, -1]
        @test charge_conservation_residual(es, one(T)) == zero(T)
        fill!(es.rho_np1, T(0.75) * M)
        fill!(es.rho_n, -T(0.75) * M)
        @test charge_conservation_residual(es, one(T)) == one(T)
        # Nonfinite inputs do not enter the exact conversion owner.
        fill!(es.rho_n, zero(T))
        fill!(es.rho_np1, T(NaN))
        @test isnan(charge_conservation_residual(es, one(T)))
        fill!(es.rho_np1, one(T))
        @test isnan(charge_conservation_residual(es, T(Inf)))
        es.J[1][1] = T(Inf)
        @test_throws ArgumentError charge_conservation_residual(es, one(T))

        # Preserve the canonical operator's existing rejection boundaries,
        # even when the complete exact diagnostic would be finite.
        g2 = FourierGrid((4, 4), (T(2pi), T(2pi)))
        es2 = EMPIC(g2, 0; n0 = zero(T))
        fill!(es2.J[1], M)
        for I in CartesianIndices(es2.rho_n)
            es2.J[2][I] = tiny * (0, 1, 0, -1)[I[2]]
            es2.rho_np1[I] = -(M * tiny) * (1, 0, -1, 0)[I[2]]
        end
        before = O.physical(es2)
        @test_throws ArgumentError charge_conservation_residual(es2, M)
        @test O.physical(es2) == before
        k = T(0.4) * M
        g3 = FourierGrid((4, 4, 4), ntuple(_ -> T(2pi) / k, 3))
        es3 = EMPIC(g3, 0; n0 = zero(T))
        fill!(es3.rho_np1, one(T))
        @test all(all(isfinite, values) for values in g3.kvec)
        @test_throws ArgumentError charge_conservation_residual(es3, one(T))
    end
end

@testset "EM dyadic digits preserve signed band boundaries" begin
    H, O = HybridPlasmaPIC, ContinuityRangeOracle
    power_of_two(e) = e >= 0 ? (BigInt(1) << e) // BigInt(1) : BigInt(1) // (BigInt(1) << -e)
    for T in (Float32, Float64), shape in ((1,), (2,), (4,), (2, 4), (4, 2, 2))
        g = FourierGrid(shape, ntuple(_ -> T(2pi), length(shape)))
        width = max(1, precision(T) - (8sizeof(Int) - leading_zeros(prod(shape) - 1)) - 2)
        tiny = O.exact(nextfloat(zero(T)))
        # Adjacent binary places exercise carry, floor and mask boundaries.
        # The oracle uses exact Gaussian-integer roots, not radix extraction.
        for exponent in (-2width - 1, -width, -1, 0, width - 1, width, width + 1, 2width + 1)
            quantum = power_of_two(exponent)
            values = Array{O.CQ}(undef, shape)
            for (j, I) in enumerate(CartesianIndices(values))
                sign = isodd(j) ? -1 : 1
                values[I] = complex(
                    sign * ((BigInt(1) << width) + j) * quantum - tiny,
                    -sign * j * quantum + tiny,
                )
            end
            before = deepcopy(values)
            forward = H._em_dyadic_transform(values, g)
            @test forward == O.dft4(values)
            @test H._em_dyadic_transform(values, g, true) == O.dft4(values; inverse = true)
            @test H._em_dyadic_transform(forward, g, true) == values
            @test values == before
        end
        for invalid in (O.Q(1) / 3, BigInt(1) // BigInt(0))
            values = fill(zero(O.CQ), shape)
            values[end] = complex(invalid)
            before = deepcopy(values)
            @test_throws ArgumentError H._em_dyadic_transform(values, g)
            @test values == before
        end
    end
end

@testset "EM continuity normalized transforms on non-power-of-two grids" begin
    O = ContinuityRangeOracle
    for T in (Float32, Float64), D = 1:3, kind in (:difference, :weighted, :zero_dt)
        shape = ntuple(d -> d == 1 ? 5 : d == 2 ? 3 : 2, D)
        g = FourierGrid(shape, ntuple(_ -> T(2pi), D))
        es = EMPIC(g, 0; n0 = zero(T))
        M = floatmax(T)
        dt = kind == :weighted ? -T(0.4) * M : kind == :zero_dt ? zero(T) : one(T)
        for I in CartesianIndices(es.rho_n)
            cosine, sine = cospi(T(2) * T(I[1] - 1) / T(5)), sinpi(T(2) * T(I[1] - 1) / T(5))
            amplitude = (kind == :difference ? T(0.75) : T(0.2)) * M
            es.rho_np1[I] = amplitude * cosine
            es.rho_n[I] = kind == :difference ? -es.rho_np1[I] : zero(T)
            es.J[1][I] = (kind == :zero_dt ? T(0.2) * M : T(2.75)) * sine
        end
        before = O.physical(es)
        got = charge_conservation_residual(es, dt)
        @test isfinite(got)
        @test O.physical(es) == before
        setprecision(BigFloat, 2048) do
            reference, scale = O.residual(
                g,
                es.rho_n,
                es.rho_np1,
                es.J,
                dt;
                transform = O.big_dft,
                convert_value = x -> Complex{BigFloat}(x),
            )
            # The budget scales with operation count and the un-cancelled
            # density/current magnitudes; it is not fitted to observed errors.
            current_bound = sum(
                maximum(abs, BigFloat.(es.J[d])) * maximum(abs, BigFloat.(g.kvec[d])) for d = 1:D
            )
            allowance =
                64 * eps(T) * length(es.rho_n) * (1 + abs(BigFloat(dt)) * current_bound / scale)
            @test allowance < BigFloat(0.01)
            lower_precision = setprecision(BigFloat, 512) do
                first(
                    O.residual(
                        g,
                        es.rho_n,
                        es.rho_np1,
                        es.J,
                        dt;
                        transform = O.big_dft,
                        convert_value = x -> Complex{BigFloat}(x),
                    ),
                )
            end
            @test abs(reference - lower_precision) < BigFloat(eps(T)) / 1024
            @test abs(BigFloat(got) - reference) <= allowance
        end
        if D == 1
            old = EMPIC1D(g, 0; n0 = zero(T))
            old.rho_n .= es.rho_n
            old.rho_np1 .= es.rho_np1
            old.Jx .= es.J[1]
            @test charge_conservation_residual(old, dt) === got
        end
    end
end

@testset "EMPIC longitudinal projection preserves extreme represented modes" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64), D = 1:3, axis = 1:D, large_domain in (false, true)
        scale =
            T === Float32 ? (large_domain ? T(1e25) : T(1e-25)) :
            (large_domain ? T(1e200) : T(1e-200))
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(d -> d == axis ? T(2pi) * scale : T(4), D))
        electrons = ParticleSet{D,T}(4^(D - 1); q = -one(T), m = one(T))
        for (p, I) in enumerate(CartesianIndices(ntuple(d -> d == axis ? 1 : 4, D)))
            for d = 1:D
                electrons.x[d][p] = T(I[d] - 1) * g.dx[d]
            end
            electrons.weight[p] = prod(g.dx)
        end
        es = EMPIC(g, length(electrons.weight); n0 = T(0.25))
        init_empic!(es, electrons)
        amplitude = setprecision(BigFloat, 256) do
            T(-inv(BigFloat(g.kvec[axis][2])))
        end
        for I in CartesianIndices(es.E[axis])
            @test es.E[axis][I] ≈ amplitude / 2 * (zero(T), one(T), zero(T), -one(T))[I[axis]] rtol =
                16eps(T) atol = zero(T)
        end
        @test all(d -> d == axis || all(iszero, es.E[d]), 1:3)

        # A represented cosine charge change with zero current has residual
        # one regardless of whether its squared wave number is representable.
        fill!(es.rho_n, zero(T))
        for I in CartesianIndices(es.rho_np1)
            es.rho_np1[I] = (one(T), zero(T), -one(T), zero(T))[I[axis]]
        end
        @test charge_conservation_residual(es, one(T)) ≈ one(T) rtol = 16eps(T)
        H._correct_current_continuity!(es, one(T))
        for I in CartesianIndices(es.J[axis])
            @test es.J[axis][I] ≈ amplitude * (zero(T), one(T), zero(T), -one(T))[I[axis]] rtol =
                16eps(T) atol = zero(T)
        end
        @test charge_conservation_residual(es, one(T)) < 64eps(T)

        # Preserve source and dt separately: target=-delta/dt itself can lose
        # range while the complete longitudinal current remains normal.
        dt, change = if large_domain
            T === Float32 ? (T(1e30), T(1e-20)) : (T(1e300), T(1e-100))
        else
            T === Float32 ? (T(1e-40), one(T)) : (T(1e-320), one(T))
        end
        es.rho_np1 .*= change
        H._correct_current_continuity!(es, dt)
        expected = setprecision(BigFloat, 256) do
            T(-BigFloat(change) / (BigFloat(dt) * BigFloat(g.kvec[axis][2])))
        end
        @test isfinite(expected) && !iszero(expected)
        for I in CartesianIndices(es.J[axis])
            @test es.J[axis][I] ≈ expected * (zero(T), one(T), zero(T), -one(T))[I[axis]] rtol =
                16eps(T) atol = zero(T)
        end
    end
    for T in (Float32, Float64), overflow in (false, true)
        kscale = T === Float32 ? (overflow ? T(1e15) : T(1e-15)) : (overflow ? T(1e150) : T(1e-150))
        amplitude =
            T === Float32 ? (overflow ? T(1e30) : T(1e-35)) : (overflow ? T(1e200) : T(1e-200))
        g = FourierGrid((4,), (T(2pi) / kscale,))
        e = ParticleSet{1,T}(0; q = -one(T), m = one(T))
        es = EMPIC(g, 0; n0 = zero(T))
        es.E[1] .= amplitude .* T[1, 0, -1, 0]
        init_empic!(es, e)
        # The complete projection removes this purely longitudinal seed,
        # including when its intermediate divergence is outside T's range.
        @test all(iszero, es.E[1])
    end
end

@testset "EMPIC null modes depend on wave-vector components" begin
    for T in (Float32, Float64), D = 1:3
        scale = T === Float32 ? T(1e25) : T(1e200)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi) * scale, D))
        es = EMPIC(g, 0; n0 = zero(T))
        for c = 1:3, I in CartesianIndices(es.Jraw[c])
            es.Jraw[c][I] = T(c) + (isodd(I[1] - 1) ? -one(T) : one(T))
        end
        HybridPlasmaPIC._correct_current_continuity!(es, one(T))
        @test all(c -> all(==(T(c)), es.J[c]), 1:3)
        for axes = 1:((1<<D)-1)
            foreach(current -> fill!(current, zero(T)), es.J)
            for I in CartesianIndices(es.rho_np1)
                parity = sum(I[d] - 1 for d = 1:D if axes & (1 << (d - 1)) != 0)
                es.rho_np1[I] = isodd(parity) ? -one(T) : one(T)
            end
            @test charge_conservation_residual(es, one(T)) == zero(T)
            es.rho_np1 .+= one(T)
            @test charge_conservation_residual(es, one(T)) ≈ one(T) rtol = 16eps(T)
        end
        if D > 1
            # Nyquist in x with nonzero ky remains represented; its Jx is
            # transverse and must not be removed as a false null mode.
            fill!(es.rho_np1, zero(T))
            foreach(current -> fill!(current, zero(T)), es.Jraw)
            pattern = similar(es.rho_n)
            for I in CartesianIndices(pattern)
                pattern[I] =
                    (isodd(I[1] - 1) ? -one(T) : one(T)) * (one(T), zero(T), -one(T), zero(T))[I[2]]
            end
            es.Jraw[1] .= pattern
            HybridPlasmaPIC._correct_current_continuity!(es, one(T))
            @test es.J[1] == pattern
            @test all(c -> all(iszero, es.J[c]), 2:3)
            foreach(current -> fill!(current, zero(T)), es.J)
            es.rho_np1 .= pattern
            @test charge_conservation_residual(es, one(T)) ≈ one(T) rtol = 16eps(T)
        end
    end
end


@testset "EMPIC continuity diagnostic retains total charge errors" begin
    for T in (Float32, Float64), D = 1:3, n in (5, 6)
        g = FourierGrid(ntuple(_ -> n, D), ntuple(_ -> T(2π), D))
        es = EMPIC(g, 0; n0 = zero(T))
        dt = T(0.125)
        tolerance = 64eps(T)

        # Integral(div J) is zero on a periodic box. A spatially uniform
        # charge change therefore has normalized continuity residual one.
        for change in (-one(T), one(T))
            fill!(es.rho_np1, change)
            @test charge_conservation_residual(es, dt) ≈ one(T) atol = tolerance
        end

        # Independent sinusoidal derivative, with and without a global defect.
        # The constant defect cannot be canceled by the periodic current.
        wave_number = T(2π) / g.L[1]
        for I in CartesianIndices(es.rho_n)
            phase = T(2) * T(I[1] - 1) / T(n)
            es.J[1][I] = sinpi(phase)
            es.rho_np1[I] = -dt * wave_number * cospi(phase)
        end
        @test charge_conservation_residual(es, dt) < tolerance
        es.rho_np1 .+= T(0.5)
        expected = T(0.5) / maximum(abs, es.rho_np1)
        rho_before = (copy(es.rho_n), copy(es.rho_np1))
        current_before = map(copy, es.J)
        @test charge_conservation_residual(es, dt) ≈ expected atol = tolerance
        @test es.rho_n == rho_before[1]
        @test es.rho_np1 == rho_before[2]
        @test es.J == current_before
        if D == 1
            legacy = EMPIC1D(g, 0; n0 = zero(T))
            legacy.rho_np1 .= es.rho_np1
            legacy.Jx .= es.J[1]
            @test charge_conservation_residual(es, dt) ≈ charge_conservation_residual(legacy, dt) atol =
                tolerance
        end

        foreach(current -> fill!(current, zero(T)), es.J)
        if iseven(n)
            # All nonempty combinations of pure axis-Nyquist modes remain
            # outside the spectral derivative range, unlike the true DC mode.
            for axes = 1:((1<<D)-1)
                for I in CartesianIndices(es.rho_n)
                    parity = sum(I[d] - 1 for d = 1:D if (axes & (1 << (d - 1))) != 0)
                    es.rho_np1[I] = isodd(parity) ? -one(T) : one(T)
                end
                @test charge_conservation_residual(es, dt) == zero(T)
                es.rho_np1 .+= one(T)
                @test charge_conservation_residual(es, dt) ≈ one(T) atol = tolerance
            end
            if D > 1
                # Nyquist in one axis with a nonzero mode in another is a
                # represented mixed mode and must still contribute.
                for I in CartesianIndices(es.rho_n)
                    nyquist = isodd(I[1] - 1) ? -one(T) : one(T)
                    es.rho_np1[I] = nyquist * cospi(T(2) * T(I[2] - 1) / T(n))
                end
                @test charge_conservation_residual(es, dt) ≈ one(T) atol = tolerance
            end
        end
    end
end

function _empic_peakfreq(series, dt)
    Nt = length(series)
    S = fft(series)
    mag = abs.(S)
    best = 2
    bm = -1.0
    for j = 2:(Nt÷2)
        if mag[j+1] > bm
            bm = mag[j+1]
            best = j
        end
    end
    a = mag[best]
    b = mag[best+1]
    c = mag[best+2]
    d = a - 2b + c
    δ = d != 0 ? 0.5 * (a - c) / d : 0.0
    return 2π * (best + δ) / (Nt * dt)
end

function _empic_gauss_residual(es, rho)
    g = es.g
    rhohat = complex.(rho)
    fft!(rhohat)
    for I in CartesianIndices(rhohat)
        # Classify the derivative null space by Fourier indices, independently
        # of the production arithmetic and the representability of k squared.
        nullmode = all(d -> I[d] == 1 || iseven(g.n[d]) && I[d] == g.n[d] ÷ 2 + 1, 1:length(g.n))
        nullmode && (rhohat[I] = zero(eltype(rhohat)))
    end
    ifft!(rhohat)
    rhorepr = real.(rhohat)
    divE = similar(rho)
    divergence!(divE, es.E, g)
    scale = max(maximum(abs, rhorepr), 1.0)
    return maximum(abs, divE .- rhorepr) / scale
end

function _empic_continuity_mode_oracle(es, dt)
    increment_hat = fft(es.rho_np1 .- es.rho_n)
    dc_increment = real(first(increment_hat)) / length(increment_hat)
    # Use Fourier indices independently of the production k² filter. These
    # are exactly the DC/pure-Nyquist combinations, not mixed nonzero modes.
    for I in CartesianIndices(increment_hat)
        if all(d -> I[d] == 1 || (iseven(es.g.n[d]) && I[d] == es.g.n[d] ÷ 2 + 1), 1:length(es.g.n))
            increment_hat[I] = 0
        end
    end
    non_dc_increment = real.(ifft(increment_hat))
    current_divergence = similar(es.rho_n)
    divergence!(current_divergence, es.J, es.g)
    non_dc_scale = maximum(abs, non_dc_increment)
    non_dc_residual = maximum(abs, non_dc_increment .+ dt .* current_divergence)
    full_scale = maximum(abs, non_dc_increment .+ dc_increment)
    return (; non_dc_scale, non_dc_residual, full_scale)
end

@testset "EMPIC constructor validation" begin
    g2 = FourierGrid((8, 6), (2π, 3π))
    es2 = EMPIC(g2, 5; n0 = 0.75, c = 4.0, mi = 100.0)
    @test es2.n0 == 0.75
    @test es2.c == 4.0
    @test es2.mi == 100.0
    @test size(es2.E[1]) == (8, 6)
    @test all(size(es2.B[c]) == (8, 6) for c = 1:3)
    @test all(length(es2.Ep[c]) == 5 for c = 1:3)
    @test all(length(es2.mide[d]) == 5 for d = 1:2)

    g3 = FourierGrid((6, 5, 4), (2π, 3π, 4π))
    es3 = EMPIC(g3, 7; n0 = 1.25)
    @test size(es3.rho_n) == (6, 5, 4)
    @test all(size(es3.J[c]) == (6, 5, 4) for c = 1:3)
    @test all(length(es3.mide[d]) == 7 for d = 1:3)

    @test_throws ArgumentError EMPIC(g2, 5; n0 = NaN)
    @test_throws ArgumentError EMPIC(g2, 5; n0 = -1.0)
    @test_throws ArgumentError EMPIC(g2, 5; c = NaN)
    @test_throws ArgumentError EMPIC(g2, 5; c = 0.0)
    @test_throws ArgumentError EMPIC(g2, 5; mi = 0.0)
    @test_throws ArgumentError EMPIC(g2, -1)

    g4 = FourierGrid((2, 2, 2, 2), (1.0, 1.0, 1.0, 1.0))
    @test_throws ArgumentError EMPIC(g4, 1)
end

@testset "EMPIC validates species and timestep before mutation" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 2π))
    e = ParticleSet{2,T}(48; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (8, 6))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, 48; n0 = 1.0, c = 5.0)
    init_empic!(es, e)

    x0 = ntuple(d -> copy(e.x[d]), 2)
    v0 = ntuple(c -> copy(e.v[c]), 3)
    E0 = ntuple(c -> copy(es.E[c]), 3)
    B0 = ntuple(c -> copy(es.B[c]), 3)
    rho_n0 = copy(es.rho_n)
    rho_np10 = copy(es.rho_np1)
    time0 = es.time[]
    step0 = es.step[]

    @test_throws ArgumentError step_empic!(es, e, NaN)
    @test_throws ArgumentError step_empic!(es, e, -0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.rho_n == rho_n0
    @test es.rho_np1 == rho_np10
    @test es.time[] == time0
    @test es.step[] == step0

    e.q = 0.0
    @test_throws ArgumentError step_empic!(es, e, 0.1)
    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.time[] == time0
    @test es.step[] == step0
end

@testset "EMPIC rejects invalid relativistic velocities before mutation" begin
    g = FourierGrid((4, 4), (2π, 2π))
    e = ParticleSet{2,Float64}(16; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (4, 4))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, 16; c = 5.0, relativistic = true)
    init_empic!(es, e)
    e.v[1][end] = 5.5
    x0 = map(copy, e.x)
    v0 = map(copy, e.v)
    E0 = map(copy, es.E)
    B0 = map(copy, es.B)
    rho_n0, rho_np10 = copy(es.rho_n), copy(es.rho_np1)
    time0, step0 = es.time[], es.step[]

    @test_throws ArgumentError step_empic!(es, e, 0.01)

    @test all(e.x[d] == x0[d] for d = 1:2)
    @test all(e.v[c] == v0[c] for c = 1:3)
    @test all(es.E[c] == E0[c] for c = 1:3)
    @test all(es.B[c] == B0[c] for c = 1:3)
    @test es.rho_n == rho_n0 && es.rho_np1 == rho_np10
    @test es.time[] == time0 && es.step[] == step0
end

@testset "EMPIC init preserves transverse field and enforces Gauss law" begin
    T = Float64
    nx, ny = 16, 8
    L = 2π
    g = FourierGrid((nx, ny), (L, L))
    N = nx * ny
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (nx, ny))
    set_density_weight!(e, 1.0, g)
    es = EMPIC(g, N; n0 = 1.0, c = 5.0)
    seed = zeros(T, nx, ny)
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        seed[i, j] = 1e-3 * cos(x)
        es.E[2][i, j] = seed[i, j]
    end

    init_empic!(es, e)
    @test maximum(abs, es.E[2] .- seed) < 1e-12
    @test _empic_gauss_residual(es, es.rho_n) < 1e-12
end

@testset "EMPIC 2D transverse dispersion ω²=ωpe²+c²k²" begin
    T = Float64
    nx, ny = 32, 4
    L = 2π
    cc = 5.0
    n0 = 1.0
    g = FourierGrid((nx, ny), (L, L))
    nppc = 80
    N = nppc * nx * ny
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    load_lattice!(e, (0.0, 0.0), g.L, (nx * nppc, ny))
    set_density_weight!(e, n0, g)
    es = EMPIC(g, N; n0 = n0, c = cc, shape = CIC())
    m = 1
    k = 2π * m / L
    for j = 1:ny, i = 1:nx
        x = (i - 1) * g.dx[1]
        es.E[2][i, j] = 1e-3 * cos(k * x)
    end
    init_empic!(es, e)
    dt = 0.01
    series = ComplexF64[]
    for _ = 1:1000
        step_empic!(es, e, dt)
        push!(series, mode_amplitude(es.E[2], g, (m, 0)))
    end
    ω = _empic_peakfreq(series, dt)
    ωth = sqrt(n0 + cc^2 * k^2)
    @test abs(ω - ωth) / ωth < 0.03
    modes = _empic_continuity_mode_oracle(es, dt)
    # Preserve the original relative bound on the non-DC modes. A nearly
    # transverse wave has tiny Δρ, so density-scale rounding in total charge
    # must instead be checked against the arithmetic that deposits that charge.
    @test modes.non_dc_residual / modes.non_dc_scale < 1e-8
    # Positive CIC accumulation has at most N terms per cell. The remaining
    # budget covers local shape/normalization arithmetic and a conservative
    # eight operations per grid point for the Fourier projection. γ_k bounds
    # accumulated relative rounding; two independently rounded snapshots enter Δρ.
    operation_budget = N + 8length(es.rho_n) + 16
    ku = operation_budget * (eps(T) / 2)
    γ = ku / (1 - ku)
    density_scale = max(n0, maximum(abs, es.ne), maximum(abs, n0 .- es.rho_n))
    allowance = 2γ * density_scale + 1e-8 * modes.non_dc_scale
    @test charge_conservation_residual(es, dt) * modes.full_scale < allowance
    @test _empic_gauss_residual(es, es.rho_np1) < 1e-8

    # A resolved total-charge defect leaves the old non-DC guarantee intact
    # but must fail the added full-residual check.
    dc_defect = sqrt(eps(T)) * density_scale
    @test dc_defect > allowance
    saved_charge = copy(es.rho_np1)
    try
        es.rho_np1 .+= dc_defect
        defective_modes = _empic_continuity_mode_oracle(es, dt)
        @test defective_modes.non_dc_residual / defective_modes.non_dc_scale < 1e-8
        @test charge_conservation_residual(es, dt) * defective_modes.full_scale > allowance
    finally
        es.rho_np1 .= saved_charge
    end
end

@testset "EMPIC 2D/3D spectral current correction conserves charge" begin
    T = Float64

    rng2 = MersenneTwister(11)
    g2 = FourierGrid((16, 12), (2π, 2π))
    e2 = ParticleSet{2,T}(1200; q = -1.0, m = 1.0)
    load_uniform!(e2, rng2, (0.0, 0.0), g2.L)
    set_density_weight!(e2, 1.0, g2)
    for p in eachindex(e2.weight)
        e2.v[1][p] = 0.2 * randn(rng2)
        e2.v[2][p] = 0.2 * randn(rng2)
        e2.v[3][p] = 0.1 * randn(rng2)
    end
    es2 = EMPIC(g2, nparticles(e2); n0 = 1.0, c = 5.0)
    init_empic!(es2, e2)
    maxres2 = 0.0
    for _ = 1:20
        step_empic!(es2, e2, 0.005)
        maxres2 = max(maxres2, charge_conservation_residual(es2, 0.005))
    end
    @test maxres2 < 1e-10
    @test _empic_gauss_residual(es2, es2.rho_np1) < 1e-10

    rng3 = MersenneTwister(13)
    g3 = FourierGrid((8, 6, 4), (2π, 2π, 2π))
    e3 = ParticleSet{3,T}(384; q = -1.0, m = 1.0)
    load_uniform!(e3, rng3, (0.0, 0.0, 0.0), g3.L)
    set_density_weight!(e3, 1.0, g3)
    for p in eachindex(e3.weight)
        e3.v[1][p] = 0.15 * randn(rng3)
        e3.v[2][p] = 0.15 * randn(rng3)
        e3.v[3][p] = 0.15 * randn(rng3)
    end
    es3 = EMPIC(g3, nparticles(e3); n0 = 1.0, c = 4.0)
    init_empic!(es3, e3)
    maxres3 = 0.0
    for _ = 1:10
        step_empic!(es3, e3, 0.003)
        maxres3 = max(maxres3, charge_conservation_residual(es3, 0.003))
    end
    @test maxres3 < 1e-10
    @test _empic_gauss_residual(es3, es3.rho_np1) < 1e-10
end

@testset "EMPIC 2D mobile-ion path" begin
    T = Float64
    g = FourierGrid((8, 6), (2π, 2π))
    N = prod(g.n)
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0, 0.0), g.L, g.n)
    load_lattice!(ions, (0.0, 0.0), g.L, g.n)
    set_density_weight!(e, 1.0, g)
    set_density_weight!(ions, 1.0, g)
    es = EMPIC(g, N; mobile = true, mi = 25.0, c = 4.0)
    @test_throws ArgumentError init_empic!(es, e)

    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.01)
    @test charge_conservation_residual(es, 0.01) < 1e-12
    @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
    @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))

    badions = ParticleSet{2,T}(N; q = 0.0, m = 25.0)
    load_lattice!(badions, (0.0, 0.0), g.L, g.n)
    set_density_weight!(badions, 1.0, g)
    @test_throws ArgumentError step_empic!(es, e, badions, 0.01)
end

@testset "EMPIC resizes electron gather buffers when the electron count grows" begin
    T = Float64
    g = FourierGrid((16, 16), (2π, 2π))
    N = 64
    e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
    ions = ParticleSet{2,T}(N; q = 1.0, m = 25.0)
    load_lattice!(e, (0.0, 0.0), g.L, (8, 8))
    load_lattice!(ions, (0.0, 0.0), g.L, (8, 8))
    es = EMPIC(g, N; mobile = true, mi = 25.0, c = 8.0)
    init_empic!(es, e, ions)
    step_empic!(es, e, ions, 0.005)
    # grow the electron population (as ionize_mcc! secondaries would); the next step must not
    # DimensionMismatch on stale-sized electron gather buffers (they are ion-symmetric now).
    xe = ParticleSet{2,T}(5; q = -1.0, m = 1.0)
    for d = 1:2, k = 1:5
        xe.x[d][k] = 0.1k
    end
    HybridPlasmaPIC.append_particles!(e, xe)
    @test HybridPlasmaPIC.nparticles(e) == N + 5
    step_empic!(es, e, ions, 0.005)
    @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))
end

@testset "EMPIC-order: step_empic! is 2nd-order in dt (leapfrog v + B priming)" begin
    # The multi-dim EM PIC primes both species' velocities and the seeded B on the first step
    # (same fix as EMPIC1D). Cold 2D longitudinal Langmuir self-convergence on Ex — rate ≈ 2.
    T = Float64
    L = 2π
    k = 2π / L
    A = 0.02
    function ex(nsteps, Tf)
        n = 16
        pd = 80
        N = pd * pd
        g = FourierGrid((n, n), (T(L), T(L)))
        e = ParticleSet{2,T}(N; q = -1.0, m = 1.0)
        ions = ParticleSet{2,T}(N; q = 1.0, m = 100.0)
        load_lattice!(e, (0.0, 0.0), g.L, (pd, pd))
        load_lattice!(ions, (0.0, 0.0), g.L, (pd, pd))
        set_density_weight!(e, 1.0, g)
        set_density_weight!(ions, 1.0, g)
        for p = 1:N
            e.x[1][p] = mod(e.x[1][p] - (A / k) * sin(k * e.x[1][p]), L)
        end
        es = EMPIC(g, N; mobile = true, mi = 100.0, c = 8.0)
        init_empic!(es, e, ions)
        for _ = 1:nsteps
            step_empic!(es, e, ions, Tf / nsteps)
        end
        abs(mode_amplitude(es.E[1], g, (1, 0)))
    end
    Tf = 1.0
    seq = (40, 80, 160, 320)
    v = [ex(ns, Tf) for ns in seq]
    r1 = log2(abs(v[1] - v[2]) / abs(v[2] - v[3]))
    r2 = log2(abs(v[2] - v[3]) / abs(v[3] - v[4]))
    @test r1 > 1.6
    @test r2 > 1.6
end

@testset "EMPIC ↔ EMPIC1D spectral equivalence incl. Nyquist current exclusion" begin
    # The continuity correction zeroes the pure-Nyquist current modes (as
    # _esirkepov_Jx! / _deposit_Jy_at_midpoints! do in EMPIC1D), so the two 1D
    # solvers must agree to roundoff on EVERY Fourier mode and Ex must never
    # accumulate a grid-Nyquist sawtooth. Pre-fix: the raw |Ĵx_Nyq| survived
    # only in EMPIC (1.5e-4 after one step here), Ex's Nyquist mode
    # random-walked to a few % of max|Ex|, and the solvers drifted apart by
    # ~2.5% within 200 steps.
    T = Float64
    n = 32
    L = 2π
    N = 1280
    dt = 0.01
    g1 = FourierGrid((n,), (L,))
    g2 = FourierGrid((n,), (L,))
    rng = MersenneTwister(1234)
    e1 = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    load_uniform!(e1, rng, (0.0,), (L,))
    set_density_weight!(e1, 1.0, g1)
    for p = 1:N
        e1.v[1][p] = 0.1 * randn(rng)
        e1.v[2][p] = 0.1 * randn(rng)     # exercise the transverse Jy channel
    end
    e2 = ParticleSet{1,T}(N; q = -1.0, m = 1.0)
    e2.x[1] .= e1.x[1]
    for c = 1:3
        e2.v[c] .= e1.v[c]
    end
    e2.weight .= e1.weight
    es1 = EMPIC1D(g1, N; n0 = 1.0, c = 5.0, shape = CIC())
    es2 = EMPIC(g2, N; n0 = 1.0, c = 5.0, shape = CIC())
    init_empic!(es1, e1)
    init_empic!(es2, e2)

    nyq(f) = abs(fft(complex.(f))[length(f)÷2+1]) / length(f)

    # one step: current spectra agree on ALL modes; Nyquist current is zero
    step_empic!(es1, e1, dt)
    step_empic!(es2, e2, dt)
    @test maximum(abs.(fft(complex.(es1.Jx)) .- fft(complex.(es2.J[1])))) / n < 1e-14
    @test maximum(abs.(fft(complex.(es1.Jy)) .- fft(complex.(es2.J[2])))) / n < 1e-14
    @test nyq(es2.J[1]) < 1e-15
    @test nyq(es2.J[2]) < 1e-15

    # 300 steps: solvers stay equivalent to roundoff and the Nyquist mode of
    # Ex stays at machine zero (measured 9e-15 / 5.5e-16 with the fix; the
    # pre-fix values 2.5e-2 / ≳1e-3 fail these bounds by many orders)
    nyqmax = 0.0
    for _ = 2:300
        step_empic!(es1, e1, dt)
        step_empic!(es2, e2, dt)
        nyqmax = max(nyqmax, nyq(es2.E[1]) / maximum(abs.(es2.E[1])))
    end
    @test maximum(abs.(es1.Ex .- es2.E[1])) / maximum(abs.(es1.Ex)) < 1e-12
    @test maximum(abs.(es1.Ey .- es2.E[2])) / maximum(abs.(es1.Ey)) < 1e-12
    @test nyqmax < 1e-13
    @test charge_conservation_residual(es2, dt) < 1e-10
end

@testset "EM PIC changing timestep preserves electron and ion half levels" begin
    # Constant E gives an exact momentum oracle, including relativistic motion.
    # Zero-weight tracers leave the imposed field unchanged through Maxwell.
    for (constructor, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
        relativistic in (false, true),
        mobile in (false, true)

        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0, D))
        e = ParticleSet{D,Float64}(1; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(1; q = 1.0, m = 16.0)
        e.weight .= 0.0
        ions.weight .= 0.0
        for d = 1:D
            e.x[d] .= 1.0
            ions.x[d] .= 2.0
        end
        es = constructor(g, 1; n0 = 0.0, c = 2.0, mobile, mi = 16.0, relativistic)
        species = mobile ? ions : nothing
        init_empic!(es, e, species)
        function set_field!()
            if es isa EMPIC1D
                es.Ex .= -0.7
                es.Ey .= 0.2
            else
                es.E[1] .= -0.7
                es.E[2] .= 0.2
            end
        end
        set_field!()
        elapsed = 0.0
        for (dt, n_sub) in ((0.1, 2), (0.2, 3), (0.05, 1))
            es isa EMPIC1D && (es.n_sub = n_sub)
            previous_dt = es.last_dt[]
            step_empic!(es, e, species, 0.0)
            @test es.last_dt[] == previous_dt
            step_empic!(es, e, species, dt)
            elapsed += dt
            electron_dt = es isa EMPIC1D ? dt / n_sub : dt
            for (s, acceleration, h) in
                ((e, (0.7, -0.2, 0.0), electron_dt), (ions, (-0.7 / 16, 0.2 / 16, 0.0), dt))
                s === ions && !mobile && continue
                u = ntuple(c -> acceleration[c] * (elapsed - h / 2), 3)
                gamma = relativistic ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0
                @test [s.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 2e-15
                if !relativistic
                    for d = 1:D
                        x0 = s === e ? 1.0 : 2.0
                        @test s.x[d][1] ≈ x0 + acceleration[d] * elapsed^2 / 2 atol = 2e-15
                    end
                end
            end
        end
        before = (deepcopy(e.x), deepcopy(e.v), es.time[], es.step[], es.last_dt[])
        for dt in (NaN, Inf, -0.1)
            @test_throws ArgumentError step_empic!(es, e, species, dt)
            @test (e.x, e.v, es.time[], es.step[], es.last_dt[]) == before
        end
        # Reinitialization takes freshly loaded physical velocities; a previous
        # run's unequal timesteps must not leak into the next run's prime.
        for c = 1:3
            e.v[c] .= 0.0
            ions.v[c] .= 0.0
        end
        init_empic!(es, e, species)
        @test es.last_dt[] == 0.0
        es isa EMPIC1D && @test es.last_sub_dt[] == 0.0
        set_field!()
        step_empic!(es, e, species, 0.08)
        u = (0.7 * 0.04, -0.2 * 0.04, 0.0)
        gamma = relativistic ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0
        @test [e.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 2e-15
    end
end

@testset "EM PIC changing timestep preserves Maxwell time staggering" begin
    function variable_dt_vacuum_error(constructor, refinement)
        g = FourierGrid((16,), (2π,))
        e = ParticleSet{1,Float64}(0; q = -1.0, m = 1.0)
        es = constructor(g, 0; n0 = 0.0, c = 2.0)
        init_empic!(es, e)
        E = es isa EMPIC1D ? es.Ey : es.E[2]
        B = es isa EMPIC1D ? es.Bz : es.B[3]
        x = range(0.0; step = g.dx[1], length = g.n[1])
        E .= sin.(x)
        for (count, dt) in ((10refinement, 0.04 / refinement), (20refinement, 0.02 / refinement))
            for _ = 1:count
                step_empic!(es, e, dt)
            end
        end
        # Analytic vacuum mode: Ey=sin(x)cos(ct), Bz=-cos(x)sin(ct)/c.
        # E is at t=0.8, while the stored B trails it by the last dt/2.
        Eerror = maximum(abs, E .- cos(1.6) .* sin.(x))
        Berror = maximum(abs, B .+ sin(2 * (0.8 - 0.01 / refinement)) / 2 .* cos.(x))
        return (Eerror, Berror)
    end
    for constructor in (EMPIC1D, EMPIC)
        errors = [variable_dt_vacuum_error(constructor, r) for r in (1, 2, 4)]
        for component = 1:2, level = 1:2
            @test log2(errors[level][component] / errors[level+1][component]) > 1.8
        end
    end
end

@testset "EM PIC changing timestep gyro-orbit convergence" begin
    function variable_dt_gyro_error(constructor, relativistic, refinement)
        g = FourierGrid((8,), (10.0,))
        e = ParticleSet{1,Float64}(1; q = -1.0, m = 1.0)
        e.weight .= 0.0
        e.x[1] .= 1.0
        e.v[1] .= 1.0
        es = constructor(g, 1; n0 = 0.0, c = 2.0, relativistic)
        init_empic!(es, e)
        if es isa EMPIC1D
            es.Bz .= 1.0
        else
            es.B[3] .= 1.0
        end
        for (count, dt) in ((10refinement, 0.04 / refinement), (20refinement, 0.02 / refinement))
            for _ = 1:count
                step_empic!(es, e, dt)
            end
        end
        omega = relativistic ? sqrt(1 - 1 / es.c^2) : 1.0
        return abs(e.x[1][1] - (1 + sin(omega * 0.8) / omega))
    end
    for constructor in (EMPIC1D, EMPIC), relativistic in (false, true)
        errors = [variable_dt_gyro_error(constructor, relativistic, r) for r in (1, 2, 4)]
        @test all(log2(errors[j] / errors[j+1]) > 1.8 for j = 1:2)
    end
end

@testset "EM PIC changing timestep retains charge and Gauss constraints" begin
    for (constructor, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3))
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> 2π, D))
        e = ParticleSet{D,Float64}(8^D; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(8^D; q = 1.0, m = 16.0)
        for s in (e, ions)
            load_lattice!(s, ntuple(_ -> 0.0, D), g.L, ntuple(_ -> 8, D))
            set_density_weight!(s, 1.0, g)
        end
        @. e.x[1] -= 0.02 * sin(e.x[1])
        e.v[1] .= 0.05
        ions.v[1] .= 0.02
        es = constructor(g, nparticles(e); mobile = true, mi = 16.0, c = 2.0)
        init_empic!(es, e, ions)
        for (dt, n_sub) in ((0.01, 1), (0.005, 3), (0.02, 2), (0.01, 1))
            es isa EMPIC1D && (es.n_sub = n_sub)
            step_empic!(es, e, ions, dt)
            sub_dt = es isa EMPIC1D ? dt / n_sub : dt
            delta_hat = fft(es.rho_np1 .- es.rho_n)
            dc = first(CartesianIndices(delta_hat))
            for I in CartesianIndices(delta_hat)
                if I != dc && all(iszero(g.kvec[d][I[d]]) for d = 1:D)
                    delta_hat[I] = 0
                end
            end
            delta = real.(ifft(delta_hat))
            divJ = similar(delta)
            if es isa EMPIC1D
                deriv!(divJ, es.Jx, g, 1)
            else
                divergence!(divJ, es.J, g)
            end
            # Density is normalized to one. Keep the DC conservation check,
            # whose deposition roundoff can dominate a tiny density increment;
            # bound absolute charge error rather than divide it by that increment.
            @test maximum(abs, delta .+ sub_dt .* divJ) < 32eps(Float64)
            @test charge_conservation_residual(es, sub_dt) * maximum(abs, delta) < 32eps(Float64)
            if es isa EMPIC1D
                rhohat = fft(es.rho_np1)
                rhohat[1] = rhohat[g.n[1]÷2+1] = 0
                rho = real.(ifft(rhohat))
                divE = similar(rho)
                deriv!(divE, es.Ex, g, 1)
                @test maximum(abs, divE .- rho) < 1e-12
            else
                @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
            end
        end
    end
end

function _empic_subcycling_fixture(
    D;
    n_sub = 1,
    shape = CIC(),
    mobile = false,
    relativistic = false,
    mi = 16.0,
)
    g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> 2π, D))
    e = ParticleSet{D,Float64}(8^D; q = -1.0, m = 1.0)
    ions = ParticleSet{D,Float64}(8^D; q = 1.0, m = mi)
    for s in (e, ions)
        load_lattice!(s, ntuple(_ -> 0.0, D), g.L, g.n)
        set_density_weight!(s, 1.0, g)
    end
    @. e.x[1] -= 0.02 * sin(e.x[1])
    @. e.v[1] = 0.1 + 0.02 * cos(e.x[1])
    @. e.v[2] = 0.01 * sin(e.x[1])
    ions.v[1] .= 0.03
    es = EMPIC(g, nparticles(e); n_sub, shape, mobile, relativistic, mi, c = 3.0)
    for I in CartesianIndices(es.E[1])
        x = (I[1] - 1) * g.dx[1]
        es.E[2][I] = 0.01 * sin(x)
        es.E[3][I] = 0.02 * cos(x)
        es.B[2][I] = 0.01 * cos(x)
        es.B[3][I] = 0.02 * sin(x)
    end
    init_empic!(es, e, mobile ? ions : nothing)
    return es, e, ions
end

@testset "EMPIC subcycling validates and resets timing state" begin
    for D = 1:3
        es, e, ions = _empic_subcycling_fixture(D; n_sub = 3, mobile = true)
        @test es.n_sub == 3
        @test_throws ArgumentError EMPIC(es.g, 1; n_sub = 0)
        @test_throws ArgumentError EMPIC(es.g, 1; n_sub = -1)
        step_empic!(es, e, ions, 0.03)
        before = (
            deepcopy(e.x),
            deepcopy(e.v),
            deepcopy(ions.x),
            deepcopy(ions.v),
            deepcopy(es.E),
            deepcopy(es.B),
            es.time[],
            es.step[],
            es.last_dt[],
            es.last_sub_dt[],
        )
        step_empic!(es, e, ions, 0.0)
        @test (
            e.x,
            e.v,
            ions.x,
            ions.v,
            es.E,
            es.B,
            es.time[],
            es.step[],
            es.last_dt[],
            es.last_sub_dt[],
        ) == before
        for n_sub in (0, -1)
            es.n_sub = n_sub
            @test_throws ArgumentError step_empic!(es, e, ions, 0.03)
            @test (
                e.x,
                e.v,
                ions.x,
                ions.v,
                es.E,
                es.B,
                es.time[],
                es.step[],
                es.last_dt[],
                es.last_sub_dt[],
            ) == before
        end
        es.n_sub = 3
        @test_throws ArgumentError step_empic!(es, e, ions, nextfloat(0.0))
        @test (
            e.x,
            e.v,
            ions.x,
            ions.v,
            es.E,
            es.B,
            es.time[],
            es.step[],
            es.last_dt[],
            es.last_sub_dt[],
        ) == before
        init_empic!(es, e, ions)
        @test (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]) == (0.0, 0, 0.0, 0.0)
    end
end

@testset "EMPIC immobile subcycling equals explicit electron steps" begin
    for D = 1:3, shape in (CIC(), TSC()), relativistic in (false, true), ns in (1, 2, 4)
        es, e, _ = _empic_subcycling_fixture(D; n_sub = ns, shape, relativistic)
        fine, efine, _ = _empic_subcycling_fixture(D; shape, relativistic)
        for (dt, n_sub) in ((0.01, ns), (0.02, 3), (0.015, ns))
            es.n_sub = n_sub
            step_empic!(es, e, dt)
            for _ = 1:n_sub
                step_empic!(fine, efine, dt / n_sub)
            end
            @test e.x == efine.x
            @test e.v == efine.v
            @test es.E == fine.E
            @test es.B == fine.B
            @test es.J == fine.J
            @test es.time[] == fine.time[]
            @test es.last_sub_dt[] == fine.last_sub_dt[]
        end
        @test es.step[] == 3
    end
end

@testset "EMPIC subcycling agrees with polarized EMPIC1D" begin
    for shape in (CIC(), TSC()), relativistic in (false, true), mobile in (false, true), ns = 1:4
        g = FourierGrid((16,), (2π,))
        e = ParticleSet{1,Float64}(256; q = -1.0, m = 1.0)
        ions = ParticleSet{1,Float64}(256; q = 1.0, m = 16.0)
        for s in (e, ions)
            load_lattice_1d!(s, 0.0, g.L[1])
            set_density_weight!(s, 1.0, g)
        end
        @. e.x[1] -= 0.02 * sin(e.x[1])
        @. e.v[2] = 0.02 * cos(e.x[1])
        e.v[1] .= relativistic ? 1.0 : 0.1
        es = EMPIC(g, nparticles(e); n_sub = ns, shape, relativistic, mobile, mi = 16.0, c = 3.0)
        ref = EMPIC1D(
            FourierGrid(g.n, g.L),
            nparticles(e);
            n_sub = ns,
            shape,
            relativistic,
            mobile,
            mi = 16.0,
            c = 3.0,
        )
        eref, iref = deepcopy(e), deepcopy(ions)
        for j = 1:g.n[1]
            x = (j - 1) * g.dx[1]
            es.E[2][j] = ref.Ey[j] = 0.01 * sin(x)
            es.B[3][j] = ref.Bz[j] = 0.02 * cos(x)
        end
        init_empic!(es, e, mobile ? ions : nothing)
        init_empic!(ref, eref, mobile ? iref : nothing)
        for (dt, n_sub) in ((0.01, ns), (0.02, 3), (0.015, 2))
            es.n_sub = ref.n_sub = n_sub
            step_empic!(es, e, mobile ? ions : nothing, dt)
            step_empic!(ref, eref, mobile ? iref : nothing, dt)
            @test all(isapprox.(e.x, eref.x; atol = 2e-13))
            @test all(isapprox.(e.v, eref.v; atol = 2e-13))
            @test all(isapprox.(ions.x, iref.x; atol = 2e-13))
            @test all(isapprox.(ions.v, iref.v; atol = 2e-13))
            @test es.E[1] ≈ ref.Ex atol = 2e-13
            @test es.E[2] ≈ ref.Ey atol = 2e-13
            @test es.B[3] ≈ ref.Bz atol = 2e-13
            @test es.time[] == ref.time[]
        end
    end
end

@testset "EMPIC multidimensional mobile subcycling conserves charge after population growth" begin
    for D in (2, 3), shape in (CIC(), TSC())
        es, e, ions = _empic_subcycling_fixture(D; n_sub = 3, shape, mobile = true)
        for (dt, ns) in ((0.01, 3), (0.02, 2), (0.005, 4))
            es.n_sub = ns
            step_empic!(es, e, ions, dt)
            @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
            delta_scale = maximum(abs, es.rho_np1 .- es.rho_n)
            @test charge_conservation_residual(es, dt / ns) * delta_scale < 128eps(Float64)
            @test isfinite(em_field_energy(es) + kinetic_energy(e) + kinetic_energy(ions))
            # An equal-position, equal-weight pair leaves the total charge
            # unchanged at insertion and exercises both lazy workspace resizes.
            for s in (e, ions)
                extra = ParticleSet{D,Float64}(2; q = s.q, m = s.m)
                for d = 1:D
                    extra.x[d] .= (0.1d, 0.2d)
                end
                extra.weight .= 0.001
                HybridPlasmaPIC.append_particles!(s, extra)
            end
        end
        step_empic!(es, e, ions, 0.005)
        @test all(length(v) == nparticles(e) for v in (es.Ep..., es.Bp..., es.mide..., es.worke))
        @test all(
            length(v) == nparticles(ions) for v in (es.Epi..., es.Bpi..., es.midi..., es.worki)
        )
        @test _empic_gauss_residual(es, es.rho_np1) < 1e-12
    end
end

@testset "EMPIC multidimensional mobile subcycling has second-order temporal convergence" begin
    function subcycling_amplitude(D, mi, transverse, refinement)
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> 2π, D))
        e = ParticleSet{D,Float64}(16^D; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(16^D; q = 1.0, m = mi)
        for s in (e, ions)
            load_lattice!(s, ntuple(_ -> 0.0, D), g.L, ntuple(_ -> 16, D))
            set_density_weight!(s, 1.0, g)
        end
        es = EMPIC(g, nparticles(e); mobile = true, mi, c = 2.0, n_sub = 3)
        if transverse
            # The polarization (1,-1,0) is orthogonal to k=(1,...,1).
            for I in CartesianIndices(es.E[1])
                phase = sum((I[d] - 1) * g.dx[d] for d = 1:D)
                es.E[1][I] = 0.001 * sin(phase)
                es.E[2][I] = -es.E[1][I]
            end
        else
            for p = 1:nparticles(e)
                phase = sum(e.x[d][p] for d = 1:D)
                e.x[1][p] -= 0.01 * sin(phase)
            end
        end
        init_empic!(es, e, ions)
        for (count, dt, ns) in
            ((10refinement, 0.04 / refinement, 2), (20refinement, 0.02 / refinement, 3))
            es.n_sub = ns
            for _ = 1:count
                step_empic!(es, e, ions, dt)
            end
        end
        return abs(mode_amplitude(es.E[1], g, ntuple(_ -> 1, D)))
    end
    for D in (2, 3), mi in (16.0, 100.0), transverse in (false, true)
        values = [subcycling_amplitude(D, mi, transverse, r) for r in (1, 2, 4)]
        @test log2(abs(values[1] - values[2]) / abs(values[2] - values[3])) > 1.7
    end
end

@testset "EMPIC multidimensional subcycling matches constant-field particle oracles" begin
    for D in (2, 3), relativistic in (false, true)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0, D))
        e = ParticleSet{D,Float64}(1; q = -1.0, m = 1.0)
        ions = ParticleSet{D,Float64}(1; q = 1.0, m = 16.0)
        for (s, x0) in ((e, 1.0), (ions, 2.0))
            s.weight .= 0.0
            for d = 1:D
                s.x[d] .= x0
            end
        end
        es = EMPIC(g, 1; n0 = 0.0, mobile = true, mi = 16.0, c = 2.0, n_sub = 2, relativistic)
        init_empic!(es, e, ions)
        force = (0.7, -0.2, 0.1)
        for c = 1:3
            es.E[c] .= -force[c]
        end
        elapsed = 0.0
        for (dt, ns) in ((0.1, 2), (0.2, 3), (0.05, 4))
            es.n_sub = ns
            step_empic!(es, e, ions, dt)
            elapsed += dt
            for (s, x0, half_dt) in ((e, 1.0, dt / (2ns)), (ions, 2.0, dt / 2))
                acceleration = ntuple(c -> -s.q / s.m * force[c], 3)
                u = ntuple(c -> acceleration[c] * (elapsed - half_dt), 3)
                gamma = relativistic ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0
                @test [s.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 5e-14
                if !relativistic
                    for d = 1:D
                        @test s.x[d][1] ≈ x0 + acceleration[d] * elapsed^2 / 2 atol = 5e-14
                    end
                end
            end
            @test es.last_sub_dt[] == dt / ns
            @test es.last_dt[] == dt
        end
    end
end

@testset "Float32 EM subcycling charge and phase constraints" begin
    for D = 1:3, mobile in (false, true), rel in (false, true), shape in (CIC(), TSC())
        g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> Float32(2π), D))
        e = ParticleSet{D,Float32}(8^D; q = -1.0f0, m = 1.0f0)
        ions = ParticleSet{D,Float32}(8^D; q = 1.0f0, m = 16.0f0)
        for s in (e, ions)
            load_lattice!(s, ntuple(_ -> 0.0f0, D), g.L, g.n)
            set_density_weight!(s, 1.0f0, g)
        end
        @. e.x[1] -= 0.02f0 * sin(e.x[1])
        e.v[1] .= 0.1f0
        @. e.v[2] = 0.02f0 * cos(e.x[1])
        es = EMPIC(
            g,
            nparticles(e);
            mobile,
            mi = 16.0f0,
            c = 2.0f0,
            relativistic = rel,
            shape,
            n_sub = 2,
        )
        init_empic!(es, e, mobile ? ions : nothing)
        for (dt, ns) in ((0.01, 2), (0.02, 3), (0.015, 1))
            es.n_sub = ns
            step_empic!(es, e, mobile ? ions : nothing, dt)
            @test es.last_dt[] === Float32(dt)
            @test es.last_sub_dt[] === Float32(dt) / ns
            @test _empic_gauss_residual(es, es.rho_np1) < 128eps(Float32)
            ds = maximum(abs, es.rho_np1 .- es.rho_n)
            @test charge_conservation_residual(es, Float32(dt) / ns) * ds < 128eps(Float32)
            @test all(all(isfinite, v) for v in (e.x..., e.v..., es.E..., es.B...))
            if rel
                @test maximum(sum(e.v[c] .^ 2 for c = 1:3)) < es.c^2
            end
        end
    end
end
@testset "Float32 EM subcycling constant-field momentum" begin
    for D = 1:3, mobile in (false, true), rel in (false, true)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0f0, D))
        e = ParticleSet{D,Float32}(1; q = -1.0f0, m = 1.0f0)
        ions = ParticleSet{D,Float32}(1; q = 1.0f0, m = 16.0f0)
        for s in (e, ions)
            s.weight .= 0.0f0
            for d = 1:D
                s.x[d] .= 1.0f0
            end
        end
        es = EMPIC(g, 1; n0 = 0.0f0, mobile, mi = 16.0f0, c = 2.0f0, relativistic = rel, n_sub = 2)
        init_empic!(es, e, mobile ? ions : nothing)
        force = (0.7f0, -0.2f0, 0.1f0)
        for c = 1:3
            es.E[c] .= -force[c]
        end
        elapsed = 0.0f0
        for (dt, ns) in ((0.1f0, 2), (0.2f0, 3), (0.05f0, 4))
            es.n_sub = ns
            step_empic!(es, e, mobile ? ions : nothing, dt)
            elapsed += dt
            for (s, h) in ((e, dt / (2ns)), (ions, dt / 2))
                s === ions && !mobile && continue
                u = ntuple(c -> -s.q / s.m * force[c] * (elapsed - h), 3)
                gamma = rel ? sqrt(1 + sum(abs2, u) / es.c^2) : 1.0f0
                @test [s.v[c][1] for c = 1:3] ≈ collect(u) ./ gamma atol = 32eps(Float32)
            end
        end
    end
end


# Four finite hardware-float factors and their sum fit exactly at this precision.
# This oracle never materializes a Rational through Base's hardware conversion.
function _empic_exact_maxwell_value(::Type{T}, E, curlB, J, dt, c) where {T}
    widest_exponent = max(exponent(floatmax(T)), -exponent(nextfloat(zero(T))))
    return setprecision(BigFloat, 8widest_exponent + 4precision(T) + 32) do
        setrounding(BigFloat, RoundNearest) do
            value = BigFloat(E) + BigFloat(dt) * (BigFloat(c)^2 * BigFloat(curlB) - BigFloat(J))
            T(value, RoundNearest)
        end
    end
end

@testset "EM Maxwell complete-expression range boundaries" begin
    for T in (Float32, Float64)
        H = HybridPlasmaPIC
        tiny = nextfloat(zero(T))
        maxval = floatmax(T)
        large_power = T === Float32 ? 80 : 600
        residual_power = T === Float32 ? -120 : -1000
        c_large = ldexp(one(T), large_power)
        dt_small = inv(c_large)
        cases = (
            # E cancels the magnetic update before the small current survives.
            (-one(T), inv(c_large), ldexp(one(T), residual_power + large_power), dt_small, c_large),
            (ldexp(one(T), residual_power), c_large, inv(c_large), c_large, inv(c_large)),
            # The force or its timestep product can overflow before cancellation.
            (zero(T), maxval, -maxval, T(0.25), one(T)),
            (-maxval, maxval / 2, -maxval / 2, T(2), one(T)),
            # Rounding the increment before adding E loses a subnormal tie.
            (tiny, floatmin(T), prevfloat(floatmin(T)), T(0.5), one(T)),
            # A truly unrepresentable final field stays nonfinite.
            (maxval, maxval, zero(T), one(T), one(T)),
        )
        for (E, curlB, J, dt, c) in cases
            c2 = c * c
            magnetic = c2 * curlB
            forcing = magnetic - J
            increment = dt * forcing
            ordinary = E + increment
            @test H._em_update_lost_range(c2, curlB, magnetic, forcing, increment, ordinary)
            coeff = H._em_exact_update_coefficients(c, dt)
            @test H._em_exact_electric_update(E, curlB, J, coeff) ===
                  _empic_exact_maxwell_value(T, E, curlB, J, dt, c)
        end
        # Neither unchanged zeros nor normal products need exact arithmetic.
        @test !H._em_update_lost_range(one(T), zero(T), zero(T), zero(T), zero(T), zero(T))
        @test !H._em_update_lost_range(T(4), T(0.5), T(2), one(T), T(0.25), T(0.75))
        @test _empic_exact_maxwell_value(T, tiny, zero(T), -tiny, T(0.5), one(T)) == 2tiny
    end
end

@testset "EM Maxwell waves across floating-point ranges" begin
    for T in (Float32, Float64)
        cases =
            T === Float32 ?
            (
                (T(1e30), T(1e-30), T(1e-30)),  # c² overflows
                (T(1e-30), T(1e30), T(1e30)),   # c² underflows
                (T(1e-20), T(1e20), T(1e20)),   # subnormal c² loses precision
                (T(1e10), T(1e-10), T(1e20)),   # c² curlB overflows
                (T(1e-10), T(1e10), T(1e-27)),  # c² curlB underflows
                (T(1e-10), T(1e-30), T(1e30)),  # forming dt*c² first would underflow
                (T(3), T(0.05), T(0.4)),
            ) :
            (
                (T(1e200), T(1e-200), T(1e-200)),
                (T(1e-200), T(1e200), T(1e200)),
                (T(1e-155), T(1e155), T(1e155)),
                (T(1e100), T(1e-100), T(1e200)),
                (T(1e-100), T(1e100), T(1e-200)),
                (T(1e-100), T(1e-200), T(1e200)),
                (T(3), T(0.05), T(0.4)),
            )
        for (State, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
            axis = 1:D,
            (c, dt, amplitude) in cases,
            sign in (-one(T), one(T))

            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
            e = ParticleSet{D,T}(0; q = -one(T), m = one(T))
            es = State(g, 0; n0 = zero(T), c)
            electric_axis, magnetic_axis = mod1(axis + 1, 3), mod1(axis + 2, 3)
            B = es isa EMPIC1D ? es.Bz : es.B[magnetic_axis]
            E = es isa EMPIC1D ? es.Ey : es.E[electric_axis]
            sine, cosine = (zero(T), one(T), zero(T), -one(T)), (one(T), zero(T), -one(T), zero(T))
            for I in CartesianIndices(B)
                B[I] = sign * amplitude * sine[I[axis]]
            end
            B0 = copy(B)
            init_empic!(es, e)
            step_empic!(es, e, dt)
            peak = _empic_exact_maxwell_value(
                T,
                zero(T),
                -sign * amplitude * g.kvec[axis][2],
                zero(T),
                dt,
                c,
            )
            @test isfinite(peak) && !iszero(peak)
            @test E ≈ [peak * cosine[I[axis]] for I in CartesianIndices(E)] rtol = 16eps(T) atol = zero(T)
            @test B == B0
            fields = es isa EMPIC1D ? (es.Ex, es.Ey, es.Bz) : (es.E..., es.B...)
            @test all(a -> all(isfinite, a), fields)
            @test es.rho_n == es.rho_np1 == zero(es.rho_n)
            @test es.time[] == dt
        end
        # The zero-field vacuum is also a solution when c² itself is not representable.
        for (State, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
            c in (cases[1][1], cases[2][1], cases[3][1])

            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
            e = ParticleSet{D,T}(0; q = -one(T), m = one(T))
            es = State(g, 0; n0 = zero(T), c)
            init_empic!(es, e)
            step_empic!(es, e, inv(c))
            fields = es isa EMPIC1D ? (es.Ex, es.Ey, es.Bz) : (es.E..., es.B...)
            @test all(a -> all(iszero, a), fields)
        end
    end
end

@testset "EM subnormal longitudinal field update" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
        tiny = nextfloat(zero(T))
        g = FourierGrid((1,), (one(T),))
        e = ParticleSet{1,T}(1; q = -one(T), m = one(T))
        e.weight[1], e.v[1][1] = 2tiny, T(0.5)
        es = State(g, 1; n0 = 2tiny, c = one(T))
        init_empic!(es, e)
        E = es isa EMPIC1D ? es.Ex : es.E[1]
        E[1] = tiny
        step_empic!(es, e, T(0.5))
        J = es isa EMPIC1D ? es.Jx : es.J[1]
        @test J[1] == -tiny
        @test es.rho_n == es.rho_np1 == T[0]
        @test E[1] == 2tiny
    end
end

@testset "EM extreme-c Maxwell phase and subcycling" begin
    for T in (Float32, Float64),
        (State, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
        large in (false, true)

        power = T === Float32 ? 80 : 600
        c = ldexp(one(T), large ? power : -power)
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2pi), D))
        e = ParticleSet{D,T}(0; q = -one(T), m = one(T))
        es = State(g, 0; n0 = zero(T), c)
        B = es isa EMPIC1D ? es.Bz : es.B[3]
        E = es isa EMPIC1D ? es.Ey : es.E[2]
        sine, cosine = (zero(T), one(T), zero(T), -one(T)), (one(T), zero(T), -one(T), zero(T))
        for I in CartesianIndices(B)
            B[I] = sine[I[1]] / c
        end
        init_empic!(es, e)
        e_amplitude, b_amplitude, previous_dt = BigFloat(0), BigFloat(1), BigFloat(0)
        for (alpha, ns) in ((T(0.125), 1), (T(0.25), 2), (T(0.0625), 1), (T(0.375), 3))
            es.n_sub = ns
            dt = alpha / c
            sub_dt = dt / ns
            for substep = 1:ns
                iszero(previous_dt) && (previous_dt = BigFloat(sub_dt))
                kick = previous_dt / 2 + BigFloat(sub_dt) / 2
                b_amplitude += BigFloat(c) * kick * BigFloat(g.kvec[1][2]) * e_amplitude
                e_amplitude -= BigFloat(sub_dt) * BigFloat(c) * BigFloat(g.kvec[1][2]) * b_amplitude
                previous_dt = BigFloat(sub_dt)
            end
            step_empic!(es, e, dt)
            @test E ≈ [T(e_amplitude) * cosine[I[1]] for I in CartesianIndices(E)] rtol = 128eps(T) atol = zero(T)
            @test c .* B ≈ [T(b_amplitude) * sine[I[1]] for I in CartesianIndices(B)] rtol = 128eps(T) atol = zero(T)
            @test all(isfinite, E) && all(isfinite, B)
            @test es.last_sub_dt[] == sub_dt
        end
    end
end

@testset "EM Maxwell squared-speed transition neighbors" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC), n in (4, 5)
        roots = (sqrt(nextfloat(zero(T))), sqrt(floatmin(T)), sqrt(floatmax(T)))
        for root in roots, c in (prevfloat(root), root, nextfloat(root))
            g = FourierGrid((n,), (T(2pi),))
            e = ParticleSet{1,T}(0; q = -one(T), m = one(T))
            es = State(g, 0; n0 = zero(T), c)
            B = es isa EMPIC1D ? es.Bz : es.B[3]
            E = es isa EMPIC1D ? es.Ey : es.E[2]
            dt, amplitude = inv(c), inv(c)
            phases = [T(2pi) * T(i - 1) / T(n) for i = 1:n]
            B .= amplitude .* sin.(phases)
            init_empic!(es, e)
            step_empic!(es, e, dt)
            peak = _empic_exact_maxwell_value(T, zero(T), -amplitude * g.kvec[1][2], zero(T), dt, c)
            @test E ≈ peak .* cos.(phases) rtol = 128eps(T) atol = zero(T)
            @test all(isfinite, E) && all(isfinite, B)
        end
    end
end

@testset "EM range repair preserves nonfinite field visibility" begin
    for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
        c = T === Float32 ? T(1e30) : T(1e200)
        g = FourierGrid((4,), (T(2pi),))
        e = ParticleSet{1,T}(0; q = -one(T), m = one(T))
        es = State(g, 0; n0 = zero(T), c)
        E = es isa EMPIC1D ? es.Ey : es.E[2]
        E[1] = T(NaN)
        init_empic!(es, e)
        # The spectral input boundary owns this rejection, before field advancement.
        @test_throws ArgumentError step_empic!(es, e, inv(c))
        @test es.step[] == 0 && iszero(es.time[])
        fields = es isa EMPIC1D ? (es.Ex, es.Ey, es.Bz) : (es.E..., es.B...)
        @test any(a -> any(!isfinite, a), fields)
    end
end

# Independent real-space discrete Fourier matrix: no FFT or production derivative
# helper is used by this oracle. Its k values include the grid's Nyquist omission.
function _empic_dft_matrices(g)
    return ntuple(length(g.n)) do axis
        n = g.n[axis]
        matrix = Matrix{BigFloat}(undef, n, n)
        for i = 1:n, j = 1:n
            matrix[i, j] =
                -sum(
                    BigFloat(g.kvec[axis][m]) * sinpi(BigFloat(2g.midx[axis][m] * (i - j)) / n) for
                    m = 1:n
                ) / n
        end
        matrix
    end
end

function _empic_dft_derivative(field, g, axis, matrices = _empic_dft_matrices(g))
    n = g.n[axis]
    matrix = matrices[axis]
    result = zeros(BigFloat, size(field))
    for I in CartesianIndices(field), j = 1:n
        J = CartesianIndex(ntuple(d -> d == axis ? j : I[d], ndims(field)))
        result[I] += matrix[I[axis], j] * BigFloat(field[J])
    end
    return result
end

function _empic_dft_curl(fields, g, matrices = _empic_dft_matrices(g))
    result = ntuple(_ -> zeros(BigFloat, g.n), 3)
    D = length(g.n)
    for c = 1:3
        plus, minus = mod1(c + 1, 3), mod1(c + 2, 3)
        plus <= D && (result[c] .+= _empic_dft_derivative(fields[minus], g, plus, matrices))
        minus <= D && (result[c] .-= _empic_dft_derivative(fields[plus], g, minus, matrices))
    end
    return result
end

function _empic_dft_step!(E, B, g, dt, ns, previous, first, matrices = _empic_dft_matrices(g))
    T = typeof(dt)
    fast = dt / ns
    if first
        derivative = _empic_dft_curl(E, g, matrices)
        for c = 1:3
            B[c] .= T.(BigFloat.(B[c]) .+ BigFloat(fast) / 2 .* derivative[c])
        end
        previous = fast
    end
    centered = map(copy, B)
    for _ = 1:ns
        derivative = _empic_dft_curl(E, g, matrices)
        for c = 1:3
            centered[c] .= T.(BigFloat.(B[c]) .- BigFloat(previous) / 2 .* derivative[c])
            B[c] .=
                T.(BigFloat.(B[c]) .- (BigFloat(previous) + BigFloat(fast)) / 2 .* derivative[c])
        end
        derivative = _empic_dft_curl(B, g, matrices)
        for c = 1:3
            E[c] .= T.(BigFloat.(E[c]) .+ BigFloat(fast) .* derivative[c])
        end
        previous = fast
    end
    return centered, fast
end

@testset "EM scaled spectral updates follow independent DFT evolution" begin
    setprecision(BigFloat, 4096) do
        for T in (Float32, Float64),
            (State, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
            axis = 1:D,
            large in (false, true),
            seed in (:electric, :magnetic),
            sign in (-1, 1)

            scale, amplitude = if T === Float32
                large ? (T(1e-10), T(1e30)) : (T(1e20), T(1e-30))
            else
                large ? (T(1e-100), T(1e300)) : (T(1e100), T(1e-300))
            end
            g = FourierGrid(
                ntuple(_ -> 4, D),
                ntuple(d -> T(2pi) * (d == axis ? scale : T(d + 1)), D),
            )
            e = ParticleSet{D,T}(0; q = -one(T), m = one(T))
            es = State(g, 0; n0 = zero(T), c = one(T))
            E = es isa EMPIC1D ? (es.Ex, es.Ey, zeros(T, g.n)) : es.E
            B = es isa EMPIC1D ? (zeros(T, g.n), zeros(T, g.n), es.Bz) : es.B
            component = seed === :electric ? mod1(axis + 1, 3) : mod1(axis + 2, 3)
            field = seed === :electric ? E[component] : B[component]
            for I in CartesianIndices(field)
                field[I] = T(sign) * amplitude * (zero(T), one(T), zero(T), -one(T))[I[axis]]
            end
            init_empic!(es, e)
            expected_E, expected_B = map(copy, E), map(copy, B)
            matrices = _empic_dft_matrices(g)
            previous = zero(T)
            for (iteration, (dt, ns)) in enumerate(((scale, 1), (scale / 2, 3), (scale, 2)))
                es.n_sub = ns
                centered, previous = _empic_dft_step!(
                    expected_E,
                    expected_B,
                    g,
                    dt,
                    ns,
                    previous,
                    iteration == 1,
                    matrices,
                )
                step_empic!(es, e, dt)
                actual_centered = es isa EMPIC1D ? (zeros(T, g.n), zeros(T, g.n), es.Bzc) : es.Bc
                for (actual, expected) in zip(
                    (E..., B..., actual_centered...),
                    (expected_E..., expected_B..., centered...),
                )
                    @test all(isfinite, actual)
                    @test actual ≈ expected rtol = 128eps(T) atol = 128eps(T) * amplitude
                end
                @test es.last_sub_dt[] == dt / ns
            end
            before = deepcopy((E, B, es.last_dt[], es.last_sub_dt[], es.time[], es.step[]))
            step_empic!(es, e, zero(T))
            @test (E, B, es.last_dt[], es.last_sub_dt[], es.time[], es.step[]) == before
            init_empic!(es, e)
            @test es.step[] == 0 && iszero(es.last_sub_dt[])
        end
    end
end

@testset "EM common curl scaling retains longitudinal and diagonal components" begin
    setprecision(BigFloat, 4096) do
        for T in (Float32, Float64),
            D = 1:3,
            axis = 1:D,
            seed in (:electric, :magnetic),
            sign in (-1, 1)
            # The large diagonal component changes no physical curl, but does
            # control the canonical curl's common FFT scale. Its transverse
            # companion must retain the bits later restored by the timestep.
            large_exp, k_exp, small_exp = T === Float32 ? (64, 100, -100) : (700, 700, -670)
            scale = ldexp(one(T), -k_exp)
            large, small = ldexp(one(T), large_exp), ldexp(T(1.3), small_exp)
            g = FourierGrid(
                ntuple(_ -> 4, D),
                ntuple(d -> T(2pi) * (d == axis ? scale : one(T)), D),
            )
            e = ParticleSet{D,T}(0; q = -one(T), m = one(T))
            es = EMPIC(g, 0; n0 = zero(T), c = one(T))
            init_empic!(es, e)
            fields = seed === :electric ? es.E : es.B
            fields[axis] .= large
            component = seed === :electric ? mod1(axis + 1, 3) : mod1(axis + 2, 3)
            for I in CartesianIndices(fields[component])
                fields[component][I] =
                    T(sign) * small * (zero(T), one(T), zero(T), -one(T))[I[axis]]
            end
            expected_E, expected_B = map(copy, es.E), map(copy, es.B)
            matrices = _empic_dft_matrices(g)
            dt = scale / 2
            centered, _ =
                _empic_dft_step!(expected_E, expected_B, g, dt, 1, zero(T), true, matrices)
            step_empic!(es, e, dt)
            for (actual, expected) in
                zip((es.E..., es.B..., es.Bc...), (expected_E..., expected_B..., centered...))
                @test all(isfinite, actual)
                @test actual ≈ expected rtol = 128eps(T) atol = 128eps(T) * small
            end
        end
    end
end

@testset "EM magnetic half intervals retain original timestep operands" begin
    setprecision(BigFloat, 4096) do
        for T in (Float32, Float64), State in (EMPIC1D, EMPIC)
            scale = ldexp(one(T), T === Float32 ? -100 : -700)
            # Subnormal timesteps (or subnormal half steps) are rejected before
            # any state changes; the smallest resolved timestep keeps the
            # original-operand half intervals.
            dt = 2floatmin(T)
            g = FourierGrid((4,), (T(2pi) * scale,))
            e = ParticleSet{1,T}(0; q = -one(T), m = one(T))
            es = State(g, 0; n0 = zero(T), c = one(T))
            E = es isa EMPIC1D ? (es.Ex, es.Ey, zeros(T, g.n)) : es.E
            B = es isa EMPIC1D ? (zeros(T, g.n), zeros(T, g.n), es.Bz) : es.B
            E[2] .= T[0, 1, 0, -1]
            init_empic!(es, e)
            before = deepcopy((E, B, es.time[], es.step[], es.last_dt[], es.last_sub_dt[]))
            for unresolved in (nextfloat(zero(T)), floatmin(T), prevfloat(dt))
                @test_throws ArgumentError step_empic!(es, e, unresolved)
                @test (E, B, es.time[], es.step[], es.last_dt[], es.last_sub_dt[]) == before
            end
            expected_E, expected_B = map(copy, E), map(copy, B)
            matrices = _empic_dft_matrices(g)
            previous = zero(T)
            for (i, next_dt) in enumerate((dt, 3dt, 2dt))
                _, previous = _empic_dft_step!(
                    expected_E,
                    expected_B,
                    g,
                    next_dt,
                    1,
                    previous,
                    i == 1,
                    matrices,
                )
                step_empic!(es, e, next_dt)
                @test B[3] ≈ expected_B[3] rtol = 16eps(T) atol = zero(T)
                @test !iszero(B[3][1])
            end
        end
        H = HybridPlasmaPIC
        for T in (Float32, Float64), backward in (false, true)
            largest, tiny = floatmax(T), nextfloat(zero(T))
            cases = (
                (-largest, largest, T(2), T(2)),
                (tiny, tiny, one(T), zero(T)),
                (largest, largest, T(2), T(2)),
            )
            for (old, derivative, before, after) in cases
                result = T[0]
                H._em_magnetic_update!(
                    result,
                    T[old],
                    T[derivative],
                    nothing,
                    before,
                    after,
                    backward,
                )
                sign = backward ? 1 : -1
                expected = T(
                    BigFloat(old) +
                    sign * (BigFloat(before) + BigFloat(after)) / 2 * BigFloat(derivative),
                    RoundNearest,
                )
                @test result[1] === expected
            end
        end
    end
end

@testset "EM scaled derivative bands preserve all signs and spectral errors" begin
    H = HybridPlasmaPIC
    setprecision(BigFloat, 4096) do
        # Odd/even transforms with mixed magnitudes and signs also exercise
        # ordinary FFT rounding, beyond the exactly represented four-point mode.
        for T in (Float32, Float64), D = 1:2, n in (5, 6), axis = 1:D
            g = FourierGrid(ntuple(_ -> n, D), ntuple(d -> T(2pi) * T(d), D))
            high, low = floatmax(T) / 128, 8floatmin(T)
            values = (high, low, -one(T), -high, 2low, one(T))
            field = reshape(T[values[mod1(i, 6)] for i = 1:prod(g.n)], g.n)
            out = similar(field)
            parts = H._em_derivative_parts!(out, field, g, axis)
            @test parts !== nothing
            lo, hi = H._em_field_exponents(field)
            klo, khi = H._em_field_exponents(g.kvec[axis])
            lower, upper, _ = H._em_scale_bounds(T, 0, 0, klo, khi, length(field))
            @test length(parts) <= cld(hi - lo + 1, upper - lower + 1)
            reference = T.(_empic_dft_derivative(field, g, axis))
            actual = [H._from_exact(T, H._em_exact_derivative(parts, i)) for i in eachindex(field)]
            @test reshape(actual, g.n) ≈ reference rtol = 64eps(T) atol = 64eps(T) * high
        end
        for T in (Float32, Float64), D = 1:3, axis = 1:D, sign in (-1, 1), cancel in (false, true)
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(d -> T(2pi) * T(d), D))
            high, low = floatmax(T), floatmin(T)
            values = T[high, low, cancel ? high : -high, -low]
            field = [T(sign) * values[I[axis]] for I in CartesianIndices(g.n)]
            out = similar(field)
            parts = H._em_derivative_parts!(out, field, g, axis)
            @test parts !== nothing
            @test length(parts) <= (T === Float32 ? 5 : 3)
            expected = _empic_dft_derivative(field, g, axis)
            for i in eachindex(field)
                actual = H._from_exact(T, H._em_exact_derivative(parts, i))
                reference = T(expected[i], RoundNearest)
                @test actual ≈ reference rtol = 64eps(T) atol = zero(T)
            end
        end
        for T in (Float32, Float64), residual in (false, true), sign in (-1, 1)
            a, b, p = T === Float32 ? (80, 60, 80) : (700, 600, 600)
            g = FourierGrid((4, 4), (T(2pi) * ldexp(one(T), -a), T(2pi) * ldexp(one(T), -b)))
            ax, ay = ldexp(one(T), p), ldexp(one(T), p + b - a)
            residual && (ay = nextfloat(ay))
            pattern = [
                (zero(T), one(T), zero(T), -one(T))[mod1(I[1] + I[2] - 1, 4)] for
                I in CartesianIndices(g.n)
            ]
            fields = (T(sign) * ax .* pattern, T(sign) * ay .* pattern, zeros(T, g.n))
            out = map(similar, fields)
            parts = H._em_curl_parts!(out, fields, g)
            @test parts !== nothing
            dt = ldexp(one(T), -a)
            expected = _empic_dft_curl(fields, g)
            for c = 1:3, i in eachindex(pattern)
                actual =
                    H._from_exact(T, H._exact_rational(dt) * H._em_exact_derivative(parts[c], i))
                reference = T(BigFloat(dt) * expected[c][i], RoundNearest)
                @test isfinite(actual)
                @test actual ≈ reference rtol = 64eps(T) atol = zero(T)
            end
        end
    end
    for T in (Float32, Float64), n in (4, 6)
        g = FourierGrid((n,), (T(2pi),))
        field = [isodd(i) ? nextfloat(zero(T)) : floatmax(T) for i = 1:n]
        out = similar(field)
        @test_throws ArgumentError deriv!(out, field, g, 1)
        @test_throws ArgumentError H._em_derivative_parts!(out, field, g, 1)
        fields = (zeros(T, n), field, zeros(T, n))
        outputs = map(similar, fields)
        @test_throws ArgumentError curl!(outputs, fields, g)
        @test_throws ArgumentError H._em_curl_parts!(outputs, fields, g)
        field[1] = T(NaN)
        @test_throws ArgumentError H._em_derivative_parts!(out, field, g, 1)
        @test_throws ArgumentError H._em_curl_parts!(outputs, fields, g)
    end
end
