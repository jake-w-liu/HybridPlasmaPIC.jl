using HybridPlasmaPIC, Test

# Independent positive-rate formula. Derive the RK4 edges from exact polygon
# vertices, then use the same binary64 source constants as the API contract.
function _bigfloat_timestep_reference(
    g;
    NB = 1,
    integrator = :rk4,
    safety = 0.8,
    d_i = 1.0,
    Ω_ci = 1.0,
    η = 0.0,
    ηH = 0.0,
)
    k = sqrt(sum(dx -> (BigFloat(π) / BigFloat(dx))^2, g.dx))
    K = k * BigFloat(d_i)
    omega = BigFloat(Ω_ci) * (K * sqrt(K^2 + 4) + K^2) / 2
    decay = BigFloat(η) * k^2 + BigFloat(ηH) * k^4
    s = BigFloat(safety)
    gyro = s * BigFloat(0.3) / BigFloat(Ω_ci)
    integrator === :cn && return gyro
    cl = integrator in (:camcl, :leapfrog)
    wave_constant = cl ? (decay > 0 ? 0.3 : 2.0) : 2.8
    wave = s * BigFloat(NB) * BigFloat(wave_constant) / omega
    iszero(decay) && return min(gyro, wave)
    if cl
        return min(gyro, wave, s * BigFloat(NB) * BigFloat(0.25) / decay)
    end
    vertices = (
        (0 // 1, 14 // 5),
        (1 // 1, 5 // 2),
        (37 // 20, 37 // 20),
        (5 // 2, 1 // 1),
        (557 // 200, 0 // 1),
    )
    polygon = minimum(1:4) do edge
        a, b = vertices[edge:edge+1]
        nx, ny = a[2] - b[2], b[1] - a[1]
        intercept = nx * a[1] + ny * a[2]
        s * BigFloat(NB) * BigFloat(Float64(intercept)) /
        (BigFloat(Float64(nx)) * decay + BigFloat(Float64(ny)) * omega)
    end
    return min(gyro, wave, polygon)
end

function _check_bigfloat_timestep(g, options; bits = 53, mode = RoundNearest)
    stored = deepcopy(options)
    before_precision, before_rounding = precision(BigFloat), rounding(BigFloat)
    result = setprecision(BigFloat, bits) do
        setrounding(BigFloat, mode) do
            value = recommended_dt(g; options...)
            @test value isa BigFloat
            @test precision(value) == bits
            @test isfinite(value) && value > 0
            @test precision(BigFloat) == bits
            @test rounding(BigFloat) == mode
            value
        end
    end
    @test precision(BigFloat) == before_precision
    @test rounding(BigFloat) == before_rounding
    @test isequal(options, stored)
    stored_bits = maximum(x -> x isa BigFloat ? precision(x) : 0, values(options))
    setprecision(BigFloat, max(2048, bits + 512, stored_bits + 512)) do
        setrounding(BigFloat, RoundNearest) do
            reference = _bigfloat_timestep_reference(g; options...)
            @test BigFloat(result) <= reference
            # Directed rounding can discard almost one target ulp at low p.
            @test BigFloat(result) >= reference * (1 - BigFloat(2)^(1 - bits))
        end
    end
    return result
end

@testset "BigFloat timestep owns its arithmetic context" begin
    g = FourierGrid((129,), (2π,))
    eta = setprecision(BigFloat, 512) do
        BigFloat(10)^200
    end
    for bits in (2, 11, 24, 53, 256),
        mode in (RoundNearest, RoundDown, RoundUp, RoundToZero, RoundFromZero),
        integrator in (:rk4, :camcl, :leapfrog, :cn)

        _check_bigfloat_timestep(g, (; η = eta, safety = 1.0, integrator); bits, mode)
    end
    for integrator in (:rk4, :camcl, :leapfrog, :cn), role in (:safety, :d_i, :Ω_ci, :η, :ηH)
        value = role === :safety ? 0.75 : 0.125
        for converted in (BigFloat(value), big(3) // 8)
            options = merge((; integrator), NamedTuple{(role,)}((converted,)))
            _check_bigfloat_timestep(g, options; bits = 11, mode = RoundUp)
        end
    end
    for integrator in (:rk4, :camcl, :leapfrog, :cn), role in (:safety, :d_i, :Ω_ci, :η, :ηH)
        options = merge((; integrator), NamedTuple{(role,)}((big(1),)))
        _check_bigfloat_timestep(g, options; bits = 11)
    end
    for integrator in (:rk4, :camcl, :leapfrog)
        _check_bigfloat_timestep(g, (; integrator, NB = big(3)); bits = 2)
        _check_bigfloat_timestep(g, (; integrator, NB = big(2)^500); bits = 24)
    end
    # NB is unused by CN and must not change its historical native result.
    @test recommended_dt(g; integrator = :cn, NB = big(3)) ===
          recommended_dt(g; integrator = :cn, NB = 3)
    for (bits, stored_precision) in ((11, 512), (256, 1024), (512, 24))
        value = setprecision(BigFloat, stored_precision) do
            BigFloat(1) + BigFloat(2)^(-stored_precision + 2)
        end
        _check_bigfloat_timestep(g, (; η = value); bits)
    end
end

@testset "Low-precision recommendation damps both production RK4 owners" begin
    for T in (Float32, Float64), evolving in (false, true), bits in (2, 11, 53)
        g = FourierGrid((129,), (T(2π),))
        eta = setprecision(BigFloat, 512) do
            BigFloat(10)^(T === Float32 ? 20 : 200)
        end
        coefficient = T(eta)
        closure = evolving ? TransportElectrons(T(5 / 3)) : IsothermalElectrons(one(T))
        st = HybridStepper(g, HybridModel(closure; η = coefficient), CIC(), 0)
        fields = st.fields
        fill!(st.fn, one(T))
        fill!(fields.ninv, one(T))
        fill!(fields.pe, one(T))
        fill!(fields.B[1], one(T))
        x = T.(0:128) .* g.dx[1]
        @. fields.B[2] = T(0.001) * cos(T(64) * x)
        @. fields.B[3] = T(0.001) * sin(T(64) * x)
        initial = sum(abs2, fields.B[2]) + sum(abs2, fields.B[3])
        step = setprecision(BigFloat, bits) do
            recommended_dt(g; η = eta, safety = 1.0)
        end
        dt = T(step)
        if evolving
            HybridPlasmaPIC._rk4_Bp!(st, dt, 1)
        else
            HybridPlasmaPIC._rk4_B!(st, dt)
        end
        measured = sqrt((sum(abs2, fields.B[2]) + sum(abs2, fields.B[3])) / initial)
        z = complex(-Float64(coefficient) * 64^2 * Float64(dt), Float64(64^2) * Float64(dt))
        polynomial = abs(1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24)
        @test isapprox(measured, polynomial; rtol = 128eps(T))
        @test measured <= 1 + 128eps(T)
    end
end

@testset "BigFloat timestep retains huge representable exponent limits" begin
    g = FourierGrid((4, 4), (2π, 4π))
    # MPFR's C-long exponent range is narrower on some supported platforms.
    # Bound wide integer scales before converting to Int; no architecture skip.
    range_scale = min(Int64(Base.MPFR.get_emax_max()), -Int64(Base.MPFR.get_emin_min())) ÷ 8
    exponent_scales =
        (Int(min(Int64(1) << 40, range_scale ÷ 2)), Int(min(Int64(1) << 50, range_scale)))
    @test 0 < first(exponent_scales) < last(exponent_scales)
    @test Int64(last(exponent_scales)) * 4 < Int64(Base.MPFR.get_emax_max())
    for exponent_scale in exponent_scales, sign in (-1, 1), role in (:d_i, :Ω_ci, :η, :ηH)
        value = setprecision(BigFloat, 256) do
            ldexp(BigFloat(1), sign * exponent_scale)
        end
        for integrator in (:rk4, :camcl, :leapfrog, :cn)
            options = merge((; integrator), NamedTuple{(role,)}((value,)))
            _check_bigfloat_timestep(g, options; bits = 53)
        end
    end
end

@testset "Timestep endpoint arithmetic encloses its supported domains" begin
    H = HybridPlasmaPIC
    B = H._TimestepBounds
    setprecision(BigFloat, 256) do
        for rational in
            (big(3) // 7, big(-3) // 7, typemin(Int) // 1, -(big(2)^257 + 1) // (big(2)^257 + 3))
            bounds = B(rational)
            @test Rational{BigInt}(bounds.lo) <= rational <= Rational{BigInt}(bounds.hi)
        end
        pi_bounds = B(π)
        setprecision(BigFloat, 1024) do
            @test pi_bounds.lo <= BigFloat(π) <= pi_bounds.hi
        end
        a = B(BigFloat(-2), BigFloat(3))
        b = B(BigFloat(1), BigFloat(4))
        @test (a + b).lo <= -1 && (a + b).hi >= 7
        @test (a - b).lo <= -6 && (a - b).hi >= 2
        @test (a * b).lo <= -8 && (a * b).hi >= 12
        @test (a / b).lo <= -2 && (a / b).hi >= 3
        @test abs2(a).lo == 0 && abs2(a).hi >= 9
        @test exp(B(-Inf)).lo == exp(B(-Inf)).hi == 0
        @test log(B(0)).lo == log(B(0)).hi == -Inf
        @test_throws ArgumentError B(BigFloat(NaN), BigFloat(1))
        @test_throws ArgumentError B(BigFloat(2), BigFloat(1))
        @test_throws ArgumentError a / B(0)
        @test_throws ArgumentError log(B(-1))
        @test_throws ArgumentError log1p(B(-2))
        @test_throws ArgumentError hypot(a, b)
        limit = exp(H._rk4_damped_log_limit(B(0), B(0)))
        result = BigFloat(limit.lo, RoundToZero; precision = 2)
        z = complex(-result, result)
        @test abs(1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24) <= 1
    end
end

@testset "BigFloat timestep restores nested and interleaved task contexts" begin
    g = FourierGrid((4,), (2π,))
    before_precision, before_rounding = precision(BigFloat), rounding(BigFloat)
    eta = BigFloat(0.125)
    tasks = map(((11, RoundUp), (53, RoundDown), (256, RoundToZero))) do configuration
        bits, mode = configuration
        @async setprecision(BigFloat, bits) do
            setrounding(BigFloat, mode) do
                yield()
                value = recommended_dt(g; η = eta)
                caught = try
                    recommended_dt(g; η = eta, safety = 2)
                    nothing
                catch error
                    error
                end
                (;
                    bits,
                    mode,
                    value,
                    caught,
                    context_precision = precision(BigFloat),
                    context_rounding = rounding(BigFloat),
                )
            end
        end
    end
    # Assert task results in the parent testset on both supported runtimes.
    for task in tasks
        result = fetch(task)
        @test precision(result.value) == result.bits
        @test result.context_precision == result.bits
        @test result.context_rounding == result.mode
        @test result.caught isa ArgumentError
        setprecision(BigFloat, 2048) do
            setrounding(BigFloat, RoundNearest) do
                @test result.value <= _bigfloat_timestep_reference(g; η = eta)
            end
        end
    end
    @test precision(BigFloat) == before_precision
    @test rounding(BigFloat) == before_rounding
end

@testset "BigFloat exponent-range checks run in an isolated process" begin
    script = joinpath(@__DIR__, "test_timestep_bigfloat_range.jl")
    project = dirname(Base.active_project())
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script`
    @test success(run(ignorestatus(command)))
end

# This Real intentionally provides value conversion without zero(Type), which
# the historical public API did not require for native parameter conversion.
struct _TimestepValueReal <: Real
    value::Float64
    _TimestepValueReal(value::Float64) = new(value)
end
Base.float(x::_TimestepValueReal) = x.value
Base.Float64(x::_TimestepValueReal) = x.value
Base.isfinite(x::_TimestepValueReal) = isfinite(x.value)
Base.:<(x::Real, y::_TimestepValueReal) = Float64(x) < y.value
Base.:<(x::_TimestepValueReal, y::Real) = x.value < Float64(y)
Base.:<(x::_TimestepValueReal, y::_TimestepValueReal) = x.value < y.value

@testset "Native timestep keeps value-based Real conversion" begin
    g = FourierGrid((4,), (2π,))
    @test_throws MethodError zero(_TimestepValueReal)
    @test recommended_dt(g; d_i = _TimestepValueReal(1.0)) === recommended_dt(g; d_i = 1.0)
end
