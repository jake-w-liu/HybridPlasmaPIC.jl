module TestBoundedMaxwell
using HybridPlasmaPIC, Test, Random, LinearAlgebra

function fill_component!(a, g, field, c, f)
    x = yee_coordinates(g, field, c)
    for I in CartesianIndices(a)
        a[I] = f(ntuple(d -> x[d][I[d]], length(x)))
    end
end
components(f) = (f.E..., f.B...)
snapshot(m) = (map(copy, components(m.fields)), m.time[], m.step[])
unchanged(m, old) =
    all(isequal.(components(m.fields), old[1])) && m.time[] == old[2] && m.step[] == old[3]
function divergence(f, g, electric)
    out = zeros(eltype(f[1]), ntuple(d -> g.n[d] - Int(electric && !g.periodic[d]), length(g.n)))
    yee_divergence!(out, f, g; field = electric ? :electric : :magnetic)
end
function discrete_energy(m, dt)
    ce = map(similar, m.fields.B)
    yee_curl!(ce, m.fields.E, m.g)
    return (
        sum(a -> sum(abs2, a), m.fields.E) + m.c^2 * sum(a -> sum(abs2, a), m.fields.B) -
        (dt * m.c)^2 * sum(a -> sum(abs2, a), ce) / 4
    ) * prod(m.g.dx) / 2
end

@testset "Yee geometry, constraints and rejection" begin
    for T in (Float32, Float64), D = 1:3, mask = 0:(2^D-1)
        n = ntuple(d -> 4 + d, D)
        bc = ntuple(d -> isodd(mask >> (d - 1)) ? :periodic : :pec, D)
        g = YeeGrid(n, ntuple(d -> T(d + 1), D); boundaries = bc)
        f = YeeFields(g)
        for c = 1:3, field in (:electric, :magnetic)
            a = field === :electric ? f.E[c] : f.B[c]
            coords = yee_coordinates(g, field, c)
            @test size(a) == length.(coords)
            for d = 1:D
                nodal = field === :electric ? c != d : c == d
                @test first(coords[d]) == (nodal ? zero(T) : g.dx[d] / 2)
                @test length(coords[d]) == n[d] + (nodal && bc[d] === :pec)
            end
            fill!(a, T(2))
        end
        apply_yee_pec!(f, g)
        for c = 1:3, field in (:electric, :magnetic)
            a = field === :electric ? f.E[c] : f.B[c]
            for I in CartesianIndices(a)
                wall = any(
                    d ->
                        bc[d] === :pec &&
                            (field === :electric ? c != d : c == d) &&
                            (I[d] == 1 || I[d] == n[d] + 1),
                    1:D,
                )
                @test a[I] == (wall ? 0 : 2)
            end
        end
        m = MaxwellVacuum(g; c = T(1.4))
        @test maxwell_cfl(m) ≈ inv(m.c * sqrt(sum(abs2, inv.(g.dx))))
        old = snapshot(m)
        for bad in (T(0), T(-1), T(NaN), T(Inf), maxwell_cfl(m), nextfloat(maxwell_cfl(m)))
            @test_throws ArgumentError step_maxwell!(m, bad)
            @test unchanged(m, old)
        end
        step_maxwell!(m, T(0.25) * maxwell_cfl(m))
        @test em_field_energy(m) === zero(T)
        @test m.step[] == 1
        @test all(a -> all(iszero, a), yee_poynting(m))
    end
    @test_throws ArgumentError YeeGrid((1,), (1.0,))
    @test_throws ArgumentError YeeGrid((typemax(Int),), (1.0,))
    @test_throws ArgumentError YeeGrid((4,), (0.0,))
    @test_throws ArgumentError YeeGrid((4,), (Inf,))
    @test_throws ArgumentError YeeGrid((4,), (nextfloat(0.0),))
    @test_throws ArgumentError YeeGrid((4,), (BigFloat(1),))
    @test_throws ArgumentError YeeGrid((4,), (1.0,); boundaries = (:open,))
    @test_throws DimensionMismatch YeeGrid((4,), (1.0,); boundaries = (:pec, :pec))
    @test_throws ArgumentError YeeGrid((4, 4, 4, 4), (1.0, 1.0, 1.0, 1.0))
    g = YeeGrid((8,), (1.0,))
    @test_throws ArgumentError yee_coordinates(g, :invalid, 1)
    @test_throws ArgumentError yee_coordinates(g, :electric, 0)
    @test_throws ArgumentError MaxwellVacuum(g; c = 0)
    @test MaxwellVacuum(g; c = floatmax(Float64)).c === floatmax(Float64)
    @test MaxwellVacuum(g; c = nextfloat(0.0)).c === nextfloat(0.0)
    m = MaxwellVacuum(g)
    for (component, index, value) in ((2, 1, 1.0), (2, 4, NaN), (1, 2, Inf))
        m.fields.E[component][index] = value
        old = snapshot(m)
        @test_throws ArgumentError step_maxwell!(m, 0.01)
        @test unchanged(m, old)
        m.fields.E[component][index] = 0
    end
    for a in m.fields.B
        fill!(a, floatmax(Float64))
        a[2:2:end] .*= -1
    end
    apply_yee_pec!(m.fields, g)
    old = snapshot(m)
    @test_throws ArgumentError step_maxwell!(m, 0.1)
    @test unchanged(m, old)
    foreach(a -> fill!(a, 0), m.fields.B)
    step_maxwell!(m, 0.01) # scratch left by the failure is safely overwritten
    @test em_field_energy(m) == 0
    for bad in (NaN, Inf, -1.0, floatmax(Float64))
        m.time[] = bad
        old = snapshot(m)
        @test_throws ArgumentError step_maxwell!(m, 0.01)
        @test all(isequal.(components(m.fields), old[1]))
        @test isequal(m.time[], old[2]) && m.step[] == old[3]
    end
    m.time[] = 0
    m.step[] = typemax(Int)
    old = snapshot(m)
    @test_throws ArgumentError step_maxwell!(m, 0.01)
    @test unchanged(m, old)
    m.step[] = 0
    resize!(m.fields.E[1], 7)
    old = snapshot(m)
    @test_throws DimensionMismatch step_maxwell!(m, 0.01)
    @test unchanged(m, old)

    gp = YeeGrid((6,), (1.0,); boundaries = (:periodic,))
    mp = MaxwellVacuum(gp)
    @test_throws ArgumentError yee_curl!(mp.fields.E, mp.fields.E, gp)
    @test_throws ArgumentError yee_divergence!(mp.fields.E[1], mp.fields.E, gp; field = :electric)
    @test_throws DimensionMismatch yee_divergence!(zeros(5), mp.fields.B, gp)
    @test_throws DimensionMismatch yee_curl!((zeros(5), zeros(5), zeros(5)), mp.fields.E, gp)
    alias_fields = YeeFields((mp.fields.E[1], mp.fields.E[1], mp.fields.E[1]), mp.fields.B)
    alias_model = MaxwellVacuum(gp, alias_fields, mp.c, mp.time, mp.step, mp._candidate, mp._cfl)
    @test_throws ArgumentError step_maxwell!(alias_model, 0.01)
    @test_throws ArgumentError apply_yee_pec!(alias_fields, gp)
    integers = ntuple(_ -> fill(3, 6), 3)
    @test_throws ArgumentError yee_curl!(integers, mp.fields.E, gp)
    @test all(a -> all(==(3), a), integers)
    @test_throws ArgumentError yee_divergence!(integers[1], mp.fields.E, gp; field = :electric)
    @test all(==(3), integers[1])
end

@testset "Compatible curls, divergences, adjoints and propagation" begin
    rng = MersenneTwister(902001)
    for T in (Float32, Float64), D = 1:3, mask = 0:(2^D-1)
        g = YeeGrid(
            ntuple(d -> 5 + d, D),
            ntuple(d -> T(1.3 + 0.2d), D);
            boundaries = ntuple(d -> isodd(mask >> (d - 1)) ? :periodic : :pec, D),
        )
        m = MaxwellVacuum(g; c = T(1.2))
        foreach(a -> randn!(rng, a), components(m.fields))
        apply_yee_pec!(m.fields, g)
        ce, cb = map(similar, m.fields.B), map(similar, m.fields.E)
        yee_curl!(ce, m.fields.E, g)
        yee_curl!(cb, m.fields.B, g; field = :magnetic)
        tol = T === Float32 ? 8e-5 : 2e-13
        @test maximum(abs, divergence(ce, g, false)) <= tol * 100
        @test maximum(abs, divergence(cb, g, true)) <= tol * 100
        lhs = sum(c -> dot(Float64.(m.fields.E[c]), Float64.(cb[c])), 1:3)
        rhs = sum(c -> dot(Float64.(ce[c]), Float64.(m.fields.B[c])), 1:3)
        @test abs(lhs - rhs) <= tol * max(1, abs(lhs), abs(rhs))
        beforeE, beforeB = divergence(m.fields.E, g, true), divergence(m.fields.B, g, false)
        dt = T(0.4) * maxwell_cfl(m)
        before_energy = discrete_energy(m, dt)
        arrays = components(m.fields)
        for _ = 1:25
            step_maxwell!(m, dt)
        end
        @test all(arrays[i] === components(m.fields)[i] for i in eachindex(arrays))
        @test maximum(abs, divergence(m.fields.E, g, true) - beforeE) <= tol * 100
        @test maximum(abs, divergence(m.fields.B, g, false) - beforeB) <= tol * 100
        @test discrete_energy(m, dt) ≈ before_energy rtol = tol
        @test em_field_energy(m) >= 0
    end
end

function wave!(m, mode, polarization, dt)
    D, T = length(m.g.n), eltype(m.fields.E[1])
    k = [d <= D ? T(2pi * mode[d]) / m.g.L[d] : zero(T) for d = 1:3]
    kd = [d <= D ? 2sin(k[d] * m.g.dx[d] / 2) / m.g.dx[d] : zero(T) for d = 1:3]
    kappa = norm(kd)
    e = cross(kd, polarization)
    e ./= norm(e)
    b = cross(kd, e) * (sqrt(1 - (m.c * dt * kappa / 2)^2) / (m.c * kappa))
    theta = 2asin(m.c * dt * kappa / 2)
    for c = 1:3
        fill_component!(
            m.fields.E[c],
            m.g,
            :electric,
            c,
            x -> e[c] * cos(sum(k[d] * x[d] for d = 1:D)),
        )
        fill_component!(
            m.fields.B[c],
            m.g,
            :magnetic,
            c,
            x -> b[c] * cos(sum(k[d] * x[d] for d = 1:D)),
        )
    end
    return k, e, b, theta
end
@testset "Periodic traveling modes: independent discrete dispersion" begin
    for T in (Float32, Float64), D = 1:3, direction in (-1, 1), harmonic in (1, 2)
        g = YeeGrid(
            ntuple(d -> 12 + d, D),
            ntuple(d -> T(1 + 0.2d), D);
            boundaries = ntuple(_ -> :periodic, D),
        )
        m = MaxwellVacuum(g; c = T(1.3))
        dt = T(0.45) * maxwell_cfl(m)
        mode = ntuple(d -> direction * (d == 1 ? harmonic : 1), D)
        k, e, b, theta = wave!(m, mode, T[0, 0, 1], dt)
        for _ = 1:23
            step_maxwell!(m, dt)
        end
        tolerance = T === Float32 ? 8e-5 : 3e-13
        for c = 1:3, field in (:electric, :magnetic)
            a = field === :electric ? m.fields.E[c] : m.fields.B[c]
            reference = similar(a)
            amplitude = field === :electric ? e[c] : b[c]
            fill_component!(
                reference,
                g,
                field,
                c,
                x -> amplitude * cos(sum(k[d] * x[d] for d = 1:D) - 23theta),
            )
            @test maximum(abs, a - reference) <= tolerance
        end
        flux = yee_poynting(m)
        @test sign(sum(flux[1])) == direction
    end
end

function cavity!(m, polarization, harmonic)
    D, T = length(m.g.n), eltype(m.fields.E[1])
    k = ntuple(
        d ->
            d == polarization ? zero(T) : T((m.g.periodic[d] ? 2pi : pi) * harmonic) / m.g.L[d],
        D,
    )
    kd = ntuple(d -> 2sin(k[d] * m.g.dx[d] / 2) / m.g.dx[d], D)
    shape = x -> prod(d -> d == polarization ? one(T) : sin(k[d] * x[d]), 1:D)
    fill_component!(m.fields.E[polarization], m.g, :electric, polarization, shape)
    apply_yee_pec!(m.fields, m.g)
    exactcurl = map(similar, m.fields.B)
    for c = 1:3
        j, q = mod1(c + 1, 3), mod1(c + 2, 3)
        d, sign = q == polarization ? (j, 1) : j == polarization ? (q, -1) : (4, 0)
        fill_component!(
            exactcurl[c],
            m.g,
            :magnetic,
            c,
            x ->
                d > D ? zero(T) :
                sign *
                kd[d] *
                prod(
                    a -> a == polarization ? one(T) : a == d ? cos(k[a] * x[a]) : sin(k[a] * x[a]),
                    1:D,
                ),
        )
    end
    # The independent analytic mode has exactly zero constrained values.
    f = YeeFields(map(copy, m.fields.E), exactcurl)
    apply_yee_pec!(f, m.g)
    return sqrt(sum(abs2, kd)), map(copy, m.fields.E), exactcurl
end
@testset "PEC/mixed cavity modes and corner constraints" begin
    for T in (Float32, Float64), D = 1:3, mask = 0:(2^D-1), c = 1:3, harmonic in (1, 2)
        D == 1 && c == 1 && continue
        g = YeeGrid(
            ntuple(d -> 7 + d, D),
            ntuple(d -> T(1.2 + 0.1d), D);
            boundaries = ntuple(d -> isodd(mask >> (d - 1)) ? :periodic : :pec, D),
        )
        m = MaxwellVacuum(g; c = T(1.1))
        dt = T(0.35) * maxwell_cfl(m)
        kappa, E0, CE0 = cavity!(m, c, harmonic)
        theta = 2asin(m.c * dt * kappa / 2)
        steps = 19
        for _ = 1:steps
            step_maxwell!(m, dt)
        end
        tolerance = T === Float32 ? 4e-5 : 7e-14
        Bfactor = -sin(steps * theta) * sqrt(1 - (m.c * dt * kappa / 2)^2) / (m.c * kappa)
        for component = 1:3
            @test maximum(abs, m.fields.E[component] - cos(steps * theta) * E0[component]) <=
                  tolerance
            @test maximum(abs, m.fields.B[component] - Bfactor * CE0[component]) <= tolerance
        end
        @test maximum(abs, divergence(m.fields.E, g, true)) <= tolerance * 100
        @test maximum(abs, divergence(m.fields.B, g, false)) <= tolerance * 100
    end
end

function cavity_error(n, steps, final_time)
    g = YeeGrid((n,), (1.0,))
    m = MaxwellVacuum(g)
    kappa, E0, CE0 = cavity!(m, 2, 2)
    dt = final_time / steps
    for _ = 1:steps
        step_maxwell!(m, dt)
    end
    continuum = maximum(abs, m.fields.E[2] - cos(2pi * final_time) * E0[2])
    semidiscrete = maximum(abs, m.fields.E[2] - cos(kappa * final_time) * E0[2])
    return continuum, semidiscrete
end
@testset "Independent spatial and temporal second order" begin
    spatial = [cavity_error(n, 20n, 0.19)[1] for n in (12, 24, 48)]
    temporal = [cavity_error(24, steps, 0.19)[2] for steps in (20, 40, 80)]
    @test all(3.8 .< spatial[1:2] ./ spatial[2:3] .< 4.2)
    @test all(3.95 .< temporal[1:2] ./ temporal[2:3] .< 4.05)
end

@testset "Analytic pulse transport and conducting reflection phase" begin
    for axis = 1:3, direction in (-1, 1), periodic in (false, true)
        D = axis
        n = ntuple(d -> d == axis ? 192 : 2, D)
        g = YeeGrid(
            n,
            ntuple(_ -> 1.0, D);
            boundaries = ntuple(d -> d == axis && !periodic ? :pec : :periodic, D),
        )
        m = MaxwellVacuum(g)
        c = mod1(axis + 1, 3)
        b = mod1(axis + 2, 3)
        x0 = direction > 0 ? 0.3 : 0.7
        sigma = 0.055
        f = x -> exp(-((x - x0) / sigma)^2)
        fill_component!(m.fields.E[c], g, :electric, c, x -> f(x[axis]))
        fill_component!(m.fields.B[b], g, :magnetic, b, x -> direction * f(x[axis]))
        apply_yee_pec!(m.fields, g)
        final_time = periodic ? 0.35 : 0.9
        steps = ceil(Int, final_time / (0.75maxwell_cfl(m)))
        dt = final_time / steps
        for _ = 1:steps
            step_maxwell!(m, dt)
        end
        function analytic(x, magnetic)
            coordinate = x[axis]
            if periodic
                return (magnetic ? direction : 1) *
                       sum(f(coordinate - direction * final_time + j) for j = -2:2)
            end
            incident = f(coordinate - direction * final_time)
            reflected = direction > 0 ? f(2 - coordinate - final_time) : f(-coordinate + final_time)
            return magnetic ? direction * (incident + reflected) : incident - reflected
        end
        for (field, component, a) in ((:electric, c, m.fields.E[c]), (:magnetic, b, m.fields.B[b]))
            expected = similar(a)
            fill_component!(expected, g, field, component, x -> analytic(x, field === :magnetic))
            @test maximum(abs, a - expected) < 0.018
            @test norm(a - expected) / norm(expected) < 0.02
        end
        if !periodic
            @test minimum(m.fields.E[c]) < -0.9
            @test sign(sum(yee_poynting(m)[axis])) == -direction
        end
    end
end
@testset "Physical energy and face/cell Poynting diagnostics" begin
    for T in (Float32, Float64), D = 1:3
        g = YeeGrid(
            ntuple(_ -> 4, D),
            ntuple(d -> T(d + 1), D);
            boundaries = ntuple(_ -> :periodic, D),
        )
        m = MaxwellVacuum(g; c = T(1.5))
        E = T[1.25, -0.5, 0.75]
        B = T[-0.25, 0.5, 1.75]
        for c = 1:3
            fill!(m.fields.E[c], E[c])
            fill!(m.fields.B[c], B[c])
        end
        expected = (sum(abs2, E) + m.c^2 * sum(abs2, B)) * prod(g.L) / 2
        @test em_field_energy(m) ≈ expected rtol = 10eps(T)
        S = m.c^2 * cross(E, B)
        @test all(all(==(S[c]), yee_poynting(m)[c]) for c = 1:3)
        flux = yee_poynting_flux(m)
        for d = 1:D
            area = prod(j -> j == d ? one(T) : g.L[j], 1:D)
            @test flux[d][1] ≈ -S[d] * area
            @test flux[d][2] ≈ S[d] * area
            @test iszero(sum(flux[d]))
        end
        # The same quadrature has zero PEC boundary flux, computed from fields.
        gp = YeeGrid(g.n, g.L)
        mp = MaxwellVacuum(gp; c = m.c)
        for c = 1:3
            fill!(mp.fields.E[c], E[c])
            fill!(mp.fields.B[c], B[c])
        end
        apply_yee_pec!(mp.fields, gp)
        @test all(face -> all(iszero, face), yee_poynting_flux(mp))
        exactsum =
            sum(a -> sum(x -> Rational{BigInt}(x)^2, a), mp.fields.E) +
            Rational{BigInt}(mp.c)^2 * sum(a -> sum(x -> Rational{BigInt}(x)^2, a), mp.fields.B)
        exact = exactsum * prod(Rational{BigInt}.(gp.dx)) / 2
        @test em_field_energy(mp) ≈ T(exact) rtol = 10eps(T)
    end
    # Cancellation of overflowing individual products retains a finite flux.
    g = YeeGrid((4,), (1.0,); boundaries = (:periodic,))
    m = MaxwellVacuum(g)
    fill!(m.fields.E[2], 1e200)
    fill!(m.fields.E[3], 1e200)
    fill!(m.fields.B[2], 1e200)
    fill!(m.fields.B[3], 1e200)
    @test all(iszero, yee_poynting(m)[1])
    @test all(face -> all(iszero, face), yee_poynting_flux(m))
end

@testset "Variable step centering and energy-range completion" begin
    for T in (Float32, Float64)
        g = YeeGrid((20,), (T(1),))
        m = MaxwellVacuum(g; c = T(1.2))
        kappa, E0, CE0 = cavity!(m, 2, 2)
        # Independent scalar oscillator update, evaluated at higher precision.
        q, p = BigFloat(1), BigFloat(0)
        for ratio in (0.2, 0.65, 0.15, 0.4, 0.3, 0.7, 0.12)
            dt = T(ratio) * maxwell_cfl(m)
            a = BigFloat(dt) * BigFloat(m.c) * BigFloat(kappa)
            p -= a * q / 2
            q += a * p
            p -= a * q / 2
            step_maxwell!(m, dt)
        end
        tolerance = T === Float32 ? 2e-6 : 3e-15
        @test maximum(abs, m.fields.E[2] - T(q) * E0[2]) < tolerance
        @test maximum(abs, m.fields.B[3] - T(p / (BigFloat(m.c) * BigFloat(kappa))) * CE0[3]) <
              tolerance
        for amplitude in (sqrt(floatmax(T) / T(8)), nextfloat(zero(T)))
            gm = YeeGrid((4,), (T(1),); boundaries = (:periodic,))
            mm = MaxwellVacuum(gm)
            fill!(mm.fields.E[1], amplitude)
            exact = Rational{BigInt}(amplitude)^2 / 2
            @test em_field_energy(mm) === T(exact)
        end
    end
end

@testset "Periodic cross-check against existing spectral vacuum EMPIC" begin
    for D = 1:3
        errors = Float64[]
        for n in (16, 32)
            shape = ntuple(d -> d == 1 ? n : 2, D)
            lengths = ntuple(_ -> 1.0, D)
            g = YeeGrid(shape, lengths; boundaries = ntuple(_ -> :periodic, D))
            m = MaxwellVacuum(g)
            spectral = EMPIC(FourierGrid(shape, lengths), 0; n0 = 0.0, c = 1.0)
            electrons = ParticleSet{D,Float64}(0; q = -1.0, m = 1.0)
            init_empic!(spectral, electrons)
            for I in CartesianIndices(spectral.E[2])
                value = sinpi(2 * (I[1] - 1) / n)
                m.fields.E[2][I] = value
                spectral.E[2][I] = value
            end
            steps = ceil(Int, 0.19 / (0.2maxwell_cfl(m)))
            dt = 0.19 / steps
            for _ = 1:steps
                step_maxwell!(m, dt)
                step_empic!(spectral, electrons, dt)
            end
            push!(errors, maximum(abs, m.fields.E[2] - spectral.E[2]))
        end
        @test 3.8 < errors[1] / errors[2] < 4.2
    end
end

@testset "Stored wave speed retains square before later rescaling" begin
    for T in (Float32, Float64)
        c = sqrt(floatmin(T)) / T(1000)
        g = YeeGrid((4,), (one(T),); boundaries = (:periodic,))
        m = MaxwellVacuum(g; c)
        amplitude = T === Float32 ? T(1e10) : T(1e100)
        fill!(m.fields.E[2], amplitude)
        fill!(m.fields.B[3], amplitude)
        expected = T(Rational{BigInt}(c)^2 * Rational{BigInt}(amplitude)^2)
        @test all(x -> isapprox(x, expected; rtol = 8eps(T), atol = zero(T)), yee_poynting(m)[1])
        length = inv(c)
        g = YeeGrid((4,), (length,); boundaries = (:periodic,))
        m = MaxwellVacuum(g; c)
        m.fields.B[3] .= T[length, -length, length, -length]
        cb = map(similar, m.fields.E)
        yee_curl!(cb, m.fields.B, g; field = :magnetic)
        dt = T === Float32 ? T(1e20) : T(1e150)
        expected =
            map(x -> T(Rational{BigInt}(dt) * Rational{BigInt}(c)^2 * Rational{BigInt}(x)), cb[2])
        step_maxwell!(m, dt)
        @test all(isapprox.(m.fields.E[2], expected; rtol = 8eps(T), atol = zero(T)))
    end
end

@testset "Curl and divergence complete overflowing differences exactly" begin
    for T in (Float32, Float64), sign in (-1, 1)
        # Differences of ±floatmax overflow; the exact quotient by dx=4 is finite.
        g = YeeGrid((2,), (T(8),); boundaries = (:periodic,))
        f = YeeFields(g)
        f.E[2] .= T(sign) .* T[floatmax(T), -floatmax(T)]
        f.E[1] .= T(sign) .* T[floatmax(T), -floatmax(T)]
        curl = map(similar, f.B)
        yee_curl!(curl, f.E, g)
        exact = T(sign) * T(Rational{BigInt}(floatmax(T)) / 2)
        @test curl[3] == [-exact, exact]
        @test all(iszero, curl[1]) && all(iszero, curl[2])
        divergence = zeros(T, 2)
        yee_divergence!(divergence, f.E, g; field = :electric)
        @test divergence == [exact, -exact]
        # An exact value outside the range is a signed infinity, not NaN.
        g1 = YeeGrid((2,), (T(1),); boundaries = (:periodic,))
        f1 = YeeFields(g1)
        f1.E[2] .= T(sign) .* T[floatmax(T), -floatmax(T)]
        curl1 = map(similar, f1.B)
        yee_curl!(curl1, f1.E, g1)
        @test isequal(curl1[3], T(sign) .* T[-Inf, Inf])
        f1.E[2][1] = T(NaN)
        @test_throws ArgumentError yee_curl!(curl1, f1.E, g1)
        @test_throws ArgumentError yee_divergence!(zeros(T, 2), f1.E, g1; field = :electric)
    end
end

@noinline allocated_step(m, dt) = @allocated step_maxwell!(m, dt)

@testset "Bitwise-deterministic, allocation-free native steps" begin
    for T in (Float32, Float64), D = 1:3, mask = 0:(2^D-1)
        g = YeeGrid(
            ntuple(d -> 5 + d, D),
            ntuple(d -> T(1 + 0.1d), D);
            boundaries = ntuple(d -> isodd(mask >> (d - 1)) ? :periodic : :pec, D),
        )
        models = ntuple(_ -> MaxwellVacuum(g; c = T(1.3)), 2)
        for m in models
            rng = MersenneTwister(4021 + D + 8mask)
            foreach(a -> randn!(rng, a), components(m.fields))
            apply_yee_pec!(m.fields, g)
        end
        dt = T(0.45) * maxwell_cfl(models[1])
        for _ = 1:5, m in models
            step_maxwell!(m, dt)
        end
        @test all(isequal.(components(models[1].fields), components(models[2].fields)))
        @test allocated_step(models[1], dt) == 0
    end
end

include("bounded_maxwell/range_oracle.jl")
include("bounded_maxwell/scaled_stage_oracle.jl")
include("bounded_maxwell/poynting_oracle.jl")
include("bounded_maxwell/cancellation_oracle.jl")
include("bounded_maxwell/accuracy_oracle.jl")

@testset "Conservative complete-range CFL certificate" begin
    for T in (Float32, Float64),
        D = 1:3,
        scale in (T(1.5) * floatmin(T), T(0.7), T(10), floatmax(T) / T(16)),
        speed in (nextfloat(zero(T)), T(0.3), T(7), floatmax(T) / T(16))

        g = YeeGrid(ntuple(_ -> 4, D), ntuple(d -> scale * T(1 + 0.125d), D))
        m = MaxwellVacuum(g; c = speed)
        exact_square = inv(Rational{BigInt}(speed)^2 * sum(dx -> inv(Rational{BigInt}(dx)^2), g.dx))
        bound = maxwell_cfl(m)
        @test isfinite(bound) && bound >= 0
        @test Rational{BigInt}(bound)^2 <= exact_square
        @test bound == floatmax(T) || Rational{BigInt}(nextfloat(bound))^2 > exact_square
        if bound > nextfloat(zero(T))
            dt = prevfloat(bound)
            step_maxwell!(m, dt)
            @test m.time[] == dt
        end
    end
    for T in (Float32, Float64)
        g = YeeGrid((4, 4, 4), ntuple(_ -> T(1.5) * floatmin(T), 3))
        m = MaxwellVacuum(g)
        @test maxwell_cfl(m) > 0
        step_maxwell!(m, g.dx[1] / 4)
        @test em_field_energy(m) === zero(T)
    end
end

end
