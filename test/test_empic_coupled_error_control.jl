module EMPICCoupledErrorControlTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC
const Q = Rational{BigInt}
cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

function control(
    ::Type{T};
    hmax = T(1 / 64),
    hmin = hmax / T(64),
    tol = T(1 / 4),
    trials = 128,
) where {T}
    return (
        pair = :projected_euler_heun_v1,
        metric = :component_max_v1,
        domain = :cic_tsc_continuous_v1,
        atol = (x = tol, v = tol, E = tol, B = tol, clock = T(128) * eps(T)),
        rtol = (x = zero(T), v = zero(T), E = zero(T), B = zero(T), clock = zero(T)),
        h_min = hmin,
        h_max = hmax,
        max_trials = trials,
        max_force_evals = 100_000,
        max_deposits = 1_000_000,
        max_fft_calls = 100_000,
        max_report_cells = 2_000_000,
    )
end
function population(
    ::Type{T},
    D;
    weighted = false,
    count = 1,
    knot = nothing,
    seam = false,
) where {T}
    rows = [
        ParticleSet{D,T}(count; q = -2, m = 3),
        ParticleSet{D,T}(1; q = 1, m = 5),
        ParticleSet{D,T}(1; q = 0, m = 1),
        ParticleSet{D,T}(0; q = 3, m = 7),
    ]
    for (s, ps) in enumerate(rows), p in eachindex(ps.weight), d = 1:D
        ps.x[d][p] = T(s == 1 ? (d == 1 ? 7 / 32 : 5 / 16) : (s == 2 ? 21 / 32 : 3 / 4))
        ps.weight[p] = weighted && s <= 2 ? (s == 1 ? T(1 / 2) / T(count) : one(T)) : zero(T)
        if s <= 2
            ps.v[d][p] = T(s == 1 ? 1 / 8 : -1 / 16)
        end
    end
    for ps in rows[1:2]
        fill!(ps.v[2], T(1 / 32))
    end
    if knot !== nothing
        fill!(rows[1].x[1], T(knot))
    end
    if seam
        fill!(rows[1].x[1], prevfloat(one(T)))
        fill!(rows[1].v[1], T(1 / 8))
    end
    return rows
end
function fixture(
    ::Type{T},
    Model;
    D = 1,
    shape = CIC(),
    weighted = false,
    proposal = control(T),
    count = 1,
    knot = nothing,
    seam = false,
    n_sub = 1,
    relativistic = false,
    initialize = true,
) where {T}
    n = D == 1 ? (4,) : (D == 2 ? (4, 4) : (4, 3, 3))
    g = FourierGrid(n, ntuple(_ -> one(T), D))
    ps = population(T, D; weighted, count, knot, seam)
    if Model === EMPIC
        fill!(ps[1].v[3], T(1 / 64))
        fill!(ps[2].v[3], T(-1 / 128))
    end
    es = Model(
        g,
        ps;
        mobile = Bool[true, true, false, true],
        c = one(T),
        shape,
        n_sub,
        species_dt = T[1/2, 3/4, 0, 1/2],
        event_budget = 256,
        error_control = proposal,
        relativistic,
    )
    if es isa EMPIC
        fill!(es.E[1], T(1 / 32))
        for I in CartesianIndices(g.n)
            es.E[2][I] = T(1 / 16) * T(cos(2pi * (I[1] - 1) / g.n[1]))
            es.B[3][I] = T(1 / 32) * T(sin(2pi * (I[1] - 1) / g.n[1]))
        end
    else
        fill!(es.Ex, T(1 / 32))
        for i in eachindex(es.Ey)
            es.Ey[i] = T(1 / 16) * T(cos(2pi * (i - 1) / g.n[1]))
            es.Bz[i] = T(1 / 32) * T(sin(2pi * (i - 1) / g.n[1]))
        end
    end
    initialize && init_empic!(es, ps)
    return es, ps
end
function eb(es)
    es isa EMPIC && return (E = map(copy, es.E), B = map(copy, es.B))
    z = zeros(eltype(es.Ex), size(es.Ex))
    return (E = (copy(es.Ex), copy(es.Ey), copy(z)), B = (copy(z), copy(z), copy(es.Bz)))
end
setfuture(es; kwargs...) =
    es isa EMPIC ? set_empic_schedule!(es; kwargs...) : set_empic1d_schedule!(es; kwargs...)
snapshot(es, ps) = H._checkpoint_state(es, ps)
physical(es, ps) =
    deepcopy((fields = eb(es), particles = [H._checkpoint_particle_state(p) for p in ps]))
function leaf_view(a)
    backing = zeros(eltype(a), length(a) + 2)
    result = view(backing, 2:length(a)+1)
    copyto!(result, a)
    return result
end

# Full-node cardinal B-splines, independent of production stencil selection.
# Only tiny n>=3 references use these finite periodic images.
cardinal(::CIC, r) = max(0.0, 1 - abs(r))
cardinal(::TSC, r) = abs(r) <= 0.5 ? 0.75 - r * r : (abs(r) < 1.5 ? (1.5 - abs(r))^2 / 2 : 0.0)
function nodal_weights(ps, p, g, shape)
    s = ntuple(d -> Float64(mod(ps.x[d][p], g.L[d]) / g.dx[d]), length(g.n))
    return [
        prod(
            d -> sum(cardinal(shape, s[d] - (I[d] - 1) + image * g.n[d]) for image = -1:1),
            1:length(g.n),
        ) for I in CartesianIndices(g.n)
    ]
end
function sources(ps, g, shape, mobile)
    rho = zeros(g.n)
    J = ntuple(_ -> zeros(g.n), 3)
    volume = prod(Float64, g.dx)
    for (s, row) in enumerate(ps), p in eachindex(row.weight)
        weights = nodal_weights(row, p, g, shape)
        charge = Float64(row.q) * Float64(row.weight[p]) / volume
        rho .+= charge .* weights
        if mobile[s]
            for c = 1:3
                J[c] .+= (charge * Float64(row.v[c][p])) .* weights
            end
        end
    end
    return rho, J
end
function sample(fields, ps, p, g, shape)
    weights = nodal_weights(ps, p, g, shape)
    return ntuple(c -> sum(Float64.(fields[c]) .* weights), 3)
end
# Direct finite DFT; no production FFT/curl/projection is called by this oracle.
function dft(values, g; inverse = false)
    result = Array{ComplexF64}(undef, g.n)
    sign = inverse ? 1 : -1
    for I in CartesianIndices(g.n)
        value = 0.0 + 0.0im
        for J in CartesianIndices(g.n)
            angle = sign * 2pi * sum(d -> (I[d] - 1) * (J[d] - 1) / g.n[d], 1:length(g.n))
            value += values[J] * cis(angle)
        end
        result[I] = inverse ? value / prod(g.n) : value
    end
    return result
end
wavevector(g, I) = ntuple(d -> d <= length(g.n) ? Float64(g.kvec[d][I[d]]) : 0.0, 3)
function curl(fields, g)
    hats = map(a -> dft(a, g), fields)
    out = map(a -> similar(a), hats)
    for I in CartesianIndices(g.n)
        value = cross3(wavevector(g, I), ntuple(c -> hats[c][I], 3))
        for c = 1:3
            out[c][I] = 1im * value[c]
        end
    end
    return map(a -> real.(dft(a, g; inverse = true)), out)
end
function corrected_current(raw, rho0, rho1, g, h)
    hats = map(a -> dft(a, g), raw)
    delta = dft(rho1 .- rho0, g)
    dc = first(CartesianIndices(g.n))
    for I in CartesianIndices(g.n)
        k = wavevector(g, I)
        k2 = sum(abs2, k)
        if k2 > 0
            longitudinal = sum(c -> k[c] * hats[c][I], 1:3)
            for c = 1:3
                hats[c][I] += k[c] * (1im * delta[I] / h - longitudinal) / k2
            end
        elseif I != dc
            for c = 1:3
                hats[c][I] = 0
            end
        end
    end
    return map(a -> real.(dft(a, g; inverse = true)), hats)
end
function gradient_lift(rho, g)
    hat = dft(rho, g)
    out = ntuple(_ -> zeros(ComplexF64, g.n), 3)
    for I in CartesianIndices(g.n)
        k = wavevector(g, I)
        k2 = sum(abs2, k)
        iszero(k2) && continue
        for c = 1:3
            out[c][I] = -1im * k[c] * hat[I] / k2
        end
    end
    return map(a -> real.(dft(a, g; inverse = true)), out)
end
function force_parts(ps, p, E, B, h)
    v = ntuple(c -> Q(ps.v[c][p]), 3)
    electric = map(Q, E)
    magnetic = map(Q, B)
    return (Q(h) * Q(ps.q) / Q(ps.m)) .* (electric .+ cross3(v, magnetic))
end
function wrapped(::Type{T}, x, L) where {T}
    value = T(mod(x, Q(L)))
    return value == L ? zero(T) : value
end
function pair_reference(es, particles, h)
    T = eltype(es.rho_n)
    g, shape, mobile = es.g, es.shape, es.collection.mobile
    start = eb(es)
    old = deepcopy(particles)
    low = deepcopy(old)
    high = deepcopy(old)
    increments = NTuple{3,Q}[]
    for (s, ps) in enumerate(old), p in eachindex(ps.weight)
        electric = T.(sample(start.E, ps, p, g, shape))
        magnetic = T.(sample(start.B, ps, p, g, shape))
        increment = mobile[s] ? force_parts(ps, p, electric, magnetic, h) : ntuple(_ -> zero(Q), 3)
        push!(increments, increment)
        mobile[s] || continue
        for c = 1:3
            low[s].v[c][p] = T(Q(ps.v[c][p]) + increment[c])
        end
        for d = 1:length(g.n)
            low[s].x[d][p] = wrapped(T, Q(ps.x[d][p]) + Q(h) * Q(ps.v[d][p]), g.L[d])
        end
    end
    rho0, J0 = sources(old, g, shape, mobile)
    rhoE, _ = sources(low, g, shape, mobile)
    rho0, rhoE = T.(rho0), T.(rhoE)
    J0 = map(a -> T.(a), J0)
    JE = corrected_current(J0, rho0, rhoE, g, Float64(h))
    curlE0, curlB0 = curl(start.E, g), curl(start.B, g)
    BE = ntuple(c -> T.(start.B[c] .- Float64(h) .* curlE0[c]), 3)
    EE = ntuple(
        c -> T.(start.E[c] .+ Float64(h) * Float64(es.c)^2 .* curlB0[c] .- Float64(h) .* JE[c]),
        3,
    )
    cursor = 0
    for (s, ps) in enumerate(low), p in eachindex(ps.weight)
        cursor += 1
        mobile[s] || continue
        electric = T.(sample(EE, ps, p, g, shape))
        magnetic = T.(sample(BE, ps, p, g, shape))
        increment = force_parts(ps, p, electric, magnetic, h)
        for c = 1:3
            high[s].v[c][p] = T(Q(old[s].v[c][p]) + (increments[cursor][c] + increment[c]) / 2)
        end
        for d = 1:length(g.n)
            high[s].x[d][p] = wrapped(
                T,
                Q(old[s].x[d][p]) + Q(h) * (Q(old[s].v[d][p]) + Q(ps.v[d][p])) / 2,
                g.L[d],
            )
        end
    end
    rhoH, _ = sources(high, g, shape, mobile)
    _, J1 = sources(low, g, shape, mobile)
    rhoH = T.(rhoH)
    J1 = map(a -> T.(a), J1)
    Jbar = ntuple(c -> T.((Float64.(J0[c]) .+ Float64.(J1[c])) ./ 2), 3)
    JH = corrected_current(Jbar, rho0, rhoH, g, Float64(h))
    curlEE, curlBE = curl(EE, g), curl(BE, g)
    BH = ntuple(c -> T.(start.B[c] .- Float64(h) / 2 .* (curlE0[c] .+ curlEE[c])), 3)
    EH = ntuple(
        c ->
            T.(
                start.E[c] .+ Float64(h) * Float64(es.c)^2 / 2 .* (curlB0[c] .+ curlBE[c]) .-
                Float64(h) .* JH[c]
            ),
        3,
    )
    return (
        low = low,
        high = high,
        E_low = EE,
        B_low = BE,
        E_high = EH,
        B_high = BH,
        J_low = JE,
        J_high = JH,
        Jraw_high = Jbar,
        rho0 = rho0,
        rho_low = rhoE,
        rho_high = rhoH,
        longitudinal = gradient_lift(Float64.(rhoH) .- Float64.(rhoE), g),
    )
end
function check_particles(actual, expected, tolerance)
    for (a, b) in zip(actual, expected)
        @test a.q == b.q && a.m == b.m && a.weight == b.weight && a.id == b.id && a.tag == b.tag
        for (u, v) in zip((a.x..., a.v...), (b.x..., b.v...))
            @test isapprox(u, v; rtol = tolerance, atol = tolerance)
        end
    end
end
function check_fields(actual, E, B, tolerance)
    if hasproperty(actual, :E)
        for c = 1:3
            @test isapprox(actual.E[c], E[c]; rtol = tolerance, atol = tolerance)
            @test isapprox(actual.B[c], B[c]; rtol = tolerance, atol = tolerance)
        end
    else
        for (a, b) in ((actual.Ex, E[1]), (actual.Ey, E[2]), (actual.Bz, B[3]))
            @test isapprox(a, b; rtol = tolerance, atol = tolerance)
        end
    end
end

@testset "Projected Euler-Heun independent coupled reference" begin
    cases = vec([
        (T, Model, 1, shape) for T in (Float32, Float64), Model in (EMPIC, EMPIC1D),
        shape in (CIC(), TSC())
    ])
    append!(cases, [(Float64, EMPIC, D, CIC()) for D in (2, 3)])
    for (T, Model, D, shape) in cases
        @testset "$T $Model D=$D $(typeof(shape))" begin
            h = T(1 / 64)
            es, ps = fixture(T, Model; D, shape, weighted = true)
            @test ps[1].v[1][1] == T(1 / 8) # Initialization is physically full time.
            original = snapshot(es, ps)
            step_empic!(es, ps, -zero(T))
            @test snapshot(es, ps) == original
            before = eb(es)
            reference = pair_reference(es, ps, h)
            step_empic!(es, ps, h)
            d = electromagnetic_diagnostics(es, ps)
            r = d.last_error_control_report
            @test d.schedule_policy === 3 && d.phase_policy === :full_time_coupled
            @test all(iszero, d.velocity_time_offsets) && iszero(d.magnetic_time_offset)
            @test es.step[] == 1 && es.last_dt[] == h && es.last_sub_dt[] == h
            @test d.source_time == Q(h) && d.public_clock_defect == Q(es.time[]) - d.source_time
            @test length(r.intervals) == 1 && only(r.intervals).dt == h
            @test r.costs.trials == r.admission_costs.trials + r.verification_costs.trials + 1
            @test r.costs.force_evals >= 4 &&
                  r.costs.fft_calls > 0 &&
                  only(r.intervals).costs.trials == 1
            interval = only(r.intervals)
            @test interval.source_before == 0 && interval.source_after == Q(h) && interval.norm <= 1
            tolerance = T === Float32 ? 1024eps(T) : 2e-10
            check_particles(ps, reference.high, tolerance)
            check_particles(interval.low.species, reference.low, tolerance)
            check_fields(interval.low.fields, reference.E_low, reference.B_low, tolerance)
            check_fields(interval.high.fields, reference.E_high, reference.B_high, tolerance)
            current = es isa EMPIC ? es.J : (es.Jx, es.Jy, zeros(T, es.g.n))
            for c = 1:3
                @test isapprox(current[c], reference.J_high[c]; rtol = tolerance, atol = tolerance)
                if es isa EMPIC
                    @test isapprox(
                        es.Jraw[c],
                        reference.Jraw_high[c];
                        rtol = tolerance,
                        atol = tolerance,
                    )
                end
            end
            @test isapprox(es.rho_np1, reference.rho_high; rtol = tolerance, atol = tolerance)
            @test maximum(a -> maximum(abs, a), reference.longitudinal) > 1e-12
            # Removing G(rhoH-rhoE) cannot satisfy this actual branch comparison.
            lowE =
                es isa EMPIC ? interval.low.fields.E :
                (interval.low.fields.Ex, interval.low.fields.Ey, zeros(T, es.g.n))
            highE =
                es isa EMPIC ? interval.high.fields.E :
                (interval.high.fields.Ex, interval.high.fields.Ey, zeros(T, es.g.n))
            for c = 1:3
                @test isapprox(
                    Float64.(highE[c]) .- Float64.(lowE[c]),
                    Float64.(reference.E_high[c]) .- Float64.(reference.E_low[c]);
                    rtol = tolerance,
                    atol = tolerance,
                )
            end
            @test (
                es isa EMPIC ? r.work_raw !== nothing && r.work_raw_transverse === nothing :
                r.work_raw === nothing && r.work_raw_transverse !== nothing
            )
            @test r.total_defect == r.field_defect + r.particle_defect + r.coupling_defect
            @test r.total_defect ==
                  r.field_after - r.field_before + r.kinetic_after - r.kinetic_before
            # Actual effective-symbol null modes are distinct from physical DC.
            hats = map(a -> dft(a, es.g), current)
            dc = first(CartesianIndices(es.g.n))
            for I in CartesianIndices(es.g.n)
                if I != dc && iszero(sum(abs2, wavevector(es.g, I)))
                    @test all(c -> abs(hats[c][I]) <= 10tolerance, 1:3)
                end
            end
            for c = 1:3
                @test isapprox(
                    sum(current[c]),
                    sum(reference.Jraw_high[c]);
                    rtol = tolerance,
                    atol = tolerance,
                )
            end
            @test any(
                c ->
                    es isa EMPIC ? es.E[c] != before.E[c] || es.B[c] != before.B[c] :
                    (c == 1 ? es.Ex != before.E[1] : es.Ey != before.E[2]),
                1:3,
            )
        end
    end
end

# Independent4096-bit evaluation of the physical-v formula. These remote checks
# do not call the production radical combiner/finisher as their reference.
function radical_reference(v0, parts0, parts1, ::Type{T}) where {T}
    return setprecision(BigFloat, 4096) do
        ntuple(3) do c
            a = BigFloat(parts0.coefficient[c]) * sqrt(BigFloat(parts0.gap))
            increment =
                parts1 === nothing ? a :
                (a + BigFloat(parts1.coefficient[c]) * sqrt(BigFloat(parts1.gap))) / 2
            T(BigFloat(v0[c]) + increment)
        end
    end
end
function physical_parts(v, E, B, q, m, h, c)
    V, e, b = map(x -> map(Q, x), (v, E, B))
    gap = (Q(c)^2 - sum(abs2, V)) / Q(c)^2
    coefficient = (Q(h) * Q(q) / Q(m)) .* (e .+ cross3(V, b) .- V .* (sum(V .* e) / Q(c)^2))
    return (coefficient = coefficient, gap = gap)
end
@testset "Original-factor and canonical radical scalar outputs" begin
    for T in (Float32, Float64)
        z = zero(T)
        v = (z, -z, z)
        E = (nextfloat(z), z, z)
        B = (z, z, z)
        row = ParticleSet{1,T}(0; q = floatmax(T), m = floatmin(T))
        factor = H._validated_boris_factor(row)
        h = Q(nextfloat(z))
        parts = H._em_lorentz_scaled_parts(v, E, B, factor, h, one(T), Val(false))
        expected = ntuple(c -> T(Q(v[c]) + h * Q(row.q) / Q(row.m) * Q(E[c])), 3)
        @test H._em_lorentz_pair_output(v, parts, nothing, one(T), Val(false)) == expected
        inert = H._em_lorentz_scaled_parts(v, B, B, factor, h, one(T), Val(true))
        @test isequal(H._em_lorentz_pair_output(v, inert, nothing, one(T), Val(true)), v)
        tiny = ParticleSet{1,T}(0; q = nextfloat(z), m = floatmax(T))
        rescued = H._em_lorentz_scaled_parts(
            v,
            (one(T), z, z),
            B,
            H._validated_boris_factor(tiny),
            Q(tiny.m),
            one(T),
            Val(false),
        )
        @test H._em_lorentz_pair_output(v, rescued, nothing, one(T), Val(false))[1] == nextfloat(z)
        # sqrt(3/4)=(3/2)*sqrt(1/3): two distinct radicands cancel exactly.
        cancel0 = (coefficient = (one(Q), one(Q), one(Q)), gap = Q(3) / 4)
        cancel1 = (coefficient = ntuple(_ -> -Q(3) / 2, 3), gap = Q(1) / 3)
        @test isequal(H._em_lorentz_pair_output(v, cancel0, cancel1, T(6), Val(true)), v)
        # c=5,v=3 has rational sqrt(g)=4/5. The parallel increment is an
        # exact midpoint at3, then just ABOVE it; raw uncanonicalized roots fail.
        v = (T(3), z, z)
        E = (one(T), z, z)
        B = (z, z, z)
        c = T(5)
        factor = H._validated_boris_factor(ParticleSet{1,T}(0))
        h = Q(125) / 64 * Q(eps(T(3))) / 2
        for electric in (E, (nextfloat(one(T)), z, z))
            own = physical_parts(v, electric, B, one(T), one(T), h, c)
            actual = H._em_lorentz_scaled_parts(v, electric, B, factor, h, c, Val(true))
            @test actual == own
            expected = electric[1] == one(T) ? T(3) : nextfloat(T(3))
            @test first(H._em_lorentz_pair_output(v, actual, nothing, c, Val(true))) == expected
        end
        E, B = (T(1 / 8), T(1 / 32), z), (z, z, T(1 / 16))
        h = Q(1) / 64
        p0 = physical_parts(v, E, B, one(T), one(T), h, c)
        low = radical_reference(v, p0, nothing, T)
        p1 = physical_parts(low, E, B, one(T), one(T), h, c)
        a0 = H._em_lorentz_scaled_parts(v, E, B, factor, h, c, Val(true))
        a1 = H._em_lorentz_scaled_parts(low, E, B, factor, h, c, Val(true))
        @test H._em_lorentz_pair_output(v, a0, a1, c, Val(true)) == radical_reference(v, p0, p1, T)
        @test_throws ArgumentError H._em_lorentz_scaled_parts(
            (one(T), z, z),
            E,
            B,
            factor,
            h,
            one(T),
            Val(true),
        )
        invalid = (coefficient = (Q(2), zero(Q), zero(Q)), gap = one(Q))
        outside = H._em_lorentz_pair_output((z, z, z), invalid, nothing, one(T), Val(true))
        @test !H._supported_relativistic_output(outside, one(T)) # Publication owns rejection.
        near = (prevfloat(one(T)), z, z)
        stationary =
            H._em_lorentz_scaled_parts(near, (z, z, z), (z, z, z), factor, h, one(T), Val(true))
        @test H._em_lorentz_pair_output(near, stationary, nothing, one(T), Val(true)) == near
    end
end

@testset "Continuous cuts, seam and pending full-time force" begin
    for T in (Float32, Float64), Model in (EMPIC, EMPIC1D), shape in (CIC(), TSC())
        knot = shape isa CIC ? T(1 / 4) : T(3 / 8)
        for seam in (false, true)
            es, ps = fixture(T, Model; shape, knot, seam)
            step_empic!(es, ps, T(1 / 64))
            d = electromagnetic_diagnostics(es, ps)
            @test d.last_deadline[1] == 0 && d.next_deadline[1] == Q(1) / 2
            interval = only(d.last_error_control_report.intervals)
            @test interval.stage1.stencil[1] == interval.lifts.low[1].stencil
            if seam
                @test ps[1].x[1][1] < T(1 / 16) && interval.lifts.high[1].winding[1][1] == 1
                @test interval.lifts.high[1].cuts[1][1].seam
            else
                knot_stencil = interval.stage0.stencil[1][1][1]
                @test knot_stencil.base == 1
                @test knot_stencil.weights ==
                      (shape isa CIC ? (one(T), zero(T)) : (T(1 / 2), T(1 / 2), zero(T)))
            end
            pending = deepcopy(d.next_deadline)
            before = deepcopy(ps[1].v)
            reference = pair_reference(es, ps, T(1 / 64))
            step_empic!(es, ps, T(1 / 64))
            d = electromagnetic_diagnostics(es, ps)
            @test d.next_deadline == pending && isempty(d.last_error_control_report.decisions)
            @test ps[1].v != before
            check_particles(ps, reference.high, T === Float32 ? 1024eps(T) : 2e-10)
        end
        # Closed new-policy restriction is distinct from old NGP support.
        @test_throws ArgumentError fixture(T, Model; shape = NGP())
        oldps = population(T, 1)
        old = Model(
            FourierGrid((4,), (one(T),)),
            oldps;
            shape = NGP(),
            mobile = Bool[true, true, false, true],
            species_dt = T[1/2, 3/4, 0, 1/2],
            event_budget = 32,
        )
        init_empic!(old, oldps)
        step_empic!(old, oldps, T(1 / 64))
        @test H._checkpoint_state(old, oldps).schema_version === 3
    end
end

@testset "Controller rejection, late rollback and complete schema4 restart" begin
    for T in (Float32, Float64), Model in (EMPIC, EMPIC1D)
        h = T(1 / 32)
        tight = control(T; hmax = h, tol = T(1 / 4096))
        es, ps = fixture(T, Model; proposal = tight)
        step_empic!(es, ps, h)
        d = electromagnetic_diagnostics(es, ps)
        report = d.last_error_control_report
        @test !isempty(report.rejections) && length(report.intervals) > 1
        @test report.costs.trials ==
              report.admission_costs.trials +
              report.verification_costs.trials +
              length(report.intervals) +
              length(report.rejections)
        @test d.actual_elapsed == sum(Q(i.dt) for i in report.intervals)
        oldreport = deepcopy(report)
        pending = deepcopy(d.next_deadline)
        # A future smaller trial budget must not invalidate the historical
        # report that actually used a larger budget and retained rejections.
        future = control(T; trials = 1)
        setfuture(es; error_control = future, species_dt = T[3/4, 1, 0, 1/2])
        d = electromagnetic_diagnostics(es, ps)
        hint = d.controller_next_h
        setfuture(
            es;
            error_control = merge(
                future,
                (h_min = min(future.h_min, hint / T(4)), h_max = hint / T(2)),
            ),
        )
        setfuture(es; error_control = future)
        d = electromagnetic_diagnostics(es, ps)
        @test d.controller_next_h == hint
        io = IOBuffer()
        serialize(io, oldreport)
        oldbytes = take!(io)
        io = IOBuffer()
        serialize(io, d.last_error_control_report)
        @test take!(io) == oldbytes && d.next_deadline == pending
        saved = H._checkpoint_state(es, ps)
        before = snapshot(es, ps)
        @test saved.schema_version === 4 &&
              saved.schedule_policy === 3 &&
              saved.phase_policy === :full_time_coupled
        @test_throws ArgumentError H._restore_checkpoint_state!(
            es,
            ps,
            merge(saved, (phase_policy = :persistent_event,)),
        )
        @test_throws ArgumentError H._restore_checkpoint_state!(
            es,
            ps,
            merge(saved, (spectral_signature = saved.spectral_signature * "x",)),
        )
        @test_throws ArgumentError H._restore_checkpoint_state!(
            es,
            ps,
            merge(
                saved,
                (
                    coupled_control = merge(
                        saved.coupled_control,
                        (h_max = Q(saved.coupled_control.h_max),),
                    ),
                ),
            ),
        )
        cycle = Any[]
        push!(cycle, cycle)
        cyclic = merge(
            saved,
            (
                coupled_phase = merge(
                    saved.coupled_phase,
                    (report = merge(saved.coupled_phase.report, (intervals = cycle,)),),
                ),
            ),
        )
        @test_throws ArgumentError H._restore_checkpoint_state!(es, ps, cyclic)
        # Both source arrays are numerically zero, so this packet differs only
        # in its original role graph, before any deepcopy could conceal it.
        @test all(iszero, saved.fields.rho_n) && all(iszero, saved.fields.rho_np1)
        bad = merge(saved, (fields = merge(saved.fields, (rho_np1 = saved.fields.rho_n,)),))
        @test_throws ArgumentError H._restore_checkpoint_state!(es, ps, bad)
        livefields =
            es isa EMPIC ?
            merge(saved.fields, (E = (saved.fields.E[1], es.E[2], saved.fields.E[3]),)) :
            merge(saved.fields, (Ey = es.Ey,))
        @test_throws ArgumentError H._restore_checkpoint_state!(
            es,
            ps,
            merge(saved, (fields = livefields,)),
        )
        @test snapshot(es, ps) == before
        @test_throws ArgumentError setfuture(
            es;
            error_control = merge(future, (domain = :single_piece,)),
        )
        @test snapshot(es, ps) == before
        mktempdir() do directory
            raw, archive = joinpath(directory, "coupled.chk"), joinpath(directory, "coupled.run")
            save_checkpoint(raw, es, ps)
            archive_run(archive, es, ps; rng_seed = 0)
            other, rows = fixture(T, Model; count = 2)
            @test other.g !== es.g && other.g.plan !== es.g.plan # Fresh target grid/plans.
            owner = rows[1]
            for row in rows
                for a in (row.x..., row.v..., row.weight)
                    fill!(a, T(NaN))
                end
                row.q, row.m = T(NaN), -one(T)
            end
            for a in (
                other isa EMPIC ?
                (other.E..., other.B..., other.J..., other.Jraw..., other.rho_n, other.rho_np1) :
                (other.Ex, other.Ey, other.Bz, other.Jx, other.Jy, other.rho_n, other.rho_np1)
            )
                fill!(a, T(NaN))
            end
            other.time[], other.step[], other.last_dt[], other.last_sub_dt[] =
                T(Inf), -1, T(Inf), T(Inf)
            load_checkpoint!(other, rows, raw)
            @test rows[1] === owner &&
                  length(rows[1].weight) == 1 &&
                  snapshot(other, rows) == before
            # A positive operation must also pay for historical admission and
            # final verification replay, unlike capture of valid saved history.
            setfuture(es; error_control = control(T))
            setfuture(other; error_control = control(T))
            step_empic!(es, ps, T(1 / 1024))
            step_empic!(other, rows, T(1 / 1024))
            @test snapshot(other, rows) == snapshot(es, ps)
            load_archive!(other, rows, archive)
            @test snapshot(other, rows) == before
        end
        # A resolved positive final edge below h_min occurs after a private
        # accepted interval. The ENTIRE external operation must roll back.
        late, rows = fixture(T, Model; proposal = control(T; hmax = T(1 / 64), hmin = T(1 / 128)))
        initial = snapshot(late, rows)
        @test_throws ArgumentError step_empic!(late, rows, T(5 / 256))
        @test snapshot(late, rows) == initial
        exhausted, rows =
            fixture(T, Model; proposal = control(T; hmax = h, tol = eps(T), trials = 1))
        initial = snapshot(exhausted, rows)
        @test_throws ArgumentError step_empic!(exhausted, rows, h)
        @test snapshot(exhausted, rows) == initial
        setfuture(exhausted; error_control = control(T; hmax = h))
        clean, cleanrows = fixture(T, Model; proposal = control(T; hmax = h))
        step_empic!(exhausted, rows, h)
        step_empic!(clean, cleanrows, h)
        @test physical(exhausted, rows) == physical(clean, cleanrows)
    end
end

@testset "Full-time physical-v relativistic caller and C1 rollback" begin
    for T in (Float32, Float64), Model in (EMPIC, EMPIC1D)
        es, ps = fixture(T, Model; relativistic = true)
        fields = es isa EMPIC ? (es.E..., es.B...) : (es.Ex, es.Ey, es.Bz)
        foreach(a -> fill!(a, zero(T)), fields)
        electric = (T(1 / 8), T(1 / 32), zero(T))
        magnetic = ntuple(_ -> zero(T), 3)
        if es isa EMPIC
            foreach(fill!, es.E, electric)
        else
            fill!(es.Ex, electric[1])
            fill!(es.Ey, electric[2])
        end
        init_empic!(es, ps)
        expected = deepcopy(ps)
        h = T(1 / 64)
        for i = 1:2, p in eachindex(ps[i].weight)
            v = ntuple(c -> ps[i].v[c][p], 3)
            p0 = physical_parts(v, electric, magnetic, ps[i].q, ps[i].m, h, es.c)
            low = radical_reference(v, p0, nothing, T)
            p1 = physical_parts(low, electric, magnetic, ps[i].q, ps[i].m, h, es.c)
            high = radical_reference(v, p0, p1, T)
            for c = 1:3
                expected[i].v[c][p] = high[c]
            end
            for d = 1:length(es.g.n)
                expected[i].x[d][p] =
                    wrapped(T, Q(ps[i].x[d][p]) + Q(h) * (Q(v[d]) + Q(low[d])) / 2, es.g.L[d])
            end
        end
        step_empic!(es, ps, h)
        check_particles(ps, expected, 64eps(T))
        setfuture(es; error_control = control(T; hmax = h, hmin = h))
        es isa EMPIC ? fill!(es.E[1], T(512)) : fill!(es.Ex, T(512))
        init_empic!(es, ps)
        before = snapshot(es, ps)
        @test_throws ArgumentError step_empic!(es, ps, h)
        @test snapshot(es, ps) == before
    end
end

@testset "Exact right edge and non-dyadic source/public clocks" begin
    for T in (Float32, Float64), Model in (EMPIC, EMPIC1D)
        h = T(1 / 64)
        es, ps = fixture(T, Model)
        setfuture(es; species_dt = T[h, 2h, 0, h])
        step_empic!(es, ps, h)
        d = electromagnetic_diagnostics(es, ps)
        @test d.last_deadline[1] == 0 && d.next_deadline[1] == Q(h)
        before = snapshot(es, ps)
        step_empic!(es, ps, zero(T))
        @test snapshot(es, ps) == before
        step_empic!(es, ps, h)
        @test electromagnetic_diagnostics(es, ps).last_deadline[1] == Q(h)
        clock, rows = fixture(
            T,
            Model;
            proposal = control(T; hmax = T(1 / 32), hmin = T(1 / 65536)),
            n_sub = 3,
        )
        step_empic!(clock, rows, T(0.1))
        d = electromagnetic_diagnostics(clock, rows)
        r = d.last_error_control_report
        @test d.source_time == Q(T(0.1)) && d.actual_elapsed == sum(Q(i.dt) for i in r.intervals)
        @test d.public_clock_defect == Q(clock.time[]) - d.source_time
        @test d.public_minus_actual == Q(clock.time[]) - d.actual_elapsed
        @test clock.step[] == 1 && all(i -> i.norm <= 1, r.intervals)
    end
end

@testset "Uninitialized schema4 seeds and ordinary view destinations" begin
    for T in (Float32, Float64), Model in (EMPIC, EMPIC1D)
        es, ps = fixture(T, Model; initialize = false)
        old_id = ps[1].id
        ps[1].id = es.collection.mobile.chunks
        @test_throws ArgumentError H._checkpoint_state(es, ps)
        ps[1].id = old_id
        saved = H._checkpoint_state(es, ps)
        @test !saved.initialized && saved.step == 0 && saved.coupled_phase.report === nothing
        @test all(iszero, values(saved.coupled_phase.costs))
        target, rows = fixture(T, Model; count = 2)
        H._restore_checkpoint_state!(target, rows, saved)
        @test !target.collection.initialized[] && snapshot(target, rows) == snapshot(es, ps)
        bad =
            merge(saved, (coupled_phase = merge(saved.coupled_phase, (initial_energy = one(Q),)),))
        before = snapshot(target, rows)
        @test_throws ArgumentError H._restore_checkpoint_state!(target, rows, bad)
        @test snapshot(target, rows) == before
        init_empic!(target, rows)
        init_empic!(es, ps)
        zero_state = H._checkpoint_state(target, rows)
        zero_before = snapshot(target, rows)
        bad = deepcopy(zero_state)
        bad.fields.rho_n[1] += one(T)
        @test_throws ArgumentError H._restore_checkpoint_state!(target, rows, bad)
        bad = deepcopy(zero_state)
        current = target isa EMPIC ? bad.fields.J[1] : bad.fields.Jx
        current[1] = one(T)
        @test_throws ArgumentError H._restore_checkpoint_state!(target, rows, bad)
        @test snapshot(target, rows) == zero_before
        step_empic!(target, rows, T(1 / 64))
        step_empic!(es, ps, T(1 / 64))
        @test snapshot(target, rows) == snapshot(es, ps)
    end
    for Model in (EMPIC, EMPIC1D)
        original = population(Float64, 1)
        rows = [
            ParticleSet{1,Float64}(
                map(leaf_view, p.x),
                map(leaf_view, p.v),
                leaf_view(p.weight),
                leaf_view(p.id),
                leaf_view(p.tag),
                p.q,
                p.m,
            ) for p in original
        ]
        es = Model(
            FourierGrid((4,), (1.0,)),
            rows;
            mobile = Bool[true, true, false, true],
            species_dt = [1 / 2, 3 / 4, 0, 1 / 2],
            event_budget = 256,
            error_control = control(Float64),
            c = 1.0,
        )
        init_empic!(es, rows)
        step_empic!(es, rows, 1 / 64)
        saved = H._checkpoint_state(es, rows)
        before = snapshot(es, rows)
        @test all(p -> p.x isa Tuple{Vector{Float64}}, saved.species)
        leaves = [(p.x..., p.v..., p.weight, p.id, p.tag) for p in rows]
        step_empic!(es, rows, 1 / 64)
        H._restore_checkpoint_state!(es, rows, saved)
        @test snapshot(es, rows) == before
        @test all(
            pair -> all(
                p -> p[1] === p[2],
                zip(pair[1], (pair[2].x..., pair[2].v..., pair[2].weight, pair[2].id, pair[2].tag)),
            ),
            zip(leaves, rows),
        )
    end
end
end
