using HybridPlasmaPIC, Test, LinearAlgebra

# Every public field array, in a fixed order, for bitwise state comparisons.
function tr_fields(f)
    return (
        f.n,
        f.ui...,
        f.B...,
        f.E...,
        f.J...,
        f.pe,
        f.gradp...,
        f.ninv,
        f.lapJ...,
        f.pforce...,
        f.cgl_dp,
    )
end

tr_state(st) = (tr_fields(st.fields)..., st.pressure_work.Eforce...)
tr_clock(st) = (st.time[], st.step[], st.last_dt[], st.fields.floor_count[])
tr_populations(ps) = ps isa ParticleSet ? (ps,) : ps
tr_particles(ps) = [(map(copy, p.x), map(copy, p.v)) for p in tr_populations(ps)]
tr_same_particles(ps, saved) =
    all(isequal(p.x, s[1]) && isequal(p.v, s[2]) for (p, s) in zip(tr_populations(ps), saved))

# A periodic plasma with flow, magnetic shear and nonuniform pressure. The
# density stays near one, far above the default floor.
function tr_fixture(
    ::Type{T},
    ::Val{D},
    closure;
    η = 0.0,
    collection = false,
    shape = CIC(),
) where {T,D}
    g = FourierGrid(ntuple(d -> d == 1 ? 8 : 4 + d, D), ntuple(_ -> T(2pi), D))
    counts = ntuple(d -> 2g.n[d], D)
    ps = ParticleSet{D,T}(prod(counts))
    load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
    set_density_weight!(ps, one(T), g)
    @. ps.v[1] = T(0.02) * sin(ps.x[1])
    @. ps.v[2] = T(0.03) * cos(ps.x[1])
    fill!(ps.v[3], T(0.01))
    particles = collection ? [ps, ParticleSet{D,T}(0; q = T(2), m = T(4))] : ps
    model = closure isa TransportElectrons ? HybridModel(closure; η) : HybridModel(closure)
    st = HybridStepper(g, model, shape, 0)
    fill!(st.fields.B[1], T(0.5))
    for I in CartesianIndices(g.n)
        x = T(I[1] - 1) * g.dx[1]
        y = D >= 2 ? T(I[2] - 1) * g.dx[2] : zero(T)
        st.fields.B[3][I] = T(0.2) * sin(x) + T(0.05) * cos(y)
        st.fields.pe[I] = T(0.4) * (one(T) + T(0.1) * sin(x) + T(0.05) * cos(y))
    end
    init!(st, particles)
    return st, particles
end

@testset "Transport electrons configuration and validation" begin
    for T in (Float32, Float64)
        c = TransportElectrons(T(5 / 3))
        @test closure_gamma(c) == T(5 / 3)
        for bad in (zero(T), one(T), -one(T), T(Inf), T(NaN))
            @test_throws ArgumentError TransportElectrons(bad)
            @test_throws ArgumentError TransportElectrons{T}(bad)
        end
        @test HybridModel(c; η = 0.02).η == 0.02
        for keywords in ((ηH = 0.01,), (de2 = 0.01,), (η = -0.1,), (η = Inf,), (η = NaN,))
            @test_throws ArgumentError HybridModel(c; keywords...)
        end
        # The ideal closure keeps its contract.
        @test_throws ArgumentError HybridModel(EvolvingElectrons(T(5 / 3)); η = 0.01)
        g = FourierGrid((8,), (T(2pi),))
        @test_throws ArgumentError CAMCLStepper(g, HybridModel(c), CIC(), 0)
        @test_throws ArgumentError HallMHDModel(c)
        @test_throws ArgumentError electron_pressure!(zeros(T, 8), ones(T, 8), c)
        @test_throws ArgumentError electron_internal_energy(ones(T, 8), c, g)
        # Only the transport closure owns the six carried-force arrays.
        transport = HybridStepper(g, HybridModel(c; η = 0.02), CIC(), 0)
        ideal = HybridStepper(g, HybridModel(EvolvingElectrons(T(5 / 3))), CIC(), 0)
        w = transport.pressure_work
        @test w.Eforce isa NTuple{3,Vector{T}} && w.saved_Eforce isa NTuple{3,Vector{T}}
        @test all(size(a) == g.n for a in (w.Eforce..., w.saved_Eforce...))
        @test ideal.pressure_work.Eforce === nothing && ideal.pressure_work.saved_Eforce === nothing
        @test HybridStepper(g, HybridModel(IsothermalElectrons(T(0.1))), CIC(), 0).pressure_work ===
              nothing
        added = Base.summarysize(transport) - Base.summarysize(ideal)
        @test 6 * prod(g.n) * sizeof(T) <= added <= 6 * prod(g.n) * sizeof(T) + 1024
        # Positional models are revalidated at the Ohm and stepping boundaries.
        raw = HybridModel(c, 0.01, 0.01, 1e-6, 0.0)
        f = HybridFields{1,T}(g.n)
        fill!(f.n, one(T))
        snapshot = map(copy, tr_fields(f))
        @test_throws ArgumentError ohms_law!(f, raw, g)
        @test all(isequal(a, b) for (a, b) in zip(tr_fields(f), snapshot))
        transport.model = raw
        @test_throws ArgumentError step!(transport, ParticleSet{1,T}(0), T(0.01))
        # Direct Ohm calls return the carried full field and check density first.
        st, ps = tr_fixture(T, Val(2), c; η = 0.04)
        f = deepcopy(st.fields)
        @test ohms_law!(f, st.model, st.g) === f.E
        @test all(isequal(a, b) for (a, b) in zip(f.E, st.fields.E))
        force = st.pressure_work.Eforce
        @test all(isequal(st.fields.E[k], force[k] .+ T(0.04) .* st.fields.J[k]) for k = 1:3)
        f.n[2] = prevfloat(T(st.model.nfloor))
        snapshot = map(copy, tr_fields(f))
        @test_throws ArgumentError ohms_law!(f, st.model, st.g)
        @test all(isequal(a, b) for (a, b) in zip(tr_fields(f), snapshot))
    end
    # A nonzero coefficient must survive conversion to the grid precision.
    g32 = FourierGrid((8,), (Float32(2pi),))
    for η in (1e-50, 1e300)
        @test_throws ArgumentError HybridStepper(
            g32,
            HybridModel(TransportElectrons(); η),
            CIC(),
            0,
        )
    end
    @test HybridStepper(g32, HybridModel(TransportElectrons(); η = 1e-40), CIC(), 0) isa
          HybridStepper
    @test HybridStepper(g32, HybridModel(TransportElectrons(); η = -0.0), CIC(), 0) isa
          HybridStepper
end

# Direct odd-grid discrete Fourier differentiation, independent of package FFTs.
function tr_dft_derivative(n)
    x = 2pi .* (0:n-1) ./ n
    modes = -(n ÷ 2):(n÷2)
    return [real(sum(im * k * cis(k * (xi - xj)) for k in modes)) / n for xi in x, xj in x]
end

# Independent semidiscrete transport rates on a 2D grid: Faraday with the full
# field E_force + ηJ and the pressure equation with Ohmic heat (γ−1)η|J|².
function tr_reference_rhs(z, n, ui, γ, η, Dx, Dy)
    dx(a) = Dx * a
    dy(a) = a * transpose(Dy)
    bx, by, bz, p = eachslice(z; dims = 3)
    j = (dy(bz), -dx(bz), dx(by) - dy(bx))
    u = ntuple(c -> ui[c] .- j[c] ./ n, 3)
    e = (
        -(u[2] .* bz - u[3] .* by) .- dx(p) ./ n .+ η .* j[1],
        -(u[3] .* bx - u[1] .* bz) .- dy(p) ./ n .+ η .* j[2],
        -(u[1] .* by - u[2] .* bx) .+ η .* j[3],
    )
    heat = (γ - 1) .* η .* (j[1] .^ 2 .+ j[2] .^ 2 .+ j[3] .^ 2)
    rate = -u[1] .* dx(p) .- u[2] .* dy(p) .- γ .* p .* (dx(u[1]) + dy(u[2])) .+ heat
    return cat(-dy(e[3]), dx(e[3]), -dx(e[2]) + dy(e[1]), rate; dims = 3), j
end

function tr_stage_rates!(st, B, pe)
    HybridPlasmaPIC._bpressure_rhs!(st.k1, st.pressure_work.kp[1], B, pe, st, 1, 1)
    return cat(map(copy, st.k1)..., copy(st.pressure_work.kp[1]); dims = 3)
end

@testset "Transport stage rates vs independent discrete DFT" begin
    γ, η = 5 / 3, 0.07
    g = FourierGrid((7, 5), (2pi, 2pi))
    st = HybridStepper(g, HybridModel(TransportElectrons(γ); η), CIC(), 0)
    x = reshape((0:6) .* g.dx[1], :, 1)
    y = reshape((0:4) .* g.dx[2], 1, :)
    n = @. 1 + 0.1 * cos(x) * cos(y)
    ui = (repeat(0.12 .* sin.(y), 7, 1), repeat(0.08 .* cos.(x), 1, 5), @. 0.04 * sin(x) * cos(y))
    z0 = cat(
        repeat(0.2 .+ 0.03 .* cos.(y), 7, 1),
        repeat(0.1 .+ 0.02 .* sin.(x), 1, 5),
        (@. 0.4 + 0.1 * sin(x) + 0.05 * cos(y)),
        (@. 0.7 + 0.04 * cos(x) + 0.03 * sin(y));
        dims = 3,
    )
    Dx, Dy = tr_dft_derivative(7), tr_dft_derivative(5)
    reference, j = tr_reference_rhs(z0, n, ui, γ, η, Dx, Dy)
    st.fn .= n
    for c = 1:3
        st.fui[c] .= ui[c]
        st.fields.B[c] .= z0[:, :, c]
    end
    st.fields.pe .= z0[:, :, 4]
    HybridPlasmaPIC._ohm_ninv!(st.fields.ninv, st.fn, st.model.nfloor, st.fields.floor_count)
    production = tr_stage_rates!(st, st.fields.B, st.fields.pe)
    @test maximum(abs, production .- reference) < 128eps(Float64)
    @test maximum(c -> maximum(abs, st.fields.J[c] .- j[c]), 1:3) < 64eps(Float64)
    # The η = 0 stage is the ideal stage, bitwise.
    st.model = HybridModel(TransportElectrons(γ); η = 0.0)
    ideal_transport = tr_stage_rates!(st, st.fields.B, st.fields.pe)
    ideal = HybridStepper(g, HybridModel(EvolvingElectrons(γ)), CIC(), 0)
    ideal.fn .= st.fn
    foreach(copyto!, ideal.fui, st.fui)
    copyto!(ideal.fields.ninv, st.fields.ninv)
    @test isequal(ideal_transport, tr_stage_rates!(ideal, st.fields.B, st.fields.pe))
    # Periodic semidiscrete Ohmic exchange: the magnetic loss from ηJ in Faraday
    # equals the internal-energy gain from the heat, both equal to ηΣ|J|²ΔV.
    dV = prod(g.dx)
    difference = production .- ideal_transport
    magnetic = sum(sum(z0[:, :, c] .* difference[:, :, c]) for c = 1:3) * dV
    internal = sum(difference[:, :, 4]) / (γ - 1) * dV
    ohmic = η * sum(j[1] .^ 2 .+ j[2] .^ 2 .+ j[3] .^ 2) * dV
    @test ohmic > 1e-3
    # Differences of O(1) stage rates carry absolute roundoff of a few eps.
    @test abs(internal - ohmic) < 1e-12
    @test abs(magnetic + ohmic) < 1e-12
end

# Circular Beltrami mode B = A(e1 cos φ + e2 sin φ), φ = k⋅x, with J = −|k|B:
# J×B = 0, |J|² is uniform and uniform pressure has no gradient. For constant
# density and ions at rest, B decays as exp(−ηk²t) and pe gains the field loss.
function tr_beltrami(::Type{T}, ::Val{D}, η, dt, steps; NB = 1) where {T,D}
    n, mode = D == 1 ? ((16,), (1,)) : D == 2 ? ((8, 8), (1, 1)) : ((6, 6, 6), (1, 1, 1))
    g = FourierGrid(n, ntuple(_ -> T(2pi), D))
    k = T[mode..., zeros(Int, 3 - D)...]
    kk = norm(k)
    khat = k ./ kk
    e1 = normalize(cross(khat, abs(khat[3]) < T(0.9) ? T[0, 0, 1] : T[1, 0, 0]))
    e2 = cross(khat, e1)
    counts = ntuple(d -> 2n[d], D)
    ps = ParticleSet{D,T}(prod(counts))
    load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
    set_density_weight!(ps, one(T), g)
    st = HybridStepper(g, HybridModel(TransportElectrons(T(5 / 3)); η), CIC(), 0)
    amplitude, pressure = T(0.3), T(0.5)
    phase = map(I -> sum(T(mode[d]) * T(I[d] - 1) * g.dx[d] for d = 1:D), CartesianIndices(n))
    for c = 1:3
        @. st.fields.B[c] = amplitude * (e1[c] * cos(phase) + e2[c] * sin(phase))
    end
    fill!(st.fields.pe, pressure)
    init!(st, ps)
    total(st) =
        energy_budget(ps, st.fields.B, st.fields.n, st.model.closure, g; pe = st.fields.pe).total
    initial = total(st)
    for _ = 1:steps
        step!(st, ps, dt; NB)
    end
    t = st.time[]
    rate = T(η) * kk^2
    A = amplitude * exp(-rate * t)
    expected_pe = pressure + (T(5 / 3) - 1) * amplitude^2 * (1 - exp(-2rate * t)) / 2
    field_error = maximum(
        maximum(abs, st.fields.B[c] .- A .* (e1[c] .* cos.(phase) .+ e2[c] .* sin.(phase))) for
        c = 1:3
    )
    return (
        field = field_error,
        pressure = maximum(abs, st.fields.pe .- expected_pe),
        heated = expected_pe - pressure,
        energy = abs(total(st) - initial),
        total = initial,
        speed = maximum(maximum(abs, v) for v in ps.v),
        force = maximum(maximum(abs, a) for a in st.pressure_work.Eforce),
        split = all(
            isequal(st.fields.E[c], st.pressure_work.Eforce[c] .+ T(η) .* st.fields.J[c]) for
            c = 1:3
        ),
        field_scale = maximum(maximum(abs, a) for a in st.fields.E),
    )
end

@testset "Transport Beltrami Ohmic heating through step!" begin
    for T in (Float32, Float64), D = 1:3, NB in (1, 2)
        # |k|² = D, so ηk² = 0.2 and RK4 truncation stays below 1e-11.
        r = tr_beltrami(T, Val(D), 0.2 / D, T(0.05), 8; NB)
        tolerance = T === Float32 ? T(64) * eps(T) : T(1e-10)
        @test r.field < tolerance * T(0.3)
        @test r.pressure < tolerance
        @test r.heated > T(0.004)
        # Ions see only E_force ≈ 0; a full-E push would accelerate them by ηJ.
        @test r.speed < T(16) * eps(T)
        @test r.force < T(64) * eps(T)
        @test r.field_scale > T(0.02)
        @test r.split
        @test r.energy < tolerance * r.total
    end
    # Fixed-moment RK4: field, heat and the total energy close at fourth order.
    results = [tr_beltrami(Float64, Val(2), 0.4, 0.4 / k, k) for k in (2, 4, 8)]
    for name in (:field, :pressure, :energy)
        errors = [getproperty(r, name) for r in results]
        @test all(ratio -> 12 < ratio < 20, errors[1:2] ./ errors[2:3])
    end
    @info "Beltrami RK4 energy errors" energy = [r.energy for r in results]
end

@testset "Transport η = 0 reduces exactly to EvolvingElectrons" begin
    for T in (Float32, Float64), D = 1:3, collection in (false, true)
        a, pa = tr_fixture(T, Val(D), TransportElectrons(T(5 / 3)); collection)
        b, pb = tr_fixture(T, Val(D), EvolvingElectrons(T(5 / 3)); collection)
        @test all(isequal(x, y) for (x, y) in zip(tr_fields(a.fields), tr_fields(b.fields)))
        @test all(isequal(x, y) for (x, y) in zip(a.pressure_work.Eforce, a.fields.E))
        for (dt, nb) in ((T(0.01), 1), (T(0.013), 2), (T(0.007), 4))
            step!(a, pa, dt; NB = nb)
            step!(b, pb, dt; NB = nb)
            @test all(isequal(x, y) for (x, y) in zip(tr_fields(a.fields), tr_fields(b.fields)))
            @test all(isequal(x, y) for (x, y) in zip(a.pressure_work.Eforce, a.fields.E))
            @test tr_same_particles(pa, tr_particles(pb))
            @test tr_clock(a) == tr_clock(b)
        end
    end
end

# NGP particles at grid nodes move less than half a cell, so the deposited
# density stays exactly uniform while the ion current is correlated with J.
function tr_energy_error(steps; multispecies = false, variable = false, η = 0.05)
    N = 33
    g = FourierGrid((N,), (2pi,))
    x = (0:N-1) .* g.dx[1]
    qs, ms, ns = multispecies ? ((1.0, 2.0), (1.0, 4.0), (0.6, 0.2)) : ((1.0,), (1.0,), (1.0,))
    populations = ParticleSet{1,Float64}[]
    for s in eachindex(qs)
        ps = ParticleSet{1,Float64}(N; q = qs[s], m = ms[s])
        ps.x[1] .= x
        set_density_weight!(ps, ns[s], g)
        @. ps.v[1] = 0.02 * cos(x)
        @. ps.v[2] = (0.03 / s) * cos(x)
        push!(populations, ps)
    end
    st = HybridStepper(g, HybridModel(TransportElectrons(); η), NGP(), N * length(qs))
    fill!(st.fields.B[1], 0.5)
    @. st.fields.B[3] = 0.2 * sin(x)
    @. st.fields.pe = 0.7 + 0.04 * sin(x)
    ps = multispecies ? populations : only(populations)
    init!(st, ps)
    closure = st.model.closure
    initial = energy_budget(ps, st.fields.B, st.fields.n, closure, g; pe = st.fields.pe).total
    duration = 0.4
    for k = 1:steps
        dt = duration / steps * (variable ? (isodd(k) ? 0.7 : 1.3) : 1.0)
        step!(st, ps, dt; NB = 2)
    end
    # Common-time kinetic energy: interpolate the staggered velocities to the
    # integer time of B and pe with one further step.
    B = map(copy, st.fields.B)
    p = copy(st.fields.pe)
    previous = [map(copy, q.v) for q in populations]
    prevdt = st.last_dt[]
    nextdt = duration / steps * (variable ? 0.7 : 1.0)
    step!(st, ps, nextdt; NB = 2)
    kinetic = sum(eachindex(populations)) do s
        q = populations[s]
        sum(1:3) do c
            vc = (nextdt .* previous[s][c] .+ prevdt .* q.v[c]) ./ (prevdt + nextdt)
            q.m * sum(q.weight .* vc .^ 2) / 2
        end
    end
    total = kinetic + magnetic_energy(B, g) + sum(p) / (5 / 3 - 1) * prod(g.dx)
    shift =
        maximum(maximum(min.(abs.(q.x[1] .- x), g.L[1] .- abs.(q.x[1] .- x))) for q in populations)
    @test shift < g.dx[1] / 2
    @test maximum(abs, st.fields.n .- sum(qs .* ns)) < 8eps(Float64)
    return abs(total - initial)
end

@testset "Transport common-time total energy converges at second order" begin
    # Predeclared windows. Fixed steps are asymptotic from 8 steps. The
    # alternating 0.7/1.3 sequence partly cancels the leading error, so its
    # ratios approach 4 from above and enter 3–5 only from 64 steps.
    for multispecies in (false, true), variable in (false, true)
        counts = variable ? (64, 128, 256) : (8, 16, 32)
        errors = [tr_energy_error(k; multispecies, variable) for k in counts]
        bounds = variable ? (3.0, 5.0) : (3.4, 4.6)
        @test all(r -> bounds[1] < r < bounds[2], errors[1:2] ./ errors[2:3])
        # Heat without ion recoil would leave |∫∫ηJi⋅J dV dt| ≈ 4e-4 here.
        @test errors[end] < (variable ? 1e-9 : 1e-7)
        @info "transport common-time energy error" multispecies variable counts errors
    end
end

# NGP lattice with one particle at each node; `index` moves by `fraction` of a
# cell in each half step. Cold electrons and B = 0 keep the force and every
# pressure rate exactly zero, so only the density contract can reject a step.
function tr_density_fixture(::Type{T}, index, fraction; nfloor = 0.5) where {T}
    g = FourierGrid((8,), (T(2pi),))
    ps = ParticleSet{1,T}(8)
    ps.x[1] .= T.(0:7) .* g.dx[1]
    set_density_weight!(ps, one(T), g)
    dt = T(0.1)
    ps.v[1][index] = fraction * g.dx[1] / (dt / 2)
    st = HybridStepper(g, HybridModel(TransportElectrons(T(5 / 3)); η = 0.05, nfloor), NGP(), 8)
    init!(st, ps)
    return st, ps, dt
end

function tr_rejected(f, st, ps, pattern)
    fields = map(copy, tr_state(st))
    identities = tr_state(st)
    particles = tr_particles(ps)
    clock = tr_clock(st)
    err = try
        f()
        nothing
    catch error
        error
    end
    @test err isa ArgumentError
    @test occursin(pattern, sprint(showerror, err))
    @test all(isequal(a, b) for (a, b) in zip(tr_state(st), fields))
    @test all(a === b for (a, b) in zip(tr_state(st), identities))
    @test tr_same_particles(ps, particles)
    @test tr_clock(st) == clock
    return err
end

@testset "Transport density admissibility and rollback" begin
    for T in (Float32, Float64)
        # The deposited density can equal nfloor; just below it is rejected.
        st, ps, _ = tr_density_fixture(T, 1, T(0.1))
        level = minimum(st.fields.n)
        exact =
            HybridStepper(st.g, HybridModel(TransportElectrons(T(5 / 3)); nfloor = level), NGP(), 8)
        fill!(exact.fields.pe, T(0.3))
        @test init!(exact, ps) === exact
        above = Float64(nextfloat(level))
        strict =
            HybridStepper(st.g, HybridModel(TransportElectrons(T(5 / 3)); nfloor = above), NGP(), 8)
        fill!(strict.fields.pe, T(0.3))
        tr_rejected(() -> init!(strict, ps), strict, ps, "density")
        # Midpoint, endpoint/predictor and empty-population violations roll back.
        st, ps, dt = tr_density_fixture(T, 3, T(0.6))
        tr_rejected(() -> step!(st, ps, dt), st, ps, "midpoint density")
        st, ps, dt = tr_density_fixture(T, 3, T(0.4))
        tr_rejected(() -> step!(st, ps, dt), st, ps, "Ohm-law density")
        retry, retryps, _ = tr_density_fixture(T, 3, T(0.4))
        step!(st, ps, dt / 4)
        step!(retry, retryps, dt / 4)
        @test all(isequal(a, b) for (a, b) in zip(tr_state(st), tr_state(retry)))
        @test tr_same_particles(ps, tr_particles(retryps))
        empty = HybridStepper(st.g, HybridModel(TransportElectrons(T(5 / 3))), CIC(), 0)
        for nothing_loaded in (ParticleSet{1,T}(0), ParticleSet{1,T}[])
            tr_rejected(() -> init!(empty, nothing_loaded), empty, nothing_loaded, "density")
        end
        # Scalar charge must be one; the one-element collection is canonical.
        charged = ParticleSet{1,T}(16; q = T(2), m = T(4))
        load_lattice_1d!(charged, zero(T), st.g.L[1])
        set_density_weight!(charged, one(T), st.g)
        target = HybridStepper(st.g, HybridModel(TransportElectrons(T(5 / 3))), CIC(), 16)
        tr_rejected(() -> init!(target, charged), target, charged, "q = 1")
        @test init!(target, [charged]) === target
        @test step!(target, [charged], T(0.01)) === target
    end
end

@testset "Transport pressure, force and storage rejection" begin
    for T in (Float32, Float64), collection in (false, true)
        # Compression empties pressure at a later subcycle; heating cannot hide it.
        st, ps = tr_fixture(T, Val(1), TransportElectrons(T(3)); η = 0.02, collection)
        for p in tr_populations(ps)
            @. p.v[1] = sin(p.x[1])
        end
        init!(st, ps)
        err = tr_rejected(() -> step!(st, ps, T(0.8); NB = 2), st, ps, "subcycle")
        control, controlps = deepcopy((st, ps))
        step!(st, ps, T(0.001); NB = 2)
        step!(control, controlps, T(0.001); NB = 2)
        @test all(isequal(a, b) for (a, b) in zip(tr_state(st), tr_state(control)))
        # A nonfinite carried force is rejected before mutation; zero dt is a no-op.
        for bad in (T(NaN), T(Inf))
            saved = st.pressure_work.Eforce[2][1]
            st.pressure_work.Eforce[2][1] = bad
            state = map(copy, tr_state(st))
            @test step!(st, ps, zero(T)) === st
            @test all(isequal(a, b) for (a, b) in zip(tr_state(st), state))
            tr_rejected(() -> step!(st, ps, T(0.001)), st, ps, "ion force")
            st.pressure_work.Eforce[2][1] = saved
        end
        # The carried force needs independent storage.
        force = st.pressure_work.Eforce
        for alias in (
            (st.fields.E[1], force[2], force[3]),
            (force[1], force[1], force[3]),
            (force[1], force[2], st.pressure_work.saved_Eforce[3]),
            (force[1], force[2], st.fields.pe),
        )
            st.pressure_work.Eforce = alias
            tr_rejected(() -> step!(st, ps, T(0.001)), st, ps, "independent storage")
            tr_rejected(() -> init!(st, ps), st, ps, "independent storage")
        end
        st.pressure_work.Eforce = force
        @test step!(st, ps, T(0.001)) === st
    end
    # Overflowing Ohmic heat is an explicit rejection, never a clipped value.
    st, ps = tr_fixture(Float32, Val(1), TransportElectrons(5.0f0 / 3); η = 1e30)
    st.fields.B[3] .*= 1.0f6
    init!(st, ps)
    tr_rejected(() -> step!(st, ps, 1.0f-3), st, ps, "rate must be finite")
end

function tr_heating_reference(γ, η, j)
    Q = Rational{BigInt}
    return (Q(γ) - 1) * Q(η) * sum(Q(c)^2 for c in j)
end

# Normwise contract: within four ulps of the exact value plus two subnormal
# quanta, or an explicit non-finite value when the exact result overflows.
function tr_heating_admissible(value::T, exact) where {T}
    Q = Rational{BigInt}
    tolerance = Q(4eps(T)) * exact + 2Q(nextfloat(zero(T)))
    isfinite(value) && return abs(Q(value) - exact) <= tolerance
    return isinf(value) && exact > Q(floatmax(T)) * (1 - Q(8eps(T)))
end

@testset "Transport Ohmic heating range vs exact rational" begin
    for T in (Float32, Float64)
        tiny = nextfloat(zero(T))
        small = sqrt(floatmin(T))
        large = sqrt(floatmax(T))
        currents = (
            (zero(T), zero(T), zero(T)),
            (one(T), -T(0.5), T(0.25)),
            (zero(T), small, -small),
            (tiny, -tiny, tiny),
            (small / 4, small / 8, zero(T)),
            (large, -large / 3, zero(T)),
            (T(2) * large, zero(T), zero(T)),
            (floatmax(T), floatmax(T), floatmax(T)),
            (T(1e-3), floatmin(T), tiny),
        )
        coefficients = (tiny, floatmin(T), T(1e-20), T(0.07), one(T), T(1e20), floatmax(T))
        gammas = (nextfloat(one(T)), T(5 / 3), T(2), T(3), T(100), floatmax(T))
        for γ in gammas, η in coefficients
            # Added to an ordinary pressure rate, and alone across the range.
            for base in (T(0.5), zero(T))
                rate = fill(base, length(currents))
                J = ntuple(c -> T[j[c] for j in currents], 3)
                HybridPlasmaPIC._add_ohmic_pressure_heating!(rate, J, γ, η)
                for (i, j) in enumerate(currents)
                    exact = tr_heating_reference(γ, η, j) + Rational{BigInt}(base)
                    @test tr_heating_admissible(rate[i], exact)
                end
            end
        end
        # Normal results whose ordinary squares overflow or underflow entirely.
        exponents = T === Float32 ? (-80, 60) : (-600, 300)
        for (j, η) in (
            ((floatmax(T), -floatmax(T), floatmax(T)), tiny),
            ((ldexp(one(T), exponents[1]), zero(T), zero(T)), ldexp(one(T), exponents[2])),
        )
            heat = zeros(T, 1)
            HybridPlasmaPIC._add_ohmic_pressure_heating!(heat, map(c -> T[c], j), T(5 / 3), η)
            exact = tr_heating_reference(T(5 / 3), η, j)
            @test Rational{BigInt}(floatmin(T)) < exact < Rational{BigInt}(floatmax(T))
            @test tr_heating_admissible(heat[1], exact)
            @test heat[1] >= floatmin(T)
        end
    end
end
