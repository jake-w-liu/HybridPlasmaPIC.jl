module ElectronIonThermalExchangeTests
using HybridPlasmaPIC, Test, Random, Serialization, LinearAlgebra
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

function fixture(
    ::Type{T};
    D = 1,
    ns = 2,
    cold = false,
    sparse = false,
    enabled = true,
    logs = T(1 / 4),
    dt = T(1 / 128),
    pressure = T(4),
    bulk = zero(T),
    uniform = false,
    epsilon = one(T),
    number_scale = one(T),
    R = Xoshiro,
    electron_mass = T(1 / 16),
    closure = EvolvingElectrons(T(5 / 3)),
    cell_width = one(T),
) where {T}
    dims = ntuple(d -> d == 1 ? 4 : 1, D)
    g = FourierGrid(dims, ntuple(d -> d == 1 ? T(4) * cell_width : one(T), D))
    species = [ParticleSet{D,T}(12; q = T(a), m = T(a + 1)) for a = 1:ns]
    for (a, p) in enumerate(species), node = 1:4, j = 1:3
        i = 3(node - 1) + j
        # NGP ownership is node+1 (periodically), unlike a floor-cell bin.
        for d = 1:D
            p.x[d][i] = (d == 1 ? T(node - 1) + T(3 / 4) : T(1 / 4)) * g.dx[d]
        end
        p.weight[i] = j == 1 ? one(T) : j == 2 && !(sparse && node == 3) ? T(3) : zero(T)
        # Weights 1:3 and deviations 3:-1 give exactly zero warm bulk.
        amplitude = T(a * (uniform ? 1 : 1 + node / 8) / 8)
        for c = 1:3
            direction = T((1, -1 / 2, 1 / 4)[c])
            p.v[c][i] =
                j == 3 ? zero(T) :
                bulk + (cold ? zero(T) : (j == 1 ? T(3) : -one(T)) * amplitude * direction)
        end
    end
    inner = HybridStepper(g, HybridModel(closure), NGP(), sum(nparticles, species))
    fill!(inner.fields.pe, pressure)
    request = enabled ? ElectronIonThermalExchange(; electron_mass, coulomb_log = logs) : nothing
    st = CollisionalHybridStepper(
        inner,
        species;
        dt,
        rng = R(913),
        coulomb_log = zero(T),
        epsilon0 = epsilon,
        weight_to_number = number_scale,
        electron_ion = request,
    )
    init!(st, species)
    return st, species
end

function bytes(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end
function snapshot(st, species)
    s = H._checkpoint_state(st, species)
    return deepcopy(merge(s, (rng = bytes(s.rng),)))
end
function nodeof(p, g, i)
    I = ntuple(d -> mod(floor(Int, p.x[d][i] / g.dx[d] + 1 / 2), g.n[d]) + 1, length(g.n))
    return LinearIndices(g.n)[I...]
end
function moments(p, g, cell)
    ids = [i for i in eachindex(p.weight) if p.weight[i] > 0 && nodeof(p, g, i) == cell]
    W, m = sum((Q(p.weight[i]) for i in ids); init = zero(Q)), Q(p.m)
    P = ntuple(c -> sum((m * Q(p.weight[i]) * Q(p.v[c][i]) for i in ids); init = zero(Q)), 3)
    mean = ntuple(c -> P[c] / (m * W), 3)
    K = sum(
        (m * Q(p.weight[i]) * sum(Q(p.v[c][i])^2 for c = 1:3) / 2 for i in ids);
        init = zero(Q),
    )
    U = sum(
        (m * Q(p.weight[i]) * sum((Q(p.v[c][i]) - mean[c])^2 for c = 1:3) / 2 for i in ids);
        init = zero(Q),
    )
    return (
        weight = W,
        mean = mean,
        kinetic = K,
        momentum = P,
        thermal = U,
        capacity = 3W / 2,
        count = length(ids),
    )
end

# Independent physical rate and a dense star solve, rather than production's
# beta elimination, exact-root normalization, or positive O(species) formula.
function star_oracle(st, cell, old)
    return setprecision(BigFloat, 256) do
        V = prod(BigFloat, st.grid.dx)
        ne, pe = BigFloat(cell.density), BigFloat(cell.pressure_before)
        Ce, Te = ne * V / (BigFloat(st.inner.model.closure.γ) - 1), pe / ne
        me, h = BigFloat(st.electron_ion.electron_mass), BigFloat(st.dt)
        n = length(old)
        A = zeros(BigFloat, n + 1, n + 1)
        rhs = zeros(BigFloat, n + 1)
        A[1, 1], rhs[1] = Ce, Ce * Te
        rates = BigFloat[]
        for a = 1:n
            s, p = old[a], st.species[a]
            C, Ti = BigFloat(s.capacity), BigFloat(s.thermal) / BigFloat(s.capacity)
            rate =
                8sqrt(2BigFloat(pi)) *
                BigFloat(st.weight_to_number) *
                (BigFloat(s.weight) / V) *
                BigFloat(p.q)^2 *
                BigFloat(st.electron_ion.logarithm[a]) / (
                    3 *
                    (4BigFloat(pi) * BigFloat(st.epsilon0))^2 *
                    me *
                    BigFloat(p.m) *
                    (Te / me + Ti / BigFloat(p.m))^(BigFloat(3) / 2)
                )
            push!(rates, rate)
            conductance = h * Ce * rate
            A[1, 1] += conductance
            A[1, a+1] = -conductance
            A[a+1, 1], A[a+1, a+1], rhs[a+1] = -conductance, C + conductance, C * Ti
        end
        return A \ rhs, rates, Ce, Te
    end
end

function check_receipts(st, species, before)
    r = electron_ion_exchange_diagnostics(st, species)
    V, gamma = prod(Q, st.grid.dx), Q(st.inner.model.closure.γ) - 1
    expected = zeros(eltype(st.inner.fields.pe), st.grid.n)
    partial = similar(expected)
    for p in species
        density!(partial, p, st.grid, NGP())
        expected .+= p.q .* partial
    end
    @test st.inner.fields.n == expected
    for cell in r.cells
        @test cell.density == expected[cell.cell]
        @test cell.pressure_before == before.pe[cell.cell]
        for (a, row) in enumerate(cell.ions)
            @test row.before == moments(before.species[a], st.grid, cell.cell)
            @test row.after == moments(species[a], st.grid, cell.cell)
            @test all(
                c ->
                    abs(row.after.momentum[c] - row.before.momentum[c]) <=
                    row.materialized.momentum_bound[c],
                1:3,
            )
            @test abs(row.after.kinetic - row.materialized.frame_energy) <=
                  row.materialized.kinetic_bound
            u = BigInt(1) // (BigInt(1) << precision(eltype(st.inner.fields.pe)))
            @test abs(row.materialized.scale^2 * row.materialized.seed_energy - row.target) <=
                  (2u + u^2) * row.target
        end
        actual = sum(row.after.kinetic - row.before.kinetic for row in cell.ions)
        expected_pressure = setprecision(BigFloat, 512) do
            typeof(cell.pressure_after)(BigFloat(Q(cell.pressure_before) - gamma * actual / V))
        end
        @test isequal(cell.pressure_after, expected_pressure)
        defect = actual + V * (Q(cell.pressure_after) - Q(cell.pressure_before)) / gamma
        @test cell.energy_defect == defect && abs(defect) <= cell.energy_bound
        @test cell.charge_defect ==
              Q(cell.density) * V -
              sum(Q(p.q) * cell.ions[a].before.weight for (a, p) in enumerate(species))
        @test iszero(cell.charge_defect)
        @test cell.entropy in (:resolved, :cold_positive, :unchanged)
        @test 0 <= cell.entropy_bounds[1] <= cell.entropy_bounds[2]
        if cell.entropy === :resolved
            setprecision(BigFloat, 256) do
                Ce = BigFloat(Q(cell.density) * V / gamma)
                entropy = Ce * log(BigFloat(cell.pressure_after) / BigFloat(cell.pressure_before))
                for row in cell.ions
                    entropy +=
                        BigFloat(row.before.capacity) *
                        log(BigFloat(row.after.thermal) / BigFloat(row.before.thermal))
                end
                @test BigFloat(cell.entropy_bounds[1]) <=
                      entropy <=
                      BigFloat(cell.entropy_bounds[2])
            end
        end
    end
    return r
end

@testset "Independent thermal rate, star targets, native payment and pressure feedback" begin
    for T in (Float32, Float64), D = 1:3, ns = 1:2
        (T === Float32 && D > 1) && continue
        st, species = fixture(T; D, ns)
        @test electron_ion_exchange_diagnostics(st, species) === nothing
        before = (species = deepcopy(species), pe = copy(st.inner.fields.pe))
        rng = bytes(st.rng)
        step!(st, species, st.dt)
        report = check_receipts(st, species, before)
        @test bytes(st.rng) == rng # Warm scaling has no thermal randomness.
        for cell in report.cells
            target, rates, Ce, Te = star_oracle(st, cell, [row.before for row in cell.ions])
            @test isapprox(BigFloat(cell.electron_target), target[1]; rtol = 32eps(T))
            if ns == 1
                setprecision(BigFloat, 256) do
                    C, Ti = BigFloat(cell.ions[1].before.capacity),
                    BigFloat(cell.ions[1].before.thermal / cell.ions[1].before.capacity)
                    equilibrium = (Ce * Te + C * Ti) / (Ce + C)
                    delta = (Ti - Te) / (1 + BigFloat(st.dt) * rates[1] * (1 + Ce / C))
                    @test isapprox(
                        target[1],
                        equilibrium - C * delta / (Ce + C);
                        rtol = BigFloat(2)^(-200),
                    )
                end
            end
            for (a, row) in enumerate(cell.ions)
                beta = BigFloat(st.dt) * Ce / BigFloat(row.before.capacity) * rates[a]
                beta /= 1 + beta
                @test BigFloat(row.beta[1]) <= beta <= BigFloat(row.beta[3])
                @test isapprox(
                    BigFloat(row.target),
                    BigFloat(row.before.capacity) * target[a+1];
                    rtol = 32eps(T),
                )
                @test row.after.thermal > 0
            end
        end
        # Four x nodes retain non-Nyquist pressure variation. With B=0 the
        # actual carried electric field must contain the new pressure force.
        gradient = similar(st.inner.fields.pe)
        deriv!(gradient, st.inner.fields.pe, st.grid, 1)
        @test maximum(abs, gradient) > 0
        expected_force = -gradient ./ st.inner.fields.n
        @test isapprox(
            st.inner.fields.E[1],
            expected_force;
            rtol = 128eps(T),
            atol = 128eps(T) * maximum(abs, expected_force),
        )
        @test any(cell -> cell.pressure_after != cell.pressure_before, report.cells)
        detached = electron_ion_exchange_diagnostics(st, species)
        empty!(detached.cells)
        @test length(electron_ion_exchange_diagnostics(st, species).cells) == 4
    end
end

@testset "Request binding and stored transport force refresh" begin
    for bad in (0, -1, Inf, NaN)
        @test_throws ArgumentError ElectronIonThermalExchange(;
            electron_mass = bad,
            coulomb_log = 1,
        )
    end
    @test_throws ArgumentError ElectronIonThermalExchange(; electron_mass = 1, coulomb_log = -1)
    @test_throws ArgumentError ElectronIonThermalExchange(;
        electron_mass = 1,
        coulomb_log = 1,
        log_bits = true,
    )
    tiny = Q(1, BigInt(2)^2000)
    @test_throws ArgumentError fixture(Float32; electron_mass = tiny)
    @test_throws ArgumentError fixture(Float32; logs = tiny)
    logs = [0.25]
    bound, bp = fixture(Float64; ns = 1, logs)
    logs[1] = 99
    @test H._checkpoint_state(bound, bp).electron_ion.logarithm == (0.25,)
    for closure in (
        TransportElectrons(5 / 3),
        ConductingElectrons(5 / 3; κ = 1 / 64),
        FieldAlignedConductingElectrons(5 / 3; κ_parallel = 1 / 64, κ_perpendicular = 1 / 64),
    )
        st, species = fixture(Float64; ns = 1, closure)
        before = (species = deepcopy(species), pe = copy(st.inner.fields.pe))
        step!(st, species, st.dt)
        check_receipts(st, species, before)
        gradient = similar(st.inner.fields.pe)
        deriv!(gradient, st.inner.fields.pe, st.grid, 1)
        @test maximum(abs, gradient) > 0
        expected_force = -gradient ./ st.inner.fields.n
        @test isapprox(
            st.inner.pressure_work.Eforce[1],
            expected_force;
            rtol = 128eps(),
            atol = 128eps() * maximum(abs, expected_force),
        )
    end
end

@testset "SI normalization carries one A and the SI permittivity" begin
    units = PlasmaUnits(8.0, 2.0, 3.0, 0.5, 0.25, 16.0)
    normalization = coulomb_normalization(units)
    st, species = fixture(
        Float64;
        ns = 1,
        epsilon = normalization.epsilon0,
        number_scale = normalization.weight_to_number,
        cell_width = 2.0,
    )
    before = (species = deepcopy(species), pe = copy(st.inner.fields.pe))
    step!(st, species, st.dt)
    cell = first(check_receipts(st, species, before).cells)
    target, rates, Ce, _ = star_oracle(st, cell, [cell.ions[1].before])
    lambda = BigFloat(st.dt) * Ce / BigFloat(cell.ions[1].before.capacity) * rates[1]
    @test BigFloat(cell.ions[1].beta[1]) <= lambda / (1 + lambda) <= BigFloat(cell.ions[1].beta[3])
    @test isapprox(BigFloat(cell.electron_target), target[1]; rtol = 32eps())
    setprecision(BigFloat, 256) do
        n0, B0, m0, e, mu, c =
            BigFloat.((units.n0, units.B0, units.mi, units.e, units.mu0, units.c))
        L, Omega, vA = sqrt(m0 / (mu * n0 * e^2)), e * B0 / m0, B0 / sqrt(mu * n0 * m0)
        epsilonSI = inv(mu * c^2)
        @test isapprox(BigFloat(normalization.weight_to_number), n0 * L^3; rtol = 16eps())
        @test isapprox(
            BigFloat(normalization.epsilon0),
            epsilonSI * m0 * L^3 * Omega^2 / e^2;
            rtol = 16eps(),
        )
        ion = cell.ions[1].before
        thetae = BigFloat(cell.pressure_before) / BigFloat(cell.density) * m0 * vA^2
        thetai = BigFloat(ion.thermal / ion.capacity) * m0 * vA^2
        me, mi = BigFloat(st.electron_ion.electron_mass) * m0, BigFloat(species[1].m) * m0
        physical_density = n0 * BigFloat(ion.weight) / prod(BigFloat, st.grid.dx)
        rateSI =
            8sqrt(2BigFloat(pi)) *
            physical_density *
            (e^2 * BigFloat(species[1].q))^2 *
            BigFloat(st.electron_ion.logarithm[1]) / (
                3 *
                (4BigFloat(pi) * epsilonSI)^2 *
                me *
                mi *
                (thetae / me + thetai / mi)^(BigFloat(3) / 2)
            )
        @test isapprox(rates[1], rateSI / Omega; rtol = 64eps())
    end
end

@testset "Cold unequal weights, sparse preflight, late rejection and RNG retry" begin
    for T in (Float32, Float64), R in (Xoshiro, MersenneTwister)
        st, species = fixture(T; ns = 1, cold = true, R)
        before = (species = deepcopy(species), pe = copy(st.inner.fields.pe))
        stream = deepcopy(st.rng)
        for _ = 1:(4*2*3)
            randn(stream, T)
        end
        step!(st, species, st.dt)
        report = check_receipts(st, species, before)
        @test bytes(st.rng) == bytes(stream)
        @test all(
            cell -> cell.entropy === :cold_positive && cell.ions[1].after.thermal > 0,
            report.cells,
        )
        @test all(c -> species[1].v[c][3:3:end] == before.species[1].v[c][3:3:end], 1:3)
        sparse, ps = fixture(T; ns = 1, cold = true, sparse = true, R)
        entry = snapshot(sparse, ps)
        @test_throws ArgumentError step!(sparse, ps, sparse.dt)
        @test isequal(snapshot(sparse, ps), entry)
        # h*bulk=2^-10 keeps the drift small; positive thermal increments are
        # below the bulk ULP. The specific post-draw rejection is checked below.
        late, lp =
            fixture(T; ns = 1, cold = true, R, bulk = T(2)^50, pressure = T(1e-6), dt = T(2)^(-60))
        entry = snapshot(late, lp)
        failure = try
            step!(late, lp, late.dt)
            nothing
        catch err
            err
        end
        @test failure isa ArgumentError
        @test occursin("positive ion thermal", sprint(showerror, failure))
        @test isequal(snapshot(late, lp), entry)
        for p in lp, v in p.v
            fill!(v, 0.0)
        end
        init!(late, lp)
        clean, cp = fixture(T; ns = 1, cold = true, R, pressure = T(1e-6), dt = T(2)^(-60))
        step!(late, lp, late.dt)
        step!(clean, cp, clean.dt)
        @test isequal(snapshot(late, lp), snapshot(clean, cp))
    end
end

@testset "Disabled and enabled-zero reduction, zero step and cold equilibrium" begin
    for T in (Float32, Float64)
        disabled, ps = fixture(T; enabled = false)
        inactive, zps = fixture(T; logs = zero(T))
        plain_particles = deepcopy(ps)
        plain = HybridStepper(
            FourierGrid(disabled.grid.n, disabled.grid.L),
            disabled.model,
            NGP(),
            sum(nparticles, plain_particles),
        )
        copyto!(plain.fields.pe, disabled.inner.fields.pe)
        init!(plain, plain_particles)
        @test H._checkpoint_state(disabled, ps).schema == 1
        @test H._checkpoint_state(inactive, zps).schema == 2
        @test_throws ArgumentError electron_ion_exchange_diagnostics(disabled, ps)
        entry = snapshot(inactive, zps)
        step!(inactive, zps, zero(T))
        @test isequal(snapshot(inactive, zps), entry)
        step!(disabled, ps, disabled.dt)
        step!(inactive, zps, inactive.dt)
        step!(plain, plain_particles, disabled.dt)
        @test H._hybrid_field_state(disabled.inner.fields) ==
              H._hybrid_field_state(inactive.inner.fields)
        @test H._hybrid_field_state(disabled.inner.fields) == H._hybrid_field_state(plain.fields)
        @test [H._checkpoint_particle_state(p) for p in ps] == [H._checkpoint_particle_state(p) for p in zps]
        @test [H._checkpoint_particle_state(p) for p in ps] == [H._checkpoint_particle_state(p) for p in plain_particles]
        @test bytes(disabled.rng) == bytes(inactive.rng)
        @test isempty(electron_ion_exchange_diagnostics(inactive, zps).cells)
        eq, ep = fixture(T; ns = 1, cold = true, pressure = zero(T))
        rng = bytes(eq.rng)
        step!(eq, ep, eq.dt)
        @test bytes(eq.rng) == rng
        @test all(
            cell -> cell.entropy === :unchanged,
            electron_ion_exchange_diagnostics(eq, ep).cells,
        )
        # This fixture has theta_i=21/512 and ne=4, hence pe=21/128.
        warm, wp = fixture(T; ns = 1, uniform = true, pressure = T(21 / 128))
        original = [map(copy, p.v) for p in wp]
        rng = bytes(warm.rng)
        step!(warm, wp, warm.dt)
        @test [p.v for p in wp] == original && bytes(warm.rng) == rng
        @test all(
            cell -> cell.entropy === :unchanged,
            electron_ion_exchange_diagnostics(warm, wp).cells,
        )
    end
end

@testset "Thermal raw/archive continuation and malformed receipts retain state" begin
    st, species = fixture(Float64; ns = 1, cold = true)
    step!(st, species, st.dt)
    saved = deepcopy(H._checkpoint_state(st, species))
    @test saved.schema == 2 && saved.electron_ion.policy == 1
    mktempdir() do dir
        path, archive = joinpath(dir, "thermal.chk"), joinpath(dir, "thermal.run")
        save_checkpoint(path, st, species)
        archive_run(archive, st, species; rng_seed = 913)
        target, tp = fixture(Float64; ns = 1, cold = true)
        identities = [map(objectid, (p.x..., p.v..., p.weight, p.id, p.tag)) for p in tp]
        load_checkpoint!(target, tp, path)
        @test isequal(snapshot(target, tp), snapshot(st, species))
        @test [map(objectid, (p.x..., p.v..., p.weight, p.id, p.tag)) for p in tp] == identities
        for _ = 1:2
            step!(st, species, st.dt)
            step!(target, tp, target.dt)
            @test isequal(snapshot(target, tp), snapshot(st, species))
        end
        load_archive!(target, tp, archive)
        baseline, baseps = fixture(Float64; ns = 1, cold = true)
        H._restore_checkpoint_state!(baseline, baseps, saved)
        @test isequal(snapshot(target, tp), snapshot(baseline, baseps))
        for damage in (:beta, :energy, :infinite, :role, :configuration, :schema)
            incoming = deepcopy(saved)
            if damage === :configuration
                incoming = merge(
                    incoming,
                    (electron_ion = merge(incoming.electron_ion, (electron_mass = 1 / 8,)),),
                )
            elseif damage === :schema
                incoming = merge(incoming, (schema = 1,))
            else
                cells = incoming.thermal_report.cells
                firstcell = cells[1]
                if damage === :beta
                    ion = firstcell.ions[1]
                    firstcell.ions[1] =
                        merge(ion, (beta = (ion.beta[1], ion.beta[2] + 1 // 8, ion.beta[3]),))
                elseif damage === :energy
                    cells[1] = merge(firstcell, (energy_defect = firstcell.energy_defect + 1,))
                elseif damage === :infinite
                    cells[1] = merge(firstcell, (energy_defect = Q(1, 0),))
                else
                    incoming = merge(
                        incoming,
                        (
                            thermal_report = merge(
                                incoming.thermal_report,
                                (cells = target.thermal_report.cells,),
                            ),
                        ),
                    )
                end
            end
            entry = snapshot(target, tp)
            @test_throws ArgumentError H._restore_checkpoint_state!(target, tp, incoming)
            @test isequal(snapshot(target, tp), entry)
        end
    end
end
end # module
