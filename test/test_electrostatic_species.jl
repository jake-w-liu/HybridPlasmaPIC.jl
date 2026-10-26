module ElectrostaticSpeciesTests
using HybridPlasmaPIC, Test, Serialization, LinearAlgebra, Random
const H = HybridPlasmaPIC

function populations(::Type{T}, D; pair = false) where {T}
    g = FourierGrid(ntuple(_ -> 8, D), ntuple(_ -> T(2π), D))
    counts = ntuple(_ -> 4, D)
    qs, ms, ns = pair ? ((-1, 1), (1, 1), (1, 1)) : ((-2, 1), (3, 5), (1, 2))
    species = ParticleSet{D,T}[]
    for k = 1:2
        ps = ParticleSet{D,T}(prod(counts); q = qs[k], m = ms[k])
        load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
        set_density_weight!(ps, T(ns[k]), g)
        ps.x[1] .+= T(0.025k) .* sin.(ps.x[1])
        for c = 1:3
            ps.v[c] .= T(0.01c) .* cos.(ps.x[1])
        end
        push!(species, ps)
    end
    return g, species
end

function snapshot(es, species)
    io = IOBuffer()
    serialize(
        io,
        (
            E = es.E,
            ne = es.ne,
            primed = es.primed[],
            last_dt = es.last_dt[],
            species = [H._checkpoint_particle_state(ps) for ps in species],
        ),
    )
    return take!(io)
end

@testset "ES collection shared charge, field and diagnostics" begin
    for T in (Float32, Float64), D = 1:3, pair in (false, true)
        g, species = populations(T, D; pair)
        es = ElectrostaticPIC(g, species)
        @test_throws ArgumentError step_espic!(es, species, T(0.01))
        init_espic!(es, species)
        expected, tmp = zeros(T, g.n), zeros(T, g.n)
        for ps in species
            density!(tmp, ps, g, CIC())
            expected .-= ps.q .* tmp
        end
        @test es.ne == expected
        diagnostics = electrostatic_diagnostics(es, species)
        @test abs(diagnostics.net_charge) <= 64eps(T) * prod(g.L)
        @test diagnostics.kinetic == kinetic_energy(species)
        @test diagnostics.momentum == total_momentum(species)
        @test diagnostics.field == field_energy(es)
        @test diagnostics.total == diagnostics.kinetic + diagnostics.field
        @test diagnostics.velocity_time_offset == 0
        @test isapprox(
            diagnostics.deposited_charge,
            diagnostics.net_charge;
            atol = 128eps(T) * prod(g.L),
        )
        p0 = total_momentum(species)
        for dt in T[0.02, 0.01, 0, 0.03]
            before = snapshot(es, species)
            step_espic!(es, species, dt)
            iszero(dt) && @test snapshot(es, species) == before
            @test all(ps -> all(a -> all(isfinite, a), (ps.x..., ps.v...)), species)
            @test isapprox(
                collect(total_momentum(species)),
                collect(p0);
                atol = 512eps(T) * prod(g.L),
            )
        end
        @test electrostatic_diagnostics(es, species).velocity_time_offset == -T(0.03) / 2
    end
end

@testset "ES collection scalar path and population permutations" begin
    for T in (Float32, Float64), D = 1:3
        g, species = populations(T, D)
        e = species[1]
        e.q, e.m = -one(T), one(T)
        scalar = ElectrostaticPIC(g, nparticles(e); n0 = one(T))
        other = deepcopy(e)
        collection = ElectrostaticPIC(
            FourierGrid(g.n, g.L),
            [other];
            background_charge = one(T),
            neutrality = :neutralize,
        )
        init_espic!(scalar, e)
        init_espic!(collection, [other])
        @test scalar.ne == collection.ne
        @test scalar.E == collection.E
        for dt in T[0.01, 0.02, 0, 0.015]
            step_espic!(scalar, e, dt)
            step_espic!(collection, [other], dt)
            @test e.x == other.x
            @test e.v == other.v
            @test scalar.ne == collection.ne
            @test scalar.E == collection.E
            @test scalar.primed[] == collection.primed[]
            @test scalar.last_dt[] == collection.last_dt[]
        end
        g, original = populations(T, D)
        neutral = ParticleSet{D,T}(3; q = 0, m = 7)
        for x in neutral.x
            x .= T[0.1, 0.3, 0.7]
        end
        neutral.v[1] .= T(0.1)
        empty = ParticleSet{D,T}(0; q = 4, m = 0.3)
        push!(original, neutral, empty)
        reversed = reverse(deepcopy(original))
        a = ElectrostaticPIC(g, original)
        b = ElectrostaticPIC(FourierGrid(g.n, g.L), reversed)
        init_espic!(a, original)
        init_espic!(b, reversed)
        x0 = copy(neutral.x[1])
        for dt in T[0.01, 0.02]
            step_espic!(a, original, dt)
            step_espic!(b, reversed, dt)
        end
        @test a.E == b.E
        @test all(p.x == q.x && p.v == q.v for (p, q) in zip(original, reverse(reversed)))
        @test neutral.x[1] ≈ x0 .+ T(0.003)
        before = snapshot(a, original)
        @test_throws ArgumentError step_espic!(a, reverse(original), T(0.01))
        @test snapshot(a, original) == before
    end
end

@testset "ES collection neutrality and empty populations" begin
    for T in (Float32, Float64), D = 1:3
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> one(T), D))
        empty = ParticleSet{D,T}[]
        es = ElectrostaticPIC(g, empty)
        init_espic!(es, empty)
        step_espic!(es, empty, T(0.1))
        @test iszero(field_energy(es))
        @test electrostatic_diagnostics(es, empty).momentum == (0, 0, 0)
        @test electrostatic_diagnostics(es, empty).kinetic == 0
        @test_throws ArgumentError ElectrostaticPIC(g, empty; background_charge = 1)
        ps = ParticleSet{D,T}(2; q = 2, m = 3)
        foreach(x -> fill!(x, T(0.3)), ps.x)
        @test_throws ArgumentError ElectrostaticPIC(g, [ps])
        balanced = ElectrostaticPIC(g, [ps]; background_charge = -4)
        init_espic!(balanced, [ps])
        @test electrostatic_diagnostics(balanced, [ps]).net_charge == 0
        neutralized = ElectrostaticPIC(g, [ps]; neutrality = :neutralize)
        init_espic!(neutralized, [ps])
        d = electrostatic_diagnostics(neutralized, [ps])
        @test d.particle_charge == 4
        @test d.neutralizing_charge == -4
        @test d.net_charge == 0
        # This is the exact integral of the represented CIC grid, whose
        # ordinary floating deposition can round relative to physical charge.
        @test d.deposited_charge ==
              -sum(Rational{BigInt}, neutralized.ne) * prod(Rational{BigInt}, g.dx)
        @test abs(d.deposited_charge - 4) <= 64eps(T) * 4
        @test_throws ArgumentError ElectrostaticPIC(g, [ps]; neutrality = :discard)
        @test_throws ArgumentError ElectrostaticPIC(g, [ps]; charge_rtol = NaN)
        @test_throws ArgumentError ElectrostaticPIC(g, [ps]; background_charge = Inf)
        @test_throws ArgumentError init_espic!(neutralized, ps)
    end
end

@testset "ES collection validation and arithmetic failures do not publish" begin
    for T in (Float32, Float64), D = 1:3
        g, species = populations(T, D)
        es = ElectrostaticPIC(g, species)
        init_espic!(es, species)
        step_espic!(es, species, T(0.01))
        last = species[end]
        for (get, put, bad) in (
            (() -> last.q, v -> (last.q = v), T(NaN)),
            (() -> last.m, v -> (last.m = v), zero(T)),
            (() -> last.weight[1], v -> (last.weight[1] = v), -one(T)),
            (() -> last.x[1][1], v -> (last.x[1][1] = v), T(Inf)),
            (() -> last.v[3][1], v -> (last.v[3][1] = v), T(NaN)),
        )
            old = get()
            put(bad)
            before = snapshot(es, species)
            @test_throws ArgumentError step_espic!(es, species, T(0.01))
            @test snapshot(es, species) == before
            @test_throws ArgumentError init_espic!(es, species)
            @test snapshot(es, species) == before
            put(old)
        end
        before = snapshot(es, species)
        for dt in (-one(T), T(Inf), T(NaN))
            @test_throws ArgumentError step_espic!(es, species, dt)
            @test snapshot(es, species) == before
        end
        @test_throws ArgumentError ElectrostaticPIC(
            g,
            [species[1], species[1]];
            neutrality = :neutralize,
        )
        aliased = deepcopy(last)
        aliased.v = (aliased.x[1], aliased.v[2], aliased.v[3])
        @test_throws ArgumentError ElectrostaticPIC(g, [aliased]; neutrality = :neutralize)
        if D == 1
            aliased = ParticleSet{1,T}(length(es.ne); q = 0)
            aliased.x = (es.ne,)
            @test_throws ArgumentError init_espic!(es, [aliased])
        end
        # A finite proposed drift overflows only in the last population. The
        # earlier population has already completed its staged kick and drift.
        last.v[1] .= floatmax(T)
        before = snapshot(es, species)
        @test_throws ArgumentError step_espic!(es, species, T(4))
        @test snapshot(es, species) == before
    end
end

@testset "ES collection checkpoint and archive preserve variable-step continuation" begin
    for T in (Float32, Float64), D = 1:3, live in (false, true)
        g, species = populations(T, D)
        es = ElectrostaticPIC(g, species)
        init_espic!(es, species)
        live && step_espic!(es, species, T(0.01))
        mktempdir() do dir
            path = joinpath(dir, "checkpoint")
            save_checkpoint(path, es, species)
            target_species = [ParticleSet{D,T}(1; q = ps.q, m = ps.m) for ps in species]
            foreach(ps -> foreach(x -> fill!(x, T(0.1)), ps.x), target_species)
            # Destination old charge need not equal the saved physical charge.
            target = ElectrostaticPIC(FourierGrid(g.n, g.L), ParticleSet{D,T}[])
            load_checkpoint!(target, target_species, path)
            @test snapshot(target, target_species) == snapshot(es, species)
            for dt in T[0.02, 0, 0.007, 0.03]
                step_espic!(es, species, dt)
                step_espic!(target, target_species, dt)
                @test snapshot(target, target_species) == snapshot(es, species)
            end
            archive = joinpath(dir, "archive")
            archive_run(archive, es, species; rng_seed = 17, normalization = "epsilon0=1")
            loaded = load_archive(archive).state
            H._restore_checkpoint_state!(target, target_species, loaded)
            @test snapshot(target, target_species) == snapshot(es, species)
            before = snapshot(target, target_species)
            bad = merge(loaded, (last_dt = -one(T),))
            @test_throws ArgumentError H._restore_checkpoint_state!(target, target_species, bad)
            @test snapshot(target, target_species) == before
            records = deepcopy(loaded.species)
            records[end] = merge(records[end], (m = zero(T),))
            @test_throws ArgumentError H._restore_checkpoint_state!(
                target,
                target_species,
                merge(loaded, (species = records,)),
            )
            @test snapshot(target, target_species) == before
            @test_throws ArgumentError H._restore_checkpoint_state!(
                target,
                [target_species[1], target_species[1]],
                loaded,
            )
            @test snapshot(target, target_species) == before
        end
    end
end

function cold_mode(dt, pair)
    T = Float64
    g = FourierGrid((64,), (2π,))
    x0 = collect(range(0.0, 2π; length = 513))[1:512] .+ 0.003
    qs, ms, ns =
        pair ? ((-1.0, 1.0), (1.0, 1.0), (1.0, 1.0)) : ((-2.0, 1.0), (3.0, 5.0), (1.0, 2.0))
    species = [ParticleSet{1,T}(length(x0); q = qs[k], m = ms[k]) for k = 1:2]
    amplitude = 1e-5
    for (k, ps) in enumerate(species)
        ps.x[1] .= x0 .+ (ps.q / ps.m) .* amplitude .* sin.(x0)
        set_density_weight!(ps, ns[k], g)
    end
    es = ElectrostaticPIC(g, species)
    init_espic!(es, species)
    initial_energy = field_energy(es)
    nsteps = round(Int, 1 / dt)
    for _ = 1:nsteps
        step_espic!(es, species, dt)
    end
    displacement = [
        2sum((ps.x[1] .- x0) .* sin.(x0)) / length(x0) / (ps.q / ps.m) / amplitude for ps in species
    ]
    omega = sqrt(sum(ns[k] * qs[k]^2 / ms[k] for k = 1:2))
    return displacement, cos(omega), species, es, initial_energy
end

function centered_energy(es, species)
    energy = field_energy(es)
    for ps in species
        Ep = ntuple(_ -> zeros(eltype(ps), nparticles(ps)), 3)
        gather_vector!(Ep, es.E, ps, es.g, CIC())
        for p in eachindex(ps.weight)
            centered = ntuple(c -> ps.v[c][p] + es.last_dt[] / 2 * ps.q / ps.m * Ep[c][p], 3)
            energy += ps.m * ps.weight[p] * sum(abs2, centered) / 2
        end
    end
    return energy
end

@testset "ES collection pair and unequal q/m cold-mode convergence" begin
    for pair in (false, true)
        coarse, exact, coarse_species, coarse_es, initial = cold_mode(1 / 20, pair)
        fine, _, fine_species, fine_es, _ = cold_mode(1 / 40, pair)
        reference, _, species, es, _ = cold_mode(1 / 320, pair)
        @test maximum(abs, reference .- exact) < 0.004
        @test norm(coarse - reference) > 3.5norm(fine - reference)
        @test maximum(abs, collect(total_momentum(species))) < 1e-13
        @test all(isfinite, electrostatic_diagnostics(es, species).momentum)
        reference_energy = centered_energy(es, species)
        coarse_error = abs(centered_energy(coarse_es, coarse_species) - reference_energy)
        fine_error = abs(centered_energy(fine_es, fine_species) - reference_energy)
        @test coarse_error > 3.5fine_error
        @test abs(reference_energy - initial) / initial < 0.005
    end
end

@testset "ES charge sum retains complete range across species" begin
    for T in (Float32, Float64)
        g = FourierGrid((2,), (one(T),))
        species = [ParticleSet{1,T}(1; q = q) for q in (one(T), -one(T), one(T))]
        for ps in species
            ps.x[1][1] = zero(T)
        end
        species[1].weight[1] = floatmax(T)
        species[2].weight[1] = floatmax(T)
        rho, scratch = zeros(T, 2), zeros(T, 2)
        state = H._charge_density_sum!(rho, scratch, species, g, CIC())
        @test state.eligible && state.lost
        @test rho == T[2, 0]
        H._charge_density_sum!(rho, scratch, reverse(species), g, CIC())
        @test rho == T[2, 0]
        # Number density underflows before the nonunit charge restores it.
        large = T === Float32 ? T(2)^100 : T(2)^800
        g = FourierGrid((2,), (large,))
        ps = ParticleSet{1,T}(1; q = large)
        ps.x[1][1] = zero(T)
        ps.weight[1] = T === Float32 ? T(2)^(-100) : T(2)^(-800)
        H._charge_density_sum!(rho, scratch, [ps], g, CIC())
        @test rho[1] == 2ps.weight[1]
        @test iszero(rho[2])
    end
end

@testset "ES field energy retains complete-expression range" begin
    for T in (Float32, Float64), State in (Electrostatic1D, ElectrostaticPIC)
        k = T === Float32 ? 100 : 800
        for sign in (-1, 1)
            g = FourierGrid((2,), (T(2)^(-sign * k),))
            es = State(g, 0; n0 = 0)
            field = es isa Electrostatic1D ? es.E : es.E[1]
            field .= T(2)^(sign * k)
            @test field_energy(es) == T(2)^(sign * k - 1)
        end
    end
end

@testset "ES collection phase, checkpoint schema and strided storage boundaries" begin
    for T in (Float32, Float64)
        g, species = populations(T, 1)
        views = [
            ParticleSet{1,T}(
                map(a -> view(a, :), ps.x),
                map(a -> view(a, :), ps.v),
                view(ps.weight, :),
                view(ps.id, :),
                view(ps.tag, :),
                ps.q,
                ps.m,
            ) for ps in species
        ]
        es = ElectrostaticPIC(g, views)
        @test_throws ArgumentError electrostatic_diagnostics(es, views)
        init_espic!(es, views)
        step_espic!(es, views, T(0.01))
        phase = es.last_dt[]
        es.last_dt[] = -one(T)
        before = snapshot(es, views)
        @test_throws ArgumentError step_espic!(es, views, T(0.01))
        @test snapshot(es, views) == before
        es.last_dt[] = phase
        saved = deepcopy(H._checkpoint_state(es, views))
        before = snapshot(es, views)
        for bad in (
            merge(saved, (schema_version = 2,)),
            merge(saved, (initialized = false,)),
            merge(saved, (neutrality = :neutralize,)),
            merge(saved, (charge_rtol = zero(T),)),
            merge(saved, (n0 = one(T),)),
            merge(saved, (species = saved.species[1:1],)),
            merge(saved, (E = (zeros(T, 3), saved.E[2], saved.E[3]),)),
        )
            @test_throws ArgumentError H._restore_checkpoint_state!(es, views, bad)
            @test snapshot(es, views) == before
        end
        H._restore_checkpoint_state!(es, views, saved)
        @test snapshot(es, views) == before
        # A heterogeneous NamedTuple vector can carry malformed/changed storage
        # records without failing inside the fixture's own setindex! conversion.
        too_small = NamedTuple[deepcopy(record) for record in saved.species]
        for i in eachindex(too_small)
            ps = too_small[i]
            too_small[i] = merge(
                ps,
                (
                    x = map(a -> a[1:1], ps.x),
                    v = map(a -> a[1:1], ps.v),
                    weight = ps.weight[1:1],
                    id = ps.id[1:1],
                    tag = ps.tag[1:1],
                ),
            )
        end
        @test_throws ArgumentError H._restore_checkpoint_state!(
            es,
            views,
            merge(saved, (species = too_small,)),
        )
        @test snapshot(es, views) == before
        underflow = ParticleSet{1,T}(0; q = nextfloat(zero(T)), m = floatmax(T))
        @test_throws ArgumentError ElectrostaticPIC(g, [underflow])
    end
end

@testset "ES strict neutrality distinguishes rounded cancellation" begin
    g = FourierGrid((4,), (1.0,))
    q = nextfloat(1.0)
    a = ParticleSet{1,Float64}(1; q = q)
    b = ParticleSet{1,Float64}(1; q = -(q * q))
    a.x[1][1] = b.x[1][1] = 0.25
    a.weight[1] = q
    @test a.q * a.weight[1] + b.q * b.weight[1] == 0
    exact = Rational{BigInt}(q)^2 - Rational{BigInt}(q * q)
    @test exact > 0
    @test_throws ArgumentError ElectrostaticPIC(g, [a, b]; charge_rtol = 0)
    es = ElectrostaticPIC(g, [a, b]; neutrality = :neutralize)
    init_espic!(es, [a, b])
    @test electrostatic_diagnostics(es, [a, b]).particle_charge == exact
    @test electrostatic_diagnostics(es, [a, b]).neutralizing_charge == -exact
end

@testset "ES collection public Poisson refresh and exact deposited charge" begin
    for T in (Float32, Float64), tolerance in (zero(T), 64eps(T))
        g = FourierGrid((8,), (one(T),))
        a, b = [ParticleSet{1,T}(17; q = q) for q in (-1, 1)]
        a.x[1] .= rand(MersenneTwister(1729), T, 17)
        b.x[1] .= mod.(a.x[1] .+ (T === Float32 ? T(1e-4) : T(1e-9)), one(T))
        es = ElectrostaticPIC(g, [a, b]; charge_rtol = tolerance)
        @test_throws ArgumentError poisson_E!(es)
        init_espic!(es, [a, b])
        step_espic!(es, [a, b], T(0.01))
        before = snapshot(es, [a, b])
        @test poisson_E!(es) === es.E
        @test snapshot(es, [a, b]) == before
        # Collection ne is derived. Refresh uses the bound populations while
        # preserving phase, rather than treating an edited buffer as physical charge.
        es.ne .+= one(T)
        @test poisson_E!(es) === es.E
        @test snapshot(es, [a, b]) == before
        b.weight[1] *= 2
        before = snapshot(es, [a, b])
        @test_throws ArgumentError poisson_E!(es)
        @test snapshot(es, [a, b]) == before
    end
    for T in (Float32, Float64)
        empty = ParticleSet{1,T}[]
        es = ElectrostaticPIC(FourierGrid((4,), (one(T),)), empty)
        init_espic!(es, empty)
        es.ne .= T[floatmax(T), floatmax(T), -floatmax(T), -floatmax(T)]
        @test electrostatic_diagnostics(es, empty).deposited_charge == 0
    end
end

@testset "ES collection replacement storage is validated before publication" begin
    g = FourierGrid((8,), (1.0,))
    species = [ParticleSet{1,Float64}(2; q = q) for q in (-1, 1)]
    for ps in species
        ps.x[1] .= [0.2, 0.6]
    end
    es = ElectrostaticPIC(g, species)
    init_espic!(es, species)
    step_espic!(es, species, 0.01)
    before = snapshot(es, species)
    views = [
        ParticleSet{1,Float64}(
            map(a -> view(a, :), ps.x),
            map(a -> view(a, :), ps.v),
            view(ps.weight, :),
            view(ps.id, :),
            view(ps.tag, :),
            ps.q,
            ps.m,
        ) for ps in species
    ]
    @test_throws ArgumentError init_espic!(es, views)
    @test snapshot(es, species) == before
    saved = deepcopy(H._checkpoint_state(es, species))
    @test_throws ArgumentError H._restore_checkpoint_state!(es, views, saved)
    @test snapshot(es, species) == before
end

# Serialization records storage types and parent links. Cross-representation
# continuation compares all physical values, while failure tests above retain
# byte snapshots to detect any publication into the unchanged representation.
physical_values(es, species) =
    (es.E, es.ne, es.primed[], es.last_dt[], [H._checkpoint_particle_state(ps) for ps in species])

@testset "ES checkpoint copies across supported source and destination storage" begin
    for T in (Float32, Float64), D = 1:3
        g, arrays = populations(T, D)
        views = [
            ParticleSet{D,T}(
                map(a -> view(a, :), ps.x),
                map(a -> view(a, :), ps.v),
                view(ps.weight, :),
                view(ps.id, :),
                view(ps.tag, :),
                ps.q,
                ps.m,
            ) for ps in deepcopy(arrays)
        ]
        a = ElectrostaticPIC(g, arrays)
        v = ElectrostaticPIC(FourierGrid(g.n, g.L), views)
        init_espic!(a, arrays)
        step_espic!(a, arrays, T(0.01))
        mktempdir() do dir
            path = joinpath(dir, "checkpoint")
            save_checkpoint(path, a, arrays)
            load_checkpoint!(v, views, path)
            @test physical_values(v, views) == physical_values(a, arrays)
            step_espic!(v, views, T(0.02))
            save_checkpoint(path, v, views)
            load_checkpoint!(a, arrays, path)
            @test physical_values(a, arrays) == physical_values(v, views)
            step_espic!(a, arrays, T(0.007))
            step_espic!(v, views, T(0.007))
            @test physical_values(a, arrays) == physical_values(v, views)
        end
    end
end

@testset "ES explicitly heterogeneous particle storage" begin
    for T in (Float32, Float64), D = 1:3
        g, species = populations(T, D)
        ps = species[2]
        species[2] = ParticleSet{D,T}(
            map(a -> view(a, :), ps.x),
            map(a -> view(a, :), ps.v),
            view(ps.weight, :),
            view(ps.id, :),
            view(ps.tag, :),
            ps.q,
            ps.m,
        )
        es = ElectrostaticPIC(g, species)
        init_espic!(es, species)
        step_espic!(es, species, T(0.01))
        before = snapshot(es, species)
        poisson_E!(es)
        @test snapshot(es, species) == before
        @test electrostatic_diagnostics(es, species).field == field_energy(es)
        saved = deepcopy(H._checkpoint_state(es, species))
        step_espic!(es, species, T(0.02))
        H._restore_checkpoint_state!(es, species, saved)
        @test snapshot(es, species) == before
        # A vector with a deliberately abstract element type retains its ability
        # to bind other supported storage after the caller resets velocities.
        arrays = [
            ParticleSet{D,T}(
                map(collect, ps.x),
                map(collect, ps.v),
                collect(ps.weight),
                collect(ps.id),
                collect(ps.tag),
                ps.q,
                ps.m,
            ) for ps in species
        ]
        foreach(ps -> foreach(a -> fill!(a, zero(T)), ps.v), arrays)
        init_espic!(es, arrays)
        step_espic!(es, arrays, T(0.01))
        @test es.primed[]
    end
end

end
