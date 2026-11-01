module EMPIC1DSpeciesTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC
Q(x) = Rational{BigInt}(x)

function fixture(::Type{T}; n = 8, shape = CIC(), n_sub = 1, relativistic = false) where {T}
    g = FourierGrid((n,), (one(T),))
    species = [ParticleSet{1,T}(k + 2; q, m = k + 1) for (k, q) in enumerate((-2, 1, 3, 0))]
    for (k, ps) in enumerate(species)
        ps.x[1] .= [T(mod(p / (nparticles(ps) + 1) + 0.07k, 1)) for p = 1:nparticles(ps)]
        fill!(ps.weight, zero(T))
        ps.weight[1] = k == 1 ? T(2) : one(T)
        if k != 3
            fill!(ps.v[1], T(0.04k))
            fill!(ps.v[2], T(0.02k))
        end
    end
    es = EMPIC1D(g, species; mobile = Bool[true, true, false, true], shape, n_sub, relativistic)
    return es, species
end

physical(es, species) = deepcopy((
    fields = H._empic_checkpoint_fields(es),
    species = [H._checkpoint_particle_state(ps) for ps in species],
))
function snapshot(es, species)
    io = IOBuffer()
    serialize(io, (fields = H._empic_transaction_fields(es),
        phase = (es.time[], es.step[], es.last_dt[], es.last_sub_dt[]),
        species = [H._checkpoint_particle_state(ps) for ps in species],
        initialized = es.collection.initialized[], owners = objectid.(es.collection.owners)))
    return take!(io)
end

@testset "Specialized planar whole-step reference" begin
    for T in (Float32, Float64)
        g = FourierGrid((1,), (one(T),))
        species = [ParticleSet{1,T}(1; q = -2, m = 3), ParticleSet{1,T}(1; q = 1, m = 5)]
        for (i, ps) in enumerate(species)
            ps.x[1][1], ps.v[1][1], ps.v[2][1] = T(i / 4), T(i / 16), T(-i / 32)
            ps.weight[1] = T(i)
        end
        es = EMPIC1D(g, species; shape = NGP())
        fill!(es.Ey, T(1 / 4))
        init_empic!(es, species)
        fill!(es.Ex, T(1 / 8))
        old = deepcopy(species)
        dt = T(1 / 32)
        step_empic!(es, species, dt)
        # Uniform fields/B=0: back-up then full Boris kick gives v0+qE*dt/(2m).
        # Expected values use original rational sources, no production mover/gather.
        for (ps, before) in zip(species, old), c = 1:2
            electric = c == 1 ? Q(T(1 / 8)) : Q(T(1 / 4))
            expected = T(Q(before.v[c][1]) + Q(before.q) * electric * Q(dt) / (2 * Q(before.m)))
            @test isapprox(ps.v[c][1], expected; rtol = 32eps(T))
        end
        for c = 1:2
            expected = T(sum(Q(ps.q) * Q(ps.weight[1]) * Q(ps.v[c][1]) for ps in species))
            current = c == 1 ? es.Jx[1] : es.Jy[1]
            field = c == 1 ? es.Ex[1] : es.Ey[1]
            initial = c == 1 ? T(1 / 8) : T(1 / 4)
            @test isapprox(current, expected; rtol = 32eps(T))
            @test isapprox(field, T(Q(initial) - Q(dt) * Q(expected)); rtol = 32eps(T))
        end
        @test all(iszero, es.Bz)
        @test all(ps -> all(iszero, ps.v[3]), species)
    end
end

# Independent four-node spectral derivative: only modes +/-k survive the
# first-derivative convention, including its even Nyquist null.
d4(a, k) = [k * (a[2] - a[4]) / 2, k * (a[3] - a[1]) / 2,
    k * (a[4] - a[2]) / 2, k * (a[1] - a[3]) / 2]
@testset "Specialized nonzero-B prime, rotation and Maxwell reference" begin
    for T in (Float32, Float64)
        ps = ParticleSet{1,T}(1; q = 1, m = 2)
        ps.x[1][1], ps.weight[1] = zero(T), zero(T)
        ps.v[1][1], ps.v[2][1] = T(1 / 8), T(1 / 16)
        es = EMPIC1D(FourierGrid((4,), (one(T),)), [ps]; shape = NGP(), c = 3)
        es.Ey .= T[0, 1 / 8, 0, -1 / 8]
        es.Bz .= T[5 / 16, 1 / 4, 3 / 16, 1 / 4]
        init_empic!(es, [ps])
        fill!(es.Ex, T(1 / 8))
        ex, ey, bz = Q(es.Ex[1]), Q(es.Ey[1]), Q(es.Bz[1])
        h, alpha = Q(T(1 / 64)), Q(ps.q) / Q(ps.m)
        vx, vy = Q(ps.v[1][1]), Q(ps.v[2][1])
        vx -= alpha * h / 2 * (ex + vy * bz)
        vy -= alpha * h / 2 * (ey - Q(ps.v[1][1]) * bz)
        ux, uy = vx + alpha * h / 2 * ex, vy + alpha * h / 2 * ey
        rotation = alpha * h / 2 * bz
        twice = 2 * rotation / (1 + rotation^2)
        vnew = (ux + (uy - ux * rotation) * twice + alpha * h / 2 * ex,
            uy - (ux + uy * rotation) * twice + alpha * h / 2 * ey)
        k = Q(imag(es.g.ik[1][2]))
        old_E, old_B = Q.(es.Ey), Q.(es.Bz)
        bnew = old_B .- h / 2 .* d4(old_E, k)
        enew = old_E .- h * Q(es.c)^2 .* d4(bnew, k)
        step_empic!(es, [ps], T(h))
        @test all(isapprox(ps.v[c][1], T(vnew[c]); rtol = 64eps(T)) for c = 1:2)
        @test isapprox(ps.x[1][1], T(h * vnew[1]); rtol = 64eps(T))
        @test isapprox(es.Bz, T.(bnew); rtol = 64eps(T), atol = 64eps(T))
        @test isapprox(es.Ey, T.(enew); rtol = 64eps(T), atol = 64eps(T))
        @test all(==(T(1 / 8)), es.Ex)
        @test all(iszero, es.Jx) && all(iszero, es.Jy)
    end
end

@testset "Specialized currents, represented continuity and synchronous phases" begin
    for T in (Float32, Float64), shape in (NGP(), CIC(), TSC())
        es, species = fixture(T; shape)
        @test es.collection.stage isa EMPIC1D
        @test_throws ArgumentError step_empic!(es, species, T(0.01))
        init_empic!(es, species)
        fixed = deepcopy(H._checkpoint_particle_state(species[3]))
        step_empic!(es, species, T(0.02))
        expected, scratch = zeros(T, es.g.n), zeros(T, es.g.n)
        for (i, ps) in enumerate(species)
            es.collection.mobile[i] || continue
            mid = deepcopy(ps)
            copyto!(mid.x[1], es.collection.midpoints[i][1])
            deposit_scalar!(scratch, mid, ps.q .* ps.weight .* ps.v[2], es.g, shape)
            expected .+= scratch ./ es.g.dx[1]
        end
        signs = T[isodd(i) ? 1 : -1 for i = 1:length(expected)]
        expected .-= signs .* (sum(signs .* expected) / length(expected))
        @test isapprox(es.Jy, expected; rtol = 256eps(T), atol = 256eps(T))
        @test H._checkpoint_particle_state(species[3]) == fixed
        d = electromagnetic_diagnostics(es, species)
        @test d.kinetic == kinetic_energy(species)
        @test d.momentum == total_momentum(species)
        @test d.field == em_field_energy(es)
        @test d.total == d.kinetic + d.field
        @test d.velocity_time_offsets == T[-es.last_sub_dt[] / 2, -es.last_sub_dt[] / 2, 0, -es.last_sub_dt[] / 2]
        if shape isa CIC
            delta = es.rho_np1 - es.rho_n
            delta .-= signs .* (sum(signs .* delta) / length(delta))
            source = sum(ps -> abs(ps.q) * sum(ps.weight), species) / es.g.L[1]
            @test maximum(abs, delta) > 64eps(T) * source
            @test charge_conservation_residual(es, es.last_sub_dt[]) < 1024eps(T)
        end
    end
    for T in (Float32, Float64)
        es, species = fixture(T; n = 7)
        init_empic!(es, species)
        fine_species = deepcopy(species)
        fine = EMPIC1D(FourierGrid(es.g.n, es.g.L), fine_species; mobile = copy(es.collection.mobile))
        init_empic!(fine, fine_species)
        for (dt, ns) in ((T(0.03), 3), (zero(T), 4), (T(0.02), 2))
            es.n_sub = ns
            before = snapshot(es, species)
            step_empic!(es, species, dt)
            if iszero(dt)
                @test snapshot(es, species) == before
            else
                for _ = 1:ns
                    step_empic!(fine, fine_species, dt / ns)
                end
                @test physical(es, species) == physical(fine, fine_species)
                @test es.time[] == fine.time[]
            end
        end
        # The count API is a separate coarse-ion model; compare n_sub=1 only.
        e, ions = ParticleSet{1,T}(3; q = -1), ParticleSet{1,T}(3; q = 1, m = 7)
        for ps in (e, ions)
            ps.x[1] .= T[0.1, 0.4, 0.8]
            ps.weight .= T[1, 2, 1]
            fill!(ps.v[1], T(0.03))
        end
        copied = deepcopy([e, ions])
        legacy = EMPIC1D(FourierGrid((8,), (one(T),)), 3; mobile = true, mi = 7)
        collection = EMPIC1D(FourierGrid((8,), (one(T),)), copied)
        init_empic!(legacy, e, ions)
        init_empic!(collection, copied)
        for dt in T[0.01, 0.015]
            step_empic!(legacy, e, ions, dt)
            step_empic!(collection, copied, dt)
            @test isapprox(legacy.Ex, collection.Ex; rtol = 128eps(T), atol = 128eps(T))
            @test isapprox(legacy.Ey, collection.Ey; rtol = 128eps(T), atol = 128eps(T))
            @test isapprox(e.x[1], copied[1].x[1]; rtol = 128eps(T), atol = 128eps(T))
        end
        empty = ParticleSet{1,T}[]
        vacuum = EMPIC1D(FourierGrid((4,), (one(T),)), empty; n_sub = 2)
        uninitialized = deepcopy(H._checkpoint_state(vacuum, empty))
        H._restore_checkpoint_state!(vacuum, empty, uninitialized)
        @test !vacuum.collection.initialized[]
        init_empic!(vacuum, empty)
        step_empic!(vacuum, empty, T(0.02))
        @test electromagnetic_diagnostics(vacuum, empty).total == zero(T)
    end
end

@testset "Specialized complete original-current rescue" begin
    for T in (Float32, Float64)
        species = [ParticleSet{1,T}(1; q) for q in (1, -1, 1)]
        for (i, ps) in enumerate(species)
            ps.x[1][1] = T(0.25)
            ps.weight[1] = i == 3 ? one(T) : floatmax(T)
            ps.v[1][1] = ps.v[2][1] = i == 3 ? T(0.25) : T(2)
        end
        es = EMPIC1D(FourierGrid((4,), (one(T),)), species; neutrality = :neutralize, shape = NGP())
        init_empic!(es, species)
        work = H._empic_collection_prepare!(es, species)
        for (ps, mid) in zip(work.staged, work.midpoints)
            copyto!(mid[1], ps.x[1])
        end
        @test H._empic1d_collection_Jx0(work.stage, work) == T(0.25)
        H._empic_collection_current!(work.stage, work)
        # Complete q*w*vy deposit is [0,1,0,0], then remove its Nyquist component.
        @test work.stage.Jy == T[1 / 4, 3 / 4, 1 / 4, -1 / 4]
    end
end

@testset "Specialized planar and late-output rejection retains public state" begin
    for T in (Float32, Float64)
        es, species = fixture(T)
        init_empic!(es, species)
        before = snapshot(es, species)
        @test_throws ArgumentError init_empic!(es, species[1])
        @test_throws ArgumentError step_empic!(es, species[1], T(0.01))
        @test snapshot(es, species) == before
        U = T === Float32 ? UInt32 : UInt64
        stored = reinterpret(T, one(U))
        species[end].v[3][1] = stored
        @test reinterpret(U, species[end].v[3][1]) == one(U)
        bad = snapshot(es, species)
        oldmode = get_zero_subnormals()
        try
            set_zero_subnormals(true)
            @test_throws ArgumentError step_empic!(es, species, T(0.01))
        finally
            set_zero_subnormals(oldmode)
        end
        @test snapshot(es, species) == bad
        species[end].v[3][1] = zero(T)
        first = species[1]
        readonly = ParticleSet{1,T}(first.x, first.v, Base.Experimental.Const(first.weight),
            first.id, first.tag, first.q, first.m)
        @test_throws ArgumentError EMPIC1D(es.g, ParticleSet{1,T}[readonly, species[2:end]...]; mobile = copy(es.collection.mobile))
        # Reuse the existing C1 boundary: first prime succeeds privately, a later
        # tiny-mass population rejects. Its zero weight keeps source charge finite.
        failed = [ParticleSet{1,T}(1; m) for m in (one(T), stored)]
        for ps in failed
            ps.x[1][1], ps.weight[1] = zero(T), zero(T)
        end
        rel = EMPIC1D(FourierGrid((4,), (T(4),)), failed; relativistic = true, c = one(T))
        init_empic!(rel, failed)
        fill!(rel.Ex, one(T))
        entry = snapshot(rel, failed)
        @test_throws ArgumentError step_empic!(rel, failed, one(T))
        @test snapshot(rel, failed) == entry
        @test rel.collection.staged[1].v != failed[1].v
        failed[end].m = one(T)
        step_empic!(rel, failed, T(0.125))
        @test rel.step[] == 1
        good, real_species = fixture(T; relativistic = true)
        fill!(good.Ey, T(0.0625))
        fill!(good.Bz, T(0.125))
        init_empic!(good, real_species)
        step_empic!(good, real_species, T(0.01))
        # Independent high-precision physical energy on stored half-time v.
        expected = setprecision(BigFloat, 256) do
            sum(BigFloat(ps.weight[p]) * BigFloat(ps.m) *
                (inv(sqrt(1 - sum(BigFloat(ps.v[c][p])^2 for c = 1:3) / BigFloat(good.c)^2)) - 1) * BigFloat(good.c)^2
                for ps in real_species for p in eachindex(ps.weight))
        end
        @test expected > 0
        @test isapprox(electromagnetic_diagnostics(good, real_species).kinetic, T(expected); rtol = 32eps(T))
    end
end

@testset "Specialized distinct checkpoint/archive and fixed strided restore" begin
    for T in (Float32, Float64)
        es, species = fixture(T; n_sub = 2)
        init_empic!(es, species)
        step_empic!(es, species, T(0.02))
        destination = [ParticleSet{1,T}(0) for _ = 1:4]
        target = EMPIC1D(FourierGrid(es.g.n, es.g.L), destination;
            mobile = copy(es.collection.mobile), c = T(0.1), relativistic = true)
        mktempdir() do dir
            path = joinpath(dir, "specialized.chk")
            save_checkpoint(path, es, species)
            load_checkpoint!(target, destination, path)
            @test physical(target, destination) == physical(es, species)
            @test !target.relativistic && target.c == es.c
            for (dt, ns) in ((T(0.03), 3), (zero(T), 1), (T(0.01), 2))
                target.n_sub = es.n_sub = ns
                step_empic!(es, species, dt)
                step_empic!(target, destination, dt)
                @test physical(target, destination) == physical(es, species)
                @test (target.time[], target.step[], target.last_dt[], target.last_sub_dt[]) ==
                    (es.time[], es.step[], es.last_dt[], es.last_sub_dt[])
            end
            archive = joinpath(dir, "specialized.run")
            archive_run(archive, es, species; rng_seed = 17)
            saved = load_archive(archive).state
            @test saved.integrator === :empic1d_collection && saved.polarization === :ex_ey_bz
            before = snapshot(target, destination)
            for changes in ((schema_version = 2,), (schema_version = true,), (polarization = :full,), (last_sub_dt = zero(T),), (relativistic = true, c = eps(T)))
                @test_throws ArgumentError H._restore_checkpoint_state!(target, destination, merge(saved, changes))
                @test snapshot(target, destination) == before
            end
            for (population, component) in ((1, 3), (3, 1))
                malformed = deepcopy(saved)
                malformed.species[population].v[component][1] = T(0.125)
                @test_throws ArgumentError H._restore_checkpoint_state!(target, destination, malformed)
                @test snapshot(target, destination) == before
            end
            generic = EMPIC(FourierGrid(es.g.n, es.g.L), destination; mobile = copy(es.collection.mobile))
            @test_throws ArgumentError load_checkpoint!(generic, destination, path)
            scalar = EMPIC1D(FourierGrid(es.g.n, es.g.L), 0)
            @test_throws ArgumentError load_checkpoint!(scalar, ParticleSet{1,T}(0; q = -1), path)
            @test_throws ArgumentError save_checkpoint(path, target, destination[1])
            swapped = merge(deepcopy(saved), (fields = merge(saved.fields, (Ex = target.Ey, Ey = target.Ex)),))
            expected = (copy(target.Ey), copy(target.Ex))
            H._restore_checkpoint_state!(target, destination, swapped)
            @test target.Ex == expected[1] && target.Ey == expected[2]
            function strided(a)
                backing = zeros(eltype(a), 2length(a) + 1)
                view = @view backing[2:2:2length(a)]
                copyto!(view, a)
                return view
            end
            views = ParticleSet{1,T}[ParticleSet{1,T}(map(strided, ps.x), map(strided, ps.v),
                strided(ps.weight), strided(ps.id), strided(ps.tag), ps.q, ps.m) for ps in species]
            viewed = EMPIC1D(FourierGrid(es.g.n, es.g.L), views; mobile = copy(es.collection.mobile))
            identities = [ps.x[1] for ps in views]
            H._restore_checkpoint_state!(viewed, views, saved)
            @test all(ps.x[1] === a for (ps, a) in zip(views, identities))
            @test physical(viewed, views) == physical(es, species)
            records = Any[saved.species...]
            last = records[end]
            records[end] = merge(last, (x = map(a -> a[1:(end-1)], last.x),
                v = map(a -> a[1:(end-1)], last.v), weight = last.weight[1:(end-1)],
                id = last.id[1:(end-1)], tag = last.tag[1:(end-1)]))
            entry = snapshot(viewed, views)
            @test_throws ArgumentError H._restore_checkpoint_state!(viewed, views, merge(saved, (species = records,)))
            @test snapshot(viewed, views) == entry
        end
    end
end
end # module
