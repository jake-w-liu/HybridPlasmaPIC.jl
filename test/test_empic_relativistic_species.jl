module EMPICRelativisticSpeciesTests
using HybridPlasmaPIC, Test, Serialization
const H = HybridPlasmaPIC

function setup(::Type{T}, D; relativistic = true, c = T(5), n_sub = 1) where {T}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
    species = [ParticleSet{D,T}(prod(g.n); q, m) for (q, m) in ((-2, 3), (1, 5), (0, 7))]
    for (i, ps) in enumerate(species)
        load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, g.n)
        set_density_weight!(ps, i == 2 ? T(2) : one(T), g)
        i == 3 && continue
        fill!(ps.v[1], T(0.03i) * c)
        fill!(ps.v[2], T(0.01i) * c)
    end
    state = EMPIC(g, species; mobile = Bool[true, true, false], relativistic, c, n_sub)
    fill!(state.E[2], T(0.02))
    fill!(state.B[3], T(0.03))
    init_empic!(state, species)
    return state, species
end

function physical(state, species)
    return deepcopy((
        fields = H._empic_checkpoint_fields(state),
        particles = [H._checkpoint_particle_state(ps) for ps in species],
    ))
end
function snapshot(state, species)
    io = IOBuffer()
    serialize(
        io,
        (
            H._empic_transaction_fields(state),
            (
                state.time[],
                state.step[],
                state.last_dt[],
                state.last_sub_dt[],
                state.c,
                state.relativistic,
            ),
            [H._checkpoint_particle_state(ps) for ps in species],
            state.collection.initialized[],
        ),
    )
    return take!(io)
end

function continuity_observation(state, sub_dt)
    g = state.g
    density = complex.(state.rho_np1 .- state.rho_n)
    g.plan * density
    dc = first(CartesianIndices(density))
    for I in CartesianIndices(density)
        I != dc && H._spectral_mode_is_null(H._spectral_wavevector(g, I)) &&
            (density[I] = zero(eltype(density)))
    end
    g.iplan * density
    divergence = similar(state.rho_n)
    divergence!(divergence, state.J, g)
    return (
        projected_density_scale = maximum(x -> abs(real(x)), density),
        absolute_residual = maximum(
            I -> abs(real(density[I]) + sub_dt * divergence[I]), eachindex(density),
        ),
        density_scale = max(maximum(abs, state.rho_n), maximum(abs, state.rho_np1)),
        current_scale = maximum(a -> maximum(abs, a), state.J),
    )
end

function reference(species, c)
    return setprecision(BigFloat, 512) do
        C = BigFloat(c)
        energy = BigFloat(0)
        momentum = [BigFloat(0) for _ = 1:3]
        for ps in species, p in eachindex(ps.weight)
            v = ntuple(d -> BigFloat(ps.v[d][p]), 3)
            gamma = C / sqrt(C^2 - sum(x -> x^2, v))
            source = BigFloat(ps.m) * BigFloat(ps.weight[p])
            energy += source * (gamma - 1) * C^2
            for d = 1:3
                momentum[d] += source * gamma * v[d]
            end
        end
        return energy, Tuple(momentum)
    end
end

@testset "Relativistic collection diagnostics retain physical original factors" begin
    for T in (Float32, Float64)
        u = nextfloat(zero(T))
        tiny = [ParticleSet{1,T}(1; m = 1) for _ = 1:3]
        for ps in tiny
            ps.x[1][1], ps.v[1][1], ps.weight[1] = zero(T), T(3), u
        end
        # c=5,v=3 gives exact gamma=5/4: 3*(25/4)*u -> 19u,
        # and 3*(15/4)*u -> 11u. Per-population rounding gives different sums.
        @test kinetic_energy_relativistic(tiny, T(5)) == 19u
        @test total_momentum_relativistic(tiny, T(5)) == (11u, zero(T), zero(T))
        # The EFT owner requires nearest/gradual arithmetic. Preserve original
        # subnormal factors through exact replay when the task enables FTZ.
        saved_flush = get_zero_subnormals()
        flushed = try
            set_zero_subnormals(true)
            (kinetic_energy_relativistic(tiny, T(5)), total_momentum_relativistic(tiny, T(5)))
        finally
            set_zero_subnormals(saved_flush)
        end
        U = T === Float32 ? UInt32 : UInt64
        expected_tiny = (reinterpret(T, U(19)), (reinterpret(T, U(11)), zero(T), zero(T)))
        @test flushed === expected_tiny
        signs = [ParticleSet{1,T}(1; m = m) for m in (floatmax(T), floatmax(T), one(T))]
        for (i, ps) in enumerate(signs)
            ps.x[1][1], ps.v[1][1] = zero(T), i == 3 ? T(3) : T(i == 1 ? 4 : -4)
        end
        # Opposite overflowing momenta cancel before the remaining 15/4 term.
        @test total_momentum_relativistic(signs, T(5)) == (T(3.75), zero(T), zero(T))
        signs[2].m = -floatmax(T)
        signs[2].v[1][1] = T(4)
        @test kinetic_energy_relativistic(signs, T(5)) == T(6.25)
        power = T === Float32 ? 60 : 500
        c = ldexp(T(5), -power)
        rescued = ParticleSet{1,T}(1; m = floatmax(T))
        rescued.x[1][1], rescued.v[1][1], rescued.weight[1] = zero(T), ldexp(T(3), -power), T(2)
        expected_energy, expected_momentum = reference([rescued], c)
        @test isfinite(T(expected_energy)) && isfinite(T(expected_momentum[1]))
        @test kinetic_energy_relativistic([rescued], c) === T(expected_energy)
        @test total_momentum_relativistic(rescued, c)[1] === T(expected_momentum[1])
        fast = ParticleSet{1,T}(1)
        fast.x[1][1], fast.v[1][1] = zero(T), T === Float32 ? T(0.9997) : T(0.99999997)
        expected_energy, expected_momentum = reference([fast], one(T))
        native_energy = kinetic_energy_relativistic(fast, one(T))
        # One particle with m=c=weight=1 isolates the scalar contribution's
        # explicit four-adjacent-spacing retention contract. Compare against
        # the wide physical reference, including its unrounded low bits.
        spacing = min(
            BigFloat(native_energy) - BigFloat(prevfloat(native_energy)),
            BigFloat(nextfloat(native_energy)) - BigFloat(native_energy),
        )
        @test abs(BigFloat(native_energy) - expected_energy) <= 4spacing
        @test kinetic_energy_relativistic([fast], one(T)) === T(expected_energy)
        @test total_momentum_relativistic([fast], one(T))[1] === T(expected_momentum[1])
        flushed_energy = try
            set_zero_subnormals(true)
            kinetic_energy_relativistic(fast, one(T))
        finally
            set_zero_subnormals(saved_flush)
        end
        @test flushed_energy === T(expected_energy)
        @test total_momentum_relativistic(ParticleSet{1,T}[], one(T)) == (zero(T), zero(T), zero(T))
        fast.v[1][1] = one(T)
        @test_throws ArgumentError total_momentum_relativistic(fast, one(T))
    end
end

@testset "Relativistic collection fast phases and canonical legacy mover" begin
    for T in (Float32, Float64), D = 1:2
        state, species = setup(T, D)
        fine, fine_species = setup(T, D)
        source_density = sum(ps -> abs(ps.q) * sum(ps.weight), species) / prod(state.g.L)
        for (dt, ns) in ((T(0.03), 3), (T(0.02), 2), (T(0.012), 4))
            state.n_sub = ns
            step_empic!(state, species, dt)
            for _ = 1:ns
                step_empic!(fine, fine_species, dt / ns)
            end
            @test physical(state, species) == physical(fine, fine_species)
            @test state.time[] == fine.time[]
            scales = continuity_observation(state, dt / ns)
            println((observation = :relativistic_collection_continuity,
                T = T, D = D, dt = dt, n_sub = ns,
                relative_residual = charge_conservation_residual(state, dt / ns), scales = scales))
            # Uniform lattices have zero physical charge increment. Normalize
            # absolute roundoff by the population sources before cancellation.
            @test scales.absolute_residual < eps(T) * source_density
        end
        diagnostic = electromagnetic_diagnostics(state, species)
        @test diagnostic.kinetic == kinetic_energy_relativistic(species, state.c)
        @test diagnostic.momentum == total_momentum_relativistic(species, state.c)
        @test diagnostic.field == em_field_energy(state)
        @test diagnostic.velocity_time_offsets ==
              [-state.last_sub_dt[] / 2, -state.last_sub_dt[] / 2, zero(T)]
        # Normalized two-species compatibility uses the same prime/kick owner.
        e, ions = deepcopy(species[1]), deepcopy(species[2])
        e.q, e.m, ions.m = -one(T), one(T), T(7)
        fill!(ions.weight, one(T))
        legacy = EMPIC(
            FourierGrid(state.g.n, state.g.L),
            nparticles(e);
            mobile = true,
            mi = T(7),
            c = T(5),
            relativistic = true,
        )
        copied = deepcopy([e, ions])
        collection = EMPIC(FourierGrid(state.g.n, state.g.L), copied; c = T(5), relativistic = true)
        init_empic!(legacy, e, ions)
        init_empic!(collection, copied)
        for dt in T[0.01, 0.02, 0, 0.007]
            step_empic!(legacy, e, ions, dt)
            step_empic!(collection, copied, dt)
            @test physical(legacy, [e, ions]) == physical(collection, copied)
        end
    end
end

@testset "Relativistic collection continuity with a resolved charge increment" begin
    for T in (Float32, Float64), D = 1:2
        state, species = setup(T, D)
        # Preserve total charge and positive weights while creating a nonuniform
        # mobile density; its changing charge resolves the relative-error norm.
        species[1].weight[1] += T(0.25)
        species[1].weight[2] -= T(0.25)
        init_empic!(state, species)
        source_density = sum(ps -> abs(ps.q) * sum(ps.weight), species) / prod(state.g.L)
        for (dt, ns) in ((T(0.03), 3), (T(0.02), 2), (T(0.012), 4))
            state.n_sub = ns
            step_empic!(state, species, dt)
            scales = continuity_observation(state, dt / ns)
            println((observation = :relativistic_collection_resolved_continuity,
                T = T, D = D, dt = dt, n_sub = ns,
                relative_residual = charge_conservation_residual(state, dt / ns), scales = scales))
            @test scales.projected_density_scale > 64eps(T) * source_density
            @test charge_conservation_residual(state, dt / ns) < 1024eps(T)
        end
    end
end

@testset "Relativistic collection C1 failure retains the complete public state" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(4),))
        species = [ParticleSet{1,T}(1; m) for m in (one(T), nextfloat(zero(T)))]
        for ps in species
            ps.x[1][1], ps.weight[1] = zero(T), zero(T)
        end
        state = EMPIC(g, species; relativistic = true, c = one(T))
        init_empic!(state, species)
        fill!(state.E[1], one(T))
        before = snapshot(state, species)
        @test_throws ArgumentError step_empic!(state, species, one(T))
        @test snapshot(state, species) == before
        @test state.collection.staged[1].v != species[1].v
        species[end].v[1][1] = one(T)
        invalid = snapshot(state, species)
        @test_throws ArgumentError init_empic!(state, species)
        @test snapshot(state, species) == invalid
    end
end

@testset "Relativistic collection checkpoint/archive restores saved mode and c" begin
    for T in (Float32, Float64), mode in (false, true)
        state, species = setup(T, 1; relativistic = mode, n_sub = 2)
        # Saved velocities can exceed the destination's old c, and nonrelativistic
        # saved velocities may exceed the saved c. Validate the saved mode/c.
        fill!(species[1].v[1], mode ? T(4) : T(6))
        init_empic!(state, species)
        destination = [ParticleSet{1,T}(0) for _ = 1:3]
        restored = EMPIC(
            FourierGrid(state.g.n, state.g.L),
            destination;
            mobile = Bool[true, true, false],
            relativistic = !mode,
            c = one(T),
        )
        mktempdir() do dir
            path = joinpath(dir, "rel.chk")
            save_checkpoint(path, state, species)
            load_checkpoint!(restored, destination, path)
            @test restored.relativistic === mode
            @test restored.c === state.c
            @test physical(restored, destination) == physical(state, species)
            for (dt, ns) in ((T(0.01), 2), (zero(T), 3), (T(0.007), 3))
                state.n_sub = restored.n_sub = ns
                step_empic!(state, species, dt)
                step_empic!(restored, destination, dt)
                @test physical(restored, destination) == physical(state, species)
                @test (
                    restored.time[],
                    restored.step[],
                    restored.last_dt[],
                    restored.last_sub_dt[],
                ) == (state.time[], state.step[], state.last_dt[], state.last_sub_dt[])
            end
            archive = joinpath(dir, "rel.run")
            archive_run(archive, state, species; rng_seed = 19)
            payload = load_archive(archive).state
            H._restore_checkpoint_state!(restored, destination, payload)
            @test physical(restored, destination) == physical(state, species)
            restored.mobile = true
            invalid = snapshot(restored, destination)
            @test_throws ArgumentError H._restore_checkpoint_state!(restored, destination, payload)
            @test snapshot(restored, destination) == invalid
            @test restored.mobile
            restored.mobile = false
            before = snapshot(restored, destination)
            malformed = merge(deepcopy(payload), (relativistic = 1,))
            @test_throws ArgumentError H._restore_checkpoint_state!(
                restored,
                destination,
                malformed,
            )
            @test snapshot(restored, destination) == before
            if !mode
                malformed = merge(deepcopy(payload), (relativistic = true,))
                @test_throws ArgumentError H._restore_checkpoint_state!(
                    restored,
                    destination,
                    malformed,
                )
                @test snapshot(restored, destination) == before
            end
        end
    end
end
end # module
