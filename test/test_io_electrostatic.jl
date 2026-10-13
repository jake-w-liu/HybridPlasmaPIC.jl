using HybridPlasmaPIC, Test, Serialization

function io_es_run(State, ::Type{T}, ::Val{D}) where {T,D}
    n = ntuple(_ -> 4, D)
    L = ntuple(_ -> T(2π), D)
    counts = ntuple(_ -> 8, D)
    grid = FourierGrid(n, L)
    electrons = ParticleSet{D,T}(prod(counts); q = -one(T), m = one(T))
    load_lattice!(electrons, ntuple(_ -> zero(T), D), L, counts)
    set_density_weight!(electrons, one(T), grid)
    electrons.weight .*= one(T) .+ T(0.03) .* cos.(electrons.x[1])
    for c = 1:3
        electrons.v[c] .= T(0.01c) .* sin.(electrons.x[1])
    end
    state = State(grid, nparticles(electrons))
    init_espic!(state, electrons)
    return state, electrons
end

function io_es_bytes(state, electrons)
    io = IOBuffer()
    serialize(io, HybridPlasmaPIC._checkpoint_state(state, electrons))
    return take!(io)
end

io_es_workspace_size(state) = state isa Electrostatic1D ? length(state.Ep) : length(state.Ep[1])

@testset "electrostatic checkpoint preserves full state and velocity phase" begin
    for T in (Float32, Float64),
        (State, D) in
        ((Electrostatic1D, 1), (ElectrostaticPIC, 1), (ElectrostaticPIC, 2), (ElectrostaticPIC, 3)),
        scenario in (:initial, :live, :reset, :edited)

        state, electrons = io_es_run(State, T, Val(D))
        if scenario != :initial
            step_espic!(state, electrons, T(0.013))
            step_espic!(state, electrons, T(0.021))
        end
        if scenario == :reset
            init_espic!(state, electrons)
        elseif scenario == :edited
            state.ne[2] += T(0.1)
            if state isa Electrostatic1D
                state.E[2] += T(0.001)
            else
                state.E[1][2] += T(0.001)
                state.E[3][2] += T(0.002)
            end
        end
        mktemp() do path, io
            close(io)
            @test save_checkpoint(path, state, electrons) == path
            restored = State(FourierGrid(state.g.n, state.g.L), 1; n0 = state.n0)
            particles = ParticleSet{D,T}(1; q = -one(T), m = one(T))
            @test load_checkpoint!(restored, particles, path) === restored
            @test io_es_bytes(restored, particles) == io_es_bytes(state, electrons)
            @test io_es_workspace_size(restored) == nparticles(electrons)
            # Include unchanged timesteps, changes of size, and an intervening no-op.
            for dt in T[0.02, 0.02, 0.007, 0.0, 0.03]
                step_espic!(state, electrons, dt)
                step_espic!(restored, particles, dt)
                @test io_es_bytes(restored, particles) == io_es_bytes(state, electrons)
            end
        end
    end
end

@testset "electrostatic checkpoint rejects incompatible destinations before mutation" begin
    state, electrons = io_es_run(ElectrostaticPIC, Float64, Val(1))
    step_espic!(state, electrons, 0.02)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, state, electrons)
        targets = (
            Electrostatic1D(state.g, 1),
            ElectrostaticPIC(state.g, 1; n0 = 2.0),
            ElectrostaticPIC(FourierGrid((8,), state.g.L), 1),
            ElectrostaticPIC(FourierGrid(state.g.n, (4π,)), 1),
            ElectrostaticPIC(FourierGrid((4, 4), (2π, 2π)), 1),
            ElectrostaticPIC(FourierGrid((4,), (Float32(2π),)), 1),
        )
        for target in targets
            T = eltype(target.ne)
            D = length(target.g.n)
            particles = ParticleSet{D,T}(1; q = -one(T), m = one(T))
            foreach(a -> fill!(a, T(0.75)), particles.x)
            before = io_es_bytes(target, particles)
            @test_throws ArgumentError load_checkpoint!(target, particles, path)
            @test io_es_bytes(target, particles) == before
            @test io_es_workspace_size(target) == 1
        end
        target = ElectrostaticPIC(state.g, 1)
        for particles in
            (ParticleSet{2,Float64}(1; q = -1.0), ParticleSet{1,Float32}(1; q = -1.0f0))
            foreach(a -> fill!(a, 0.75), particles.x)
            before = io_es_bytes(target, particles)
            @test_throws ArgumentError load_checkpoint!(target, particles, path)
            @test io_es_bytes(target, particles) == before
            @test io_es_workspace_size(target) == 1
        end
    end
end

@testset "electrostatic checkpoint validates field, species and phase data" begin
    for State in (Electrostatic1D, ElectrostaticPIC)
        state, electrons = io_es_run(State, Float64, Val(1))
        step_espic!(state, electrons, 0.02)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, state, electrons)
            saved = deserialize(path)
            target = State(state.g, 1)
            particles = ParticleSet{1,Float64}(1; q = -1.0)
            fill!(particles.x[1], 0.75)
            before = io_es_bytes(target, particles)
            malformed_E = State === Electrostatic1D ? (saved.E,) : saved.E[1]
            wrong_shape_E =
                State === Electrostatic1D ? zeros(3) : (zeros(3), saved.E[2], saved.E[3])
            wrong_type_E =
                State === Electrostatic1D ? Float32.(saved.E) :
                (Float32.(saved.E[1]), saved.E[2], saved.E[3])
            missing_phase = (; (k => v for (k, v) in pairs(saved) if k != :primed)...)
            malformed = (
                merge(saved, (integrator = :unsupported,)),
                merge(saved, (E = malformed_E,)),
                merge(saved, (E = wrong_shape_E,)),
                merge(saved, (E = wrong_type_E,)),
                merge(saved, (ne = zeros(3),)),
                merge(saved, (ne = Float32.(saved.ne),)),
                merge(saved, (q = 1.0,)),
                merge(saved, (m = 2.0,)),
                merge(saved, (q = -1,)),
                merge(saved, (n0 = NaN,)),
                merge(saved, (n0 = 1.0f0,)),
                merge(saved, (primed = false,)),
                merge(saved, (primed = 1,)),
                merge(saved, (last_dt = 0.0,)),
                merge(saved, (last_dt = NaN,)),
                merge(saved, (last_dt = -0.01,)),
                merge(saved, (last_dt = 0.02f0,)),
                missing_phase,
            )
            for bad in malformed
                serialize(path, bad)
                @test_throws ArgumentError load_checkpoint!(target, particles, path)
                @test io_es_bytes(target, particles) == before
                @test io_es_workspace_size(target) == 1
            end
            # Invalid electron species cannot replace an already valid checkpoint.
            serialize(path, saved)
            electrons.q = 1.0
            @test_throws ArgumentError save_checkpoint(path, state, electrons)
            @test deserialize(path) == saved
        end
    end
end

@testset "electrostatic checkpoint supports empty and view-backed particles" begin
    for State in (Electrostatic1D, ElectrostaticPIC)
        state, electrons = io_es_run(State, Float64, Val(1))
        particles = deepcopy(electrons)
        viewed = ParticleSet{1,Float64}(
            map(a -> view(a, :), particles.x),
            map(a -> view(a, :), particles.v),
            view(particles.weight, :),
            view(particles.id, :),
            view(particles.tag, :),
            particles.q,
            particles.m,
        )
        for a in (particles.x..., particles.v..., particles.weight)
            fill!(a, -0.75)
        end
        fill!(particles.id, 0)
        fill!(particles.tag, 7)
        target = State(state.g, 1)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, state, electrons)
            @test load_checkpoint!(target, viewed, path) === target
            @test io_es_bytes(target, particles) == io_es_bytes(state, electrons)
            short_views = ParticleSet{1,Float64}(
                map(a -> view(a, 1:1), particles.x),
                map(a -> view(a, 1:1), particles.v),
                view(particles.weight, 1:1),
                view(particles.id, 1:1),
                view(particles.tag, 1:1),
                particles.q,
                particles.m,
            )
            blocked = State(state.g, 1)
            before = io_es_bytes(blocked, short_views)
            @test_throws ArgumentError load_checkpoint!(blocked, short_views, path)
            @test io_es_bytes(blocked, short_views) == before
            @test io_es_workspace_size(blocked) == 1
            # A zero-particle state remains restartable without inventing particles.
            empty_particles = ParticleSet{1,Float64}(0; q = -1.0)
            init_espic!(state, empty_particles)
            save_checkpoint(path, state, empty_particles)
            load_checkpoint!(target, particles, path)
            @test nparticles(particles) == 0
            @test io_es_workspace_size(target) == 0
            @test io_es_bytes(target, particles) == io_es_bytes(state, empty_particles)
            @test step_espic!(target, particles, 0.02) === target
        end
    end
end
