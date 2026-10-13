# IO-001 — checkpoint/restart. A run interrupted and restarted from a checkpoint
# must continue bitwise-identically to the uninterrupted run (stepping is
# deterministic: no RNG, serial loops).

using HybridPlasmaPIC, Test, Random, Serialization

struct CheckpointArrayElectrons <: HybridPlasmaPIC.ElectronClosure
    temperature::Vector{Float64}
    alias::Vector{Float64}
end

function HybridPlasmaPIC.electron_pressure!(pe, n, c::CheckpointArrayElectrons)
    @. pe = c.temperature[1] * n
    return pe
end

mutable struct CheckpointParameterizedShape <: ShapeFunction
    parameters::Vector{Float64}
end

function build_camcl_checkpoint_run(::Type{T}, ::Val{D}) where {T,D}
    n = ntuple(d -> d == 1 ? 8 : 4, D)
    L = ntuple(_ -> T(2π), D)
    counts = ntuple(d -> 2n[d], D)
    g = FourierGrid(n, L)
    ps = ParticleSet{D,T}(prod(counts))
    load_lattice!(ps, ntuple(_ -> zero(T), D), L, counts)
    set_density_weight!(ps, one(T), g)
    ps.v[1] .= T(0.03) .* sin.(ps.x[1])
    ps.v[2] .= T(0.01) .* cos.(ps.x[1])
    st = CAMCLStepper(g, HybridModel(IsothermalElectrons(T(0.1))), CIC(), nparticles(ps))
    fill!(st.fields.B[1], one(T))
    for I in CartesianIndices(n)
        x = T(I[1] - 1) * g.dx[1]
        st.fields.B[2][I] = T(0.05) * sin(x)
        st.fields.B[3][I] = T(0.02) * cos(x)
    end
    init_camcl!(st, ps)
    return ps, st
end

@testset "CAM-CL checkpoint preserves persistent field and particle phases" begin
    for T in (Float32, Float64),
        D in (1, 2, 3),
        scenario in (:initial, :live, :reset, :edited, :inactive)

        ps, st = build_camcl_checkpoint_run(T, Val(D))
        dt = T(0.02)
        if scenario != :initial
            for _ = 1:3
                step_camcl!(st, ps, dt; NB = 2)
            end
            @test st.cl_live[]
        end
        if scenario == :reset
            init_camcl!(st, ps)
            @test st.cl_h[] > zero(T)
        elseif scenario == :edited
            st.fields.B[2][1] += T(0.001)
            @test st.fields.B[2] != st.Bout[2]
        elseif scenario == :inactive
            st.cl_live[] = false
        end
        mktemp() do path, io
            close(io)
            @test save_checkpoint(path, st, ps) == path
            grid = FourierGrid(st.g.n, st.g.L)
            restored = CAMCLStepper(grid, st.model, st.shape, 1)
            particles = ParticleSet{D,T}(1)
            @test load_checkpoint!(restored, particles, path) === restored
            @test HybridPlasmaPIC._checkpoint_state(restored, particles) ==
                  HybridPlasmaPIC._checkpoint_state(st, ps)
            @test length(restored.work) == nparticles(ps)
            # Equal substeps first exercise pair persistence; changes afterward
            # exercise the same rebootstrap decisions as uninterrupted evolution.
            for (next_dt, NB) in ((dt, 2), (dt / 2, 2), (dt / 2, 4), (dt, 4))
                step_camcl!(st, ps, next_dt; NB)
                step_camcl!(restored, particles, next_dt; NB)
                @test HybridPlasmaPIC._checkpoint_state(restored, particles) ==
                      HybridPlasmaPIC._checkpoint_state(st, ps)
            end
        end
    end
end

@testset "CAM-CL checkpoint validates persistent state before mutation" begin
    ps, st = build_camcl_checkpoint_run(Float64, Val(1))
    step_camcl!(st, ps, 0.02; NB = 2)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, st, ps)
        state = deserialize(path)
        destination = CAMCLStepper(st.g, HybridModel(IsothermalElectrons(2.0)), CIC(), 1)
        particles = ParticleSet{1,Float64}(1)
        fill!(particles.x[1], 0.75)
        before = deepcopy(HybridPlasmaPIC._checkpoint_state(destination, particles))
        missing_phase = (; (k => v for (k, v) in pairs(state) if k != :cl_h)...)
        malformed = (
            merge(state, (Blead = (state.Blead[1], state.Blead[2]),)),
            merge(state, (Blag = ([1.0], state.Blag[2], state.Blag[3]),)),
            merge(state, (Bout = map(x -> Float32.(x), state.Bout),)),
            merge(state, (cl_h = NaN,)),
            merge(state, (cl_h = -0.01,)),
            merge(state, (cl_h = 0.01f0,)),
            merge(state, (cl_h = 0.0, cl_live = true)),
            merge(state, (cl_live = 1,)),
            missing_phase,
        )
        for bad in malformed
            serialize(path, bad)
            @test_throws ArgumentError load_checkpoint!(destination, particles, path)
            @test HybridPlasmaPIC._checkpoint_state(destination, particles) == before
            @test length(destination.work) == 1
        end
    end
end

@testset "checkpoints reject switching the field integrator" begin
    ps, camcl = build_camcl_checkpoint_run(Float64, Val(1))
    hybrid = HybridStepper(camcl.g, camcl.model, camcl.shape, nparticles(ps))
    for (source, target) in ((camcl, hybrid), (hybrid, camcl))
        before = deepcopy(HybridPlasmaPIC._checkpoint_state(target, ps))
        mktemp() do path, io
            close(io)
            save_checkpoint(path, source, ps)
            @test_throws ArgumentError load_checkpoint!(target, ps, path)
            @test HybridPlasmaPIC._checkpoint_state(target, ps) == before
        end
    end
end

function build_run(seed)
    T = Float64
    n = 16
    L = 2π
    g = FourierGrid((n,), (L,))
    nppc = 50
    N = nppc * n
    ps = ParticleSet{1,T}(N)
    load_lattice_1d!(ps, 0.0, L)
    set_density_weight!(ps, 1.0, g)
    load_quiet_velocities!(ps, MersenneTwister(seed), (0.0, 0.0, 0.0), (0.1, 0.1, 0.1))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), N)
    fill!(st.fields.B[1], 1.0)
    x = [(i - 1) * g.dx[1] for i = 1:n]
    st.fields.B[2] .= 0.01 .* cos.(x)
    init!(st, ps)
    return ps, st
end

@testset "IO-001 checkpoint/restart bitwise" begin
    dt = 0.02
    NB = 2

    # reference: 20 uninterrupted steps
    psr, str = build_run(42)
    for _ = 1:20
        step!(str, psr, dt; NB)
    end

    # restart: 10 steps, checkpoint, restore into fresh state, 10 more steps
    psa, sta = build_run(42)
    for _ = 1:10
        step!(sta, psa, dt; NB)
    end
    path = tempname()
    save_checkpoint(path, sta, psa)

    psb, stb = build_run(42)           # any state; load overwrites it
    load_checkpoint!(stb, psb, path)
    @test stb.step[] == 10             # restored counters
    for _ = 1:10
        step!(stb, psb, dt; NB)
    end
    rm(path; force = true)

    # bitwise identical to the uninterrupted run
    @test psb.x[1] == psr.x[1]
    for c = 1:3
        @test psb.v[c] == psr.v[c]
        @test stb.fields.B[c] == str.fields.B[c]
        @test stb.fields.E[c] == str.fields.E[c]
    end
    @test psb.id == psr.id
    @test stb.step[] == str.step[] == 20
end

@testset "load_checkpoint! resizes particle workspaces" begin
    g = FourierGrid((8,), (1.0,))
    sta = HybridStepper(g, HybridModel(IsothermalElectrons(0.1)), CIC(), 6)
    psa = ParticleSet{1,Float64}(6)
    load_lattice_1d!(psa, 0.0, 1.0)
    set_density_weight!(psa, 1.0, g)
    init!(sta, psa)

    path = tempname()
    save_checkpoint(path, sta, psa)

    stb = HybridStepper(g, HybridModel(IsothermalElectrons(0.1)), CIC(), 1)
    psb = ParticleSet{1,Float64}(1)
    load_checkpoint!(stb, psb, path)
    @test nparticles(psb) == 6
    @test length(stb.work) == 6
    @test all(length(stb.Ep[c]) == 6 for c = 1:3)
    @test all(length(stb.Bp[c]) == 6 for c = 1:3)
    @test length(stb.xmid[1]) == 6
    @test step!(stb, psb, 0.01) === stb
    rm(path; force = true)
end

@testset "checkpoint restores saved physics configuration" begin
    for (source_closure, target_closure) in (
        (IsothermalElectrons(0.1), IsothermalElectrons(2.0)),
        (PolytropicElectrons(0.1, 1.0, 5 / 3), PolytropicElectrons(2.0, 2.0, 1.4)),
        (CGLElectrons(0.1, 0.1, 1.0, 1.0), CGLElectrons(2.0, 3.0, 1.0, 1.0)),
        (CheckpointArrayElectrons([0.1], [0.1]), CheckpointArrayElectrons([2.0], [2.0])),
    )
        g = FourierGrid((8,), (1.0,))
        model = HybridModel(source_closure; η = 0.02, ηH = 0.001, nfloor = 1e-5, de2 = 0.1)
        a = HybridStepper(g, model, CIC(), 64)
        p = ParticleSet{1,Float64}(64)
        load_lattice_1d!(p, 0.0, 1.0)
        set_density_weight!(p, 1.0, g)
        p.weight .*= 1 .+ 0.1cos.(2π .* p.x[1])
        fill!(a.fields.B[1], 1.0)
        init!(a, p)
        step!(a, p, 0.005)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, a, p)
            # Different resistivity/inertia/floor must not survive restoration.
            b = HybridStepper(g, HybridModel(target_closure), CIC(), 1)
            q = ParticleSet{1,Float64}(1)
            load_checkpoint!(b, q, path)
            @test b.model.η == a.model.η
            @test b.model.ηH == a.model.ηH
            @test b.model.nfloor == a.model.nfloor
            @test b.model.de2 == a.model.de2
            step!(a, p, 0.005)
            step!(b, q, 0.005)
            @test p.x == q.x
            @test p.v == q.v
            @test a.fields.B == b.fields.B
            @test a.fields.E == b.fields.E
        end
    end
end

@testset "checkpoint owns custom mutable configuration snapshots" begin
    temperature = [0.1]
    closure = CheckpointArrayElectrons(temperature, temperature)
    shape = CheckpointParameterizedShape([1.0, 2.0])
    g = FourierGrid((8,), (1.0,))
    st = HybridStepper(g, HybridModel(closure), shape, 2)
    ps = ParticleSet{1,Float64}(2)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, st, ps)
        temperature[1] = 9.0
        shape.parameters[1] = 8.0
        load_checkpoint!(st, ps, path)
        @test st.model.closure.temperature == [0.1]
        @test st.model.closure.temperature === st.model.closure.alias
        @test st.model.closure.temperature !== temperature
        @test st.shape.parameters == [1.0, 2.0]
        @test st.shape !== shape
    end
end

@testset "checkpoint preserves variable-step velocity time level" begin
    reference_particles, reference = build_run(7)
    particles, stepper = build_run(7)
    steps = (0.01, 0.015, 0.005, 0.012)
    for dt in steps
        step!(reference, reference_particles, dt; NB = 2)
    end
    for dt in steps[1:2]
        step!(stepper, particles, dt; NB = 2)
    end
    mktemp() do path, io
        close(io)
        save_checkpoint(path, stepper, particles)
        restored_particles, restored = build_run(8)
        load_checkpoint!(restored, restored_particles, path)
        @test restored.last_dt[] == steps[2]
        for dt in steps[3:4]
            step!(restored, restored_particles, dt; NB = 2)
        end
        @test restored_particles.x == reference_particles.x
        @test restored_particles.v == reference_particles.v
        @test restored.fields.B == reference.fields.B
        @test restored.fields.E == reference.fields.E
        @test restored.last_dt[] == steps[end]
    end
end

@testset "checkpoint restores writable views without resizing" begin
    source_particles, source = build_run(7)
    step!(source, source_particles, 0.01)
    destination_particles, destination = build_run(8)
    viewed = ParticleSet{1,Float64}(
        map(a -> view(a, :), destination_particles.x),
        map(a -> view(a, :), destination_particles.v),
        view(destination_particles.weight, :),
        view(destination_particles.id, :),
        view(destination_particles.tag, :),
        destination_particles.q,
        destination_particles.m,
    )
    mktemp() do path, io
        close(io)
        save_checkpoint(path, source, source_particles)
        @test load_checkpoint!(destination, viewed, path) === destination
        @test destination_particles.x == source_particles.x
        @test destination_particles.v == source_particles.v
        @test destination_particles.weight == source_particles.weight
        @test destination_particles.id == source_particles.id
        @test destination_particles.tag == source_particles.tag
        step!(source, source_particles, 0.015)
        step!(destination, viewed, 0.015)
        @test destination_particles.x == source_particles.x
        @test destination_particles.v == source_particles.v
        @test destination.fields.B == source.fields.B
        @test destination.fields.E == source.fields.E
    end
end

@testset "checkpoint preflights every destination resize before mutation" begin
    source_particles, source = build_run(2)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, source, source_particles)
        for view_field in (:x, :v, :weight, :id, :tag)
            model = HybridModel(IsothermalElectrons(2.0))
            destination = HybridStepper(source.g, model, CIC(), 1)
            p = ParticleSet{1,Float64}(1)
            p.x[1][1] = 0.75
            particles = ParticleSet{1,Float64}(
                view_field == :x ? map(a -> view(a, :), p.x) : p.x,
                view_field == :v ? map(a -> view(a, :), p.v) : p.v,
                view_field == :weight ? view(p.weight, :) : p.weight,
                view_field == :id ? view(p.id, :) : p.id,
                view_field == :tag ? view(p.tag, :) : p.tag,
                p.q,
                p.m,
            )
            @test_throws ArgumentError load_checkpoint!(destination, particles, path)
            @test particles.x[1] == [0.75]
            @test particles.v == ([0.0], [0.0], [0.0])
            @test particles.weight == [1.0]
            @test particles.id == UInt64[1]
            @test particles.tag == UInt32[0]
            @test destination.model === model
            @test length(destination.work) == 1
            @test destination.step[] == 0
            @test destination.last_dt[] == 0.0
        end
    end
end

@testset "checkpoint rejects incompatible configuration types before mutation" begin
    ps, st = build_run(2)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, st, ps)
        for (model, shape) in (
            (HybridModel(IsothermalElectrons(0.5)), TSC()),
            (HybridModel(PolytropicElectrons(0.5, 1.0, 5 / 3)), CIC()),
        )
            target = HybridStepper(st.g, model, shape, 1)
            particles = ParticleSet{1,Float64}(1)
            particles.x[1][1] = 0.75
            target.fields.B[1] .= 7.0
            @test_throws ArgumentError load_checkpoint!(target, particles, path)
            @test particles.x[1] == [0.75]
            @test all(==(7.0), target.fields.B[1])
            @test target.model === model
            @test target.shape === shape
            @test length(target.work) == 1
        end
    end
end

@testset "load_checkpoint! rejects eltype mismatch" begin
    # Regression: only D and grid were validated; a Float64 checkpoint loaded into
    # a Float32 stepper was silently converted, breaking the bitwise guarantee.
    ga = FourierGrid((16,), (1.0,))
    sta = HybridStepper(ga, HybridModel(IsothermalElectrons(0.5)), CIC(), 32)
    psa = ParticleSet{1,Float64}(32)
    psa.x[1] .= 0.0
    path = tempname()
    save_checkpoint(path, sta, psa)
    gb = FourierGrid((16,), (1.0f0,))
    stb = HybridStepper(gb, HybridModel(IsothermalElectrons(0.5f0)), CIC(), 32)
    psb = ParticleSet{1,Float32}(32)
    @test_throws ArgumentError load_checkpoint!(stb, psb, path)
    rm(path; force = true)
end

@testset "load_checkpoint! rejects box-length mismatch" begin
    ga = FourierGrid((16,), (1.0,))
    sta = HybridStepper(ga, HybridModel(IsothermalElectrons(0.5)), CIC(), 32)
    psa = ParticleSet{1,Float64}(32)
    psa.x[1] .= 0.0
    path = tempname()
    save_checkpoint(path, sta, psa)
    gb = FourierGrid((16,), (2.0,))
    stb = HybridStepper(gb, HybridModel(IsothermalElectrons(0.5)), CIC(), 32)
    psb = ParticleSet{1,Float64}(32)
    @test_throws ArgumentError load_checkpoint!(stb, psb, path)
    rm(path; force = true)
end

@testset "load_checkpoint! rejects checkpoints missing box lengths" begin
    g = FourierGrid((16,), (1.0,))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), 32)
    ps = ParticleSet{1,Float64}(32)
    legacy = (
        integrator = :rk4,
        D = 1,
        T = Float64,
        ncell = st.g.n,
        model = st.model,
        shape = st.shape,
        x = ps.x,
        v = ps.v,
        weight = ps.weight,
        id = ps.id,
        tag = ps.tag,
        q = ps.q,
        m = ps.m,
        B = st.fields.B,
        E = st.fields.E,
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
    )
    path = tempname()
    Serialization.serialize(path, legacy)
    @test_throws ArgumentError load_checkpoint!(st, ps, path)
    rm(path; force = true)
end

@testset "load_checkpoint! rejects inconsistent particle array lengths" begin
    g = FourierGrid((16,), (1.0,))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), 2)
    ps = ParticleSet{1,Float64}(2)
    fill!(ps.x[1], 0.25)
    fill!(st.fields.B[1], 7.0)
    legacy = (
        integrator = :rk4,
        D = 1,
        T = Float64,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        shape = st.shape,
        x = ([0.1],),
        v = ([0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
        weight = [1.0, 1.0],
        id = UInt64[1, 2],
        tag = UInt32[0, 0],
        q = ps.q,
        m = ps.m,
        B = (copy(st.fields.B[1]), copy(st.fields.B[2]), copy(st.fields.B[3])),
        E = (copy(st.fields.E[1]), copy(st.fields.E[2]), copy(st.fields.E[3])),
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
    )
    path = tempname()
    Serialization.serialize(path, legacy)
    @test_throws ArgumentError load_checkpoint!(st, ps, path)
    @test length(ps.x[1]) == 2
    @test ps.x[1] == fill(0.25, 2)
    @test all(==(7.0), st.fields.B[1])
    rm(path; force = true)
end

@testset "load_checkpoint! rejects truncated field arrays" begin
    g = FourierGrid((16,), (1.0,))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), 2)
    ps = ParticleSet{1,Float64}(2)
    fill!(st.fields.B[1], 9.0)
    fill!(st.fields.E[1], 8.0)
    legacy = (
        integrator = :rk4,
        D = 1,
        T = Float64,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        shape = st.shape,
        x = ([0.1, 0.2],),
        v = ([0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
        weight = [1.0, 1.0],
        id = UInt64[1, 2],
        tag = UInt32[0, 0],
        q = ps.q,
        m = ps.m,
        B = ([1.5], copy(st.fields.B[2]), copy(st.fields.B[3])),
        E = ([4.5], copy(st.fields.E[2]), copy(st.fields.E[3])),
        time = st.time[],
        step = st.step[],
        last_dt = st.last_dt[],
    )
    path = tempname()
    Serialization.serialize(path, legacy)
    @test_throws ArgumentError load_checkpoint!(st, ps, path)
    @test all(==(9.0), st.fields.B[1])
    @test all(==(8.0), st.fields.E[1])
    rm(path; force = true)
end

@testset "load_checkpoint! rejects scalar/type conversion before mutation" begin
    g = FourierGrid((16,), (1.0,))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), 2)
    ps = ParticleSet{1,Float64}(2; q = 2.0, m = 3.0)
    fill!(ps.x[1], 0.25)
    fill!(st.fields.B[1], 7.0)
    st.step[] = 3
    st.time[] = 0.75
    base = (
        integrator = :rk4,
        D = 1,
        T = Float64,
        ncell = st.g.n,
        L = st.g.L,
        model = st.model,
        shape = st.shape,
        x = ([0.1, 0.2],),
        v = ([1.0, 1.0], [2.0, 2.0], [3.0, 3.0]),
        weight = [4.0, 4.0],
        id = UInt64[10, 11],
        tag = UInt32[1, 1],
        q = 2.0,
        m = 3.0,
        B = (fill(8.0, size(st.fields.B[1])), copy(st.fields.B[2]), copy(st.fields.B[3])),
        E = (copy(st.fields.E[1]), copy(st.fields.E[2]), copy(st.fields.E[3])),
        time = 1.25,
        step = 4,
        last_dt = 0.25,
    )

    function assert_rejects_without_mutation(bad)
        fill!(ps.x[1], 0.25)
        ps.q = 2.0
        ps.m = 3.0
        fill!(st.fields.B[1], 7.0)
        st.step[] = 3
        st.time[] = 0.75
        path = tempname()
        try
            Serialization.serialize(path, bad)
            @test_throws ArgumentError load_checkpoint!(st, ps, path)
            @test ps.x[1] == fill(0.25, 2)
            @test ps.q == 2.0
            @test ps.m == 3.0
            @test all(==(7.0), st.fields.B[1])
            @test st.step[] == 3
            @test st.time[] == 0.75
        finally
            rm(path; force = true)
        end
    end

    assert_rejects_without_mutation(merge(base, (step = "bad",)))
    assert_rejects_without_mutation(merge(base, (q = 2,)))
    assert_rejects_without_mutation(merge(base, (last_dt = 0.0,)))
    assert_rejects_without_mutation(merge(base, (last_dt = NaN,)))
    assert_rejects_without_mutation(merge(base, (last_dt = Inf,)))
    assert_rejects_without_mutation(merge(base, (last_dt = -0.25,)))
    assert_rejects_without_mutation(merge(base, (last_dt = 0.25f0,)))
    assert_rejects_without_mutation(merge(base, (step = 0,)))
    assert_rejects_without_mutation(merge(base, (step = -1,)))
    assert_rejects_without_mutation(merge(base, (x = (Float32[0.1, 0.2],),)))
    assert_rejects_without_mutation(merge(base, (B = (Float32.(base.B[1]), base.B[2], base.B[3]),)))
end

@testset "load_checkpoint! rejects invalid container shape before mutation" begin
    g = FourierGrid((16,), (1.0,))
    st = HybridStepper(g, HybridModel(IsothermalElectrons(0.5)), CIC(), 2)
    ps = ParticleSet{1,Float64}(2)
    fill!(ps.x[1], 0.25)
    fill!(st.fields.B[1], 7.0)
    st.step[] = 3
    path = tempname()
    Serialization.serialize(path, 7)
    @test_throws ArgumentError load_checkpoint!(st, ps, path)
    @test ps.x[1] == fill(0.25, 2)
    @test all(==(7.0), st.fields.B[1])
    @test st.step[] == 3
    rm(path; force = true)
end

include("test_io_electrostatic.jl")
include("test_io_electromagnetic.jl")
include("test_io_hall_mhd.jl")
include("test_io_collections.jl")
include("test_io_particle_storage.jl")
include("test_io_clock.jl")
include("test_io_configuration_ownership.jl")
include("test_io_publication.jl")
