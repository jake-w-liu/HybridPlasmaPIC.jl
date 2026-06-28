# IO-001 — checkpoint/restart. A run interrupted and restarted from a checkpoint
# must continue bitwise-identically to the uninterrupted run (stepping is
# deterministic: no RNG, serial loops).

using HybridPlasmaPIC, Test, Random, Serialization

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
        D = 1,
        T = Float64,
        ncell = st.g.n,
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
        D = 1,
        T = Float64,
        ncell = st.g.n,
        L = st.g.L,
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
        D = 1,
        T = Float64,
        ncell = st.g.n,
        L = st.g.L,
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
        D = 1,
        T = Float64,
        ncell = st.g.n,
        L = st.g.L,
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
    assert_rejects_without_mutation(merge(base, (x = (Float32[0.1, 0.2],),)))
    assert_rejects_without_mutation(
        merge(base, (B = (Float32.(base.B[1]), base.B[2], base.B[3]),)),
    )
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
