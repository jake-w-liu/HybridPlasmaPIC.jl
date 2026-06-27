# IO-001 — checkpoint/restart. A run interrupted and restarted from a checkpoint
# must continue bitwise-identically to the uninterrupted run (stepping is
# deterministic: no RNG, serial loops).

using HybridPlasmaPIC, Test, Random

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
