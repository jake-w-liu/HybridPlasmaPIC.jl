using HybridPlasmaPIC, Test
using Random
using Serialization

# Build a small but real 1D3V hybrid run: grid, model, stepper, particles.
function _build_run(; N = 64, n = (16,), L = (4.0,), seed = 2026)
    g = FourierGrid(n, L)
    model = HybridModel(IsothermalElectrons(1.0); η = 0.0, nfloor = 1e-6)
    st = HybridStepper(g, model, CIC(), N)
    ps = ParticleSet{1,Float64}(N; q = 1.0, m = 1.0)
    rng = MersenneTwister(seed)
    load_lattice_1d!(ps, 0.0, L[1])
    load_maxwellian!(ps, rng, (0.5, 0.0, 0.0), (0.1, 0.1, 0.1))
    set_density_weight!(ps, 1.0, g)
    # A nonzero guide field so B and E are nontrivial.
    fill!(st.fields.B[3], 0.7)
    init!(st, ps)
    return g, model, st, ps
end

@testset "operators_match (Phase-1 migration check)" begin
    for n in ((16,), (32,))
        g = FourierGrid(n, (3.0,))
        @test operators_match(g) === true
    end
    g2 = FourierGrid((8, 8), (2.0, 3.0))
    @test operators_match(g2) === true
    g3 = FourierGrid((4, 5, 6), (1.0, 2.0, 3.0))
    @test operators_match(g3) === true

    # Independently confirm it is genuinely measuring identity to ~1e-14:
    # compute the max difference by hand on a fresh field.
    g = FourierGrid((24,), (5.0,))
    f = rand(MersenneTwister(7), Float64, 24)
    a = HybridPlasmaPIC.deriv(f, g, 1)
    g2b = FourierGrid((24,), (5.0,))
    b = SpectralOperators.deriv!(similar(f), f, g2b, 1)
    @test maximum(abs.(a .- b)) <= 1e-14
end

@testset "sample_particles (sampled dumps, §6)" begin
    _, _, _, ps = _build_run(; N = 100)
    s = sample_particles(ps, 10)
    # cld(100,10) = 10 entries: indices 1,11,21,...,91.
    @test length(s.index) == 10
    @test s.index == collect(1:10:100)
    @test length(s.x[1]) == 10
    @test all(length(v) == 10 for v in s.v)
    # Values match the source particles at the sampled indices.
    @test s.x[1] == ps.x[1][s.index]
    @test s.v[2] == ps.v[2][s.index]
    # Copy semantics: mutating the dump must not touch ps.
    s.x[1][1] = -999.0
    @test ps.x[1][1] != -999.0

    # ~N/10 for a non-divisor count too, and stride=1 returns all.
    _, _, _, ps2 = _build_run(; N = 64)
    @test length(sample_particles(ps2, 10).index) == cld(64, 10)
    @test length(sample_particles(ps2, 1).index) == 64
    @test_throws ArgumentError sample_particles(ps2, 0)
end

@testset "archive_run / load_archive round-trip" begin
    g, model, st, ps = _build_run()
    # Advance a couple of steps so time/step/fields are nonzero.
    step!(st, ps, 0.05)
    step!(st, ps, 0.05)

    seed = 31415
    mktemp() do path, io
        close(io)
        ret = archive_run(
            path,
            st,
            ps;
            rng_seed = seed,
            filter_desc = "exp_filter(36,8)",
            diagnostic_desc = "energy",
            rank_layout = "mpi;ranks=2;dims=(2);coords=(1);periodic=(true)",
        )
        @test ret == path

        a = load_archive(path)

        # Metadata provenance.
        @test a.meta isa RunMetadata
        @test a.meta.julia_version == string(VERSION)
        @test !isempty(a.meta.git_commit)
        @test a.meta.rng_seed == seed
        @test a.meta.normalization == "Omega_ci"
        @test a.meta.boundary_desc == "periodic"
        @test a.meta.filter_desc == "exp_filter(36,8)"
        @test a.meta.diagnostic_desc == "energy"
        @test a.meta.rank_layout == "mpi;ranks=2;dims=(2);coords=(1);periodic=(true)"

        # State round-trips exactly.
        @test a.state.D == 1
        @test a.state.ncell == g.n
        @test a.state.L == g.L
        @test a.state.step == st.step[]
        @test a.state.time == st.time[]
        @test a.state.q == ps.q && a.state.m == ps.m
        @test a.state.x == ps.x
        @test a.state.v == ps.v
        @test a.state.weight == ps.weight
        @test a.state.id == ps.id
        @test a.state.tag == ps.tag
        @test a.state.B == st.fields.B
        @test a.state.E == st.fields.E
    end
end

@testset "archive_run validates rng_seed range" begin
    _, _, st, ps = _build_run(; N = 16)
    mktemp() do path, io
        close(io)
        @test_throws ArgumentError archive_run(path, st, ps; rng_seed = big(typemax(Int)) + 1)
    end
end

@testset "archive corruption detected (checksum via load_run)" begin
    g, model, st, ps = _build_run(; N = 32)
    mktemp() do path, io
        close(io)
        archive_run(path, st, ps; rng_seed = 1)
        # load_archive must reject a tampered checksum (it delegates to load_run).
        cont = deserialize(path)
        bad = (
            schema = cont.schema,
            meta = cont.meta,
            state = cont.state,
            checksum = cont.checksum == repeat("0", 64) ? repeat("1", 64) : repeat("0", 64),
        )
        serialize(path, bad)
        @test_throws ErrorException load_archive(path)
    end
end
