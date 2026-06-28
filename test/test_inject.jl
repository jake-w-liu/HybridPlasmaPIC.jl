# §11.4 flux-weighted injection: the sampler reproduces the analytic inward-flux
# moments, and inject_face_1d! delivers the correct number flux into a box.

using HybridPlasmaPIC, Test, Random, Statistics

# independent oracle: ⟨v_n⟩ of p(s) ∝ s·exp(−(s−a)²/2σ²) by quadrature
function flux_mean_oracle(a, σ)
    N = 400_000
    smax = a + 14σ
    ds = smax / N
    num = 0.0
    den = 0.0
    for i = 1:N
        s = (i - 0.5) * ds
        wt = s * exp(-(s - a)^2 / (2σ^2))
        num += s * wt * ds
        den += wt * ds
    end
    return num / den
end

@testset "flux sampler reproduces inward-flux moments (LOAD)" begin
    rng = MersenneTwister(1)
    M = 400_000   # SE ≪ rtol for this M
    for (a, σ) in ((0.0, 1.0), (1.5, 1.0), (3.0, 0.7))
        ss = [flux_speed(rng, a, σ) for _ = 1:M]
        @test isapprox(mean(ss), flux_mean_oracle(a, σ); rtol = 0.01)
        if a == 0
            @test isapprox(mean(abs2, ss), 2σ^2; rtol = 0.01)   # Rayleigh ⟨s²⟩ = 2σ²
        end
    end
end

@testset "flux sampler handles the cold-beam limit" begin
    rng = MersenneTwister(2)
    @test flux_per_density(2.0, 0.0) == 2.0
    @test flux_per_density(0.0, 0.0) == 0.0
    @test flux_per_density(-2.0, 0.0) == 0.0
    @test flux_speed(rng, 2.0, 0.0) == 2.0
    @test flux_speed(rng, 0.0, 0.0) == 0.0

    T = Float64
    ps = ParticleSet{1,T}(0)
    acc = Ref(0.0)
    nid = Ref(UInt64(1))
    ninj = inject_face_1d!(ps, rng, 0.0, +1, 1.0, 2.0, 0.0, (0.0, 0.0), 0.3, 1.0, 0.5, acc, nid)
    @test ninj == 4
    @test all(==(2.0), ps.v[1])
    @test all(>=(0.0), ps.x[1])
    @test all(<=(2.0), ps.x[1])   # swept slab face_x .. face_x + a*dt
end

@testset "inject_face_1d! validates inputs before mutation" begin
    rng = MersenneTwister(3)

    function snapshot()
        ps = ParticleSet{1,Float64}(0)
        acc = Ref(10.0)
        nid = Ref(UInt64(1))
        return ps, acc, nid
    end

    for kwargs in (
        (; face_x = NaN, inward = +1, n0 = 1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = 0.3, dt = 1.0, w = 1.0),
        (; face_x = 0.0, inward = 0, n0 = 1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = 0.3, dt = 1.0, w = 1.0),
        (; face_x = 0.0, inward = 2, n0 = 1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = 0.3, dt = 1.0, w = 1.0),
        (; face_x = 0.0, inward = +1, n0 = -1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = 0.3, dt = 1.0, w = 1.0),
        (; face_x = 0.0, inward = +1, n0 = 1.0, a = 2.0, σ = 0.5, ut = (NaN, 0.0), σt = 0.3, dt = 1.0, w = 1.0),
        (; face_x = 0.0, inward = +1, n0 = 1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = NaN, dt = 1.0, w = 1.0),
        (; face_x = 0.0, inward = +1, n0 = 1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = 0.3, dt = NaN, w = 1.0),
        (; face_x = 0.0, inward = +1, n0 = 1.0, a = 2.0, σ = 0.5, ut = (0.0, 0.0), σt = 0.3, dt = 1.0, w = 0.0),
    )
        ps, acc, nid = snapshot()
        @test_throws ArgumentError inject_face_1d!(
            ps,
            rng,
            kwargs.face_x,
            kwargs.inward,
            kwargs.n0,
            kwargs.a,
            kwargs.σ,
            kwargs.ut,
            kwargs.σt,
            kwargs.dt,
            kwargs.w,
            acc,
            nid,
        )
        @test nparticles(ps) == 0
        @test isempty(ps.x[1]) && isempty(ps.v[1]) && isempty(ps.v[2]) && isempty(ps.v[3])
        @test isempty(ps.weight) && isempty(ps.id) && isempty(ps.tag)
        @test acc[] == 10.0
        @test nid[] == UInt64(1)
    end
end

@testset "inject_face_1d! delivers the target number flux" begin
    T = Float64
    n0 = 1.0
    a = 2.0
    σ = 0.5
    w = 0.01
    dt = 0.02
    ps = ParticleSet{1,T}(0)                     # empty box, inject at x=0 toward +x
    acc = Ref(0.0)
    nid = Ref(UInt64(1))
    rng = MersenneTwister(7)
    nsteps = 2000
    for _ = 1:nsteps
        inject_face_1d!(ps, rng, 0.0, +1, n0, a, σ, (0.0, 0.0), 0.3, dt, w, acc, nid)
    end
    Ttot = nsteps * dt
    flux = n0 * flux_per_density(a, σ)           # particles per unit area per time
    @test isapprox(sum(ps.weight), flux * Ttot; rtol = 0.02)        # number flux
    @test isapprox(mean(ps.v[1]), flux_mean_oracle(a, σ); rtol = 0.02)  # flux-weighted ⟨v_n⟩
    @test all(>(0), ps.v[1])                     # all injected moving inward
    @test allunique(ps.id)                       # unique ids
end
