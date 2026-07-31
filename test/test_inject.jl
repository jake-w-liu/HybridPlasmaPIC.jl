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

function negative_flux_oracle(a, σ)
    b = BigFloat(-a / σ)
    xmax = max(BigFloat(2), BigFloat(12) / b)
    N = 400_000
    dx = xmax / N
    den = zero(BigFloat)
    num = zero(BigFloat)
    for i = 1:N
        x = (BigFloat(i) - big"0.5") * dx
        wt = x * exp(-b * x - x^2 / 2)
        den += wt
        num += x * wt
    end
    shape = den * dx
    flux = BigFloat(σ) / sqrt(BigFloat(2) * BigFloat(π)) * exp(-(b^2) / 2) * shape
    mean_speed = BigFloat(σ) * num / den
    return Float64(flux), Float64(mean_speed)
end

function negative_shape_asymptotic_oracle(b::BigFloat)
    invb2 = inv(b)^2
    term = invb2
    total = term
    for n = 1:10_000
        nextterm = term * BigFloat(2n + 1) * invb2
        nextterm >= term && break
        total += isodd(n) ? -nextterm : nextterm
        nextterm <= eps(BigFloat) * abs(total) && break
        term = nextterm
    end
    return total
end

function negative_tail_flux_oracle(σ, b)
    return setprecision(BigFloat, 1024) do
        σbig = BigFloat(σ)
        bbig = BigFloat(b)
        σbig / sqrt(BigFloat(2) * BigFloat(π)) *
        exp(-(bbig^2) / 2) *
        negative_shape_asymptotic_oracle(bbig)
    end
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

@testset "flux sampler resolves the strong outward-drift tail" begin
    rng = MersenneTwister(11)
    for (a, σ) in ((-0.5, 0.1), (-1.5, 0.1))
        flux_ref, mean_ref = negative_flux_oracle(a, σ)
        @test isapprox(flux_per_density(a, σ), flux_ref; rtol = 2e-5)
        samples = [flux_speed(rng, a, σ) for _ = 1:100_000]
        @test all(>(0), samples)
        @test isapprox(mean(samples), mean_ref; rtol = 0.01)
    end
    @test flux_per_density(-1.5, 0.1) > 0
    # the reachable, positive-flux regime is unchanged and stays strictly positive
    for _ = 1:1000
        @test flux_speed(rng, 2.0, 1.0) > 0
    end
end

@testset "negative-tail flux survives Gaussian underflow" begin
    for (T, σ, tail_points) in (
        (Float64, 1.0e300, (30.0, 38.0, 40.0, 50.0, 55.0)),
        (Float32, 1.0f30, (10.0f0, 14.0f0, 15.0f0, 18.0f0, 21.0f0)),
    )
        tiny = nextfloat(zero(T))
        for b in tail_points
            a = -b * σ
            expected = T(negative_tail_flux_oracle(σ, -a / σ))
            actual = flux_per_density(a, σ)
            if iszero(expected)
                @test iszero(actual)
            elseif expected < floatmin(T)
                @test abs(actual - expected) <= 2tiny
            else
                @test isapprox(actual, expected; rtol = 16eps(T))
            end
        end

        gaussian_normal_limit = sqrt(-T(2) * log(floatmin(T)))
        for b in (prevfloat(gaussian_normal_limit), nextfloat(gaussian_normal_limit))
            a = -b * σ
            expected = T(negative_tail_flux_oracle(σ, -a / σ))
            actual = flux_per_density(a, σ)
            # The tail is conditioned by exp(-b²/2): one rounding error in
            # b² is amplified by O(b²), even when the returned flux is normal.
            @test isapprox(actual, expected; rtol = 4b^2 * eps(T), atol = 2tiny)
        end

        @test flux_per_density(-floatmax(T), one(T)) == zero(T)
        @test flux_per_density(floatmax(T), one(T)) == floatmax(T)
        @test flux_per_density(-one(T), tiny) == zero(T)
        @test flux_per_density(one(T), tiny) == one(T)
        positive_scale = T === Float32 ? T(1.0e30) : T(1.0e300)
        @test flux_per_density(positive_scale, positive_scale) / positive_scale ≈
              flux_per_density(one(T), one(T)) rtol = 4eps(T)
    end

    @test flux_per_density(-4.0e301, 1.0e300) == 9.128344722912974e-52
end

@testset "flux arithmetic is scale-stable near the Float64 range" begin
    scale = 1.0e200
    unit_flux = flux_per_density(1.0, 1.0)
    scaled_flux = flux_per_density(scale, scale)
    @test isfinite(scaled_flux)
    @test scaled_flux / scale ≈ unit_flux rtol = 2eps(Float64)

    unit_speed = flux_speed(MersenneTwister(19), 1.0, 1.0)
    scaled_speed = flux_speed(MersenneTwister(19), scale, scale)
    @test isfinite(scaled_speed)
    @test scaled_speed / scale ≈ unit_speed rtol = 2eps(Float64)

    # n0*flux overflows before a subnormal dt brings the complete metered
    # increment back to O(1).
    ps = ParticleSet{1,Float64}(0)
    acc = Ref(0.5)
    nextid = Ref(UInt64(1))
    ninj = inject_face_1d!(
        ps,
        MersenneTwister(20),
        0.0,
        +1,
        1.0e200,
        1.0e120,
        0.0,
        (0.0, 0.0),
        0.0,
        1.0e-320,
        1.0,
        acc,
        nextid,
    )
    @test ninj == 1
    @test nparticles(ps) == 1
    @test isfinite(acc[]) && 0.0 <= acc[] < 1.0
    @test nextid[] == UInt64(2)
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
        (;
            face_x = NaN,
            inward = +1,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = 0.3,
            dt = 1.0,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = 0,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = 0.3,
            dt = 1.0,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = 2,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = 0.3,
            dt = 1.0,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = +1,
            n0 = -1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = 0.3,
            dt = 1.0,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = +1,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (NaN, 0.0),
            σt = 0.3,
            dt = 1.0,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = +1,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = NaN,
            dt = 1.0,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = +1,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = 0.3,
            dt = NaN,
            w = 1.0,
        ),
        (;
            face_x = 0.0,
            inward = +1,
            n0 = 1.0,
            a = 2.0,
            σ = 0.5,
            ut = (0.0, 0.0),
            σt = 0.3,
            dt = 1.0,
            w = 0.0,
        ),
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

@testset "inject_face_1d! rejects invalid metering state before mutation" begin
    rng = MersenneTwister(13)
    for acc0 in (-1.0, NaN, Inf)
        ps = ParticleSet{1,Float64}(0)
        acc = Ref(acc0)
        nid = Ref(UInt64(1))
        @test_throws ArgumentError inject_face_1d!(
            ps,
            rng,
            0.0,
            +1,
            0.0,
            1.0,
            0.0,
            (0.0, 0.0),
            0.0,
            1.0,
            1.0,
            acc,
            nid,
        )
        @test nparticles(ps) == 0
        @test isequal(acc[], acc0)
        @test nid[] == UInt64(1)
    end

    for nid0 in (UInt64(0), typemax(UInt64))
        ps = ParticleSet{1,Float64}(0)
        acc = Ref(0.0)
        nid = Ref(nid0)
        @test_throws ArgumentError inject_face_1d!(
            ps,
            rng,
            0.0,
            +1,
            1.0,
            1.0,
            0.0,
            (0.0, 0.0),
            0.0,
            1.0,
            1.0,
            acc,
            nid,
        )
        @test nparticles(ps) == 0
        @test acc[] == 0.0
        @test nid[] == nid0
    end

    ps = ParticleSet{1,Float64}(0)
    acc = Ref(0.0)
    nid = Ref(UInt64(1))
    @test_throws ArgumentError inject_face_1d!(
        ps,
        rng,
        0.0,
        +1,
        floatmax(Float64),
        1.0,
        0.0,
        (0.0, 0.0),
        0.0,
        2.0,
        1.0,
        acc,
        nid,
    )
    @test nparticles(ps) == 0
    @test acc[] == 0.0
    @test nid[] == UInt64(1)
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
