module LegacyCurrentModeTests
using HybridPlasmaPIC, Test, Random
const H = HybridPlasmaPIC
const Q = Rational{BigInt}

original(source, k, dt) = -(source / dt) / complex(zero(k), k)
function exact(source::Complex{T}, k, dt) where {T}
    denominator = Q(dt) * Q(k)
    return complex(T(-Q(imag(source)) / denominator), T(Q(real(source)) / denominator))
end
cost(source, k, dt) = @allocated H._spectral_longitudinal_current(source, k, dt)

@testset "Legacy current retains complete modal quotients" begin
    for T in (Float32, Float64)
        shift = T === Float32 ? 100 : 900
        large, small = ldexp(one(T), shift), ldexp(one(T), -shift)
        qmin = nextfloat(zero(T))
        cases = (
            (complex(one(T), T(0.75)), large, floatmin(T) / 128),
            (complex(small, -small / 2), small, large),
            (complex(floatmax(T) / 8, floatmax(T) / 16), large, large),
            (complex(small, small / 2), small, small),
            (complex(floatmax(T), one(T)), small, floatmin(T)),
            (complex(qmin, zero(T)), one(T), T(2)),
            (complex(qmin, zero(T)), one(T), prevfloat(T(2))),
            (complex(qmin, zero(T)), one(T), nextfloat(T(2))),
            (zero(Complex{T}), large, floatmin(T) / 128),
        )
        for (source, k, dt) in cases, sign_source in (-1, 1), sign_k in (-1, 1)
            z, wave = sign_source * source, sign_k * k
            @test H._spectral_longitudinal_current(z, wave, dt) == exact(z, wave, dt)
            rotated = complex(-imag(z), real(z))
            @test H._spectral_longitudinal_current(rotated, wave, dt) == exact(rotated, wave, dt)
        end
        # The first two rows reject the old independently narrowed quotients.
        for (source, k, dt) in cases[1:2]
            @test !isequal(original(source, k, dt), exact(source, k, dt))
        end
        results = Vector{Complex{T}}(undef, 64)
        Threads.@threads for index in eachindex(results)
            source, k, dt = cases[mod1(index, length(cases))]
            results[index] = H._spectral_longitudinal_current(source, k, dt)
        end
        for index in eachindex(results)
            @test results[index] == exact(cases[mod1(index, length(cases))]...)
        end
    end
end

@testset "Legacy current preserves ordinary and ineligible arithmetic" begin
    rng = MersenneTwister(820917)
    for T in (Float32, Float64)
        for _ = 1:256
            source = complex(T(2rand(rng) - 1), T(2rand(rng) - 1))
            k = T(0.5 + rand(rng)) * (rand(rng, Bool) ? one(T) : -one(T))
            dt = T(0.1 + rand(rng))
            @test isequal(H._spectral_longitudinal_current(source, k, dt), original(source, k, dt))
            cost(source, k, dt)
            @test cost(source, k, dt) == 0
        end
        for source in (complex(T(NaN), zero(T)), complex(T(Inf), one(T)), complex(one(T), T(-Inf))),
            k in (one(T), -one(T)),
            dt in (one(T), T(Inf), zero(T))

            @test isequal(H._spectral_longitudinal_current(source, k, dt), original(source, k, dt))
        end
        for k in (T(Inf), T(-Inf), T(NaN)), dt in (one(T), T(Inf), zero(T))
            source = complex(one(T), T(0.5))
            @test isequal(H._spectral_longitudinal_current(source, k, dt), original(source, k, dt))
        end
    end
end

@testset "Legacy current completes finite four-point Fourier modes" begin
    for T in (Float32, Float64),
        profile in (:ordinary, :overflow, :underflow),
        sign_source in (-1, 1),
        sine in (false, true)

        shift = T === Float32 ? 100 : 900
        wave =
            profile === :ordinary ? one(T) : ldexp(one(T), profile === :overflow ? shift : -shift)
        dt =
            profile === :ordinary ? T(0.125) :
            profile === :overflow ? floatmin(T) / 128 : ldexp(one(T), shift)
        amplitude = profile === :underflow ? ldexp(one(T), -shift) : one(T)
        g = FourierGrid((4,), (T(2pi) / wave,))
        es = EMPIC1D(g, 0)
        e = ParticleSet{1,T}(0; q = -one(T))
        shape = sine ? T[0, 1, 0, -1] : T[1, 0, -1, 0]
        es.rho_np1 .= sign_source * amplitude .* shape
        before = (copy(es.rho_n), copy(es.rho_np1), copy(es.Ex), copy(es.Ey), copy(es.Bz))
        scale = T(sign_source * Q(amplitude) / (Q(dt) * Q(g.kvec[1][2])))
        expected = scale .* (sine ? T[1, 0, -1, 0] : T[0, -1, 0, 1])
        @test isfinite(scale) && !iszero(scale)
        @test H._esirkepov_Jx!(es, e, nothing, dt) == expected
        @test isequal(before, (es.rho_n, es.rho_np1, es.Ex, es.Ey, es.Bz))
    end
end

@testset "Legacy current retains physical DC and excludes Nyquist sources" begin
    for T in (Float32, Float64), mobile in (false, true)
        g = FourierGrid((4,), (T(4),))
        es = EMPIC1D(g, 1; mobile)
        e = ParticleSet{1,T}(1; q = -one(T))
        e.weight[1] = T(2)
        e.v[1][1] = T(0.5)
        ions = mobile ? ParticleSet{1,T}(1; q = one(T), m = es.mi) : nothing
        if mobile
            ions.weight[1] = one(T)
            ions.v[1][1] = T(2)
        end
        expected = mobile ? T(0.25) : T(-0.25)
        for source in (ones(T, 4), T[1, -1, 1, -1], zeros(T, 4))
            es.rho_np1 .= source
            @test all(==(expected), H._esirkepov_Jx!(es, e, ions, T(0.125)))
        end
    end
end
end
