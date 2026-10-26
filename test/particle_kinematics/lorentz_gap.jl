# Included in ParticleKinematicsRangeTests so all references are the independent
# Cayley/exact-gap oracles above; no production arithmetic is used by the oracle.
@testset "Relativistic kinematics resolve an ill-conditioned Lorentz gap" begin
    for T in (Float16, Float32, Float64)
        z = (zero(T), zero(T), zero(T))
        k = T === Float16 ? 2 : T === Float32 ? 2 : 3
        witness = (prevfloat(T(0.6), k), T(0.8), zero(T))
        # Include non-boundary high-speed cases: a repair limited to a few eps
        # around c misses the same conditioning mechanism farther from c.
        speeds = (T(0.75), T(0.9), T(0.99))
        states = (witness, ((s * T(0.6), s * T(0.8), zero(T)) for s in speeds)...)
        for raw in states,
            order in ((1, 2, 3), (2, 3, 1), (3, 1, 2)),
            sign in (-one(T), one(T)),
            e in (-8, 0, 8)

            c = ldexp(one(T), e)
            v = ntuple(i -> sign * ldexp(raw[order[i]], e), 3)
            ratios = v ./ c
            β2 = sum(abs2, ratios)
            # The existing accepted domain requires the rounded comparison as
            # well as a physically positive exact gap. Preserve that domain.
            β2 < one(T) && sum(x -> Q(x)^2, v) < Q(c)^2 || continue
            @test !H._lorentz_gap_in_range(β2)
            B = ntuple(i -> order[i] == 3 ? inv(sqrt(eps(T))) : zero(T), 3)
            qm, dt = one(T), T(2)
            reference, scale = boris_rel_reference(v, z, B, qm, dt, c)
            actual = H._boris_rel(v..., z..., B..., qm, dt, c)
            @test rel_ok(actual, reference, scale)
            for h in (T(0.125), T(2))
                reference, scale = rel_prime_reference(v, z, B, qm, h, c)
                vv = ntuple(i -> [v[i]], 3)
                H._prime_leapfrog!(
                    vv,
                    ntuple(i -> [z[i]], 3),
                    ntuple(i -> [B[i]], 3),
                    qm,
                    h,
                    1,
                    true,
                    c,
                )
                @test rel_ok(map(only, vv), reference, scale)
            end
        end
        # Admission-boundary neighbors and ordinary native identity.
        for β2 in (prevfloat(T(0.5)), T(0.5), nextfloat(T(0.5)))
            @test H._lorentz_gap_in_range(β2) == (β2 <= T(0.5))
        end
        for speed in (prevfloat(T(0.5)), T(0.5), nextfloat(T(0.5)))
            v = (speed, zero(T), zero(T))
            E, B, qm, dt, c = z, (zero(T), zero(T), one(T)), one(T), T(0.125), one(T)
            @test isequal(
                H._boris_rel(v..., E..., B..., qm, dt, c),
                old_boris_rel(v..., E..., B..., qm, dt, c),
            )
        end
    end
end

@testset "Relativistic batch owners cannot bypass Lorentz-gap admission" begin
    for T in (Float32, Float64)
        k = T === Float32 ? 2 : 3
        fast = (prevfloat(T(0.6), k), T(0.8), zero(T))
        slow = (T(0.125), -T(0.0625), T(0.03125))
        z = (zero(T), zero(T), zero(T))
        B = (zero(T), zero(T), inv(sqrt(eps(T))))
        c, dt = one(T), T(2)
        for D in (1, 3), N in (1, 17), index in unique((1, N))
            ps = ParticleSet{D,T}(N)
            foreach(x -> fill!(x, T(0.25)), ps.x)
            for d = 1:3
                fill!(ps.v[d], slow[d])
                ps.v[d][index] = fast[d]
            end
            Ep, Bp = ntuple(i -> fill(z[i], N), 3), ntuple(i -> fill(B[i], N), 3)
            scales = H._boris_batch_scales(ps.v, Ep, Bp, dt / 2, 1:N)
            @test !H._boris_rel_batch_certificate(dt / 2, one(T), dt, c, scales...)
            g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
            es =
                D == 1 ? EMPIC1D(g, N; c, relativistic = true) : EMPIC(g, N; c, relativistic = true)
            before = map(copy, ps.v)
            if D == 1
                H._push_velocities!(es, ps, Ep[1], Ep[2], Bp[3], dt)
            else
                H._push_velocities!(es, ps, Ep, Bp, dt)
            end
            reference, scale = boris_rel_reference(fast, z, B, one(T), dt, c)
            @test rel_ok(ntuple(i -> ps.v[i][index], 3), reference, scale)
            for p in setdiff(1:N, (index,))
                @test isequal(
                    ntuple(i -> ps.v[i][p], 3),
                    old_boris_rel(slow..., z..., B..., one(T), dt, c),
                )
            end
            for planar in (false, true)
                vv = map(copy, before)
                if planar
                    H._prime_empic1d!(vv, Ep[1], Ep[2], Bp[3], one(T), dt, N, true, c)
                else
                    H._prime_leapfrog!(vv, Ep, Bp, one(T), dt, N, true, c)
                end
                reference, scale = rel_prime_reference(fast, z, B, one(T), dt, c)
                @test rel_ok(ntuple(i -> vv[i][index], 3), reference, scale)
            end
        end
    end
end

@testset "Public EMPIC steps retain exact-gap particle dynamics" begin
    # Zero-weight tracer particles feel the uniform seeded magnetic field but
    # source no current. This makes the full init + prime + step sequence
    # independently reproducible without replacing any production operator.
    for T in (Float32, Float64), D in (1, 3), speed in (T(0.9), T(0.99))
        dt, c, qm = T(0.125), one(T), -one(T)
        v = (T(0.6) * speed, T(0.8) * speed, zero(T))
        z = (zero(T), zero(T), zero(T))
        B = (zero(T), zero(T), inv(sqrt(eps(T))))
        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(4), D))
        es =
            D == 1 ? EMPIC1D(g, 1; n0 = zero(T), c, relativistic = true) :
            EMPIC(g, 1; n0 = zero(T), c, relativistic = true)
        ps = ParticleSet{D,T}(1; q = qm)
        foreach(x -> fill!(x, T(0.25)), ps.x)
        fill!(ps.weight, zero(T))
        for i = 1:3
            ps.v[i][1] = v[i]
        end
        if D == 1
            fill!(es.Bz, B[3])
        else
            for i = 1:3
                fill!(es.B[i], B[i])
            end
        end
        init_empic!(es, ps)
        reference, _ = rel_prime_reference(v, z, B, qm, dt / 2, c)
        primed = T.(reference)
        @test sum(abs2, primed) < one(T)
        reference, scale = boris_rel_reference(primed, z, B, qm, dt, c)
        step_empic!(es, ps, dt)
        @test rel_ok(map(only, ps.v), reference, scale)
        @test es.time[] == dt && es.step[] == 1
        @test ps.weight == T[0]
        @test all(a -> all(isfinite, a), (ps.x..., ps.v...))
    end
end

include("lorentz_factor.jl")
