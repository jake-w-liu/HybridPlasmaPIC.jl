using HybridPlasmaPIC, Test

@testset "Radiation input conversion validates before mutation" begin
    for T in (Float16, Float32, Float64), D = 1:3
        ps = ParticleSet{D,T}(2)
        foreach(x -> fill!(x, T(0.25)), ps.x)
        ps.v[1] .= T(0.5)
        before = deepcopy(ps)
        zero_field = (zero(T), zero(T), zero(T))
        wide = BigFloat(floatmax(T)) * 2
        bad_values = (NaN, Inf, -Inf, wide)
        for bad in bad_values, component = 1:3, field in (:E, :B), noop in (:none, :K, :dt)
            values = ntuple(i -> i == component ? bad : zero(bad), 3)
            E, B = field === :E ? (values, zero_field) : (zero_field, values)
            dt = noop === :dt ? zero(T) : one(T)
            K = noop === :K ? zero(T) : one(T)
            @test_throws ArgumentError apply_radiation_reaction!(ps, E, B, dt; K)
            @test ps.v == before.v
        end
        # Check source signs before narrowing: a negative subnormal source may
        # otherwise become signed zero in the destination precision.
        negative_tiny = -BigFloat(nextfloat(zero(T))) / 2
        for bad in (NaN, Inf, -Inf, wide, -one(T), negative_tiny)
            @test_throws ArgumentError apply_radiation_reaction!(
                ps,
                zero_field,
                zero_field,
                one(T);
                K = bad,
            )
            @test ps.v == before.v
            @test_throws ArgumentError apply_radiation_reaction!(
                ps,
                zero_field,
                zero_field,
                bad;
                K = one(T),
            )
            @test ps.v == before.v
        end
        for (K, dt) in ((zero(T), one(T)), (one(T), zero(T)))
            @test apply_radiation_reaction!(ps, zero_field, zero_field, dt; K) === ps
            @test ps.v == before.v
            # The established no-op contract does not inspect particle velocities.
            ps.v[1][end] = T(NaN)
            @test apply_radiation_reaction!(ps, zero_field, zero_field, dt; K) === ps
            @test isnan(ps.v[1][end])
            ps.v[1][end] = before.v[1][end]
        end
        # Wider finite inputs use the same validated conversion as the pusher.
        ordinary = deepcopy(ps)
        wide_inputs = deepcopy(ps)
        E = (T(0.3), T(0.2), T(-0.1))
        B = (T(0.1), T(-0.2), T(0.4))
        apply_radiation_reaction!(ordinary, E, B, T(0.125); K = T(0.25))
        apply_radiation_reaction!(
            wide_inputs,
            BigFloat.(E),
            BigFloat.(B),
            BigFloat(T(0.125));
            K = BigFloat(T(0.25)),
        )
        @test wide_inputs.v == ordinary.v
        @test all(isfinite, Iterators.flatten(wide_inputs.v))
        @test ps.x == before.x
        @test ps.weight == before.weight
        @test ps.id == before.id
        @test ps.tag == before.tag
    end
end
