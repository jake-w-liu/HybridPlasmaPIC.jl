@testset "Scaled compensated Lorentz factor resolves stored operands" begin
    for T in (Float16, Float32, Float64)
        rng = Xoshiro(0x6c6f72656e747a)
        admitted = 0
        for _ = 1:2000
            direction = randn(rng, 3)
            direction ./= norm(direction)
            k = rand(rng, 2:(precision(T)-2))
            speed = one(T) - ldexp(one(T), -k)
            span = T === Float16 ? 8 : T === Float32 ? 60 : 500
            c = ldexp(one(T), rand(rng, (-span):span))
            v = Tuple(T.(direction) .* speed .* c)
            β = v ./ c
            β2 = β[1] * β[1] + β[2] * β[2] + β[3] * β[3]
            gap = (Q(c)^2 - sum(x -> Q(x)^2, v)) / Q(c)^2
            β2 < one(T) && gap > 0 || continue
            factor, resolved = H._conditioned_lorentz_factor(v, c, β2)
            @test resolved
            reference = H._sqrt_from_exact(T, inv(gap))
            @test isfinite(factor) && isfinite(reference)
            @test abs(Q(factor) - Q(reference)) <= Q(eps(reference))
            admitted += 1
        end
        @test admitted >= 1900
    end
    # A rounded subluminal comparison may hide a negative exact gap. It must
    # still reach the exact rejection policy rather than an artificial factor.
    for (T, raw) in (
        (Float32, (-0.86362565, -0.49116498, 0.113612115)),
        (Float64, (-0.474852586251129, 0.8763033506698198, -0.0812862776577999)),
    )
        v = T.(raw)
        z = ntuple(_ -> zero(T), 3)
        @test sum(x -> Q(x)^2, v) > 1
        @test sum(abs2, v) < one(T)
        @test_throws ArgumentError H._boris_rel(v..., z..., z..., one(T), T(0.125), one(T))
    end
end
