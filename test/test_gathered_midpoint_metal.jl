# Included only after the optional-backend owner has established a functional
# device. The external audit also runs this file with an actual Metal device.
using HybridPlasmaPIC, Test
import Metal

@testset "Metal gathered midpoint storage preflight" begin
    @test !isdefined(Metal, :functional) || getproperty(Metal, :functional)()
    for D = 1:3, N in (0, 4)
        host = ParticleSet{D,Float32}(N)
        for d = 1:D
            fill!(host.x[d], Float32(d))
        end
        for c = 1:3
            fill!(host.v[c], Float32(2c))
        end
        ps = copy_particles_to_backend(Val(:metal), host)
        E = ntuple(_ -> Metal.MtlArray(zeros(Float32, N)), 3)
        B = ntuple(_ -> Metal.MtlArray(zeros(Float32, N)), 3)
        midpoint = ntuple(_ -> Metal.MtlArray(fill(99.0f0, N)), D)
        function rejected_unchanged(mid, fields = E)
            arrays = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, fields..., B..., mid...)
            before = map(Array, arrays)
            @test_throws ArgumentError push_gathered!(ps, fields, B, 0.5f0; xmid = mid)
            Metal.synchronize()
            @test isequal(map(Array, arrays), before)
        end
        for source in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, E..., B...), d = 1:D
            rejected_unchanged(Base.setindex(midpoint, source, d))
        end
        if D >= 2
            rejected_unchanged(Base.setindex(midpoint, midpoint[1], 2))
        end
        if N > 0
            backing = Metal.MtlArray(zeros(Float32, N + 1))
            field, output = view(backing, 1:N), view(backing, 2:N+1)
            @test field !== output && Base.mightalias(field, output)
            rejected_unchanged(Base.setindex(midpoint, output, 1), (field, E[2], E[3]))
        end
        shared = (E[1], E[1], E[1])
        @test push_gathered!(ps, shared, shared, 0.5f0; xmid = midpoint) === ps
        Metal.synchronize()
        for d = 1:D
            @test all(==(Float32(2d)), Array(ps.x[d]))
            @test all(==(Float32(3d / 2)), Array(midpoint[d]))
        end
        @test Array(ps.weight) == host.weight
        @test Array(ps.id) == host.id
        @test Array(ps.tag) == host.tag
        @test all(iszero, Array(E[1]))
    end
end
