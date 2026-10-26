using HybridPlasmaPIC, Test
import Metal

@testset "Metal gathered inputs reject written-particle aliases" begin
    @test Metal.functional()
    for D = 1:3, N in (0, 4)
        host = ParticleSet{D,Float32}(N)
        foreach(x -> fill!(x, 1.0f0), host.x)
        foreach(v -> fill!(v, 0.25f0), host.v)
        ps = copy_particles_to_backend(Val(:metal), host)
        E = ntuple(_ -> Metal.MtlArray(zeros(Float32, N)), 3)
        B = ntuple(_ -> Metal.MtlArray(zeros(Float32, N)), 3)
        midpoint = ntuple(_ -> Metal.MtlArray(zeros(Float32, N)), D)
        function rejected_unchanged(fields, magnetic, mid, dt)
            arrays = (
                ps.x...,
                ps.v...,
                ps.weight,
                ps.id,
                ps.tag,
                fields...,
                magnetic...,
                (mid === nothing ? () : mid)...,
            )
            before = map(Array, arrays)
            @test_throws ArgumentError push_gathered!(ps, fields, magnetic, dt; xmid = mid)
            Metal.synchronize()
            @test isequal(map(Array, arrays), before)
        end
        for mid in (nothing, midpoint), source in (ps.x..., ps.v...), c = 1:3, dt in (0.0f0, 0.5f0)
            rejected_unchanged(Base.setindex(E, source, c), B, mid, dt)
            rejected_unchanged(E, Base.setindex(B, source, c), mid, dt)
            if N > 0
                # Contiguous views remain native MtlArray inputs to the device method.
                wrapped = view(source, 1:N)
                @test wrapped isa Metal.MtlArray{Float32,1}
                @test wrapped !== source && Base.dataids(wrapped) == Base.dataids(source)
                @test HybridPlasmaPIC._particle_arrays_mightalias(wrapped, source)
                rejected_unchanged(Base.setindex(E, wrapped, c), B, mid, dt)
                rejected_unchanged(E, Base.setindex(B, wrapped, c), mid, dt)
            end
        end
        # The weights are read-only during the push, so sharing them remains valid.
        inputs = (ps.weight, ps.weight, ps.weight)
        reference = copy_particles_to_backend(Val(:metal), host)
        frozen = ntuple(_ -> Metal.MtlArray(ones(Float32, N)), 3)
        reference_midpoint = ntuple(_ -> Metal.MtlArray(zeros(Float32, N)), D)
        push_gathered!(reference, frozen, frozen, 0.5f0; xmid = reference_midpoint)
        @test push_gathered!(ps, inputs, inputs, 0.5f0; xmid = midpoint) === ps
        Metal.synchronize()
        @test isequal(
            map(Array, (ps.x..., ps.v..., midpoint...)),
            map(Array, (reference.x..., reference.v..., reference_midpoint...)),
        )
        @test Array(ps.weight) == host.weight
    end
end
