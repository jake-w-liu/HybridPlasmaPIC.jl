using HybridPlasmaPIC, Test

function gms_fixture(::Type{T}, ::Val{D}, N) where {T,D}
    ps = ParticleSet{D,T}(N)
    for d = 1:D
        fill!(ps.x[d], T(d))
    end
    for c = 1:3
        fill!(ps.v[c], T(2c))
    end
    E = ntuple(_ -> zeros(T, N), 3)
    B = ntuple(_ -> zeros(T, N), 3)
    midpoint = ntuple(_ -> fill(T(99), N), D)
    return ps, E, B, midpoint
end

function gms_rejected_unchanged(ps, E, B, midpoint, dt)
    arrays = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, E..., B..., midpoint...)
    before = map(copy, arrays)
    @test_throws ArgumentError push_gathered!(ps, E, B, dt; xmid = midpoint)
    @test all(isequal(a, saved) for (a, saved) in zip(arrays, before))
end

@testset "Gathered midpoint aliases are rejected before publication" begin
    for T in (Float16, Float32, Float64), D = 1:3, N in (0, 1, 4), dt in (0, 0.5)
        ps, E, B, midpoint = gms_fixture(T, Val(D), N)
        # Include weights and integer metadata: a push must not overwrite them.
        for source in (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, E..., B...)
            for d = 1:D
                bad = Base.setindex(midpoint, source, d)
                gms_rejected_unchanged(ps, E, B, bad, dt)
            end
        end
        for d = 1:D
            bad = Base.setindex(midpoint, view(ps.x[d], :), d)
            gms_rejected_unchanged(ps, E, B, bad, dt)
        end
        if D >= 2
            bad = Base.setindex(midpoint, midpoint[1], 2)
            gms_rejected_unchanged(ps, E, B, bad, dt)
        end
        if N == 4
            backing = zeros(T, N + 1)
            field = view(backing, 1:N)
            output = view(backing, 2:N+1)
            @test field !== output && Base.mightalias(field, output)
            gms_rejected_unchanged(
                ps,
                (field, E[2], E[3]),
                B,
                Base.setindex(midpoint, output, 1),
                dt,
            )
            if D >= 2
                bad = Base.setindex(Base.setindex(midpoint, field, 1), output, 2)
                gms_rejected_unchanged(ps, E, B, bad, dt)
            end
        end
    end
end

@testset "Disjoint midpoint outputs preserve the gathered-push contract" begin
    for T in (Float16, Float32, Float64), D = 1:3, N in (0, 1, 4)
        ps, E, B, midpoint = gms_fixture(T, Val(D), N)
        shared = zeros(T, N)
        readonly = (shared, shared, shared)
        @test push_gathered!(ps, readonly, readonly, T(0.5); xmid = midpoint) === ps
        for d = 1:D
            @test all(==(T(2d)), ps.x[d])
            @test all(==(T(3d / 2)), midpoint[d])
        end
        @test all(iszero, shared)
        @test all(isone, ps.weight)
        @test ps.id == UInt64.(1:N)
        @test all(iszero, ps.tag)
        if N > 0
            # One parent allocation, provably disjoint element ranges. The
            # field and midpoint components are distinct views, not copies.
            ps, E, B, _ = gms_fixture(T, Val(D), N)
            storage = zeros(T, (D + 1) * N)
            field = view(storage, 1:N)
            midpoint = ntuple(d -> view(storage, d*N+1:(d+1)*N), D)
            @test all(a -> !Base.mightalias(a, field), midpoint)
            @test push_gathered!(ps, (field, field, field), B, T(0.5); xmid = midpoint) === ps
            for d = 1:D
                @test all(==(T(2d)), ps.x[d])
                @test all(==(T(3d / 2)), midpoint[d])
            end
            @test all(iszero, field)
        end
    end
end

@testset "Gathered midpoint components must be vectors" begin
    ps, E, B, _ = gms_fixture(Float64, Val(1), 1)
    for bad in ((99.0,), (fill(99.0, 1, 1),))
        before = map(copy, (ps.x..., ps.v...))
        @test_throws ArgumentError push_gathered!(ps, E, B, 0.5; xmid = bad)
        @test isequal((ps.x..., ps.v...), before)
    end
end
