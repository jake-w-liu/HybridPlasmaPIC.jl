using HybridPlasmaPIC, Test

# The midpoint suite owns the common initialized storage fixture.
function gfs_rejected_unchanged(ps, E, B, mid, dt, label, component)
    arrays =
        (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, E..., B..., (mid === nothing ? () : mid)...)
    before = map(copy, arrays)
    error = try
        push_gathered!(ps, E, B, dt; xmid = mid)
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test error isa ArgumentError && occursin("$label[$component]", sprint(showerror, error))
    @test all(isequal(array, saved) for (array, saved) in zip(arrays, before))
end

@testset "Gathered field inputs cannot overlap written particle storage" begin
    for T in (Float16, Float32, Float64), D = 1:3, N in (0, 1, 4), dt in (0, 0.5)
        ps, E, B, midpoint = gms_fixture(T, Val(D), N)
        for mid in (nothing, midpoint), source in (ps.x..., ps.v...), c = 1:3
            gfs_rejected_unchanged(ps, Base.setindex(E, source, c), B, mid, dt, "E", c)
            gfs_rejected_unchanged(ps, E, Base.setindex(B, source, c), mid, dt, "B", c)
            wrapped = view(source, :)
            @test wrapped !== source
            @test HybridPlasmaPIC._particle_arrays_mightalias(wrapped, source)
            gfs_rejected_unchanged(ps, Base.setindex(E, wrapped, c), B, mid, dt, "E", c)
            gfs_rejected_unchanged(ps, E, Base.setindex(B, wrapped, c), mid, dt, "B", c)
            if N > 0
                reversed = view(source, N:-1:1)
                @test reversed !== source && Base.mightalias(reversed, source)
                gfs_rejected_unchanged(ps, Base.setindex(E, reversed, c), B, mid, dt, "E", c)
                gfs_rejected_unchanged(ps, E, Base.setindex(B, reversed, c), mid, dt, "B", c)
            end
        end
    end
end

@testset "Reversed gathered input has an independent fixed-field result" begin
    for T in (Float32, Float64), D = 1:3
        ps = ParticleSet{D,T}(3)
        foreach(x -> fill!(x, zero(T)), ps.x)
        ps.v[1] .= T[1, 2, 3]
        E = (view(ps.v[1], 3:-1:1), zeros(T, 3), zeros(T, 3))
        B = ntuple(_ -> zeros(T, 3), 3)
        reference = deepcopy(ps)
        push_gathered!(reference, map(copy, E), B, T(0.5))
        @test reference.v[1] == T[2.5, 3, 3.5]
        gfs_rejected_unchanged(ps, E, B, nothing, T(0.5), "E", 1)
    end
end

@testset "Read-only particle weights and disjoint views remain valid inputs" begin
    for T in (Float16, Float32, Float64), D = 1:3, N in (0, 1, 4)
        for readonly in (:weight, :id, :tag)
            ps, _, _, mid = gms_fixture(T, Val(D), N)
            shared = getfield(ps, readonly)
            fields = (shared, shared, shared)
            reference, reference_mid = deepcopy(ps), map(copy, mid)
            frozen = map(copy, fields)
            push_gathered!(reference, frozen, frozen, T(0.5); xmid = reference_mid)
            @test push_gathered!(ps, fields, fields, T(0.5); xmid = mid) === ps
            @test isequal(
                (ps.x..., ps.v..., mid...),
                (reference.x..., reference.v..., reference_mid...),
            )
            @test all(isone, ps.weight)
            @test ps.id == UInt64.(1:N)
            @test all(iszero, ps.tag)
        end
        if N > 0
            host, _, _, mid = gms_fixture(T, Val(D), N)
            backing = zeros(T, 4N)
            velocities = ntuple(c -> view(backing, (c-1)*N+1:c*N), 3)
            for c = 1:3
                fill!(velocities[c], T(c / 16))
            end
            ps =
                ParticleSet{D,T}(host.x, velocities, host.weight, host.id, host.tag, host.q, host.m)
            input = view(backing, 3N+1:4N)
            fill!(input, T(1 / 32))
            @test all(v -> !Base.mightalias(input, v), ps.v)
            E, B = (input, input, input), ntuple(_ -> zeros(T, N), 3)
            reference, reference_mid = deepcopy(ps), map(copy, mid)
            push_gathered!(reference, map(copy, E), B, T(0.5); xmid = reference_mid)
            @test push_gathered!(ps, E, B, T(0.5); xmid = mid) === ps
            @test isequal(
                (ps.x..., ps.v..., mid...),
                (reference.x..., reference.v..., reference_mid...),
            )
            @test all(==(T(1 / 32)), input)
        end
    end
end
