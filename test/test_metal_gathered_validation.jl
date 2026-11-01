module MetalGatheredValidationTests
using HybridPlasmaPIC, Test

function population(D, count)
    ps = ParticleSet{D,Float32}(count; q = 1.5f0, m = 2.0f0)
    for d = 1:D, p = 1:count
        ps.x[d][p] = Float32(mod(p, 11) + d) / 8
    end
    for c = 1:3, p = 1:count
        ps.v[c][p] = Float32(mod(p + 3c, 7) - 3) / 16
    end
    return ps
end
fields(count) = (
    ntuple(c -> [Float32(mod(p + c, 7) - 3) / 32 for p = 1:count], 3),
    ntuple(c -> [Float32(mod(p + 2c, 5) - 2) / 64 for p = 1:count], 3),
)
function observed_arrays(ps, E, B, mid)
    return (ps.x..., ps.v..., ps.weight, ps.id, ps.tag, E..., B..., (mid === nothing ? () : mid)...)
end
function rejects_unchanged(ps, E, B, mid, dt, error = ArgumentError)
    arrays = observed_arrays(ps, E, B, mid)
    before = map(Array, arrays)
    @test_throws error push_gathered!(ps, E, B, dt; xmid = mid)
    @test isequal(map(Array, arrays), before)
    return nothing
end

@testset "Gathered finite-field admission precedes CPU publication" begin
    for component = 1:6, bad in (NaN32, Inf32, -Inf32), index in (1, 257)
        ps = population(2, 257)
        E, B = fields(257)
        (E..., B...)[component][index] = bad
        mid = ntuple(_ -> fill(-1.0f0, 257), 2)
        rejects_unchanged(ps, E, B, mid, 0.125f0)
    end
    for count in (0, 1), dt in (NaN32, Inf32, 1.0e40)
        ps = population(1, count)
        E, B = fields(count)
        rejects_unchanged(ps, E, B, nothing, dt)
    end
end

if Sys.isapple() && Base.find_package("Metal") !== nothing
    import Metal
    if Metal.functional()
        prepare_gpu_backend!(Val(:metal))
        disallow_scalar_indexing!(Val(:metal))

        # Freeze the current parent broadcast schedule, using its canonical
        # scalar owners rather than a second implementation of Boris arithmetic.
        function parent_broadcast!(ps, E, B, dt; xmid = nothing)
            vx, vy, vz = ps.v
            Ex, Ey, Ez = E
            Bx, By, Bz = B
            qm = ps.q / ps.m
            kicked = HybridPlasmaPIC._boris_kick_native.(vx, vy, vz, Ex, Ey, Ez, Bx, By, Bz, qm, dt)
            vx .= getindex.(kicked, 1)
            vy .= getindex.(kicked, 2)
            vz .= getindex.(kicked, 3)
            for d = 1:length(ps.x)
                if xmid !== nothing
                    xmid[d] .=
                        HybridPlasmaPIC._particle_coordinate_native.(ps.x[d], ps.v[d], dt / 2.0f0)
                end
                ps.x[d] .= HybridPlasmaPIC._particle_coordinate_native.(ps.x[d], ps.v[d], dt)
            end
            return ps
        end
        function device_population(D, count, alias = :none)
            ps = copy_particles_to_backend(Val(:metal), population(D, count))
            if alias == :velocity
                ps.v = (ps.v[1], ps.v[1], ps.v[3])
            elseif alias == :position
                ps.x = Base.setindex(ps.x, ps.v[1], 1)
            elseif alias == :overlap
                backing = Metal.MtlArray([Float32(mod(p, 7) - 3) / 16 for p = 1:(count+1)])
                ps.v = (view(backing, 1:count), view(backing, 2:(count+1)), ps.v[3])
            elseif alias == :weight
                ps.weight = ps.v[1]
            end
            return ps
        end
        bits(arrays) = map(a -> reinterpret(UInt32, Array(a)), arrays)

        @testset "Metal gathered fields reject before every public publication" begin
            for D = 1:3, count in (1, 257), component = 1:6, bad in (NaN32, Inf32, -Inf32)
                ps = device_population(D, count)
                E, B = fields(count)
                (E..., B...)[component][end] = bad
                Eg, Bg = map(Metal.MtlArray, E), map(Metal.MtlArray, B)
                mid = ntuple(_ -> Metal.MtlArray(fill(-1.0f0, count)), D)
                rejects_unchanged(ps, Eg, Bg, mid, 0.125f0)
            end
            # Exercise the ordered overlapping-output route without a midpoint.
            ps = device_population(2, 257, :overlap)
            E, B = fields(257)
            E[2][end] = NaN32
            rejects_unchanged(ps, map(Metal.MtlArray, E), map(Metal.MtlArray, B), nothing, 0.0f0)
            E[2][end] = 0.0f0
            @test push_gathered!(ps, map(Metal.MtlArray, E), map(Metal.MtlArray, B), 0.125f0) === ps
            Metal.synchronize()
            @test all(a -> all(isfinite, Array(a)), (ps.x..., ps.v...))
        end

        @testset "Metal gathered storage admission includes empty populations" begin
            for count in (0, 2)
                ps = device_population(2, count)
                E, B = fields(count)
                Eg, Bg = map(Metal.MtlArray, E), map(Metal.MtlArray, B)
                mid = ntuple(_ -> Metal.MtlArray(fill(-1.0f0, count)), 2)
                for c = 1:3
                    rejects_unchanged(ps, Base.setindex(Eg, E[c], c), Bg, mid, 0.125f0)
                    rejects_unchanged(ps, Eg, Base.setindex(Bg, B[c], c), mid, 0.125f0)
                    wrongtype = Metal.MtlArray(zeros(Int32, count))
                    rejects_unchanged(ps, Base.setindex(Eg, wrongtype, c), Bg, mid, 0.125f0)
                end
                rejects_unchanged(ps, Eg, Bg, ntuple(_ -> zeros(Float32, count), 2), 0.125f0)
                short = Metal.MtlArray(zeros(Float32, count + 1))
                rejects_unchanged(
                    ps,
                    Base.setindex(Eg, short, 1),
                    Bg,
                    mid,
                    0.125f0,
                    DimensionMismatch,
                )
            end
        end

        @testset "Metal gathered scheduling retains native parent bits" begin
            for D = 1:3,
                count in (0, 1, 257),
                midpoint in (false, true),
                dt in (-0.125f0, 0.0f0, 0.125f0)

                ps, reference = device_population(D, count), device_population(D, count)
                E, B = map(group -> map(Metal.MtlArray, group), fields(count))
                before_fields = map(Array, (E..., B...))
                mid = midpoint ? ntuple(_ -> Metal.MtlArray(fill(-1.0f0, count)), D) : nothing
                rmid = midpoint ? ntuple(_ -> Metal.MtlArray(fill(-1.0f0, count)), D) : nothing
                parent_broadcast!(reference, E, B, dt; xmid = rmid)
                @test push_gathered!(ps, E, B, dt; xmid = mid) === ps
                @test bits((ps.x..., ps.v..., (mid === nothing ? () : mid)...)) ==
                      bits((reference.x..., reference.v..., (rmid === nothing ? () : rmid)...))
                @test map(Array, (E..., B...)) == before_fields
            end
            for alias in (:velocity, :position, :overlap, :weight)
                ps, reference = device_population(2, 257, alias), device_population(2, 257, alias)
                E, B = map(group -> map(Metal.MtlArray, group), fields(257))
                mid = ntuple(_ -> Metal.MtlArray(fill(-1.0f0, 257)), 2)
                rmid = ntuple(_ -> Metal.MtlArray(fill(-1.0f0, 257)), 2)
                parent_broadcast!(reference, E, B, 0.125f0; xmid = rmid)
                push_gathered!(ps, E, B, 0.125f0; xmid = mid)
                @test bits((ps.x..., ps.v..., ps.weight, mid...)) ==
                      bits((reference.x..., reference.v..., reference.weight, rmid...))
            end
        end

        # A counter detects stores even when dt=0 leaves the stored value equal.
        struct CountedDestination{A,C}
            values::A
            stores::C
        end
        Metal.Adapt.adapt_structure(to, dest::CountedDestination) = CountedDestination(
            Metal.Adapt.adapt(to, dest.values),
            Metal.Adapt.adapt(to, dest.stores),
        )
        @inline Base.getindex(dest::CountedDestination, index) = @inbounds dest.values[index]
        @inline function Base.setindex!(dest::CountedDestination, value, index)
            Metal.atomic_fetch_add_explicit(pointer(dest.stores, 1), UInt32(1))
            @inbounds dest.values[index] = value
            return dest
        end
        @eval @testset "Metal barrier and guarded commits perform no rejected stores" begin
            ext = Base.get_extension(HybridPlasmaPIC, :HybridPlasmaPICMetalExt)
            for count in (1, 257), midpoint in (false, true), invalid in (false, true)
                cpu = population(2, count)
                stores = Metal.zeros(UInt32, 1)
                counted(a) = CountedDestination(Metal.MtlArray(a), stores)
                x, v = map(counted, cpu.x), map(counted, cpu.v)
                mid = midpoint ? ntuple(_ -> counted(fill(-1.0f0, count)), 2) : nothing
                E, B = fields(count)
                invalid && (B[3][end] = NaN32)
                Eg, Bg = map(Metal.MtlArray, E), map(Metal.MtlArray, B)
                flag = Metal.MtlArray{UInt32,1,Metal.SharedStorage}(UInt32[0])
                # Three threads force repeated strides and inactive lanes while
                # the entire scan remains inside the required one-group grid.
                Metal.@metal threads = 3 groups = 1 ext._push_gathered_single_group_kernel!(
                    x,
                    v,
                    Eg,
                    Bg,
                    0.75f0,
                    0.0f0,
                    mid,
                    flag,
                    UInt64(count),
                )
                @test only(Array(flag)) == UInt32(invalid)
                expected_stores = invalid ? 0 : count * (5 + (midpoint ? 2 : 0))
                @test only(Array(stores)) == UInt32(expected_stores)
                @test bits(map(a -> a.values, x)) == bits(cpu.x)
                @test bits(map(a -> a.values, v)) == bits(cpu.v)
                if midpoint
                    expected = invalid ? ntuple(_ -> fill(-1.0f0, count), 2) : cpu.x
                    @test bits(map(a -> a.values, mid)) == bits(expected)
                end
            end
            for status in (UInt32(0), UInt32(1), typemax(UInt32))
                values = Metal.MtlArray(fill(-0.0f0, 7))
                stores = Metal.zeros(UInt32, 1)
                dest = CountedDestination(values, stores)
                flag = Metal.MtlArray([status])
                source = Metal.MtlArray(fill(1.0f0, 7))
                expression =
                    Base.Broadcast.instantiate(Base.Broadcast.broadcasted(identity, source))
                Metal.@metal threads = 3 groups = 1 ext._gathered_commit_kernel!(
                    dest,
                    expression,
                    flag,
                    UInt64(7),
                )
                @test only(Array(stores)) == (iszero(status) ? UInt32(7) : UInt32(0))
                @test bits((values,)) == bits((fill(iszero(status) ? 1.0f0 : -0.0f0, 7),))
            end
        end
    else
        @test_skip "Metal gathered validation requires a functional remote Metal device"
    end
end
end
