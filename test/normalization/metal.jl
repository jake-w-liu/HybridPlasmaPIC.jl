module NormalizationMetalTests
using HybridPlasmaPIC, Test
const H = HybridPlasmaPIC
if Sys.isapple() && Base.find_package("Metal") !== nothing
    import Metal
    if Metal.functional()
        Metal.allowscalar(false)
        @testset "Native device singleton normalization at IEEE boundaries" begin
            for (T, halfpower) in ((Float16, 5), (Float32, 50)),
                storage in (Metal.PrivateStorage, Metal.SharedStorage),
                direction in (-1, 1),
                inverse in (false, true)

                @testset "$T $storage sign=$direction inverse=$inverse" begin
                    B = ldexp(one(T), direction * halfpower)
                    u = PlasmaUnits{T}(1, B, 1, 1, 1, 1)
                    values = T[
                        0,
                        -0.0,
                        nextfloat(zero(T)),
                        -nextfloat(zero(T)),
                        nextfloat(nextfloat(zero(T))),
                        prevfloat(floatmin(T)),
                        floatmin(T),
                        -floatmin(T),
                        floatmax(T),
                        -floatmax(T),
                        Inf,
                        -Inf,
                        NaN,
                    ]
                    for value in values
                        singleton = T[value]
                        source = Metal.MtlArray{T,1,storage}(singleton)
                        shift = (inverse ? -1 : 1) * direction * 2halfpower
                        # A positive power-of-two physical scale makes ldexp an independent
                        # one-round reference, including signs, subnormals and overflow endpoints.
                        expected = ldexp.(singleton, shift)
                        actual =
                            inverse ? to_normalized(source, :pressure, u) :
                            to_SI(source, :pressure, u)
                        @test typeof(actual) === typeof(source)
                        @test isequal(Array(actual), expected)
                        @test isequal(Array(source), singleton)
                    end
                end
            end
        end
        @testset "Metal normalization storage and layout" begin
            for storage in (Metal.PrivateStorage, Metal.SharedStorage), ordinary in (false, true)
                B = ordinary ? 2.0f0 : ldexp(1.0f0, 100)
                u = PlasmaUnits{Float32}(1, B, 1, 1, 1, 1)
                x = ldexp(1.0f0, -100)
                base = Metal.MtlArray{Float32,1,storage}(Float32[x, -x, 0, -0.0])
                matrix = reshape(base, 2, 2)
                sources = (
                    base,
                    matrix,
                    view(base, 1:3),
                    view(base, 1:2:4),
                    transpose(matrix),
                    view(matrix, :, 1:2:2),
                    Metal.MtlArray{Float32,1,storage}(Float32[]),
                    Metal.MtlArray{Float32,2,storage}(zeros(Float32, 0, 3)),
                    Metal.MtlArray{Float32,0,storage}(fill(x)),
                )
                for (index, source) in enumerate(sources), inverse in (false, true)
                    @testset "storage=$storage ordinary=$ordinary index=$index inverse=$inverse" begin
                        before = Array(source)
                        expected =
                            inverse ? to_normalized(before, :pressure, u) :
                            to_SI(before, :pressure, u)
                        actual =
                            inverse ? to_normalized(source, :pressure, u) :
                            to_SI(source, :pressure, u)
                        host = actual isa AbstractArray ? Array(actual) : actual
                        # Zero-dimensional conversion follows Julia's scalar broadcast contract.
                        observed = host isa AbstractArray && ndims(host) == 0 ? host[] : host
                        @test isequal(observed, expected)
                        @test isequal(Array(source), before)
                        @test !(actual isa Array)
                        @test ndims(source) == 0 ? actual isa Float32 :
                              actual isa Metal.MtlArray{Float32,ndims(source),storage}
                    end
                end
            end
        end
        @testset "Native normalization flag boundaries and lifetime" begin
            for T in (Float16, Float32), storage in (Metal.PrivateStorage, Metal.SharedStorage)
                @testset "$T $storage" begin
                    for count in (0, 1, 2, 31, 32, 33, 255, 256, 257, 1023, 1024, 1025, 4097, 65537)
                        input = fill(one(T), count)
                        output = copy(input)
                        a = Metal.MtlArray{T,1,storage}(input)
                        b = Metal.MtlArray{T,1,storage}(output)
                        @test !H._normalization_device_needs_complete(a, b)
                        for index in unique(
                            filter(
                                i -> 1 <= i <= count,
                                (1, count, div(count, 2) + 1, 32, 33, 256, 257, 1024, 1025),
                            ),
                        )
                            output[index] = zero(T)
                            copyto!(b, output)
                            @test H._normalization_device_needs_complete(a, b)
                            output[index] = one(T)
                            copyto!(b, output)
                            @test !H._normalization_device_needs_complete(a, b)
                        end
                        if count > 0
                            fill!(b, T(Inf))
                            @test H._normalization_device_needs_complete(a, b)
                            # All threads may set the flag concurrently; a later call must start clean.
                            fill!(b, one(T))
                            @test !H._normalization_device_needs_complete(a, b)
                        end
                        @test Array(a) == input
                    end
                    input = T[0, -0.0, Inf, -Inf, NaN]
                    a = Metal.MtlArray{T,1,storage}(input)
                    b = Metal.MtlArray{T,1,storage}(fill(T(NaN), 5))
                    @test !H._normalization_device_needs_complete(a, b)
                    @test isequal(Array(a), input)
                    power = T === Float16 ? 10 : 100
                    u = PlasmaUnits{T}(1, ldexp(one(T), power), 1, 1, 1, 1)
                    x = Metal.MtlArray{T,1,storage}(fill(ldexp(one(T), -power), 129))
                    converted = to_SI(x, :pressure, u)
                    GC.gc(true)
                    @test Array(converted) == fill(ldexp(one(T), power), 129)
                    @test Array(x) == fill(ldexp(one(T), -power), 129)
                end
            end
        end
        @testset "Native complex normalization preserves component range and signs" begin
            for T in (Float16, Float32),
                storage in (Metal.PrivateStorage, Metal.SharedStorage),
                halfpower in (T === Float16 ? (5, 10) : (50, 100)),
                direction in (-1, 1),
                inverse in (false, true)

                u = PlasmaUnits{T}(1, ldexp(one(T), direction * halfpower), 1, 1, 1, 1)
                small = nextfloat(zero(T))
                normal = floatmin(T)
                values = Complex{T}[
                    complex(small, -small),
                    complex(normal, -normal),
                    complex(one(T), -one(T)),
                    complex(zero(T), -zero(T)),
                    complex(T(Inf), T(-Inf)),
                    complex(T(NaN), one(T)),
                ]
                shift = (inverse ? -1 : 1) * direction * 2halfpower
                expected = [complex(ldexp(real(z), shift), ldexp(imag(z), shift)) for z in values]
                source = Metal.MtlArray{Complex{T},1,storage}(values)
                actual = inverse ? to_normalized(source, :pressure, u) : to_SI(source, :pressure, u)
                @test typeof(actual) === typeof(source)
                @test isequal(Array(actual), expected)
                @test isequal(Array(source), values)
            end
        end
    end
end
end
