using HybridPlasmaPIC, Test

const _MEMORY_COUNTER_FIELDS =
    (:total_bytes, :free_bytes, :used_bytes, :cached_bytes, :reserved_bytes)

function _construct_memory_status(style, counters; note = "")
    if style === :wrapper
        return BackendMemoryStatus(:cpu, true, false; counters..., note)
    elseif style === :keyword
        return BackendMemoryStatus(;
            backend = :cpu,
            device_available = true,
            pool_supported = false,
            counters...,
            note,
        )
    end
    return BackendMemoryStatus(:cpu, true, false, values(counters)..., note)
end

@testset "Backend memory counters share a constructor invariant" begin
    valid = (
        nothing,
        false,
        true,
        Int8(7),
        UInt8(9),
        0,
        1,
        typemax(Int),
        Int128(typemax(Int)),
        UInt128(typemax(Int)),
        big(17),
        big(typemax(Int)),
    )
    invalid = (
        -1,
        typemin(Int),
        Int128(typemax(Int)) + 1,
        UInt128(typemax(Int)) + 1,
        big(typemax(Int)) + 1,
        big(-1),
        -1.0,
        0.0,
        1.0,
        0.5,
        NaN,
        Inf,
        -Inf,
        1 // 1,
        1 // 2,
        1 + 0im,
        missing,
        "1",
    )
    for style in (:wrapper, :keyword, :positional), field in _MEMORY_COUNTER_FIELDS
        for value in valid
            counters = NamedTuple{_MEMORY_COUNTER_FIELDS}(
                ntuple(i -> _MEMORY_COUNTER_FIELDS[i] === field ? value : nothing, 5),
            )
            status = _construct_memory_status(style, counters)
            expected = value === nothing ? nothing : Int(value)
            @test getproperty(status, field) === expected
            @test status.backend === :cpu
            @test status.device_available && !status.pool_supported
            @test status.note == ""
            @test all(
                getproperty(status, other) === nothing for
                other in _MEMORY_COUNTER_FIELDS if other !== field
            )
        end
        for value in invalid
            counters = NamedTuple{_MEMORY_COUNTER_FIELDS}(
                ntuple(i -> _MEMORY_COUNTER_FIELDS[i] === field ? value : nothing, 5),
            )
            @test_throws ArgumentError _construct_memory_status(style, counters)
        end
    end
end

@testset "Backend memory defaults and metadata conversions" begin
    defaults = (
        BackendMemoryStatus(:cpu, true, false),
        BackendMemoryStatus(; backend = :cpu, device_available = true, pool_supported = false),
    )
    for status in defaults
        @test all(getproperty(status, field) === nothing for field in _MEMORY_COUNTER_FIELDS)
        @test status.note == ""
        @test memory_pressure(status) === nothing
    end
    counters = (
        total_bytes = UInt16(100),
        free_bytes = Int8(25),
        used_bytes = big(75),
        cached_bytes = nothing,
        reserved_bytes = Int128(125),
    )
    for style in (:wrapper, :keyword, :positional)
        status = _construct_memory_status(style, counters; note = SubString("telemetry", 2, 5))
        @test (
            status.total_bytes,
            status.free_bytes,
            status.used_bytes,
            status.cached_bytes,
            status.reserved_bytes,
        ) === (100, 25, 75, nothing, 125)
        @test status.note == "elem"
        @test memory_pressure(status) === 0.75
    end
    for status in (
        BackendMemoryStatus(:cpu, 1, 0, nothing, nothing, nothing, nothing, nothing, ""),
        BackendMemoryStatus(; backend = :cpu, device_available = 1, pool_supported = 0),
    )
        @test status.device_available === true
        @test status.pool_supported === false
    end
    @test_throws InexactError BackendMemoryStatus(
        :cpu,
        2,
        false,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        "",
    )
    @test_throws InexactError BackendMemoryStatus(;
        backend = :cpu,
        device_available = 2,
        pool_supported = false,
    )
    @test_throws MethodError BackendMemoryStatus("cpu", true, false)
    @test_throws MethodError BackendMemoryStatus(;
        backend = "cpu",
        device_available = true,
        pool_supported = false,
    )
end

@testset "Independent memory samples and pressure" begin
    for style in (:wrapper, :keyword, :positional)
        for (total, free, used, expected) in (
            (nothing, 2, 3, nothing),
            (0, 2, 3, nothing),
            (100, nothing, nothing, nothing),
            (100, 75, nothing, 0.25),
            (100, 150, nothing, 0.0),
            (100, 150, 125, 1.0),
            (100, 150, 25, 0.25),
            (typemax(Int), 0, typemax(Int), 1.0),
        )
            counters = (
                total_bytes = total,
                free_bytes = free,
                used_bytes = used,
                cached_bytes = typemax(Int),
                reserved_bytes = 0,
            )
            status = _construct_memory_status(style, counters)
            @test memory_pressure(status) === expected
            @test status.free_bytes === free
            @test status.used_bytes === used
            @test status.cached_bytes === typemax(Int)
            @test status.reserved_bytes === 0
        end
    end
end
