using HybridPlasmaPIC, Test, Serialization

# This wrapper's iteration and stored dictionary expose the same objects. Only
# stored edges count as repeated serialization of an immutable view.
struct IOSerializationDict{T} <: AbstractDict{Symbol,T}
    data::Dict{Symbol,T}
end
Base.length(d::IOSerializationDict) = length(d.data)
Base.iterate(d::IOSerializationDict, state...) = iterate(d.data, state...)
Base.getindex(d::IOSerializationDict, key) = d.data[key]

struct IOSerializationPointer
    pointer::Ptr{Cvoid}
end

io_serialization_copy_value(value::BigInt) = deepcopy(value)
io_serialization_copy_value(value::Regex) =
    Regex(value.pattern, value.compile_options, value.match_options)

function io_serialization_invalid(::Type{T}) where {T}
    a = T[0.1, 0.2, 0.3, 0.4]
    v = view(a, 2:3)
    index = [2, 3]
    indexed = view(a, index)
    wrapper = IOOwnershipArray(v, Ref(zero(T)), nothing)
    dictionary = Dict(:temperature => T(0.1))
    identity_dictionary = IdDict{Any,Any}(:temperature => T(0.1))
    members = Set(T[0.1])
    expression = Expr(:tuple, copy(a))
    sibling = Expr(:tuple)
    sibling.args = expression.args
    integer = big(12)
    pattern = Regex("parameter")
    pointer = Ptr{Cvoid}(UInt(0x1234)) # Inspected only, never dereferenced.
    return (
        (v, v),
        (a, v),
        (v, a),
        (view(a, 1:2), view(a, 2:3)),
        (a, reshape(a, 2, 2)),
        (index, indexed),
        (indexed, index),
        (indexed, view(copy(a), index)),
        view(reshape(T.(1:16), 4, 4), index, index),
        (wrapper, wrapper),
        (Ref(wrapper), wrapper),
        ((view = v,), (view = v,)),
        (dictionary, dictionary.slots),
        (dictionary.keys, dictionary),
        (dictionary, view(dictionary.vals, :)),
        (identity_dictionary, identity_dictionary.ht),
        (members, members.dict.vals),
        (expression, expression.args),
        (expression.args, expression),
        (expression, view(expression.args, :)),
        (expression, sibling),
        (sibling, expression),
        (integer, integer),
        (pattern, pattern),
        pointer,
        Ref(pointer),
        [pointer],
        IOSerializationPointer(pointer),
        [IOSerializationPointer(pointer)],
        view([pointer, C_NULL], 1:1),
    )
end

@testset "Mutable value codecs preserve shared containing objects" begin
    source, args, config = io_ownership_run(HybridStepper, Float64, Val(1), false)
    for value in (big(12), Regex("parameter"))
        for extra in (
            Any[value, value],
            (Ref(value), Ref(value)),
            ((parameter = value,), (parameter = value,)),
            Dict(:first => value, :second => value),
            Expr(:tuple, value, value),
        )
            config.extra = extra
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(source, args),
                source,
            )
        end
        for container in (
            Ref(io_serialization_copy_value(value)),
            [io_serialization_copy_value(value)],
            Dict(:parameter => io_serialization_copy_value(value)),
        )
            config.extra = (container, container)
            mktempdir() do directory
                checkpoint, archive =
                    joinpath(directory, "checkpoint"), joinpath(directory, "archive")
                save_checkpoint(checkpoint, source, args...)
                archive_run(archive, source, args...; rng_seed = 12)
                for state in (deserialize(checkpoint), load_archive(archive).state)
                    first, second = io_ownership_config(state).extra
                    @test first === second
                    @test first !== container
                    getvalue(c) = c isa Ref ? c[] : c isa AbstractDict ? c[:parameter] : c[1]
                    if value isa BigInt
                        Base.GMP.MPZ.add!(getvalue(first), big(1))
                        @test getvalue(second) == 13
                    else
                        getvalue(first).match_options |= Base.PCRE.ANCHORED
                        @test match(getvalue(second), "xparameter") === nothing
                    end
                    @test getvalue(container) == value
                end
            end
        end
        config.extra = (io_serialization_copy_value(value), io_serialization_copy_value(value))
        @test config.extra[1] !== config.extra[2]
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(io_ownership_state(source, args), source),
        )
    end
end

@testset "Independent view indices retain logical mutation behavior" begin
    for T in (Float32, Float64), mutable_indices in (false, true)
        source, args, config = io_ownership_run(HybridStepper, T, Val(1), false)
        original = if mutable_indices
            view(reshape(T.(1:16), 4, 4), [1, 2], [1, 2])
        else
            view(reshape(T.(1:16), 4, 4), 1:2, 1:2)
        end
        config.extra = original
        mktempdir() do directory
            checkpoint, archive = joinpath(directory, "checkpoint"), joinpath(directory, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 13)
            restored = map((deserialize(checkpoint), load_archive(archive).state)) do state
                value = io_ownership_config(state).extra
                @test value == original
                value
            end
            if mutable_indices
                for value in (original, restored...)
                    @test parentindices(value)[1] !== parentindices(value)[2]
                    parentindices(value)[1][1] = 2
                end
                @test all(value -> value == original, restored)
            end
            for value in (original, restored...)
                value[1, 2] = T(99)
            end
            @test all(value -> value == original, restored)
        end
    end
end

@testset "Configuration serialization rejects detached storage before publication or restore" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES, D = 1:3
        State === EMPIC1D && D != 1 && continue
        source, args, config = io_ownership_run(State, T, Val(D), variant)
        target, destination, target_config = io_ownership_run(State, T, Val(D), variant)
        target_before = io_ownership_bytes((
            io_ownership_state(target, destination),
            io_ownership_storage(target, destination),
        ))
        incoming = deepcopy(io_ownership_state(source, args))
        mktempdir() do directory
            checkpoint = joinpath(directory, "checkpoint")
            archive = joinpath(directory, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 10)
            checkpoint_before, archive_before = read(checkpoint), read(archive)
            for invalid in io_serialization_invalid(T)
                config.extra = invalid
                @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                @test read(checkpoint) == checkpoint_before
                @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 10)
                @test read(archive) == archive_before
                io_ownership_config(incoming).extra = invalid
                @test_throws ArgumentError HybridPlasmaPIC._restore_checkpoint_state!(
                    target,
                    destination...,
                    incoming,
                )
                @test io_ownership_config(target) === target_config
                @test io_ownership_bytes((
                    io_ownership_state(target, destination),
                    io_ownership_storage(target, destination),
                )) == target_before
            end
            # Pointer arrays retain their addresses in Julia's raw array encoding;
            # public restore must reject such incoming storage before mutation.
            io_ownership_config(incoming).extra = [Ptr{Cvoid}(UInt(0x1234))]
            malformed = joinpath(directory, "malformed")
            serialize(malformed, incoming)
            @test_throws ArgumentError load_checkpoint!(target, destination..., malformed)
            @test io_ownership_config(target) === target_config
            @test io_ownership_bytes((
                io_ownership_state(target, destination),
                io_ownership_storage(target, destination),
            )) == target_before
        end
    end
end

@testset "Serialization edge counts distinguish wrappers, cycles and disjoint storage" begin
    for T in (Float32, Float64)
        source, args, config = io_ownership_run(HybridStepper, T, Val(1), false)
        a = T[0.1, 0.2, 0.3, 0.4]
        v = view(a, 2:3)
        wrapper = IOOwnershipArray(v, Ref(zero(T)), nothing)
        mutable_wrapper = Ref(wrapper)
        cycle = Any[]
        push!(cycle, cycle, a, a)
        bits = BitArray([true, false, true, false])
        bit_views = (view(bits, 1:2), view(bits, 3:4), bits.chunks)
        members = Set(T[0.1])
        expression = Expr(:tuple, a)
        positive = (
            (a, a),
            v,
            view(a, [2, 3]),
            view(reshape(T.(1:16), 4, 4), [1, 2], [1, 2]),
            view(reshape(T.(1:16), 4, 4), 1:2, 1:2),
            (view(a, 1:2), view(a, 3:4)),
            (view(T[], :), view(T[], :)),
            wrapper,
            IOSerializationDict(Dict(:parameters => v)),
            (mutable_wrapper, mutable_wrapper),
            cycle,
            bit_views,
            (members, members.dict),
            expression,
            (expression, expression),
            (expression, a),
            (expression, Expr(:tuple, a)),
            Expr(:tuple, v),
            (C_NULL, Ref(C_NULL), [C_NULL], IOSerializationPointer(C_NULL)),
            view([Ptr{Cvoid}(UInt(0x1234)), C_NULL], 2:2),
            (big(12), BigFloat(0.125), r"parameter"),
        )
        mktempdir() do directory
            checkpoint = joinpath(directory, "checkpoint")
            for extra in positive
                config.extra = extra
                @test isnothing(
                    HybridPlasmaPIC._validate_checkpoint_state(
                        io_ownership_state(source, args),
                        source,
                    ),
                )
                save_checkpoint(checkpoint, source, args...)
                saved = deserialize(checkpoint)
                @test isnothing(HybridPlasmaPIC._validate_checkpoint_state(saved, source))
            end
        end
        # Metadata-first and logical-first visits must give the same result.
        for extra in ((wrapper, v), (v, wrapper), (a, v), (v, a))
            config.extra = extra
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(source, args),
                source,
            )
        end
        config.extra = IOSerializationDict(Dict(:first => v, :second => v))
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(source, args),
            source,
        )
        # Unused table capacity is not serialized and cannot introduce a native
        # resource or an array relationship into the logical configuration.
        dictionary = Dict{Symbol,Any}(:old => Ref(1))
        empty!(dictionary)
        dictionary.vals[1] = Ptr{Cvoid}(UInt(0x1234))
        config.extra = dictionary
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(io_ownership_state(source, args), source),
        )
        # Immutable fan-out is counted up to two, not expanded exponentially.
        nested = (a,)
        for _ = 1:12
            nested = (nested, nested)
        end
        config.extra = nested
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(io_ownership_state(source, args), source),
        )
    end
end

@testset "Independent configuration preserves mutable continuation and safe value codecs" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES, D = 1:3
        State === EMPIC1D && D != 1 && continue
        source, args, config = io_ownership_run(State, T, Val(D), variant)
        parameters = T[0.1, 0.2, 0.3, 0.4]
        if config isa IOOwnershipClosure
            config.temperature = view(parameters, 1:2)
        else
            config.parameters = view(parameters, 1:2)
        end
        shared = T[0.5, 0.6]
        cycle = Any[]
        push!(cycle, cycle, shared)
        config.extra = (
            other = view(parameters, 3:4),
            first = shared,
            second = shared,
            cycle = cycle,
            values = (big(12), BigFloat(0.125), r"parameter"),
        )
        io_ownership_initialize!(source, args)
        io_ownership_step!(source, args, T(0.001))
        target, destination, _ = io_ownership_run(State, T, Val(D), variant; count = 1)
        mktempdir() do directory
            checkpoint, archive = joinpath(directory, "checkpoint"), joinpath(directory, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 11)
            load_checkpoint!(target, destination..., checkpoint)
            archived = load_archive(archive).state
            io_ownership_compare_state(archived, io_ownership_state(source, args))
            for restored in (io_ownership_config(target), io_ownership_config(archived))
                @test restored.extra.first === restored.extra.second
                @test restored.extra.cycle[1] === restored.extra.cycle
                @test restored.extra.cycle[2] === restored.extra.first
                @test restored.extra.other == T[0.3, 0.4]
                @test restored.extra.values[1] == big(12)
                @test restored.extra.values[2] == BigFloat(0.125)
                @test match(restored.extra.values[3], "parameter").match == "parameter"
            end
            for c in (config, io_ownership_config(target))
                parameter = c isa IOOwnershipClosure ? c.temperature : c.parameters
                parameter[1] = T(0.15)
                c.extra.first[1] = T(0.7)
                @test c.extra.second[1] == T(0.7)
                @test c.extra.other == T[0.3, 0.4]
            end
            io_ownership_step!(source, args, T(0.002))
            io_ownership_step!(target, destination, T(0.002))
            io_ownership_compare_state(
                io_ownership_state(source, args),
                io_ownership_state(target, destination),
            )
        end
    end
end

@testset "Independent native plans reject before publication or public restore" begin
    for T in (Float32, Float64),
        State in (HybridStepper, CAMCLStepper, EMPIC1D, EMPIC, HallMHDState)

        source, args, config = io_ownership_run(State, T, Val(1), false)
        target, destination, target_config = io_ownership_run(State, T, Val(1), false)
        target_before = io_ownership_bytes((
            io_ownership_state(target, destination),
            io_ownership_storage(target, destination),
        ))
        independent = FourierGrid((4,), (T(2π),))
        mktempdir() do directory
            checkpoint, archive = joinpath(directory, "checkpoint"), joinpath(directory, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 12)
            before, archive_before = read(checkpoint), read(archive)
            for plan in (independent.plan, independent.iplan)
                config.extra = plan
                @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                @test read(checkpoint) == before
                @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 12)
                @test read(archive) == archive_before
                # Legacy serialized plans retain their type but lose the handle.
                # Reject that identifiable invalid resource too; never execute it.
                malformed = joinpath(directory, "malformed")
                serialize(malformed, io_ownership_state(source, args))
                @test_throws ArgumentError load_checkpoint!(target, destination..., malformed)
                @test io_ownership_config(target) === target_config
                @test io_ownership_bytes((
                    io_ownership_state(target, destination),
                    io_ownership_storage(target, destination),
                )) == target_before
            end
        end
    end
end

@testset "Reconstructed configuration buffers cannot own evolving particle storage" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES, D = 1:3
        (State === HallMHDState || (State === EMPIC1D && D != 1)) && continue
        dictionary = Dict(:parameter => T(0.1))
        source, args, config =
            io_ownership_run(State, T, Val(D), variant; count = length(dictionary.vals))
        particles = first(io_ownership_particles(args))
        dictionary.vals .= particles.x[1]
        replacement = ParticleSet{D,T}(
            ntuple(
                d ->
                    d == 1 ? dictionary.vals :
                    copyto!(Memory{T}(undef, length(dictionary.vals)), particles.x[d]),
                D,
            ),
            particles.v,
            particles.weight,
            particles.id,
            particles.tag,
            particles.q,
            particles.m,
        )
        if first(args) isa AbstractVector
            populations = ParticleSet{D,T}[ps for ps in first(args)]
            populations[1] = replacement
            args = (populations,)
        else
            args = (replacement, Base.tail(args)...)
        end
        # Dict uses Memory on supported Julia versions. An Array assignment to
        # this field would convert/copy and would not reproduce the lost link.
        @test dictionary.vals === first(io_ownership_particles(args)).x[1]
        slot = only(findall(!iszero, dictionary.slots))
        original = dictionary[:parameter]
        replacement.x[1][slot] += one(T)
        @test dictionary[:parameter] == original + one(T)
        replacement.x[1][slot] = original
        target, destination, target_config = io_ownership_run(State, T, Val(D), variant)
        target_before = io_ownership_bytes((
            io_ownership_state(target, destination),
            io_ownership_storage(target, destination),
        ))
        mktempdir() do directory
            checkpoint, archive = joinpath(directory, "checkpoint"), joinpath(directory, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 13)
            checkpoint_before, archive_before = read(checkpoint), read(archive)
            config.extra = dictionary
            incoming = io_ownership_state(source, args)
            before = io_ownership_bytes(incoming)
            @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
            @test read(checkpoint) == checkpoint_before
            @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 13)
            @test read(archive) == archive_before
            @test_throws ArgumentError HybridPlasmaPIC._restore_checkpoint_state!(
                target,
                destination...,
                incoming,
            )
            @test io_ownership_bytes(incoming) == before
            @test io_ownership_config(target) === target_config
            @test io_ownership_bytes((
                io_ownership_state(target, destination),
                io_ownership_storage(target, destination),
            )) == target_before
        end
    end
end
