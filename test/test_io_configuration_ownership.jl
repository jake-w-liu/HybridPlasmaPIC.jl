using HybridPlasmaPIC, Test, Serialization

mutable struct IOOwnershipClosure{T} <: HybridPlasmaPIC.ElectronClosure
    temperature::AbstractVector{T}
    extra::Any
end

function HybridPlasmaPIC.electron_pressure!(pe, n, c::IOOwnershipClosure)
    temperature = c.temperature[1]
    @. pe = temperature * n
    return pe
end

mutable struct IOOwnershipShape{T} <: ShapeFunction
    parameters::AbstractVector{T}
    extra::Any
end
HybridPlasmaPIC.width(::IOOwnershipShape) = 2
HybridPlasmaPIC._stencil1d(::IOOwnershipShape, s) = HybridPlasmaPIC._stencil1d(CIC(), s)

struct IOOwnershipString{T} <: AbstractString
    text::String
    parameters::AbstractArray{T}
end
Base.ncodeunits(s::IOOwnershipString) = ncodeunits(s.text)
Base.codeunit(s::IOOwnershipString) = codeunit(s.text)
Base.codeunit(s::IOOwnershipString, i::Int) = codeunit(s.text, i)
Base.isvalid(s::IOOwnershipString, i::Int) = isvalid(s.text, i)
Base.iterate(s::IOOwnershipString, i::Int = 1) = iterate(s.text, i)

struct IOOwnershipArray{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    offset::Base.RefValue{T}
    extra::Any
end
Base.size(a::IOOwnershipArray) = size(a.data)
Base.IndexStyle(::Type{<:IOOwnershipArray}) = IndexLinear()
Base.parent(a::IOOwnershipArray) = a.data
Base.dataids(a::IOOwnershipArray) = Base.dataids(a.data)
# Preserve the declared parent's precise view overlap checks in this wrapper.
Base.mightalias(a::IOOwnershipArray, b::AbstractArray) = Base.mightalias(a.data, b)
Base.mightalias(a::AbstractArray, b::IOOwnershipArray) = Base.mightalias(a, b.data)
Base.mightalias(a::IOOwnershipArray, b::IOOwnershipArray) = Base.mightalias(a.data, b.data)
Base.getindex(a::IOOwnershipArray, i::Int) = a.data[i] + a.offset[]

struct IOOwnershipTrackedArray{T} <: AbstractVector{T}
    data::Vector{T}
    members::Set{T}
end
Base.size(a::IOOwnershipTrackedArray) = size(a.data)
Base.IndexStyle(::Type{<:IOOwnershipTrackedArray}) = IndexLinear()
Base.parent(a::IOOwnershipTrackedArray) = a.data
Base.dataids(a::IOOwnershipTrackedArray) = Base.dataids(a.data)
Base.getindex(a::IOOwnershipTrackedArray, i::Int) = a.data[i]
function Base.setindex!(a::IOOwnershipTrackedArray, value, i::Int)
    a.data[i] = value
    empty!(a.members)
    union!(a.members, a.data)
    return value
end

struct IOOwnershipDict{T} <: AbstractDict{Symbol,T}
    data::Dict{Symbol,T}
    offset::Base.RefValue{T}
    extra::Any
end
Base.length(d::IOOwnershipDict) = length(d.data)
Base.getindex(d::IOOwnershipDict, key) = d.data[key] + d.offset[]
function Base.iterate(d::IOOwnershipDict, state...)
    item = iterate(d.data, state...)
    item === nothing && return nothing
    pair, next = item
    return (first(pair) => last(pair) + d.offset[], next)
end

struct IOOwnershipSet{T} <: AbstractSet{T}
    data::Set{T}
    offset::Base.RefValue{T}
    extra::Any
end
Base.length(s::IOOwnershipSet) = length(s.data)
Base.in(value, s::IOOwnershipSet) = any(isequal(value), s)
function Base.iterate(s::IOOwnershipSet, state...)
    item = iterate(s.data, state...)
    item === nothing && return nothing
    value, next = item
    return (value + s.offset[], next)
end

function io_ownership_wrappers(::Type{T}, offset, extra = nothing) where {T}
    array = IOOwnershipArray(T[0.1, 0.2], offset, extra)
    return (
        array,
        view(array, 1:1),
        reshape(view(array, 1:1), 1, 1),
        IOOwnershipDict(Dict(:temperature => T(0.1)), offset, extra),
        IOOwnershipSet(Set(T[0.1]), offset, extra),
    )
end

function io_ownership_run(State, ::Type{T}, ::Val{D}, variant; count = 3^D) where {T,D}
    g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> T(2π), D))
    closure = IOOwnershipClosure(T[0.1], nothing)
    shape = IOOwnershipShape(T[0.1], nothing)
    if State === HallMHDState
        st = State(g, HallMHDModel(closure))
        fill!(st.fields.n, one(T))
        fill!(st.fields.B[1], one(T))
        return st, (), closure
    end
    function population(q, m, density)
        ps = ParticleSet{D,T}(count; q, m)
        if count == prod(g.n)
            # Continuation fixtures resolve every cell; smaller counts are
            # destination placeholders that are replaced before stepping.
            load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, g.n)
        else
            for d = 1:D
                ps.x[d] .= T.(1:count) .* T(0.1d)
            end
        end
        count > 0 && set_density_weight!(ps, T(density), g)
        foreach(a -> fill!(a, T(0.001)), ps.v)
        ps.v[1] .= T(0.001) .* (one(T) .+ sin.(ps.x[1]))
        return ps
    end
    if State in (HybridStepper, CAMCLStepper)
        st = State(g, HybridModel(closure), shape, count)
        fill!(st.fields.B[1], one(T))
        ps = population(one(T), one(T), variant ? 0.5 : 1)
        input = variant ? [ps, population(T(2), T(4), 0.25)] : ps
        return st, (input,), closure
    end
    st = State(g, count; shape, mobile = variant, mi = T(16), c = T(3))
    electrons = population(-one(T), one(T), 1)
    ions = variant ? population(one(T), T(16), 1) : nothing
    return st, (electrons, ions), shape
end

io_ownership_initialize!(st::HybridStepper, args) = init!(st, args...)
io_ownership_initialize!(st::CAMCLStepper, args) = init_camcl!(st, args...)
io_ownership_initialize!(st::Union{EMPIC1D,EMPIC}, args) = init_empic!(st, args...)
io_ownership_initialize!(st::HallMHDState, args) = nothing
io_ownership_step!(st::HybridStepper, args, dt) = step!(st, args..., dt; NB = 2)
io_ownership_step!(st::CAMCLStepper, args, dt) = step_camcl!(st, args..., dt; NB = 2)
io_ownership_step!(st::Union{EMPIC1D,EMPIC}, args, dt) = step_empic!(st, args..., dt)
io_ownership_step!(st::HallMHDState, args, dt) = step_hall_mhd!(st, dt)

function io_ownership_particles(args)
    isempty(args) && return ()
    first(args) isa AbstractVector && return first(args)
    return filter(!isnothing, args)
end

io_ownership_config(st::Union{HybridStepper,CAMCLStepper,HallMHDState}) = st.model.closure
io_ownership_config(st::Union{EMPIC1D,EMPIC}) = st.shape
io_ownership_config(s::NamedTuple) = haskey(s, :model) ? s.model.closure : s.shape
io_ownership_state(st, args) = HybridPlasmaPIC._checkpoint_state(st, args...)
function io_ownership_bytes(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

function io_ownership_compare_state(a, b)
    @test keys(a) == keys(b)
    for name in keys(a)
        name in (:model, :shape) && continue
        @test isequal(getproperty(a, name), getproperty(b, name))
    end
end

function io_ownership_validation_error(s, st)
    try
        HybridPlasmaPIC._validate_checkpoint_state(s, st)
    catch err
        return err
    end
    return nothing
end

const IO_OWNERSHIP_FAMILIES = (
    (HybridStepper, false),
    (HybridStepper, true),
    (CAMCLStepper, false),
    (CAMCLStepper, true),
    (EMPIC1D, false),
    (EMPIC1D, true),
    (EMPIC, false),
    (EMPIC, true),
    (HallMHDState, false),
)

@testset "CAMCL configuration validation retains magnetic-field error order" begin
    for T in (Float32, Float64), D = 1:3, collection in (false, true)
        st, args, _ = io_ownership_run(CAMCLStepper, T, Val(D), collection)
        saved = io_ownership_state(st, args)
        for name in (:Blead, :Blag, :Bout)
            malformed = merge(saved, NamedTuple{(name,)}((nothing,)))
            err = io_ownership_validation_error(malformed, st)
            @test err isa ArgumentError
            @test sprint(showerror, err) ==
                  "ArgumentError: checkpoint $name fields must be stored as a tuple"
            malformed = merge(
                saved,
                NamedTuple{(name,)}(((
                    zeros(T, ntuple(_ -> 2, D)),
                    getproperty(saved, name)[2],
                    getproperty(saved, name)[3],
                ),)),
            )
            err = io_ownership_validation_error(malformed, st)
            @test err isa ArgumentError
            @test occursin("checkpoint $(name)1 axes", sprint(showerror, err))
        end
        malformed = merge(saved, (Blead = nothing, Blag = nothing, Bout = nothing))
        err = io_ownership_validation_error(malformed, st)
        @test err isa ArgumentError
        @test sprint(showerror, err) ==
              "ArgumentError: checkpoint Blead fields must be stored as a tuple"
        missing = (; (k => v for (k, v) in pairs(malformed) if k != :cl_h)...)
        err = io_ownership_validation_error(missing, st)
        @test err isa ArgumentError
        @test sprint(showerror, err) == "ArgumentError: checkpoint is missing required fields: cl_h"
    end
end

function io_ownership_storage(st, args)
    arrays = Any[st.g.cbuf, st.g.tbuf, st.g.abuf, st.g.sbuf]
    for ps in io_ownership_particles(args)
        append!(arrays, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
    end
    if st isa Union{HybridStepper,CAMCLStepper}
        append!(arrays, (st.fields.B..., st.fields.E..., st.fields.n, st.fields.pe))
        append!(arrays, (st.fn, st.fui..., st.Ep..., st.Bp..., st.vpred..., st.xmid..., st.work))
        append!(arrays, (st.time, st.step, st.last_dt, st.fields.floor_count))
        if st isa CAMCLStepper
            append!(arrays, (st.Blead..., st.Blag..., st.Bout..., st.cl_h, st.cl_live))
        end
    elseif st isa HallMHDState
        append!(
            arrays,
            (st.fields.n, st.fields.ui..., st.fields.B..., st.fields.E..., st.fields.J...),
        )
        append!(arrays, (st.stage.n, st.stage.lapJ..., st.rhs_n, st.n0, st.ntmp, st.kn...))
        append!(arrays, (st.time, st.step, st.fields.floor_count))
    elseif st isa EMPIC1D
        append!(arrays, (st.Ex, st.Ey, st.Bz, st.Jx, st.Jy, st.rho_n, st.rho_np1))
        append!(arrays, (st.Bzc, st.ne, st.Exp, st.Exi, st.work, st.worki))
        append!(arrays, (st.time, st.step, st.last_dt, st.last_sub_dt))
    else
        append!(arrays, (st.E..., st.B..., st.J..., st.rho_n, st.rho_np1))
        append!(arrays, (st.Bc..., st.Jraw..., st.ne, st.Ep..., st.Epi..., st.mide..., st.worke))
        append!(arrays, (st.time, st.step, st.last_dt, st.last_sub_dt))
    end
    return arrays
end

@testset "Checkpoint configuration rejects evolving storage ownership" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES, D = 1:3
        State === EMPIC1D && D != 1 && continue
        st, args, config = io_ownership_run(State, T, Val(D), variant)
        storage = io_ownership_storage(st, args)
        for value in storage
            config.extra = value
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(st, args),
                st,
            )
        end
        for value in (Core.svec(first(storage)), IOOwnershipString("parameters", first(storage)))
            config.extra = value
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(st, args),
                st,
            )
        end
        # Cycles must terminate; a nested reference must not hide the ownership link.
        config.extra = Any[]
        push!(config.extra, config.extra, Dict(:nested => Ref(first(storage))))
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(st, args),
            st,
        )
        config.extra = nothing
        if st isa Union{HybridStepper,CAMCLStepper}
            st.shape.extra = st.fields.B[1]
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(st, args),
                st,
            )
            st.shape.extra = nothing
        end
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, st, args...)
            archive_run(archive, st, args...; rng_seed = 7)
            before_checkpoint, before_archive = read(checkpoint), read(archive)
            # Include both represented and omitted storage at the public I/O boundary.
            for value in (first(storage), last(storage), st, st.g.plan, st.g.iplan, st.g.iplan.p)
                config.extra = value
                @test_throws ArgumentError save_checkpoint(checkpoint, st, args...)
                @test read(checkpoint) == before_checkpoint
                @test_throws ArgumentError archive_run(archive, st, args...; rng_seed = 7)
                @test read(archive) == before_archive
            end
        end
    end
end

@testset "Checkpoint configuration detects wrappers before serialization erases links" begin
    for T in (Float32, Float64), State in (HybridStepper, CAMCLStepper), count in (0, 3)
        st, args, config = io_ownership_run(State, T, Val(1), false; count)
        ps = only(args)
        for value in
            (view(ps.x[1], :), reshape(view(ps.x[1], :), 1, :), reinterpret(UInt8, ps.x[1]))
            config.extra = value
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(st, args),
                st,
            )
        end
        if count == 0
            # Reshaping an empty Array creates a separate array header; later
            # vector growth detaches its storage and leaves this Matrix empty.
            config.extra = reshape(ps.x[1], 1, :)
            @test isnothing(
                HybridPlasmaPIC._validate_checkpoint_state(io_ownership_state(st, args), st),
            )
        end
        config.extra = nothing
        config.temperature = ps.x[1]
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(st, args),
            st,
        )
    end
end

@testset "Incoming configuration links reject before destination mutation" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES
        source, args, _ = io_ownership_run(State, T, Val(1), variant)
        target, destination, _ = io_ownership_run(State, T, Val(1), variant; count = 1)
        before = io_ownership_bytes(io_ownership_state(target, destination))
        original_config = io_ownership_config(target)
        buffers = io_ownership_storage(target, destination)
        contents = map(a -> a isa AbstractArray ? copy(a) : a[], buffers)
        saved = deepcopy(io_ownership_state(source, args))
        linked = if haskey(saved, :species)
            last(saved.species).x[1]
        elseif haskey(saved, :electrons)
            (saved.mobile ? saved.ions : saved.electrons).x[1]
        elseif haskey(saved, :x)
            saved.x[1]
        else
            saved.fields.n
        end
        for value in (linked, Core.svec(linked), IOOwnershipString("parameters", linked))
            io_ownership_config(saved).extra = value
            mktemp() do path, io
                close(io)
                serialize(path, saved)
                @test_throws ArgumentError load_checkpoint!(target, destination..., path)
                @test io_ownership_config(target) === original_config
                @test io_ownership_bytes(io_ownership_state(target, destination)) == before
                @test all(
                    a === b for (a, b) in zip(buffers, io_ownership_storage(target, destination))
                )
                @test all(
                    isequal(a isa AbstractArray ? a : a[], b) for (a, b) in zip(buffers, contents)
                )
            end
        end
    end
end

@testset "Independent configuration storage and disjoint particle views remain supported" begin
    for T in (Float32, Float64), State in (HybridStepper, CAMCLStepper)
        st, args, config = io_ownership_run(State, T, Val(1), false)
        ps = only(args)
        backing = T[ps.x[1]..., 0.1, 0.2]
        particles =
            ParticleSet{1,T}((view(backing, 1:3),), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        config.temperature = view(backing, 4:5)
        config.extra = Set([Ref(T[0.3])])
        @test !HybridPlasmaPIC._particle_arrays_mightalias(config.temperature, particles.x[1])
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(io_ownership_state(st, (particles,)), st),
        )
        mktempdir() do dir
            path = joinpath(dir, "checkpoint")
            save_checkpoint(path, st, particles)
            target, destination, _ = io_ownership_run(State, T, Val(1), false)
            load_checkpoint!(target, destination..., path)
            @test io_ownership_config(target).temperature == T[0.1, 0.2]
            @test only(io_ownership_config(target).extra)[] == T[0.3]
        end
        # Container keys are part of configuration too, and can hold a state link.
        config.extra = IdDict(particles.x[1] => nothing)
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(st, (particles,)),
            st,
        )
    end
end

@testset "Independent mutable configuration retains cycles, aliases and continuation" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES, D = 1:3
        State === EMPIC1D && D != 1 && continue
        source, args, config = io_ownership_run(State, T, Val(D), variant)
        parameters = config isa IOOwnershipClosure ? config.temperature : config.parameters
        cycle = Any[]
        push!(cycle, cycle)
        config.extra = Dict{Symbol,Any}(
            :cycle => cycle,
            :first => parameters,
            :second => parameters,
            :nested => Ref(parameters),
            :simple => Core.svec(parameters, parameters),
            :label => IOOwnershipString("parameters", parameters),
            :geometry => (source.g.kfull, source.g.kvec, source.g.midx, source.g.ik),
        )
        if source isa Union{HybridStepper,CAMCLStepper}
            source.shape.parameters = parameters
        end
        io_ownership_initialize!(source, args)
        io_ownership_step!(source, args, T(0.001))
        if source isa Union{HybridStepper,CAMCLStepper}
            @test minimum(source.fields.n) > T(0.9)
        end
        target, destination, _ = io_ownership_run(State, T, Val(D), variant; count = 1)
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 7)
            archived = load_archive(archive).state
            io_ownership_compare_state(archived, deserialize(checkpoint))
            archived_config = io_ownership_config(archived)
            @test archived_config.extra[:first] === archived_config.extra[:second]
            @test archived_config.extra[:first] == parameters
            @test archived_config.extra[:simple][1] === archived_config.extra[:first]
            @test archived_config.extra[:label].parameters === archived_config.extra[:first]
            @test archived_config.extra[:cycle][1] === archived_config.extra[:cycle]
            load_checkpoint!(target, destination..., checkpoint)
            restored = io_ownership_config(target)
            @test restored.extra[:first] === restored.extra[:second] === restored.extra[:nested][]
            @test restored.extra[:first] !== parameters
            @test restored.extra[:first] == parameters
            @test restored.extra[:simple][1] ===
                  restored.extra[:simple][2] ===
                  restored.extra[:first]
            @test restored.extra[:label].parameters === restored.extra[:first]
            @test String(restored.extra[:label]) == "parameters"
            @test restored.extra[:cycle][1] === restored.extra[:cycle]
            if target isa Union{HybridStepper,CAMCLStepper}
                @test target.shape.parameters === restored.temperature
            end
            for dt in T[0.001, 0.0007]
                io_ownership_step!(source, args, dt)
                io_ownership_step!(target, destination, dt)
                if source isa Union{HybridStepper,CAMCLStepper}
                    @test minimum(source.fields.n) > T(0.9)
                end
                # Compare physical/phase state separately from cyclic configuration graphs.
                io_ownership_compare_state(
                    io_ownership_state(source, args),
                    io_ownership_state(target, destination),
                )
            end
        end
    end
end


@testset "Custom container metadata respects logical storage and retained ownership" begin
    for T in (Float32, Float64), State in (HybridStepper, CAMCLStepper, EMPIC, HallMHDState)
        source, args, config = io_ownership_run(State, T, Val(1), false)
        # Custom array getindex and Dict/Set interfaces use their stored offset.
        # Their primitive logical entries alone cannot expose this clock link.
        for wrapper in io_ownership_wrappers(T, source.time)
            config.extra = wrapper
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(source, args),
                source,
            )
        end
        config.extra = nothing
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 7)
            before_checkpoint, before_archive = read(checkpoint), read(archive)
            for wrapper in io_ownership_wrappers(T, source.time)
                config.extra = wrapper
                @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                @test read(checkpoint) == before_checkpoint
                @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 7)
                @test read(archive) == before_archive
            end
        end
        config.extra = nothing
        target, destination, _ = io_ownership_run(State, T, Val(1), false; count = 1)
        before = io_ownership_bytes(io_ownership_state(target, destination))
        original_config = io_ownership_config(target)
        saved = deepcopy(io_ownership_state(source, args))
        linked =
            haskey(saved, :x) ? saved.x[1] :
            haskey(saved, :electrons) ? saved.electrons.x[1] : saved.fields.n
        for wrapper in io_ownership_wrappers(T, Ref(zero(T)), Ref(linked))
            io_ownership_config(saved).extra = wrapper
            mktemp() do path, io
                close(io)
                serialize(path, saved)
                @test_throws ArgumentError load_checkpoint!(target, destination..., path)
                @test io_ownership_config(target) === original_config
                @test io_ownership_bytes(io_ownership_state(target, destination)) == before
            end
        end
    end
end

@testset "Wrapper parent metadata and logical visitation are order independent" begin
    for T in (Float32, Float64)
        source, args, config = io_ownership_run(HybridStepper, T, Val(1), false)
        ps = only(args)
        independent = T[0.1, 0.2, 0.3]
        wrapper = IOOwnershipArray(independent, Ref(zero(T)), nothing)
        view_of_wrapper = view(wrapper, 1:1)
        # LIFO traversal sees the view's parent as metadata before it later sees
        # the explicitly stored parent as logical data. The reverse order agrees.
        for roots in ((wrapper, view_of_wrapper), (view_of_wrapper, wrapper))
            visits = IdDict{Any,Int}()
            HybridPlasmaPIC._visit_checkpoint_objects(roots) do value
                visits[value] = get(visits, value, 0) + 1
            end
            @test get(visits, wrapper, 0) == 1
            @test get(visits, view_of_wrapper, 0) == 1
            @test !haskey(visits, independent)
        end
        for roots in ((ps.x[1], view(ps.x[1], 1:1)), (view(ps.x[1], 1:1), ps.x[1]))
            config.extra = roots
            @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(source, args),
                source,
            )
        end
        # Reaching only a disjoint view keeps the parent's other elements outside
        # the logical payload, while wrapper fields still participate in checks.
        payload = Any[ps.x[1], T[0.1]]
        config.extra = view(payload, 2:2)
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(io_ownership_state(source, args), source),
        )
        config.extra = (payload, view(payload, 2:2))
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(source, args),
            source,
        )
        backing = T[ps.x[1]..., 0.1, 0.2]
        particles =
            ParticleSet{1,T}((view(backing, 1:3),), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        wrapper = IOOwnershipArray(view(backing, 4:5), Ref(zero(T)), nothing)
        config.temperature = wrapper
        cycle = Any[]
        push!(cycle, cycle)
        offset = Ref(zero(T))
        wrapped = io_ownership_wrappers(T, offset, cycle)
        push!(cycle, wrapped)
        config.extra = wrapped
        visits = IdDict{Any,Int}()
        HybridPlasmaPIC._visit_checkpoint_objects(wrapped) do value
            visits[value] = get(visits, value, 0) + 1
        end
        @test all(==(1), values(visits))
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(source, (particles,)),
                source,
            ),
        )
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, particles)
            archive_run(archive, source, particles; rng_seed = 7)
            target, destination, _ = io_ownership_run(HybridStepper, T, Val(1), false)
            load_checkpoint!(target, destination..., checkpoint)
            restored = io_ownership_config(target)
            @test collect(restored.temperature) == T[0.1, 0.2]
            @test restored.extra[1].offset === restored.extra[4].offset === restored.extra[5].offset
            @test restored.extra[1].extra[1] === restored.extra[1].extra
            @test restored.extra[1].extra[2][1].offset === restored.extra[1].offset
            archived = io_ownership_config(load_archive(archive).state)
            @test archived.extra[1].offset === archived.extra[4].offset === archived.extra[5].offset
        end
        # A clock link in the same parent wrapper is detected even through a view.
        wrapper = IOOwnershipArray(view(backing, 4:5), source.time, nothing)
        config.temperature = view(wrapper, 1:1)
        config.extra = nothing
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(source, (particles,)),
            source,
        )
    end
end

@testset "Set dictionary ownership follows stored identity without exposing capacity" begin
    for T in (Float32, Float64), (State, variant) in IO_OWNERSHIP_FAMILIES
        State === HallMHDState && continue
        source, original_args, config = io_ownership_run(State, T, Val(1), variant)
        ps = first(io_ownership_particles(original_args))
        x = IOOwnershipTrackedArray(copy(ps.x[1]), Set(ps.x[1]))
        particles = ParticleSet{1,T}((x,), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        args = if first(original_args) isa AbstractVector
            ([particles; first(original_args)[2:end]],)
        else
            (particles, Base.tail(original_args)...)
        end
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 7)
            before_checkpoint, before_archive = read(checkpoint), read(archive)
            config.extra = x.members.dict
            @test config.extra === particles.x[1].members.dict
            particles.x[1][1] = T(0.25)
            @test haskey(config.extra, T(0.25))
            @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
            @test read(checkpoint) == before_checkpoint
            @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 7)
            @test read(archive) == before_archive

            # Generic Set serialization retains the dictionary link in a payload.
            saved = io_ownership_state(source, args)
            serialize(checkpoint, saved)
            incoming = deserialize(checkpoint)
            record =
                haskey(incoming, :species) ? first(incoming.species) :
                haskey(incoming, :electrons) ? incoming.electrons : incoming
            @test io_ownership_config(incoming).extra === record.x[1].members.dict
            target, destination, original_config = io_ownership_run(State, T, Val(1), variant)
            before = io_ownership_bytes(io_ownership_state(target, destination))
            @test_throws ArgumentError load_checkpoint!(target, destination..., checkpoint)
            @test io_ownership_config(target) === original_config
            @test io_ownership_bytes(io_ownership_state(target, destination)) == before
        end
    end
    for T in (Float32, Float64)
        source, args, config = io_ownership_run(HybridStepper, T, Val(1), false)
        members = Set(T[0.1, 0.2])
        config.extra = (members = members, dictionary = members.dict)
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 7)
            target, destination, _ = io_ownership_run(HybridStepper, T, Val(1), false)
            load_checkpoint!(target, destination..., checkpoint)
            for restored in
                (io_ownership_config(target), io_ownership_config(load_archive(archive).state))
                @test restored.extra.members.dict === restored.extra.dictionary
                push!(restored.extra.members, T(0.3))
                @test haskey(restored.extra.dictionary, T(0.3))
            end
        end
    end
end

@testset "Immutable metadata equality does not imply evolving storage ownership" begin
    for T in (Float32, Float64), State in (HybridStepper, CAMCLStepper, EMPIC1D, EMPIC)
        source, original_args, config = io_ownership_run(State, T, Val(1), false)
        ps = first(original_args)
        metadata = (label = "ion",)
        offset = Ref(zero(T))
        x = IOOwnershipArray(copy(ps.x[1]), offset, metadata)
        particles = ParticleSet{1,T}((x,), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        args = (particles, Base.tail(original_args)...)
        config.extra = (label = String(UInt8[0x69, 0x6f, 0x6e]),)
        @test config.extra === metadata
        @test !ismutable(metadata)
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 7)
            target, destination, _ = io_ownership_run(State, T, Val(1), false)
            load_checkpoint!(target, destination..., checkpoint)
            @test io_ownership_config(target).extra == metadata
            @test io_ownership_config(load_archive(archive).state).extra == metadata
            # Immutable parents still expose mutable children that belong to state.
            for child in (offset, particles.x[1], source.g.iplan)
                config.extra = (metadata = metadata, nested = (child = child,))
                @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 7)
            end
        end
    end
end

@testset "Packed backing storage does not expand logical bit views" begin
    for T in (Float32, Float64), State in (HybridStepper, CAMCLStepper, EMPIC1D, EMPIC)
        source, original_args, config = io_ownership_run(State, T, Val(1), false)
        ps = first(original_args)
        bits = trues(6)
        owned_view, independent_view = view(bits, 1:3), view(bits, 4:6)
        @test !Base.mightalias(owned_view, independent_view)
        independent_view[1] = false
        @test all(owned_view)
        x = IOOwnershipArray(copy(ps.x[1]), Ref(zero(T)), owned_view)
        particles = ParticleSet{1,T}((x,), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        args = (particles, Base.tail(original_args)...)
        config.extra = independent_view
        mktempdir() do dir
            checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 7)
            target, destination, original_config = io_ownership_run(State, T, Val(1), false)
            load_checkpoint!(target, destination..., checkpoint)
            @test io_ownership_config(target).extra == independent_view
            @test io_ownership_config(load_archive(archive).state).extra == independent_view
            before_checkpoint, before_archive = read(checkpoint), read(archive)
            before = io_ownership_bytes(io_ownership_state(target, destination))
            original_config = io_ownership_config(target)
            for linked in (view(bits, 2:4), bits, bits.chunks, view(bits.chunks, 1:1))
                config.extra = linked
                @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                @test read(checkpoint) == before_checkpoint
                @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 7)
                @test read(archive) == before_archive
            end
            # Array-backed chunk views are copied by Serialization; only these
            # three relationships remain available for incoming ownership checks.
            for linked in (view(bits, 2:4), bits, bits.chunks)
                config.extra = linked
                incoming = joinpath(dir, "incoming")
                serialize(incoming, io_ownership_state(source, args))
                saved = deserialize(incoming)
                record = haskey(saved, :electrons) ? saved.electrons : saved
                saved_bits = parent(record.x[1].extra)
                saved_extra = io_ownership_config(saved).extra
                if linked isa SubArray
                    @test parent(saved_extra) === saved_bits
                elseif linked isa BitArray
                    @test saved_extra === saved_bits
                else
                    @test saved_extra === saved_bits.chunks
                end
                @test_throws ArgumentError load_checkpoint!(target, destination..., incoming)
                @test io_ownership_config(target) === original_config
                @test io_ownership_bytes(io_ownership_state(target, destination)) == before
            end
            config.extra = copy(bits)
            save_checkpoint(checkpoint, source, args...)
            load_checkpoint!(target, destination..., checkpoint)
            @test io_ownership_config(target).extra == bits
        end
        # Whole packed arrays also own their backing buffer.
        x = IOOwnershipArray(copy(ps.x[1]), Ref(zero(T)), bits)
        particles = ParticleSet{1,T}((x,), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        config.extra = bits.chunks
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_state(
            io_ownership_state(source, (particles, Base.tail(original_args)...)),
            source,
        )
    end
    bits = trues(6)
    slice = view(bits, 1:3)
    for roots in ((bits.chunks, slice), (slice, bits.chunks))
        visits, backing = IdDict{Any,Int}(), IdDict{Any,Int}()
        remember_backing(value) = (backing[value] = get(backing, value, 0) + 1)
        HybridPlasmaPIC._visit_checkpoint_objects(roots; backing = remember_backing) do value
            visits[value] = get(visits, value, 0) + 1
        end
        @test get(visits, slice, 0) == 1
        @test get(visits, bits.chunks, 0) == 1
        @test !haskey(visits, bits)
        @test get(backing, bits.chunks, 0) == 1
    end
    for (roots, backing_visits) in (((bits, slice), 2), ((slice, bits), 1))
        visits, backing = IdDict{Any,Int}(), IdDict{Any,Int}()
        remember_backing(value) = (backing[value] = get(backing, value, 0) + 1)
        HybridPlasmaPIC._visit_checkpoint_objects(roots; backing = remember_backing) do value
            visits[value] = get(visits, value, 0) + 1
        end
        @test get(visits, bits, 0) == 1
        @test get(visits, slice, 0) == 1
        @test !haskey(visits, bits.chunks)
        @test get(backing, bits.chunks, 0) == backing_visits
    end
end


@testset "Dictionary backing ownership keeps unused capacity opaque" begin
    for T in (Float32, Float64), State in (HybridStepper, CAMCLStepper, EMPIC1D, EMPIC)
        source, original_args, config = io_ownership_run(State, T, Val(1), false)
        ps = first(original_args)
        tables = (
            Dict{Symbol,Any}(:parameter => Ref(T(0.1))),
            IdDict{Any,Any}(:parameter => Ref(T(0.1))),
            Set(T[0.1, 0.2]),
        )
        for table in tables
            x = IOOwnershipArray(copy(ps.x[1]), Ref(zero(T)), table)
            particles = ParticleSet{1,T}((x,), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
            args = (particles, Base.tail(original_args)...)
            dictionary = table isa Set ? table.dict : table
            buffers =
                dictionary isa Dict ? (dictionary.slots, dictionary.keys, dictionary.vals) :
                (dictionary.ht,)
            config.extra = nothing
            mktempdir() do dir
                checkpoint, archive = joinpath(dir, "checkpoint"), joinpath(dir, "archive")
                save_checkpoint(checkpoint, source, args...)
                archive_run(archive, source, args...; rng_seed = 7)
                before_checkpoint, before_archive = read(checkpoint), read(archive)
                for buffer in buffers, linked in (buffer, view(buffer, :))
                    config.extra = linked
                    @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                    @test read(checkpoint) == before_checkpoint
                    @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 7)
                    @test read(archive) == before_archive
                end
                independent = T[0.1, 0.2, 0.3, 0.4]
                parameters = view(independent, 3:4)
                config.extra = Dict(:parameters => parameters, :other => view(independent, 1:2))
                @test !Base.mightalias(config.extra[:parameters], config.extra[:other])
                save_checkpoint(checkpoint, source, args...)
                archive_run(archive, source, args...; rng_seed = 7)
                target, destination, _ = io_ownership_run(State, T, Val(1), false)
                load_checkpoint!(target, destination..., checkpoint)
                @test io_ownership_config(target).extra[:parameters] == T[0.3, 0.4]
                @test io_ownership_config(load_archive(archive).state).extra[:other] == T[0.1, 0.2]
            end
        end
        # An unused private table slot is outside the logical dictionary payload.
        # Register the allocation without following stale or unassigned elements.
        parameter = T[0.1]
        dictionary = Dict{Symbol,Any}(:old => parameter)
        empty!(dictionary)
        dictionary.vals[1] = parameter
        @test isempty(dictionary)
        @test isempty(collect(dictionary))
        x = IOOwnershipArray(copy(ps.x[1]), Ref(zero(T)), dictionary)
        particles = ParticleSet{1,T}((x,), ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
        config.extra = parameter
        @test isnothing(
            HybridPlasmaPIC._validate_checkpoint_state(
                io_ownership_state(source, (particles, Base.tail(original_args)...)),
                source,
            ),
        )
    end
    for dictionary in (Dict{Symbol,Any}(:value => Ref(1)), IdDict{Any,Any}(:value => Ref(1)))
        buffer = dictionary isa Dict ? dictionary.vals : dictionary.ht
        slice = view(buffer, :)
        for roots in ((buffer, slice, dictionary), (dictionary, slice, buffer))
            visits, backing = IdDict{Any,Int}(), IdDict{Any,Int}()
            remember_backing(value) = (backing[value] = get(backing, value, 0) + 1)
            HybridPlasmaPIC._visit_checkpoint_objects(roots; backing = remember_backing) do value
                visits[value] = get(visits, value, 0) + 1
            end
            @test get(visits, dictionary, 0) == 1
            @test get(visits, buffer, 0) == 1
            @test get(visits, slice, 0) == 1
            @test get(backing, buffer, 0) == 1
        end
    end
end
