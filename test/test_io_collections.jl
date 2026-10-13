using HybridPlasmaPIC, Test, Serialization

io_multi_init!(st::HybridStepper, species) = init!(st, species)
io_multi_init!(st::CAMCLStepper, species) = init_camcl!(st, species)
io_multi_step!(st::HybridStepper, species, dt; NB = 2) = step!(st, species, dt; NB)
io_multi_step!(st::CAMCLStepper, species, dt; NB = 2) = step_camcl!(st, species, dt; NB)

function io_multi_run(
    State,
    ::Type{T},
    ::Val{D};
    shape = CIC(),
    closure = IsothermalElectrons(T(0.1)),
) where {T,D}
    g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> T(2π), D))
    species = ParticleSet{D,T}[]
    for s = 1:3
        # Multiple particles per cell keep NGP occupied as particles cross assignment seams.
        counts = ntuple(_ -> s == 3 ? 0 : 3(s + 1), D)
        ps = ParticleSet{D,T}(prod(counts); q = T(s), m = T(s^2))
        if nparticles(ps) > 0
            load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
            set_density_weight!(ps, T(0.4 / s), g)
            for c = 1:3
                ps.v[c] .= T(0.002c * s) .* cos.(ps.x[1])
            end
            ps.id .= UInt64(1000s) .+ UInt64.(1:nparticles(ps))
            ps.tag .= UInt32(s)
        end
        push!(species, ps)
    end
    st = State(g, HybridModel(closure; η = 0.001, ηH = 0.0001, de2 = 0.01), shape, 0)
    fill!(st.fields.B[1], one(T))
    for I in CartesianIndices(g.n)
        st.fields.B[2][I] = T(0.005) * cos(T(I[1] - 1) * g.dx[1])
    end
    io_multi_init!(st, species)
    return st, species
end

# Read persistent solver and population state independently of the I/O helpers.
function io_multi_bytes(st, species)
    cl = st isa CAMCLStepper ? (st.Blead, st.Blag, st.Bout, st.cl_h[], st.cl_live[]) : ()
    populations = [
        (
            map(collect, ps.x),
            map(collect, ps.v),
            collect(ps.weight),
            collect(ps.id),
            collect(ps.tag),
            ps.q,
            ps.m,
        ) for ps in species
    ]
    state = (
        st.model,
        st.shape,
        st.fields.B,
        st.fields.E,
        st.time[],
        st.step[],
        st.last_dt[],
        cl,
        populations,
    )
    io = IOBuffer()
    serialize(io, state)
    return take!(io)
end

io_multi_buffers(st) = (st.Ep..., st.Bp..., st.vpred..., st.xmid..., st.work)
io_multi_arrays(ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)

function io_multi_reject_unchanged(st, species, path)
    before = io_multi_bytes(st, species)
    buffers = io_multi_buffers(st)
    arrays = map(io_multi_arrays, species)
    @test_throws ArgumentError load_checkpoint!(st, species, path)
    @test io_multi_bytes(st, species) == before
    @test all(a === b for (a, b) in zip(buffers, io_multi_buffers(st)))
    @test all(
        all(a === b for (a, b) in zip(old, io_multi_arrays(ps))) for
        (old, ps) in zip(arrays, species)
    )
end

function io_multi_views(ps::ParticleSet{D,T}, range = 1:nparticles(ps)) where {D,T}
    return ParticleSet{D,T}(
        map(a -> view(a, range), ps.x),
        map(a -> view(a, range), ps.v),
        view(ps.weight, range),
        view(ps.id, range),
        view(ps.tag, range),
        ps.q,
        ps.m,
    )
end

@testset "Hybrid collection checkpoints preserve ordered species and integrator phase" begin
    for T in (Float32, Float64),
        D = 1:3,
        shape in (NGP(), CIC(), TSC()),
        State in (HybridStepper, CAMCLStepper),
        phase in (:initial, :live, :reset, :edited)

        st, species = io_multi_run(State, T, Val(D); shape)
        if phase !== :initial
            io_multi_step!(st, species, T(0.002); NB = 2)
            io_multi_step!(st, species, T(0.001); NB = 3)
        end
        phase === :reset && io_multi_init!(st, species)
        phase === :edited && (st.fields.B[2][1] += T(0.001))
        mktemp() do path, io
            close(io)
            @test save_checkpoint(path, st, species) == path
            target = State(
                FourierGrid(st.g.n, st.g.L),
                HybridModel(IsothermalElectrons(T(0.2))),
                shape,
                1,
            )
            destinations = [ParticleSet{D,T}(1; q = T(4), m = T(7)) for _ in species]
            objects = copy(destinations)
            @test load_checkpoint!(target, destinations, path) === target
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
            @test all(a === b for (a, b) in zip(objects, destinations))
            @test all(length(a) == sum(nparticles, species) for a in io_multi_buffers(target))
            for (dt, NB) in ((T(0.0015), 4), (T(0.0008), 2), (T(0.0008), 3))
                io_multi_step!(st, species, dt; NB)
                io_multi_step!(target, destinations, dt; NB)
                @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
            end
        end
    end
end

@testset "Collection restart accepts independent views and replaces placeholder values" begin
    for State in (HybridStepper, CAMCLStepper), T in (Float32, Float64), D = 1:3
        st, species = io_multi_run(State, T, Val(D))
        io_multi_step!(st, species, T(0.001))
        target, backing = io_multi_run(State, T, Val(D))
        destinations = [io_multi_views(ps) for ps in backing]
        arrays = map(io_multi_arrays, destinations)
        for ps in destinations
            ps.q = T(NaN)
            ps.m = zero(T)
            foreach(a -> fill!(a, T(NaN)), (ps.x..., ps.v..., ps.weight))
        end
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, [io_multi_views(ps) for ps in species])
            load_checkpoint!(target, destinations, path)
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
            @test all(
                all(a === b for (a, b) in zip(old, io_multi_arrays(ps))) for
                (old, ps) in zip(arrays, destinations)
            )
            io_multi_step!(st, species, T(0.002); NB = 3)
            io_multi_step!(target, destinations, T(0.002); NB = 3)
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
        end
    end
end

@testset "Collection checkpoint rejects later invalid populations before mutation" begin
    for State in (HybridStepper, CAMCLStepper)
        st, species = io_multi_run(State, Float64, Val(1))
        io_multi_step!(st, species, 0.001)
        target, destinations = io_multi_run(State, Float64, Val(1))
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, species)
            saved = deserialize(path)
            record = saved.species[2]
            bad_records = (
                merge(record, (q = -1.0,)),
                merge(record, (q = NaN,)),
                merge(record, (m = 0.0,)),
                merge(record, (m = 2.0f0,)),
                merge(record, (q = floatmax(Float64), m = eps(Float64))),
                merge(record, (x = (fill(Inf, length(record.weight)),),)),
                merge(record, (v = (record.v[1], record.v[2], fill(NaN, length(record.weight))),)),
                merge(record, (weight = fill(-1.0, length(record.weight)),)),
                merge(record, (id = UInt64[1],)),
                merge(record, (v = (record.v[1], record.v[1], record.v[3]),)),
                saved.species[1],
                nothing,
            )
            for record in bad_records
                serialize(
                    path,
                    merge(saved, (species = Any[saved.species[1], record, saved.species[3]],)),
                )
                io_multi_reject_unchanged(target, destinations, path)
            end
            for bad in (
                merge(saved, (particle_layout = :unknown,)),
                merge(saved, (species = (),)),
                merge(saved, (last_dt = 0.0,)),
                merge(saved, (integrator = :unsupported,)),
                (; (k => v for (k, v) in pairs(saved) if k != :species)...),
                (; (k => v for (k, v) in pairs(saved) if k != :particle_layout)...),
            )
                serialize(path, bad)
                io_multi_reject_unchanged(target, destinations, path)
            end
            save_checkpoint(path, st, species)
            bytes = read(path)
            species[2].q = -1.0
            @test_throws ArgumentError save_checkpoint(path, st, species)
            @test read(path) == bytes
            @test_throws ArgumentError archive_run(path, st, species; rng_seed = 5)
            @test read(path) == bytes
        end
    end
end

@testset "Collection checkpoint preflights count, layout, storage and configuration" begin
    for State in (HybridStepper, CAMCLStepper)
        st, species = io_multi_run(State, Float64, Val(1))
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, species)
            target, destinations = io_multi_run(State, Float64, Val(1))
            target.model = HybridModel(IsothermalElectrons(0.7))
            fill!(target.fields.B[1], 2.0)
            destinations[1].x[1] .+= 0.1
            destinations[1].v[1] .+= 0.2
            destinations[1].q = 3.0
            io_multi_reject_unchanged(target, destinations[1:2], path)
            io_multi_reject_unchanged(
                target,
                [destinations[1], destinations[1], destinations[3]],
                path,
            )
            empty = ParticleSet{1,Float64}(0)
            io_multi_reject_unchanged(target, [empty, empty, destinations[3]], path)
            other = ParticleSet{1,Float64}(0)
            shared_empty = ParticleSet{1,Float64}(
                empty.x,
                other.v,
                other.weight,
                other.id,
                other.tag,
                other.q,
                other.m,
            )
            io_multi_reject_unchanged(target, [empty, shared_empty, destinations[3]], path)
            within_empty = ParticleSet{1,Float64}(
                empty.x,
                (empty.x[1], empty.v[2], empty.v[3]),
                empty.weight,
                empty.id,
                empty.tag,
                empty.q,
                empty.m,
            )
            io_multi_reject_unchanged(
                target,
                [destinations[1], within_empty, destinations[3]],
                path,
            )
            late_view = [destinations[1], io_multi_views(destinations[2], 1:1), destinations[3]]
            io_multi_reject_unchanged(target, late_view, path)
            wrong_dimension = [ParticleSet{2,Float64}(1) for _ in species]
            io_multi_reject_unchanged(target, wrong_dimension, path)
            wrong_precision = [ParticleSet{1,Float32}(1) for _ in species]
            io_multi_reject_unchanged(target, wrong_precision, path)
            for incompatible in (
                State(FourierGrid((4,), st.g.L), st.model, CIC(), 1),
                State(FourierGrid(st.g.n, (4π,)), st.model, CIC(), 1),
                State(st.g, HybridModel(PolytropicElectrons(0.1, 1.0, 5 / 3)), CIC(), 1),
                State(st.g, st.model, TSC(), 1),
                (State === HybridStepper ? CAMCLStepper : HybridStepper)(st.g, st.model, CIC(), 1),
            )
                io_multi_reject_unchanged(incompatible, destinations, path)
            end
            # Even a one-population collection has a distinct checkpoint layout.
            save_checkpoint(path, st, species[1:1])
            before = io_multi_bytes(target, destinations)
            @test_throws ArgumentError load_checkpoint!(target, destinations[1], path)
            @test io_multi_bytes(target, destinations) == before
            save_checkpoint(path, st, species[1])
            io_multi_reject_unchanged(target, destinations[1:1], path)
        end
    end
end

@testset "Typed empty collections preserve shared fields and phase on restart" begin
    for State in (HybridStepper, CAMCLStepper), T in (Float32, Float64), D = 1:3
        st, _ = io_multi_run(State, T, Val(D))
        species = ParticleSet{D,T}[]
        fill!(st.fields.B[2], T(0.005))
        io_multi_init!(st, species)
        io_multi_step!(st, species, T(0.001))
        target = State(FourierGrid(st.g.n, st.g.L), st.model, st.shape, 1)
        destinations = ParticleSet{D,T}[]
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, species)
            load_checkpoint!(target, destinations, path)
            @test isempty(destinations)
            @test all(isempty, io_multi_buffers(target))
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
            io_multi_step!(st, species, T(0.002); NB = 3)
            io_multi_step!(target, destinations, T(0.002); NB = 3)
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
        end
    end
end

struct IOMultiArrayClosure <: HybridPlasmaPIC.ElectronClosure
    temperature::Vector{Float64}
    alias::Vector{Float64}
end
function HybridPlasmaPIC.electron_pressure!(pe, n, closure::IOMultiArrayClosure)
    @. pe = closure.temperature[1] * n
    return pe
end

@testset "Collection checkpoints own saved custom configuration" begin
    for State in (HybridStepper, CAMCLStepper)
        values = [0.1]
        st, species =
            io_multi_run(State, Float64, Val(1); closure = IOMultiArrayClosure(values, values))
        expected = io_multi_bytes(st, species)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, species)
            values[1] = 0.7
            replacement = [0.2]
            target, destinations = io_multi_run(
                State,
                Float64,
                Val(1);
                closure = IOMultiArrayClosure(replacement, replacement),
            )
            load_checkpoint!(target, destinations, path)
            @test io_multi_bytes(target, destinations) == expected
            @test target.model.closure.temperature === target.model.closure.alias
            @test target.model.closure.temperature !== values
            values[1] = 0.1
            io_multi_step!(st, species, 0.001)
            io_multi_step!(target, destinations, 0.001)
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
        end
    end
end

@testset "Hybrid collection archives preserve canonical species and shared phase" begin
    for State in (HybridStepper, CAMCLStepper), T in (Float32, Float64), D = 1:3
        st, species = io_multi_run(State, T, Val(D))
        io_multi_step!(st, species, T(0.001))
        mktempdir() do dir
            checkpoint = joinpath(dir, "checkpoint.bin")
            archive = joinpath(dir, "archive.bin")
            save_checkpoint(checkpoint, st, species)
            @test archive_run(archive, st, species; rng_seed = 21) == archive
            loaded = load_archive(archive)
            @test loaded.meta.rng_seed == 21
            @test loaded.state == deserialize(checkpoint)
            target, destinations = io_multi_run(State, T, Val(D))
            HybridPlasmaPIC._restore_checkpoint_state!(target, destinations, loaded.state)
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
            io_multi_step!(st, species, T(0.002); NB = 3)
            io_multi_step!(target, destinations, T(0.002); NB = 3)
            @test io_multi_bytes(target, destinations) == io_multi_bytes(st, species)
            packed = deserialize(archive)
            corrupt = copy(packed.state)
            corrupt[end] ⊻= 0x01
            serialize(archive, merge(packed, (state = corrupt,)))
            @test_throws ErrorException load_archive(archive)
        end
    end
end
