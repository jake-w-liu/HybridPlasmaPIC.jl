using HybridPlasmaPIC, Test, Serialization

function io_clock_run(State, ::Type{T}, variant; target = false) where {T}
    g = FourierGrid((3,), (T(2π),))
    closure = IsothermalElectrons(T(target ? 0.2 : 0.1))
    if State === HallMHDState
        st = State(g, HallMHDModel(closure))
        fill!(st.fields.n, one(T))
        fill!(st.fields.B[1], one(T))
        return st, ()
    end
    count = target ? 1 : 3
    hybrid = State in (HybridStepper, CAMCLStepper)
    function population(q, m, density)
        ps = ParticleSet{1,T}(count; q, m)
        load_lattice!(ps, (zero(T),), g.L, (count,))
        set_density_weight!(ps, T(density), g)
        foreach(a -> fill!(a, T(0.001)), ps.v)
        return ps
    end
    if hybrid
        st = State(g, HybridModel(closure), CIC(), count)
        fill!(st.fields.B[1], one(T))
        first = population(one(T), one(T), variant ? 0.5 : 1)
        particles = variant ? [first, population(T(2), T(4), 0.25)] : first
        State === HybridStepper ? init!(st, particles) : init_camcl!(st, particles)
        return st, (particles,)
    end
    st = State(g, count; mobile = variant, mi = T(16), c = T(target ? 4 : 3), n_sub = 2)
    electrons = population(-one(T), one(T), 1)
    ions = variant ? population(one(T), T(16), 1) : nothing
    init_empic!(st, electrons, ions)
    return st, (electrons, ions)
end

io_clock_step!(st::HybridStepper, args, dt) = step!(st, args..., dt; NB = 2)
io_clock_step!(st::CAMCLStepper, args, dt) = step_camcl!(st, args..., dt; NB = 2)
io_clock_step!(st::Union{EMPIC1D,EMPIC}, args, dt) = step_empic!(st, args..., dt)
io_clock_step!(st::HallMHDState, args, dt) = step_hall_mhd!(st, dt)

io_clock_buffers(st::Union{HybridStepper,CAMCLStepper}) =
    (st.Ep..., st.Bp..., st.vpred..., st.xmid..., st.work)
io_clock_buffers(st::EMPIC1D) = (st.Exp, st.Eyp, st.Bzp, st.work, st.Exi, st.Eyi, st.Bzi, st.worki)
io_clock_buffers(st::EMPIC) =
    (st.Ep..., st.Bp..., st.mide..., st.worke, st.Epi..., st.Bpi..., st.midi..., st.worki)
io_clock_buffers(st::HallMHDState) = ()

function io_clock_bytes(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

io_clock_state(st, args) = HybridPlasmaPIC._checkpoint_state(st, args...)
io_clock_snapshot(st, args) = io_clock_bytes((io_clock_state(st, args), io_clock_buffers(st)))

const IO_CLOCK_FAMILIES = (
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

@testset "Invalid checkpoint clocks reject before writing or restoring" begin
    for T in (Float32, Float64), (State, variant) in IO_CLOCK_FAMILIES, live in (false, true)
        source, args = io_clock_run(State, T, variant)
        live && io_clock_step!(source, args, T(0.001))
        target, destination = io_clock_run(State, T, variant; target = true)
        before = io_clock_snapshot(target, destination)
        buffers = io_clock_buffers(target)
        mktempdir() do dir
            checkpoint = joinpath(dir, "checkpoint.bin")
            archive = joinpath(dir, "archive.bin")
            malformed = joinpath(dir, "malformed.bin")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 9)
            checkpoint_bytes, archive_bytes = read(checkpoint), read(archive)
            valid = deserialize(checkpoint)
            for time in T[NaN, Inf, -Inf, -1]
                source.time[] = time
                @test_throws ArgumentError save_checkpoint(checkpoint, source, args...)
                @test read(checkpoint) == checkpoint_bytes
                @test_throws ArgumentError archive_run(archive, source, args...; rng_seed = 9)
                @test read(archive) == archive_bytes
                serialize(malformed, merge(valid, (; time)))
                @test_throws ArgumentError load_checkpoint!(target, destination..., malformed)
                @test io_clock_snapshot(target, destination) == before
                @test all(a === b for (a, b) in zip(buffers, io_clock_buffers(target)))
            end
        end
    end
end

@testset "Checkpoint clock origins retain exact continuation and archive parity" begin
    for T in (Float32, Float64),
        (State, variant) in IO_CLOCK_FAMILIES,
        live in (false, true),
        origin in (T(-0.0), T(2))

        source, args = io_clock_run(State, T, variant)
        source.time[] = origin
        live && io_clock_step!(source, args, T(0.001))
        target, destination = io_clock_run(State, T, variant; target = true)
        mktempdir() do dir
            checkpoint = joinpath(dir, "checkpoint.bin")
            archive = joinpath(dir, "archive.bin")
            save_checkpoint(checkpoint, source, args...)
            archive_run(archive, source, args...; rng_seed = 9)
            @test io_clock_bytes(load_archive(archive).state) == read(checkpoint)
            @test load_checkpoint!(target, destination..., checkpoint) === target
            @test isequal(target.time[], source.time[])
            @test io_clock_bytes(io_clock_state(target, destination)) ==
                  io_clock_bytes(io_clock_state(source, args))
            for dt in T[0.001, 0.0007]
                io_clock_step!(source, args, dt)
                io_clock_step!(target, destination, dt)
                @test target.time[] === source.time[]
                @test io_clock_bytes(io_clock_state(target, destination)) ==
                      io_clock_bytes(io_clock_state(source, args))
            end
        end
    end
end

@testset "Checkpoint clock boundary contract" begin
    for T in (Float32, Float64)
        # A saved terminal state need not be able to advance its counter.
        @test HybridPlasmaPIC._validate_checkpoint_clock_state(
            (time = floatmax(T), step = typemax(Int)),
            T,
        ) === nothing
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_clock_state(
            (time = one(T), step = -1),
            T,
        )
        @test_throws ArgumentError HybridPlasmaPIC._validate_checkpoint_clock_state(
            (time = 1, step = 0),
            T,
        )
    end
end
