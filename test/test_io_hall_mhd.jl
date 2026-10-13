using HybridPlasmaPIC, Test, Serialization

function io_hall_run(::Type{T}, ::Val{D}, closure; ηH = 0.003) where {T,D}
    g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> T(2π), D))
    model = HallMHDModel(closure; Ti = 0.07, η = 0.002, ηH, nfloor = 1e-5)
    st = HallMHDState(g, model)
    for I in CartesianIndices(g.n)
        x = T(I[1] - 1) * g.dx[1]
        st.fields.n[I] = one(T) + T(0.02) * cos(x)
        for c = 1:3
            st.fields.ui[c][I] = T(0.004c) * sin(x)
            st.fields.B[c][I] = T(c == 1) + T(0.003c) * cos(x)
        end
    end
    return st
end

function io_hall_closures(::Type{T}) where {T}
    return (
        IsothermalElectrons(T(0.12)),
        PolytropicElectrons(T(0.12), one(T), T(5 / 3)),
        CGLElectrons(T(0.10), T(0.12), one(T), one(T)),
    )
end

# Inspect every public field directly, independently of checkpoint capture.
function io_hall_bytes(st)
    values = map(fieldnames(typeof(st.fields))) do name
        value = getfield(st.fields, name)
        value isa Ref ? value[] : value
    end
    io = IOBuffer()
    serialize(io, (st.model, values, st.time[], st.step[]))
    return take!(io)
end

io_hall_references(st) = map(name -> getfield(st.fields, name), fieldnames(typeof(st.fields)))

function io_hall_same_references(st, references)
    return all(
        getfield(st.fields, name) === reference for
        (name, reference) in zip(fieldnames(typeof(st.fields)), references)
    )
end

function io_hall_reject_unchanged(st, path)
    before = io_hall_bytes(st)
    references = io_hall_references(st)
    model = st.model
    stage_lapJ = deepcopy(st.stage.lapJ)
    @test_throws ArgumentError load_checkpoint!(st, path)
    @test io_hall_bytes(st) == before
    @test st.model === model
    @test io_hall_same_references(st, references)
    @test st.stage.lapJ == stage_lapJ
end

@testset "Hall-MHD checkpoint retains inactive hyperresistive stage history" begin
    for T in (Float32, Float64), D = 1:3, closure in io_hall_closures(T)
        st = io_hall_run(T, Val(D), closure; ηH = 0.02)
        step_hall_mhd!(st, T(0.001))
        st.model = HallMHDModel(closure; Ti = st.model.Ti, η = st.model.η, ηH = 0.0)
        for a in st.fields.lapJ
            fill!(a, T(0.125))
        end
        @test st.fields.lapJ != st.stage.lapJ
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st)
            target = io_hall_run(T, Val(D), closure; ηH = 0.0)
            stage_references = target.stage.lapJ
            load_checkpoint!(target, path)
            @test io_hall_bytes(target) == io_hall_bytes(st)
            @test target.stage.lapJ == st.stage.lapJ
            @test target.stage.lapJ === stage_references
            for dt in (T(0.002), T(0.001))
                step_hall_mhd!(st, dt)
                step_hall_mhd!(target, dt)
                @test io_hall_bytes(target) == io_hall_bytes(st)
            end
        end
    end
end

@testset "Hall-MHD checkpoint preserves fluid, derived fields, model and clock" begin
    for T in (Float32, Float64),
        D = 1:3,
        closure in io_hall_closures(T),
        phase in (:initial, :ohm, :live),
        ηH in (0.0, 0.003)

        st = io_hall_run(T, Val(D), closure; ηH)
        phase === :initial || hall_mhd_ohms_law!(st)
        if phase === :live
            step_hall_mhd!(st, T(0.001); project_b = false)
            step_hall_mhd!(st, T(0.002))
        end
        mktemp() do path, io
            close(io)
            @test save_checkpoint(path, st) == path
            target = HallMHDState(
                FourierGrid(st.g.n, st.g.L),
                HallMHDModel(closure; Ti = 0.2, η = 0.0, ηH = 0.01, nfloor = 1e-4),
            )
            references = io_hall_references(target)
            @test load_checkpoint!(target, path) === target
            @test io_hall_bytes(target) == io_hall_bytes(st)
            @test io_hall_same_references(target, references)
            @test magnetic_energy(target.fields.B, target.g) === magnetic_energy(st.fields.B, st.g)
            @test jdotE_density(target.fields.J, target.fields.E) ==
                  jdotE_density(st.fields.J, st.fields.E)
            for (dt, project_b) in ((T(0.0015), false), (T(0.0007), true), (T(0.0007), false))
                step_hall_mhd!(st, dt; project_b)
                step_hall_mhd!(target, dt; project_b)
                @test io_hall_bytes(target) == io_hall_bytes(st)
            end
        end
    end
end

@testset "Hall-MHD checkpoint preserves unevaluated and externally edited fields" begin
    for T in (Float32, Float64), D = 1:3, closure in io_hall_closures(T)
        g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> T(2π), D))
        st = HallMHDState(g, HallMHDModel(closure))
        # Constructor density is zero: preserving an initialization snapshot
        # must not evaluate Ohm's law or demand an already-steppable density.
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st)
            target = io_hall_run(T, Val(D), closure)
            load_checkpoint!(target, path)
            @test io_hall_bytes(target) == io_hall_bytes(st)
            for name in fieldnames(typeof(st.fields))
                value = getfield(st.fields, name)
                if value isa Ref
                    value[] = length(st.fields.n)
                else
                    for a in (value isa Tuple ? value : (value,))
                        fill!(a, T(0.125))
                    end
                end
            end
            save_checkpoint(path, st)
            saved = deserialize(path)
            # Named fields are matched by name, including reordered containers.
            reordered =
                (; (k => getproperty(saved.fields, k) for k in reverse(keys(saved.fields)))...)
            serialize(path, merge(saved, (fields = reordered,)))
            load_checkpoint!(target, path)
            @test io_hall_bytes(target) == io_hall_bytes(st)
        end
    end
end

@testset "Hall-MHD checkpoint rejects incompatible destinations before mutation" begin
    st = io_hall_run(Float64, Val(1), IsothermalElectrons(0.12))
    step_hall_mhd!(st, 0.001)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, st)
        for target in (
            HallMHDState(FourierGrid((4,), (2π,)), st.model),
            HallMHDState(FourierGrid((3,), (4π,)), st.model),
            HallMHDState(FourierGrid((3, 3), (2π, 2π)), st.model),
            HallMHDState(FourierGrid((3,), (Float32(2π),)), st.model),
            io_hall_run(Float64, Val(1), PolytropicElectrons(0.12, 1.0, 5 / 3)),
            io_hall_run(Float64, Val(1), IsothermalElectrons(0.12f0)),
        )
            io_hall_reject_unchanged(target, path)
        end
    end
end

@testset "Hall-MHD checkpoint validates the complete payload before mutation" begin
    for D = 1:3, closure in io_hall_closures(Float64)
        st = io_hall_run(Float64, Val(D), closure)
        step_hall_mhd!(st, 0.001)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st)
            saved = deserialize(path)
            target = io_hall_run(Float64, Val(D), closure; ηH = 0.004)
            step_hall_mhd!(target, 0.002)
            wrong_model = HallMHDModel(closure, -1.0, 0.0, 0.0, 1e-5)
            malformed = Any[
                nothing,
                merge(saved, (integrator = :rk4,)),
                merge(saved, (D = D + 1,)),
                merge(saved, (T = Float32,)),
                merge(saved, (model = nothing,)),
                merge(saved, (model = wrong_model,)),
                merge(saved, (time = 0.0f0,)),
                merge(saved, (time = NaN,)),
                merge(saved, (time = Inf,)),
                merge(saved, (time = -1.0,)),
                merge(saved, (step = -1,)),
                merge(saved, (step = 1.0,)),
                merge(saved, (floor_count = 0.0,)),
                merge(saved, (floor_count = -1,)),
                merge(saved, (floor_count = length(st.fields.n) + 1,)),
                merge(saved, (fields = nothing,)),
                merge(saved, (stage_lapJ = (),)),
                merge(saved, (stage_lapJ = ntuple(_ -> zeros(2), 3),)),
                merge(saved, (fields = merge(saved.fields, (n = zeros(Float32, st.g.n),)),)),
                merge(saved, (fields = merge(saved.fields, (ui = (saved.fields.ui[1],),)),)),
                merge(saved, (fields = merge(saved.fields, (gradp = (),)),)),
                merge(saved, (fields = merge(saved.fields, (pforce = [],)),)),
                merge(saved, (fields = merge(saved.fields, (cgl_dp = zeros(2),)),)),
                (; (k => v for (k, v) in pairs(saved) if k != :model)...),
                (; (k => v for (k, v) in pairs(saved) if k != :time)...),
            ]
            for bad in malformed
                serialize(path, bad)
                io_hall_reject_unchanged(target, path)
            end
            # A failed capture must preserve the existing file as well.
            save_checkpoint(path, st)
            bytes = read(path)
            st.time[] = NaN
            @test_throws ArgumentError save_checkpoint(path, st)
            @test read(path) == bytes
            @test_throws ArgumentError archive_run(path, st; rng_seed = 2)
            @test read(path) == bytes
        end
    end
end

mutable struct IOHallMutableClosure <: HybridPlasmaPIC.ElectronClosure
    temperature::Vector{Float64}
end
function HybridPlasmaPIC.electron_pressure!(pe, n, closure::IOHallMutableClosure)
    @. pe = closure.temperature[1] * n
    return pe
end

@testset "Hall-MHD checkpoint owns mutable closure configuration" begin
    closure = IOHallMutableClosure([0.12])
    st = io_hall_run(Float64, Val(1), closure)
    step_hall_mhd!(st, 0.001)
    expected = io_hall_bytes(st)
    mktemp() do path, io
        close(io)
        save_checkpoint(path, st)
        closure.temperature[1] = 0.7
        target = io_hall_run(Float64, Val(1), IOHallMutableClosure([0.2]))
        load_checkpoint!(target, path)
        @test io_hall_bytes(target) == expected
        @test target.model.closure.temperature !== closure.temperature
        closure.temperature[1] = 0.12
        step_hall_mhd!(st, 0.002)
        step_hall_mhd!(target, 0.002)
        @test io_hall_bytes(target) == io_hall_bytes(st)
    end
end

@testset "Hall-MHD archives preserve the canonical restart state and checksum" begin
    for T in (Float32, Float64), D = 1:3, closure in io_hall_closures(T)
        st = io_hall_run(T, Val(D), closure)
        step_hall_mhd!(st, T(0.001))
        mktempdir() do dir
            checkpoint = joinpath(dir, "checkpoint.bin")
            archive = joinpath(dir, "archive.bin")
            save_checkpoint(checkpoint, st)
            @test archive_run(archive, st; rng_seed = 7, diagnostic_desc = "Hall fields") == archive
            loaded = load_archive(archive)
            @test loaded.meta.rng_seed == 7
            @test loaded.state == deserialize(checkpoint)
            target = io_hall_run(T, Val(D), closure)
            HybridPlasmaPIC._restore_checkpoint_state!(target, loaded.state)
            @test io_hall_bytes(target) == io_hall_bytes(st)
            step_hall_mhd!(st, T(0.002); project_b = false)
            step_hall_mhd!(target, T(0.002); project_b = false)
            @test io_hall_bytes(target) == io_hall_bytes(st)
            packed = deserialize(archive)
            corrupt = copy(packed.state)
            corrupt[end] ⊻= 0x01
            serialize(archive, merge(packed, (state = corrupt,)))
            @test_throws ErrorException load_archive(archive)
        end
    end
end
