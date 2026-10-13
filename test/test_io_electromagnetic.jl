using HybridPlasmaPIC, Test, Serialization

function io_em_run(
    State,
    ::Type{T},
    ::Val{D};
    shape = CIC(),
    mobile = false,
    relativistic = false,
) where {T,D}
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
    function population(count, charge, mass)
        counts = ntuple(_ -> count, D)
        ps = ParticleSet{D,T}(prod(counts); q = charge, m = mass)
        load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, counts)
        set_density_weight!(ps, T(0.7), g)
        ps.weight .*= one(T) .+ T(0.01) .* cos.(ps.x[1])
        for c = 1:3
            ps.v[c] .= T(0.01c) .* cos.(ps.x[1])
        end
        return ps
    end
    electrons = population(2, -one(T), one(T))
    ions = mobile ? population(3, one(T), T(16)) : nothing
    st = State(
        g,
        nparticles(electrons);
        n0 = T(0.7),
        c = T(3),
        mi = T(16),
        shape,
        mobile,
        relativistic,
        n_sub = 2,
    )
    for I in CartesianIndices(g.n)
        x = T(I[1] - 1) * g.dx[1]
        if st isa EMPIC1D
            st.Ey[I] = T(0.01) * sin(x)
            st.Bz[I] = T(0.02) * cos(x)
        else
            st.E[2][I] = T(0.01) * sin(x)
            st.B[3][I] = T(0.02) * cos(x)
        end
    end
    init_empic!(st, electrons, ions)
    return st, electrons, ions
end

function io_em_serialized(value)
    io = IOBuffer()
    serialize(io, value)
    return take!(io)
end

io_em_particle_values(ps) =
    ps === nothing ? nothing : (ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m)
io_em_fields(st::EMPIC1D) = (st.Ex, st.Ey, st.Bz, st.Jx, st.Jy, st.rho_n, st.rho_np1)
io_em_fields(st::EMPIC) = (st.E, st.B, st.J, st.rho_n, st.rho_np1)

function io_em_bytes(st, electrons, ions = nothing)
    # Read the solver directly so omission from checkpoint capture cannot make
    # the restart comparison pass merely by omitting the same state on both sides.
    return io_em_serialized((
        typeof(st),
        st.g.n,
        st.g.L,
        st.n0,
        st.c,
        st.shape,
        st.relativistic,
        st.mobile,
        st.mi,
        st.n_sub,
        io_em_particle_values(electrons),
        io_em_particle_values(ions),
        io_em_fields(st),
        st.time[],
        st.step[],
        st.last_dt[],
        st.last_sub_dt[],
    ))
end

function io_em_workspaces(st::EMPIC1D)
    return (
        map(length, (st.Exp, st.Eyp, st.Bzp, st.work)),
        map(length, (st.Exi, st.Eyi, st.Bzi, st.worki)),
    )
end

function io_em_workspaces(st::EMPIC)
    return (
        map(length, (st.Ep..., st.Bp..., st.mide..., st.worke)),
        map(length, (st.Epi..., st.Bpi..., st.midi..., st.worki)),
    )
end

function io_em_reject_unchanged(st, electrons, ions, path)
    before = io_em_bytes(st, electrons, ions)
    workspaces = io_em_workspaces(st)
    @test_throws ArgumentError load_checkpoint!(st, electrons, ions, path)
    @test io_em_bytes(st, electrons, ions) == before
    @test io_em_workspaces(st) == workspaces
end

@testset "EM checkpoint preserves species, diagnostics, configuration and phase" begin
    for T in (Float32, Float64),
        (State, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
        shape in (NGP(), CIC(), TSC()),
        mobile in (false, true),
        relativistic in (false, true),
        scenario in (:initial, :live, :reset)

        st, electrons, ions = io_em_run(State, T, Val(D); shape, mobile, relativistic)
        if scenario != :initial
            step_empic!(st, electrons, ions, T(0.008))
            st.n_sub = 3
            step_empic!(st, electrons, ions, T(0.005))
        end
        scenario == :reset && init_empic!(st, electrons, ions)
        # A setting changed between steps must not rewrite the saved fast half level.
        st.n_sub = 5
        mktempdir() do dir
            path = joinpath(dir, "checkpoint.bin")
            @test save_checkpoint(path, st, electrons, ions) == path
            restored = State(
                FourierGrid(st.g.n, st.g.L),
                1;
                n0 = T(0.2),
                c = T(5),
                mi = T(25),
                shape = typeof(shape)(),
                mobile,
                relativistic = !relativistic,
            )
            e2 = ParticleSet{D,T}(1; q = -one(T))
            i2 = mobile ? ParticleSet{D,T}(1; q = one(T), m = T(25)) : nothing
            @test load_checkpoint!(restored, e2, i2, path) === restored
            @test io_em_bytes(restored, e2, i2) == io_em_bytes(st, electrons, ions)
            @test all(==(nparticles(electrons)), io_em_workspaces(restored)[1])
            @test all(==(mobile ? nparticles(ions) : 0), io_em_workspaces(restored)[2])
            @test em_field_energy(restored) === em_field_energy(st)
            diagnostic_dt = st.step[] == 0 ? T(0.001) : st.last_sub_dt[]
            @test charge_conservation_residual(restored, diagnostic_dt) ===
                  charge_conservation_residual(st, diagnostic_dt)
            if scenario == :live && shape isa CIC
                archive = joinpath(dir, "archive.bin")
                @test archive_run(archive, st, electrons, ions; rng_seed = 19) == archive
                archived = load_archive(archive)
                @test io_em_serialized(archived.state) == io_em_serialized(deserialize(path))
                @test archived.meta.rng_seed == 19
            end
            for (dt, ns) in ((T(0.007), 5), (T(0.003), 1), (zero(T), 4), (T(0.009), 4))
                st.n_sub = restored.n_sub = ns
                step_empic!(st, electrons, ions, dt)
                step_empic!(restored, e2, i2, dt)
                @test io_em_bytes(restored, e2, i2) == io_em_bytes(st, electrons, ions)
            end
        end
    end
end

@testset "EM checkpoint preserves externally edited diagnostic and field arrays" begin
    for State in (EMPIC1D, EMPIC)
        st, electrons, ions = io_em_run(State, Float64, Val(1); mobile = true)
        step_empic!(st, electrons, ions, 0.01)
        for arrays in io_em_fields(st)
            for a in (arrays isa Tuple ? arrays : (arrays,))
                a[1] += 0.125
            end
        end
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, electrons, ions)
            saved = deserialize(path)
            reordered = (;
                extra = :unused,
                (
                    name => getproperty(saved.fields, name) for name in reverse(keys(saved.fields))
                )...,
            )
            serialize(path, merge(saved, (fields = reordered,)))
            restored, e2, i2 = io_em_run(State, Float64, Val(1); mobile = true)
            load_checkpoint!(restored, e2, i2, path)
            @test io_em_bytes(restored, e2, i2) == io_em_bytes(st, electrons, ions)
            @test charge_conservation_residual(restored, st.last_sub_dt[]) ===
                  charge_conservation_residual(st, st.last_sub_dt[])
            step_empic!(st, electrons, ions, 0.006)
            step_empic!(restored, e2, i2, 0.006)
            @test io_em_bytes(restored, e2, i2) == io_em_bytes(st, electrons, ions)
        end
    end
end

@testset "EM checkpoint rejects incompatible destinations before mutation" begin
    for State in (EMPIC1D, EMPIC)
        st, electrons, ions = io_em_run(State, Float64, Val(1); mobile = true)
        step_empic!(st, electrons, ions, 0.01)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, electrons, ions)
            Other = State === EMPIC1D ? EMPIC : EMPIC1D
            targets = (
                Other(st.g, 1; mobile = true, mi = 16.0),
                State(st.g, 1; mobile = false, mi = 16.0),
                State(st.g, 1; mobile = true, mi = 16.0, shape = TSC()),
                State(FourierGrid((8,), st.g.L), 1; mobile = true, mi = 16.0),
                State(FourierGrid(st.g.n, (4π,)), 1; mobile = true, mi = 16.0),
                State(FourierGrid(st.g.n, (Float32(2π),)), 1; mobile = true, mi = 16.0),
                EMPIC(FourierGrid((4, 4), (2π, 2π)), 1; mobile = true, mi = 16.0),
            )
            for target in targets
                D, T = length(target.g.n), typeof(target.n0)
                e2 = ParticleSet{D,T}(1; q = -one(T))
                i2 = target.mobile ? ParticleSet{D,T}(1; m = T(16)) : nothing
                foreach(a -> fill!(a, T(0.3)), e2.x)
                i2 === nothing || foreach(a -> fill!(a, T(0.7)), i2.x)
                io_em_reject_unchanged(target, e2, i2, path)
            end
            target, e2, i2 = io_em_run(State, Float64, Val(1); mobile = true)
            io_em_reject_unchanged(target, e2, nothing, path)
            for bad in (ParticleSet{2,Float64}(1), ParticleSet{1,Float32}(1))
                foreach(a -> fill!(a, 0.4), bad.x)
                io_em_reject_unchanged(target, bad, i2, path)
                io_em_reject_unchanged(target, e2, bad, path)
            end
            before = io_em_bytes(target, e2, i2)
            @test_throws ArgumentError load_checkpoint!(target, e2, path)
            @test io_em_bytes(target, e2, i2) == before
        end
    end
end

@testset "EM checkpoint validates both species, fields, settings and phase" begin
    for State in (EMPIC1D, EMPIC)
        st, electrons, ions = io_em_run(State, Float64, Val(1); mobile = true, relativistic = true)
        step_empic!(st, electrons, ions, 0.01)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, electrons, ions)
            saved = deserialize(path)
            target = State(st.g, 1; mobile = true, mi = 25.0, n0 = 2.0, c = 6.0)
            e2 = ParticleSet{1,Float64}(1; q = -1.0)
            i2 = ParticleSet{1,Float64}(1; m = 25.0)
            e2.x[1][1] = 0.3
            i2.x[1][1] = 0.7
            field_name = first(keys(HybridPlasmaPIC._empic_checkpoint_fields(st)))
            bad_fields = merge(saved.fields, NamedTuple{(field_name,)}((zeros(Float32, 3),)))
            fast_ion = merge(
                saved.ions,
                (v = (fill(saved.c, length(saved.ions.weight)), saved.ions.v[2], saved.ions.v[3]),),
            )
            malformed = (
                merge(saved, (integrator = :unsupported,)),
                merge(saved, (n0 = NaN,)),
                merge(saved, (n0 = -1.0,)),
                merge(saved, (c = 0.0,)),
                merge(saved, (c = 3.0f0,)),
                merge(saved, (mi = -1.0,)),
                merge(saved, (mi = 16.0f0,)),
                merge(saved, (mobile = 1,)),
                merge(saved, (relativistic = 1,)),
                merge(saved, (n_sub = 0,)),
                merge(saved, (n_sub = 2.0,)),
                merge(saved, (time = NaN,)),
                merge(saved, (step = -1,)),
                merge(saved, (last_dt = 0.0,)),
                merge(saved, (last_sub_dt = 0.0,)),
                merge(saved, (last_sub_dt = 0.1,)),
                merge(saved, (last_sub_dt = NaN,)),
                merge(saved, (last_sub_dt = 0.005f0,)),
                merge(saved, (ions = nothing,)),
                merge(saved, (ions = fast_ion,)),
                merge(saved, (electrons = merge(saved.electrons, (q = 1.0,)),)),
                merge(saved, (electrons = merge(saved.electrons, (m = 2.0,)),)),
                merge(saved, (ions = merge(saved.ions, (q = -1.0,)),)),
                merge(saved, (ions = merge(saved.ions, (m = 25.0,)),)),
                merge(saved, (ions = merge(saved.ions, (m = 16,)),)),
                merge(saved, (ions = merge(saved.ions, (weight = [1.0],)),)),
                merge(saved, (ions = merge(saved.ions, (tag = Float64.(saved.ions.tag),)),)),
                merge(saved, (fields = bad_fields,)),
                merge(saved, (fields = merge(saved.fields, (rho_np1 = zeros(3),)),)),
                (; (k => v for (k, v) in pairs(saved) if k != :last_sub_dt)...),
            )
            for bad in malformed
                serialize(path, bad)
                io_em_reject_unchanged(target, e2, i2, path)
            end
            # Failed captures leave existing checkpoint/archive files untouched.
            serialize(path, saved)
            bytes = read(path)
            @test_throws ArgumentError save_checkpoint(path, st, electrons)
            @test read(path) == bytes
            ions.m = 25.0
            @test_throws ArgumentError save_checkpoint(path, st, electrons, ions)
            @test read(path) == bytes
            @test_throws ArgumentError archive_run(path, st, electrons, ions; rng_seed = 4)
            @test read(path) == bytes
        end
    end
end

function io_em_views(ps::ParticleSet{D,T}; short = false) where {D,T}
    indices = short ? (1:1) : Colon()
    return ParticleSet{D,T}(
        map(a -> view(a, indices), ps.x),
        map(a -> view(a, indices), ps.v),
        view(ps.weight, indices),
        view(ps.id, indices),
        view(ps.tag, indices),
        ps.q,
        ps.m,
    )
end

@testset "EM checkpoint preflights both destinations and supports independent views" begin
    for State in (EMPIC1D, EMPIC)
        st, electrons, ions = io_em_run(State, Float64, Val(1); mobile = true)
        step_empic!(st, electrons, ions, 0.01)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, electrons, ions)
            target, e2, i2 = io_em_run(State, Float64, Val(1); mobile = true)
            for ps in (e2, i2), a in (ps.x..., ps.v..., ps.weight)
                fill!(a, 0.2)
            end
            ev, iv = io_em_views(e2), io_em_views(i2)
            @test load_checkpoint!(target, ev, iv, path) === target
            @test io_em_bytes(target, e2, i2) == io_em_bytes(st, electrons, ions)
            step_empic!(target, ev, iv, 0.005)
            step_empic!(st, electrons, ions, 0.005)
            @test io_em_bytes(target, e2, i2) == io_em_bytes(st, electrons, ions)
            # Invalid ion resizing must be rejected before electron resizing/configuration.
            target = State(st.g, 1; mobile = true, mi = 25.0, c = 6.0)
            tiny_e = ParticleSet{1,Float64}(1; q = -1.0)
            tiny_e.x[1][1] = 0.3
            io_em_reject_unchanged(target, tiny_e, io_em_views(i2; short = true), path)
            io_em_reject_unchanged(target, io_em_views(e2; short = true), i2, path)
            for name in (:x, :v, :weight, :id, :tag)
                base = ParticleSet{1,Float64}(1; m = 16.0)
                base.x[1][1] = 0.7
                alias = ParticleSet{1,Float64}(
                    name == :x ? tiny_e.x : base.x,
                    name == :v ? tiny_e.v : base.v,
                    name == :weight ? tiny_e.weight : base.weight,
                    name == :id ? tiny_e.id : base.id,
                    name == :tag ? tiny_e.tag : base.tag,
                    1.0,
                    16.0,
                )
                io_em_reject_unchanged(target, tiny_e, alias, path)
            end
        end
    end
end

mutable struct IOEMCheckpointShape <: HybridPlasmaPIC.ShapeFunction
    parameters::Vector{Float64}
end

@testset "EM checkpoint owns mutable shape configuration" begin
    for State in (EMPIC1D, EMPIC)
        shape = IOEMCheckpointShape([1.0, 2.0])
        st = State(FourierGrid((4,), (2π,)), 1; shape)
        electrons = ParticleSet{1,Float64}(1; q = -1.0)
        electrons.x[1][1] = 0.3
        mktemp() do path, io
            close(io)
            save_checkpoint(path, st, electrons)
            shape.parameters[1] = 8.0
            load_checkpoint!(st, electrons, path)
            @test st.shape !== shape
            @test st.shape.parameters == [1.0, 2.0]
        end
    end
end

@testset "EM checkpoint electron-only API, ion presence and empty populations" begin
    for State in (EMPIC1D, EMPIC), mobile in (false, true)
        st = State(FourierGrid((4,), (2π,)), 0; mobile, mi = 16.0)
        electrons = ParticleSet{1,Float64}(0; q = -1.0)
        ions = mobile ? ParticleSet{1,Float64}(0; m = 16.0) : nothing
        init_empic!(st, electrons, ions)
        mktempdir() do dir
            path = joinpath(dir, "checkpoint.bin")
            mobile ? save_checkpoint(path, st, electrons, ions) :
            save_checkpoint(path, st, electrons)
            target = State(st.g, 1; mobile, mi = 25.0)
            e2 = ParticleSet{1,Float64}(1; q = -1.0)
            i2 = mobile ? ParticleSet{1,Float64}(1; m = 25.0) : nothing
            mobile ? load_checkpoint!(target, e2, i2, path) : load_checkpoint!(target, e2, path)
            @test io_em_bytes(target, e2, i2) == io_em_bytes(st, electrons, ions)
            @test all(iszero, io_em_workspaces(target)[1])
            @test all(iszero, io_em_workspaces(target)[2])
            step_empic!(st, electrons, ions, 0.01)
            step_empic!(target, e2, i2, 0.01)
            @test io_em_bytes(target, e2, i2) == io_em_bytes(st, electrons, ions)
            if !mobile
                extra = ParticleSet{1,Float64}(0; m = 16.0)
                io_em_reject_unchanged(target, e2, extra, path)
                bytes = read(path)
                @test_throws ArgumentError save_checkpoint(path, st, electrons, extra)
                @test read(path) == bytes
                archive = joinpath(dir, "archive.bin")
                @test archive_run(archive, st, electrons; rng_seed = 5) == archive
                @test load_archive(archive).state.ions === nothing
            end
        end
    end
end
