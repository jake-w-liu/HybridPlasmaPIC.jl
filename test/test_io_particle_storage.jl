using HybridPlasmaPIC, Test, Serialization

function io_storage_state(
    State,
    ::Type{T},
    ::Val{D},
    count;
    mobile = false,
    target = false,
) where {T,D}
    g = FourierGrid(ntuple(_ -> 3, D), ntuple(_ -> T(2π), D))
    if State in (HybridStepper, CAMCLStepper)
        model = HybridModel(IsothermalElectrons(T(target ? 0.2 : 0.1)))
        return State(g, model, CIC(), count)
    elseif State in (EMPIC1D, EMPIC)
        return State(g, count; mobile, mi = T(16), c = T(target ? 4 : 3))
    end
    return State(g, count)
end

function io_storage_particles(::Type{T}, ::Val{D}, count; ion = false, hybrid = false) where {T,D}
    ps = ParticleSet{D,T}(count; q = ion || hybrid ? one(T) : -one(T), m = ion ? T(16) : one(T))
    for d = 1:D
        ps.x[d] .= T.(1:count) .* T(0.1d)
    end
    for c = 1:3
        ps.v[c] .= T.(1:count) .* T(0.01c)
    end
    ps.weight .= T(0.3)
    ps.id .= UInt64.(11:(10+count))
    ps.tag .= UInt32(7)
    return ps
end

function io_storage_alias(ps::ParticleSet{D,T}, kind) where {D,T}
    x, v, weight, id, tag = ps.x, ps.v, ps.weight, ps.id, ps.tag
    if kind === :position_velocity
        x = (v[1], Base.tail(x)...)
    elseif kind === :position_position
        x = ntuple(d -> d == D ? x[1] : x[d], D)
    elseif kind === :velocity_velocity
        v = (v[1], v[1], v[3])
    elseif kind === :weight_velocity
        weight = v[3]
    elseif kind === :weight_provenance
        weight = reinterpret(T, T === Float64 ? id : tag)
    elseif kind === :id_tag
        tag = view(reinterpret(UInt32, id), 1:length(tag))
    elseif kind === :indices
        weight = view(ones(T, 3), id)
    else
        error("unknown storage test case")
    end
    return ParticleSet{D,T}(x, v, weight, id, tag, ps.q, ps.m)
end

io_storage_arrays(ps) = (ps.x..., ps.v..., ps.weight, ps.id, ps.tag)
io_storage_buffers(st::Union{HybridStepper,CAMCLStepper}) =
    (st.Ep..., st.Bp..., st.vpred..., st.xmid..., st.work)
io_storage_buffers(st::Electrostatic1D) = (st.Ep,)
io_storage_buffers(st::ElectrostaticPIC) = st.Ep
io_storage_buffers(st::EMPIC1D) =
    (st.Exp, st.Eyp, st.Bzp, st.work, st.Exi, st.Eyi, st.Bzi, st.worki)
io_storage_buffers(st::EMPIC) =
    (st.Ep..., st.Bp..., st.mide..., st.worke, st.Epi..., st.Bpi..., st.midi..., st.worki)

function io_storage_snapshot(st, populations...)
    io = IOBuffer()
    serialize(io, (HybridPlasmaPIC._checkpoint_state(st, populations...), io_storage_buffers(st)))
    return take!(io)
end

function io_storage_reject_unchanged(st, path, populations...)
    before = io_storage_snapshot(st, populations...)
    buffers = io_storage_buffers(st)
    arrays = map(io_storage_arrays, populations)
    @test_throws r"particle components must not alias" load_checkpoint!(st, populations..., path)
    @test io_storage_snapshot(st, populations...) == before
    @test all(a === b for (a, b) in zip(buffers, io_storage_buffers(st)))
    @test all(
        all(a === b for (a, b) in zip(old, io_storage_arrays(ps))) for
        (old, ps) in zip(arrays, populations)
    )
end

const IO_STORAGE_FAMILIES = (
    (HybridStepper, 1),
    (HybridStepper, 2),
    (HybridStepper, 3),
    (CAMCLStepper, 1),
    (CAMCLStepper, 2),
    (CAMCLStepper, 3),
    (Electrostatic1D, 1),
    (ElectrostaticPIC, 1),
    (ElectrostaticPIC, 2),
    (ElectrostaticPIC, 3),
    (EMPIC1D, 1),
    (EMPIC, 1),
    (EMPIC, 2),
    (EMPIC, 3),
)

@testset "Particle storage predicate checks every component pair" begin
    for T in (Float32, Float64), D = 1:3, count in (0, 2)
        types = (ntuple(_ -> T, D + 4)..., UInt64, UInt32)
        # Byte-backed arrays let every pair share storage even across element types.
        function population(first = 0, second = 0)
            backing = [zeros(UInt8, count * sizeof(UInt64)) for _ in types]
            first == 0 || (backing[second] = backing[first])
            arrays = map(types, Tuple(backing)) do element, storage
                reinterpret(element, view(storage, 1:(count*sizeof(element))))
            end
            return ParticleSet{D,T}(
                ntuple(d -> arrays[d], D),
                ntuple(c -> arrays[D+c], 3),
                arrays[D+4],
                arrays[D+5],
                arrays[D+6],
                one(T),
                one(T),
            )
        end
        @test !HybridPlasmaPIC._particle_components_mightalias(population())
        for second = 2:length(types), first = 1:(second-1)
            @test HybridPlasmaPIC._particle_components_mightalias(population(first, second))
        end
    end
end

@testset "Checkpoint preflight rejects overlapping components before mutation" begin
    for T in (Float32, Float64), (State, D) in IO_STORAGE_FAMILIES, count in (0, 3)
        hybrid = State in (HybridStepper, CAMCLStepper)
        source = io_storage_state(State, T, Val(D), count)
        particles = io_storage_particles(T, Val(D), count; hybrid)
        target = io_storage_state(State, T, Val(D), 1; target = true)
        kinds =
            [:position_velocity, :velocity_velocity, :weight_velocity, :weight_provenance, :id_tag]
        D > 1 && push!(kinds, :position_position)
        count == 0 && push!(kinds, :indices)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, source, particles)
            for kind in kinds, destination_count in unique((min(count, 1), count))
                destination = io_storage_alias(
                    io_storage_particles(T, Val(D), destination_count; hybrid),
                    kind,
                )
                @test HybridPlasmaPIC._particle_components_mightalias(destination)
                io_storage_reject_unchanged(target, path, destination)
            end
        end
    end
end

@testset "Mobile checkpoint preflights ion components before resizing electrons" begin
    for T in (Float32, Float64),
        (State, D) in ((EMPIC1D, 1), (EMPIC, 1), (EMPIC, 2), (EMPIC, 3)),
        count in (0, 3)

        source = io_storage_state(State, T, Val(D), 2; mobile = true)
        electrons = io_storage_particles(T, Val(D), 2)
        ions = io_storage_particles(T, Val(D), count; ion = true)
        target = io_storage_state(State, T, Val(D), 1; mobile = true, target = true)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, source, electrons, ions)
            for kind in (:position_velocity, :weight_provenance, :id_tag)
                destination_e = io_storage_particles(T, Val(D), 1)
                destination_i =
                    io_storage_alias(io_storage_particles(T, Val(D), count; ion = true), kind)
                io_storage_reject_unchanged(target, path, destination_e, destination_i)
            end
        end
    end
end

function io_storage_disjoint_views(ps::ParticleSet{D,T}) where {D,T}
    count = nparticles(ps)
    storage = fill(T(-7), (D + 4) * (count + 1))
    # One sentinel separates each disjoint view; all share a single writable parent.
    component(i) = view(storage, ((i-1)*(count+1)+1):((i-1)*(count+1)+count))
    x = ntuple(component, D)
    v = ntuple(c -> component(D + c), 3)
    destination = ParticleSet{D,T}(
        x,
        v,
        component(D + 4),
        view(copy(ps.id), :),
        view(copy(ps.tag), :),
        zero(T),
        T(NaN),
    )
    return destination, storage
end

io_storage_init!(st::HybridStepper, ps) = init!(st, ps)
io_storage_init!(st::CAMCLStepper, ps) = init_camcl!(st, ps)
io_storage_init!(st::Union{Electrostatic1D,ElectrostaticPIC}, ps) = init_espic!(st, ps)
io_storage_init!(st::Union{EMPIC1D,EMPIC}, populations...) = init_empic!(st, populations...)
io_storage_step!(st::HybridStepper, ps, dt) = step!(st, ps, dt; NB = 2)
io_storage_step!(st::CAMCLStepper, ps, dt) = step_camcl!(st, ps, dt; NB = 2)
io_storage_step!(st::Union{Electrostatic1D,ElectrostaticPIC}, ps, dt) = step_espic!(st, ps, dt)
io_storage_step!(st::Union{EMPIC1D,EMPIC}, ps, dt) = step_empic!(st, ps, dt)

@testset "Checkpoint restores independent writable component views in place" begin
    for T in (Float32, Float64), (State, D) in IO_STORAGE_FAMILIES
        hybrid = State in (HybridStepper, CAMCLStepper)
        count = 3^D
        source = io_storage_state(State, T, Val(D), count)
        particles = io_storage_particles(T, Val(D), count; hybrid)
        load_lattice!(particles, ntuple(_ -> zero(T), D), source.g.L, source.g.n)
        set_density_weight!(particles, one(T), source.g)
        foreach(a -> fill!(a, T(0.001)), particles.v)
        io_storage_init!(source, particles)
        io_storage_step!(source, particles, T(0.001))
        target = io_storage_state(State, T, Val(D), 1; target = true)
        destination, backing = io_storage_disjoint_views(particles)
        arrays = io_storage_arrays(destination)
        @test !HybridPlasmaPIC._particle_components_mightalias(destination)
        mktemp() do path, io
            close(io)
            save_checkpoint(path, source, particles)
            @test load_checkpoint!(target, destination, path) === target
            @test isequal(
                HybridPlasmaPIC._checkpoint_state(target, destination),
                HybridPlasmaPIC._checkpoint_state(source, particles),
            )
            @test all(a === b for (a, b) in zip(arrays, io_storage_arrays(destination)))
            @test all(==(T(-7)), backing[(count+1):(count+1):end])
            for dt in T[0.001, 0.0007]
                io_storage_step!(source, particles, dt)
                io_storage_step!(target, destination, dt)
                @test isequal(
                    HybridPlasmaPIC._checkpoint_state(target, destination),
                    HybridPlasmaPIC._checkpoint_state(source, particles),
                )
            end
        end
    end
end
