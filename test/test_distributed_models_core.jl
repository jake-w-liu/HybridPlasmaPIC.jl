using HybridPlasmaPIC, PencilArrays, PencilFFTs, MPI, Test

# These are actual physical-pencil models, distinct from the replicated MPI
# reference. A dyadic charged CIC load gives an independent initial density
# oracle; vacuum Maxwell additionally has an analytic single-mode KDK oracle.
function _pencil_charged_species(comm)
    ps = ParticleSet{3,Float64}(MPI.Comm_rank(comm) == 0 ? 64 : 0; q = 1, m = 1)
    pattern = (1, 0, -1, 0)
    if MPI.Comm_rank(comm) == 0
        p = 0
        for k = 0:3, j = 0:3, i = 0:3
            p += 1
            ps.x[1][p], ps.x[2][p], ps.x[3][p] = (
                (i + 7 / 8) / 4,
                (j + 7 / 8) / 4,
                (k + 7 / 8) / 4,
            )
            for c = 1:3
                ps.v[c][p] = 2
            end
            ps.weight[p] = (1 + pattern[i+1] / 8) / 64
            ps.id[p], ps.tag[p] = UInt64(p), UInt32(17 + p)
        end
    end
    return [ps]
end

Base.@nospecializeinfer function _pencil_fill_mode!(field, @nospecialize(value))
    values = PencilArrays.global_view(field)
    for I in CartesianIndices(values)
        values[I] = value(Tuple(I) .- 1)
    end
    return field
end

function _pencil_initial_density(density, background)
    pattern = (1, 0, -1, 0)
    values = PencilArrays.global_view(density)
    for I in CartesianIndices(values)
        # CIC weights are 1/8 and 7/8 on each axis. Uniform y/z sums are one;
        # x receives its own marker and the preceding cell's marker.
        expected = background + (pattern[I[1]] + 7 * pattern[mod1(I[1] - 1, 4)]) / 64
        @test values[I] == expected
    end
end

Base.@nospecializeinfer function _pencil_boundary(
    @nospecialize(fields),
    @nospecialize(species),
    @nospecialize(phase),
    @nospecialize(caches),
)
    particles = map(species) do ps
        (map(copy, ps.x), map(copy, ps.v), copy(ps.weight), copy(ps.id), copy(ps.tag), ps.q, ps.m)
    end
    return (map(a -> copy(parent(a)), fields), particles, map(r -> r[], phase), map(copy, caches))
end

Base.@nospecializeinfer function _pencil_inventory(@nospecialize(species), field, comm)
    counts = zeros(Int, 64)
    owned = PencilArrays.range_local(field, PencilArrays.LogicalOrder())
    pattern = (1, 0, -1, 0)
    for ps in species
        @test ps.q == ps.m == 1
        for p in eachindex(ps.weight)
            valid_id = 1 <= ps.id[p] <= 64
            @test valid_id
            if valid_id
                id = Int(ps.id[p])
                counts[id] += 1
                @test ps.tag[p] == UInt32(17 + id)
                @test ps.weight[p] == (1 + pattern[mod1(id, 4)] / 8) / 64
            end
            for d = 1:3
                @test 0 <= ps.x[d][p] < 1
                @test floor(Int, 4 * ps.x[d][p]) + 1 in owned[d]
                @test isfinite(ps.v[d][p])
            end
        end
    end
    @test MPI.Allreduce(counts, +, comm) == ones(Int, 64)
    return nothing
end

_pencil_magnitude(fields, comm) =
    MPI.Allreduce(maximum(a -> maximum(abs, parent(a)), fields), MPI.MAX, comm)

# Only the deliberately injected validation boundary captures expected errors;
# unexpected runtime failures propagate instead of becoming a successful check.
Base.@nospecializeinfer function _pencil_reject_step(
    @nospecialize(advance),
    @nospecialize(fields),
    @nospecialize(species),
    @nospecialize(phase),
    @nospecialize(caches),
    comm,
)
    before = _pencil_boundary(fields, species, phase, caches)
    origin = MPI.Comm_rank(comm) == 1
    caught = try
        advance(origin ? -1 / 32 : 1 / 32)
        nothing
    catch err
        err isa Union{ArgumentError,ErrorException} || rethrow()
        err
    end
    @test caught isa (origin ? ArgumentError : ErrorException)
    if !origin && caught isa ErrorException
        @test occursin("failed on another rank", sprint(showerror, caught))
    end
    @test isequal(_pencil_boundary(fields, species, phase, caches), before)
    @test MPI.Allreduce(1, +, comm) == MPI.Comm_size(comm)
    return nothing
end

Base.@nospecializeinfer function _pencil_checkpoint_continuation(
    @nospecialize(advance),
    @nospecialize(model),
    @nospecialize(fields),
    @nospecialize(species),
    @nospecialize(phase),
    @nospecialize(caches),
    comm,
)
    saved = _pencil_boundary(fields, species, phase, caches)
    owners, parents = copy(species), map(parent, fields)
    directory = MPI.bcast(MPI.Comm_rank(comm) == 0 ? mktempdir() : "", 0, comm)
    try
        manifest = save_mpi_checkpoint(directory, model)
        @test isfile(manifest)
        @test advance(3 / 128) === model
        expected = _pencil_boundary(fields, species, phase, caches)
        # Valid destination representation survives while replaceable values,
        # local lengths and cached counts deliberately cease to be ready state.
        foreach(a -> fill!(parent(a), NaN), fields)
        model.time[], model.step[] = NaN, -1
        for ps in species
            foreach(empty!, (ps.x..., ps.v..., ps.weight, ps.id, ps.tag))
            ps.q, ps.m = NaN, NaN
        end
        if !isempty(species)
            model.initialized[], model.primed[], model.last_dt[] = false, true, NaN
            hasproperty(model, :_population_counts) && fill!(model._population_counts, -1)
        end
        @test load_mpi_checkpoint!(model, directory) === model
        @test isequal(_pencil_boundary(fields, species, phase, caches), saved)
        @test all(i -> parent(fields[i]) === parents[i], eachindex(fields))
        @test model.time === phase[end-1] && model.step === phase[end]
        if !isempty(species)
            @test all(
                s -> model.owners[s] === owners[s] && species[s] === owners[s],
                eachindex(species),
            )
            @test model.initialized === phase[1] && model.primed === phase[2]
            @test model.last_dt === phase[3] && model._species_scalars === caches[1]
            if hasproperty(model, :_population_counts)
                @test model._population_counts === caches[2]
            end
        end
        @test advance(3 / 128) === model
        @test isequal(_pencil_boundary(fields, species, phase, caches), expected)
        @test model.time[] == 7 / 128 && model.step[] == 2
        @test all(a -> all(isfinite, parent(a)), fields)
    finally
        MPI.Barrier(comm)
        MPI.Comm_rank(comm) == 0 && rm(directory; recursive = true)
    end
    return nothing
end

Base.@nospecializeinfer function _pencil_particle_lifecycle(
    @nospecialize(advance),
    @nospecialize(diagnose),
    @nospecialize(model),
    @nospecialize(fields),
    @nospecialize(species),
    @nospecialize(caches),
    comm,
)
    phase = (model.initialized, model.primed, model.last_dt, model.time, model.step)
    @test model.initialized[] && !model.primed[]
    @test model.time[] == model.last_dt[] == model.step[] == 0
    @test MPI.Allreduce(Int(!isempty(species[1].weight)), +, comm) == MPI.Comm_size(comm)
    _pencil_inventory(species, fields[1], comm)
    owners = copy(species)
    before = _pencil_boundary(fields, species, phase, caches)
    @test advance(0.0) === model
    @test isequal(_pencil_boundary(fields, species, phase, caches), before)
    _pencil_reject_step(advance, fields, species, phase, caches, comm)
    identifiers = Set(species[1].id)
    @test advance(1 / 32) === model
    @test model.primed[] && model.last_dt[] == model.time[] == 1 / 32
    @test model.step[] == 1
    @test !isequal(_pencil_boundary(fields, species, phase, caches), before)
    @test all(s -> species[s] === owners[s] && model.owners[s] === owners[s], eachindex(species))
    changed_fields = any(i -> !isequal(parent(fields[i]), before[1][i]), eachindex(fields))
    accelerated = any(ps -> any(a -> any(!=(2), a), ps.v), species)
    @test MPI.Allreduce(Int(changed_fields), MPI.MAX, comm) == 1
    @test MPI.Allreduce(Int(accelerated), MPI.MAX, comm) == 1
    @test all(a -> all(isfinite, parent(a)), fields)
    @test MPI.Allreduce(Int(Set(species[1].id) != identifiers), MPI.MAX, comm) == 1
    _pencil_inventory(species, fields[1], comm)
    diagnostics = diagnose()
    @test diagnostics.particle_count == 64 && diagnostics.particle_charge == 1
    @test diagnostics.time == 1 / 32 && diagnostics.step == 1
    @test diagnostics.velocity_time_offset == -1 / 64
    @test isfinite(diagnostics.kinetic) && isfinite(diagnostics.total)
    if hasproperty(diagnostics, :physical_neutrality_residual)
        @test diagnostics.physical_neutrality_residual == 0
    end
    _pencil_checkpoint_continuation(advance, model, fields, species, phase, caches, comm)
    _pencil_inventory(species, fields[1], comm)
    @test all(s -> species[s] === owners[s] && model.owners[s] === owners[s], eachindex(species))
    return nothing
end

@testset "Actual-pencil charged model lifecycles" begin
    ensure_mpi_initialized!()
    comm = MPI.COMM_WORLD
    @test MPI.Comm_size(comm) in (2, 4, 8)
    @testset "Electrostatic" begin
        g = distributed_fourier_grid((4, 4, 4), (1.0, 1.0, 1.0); comm)
        species = _pencil_charged_species(comm)
        owner = species[1]
        model = distributed_electrostatic_pic(
            g,
            species;
            grid_dc = :project,
            background_charge = -1.0,
            charge_rtol = 0.0,
        )
        @test init_espic!(model, species) === model
        @test species[1] === owner && model.owners[1] === owner
        _pencil_initial_density(model.rho, 0.0)
        @test _pencil_magnitude(model.E, comm) > 0
        _pencil_particle_lifecycle(
            h -> step_espic!(model, species, h),
            () -> electrostatic_diagnostics(model, species),
            model,
            (model.E..., model.rho),
            species,
            (model._species_scalars,),
            comm,
        )
    end
    @testset "Electromagnetic" begin
        g = distributed_fourier_grid((4, 4, 4), (1.0, 1.0, 1.0); comm)
        species = _pencil_charged_species(comm)
        owner = species[1]
        model = distributed_empic(
            g,
            species;
            grid_dc = :project,
            c = 1.0,
            background_charge = -1.0,
            charge_rtol = 0.0,
        )
        _pencil_fill_mode!(model.E[2], i -> cospi(i[1] / 2) / 64)
        _pencil_fill_mode!(model.B[3], i -> 1 / 8 + sinpi(i[1] / 2) / 64)
        @test init_empic!(model, species) === model
        @test species[1] === owner && model.owners[1] === owner
        _pencil_initial_density(model.rho_n, 0.0)
        _pencil_initial_density(model.rho_np1, 0.0)
        @test _pencil_magnitude(model.J, comm) == 0
        @test 1 / 32 < maxwell_cfl(model)
        fields = (model.E..., model.B..., model.J..., model.rho_n, model.rho_np1)
        _pencil_particle_lifecycle(
            h -> step_empic!(model, species, h),
            () -> electromagnetic_diagnostics(model, species),
            model,
            fields,
            species,
            (model._species_scalars, model._population_counts),
            comm,
        )
        @test _pencil_magnitude(model.J, comm) > 0
    end
    @testset "Isothermal hybrid" begin
        g = distributed_fourier_grid((4, 4, 4), (1.0, 1.0, 1.0); comm)
        species = _pencil_charged_species(comm)
        owner = species[1]
        physics = HybridModel(IsothermalElectrons(1 / 8); nfloor = 1 / 16)
        model = distributed_hybrid(g, species, physics)
        f = model.fields
        _pencil_fill_mode!(f.B[3], i -> 1 / 8 + sinpi(i[1] / 2) / 64)
        @test init!(model, species) === model
        @test species[1] === owner && model.owners[1] === owner
        _pencil_initial_density(f.n, 1.0)
        @test parent(f.pe) == parent(f.n) ./ 8
        @test _pencil_magnitude(f.gradp, comm) > 0
        fields = (f.B..., f.E..., f.n, f.ui..., f.pe, f.J..., f.gradp..., f.ninv)
        _pencil_particle_lifecycle(
            h -> step!(model, species, h; NB = 2),
            () -> hybrid_diagnostics(model, species),
            model,
            fields,
            species,
            (model._species_scalars, model._population_counts),
            comm,
        )
        @test all(value -> value >= model.nfloor, parent(f.n))
    end
end

@testset "Actual-pencil vacuum Maxwell mode and continuation" begin
    comm = MPI.COMM_WORLD
    g = distributed_fourier_grid((4, 4, 4), (1.0, 1.0, 1.0); comm)
    model = distributed_maxwell_vacuum(g; c = 1.0)
    e, b, h, k = 1 / 16, 1 / 32, 1 / 32, 2pi
    _pencil_fill_mode!(model.fields.E[2], i -> e * cospi(i[1] / 2))
    _pencil_fill_mode!(model.fields.B[3], i -> b * sinpi(i[1] / 2))
    fields, phase = (model.fields.E..., model.fields.B...), (model.time, model.step)
    advance(dt) = step_maxwell!(model, dt)
    _pencil_reject_step(advance, fields, (), phase, (), comm)
    @test h < maxwell_cfl(model)
    @test advance(h) === model
    # For E_y=e*cos(kx), B_z=b*sin(kx), curls close on these two
    # amplitudes. Allow FFT/three-stage normwise roundoff, not time error.
    bhalf = b + h * k * e / 2
    enew = e - h * k * bhalf
    bnew = bhalf + h * k * enew / 2
    tolerance = 128 * eps(Float64) * (abs(e) + abs(b) + abs(h * k * e))
    E = map(PencilArrays.global_view, model.fields.E)
    B = map(PencilArrays.global_view, model.fields.B)
    for I in CartesianIndices(E[1])
        @test abs(E[2][I] - enew * cospi((I[1] - 1) / 2)) <= tolerance
        @test abs(B[3][I] - bnew * sinpi((I[1] - 1) / 2)) <= tolerance
        @test abs(E[1][I]) <= tolerance && abs(E[3][I]) <= tolerance
        @test abs(B[1][I]) <= tolerance && abs(B[2][I]) <= tolerance
    end
    diagnostics = distributed_maxwell_diagnostics(model)
    @test diagnostics.time == h && diagnostics.step == 1
    @test isfinite(diagnostics.energy) && diagnostics.energy > 0
    @test diagnostics.max_divE <= tolerance && diagnostics.max_divB <= tolerance
    _pencil_checkpoint_continuation(advance, model, fields, (), phase, (), comm)
end
