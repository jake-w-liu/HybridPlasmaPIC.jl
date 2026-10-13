using HybridPlasmaPIC, Test
import MPI

function mpi_ion_fixture(::Type{T}, ::Val{D}, shape, closure; de2 = zero(T)) where {T,D}
    g = FourierGrid(ntuple(d -> d == 1 ? 8 : 4, D), ntuple(_ -> T(2π), D))
    species = [ParticleSet{D,T}(s == 3 ? 0 : s * prod(g.n); q = T(s), m = T(s^2)) for s = 1:3]
    for (s, ps) in enumerate(species)
        s == 3 && continue
        p = 0
        for I in CartesianIndices(g.n), j = 1:s
            p += 1
            for d = 1:D
                fraction = d == 1 ? (T(j) - T(0.5)) / T(s) : T(0.5)
                ps.x[d][p] = (T(I[d] - 1) + fraction) * g.dx[d]
            end
            phase = ps.x[1][p]
            ps.weight[p] = T(0.5 / s^2) * prod(g.dx) * (one(T) + T(0.02) * cos(phase))
            for c = 1:3
                ps.v[c][p] = T(0.015 * (s - c)) + T(0.002) * sin(T(c) * phase)
            end
            ps.id[p] = UInt64(p)
            ps.tag[p] = UInt32(100s + p)
        end
    end
    st = HybridStepper(g, HybridModel(closure; nfloor = T(0.01), de2), shape, 0)
    fill!(st.fields.B[1], T(0.8))
    for I in CartesianIndices(g.n)
        phase = T(I[1] - 1) * g.dx[1]
        st.fields.B[2][I] = T(0.002) * cos(phase)
        st.fields.B[3][I] = T(0.001) * sin(phase)
    end
    return st, species
end

function mpi_ion_subset(ps::ParticleSet{D,T}, rank, size) where {D,T}
    indices = [p for p in eachindex(ps.weight) if mod(p - 1, size) == rank]
    return ParticleSet{D,T}(
        map(a -> a[indices], ps.x),
        map(a -> a[indices], ps.v),
        ps.weight[indices],
        ps.id[indices],
        ps.tag[indices],
        ps.q,
        ps.m,
    )
end

function mpi_ion_state(st, species)
    fields = map(fieldnames(typeof(st.fields))) do name
        value = getfield(st.fields, name)
        value isa Ref ? value[] : value
    end
    particles = [(ps.x, ps.v, ps.weight, ps.id, ps.tag, ps.q, ps.m) for ps in species]
    return deepcopy((
        fields,
        particles,
        st.time[],
        st.step[],
        st.last_dt[],
        st.Ep,
        st.Bp,
        st.xmid,
        st.vpred,
        st.work,
    ))
end

function mpi_ion_compare(st, species, reference, full, ctx; migrated = false)
    T = eltype(st.fields.n)
    # Collective summation changes arithmetic order; this bound is small
    # compared with the imposed density/velocity/field perturbations.
    tolerance = T(256) * eps(T)
    for name in (:n, :ui, :E, :B, :J, :pe, :gradp, :ninv)
        actual, expected = getproperty(st.fields, name), getproperty(reference.fields, name)
        aa = actual isa Tuple ? actual : (actual,)
        ee = expected isa Tuple ? expected : (expected,)
        @test all(isapprox(a, e; rtol = tolerance, atol = tolerance) for (a, e) in zip(aa, ee))
    end
    @test (st.time[], st.step[], st.last_dt[]) ==
          (reference.time[], reference.step[], reference.last_dt[])
    for (ps, original) in zip(species, full)
        @test (ps.q, ps.m) == (original.q, original.m)
        @test MPI.Allreduce(nparticles(ps), +, ctx.comm) == nparticles(original)
        @test MPI.Allreduce(sum(ps.id), +, ctx.comm) == sum(original.id)
        if !isempty(original.id)
            occurrences = zeros(Int, nparticles(original))
            for id in ps.id
                occurrences[Int(id)] += 1
            end
            @test MPI.Allreduce(occurrences, +, ctx.comm) == ones(Int, nparticles(original))
        end
        for (p, id) in enumerate(ps.id)
            index = Int(id)
            @test ps.weight[p] == original.weight[index]
            @test ps.tag[p] == original.tag[index]
            @test all(
                isapprox(ps.x[d][p], original.x[d][index]; rtol = tolerance, atol = tolerance) for
                d in eachindex(ps.x)
            )
            @test all(
                isapprox(ps.v[c][p], original.v[c][index]; rtol = tolerance, atol = tolerance) for
                c = 1:3
            )
            if migrated
                position = map(a -> a[p], ps.x)
                @test rank_of_position(position, st.g, ctx.layout) == ctx.logical_rank
            end
        end
    end
end

@testset "MPI ion collections" begin
    ensure_mpi_initialized!()
    size, rank = MPI.Comm_size(MPI.COMM_WORLD), MPI.Comm_rank(MPI.COMM_WORLD)
    status = GPUAwareMPIStatus(false, false, false, :test, "CPU collection regression")
    for T in (Float32, Float64), D = 1:3
        layout =
            LogicalRankLayout(ntuple(d -> d == 1 ? size : 1, D); periodic = ntuple(_ -> true, D))
        ctx = create_cartesian_communicator(layout)
        try
            @testset "serial agreement $T $D dimensions" begin
                configurations = (
                    (NGP(), IsothermalElectrons(T(0.03)), zero(T)),
                    (CIC(), IsothermalElectrons(T(0.03)), zero(T)),
                    (TSC(), IsothermalElectrons(T(0.03)), T(0.1)),
                    (CIC(), PolytropicElectrons(T(0.03), one(T), T(5 / 3)), zero(T)),
                    (CIC(), CGLElectrons(T(0.02), T(0.03), one(T), one(T)), T(0.1)),
                )
                for (shape, closure, de2) in configurations
                    reference, full = mpi_ion_fixture(T, Val(D), shape, closure; de2)
                    st = deepcopy(reference)
                    species = [mpi_ion_subset(ps, rank, size) for ps in full]
                    init!(reference, full)
                    mpi_init!(st, species, ctx; gpu_status = status)
                    mpi_ion_compare(st, species, reference, full, ctx)
                    for (i, dt) in enumerate(T.((0.001, 0.0007, 0.0013)))
                        step!(reference, full, dt; NB = i)
                        mpi_step!(
                            st,
                            species,
                            ctx,
                            dt;
                            NB = i,
                            gpu_status = status,
                            migrate_particles = isodd(i),
                        )
                        mpi_ion_compare(st, species, reference, full, ctx; migrated = isodd(i))
                    end
                end
            end

            @testset "independent global moments $T $D dimensions" begin
                st, _ = mpi_ion_fixture(T, Val(D), NGP(), IsothermalElectrons(zero(T)))
                species = [ParticleSet{D,T}(1; q = T(s), m = T(s^2)) for s = 1:2]
                for (s, ps) in enumerate(species)
                    for d = 1:D
                        ps.x[d][1] = T(0.1) * st.g.dx[d]
                    end
                    ps.weight[1] = T(0.125 * s / size) * prod(st.g.dx)
                    for c = 1:3
                        ps.v[c][1] = T((rank + 1) * s * c / 16)
                    end
                end
                for nf in T.((0.5, 0.875)), scratch in (nothing, zeros(T, 2))
                    mpi_compute_moments!(
                        st.fields,
                        species,
                        st.g,
                        NGP(),
                        nf,
                        ctx;
                        work = scratch,
                        gpu_status = status,
                    )
                    expected_n = zeros(T, st.g.n)
                    expected_n[1] = T(5 / 8)
                    @test isapprox(st.fields.n, expected_n; rtol = 32eps(T), atol = 32eps(T))
                    for c = 1:3
                        expected = zeros(T, st.g.n)
                        # Σ_rank Σ_s q_s n_s v_s / max(Σ q_s n_s, floor).
                        expected[1] = T(9c * (size + 1) / 256) / max(T(5 / 8), nf)
                        @test isapprox(st.fields.ui[c], expected; rtol = 32eps(T), atol = 32eps(T))
                    end
                end
                # An entire rank may have no ions; it still joins every reduction.
                if rank != 0
                    species = [ParticleSet{D,T}(0; q = T(s), m = T(s^2)) for s = 1:2]
                end
                mpi_init!(st, species, ctx; gpu_status = status)
                mpi_step!(st, species, ctx, T(0.0001); gpu_status = status)
                @test MPI.Allreduce(sum(nparticles, species), +, ctx.comm) == 2
                for empty in (ParticleSet{D,T}[], [ParticleSet{D,T}(0; q = T(2), m = T(4))])
                    mpi_init!(st, empty, ctx; gpu_status = status)
                    mpi_step!(st, empty, ctx, T(0.001); gpu_status = status)
                    @test all(iszero, st.fields.n)
                    @test all(a -> all(iszero, a), st.fields.ui)
                    @test st.step[] == 1
                end
            end
        finally
            free_mpi_communicator!(ctx)
        end
    end

    ctx = create_cartesian_communicator(LogicalRankLayout((size,); periodic = (true,)))
    try
        @testset "collection validation is collective and precedes mutation" begin
            # Model parameters are Float64; validate the floor after conversion
            # to particle/grid precision, before any output or workspace changes.
            for floor in (1e-100, 1e100), collection in (false, true), operation in (:init, :step)
                st, full = mpi_ion_fixture(Float32, Val(1), CIC(), IsothermalElectrons(0.03f0))
                species = [mpi_ion_subset(ps, rank, size) for ps in full]
                populations = collection ? species : species[1]
                mpi_init!(st, populations, ctx; gpu_status = status)
                if !collection || rank == 0
                    st.model = HybridModel(st.model.closure; nfloor = floor)
                end
                before = mpi_ion_state(st, species)
                if operation == :init
                    @test_throws ArgumentError mpi_init!(st, populations, ctx; gpu_status = status)
                else
                    @test_throws ArgumentError mpi_step!(
                        st,
                        populations,
                        ctx,
                        0.001f0;
                        gpu_status = status,
                    )
                end
                @test isequal(before, mpi_ion_state(st, species))
            end
            for defect in (
                :charge,
                :position,
                :velocity,
                :weight,
                :length,
                :alias,
                :clock,
                :counter,
                :previous_dt,
            )
                st, full = mpi_ion_fixture(Float64, Val(1), CIC(), IsothermalElectrons(0.03))
                species = [mpi_ion_subset(ps, rank, size) for ps in full]
                mpi_init!(st, species, ctx; gpu_status = status)
                if rank == size - 1
                    ps = species[2]
                    defect == :charge && (ps.q = -1.0)
                    defect == :position && (ps.x[1][1] = Inf)
                    defect == :velocity && (ps.v[3][1] = NaN)
                    defect == :weight && (ps.weight[1] = -1.0)
                    defect == :length && pop!(ps.id)
                    defect == :alias && (species[2] = species[1])
                    defect == :clock && (st.time[] = Inf)
                    defect == :counter && (st.step[] = typemax(Int))
                    defect == :previous_dt && (st.last_dt[] = -1.0)
                end
                before = mpi_ion_state(st, species)
                @test_throws Union{ArgumentError,DimensionMismatch} mpi_step!(
                    st,
                    species,
                    ctx,
                    0.001;
                    gpu_status = status,
                )
                @test isequal(before, mpi_ion_state(st, species))
            end
            st, species = mpi_ion_fixture(Float64, Val(1), CIC(), IsothermalElectrons(0.03))
            species = [mpi_ion_subset(ps, rank, size) for ps in species]
            for bad in (-0.1, Inf, NaN)
                before = mpi_ion_state(st, species)
                @test_throws ArgumentError mpi_step!(
                    st,
                    species,
                    ctx,
                    rank == 0 ? bad : 0.001;
                    gpu_status = status,
                )
                @test isequal(before, mpi_ion_state(st, species))
            end
            before = mpi_ion_state(st, species)
            @test_throws ArgumentError mpi_step!(
                st,
                species,
                ctx,
                0.001;
                NB = rank == 0 ? 0 : 1,
                gpu_status = status,
            )
            @test isequal(before, mpi_ion_state(st, species))
            # All ranks agree on a zero step; malformed particles remain untouched.
            species[2].q = -1.0
            before = mpi_ion_state(st, species)
            mpi_step!(st, species, ctx, 0.0; gpu_status = status)
            @test isequal(before, mpi_ion_state(st, species))
            @test_throws ArgumentError mpi_init!(st, species, ctx; gpu_status = status)
            @test isequal(before, mpi_ion_state(st, species))
            species[2].q = 2.0
            before = mpi_ion_state(st, species)
            @test_throws Union{ArgumentError,DimensionMismatch} mpi_compute_moments!(
                st.fields,
                species,
                st.g,
                CIC(),
                0.1,
                ctx;
                work = zeros(rank == 0 ? 1 : sum(nparticles, species)),
                gpu_status = status,
            )
            @test isequal(before, mpi_ion_state(st, species))

            for floor in (0.0, Inf, NaN)
                before = mpi_ion_state(st, species)
                @test_throws ArgumentError mpi_compute_moments!(
                    st.fields,
                    species,
                    st.g,
                    CIC(),
                    rank == 0 ? floor : 0.1,
                    ctx;
                    gpu_status = status,
                )
                @test isequal(before, mpi_ion_state(st, species))
            end

            for alias in (:particle_work, :output_work, :outputs, :output_particle)
                st, full = mpi_ion_fixture(Float64, Val(1), CIC(), IsothermalElectrons(0.03))
                # Keep eight particles per rank so output and work extents match.
                species = [full[1]]
                work = zeros(8)
                if rank == 0
                    alias == :particle_work && (work = species[1].weight)
                    alias == :output_work && (work = st.fields.n)
                    alias == :outputs &&
                        (st.fields.ui = (st.fields.n, st.fields.ui[2], st.fields.ui[3]))
                    alias == :output_particle && (st.fields.n = species[1].weight)
                end
                before = mpi_ion_state(st, species)
                @test_throws ArgumentError mpi_compute_moments!(
                    st.fields,
                    species,
                    st.g,
                    CIC(),
                    0.1,
                    ctx;
                    work,
                    gpu_status = status,
                )
                @test isequal(before, mpi_ion_state(st, species))
            end

            @testset "views retain the nonmigrating evolution path" begin
                st, full = mpi_ion_fixture(Float64, Val(1), CIC(), IsothermalElectrons(0.03))
                reference = deepcopy(st)
                species = [mpi_ion_subset(ps, rank, size) for ps in full]
                viewed = map(species) do ps
                    ParticleSet{1,Float64}(
                        map(a -> view(a, :), ps.x),
                        map(a -> view(a, :), ps.v),
                        view(ps.weight, :),
                        view(ps.id, :),
                        view(ps.tag, :),
                        ps.q,
                        ps.m,
                    )
                end
                init!(reference, full)
                mpi_init!(st, viewed, ctx; gpu_status = status)
                before = mpi_ion_state(st, viewed)
                @test_throws ArgumentError mpi_step!(st, viewed, ctx, 0.001; gpu_status = status)
                @test isequal(before, mpi_ion_state(st, viewed))
                step!(reference, full, 0.001)
                mpi_step!(st, viewed, ctx, 0.001; gpu_status = status, migrate_particles = false)
                mpi_ion_compare(st, viewed, reference, full, ctx)
            end

            if size > 1
                for mismatch in
                    (:count, :order, :mass, :step, :time, :last_dt, :dt, :zero_dt, :NB, :migration)
                    st, full = mpi_ion_fixture(Float64, Val(1), CIC(), IsothermalElectrons(0.03))
                    species = [mpi_ion_subset(ps, rank, size) for ps in full]
                    mpi_init!(st, species, ctx; gpu_status = status)
                    if rank == 0
                        mismatch == :count && pop!(species)
                        mismatch == :order && reverse!(species)
                        mismatch == :mass && (species[2].m = 5.0)
                        mismatch == :step && (st.step[] = 1)
                        mismatch == :time && (st.time[] = 0.1)
                        mismatch == :last_dt && (st.last_dt[] = 0.01)
                    end
                    dt =
                        rank == 0 && mismatch == :dt ? 0.002 :
                        rank == 0 && mismatch == :zero_dt ? 0.0 : 0.001
                    NB = rank == 0 && mismatch == :NB ? 2 : 1
                    migrate = !(rank == 0 && mismatch == :migration)
                    before = mpi_ion_state(st, species)
                    @test_throws ArgumentError mpi_step!(
                        st,
                        species,
                        ctx,
                        dt;
                        NB,
                        gpu_status = status,
                        migrate_particles = migrate,
                    )
                    @test isequal(before, mpi_ion_state(st, species))
                end
            end
        end
    finally
        free_mpi_communicator!(ctx)
    end
end
