using HybridPlasmaPIC, Test, Serialization
import MPI

function clock_snapshot(st, ps)
    # Include particle buffers and every solver workspace; omit the FFT plan,
    # whose foreign pointers are unrelated to the simulation clock.
    fields = (; (name => getfield(st, name) for name in fieldnames(typeof(st)) if name != :g)...)
    io = IOBuffer()
    serialize(io, (fields, ps))
    return take!(io)
end

function clock_case(family, D, T; n_sub = 1, populated = false)
    g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> T(2π), D))
    if family === HallMHDState
        st = HallMHDState(g, HallMHDModel(IsothermalElectrons(zero(T))))
        st.fields.n .= one(T)
        if populated
            for I in CartesianIndices(g.n)
                st.fields.n[I] += T(0.1) * cos(T(I[1] - 1) * g.dx[1])
            end
            st.fields.ui[1] .= T(0.2)
        end
        return st, nothing, dt -> step_hall_mhd!(st, dt)
    end
    electron = family === EMPIC1D || family === EMPIC
    np = populated ? prod(g.n) : 0
    ps = ParticleSet{D,T}(np; q = electron ? -one(T) : one(T))
    if populated
        load_lattice!(ps, ntuple(_ -> zero(T), D), g.L, g.n)
        set_density_weight!(ps, one(T), g)
        ps.v[1] .= T(0.2)
    end
    if electron
        st = family(g, np; n0 = populated ? one(T) : zero(T), n_sub)
        init_empic!(st, ps)
        return st, ps, dt -> step_empic!(st, ps, dt)
    end
    st = family(g, HybridModel(IsothermalElectrons(zero(T))), CIC(), np)
    initialize, advance = family === HybridStepper ? (init!, step!) : (init_camcl!, step_camcl!)
    initialize(st, ps)
    return st, ps, dt -> advance(st, ps, dt; NB = 2)
end

@testset "Clock rejection preserves moving populations and evolving fields" begin
    for T in (Float32, Float64),
        (family, D) in
        ((HybridStepper, 2), (CAMCLStepper, 2), (HallMHDState, 2), (EMPIC1D, 1), (EMPIC, 2))

        st, ps, advance = clock_case(family, D, T; populated = true)
        st.step[] = typemax(Int)
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(T(0.01))
        @test clock_snapshot(st, ps) == before
        st.step[] = 0
        st.time[] = floatmax(T)
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(floatmax(T))
        @test clock_snapshot(st, ps) == before
        st.time[] = zero(T)
        before = clock_snapshot(st, ps)
        advance(T(0.01))
        @test clock_snapshot(st, ps) != before
        @test st.time[] == T(0.01)
    end
end

@testset "Solver clocks reject overflow before state mutation" begin
    for T in (Float32, Float64),
        D = 1:3,
        family in (HybridStepper, CAMCLStepper, HallMHDState, EMPIC)

        st, ps, advance = clock_case(family, D, T)
        @test advance(floatmax(T)) === st
        @test st.time[] == floatmax(T)
        @test st.step[] == 1
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(floatmax(T))
        @test clock_snapshot(st, ps) == before
        @test advance(zero(T)) === st
        @test clock_snapshot(st, ps) == before

        st.time[] = zero(T)
        st.step[] = typemax(Int) - 1
        @test advance(T(0.01)) === st
        @test st.step[] == typemax(Int)
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(T(0.01))
        @test clock_snapshot(st, ps) == before
        @test advance(zero(T)) === st
        @test clock_snapshot(st, ps) == before

        for invalid in (T(Inf), T(NaN), -one(T))
            st.time[] = invalid
            st.step[] = 1
            before = clock_snapshot(st, ps)
            @test_throws ArgumentError advance(T(0.01))
            @test clock_snapshot(st, ps) == before
        end
        st.time[] = zero(T)
        st.step[] = -1
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(T(0.01))
        @test clock_snapshot(st, ps) == before
    end
end

@testset "EM clocks preflight the actual sequence of substep additions" begin
    for T in (Float32, Float64), family in (EMPIC1D, EMPIC)
        st, ps, advance = clock_case(family, 1, T)
        advance(floatmax(T))
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(floatmax(T))
        @test clock_snapshot(st, ps) == before

        st, ps, advance = clock_case(family, 1, T; n_sub = 3)
        st.time[] = prevfloat(floatmax(T), 2)
        # Three separately rounded increments overflow, even though the single
        # full-step addition rounds to floatmax(T). Checking only time+dt misses it.
        dt = T(2.4) * eps(floatmax(T))
        @test isfinite(st.time[] + dt)
        before = clock_snapshot(st, ps)
        @test_throws ArgumentError advance(dt)
        @test clock_snapshot(st, ps) == before

        for ns in (1, 3, 7)
            st, ps, advance = clock_case(family, 1, T; n_sub = ns)
            expected = zero(T)
            for dt in T.((0.1, 0.03, 0.27, 0.005))
                # The public clock retains the established fast-step rounding.
                for _ = 1:ns
                    expected += dt / ns
                end
                advance(dt)
                @test st.time[] === expected
            end
            @test st.step[] == 4
            st.step[] = typemax(Int)
            before = clock_snapshot(st, ps)
            @test_throws ArgumentError advance(T(0.01))
            @test clock_snapshot(st, ps) == before
        end
    end
end

@testset "MPI clock rejects invalid advancement before collectives or mutation" begin
    ensure_mpi_initialized!()
    ctx = create_cartesian_communicator(LogicalRankLayout((1,)); comm = MPI.COMM_SELF)
    status = GPUAwareMPIStatus(false, false, false, :test, "CPU clock test")
    try
        for T in (Float32, Float64)
            st, ps, _ = clock_case(HybridStepper, 1, T)
            mpi_init!(st, ps, ctx; gpu_status = status)
            advance(dt) = mpi_step!(st, ps, ctx, dt; NB = 2, gpu_status = status)
            advance(floatmax(T))
            before = clock_snapshot(st, ps)
            @test_throws ArgumentError advance(floatmax(T))
            @test clock_snapshot(st, ps) == before
            st.time[] = zero(T)
            st.step[] = typemax(Int)
            before = clock_snapshot(st, ps)
            @test_throws ArgumentError advance(T(0.01))
            @test clock_snapshot(st, ps) == before
        end
    finally
        free_mpi_communicator!(ctx)
    end
end
