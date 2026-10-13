using HybridPlasmaPIC, Test, Random
import MPI

@testset "Hybrid steppers preserve half-time levels across timestep changes" begin
    for (ctor, initialize!, advance!) in
        ((HybridStepper, init!, step!), (CAMCLStepper, init_camcl!, step_camcl!)),
        D = 1:3

        g = FourierGrid(ntuple(_ -> 4, D), ntuple(_ -> 10.0, D))
        ps = ParticleSet{D,Float64}(1)
        ps.weight .= 0.0
        for d = 1:D
            ps.x[d] .= 1.0
        end
        st = ctor(g, HybridModel(IsothermalElectrons(0.0)), CIC(), 1)
        initialize!(st, ps)
        acceleration = (0.7, -0.3, 0.2)
        elapsed = 0.0
        for dt in (0.1, 0.2, 0.05, 0.08)
            previous_dt = st.last_dt[]
            advance!(st, ps, 0.0; NB = 2)
            @test st.last_dt[] == previous_dt
            for c = 1:3
                fill!(st.fields.E[c], acceleration[c])
            end
            advance!(st, ps, dt; NB = 2)
            elapsed += dt
            for d = 1:D
                @test ps.x[d][1] ≈ 1 + acceleration[d] * elapsed^2 / 2 atol = 2e-15
            end
            for c = 1:3
                @test ps.v[c][1] ≈ acceleration[c] * (elapsed - dt / 2) atol = 2e-15
            end
            @test st.time[] ≈ elapsed
            @test st.last_dt[] == dt
        end
        state = (deepcopy(ps.x), deepcopy(ps.v), st.last_dt[], st.time[], st.step[])
        for dt in (-0.1, NaN, Inf)
            @test_throws ArgumentError advance!(st, ps, dt; NB = 2)
            @test (ps.x, ps.v, st.last_dt[], st.time[], st.step[]) == state
        end
        foreach(v -> fill!(v, 0.0), ps.v)
        initialize!(st, ps)
        @test st.last_dt[] == 0.0
        for c = 1:3
            fill!(st.fields.E[c], acceleration[c])
        end
        advance!(st, ps, 0.07; NB = 2)
        for d = 1:D
            @test ps.x[d][1] ≈ state[1][d][1] + acceleration[d] * 0.07^2 / 2 atol = 2e-15
        end
    end
end

@testset "Variable-step hybrid gyromotion converges to analytic orbit" begin
    for (ctor, initialize!, advance!) in
        ((HybridStepper, init!, step!), (CAMCLStepper, init_camcl!, step_camcl!))
        function orbit_error(refinement)
            g = FourierGrid((4, 4), (10.0, 10.0))
            ps = ParticleSet{2,Float64}(1)
            ps.x[1] .= 1.0
            ps.x[2] .= 1.0
            ps.v[1] .= 0.3
            ps.v[2] .= -0.2
            ps.weight .= 0.0
            st = ctor(g, HybridModel(IsothermalElectrons(0.0)), CIC(), 1)
            fill!(st.fields.B[3], 1.0)
            initialize!(st, ps)
            for _ = 1:(10refinement), dt in (0.04 / refinement, 0.02 / refinement)
                advance!(st, ps, dt; NB = 2)
            end
            t = 0.6
            exact = (1 + 0.3sin(t) - 0.2(1 - cos(t)), 1 - 0.2sin(t) - 0.3(1 - cos(t)))
            return hypot(ps.x[1][1] - exact[1], ps.x[2][1] - exact[2])
        end
        errors = orbit_error.((1, 2, 4))
        @test log2(errors[1] / errors[2]) > 1.8
        @test log2(errors[2] / errors[3]) > 1.8
    end
end

@testset "MPI variable steps and checkpoint preserve serial trajectory" begin
    ensure_mpi_initialized!()
    ctx = create_cartesian_communicator(LogicalRankLayout((1,)); comm = MPI.COMM_SELF)
    try
        g = FourierGrid((8,), (2π,))
        ps = ParticleSet{1,Float64}(128)
        load_lattice_1d!(ps, 0.0, 2π)
        set_density_weight!(ps, 1.0, g)
        @. ps.v[1] = 0.01sin(ps.x[1])
        pm = deepcopy(ps)
        st = HybridStepper(g, HybridModel(IsothermalElectrons(0.3)), CIC(), nparticles(ps))
        sm = HybridStepper(g, st.model, CIC(), nparticles(pm))
        init!(st, ps)
        mpi_init!(sm, pm, ctx)
        for dt in (0.01, 0.005, 0.007, 0.0)
            step!(st, ps, dt; NB = 2)
            mpi_step!(sm, pm, ctx, dt; NB = 2)
            @test pm.x == ps.x
            @test pm.v == ps.v
            @test sm.fields.E == st.fields.E
            @test sm.fields.B == st.fields.B
            @test sm.last_dt[] == st.last_dt[]
        end
        mktempdir() do dir
            save_mpi_checkpoint(dir, sm, pm, ctx)
            restored = HybridStepper(g, st.model, CIC(), 0)
            pr = ParticleSet{1,Float64}(0)
            load_mpi_checkpoint!(restored, pr, dir, ctx)
            @test restored.last_dt[] == sm.last_dt[]
            for dt in (0.003, 0.006)
                mpi_step!(restored, pr, ctx, dt; NB = 2)
                step!(st, ps, dt; NB = 2)
            end
            @test pr.x == ps.x
            @test pr.v == ps.v
            @test restored.fields.E == st.fields.E
        end
        mpi_init!(sm, pm, ctx)
        @test sm.last_dt[] == 0.0
    finally
        free_mpi_communicator!(ctx)
    end
end

@testset "Variable-step coupled hybrid field converges at second order" begin
    for (ctor, initialize!, advance!) in
        ((HybridStepper, init!, step!), (CAMCLStepper, init_camcl!, step_camcl!))
        function field_at_end(refinement)
            g = FourierGrid((8,), (2π,))
            ps = ParticleSet{1,Float64}(1024)
            load_lattice_1d!(ps, 0.0, 2π)
            set_density_weight!(ps, 1.0, g)
            @. ps.v[1] = 0.01sin(ps.x[1])
            st = ctor(g, HybridModel(IsothermalElectrons(0.3)), CIC(), nparticles(ps))
            fill!(st.fields.B[1], 1.0)
            x = [i * g.dx[1] for i = 0:7]
            @. st.fields.B[2] = 0.02cos(x)
            @. st.fields.B[3] = 0.02sin(x)
            initialize!(st, ps)
            for _ = 1:(10refinement), dt in (0.02 / refinement, 0.01 / refinement)
                advance!(st, ps, dt; NB = 2)
            end
            return st.fields.B[2]
        end
        fields = field_at_end.((1, 2, 4))
        coarse = sqrt(sum(abs2, fields[1] .- fields[2]))
        fine = sqrt(sum(abs2, fields[2] .- fields[3]))
        @test log2(coarse / fine) > 1.7
    end
end
