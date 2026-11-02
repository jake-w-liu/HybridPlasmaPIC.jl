using Test
using HybridPlasmaPIC

function _implicit_yee_em_snapshot(state, species)
    saved = HybridPlasmaPIC._checkpoint_state(state, species)
    report = saved.last_report
    values =
        report === nothing ? nothing :
        NamedTuple{fieldnames(typeof(report))}(
            ntuple(i -> getfield(report, i), fieldcount(typeof(report))),
        )
    return merge(saved, (last_report = values,))
end

# Neutral two-species lattice: 4n co-located particles per species, electrons
# sinusoidally perturbed around the uniform ion lattice.
function _implicit_yee_plasma(::Type{T}, n, amp; ion_mass = T(100)) where {T}
    L, np = T(n), 4n
    xs = [(T(i) - T(1) / 2) * L / np for i = 1:np]
    electrons = ParticleSet{1,T}(np; q = -one(T), m = one(T))
    ions = ParticleSet{1,T}(np; q = one(T), m = ion_mass)
    electrons.x[1] .= mod.(xs .+ amp .* sin.(2 * T(pi) .* xs ./ L), L)
    ions.x[1] .= xs
    return [electrons, ions]
end

@testset "Compatible periodic Yee implicit EM" begin
    Q = HybridPlasmaPIC._exact_rational
    for T in (Float32, Float64)
        gate = Q(256eps(T))
        @testset "analytic uniform DC step $T" begin
            g = YeeGrid((16,), (T(16),); boundaries = (:periodic,))
            np = 64
            xs = [(T(i) - T(1) / 2) / T(4) for i = 1:np]
            species = [
                ParticleSet{1,T}(np; q = -one(T), m = one(T)),
                ParticleSet{1,T}(np; q = one(T), m = one(T)),
            ]
            for ps in species
                ps.x[1] .= xs
            end
            st = ImplicitYeeEMPICStepper(g, species; shape = YeeTriangularCloud(one(T)), c = T(10))
            st.fields.E[1] .= T(1 / 8)
            init!(st, species)
            step!(st, species, T(1 / 8))
            atol = 32 * (st.last_report.tolerance + eps(T))
            # omega^2 = 128/16 = 8, omega^2*dt^2/4 = 1/32: E1 = E0(1-1/32)/(1+1/32),
            # v1 = (q/m) dt E0/(1+1/32).
            @test isapprox(
                st.fields.E[1],
                fill(T(1 / 8) * T(31 / 33), 16);
                atol = atol,
                rtol = zero(T),
            )
            for ps in species
                @test isapprox(ps.v[1], fill(ps.q / 66, np); atol = atol, rtol = zero(T))
            end
        end
        @testset "nonuniform published Ex $T" begin
            for n in (16, 15), (amp, dt) in ((T(0.05), T(0.1)), (T(0.3), T(0.1)), (T(0.3), T(0.5)))

                g = YeeGrid((n,), (T(n),); boundaries = (:periodic,))
                species = _implicit_yee_plasma(T, n, amp)
                st = ImplicitYeeEMPICStepper(
                    g,
                    species;
                    shape = YeeTriangularCloud(one(T)),
                    c = T(10),
                )
                init!(st, species)
                for k = 1:3
                    step!(st, species, dt)
                    r = st.last_report
                    @test r.gauss <= gate
                    @test r.continuity <= gate
                    @test r.energy_change == sum(r.defects)
                    @test iszero(r.defects[2])
                    @test iszero(r.curl_adjoint_work)
                    @test r.defects[1] == r.particle_residual_work
                    @test r.defects[4] == r.field_residual_work + r.curl_adjoint_work
                    E1 = st.fields.E[1]
                    @test maximum(
                        i -> abs(Q(E1[i]) - Q(E1[i == 1 ? n : i - 1]) - st.charge_np1[i]),
                        1:n,
                    ) / r.charge_scale <= gate
                end
            end
        end
        @testset "persistence $T" begin
            g = YeeGrid((16,), (T(16),); boundaries = (:periodic,))
            species = _implicit_yee_plasma(T, 16, T(0.3))
            st = ImplicitYeeEMPICStepper(g, species; shape = YeeTriangularCloud(one(T)), c = T(10))
            init!(st, species)
            step!(st, species, T(0.1))
            step!(st, species, T(0.1))
            mktempdir() do directory
                checkpoint = joinpath(directory, "pyee.chk")
                archive = joinpath(directory, "pyee.run")
                save_checkpoint(checkpoint, st, species)
                saved = _implicit_yee_em_snapshot(st, species)
                step!(st, species, T(0.1))
                load_checkpoint!(st, species, checkpoint)
                @test isequal(_implicit_yee_em_snapshot(st, species), saved)
                archive_run(archive, st, species; rng_seed = 0)
                step!(st, species, T(0.1))
                load_archive!(st, species, archive)
                @test isequal(_implicit_yee_em_snapshot(st, species), saved)
            end
        end
        @testset "failure preservation $T" begin
            g = YeeGrid((16,), (T(16),); boundaries = (:periodic,))
            species = _implicit_yee_plasma(T, 16, T(0.3))
            failed = ImplicitYeeEMPICStepper(
                g,
                species;
                shape = YeeTriangularCloud(one(T)),
                c = T(10),
                max_newton = 0,
            )
            init!(failed, species)
            before = _implicit_yee_em_snapshot(failed, species)
            @test_throws ImplicitEMConvergenceError step!(failed, species, T(0.1))
            @test isequal(_implicit_yee_em_snapshot(failed, species), before)
            @test step!(failed, species, zero(T)) === failed
            @test step!(failed, species, -zero(T)) === failed
            @test isequal(_implicit_yee_em_snapshot(failed, species), before)
            @test_throws ArgumentError step!(failed, species, -one(T))
            @test isequal(_implicit_yee_em_snapshot(failed, species), before)
        end
        @testset "constructor rejections $T" begin
            grid() = YeeGrid((16,), (T(16),); boundaries = (:periodic,))
            shape1 = YeeTriangularCloud(one(T))
            species = _implicit_yee_plasma(T, 16, zero(T))
            @test_throws ArgumentError ImplicitYeeEMPICStepper(
                YeeGrid((16,), (T(16),)),
                species;
                shape = shape1,
                c = T(10),
            )
            @test_throws ArgumentError ImplicitYeeEMPICStepper(
                YeeGrid((3,), (one(T),); boundaries = (:periodic,)),
                species;
                shape = shape1,
                c = T(10),
            )
            @test_throws ArgumentError ImplicitYeeEMPICStepper(
                grid(),
                species;
                shape = YeeTriangularCloud(T(9)),
                c = T(10),
            )
            electrons = [ParticleSet{1,T}(4; q = -one(T), m = one(T))]
            electrons[1].x[1] .= T.(1:4)
            @test_throws ArgumentError ImplicitYeeEMPICStepper(
                grid(),
                electrons;
                shape = shape1,
                c = T(10),
            )
            boundary_species = _implicit_yee_plasma(T, 16, zero(T))
            boundary_species[1].x[1][1] = T(16)
            @test_throws ArgumentError ImplicitYeeEMPICStepper(
                grid(),
                boundary_species;
                shape = shape1,
                c = T(10),
            )
        end
        @testset "diagnostics $T" begin
            g = YeeGrid((16,), (T(16),); boundaries = (:periodic,))
            species = _implicit_yee_plasma(T, 16, T(0.3))
            st = ImplicitYeeEMPICStepper(g, species; shape = YeeTriangularCloud(one(T)), c = T(10))
            init!(st, species)
            step!(st, species, T(0.1))
            diag = implicit_yee_em_diagnostics(st, species)
            @test diag.total isa T
            @test diag.gauss <= gate
            @test diag.current_time_offset == -T(0.1) / 2
        end
    end
end
