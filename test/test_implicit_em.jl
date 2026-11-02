using Test
using HybridPlasmaPIC

# Check all scalar report values after detachment; custom records containing
# BigInts have no package-defined value equality, so retain their named fields.
function _implicit_em_snapshot(state, species)
    saved = HybridPlasmaPIC._checkpoint_state(state, species)
    report = saved.last_report
    values =
        report === nothing ? nothing :
        NamedTuple{fieldnames(typeof(report))}(
            ntuple(i -> getfield(report, i), fieldcount(typeof(report))),
        )
    return merge(saved, (last_report = values,))
end

@testset "Synchronous coupled implicit EM essentials" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64)
        @testset "$T" begin
            g = FourierGrid((5,), (T(5),))
            populations = typeof(ParticleSet{1,T}(0))[]
            state = ImplicitEMPICStepper(g, CIC(), populations; c = T(2))
            wave = [cospi(T(2) * T(i - 1) / T(5)) for i = 1:5]
            sine = [sinpi(T(2) * T(i - 1) / T(5)) for i = 1:5]
            state.fields.E[2] .= wave
            state.fields.E[3] .= wave ./ T(2)
            init!(state, populations)
            before = _implicit_em_snapshot(state, populations)
            @test step!(state, populations, -zero(T)) === state
            @test isequal(_implicit_em_snapshot(state, populations), before)
            @test_throws ArgumentError step!(state, populations, -one(T))
            @test isequal(_implicit_em_snapshot(state, populations), before)
            h = T(1 / 8)
            step!(state, populations, h)
            # Independent exact-in-time midpoint response of this single
            # represented Fourier mode, for both transverse polarizations.
            k = g.kvec[1][2]
            denominator = one(T) + (state.c * k * h / T(2))^2
            ratio = (one(T) - (state.c * k * h / T(2))^2) / denominator
            allowance = T(32) * (state.last_report.tolerance + eps(T))
            @test isapprox(state.fields.E[2], ratio .* wave; atol = allowance, rtol = zero(T))
            @test isapprox(
                state.fields.E[3],
                ratio .* wave ./ T(2);
                atol = allowance,
                rtol = zero(T),
            )
            @test isapprox(
                state.fields.B[3],
                (h * k / denominator) .* sine;
                atol = allowance,
                rtol = zero(T),
            )
            @test isapprox(
                state.fields.B[2],
                -(h * k / denominator) .* sine ./ T(2);
                atol = allowance,
                rtol = zero(T),
            )
            r = state.last_report
            @test r.energy_change == sum(r.defects)
            bad_report = typeof(r)(
                map(
                    name -> name == :time ? h / T(2) : getproperty(r, name),
                    fieldnames(typeof(r)),
                )...,
            )
            @test_throws ArgumentError H._implicit_em_phase(state, true, h / T(2), 1, h, bad_report)
            mktempdir() do directory
                checkpoint = joinpath(directory, "implicit.chk")
                archive = joinpath(directory, "implicit.run")
                save_checkpoint(checkpoint, state, populations)
                saved = _implicit_em_snapshot(state, populations)
                step!(state, populations, h)
                load_checkpoint!(state, populations, checkpoint)
                @test isequal(_implicit_em_snapshot(state, populations), saved)
                archive_run(archive, state, populations; rng_seed = 0)
                step!(state, populations, h)
                load_archive!(state, populations, archive)
                @test isequal(_implicit_em_snapshot(state, populations), saved)
            end
            # A genuinely coupled neutral two-species plasma, with uniform
            # density and a self-consistent longitudinal E/current response.
            g2 = FourierGrid((5,), (T(5),))
            species = [ParticleSet{1,T}(5; q = T(q), m = one(T)) for q in (-1, 1)]
            for ps in species
                ps.x[1] .= T.(0:4)
            end
            plasma = ImplicitEMPICStepper(g2, CIC(), species; c = one(T))
            fill!(plasma.fields.E[1], T(1 / 8))
            init!(plasma, species)
            step!(plasma, species, h)
            denominator = one(T) + h^2 / T(2)
            electric = T(1 / 8) * (one(T) - h^2 / T(2)) / denominator
            allowance = T(32) * (plasma.last_report.tolerance + eps(T))
            @test isapprox(plasma.fields.E[1], fill(electric, 5); atol = allowance, rtol = zero(T))
            for ps in species
                velocity = ps.q * h * T(1 / 8) / denominator
                @test isapprox(ps.v[1], fill(velocity, 5); atol = allowance, rtol = zero(T))
            end
            @test implicit_em_diagnostics(plasma).total isa T
            # Exhausting the declared solver bound preserves the public phase
            # and leaves the same instance usable for a validated zero call.
            failed =
                ImplicitEMPICStepper(FourierGrid((5,), (T(5),)), CIC(), populations; max_newton = 0)
            failed.fields.E[2] .= wave
            init!(failed, populations)
            before = _implicit_em_snapshot(failed, populations)
            @test_throws ImplicitEMConvergenceError step!(failed, populations, h)
            @test isequal(_implicit_em_snapshot(failed, populations), before)
            @test step!(failed, populations, zero(T)) === failed
        end
    end
end

@testset "Gauss- and divergence-consistent published fields" begin
    for T in (Float32, Float64)
        gate = T(256) * eps(T)
        @testset "nonuniform 1D $T" begin
            for (amp, dt) in ((T(0.05), T(0.1)), (T(0.3), T(0.1)), (T(0.3), T(0.5)))
                g = FourierGrid((15,), (T(15),))
                np, L = 60, T(15)
                xs = [(T(i) - T(1) / 2) * L / np for i = 1:np]
                electrons = ParticleSet{1,T}(np; q = -one(T), m = one(T))
                ions = ParticleSet{1,T}(np; q = one(T), m = T(100))
                electrons.x[1] .= mod.(xs .+ amp .* sin.(2 * T(pi) .* xs ./ L), L)
                ions.x[1] .= xs
                species = [electrons, ions]
                state = ImplicitEMPICStepper(g, CIC(), species; c = T(10))
                init!(state, species)
                for k = 1:3
                    step!(state, species, dt)
                    r = state.last_report
                    @test r.continuity <= gate
                    @test r.gauss <= gate
                    @test r.magnetic_divergence <= gate
                    @test r.energy_change == sum(r.defects)
                end
            end
        end
        @testset "2D vacuum $T" begin
            g = FourierGrid((7, 7), (T(7), T(7)))
            species = typeof(ParticleSet{2,T}(0))[]
            state = ImplicitEMPICStepper(g, CIC(), species; c = T(2))
            for i = 1:7, j = 1:7
                state.fields.E[3][i, j] = cospi(2 * T(i - 1) / 7 + 2 * T(j - 1) / 7)
            end
            init!(state, species)
            for k = 1:4
                step!(state, species, T(1 / 8))
                r = state.last_report
                @test r.gauss <= gate
                @test r.magnetic_divergence <= gate
            end
        end
    end
end
