using Test
using LinearAlgebra
using HybridPlasmaPIC

# Compare every saved report field by value after detaching its BigInt scalars.
function _static_composite_snapshot(state, species)
    saved = HybridPlasmaPIC._checkpoint_state(state, species)
    report = saved.last_report
    values = report === nothing ? nothing : NamedTuple{fieldnames(typeof(report))}(
        ntuple(i -> getfield(report, i), fieldcount(typeof(report))))
    return merge(saved, (last_report = values,))
end

@testset "Static composite electrostatic essentials" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64)
        @testset "$T" begin
            mesh = StaticAMRHierarchy1D(4, T(8), ((1, 2, 2),))
            @test mesh.vertices == [0, 2, 3, 4, 6, 8]
            @test mesh.widths == [2, 1, 1, 2, 2]
            @test sum(mesh.dual_volumes) == 8
            @test mesh.patches[1].valid == [true, false, true, true]
            @test mesh.patches[1].covered == [false, true, false, false]
            @test_throws ArgumentError StaticAMRHierarchy1D(4, T(8), ((1, 2, 2), (1, 2, 3)))
            @test_throws ArgumentError StaticAMRHierarchy1D(4, T(8), ((1, 2, 2), (2, 1, 1)))
            particle = ParticleSet{1,T}(1)
            particle.x[1][1], particle.v[1][1] = T(7 / 2), T(3)
            species = [particle, ParticleSet{1,T}(0)]
            state = StaticCompositeElectrostaticPIC(mesh, species; shape = YeeTriangularCloud(T(1 / 4)))
            init!(state, species)
            # The symmetric cloud lies wholly in [3,4]. Assemble the mixed
            # Neumann matrix independently from the leaf incidence and widths.
            source = T[0, 0, 1 / 2, 1 / 2, 0, -1]
            matrix = zeros(T, 6, 6)
            for cell in eachindex(mesh.widths)
                conductance = inv(T(mesh.widths[cell]))
                matrix[cell, cell] += conductance
                matrix[cell + 1, cell + 1] += conductance
                matrix[cell, cell + 1] -= conductance
                matrix[cell + 1, cell] -= conductance
            end
            potential = vcat(zero(T), matrix[2:end, 2:end] \ source[2:end])
            @test isapprox(state.phi, potential; rtol = T(16) * eps(T), atol = T(16) * eps(T))
            @test state.phi == T[0, 0, 0, -1 / 2, -5 / 2, -9 / 2]
            @test state.E == T[0, 0, 1 / 2, 1, 1]
            @test patch_electric_field(state, 1) == T[0, 1 / 4, 1, 1]
            diagnostics = composite_electrostatic_diagnostics(state, species)
            @test diagnostics.nodal_source == source
            @test diagnostics.particle_charge == 1 && diagnostics.right_plate_charge == -1
            before = _static_composite_snapshot(state, species)
            @test step!(state, species, -zero(T)) === state
            @test isequal(_static_composite_snapshot(state, species), before)
            @test_throws ArgumentError step!(state, species, -one(T))
            @test isequal(_static_composite_snapshot(state, species), before)
            position_leaf = particle.x[1]
            h = T(1 / 4)
            step!(state, species, h)
            # Both endpoint clouds lie within single cells. Exact Verlet gives
            # vhalf=49/16, x1=273/64 and endpoint E=111/128 after crossing x=4.
            @test particle.x[1][1] == T(273 / 64)
            @test particle.v[1][1] == T(3247 / 1024)
            @test particle.x[1] === position_leaf
            @test state.E == T[0, 0, 0, 111 / 128, 1]
            @test state.last_impulse == [0, 0, 1 // 2, 17 // 128, 0]
            @test state.last_report[].interface_crossings == [0, 1]
            report = state.last_report[]
            @test report.field_solve_defect == 0
            @test report.total_energy_defect == report.field_solve_defect + report.particle_work_defect
            saved = _static_composite_snapshot(state, species)
            bad = merge(H._checkpoint_state(state, species), (time = h / T(2),))
            @test_throws ArgumentError H._restore_checkpoint_state!(state, species, bad)
            @test isequal(_static_composite_snapshot(state, species), saved)
            mktempdir() do directory
                checkpoint, archive = joinpath(directory, "composite.chk"), joinpath(directory, "composite.run")
                save_checkpoint(checkpoint, state, species)
                step!(state, species, h)
                particle.x[1][1], particle.v[1][1], particle.weight[1] = T(NaN), T(NaN), T(NaN)
                particle.q, particle.m = T(NaN), -one(T)
                fill!(state.E, T(NaN)); fill!(state.phi, T(NaN))
                load_checkpoint!(state, species, checkpoint)
                @test isequal(_static_composite_snapshot(state, species), saved)
                archive_run(archive, state, species; rng_seed = 0)
                step!(state, species, h)
                load_archive!(state, species, archive)
                @test isequal(_static_composite_snapshot(state, species), saved)
            end
            # A positive subnormal duration keeps an exact half duration, even
            # though that half cannot itself be stored in the field precision.
            empty_species = [ParticleSet{1,T}(0)]
            vacuum = StaticCompositeElectrostaticPIC(mesh, empty_species; shape = YeeTriangularCloud(T(1 / 4)))
            init!(vacuum, empty_species)
            tiny = nextfloat(zero(T))
            step!(vacuum, empty_species, tiny)
            @test vacuum.time[] == tiny && vacuum.last_dt[] == tiny && vacuum.step[] == 1
        end
    end
end
