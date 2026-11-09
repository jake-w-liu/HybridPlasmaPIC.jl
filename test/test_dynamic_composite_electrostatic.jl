using Test
using HybridPlasmaPIC

# Immutable reports contain BigInts; compare their fields by numeric value.
function _dynamic_composite_snapshot(st, species)
    s = HybridPlasmaPIC._checkpoint_state(st, species)
    values(r) =
        r === nothing ? nothing :
        NamedTuple{fieldnames(typeof(r))}(ntuple(i -> getfield(r, i), fieldcount(typeof(r))))
    epoch(e) = merge(e, (last_report = values(e.last_report),))
    history = s.history === nothing ? nothing : merge(s.history, (state = epoch(s.history.state),))
    return merge(
        s,
        (
            current = epoch(s.current),
            history = history,
            last_regrid_report = values(s.last_regrid_report),
        ),
    )
end

@testset "Dynamic composite electrostatic essentials" begin
    H = HybridPlasmaPIC
    for T in (Float32, Float64)
        @testset "$T" begin
            base = StaticAMRHierarchy1D(4, T(8))
            fine = StaticAMRHierarchy1D(4, T(8), ((1, 2, 2),))
            p = ParticleSet{1,T}(1)
            p.x[1][1], p.v[1][1] = T(7 / 2), T(3)
            species = [p, ParticleSet{1,T}(0)]
            st = DynamicCompositeElectrostaticPIC(
                base,
                species;
                shape = YeeTriangularCloud(T(1 / 4)),
            )
            init!(st, species)
            @test st.E == T[0, 1/4, 1, 1]
            @test st.phi == T[0, 0, -1/2, -5/2, -9/2]
            leaves = (p.x[1], p.v..., p.weight, p.id, p.tag)
            regrid!(st, species, fine)
            # Symmetric cloud at 7/2: nodal source is (0,0,1/2,1/2,0,-1).
            # Half sum(width*E^2) changes from 33/16 to 17/8, hence +1/16.
            @test st.E == T[0, 0, 1/2, 1, 1]
            @test st.phi == T[0, 0, 0, -1/2, -5/2, -9/2]
            @test st.mesh_work == 1 // 16 && iszero(st.physical_defect)
            @test st.time == 0 && st.step == 0 && st.regrid_count == 1
            @test p.x[1][1] == T(7 / 2) && p.v[1][1] == T(3)
            E, report = st.E, st.last_regrid_report
            @test regrid!(st, species, fine) === st
            @test st.E === E && st.last_regrid_report === report && st.regrid_count == 1
            before = _dynamic_composite_snapshot(st, species)
            @test step!(st, species, -zero(T)) === st
            @test_throws ArgumentError regrid!(st, species, StaticAMRHierarchy1D(4, T(16)))
            @test isequal(_dynamic_composite_snapshot(st, species), before)
            h = T(1 / 4)
            step!(st, species, h)
            # Exact KDK: vhalf=49/16, x1=273/64, endpoint field=111/128.
            @test p.x[1][1] == T(273 / 64) && p.v[1][1] == T(3247 / 1024)
            @test all(
                pair -> pair[1] === pair[2],
                zip(leaves, (p.x[1], p.v..., p.weight, p.id, p.tag)),
            )
            regrid!(st, species, base)
            @test st.step == 1 && st.time == h && st._state.step[] == 0
            @test st.E == T[0, 0, 111/128, 1]
            d = composite_electrostatic_diagnostics(st, species)
            @test iszero(d.accounting_defect) &&
                  d.particle_charge == 1 &&
                  d.right_plate_charge == -1
            @test d.last_physical_step.state.configuration.refinements == fine.refinements
            @test d.last_physical_step.global_step == 1 && d.last_physical_step.state.step == 1
            # Repeated regrid with no physical step retains the preceding OLD mesh history.
            history = deepcopy(d.last_physical_step.state.last_impulse)
            regrid!(st, species, fine)
            d = composite_electrostatic_diagnostics(st, species)
            @test d.last_physical_step.state.last_impulse == history
            @test d.last_physical_step.state.configuration.refinements == fine.refinements
            saved = _dynamic_composite_snapshot(st, species)
            @test_throws ArgumentError step!(st, species, h / T(2))
            @test isequal(_dynamic_composite_snapshot(st, species), saved)
            raw = H._checkpoint_state(st, species)
            bad_origin = merge(raw, (history = merge(raw.history, (epoch_start_time = h / T(2),)),))
            @test_throws ArgumentError H._restore_checkpoint_state!(st, species, bad_origin)
            aliased = merge(raw, (current = merge(raw.current, (E = st.E,)),))
            @test_throws ArgumentError H._restore_checkpoint_state!(st, species, aliased)
            @test isequal(_dynamic_composite_snapshot(st, species), saved)
            mktempdir() do directory
                checkpoint, archive =
                    joinpath(directory, "dynamic.chk"), joinpath(directory, "dynamic.run")
                save_checkpoint(checkpoint, st, species)
                archive_run(archive, st, species; rng_seed = 0)
                # A different topology/count and stale old values must be replaced, not reinitialized.
                other = [ParticleSet{1,T}(2), ParticleSet{1,T}(0)]
                other[1].x[1] .= T.((1, 2))
                destination = DynamicCompositeElectrostaticPIC(
                    base,
                    other;
                    shape = YeeTriangularCloud(T(1 / 4)),
                )
                init!(destination, other)
                other[1].x[1] .= T(NaN)
                other[1].weight .= T(NaN)
                other[1].q, other[1].m = T(NaN), -one(T)
                fill!(destination.E, T(NaN))
                destination.time, destination.step, destination.last_dt = T(Inf), -1, T(Inf)
                owner = other[1]
                load_checkpoint!(destination, other, checkpoint)
                @test other[1] === owner && length(other[1].weight) == 1
                @test isequal(_dynamic_composite_snapshot(destination, other), saved)
                step!(st, species, h)
                step!(destination, other, h)
                continued = _dynamic_composite_snapshot(st, species)
                @test isequal(_dynamic_composite_snapshot(destination, other), continued)
                @test st.step == 2 && st.time == 2h && st.epoch_start_step == 1
                @test st._history === nothing
                load_archive!(destination, other, archive)
                @test isequal(_dynamic_composite_snapshot(destination, other), saved)
                step!(destination, other, h)
                @test isequal(_dynamic_composite_snapshot(destination, other), continued)
            end
            # Independent curvature rule on a base field: both interior cells tag.
            tagged = [ParticleSet{1,T}(1)]
            tagged[1].x[1][1] = T(7 / 2)
            adaptive =
                DynamicCompositeElectrostaticPIC(base, tagged; shape = YeeTriangularCloud(T(1 / 4)))
            init!(adaptive, tagged)
            adapt!(adaptive, tagged; threshold = T(1 / 4), buffer = 0)
            @test adaptive.mesh.refinements == ((1, 2, 3),)
            @test adaptive.last_regrid_report.selection == :base_electric_curvature
            adapt!(adaptive, tagged; threshold = floatmax(T), buffer = 0)
            @test isempty(adaptive.mesh.refinements)
            @test iszero(adaptive.mesh_work) && iszero(adaptive.physical_defect)
            @test iszero(composite_electrostatic_diagnostics(adaptive, tagged).accounting_defect)
        end
    end
end
