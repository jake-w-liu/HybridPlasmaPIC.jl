# §9.4 particle provenance event log: id-keyed records reduce to the full
# §9.4 retain-list (id, species, source region, injection batch, injection time,
# first/last crossing, crossing count, reflection flag, max kinetic energy), and
# the log is invariant under particle reordering (sort/migration).

using HybridPlasmaPIC, Test

@testset "§9.4 provenance event log reduces to the full record" begin
    T = Float64
    log = ParticleProvenanceLog()
    id1 = global_particle_id(1; species = 0)
    id2 = global_particle_id(2; species = 3)
    record_injection!(log, id1, 0.0; source_region = 7, batch = 2)
    record_injection!(log, id2, 0.5; source_region = 4, batch = 2)
    record_crossing!(log, id1, 1.0)
    record_crossing!(log, id1, 3.0)
    record_crossing!(log, id1, 2.0)
    record_reflection!(log, id1, 1.5)

    s1 = provenance_summary(log, id1)
    @test s1.id == id1
    @test s1.species == 0
    @test s1.injection_time == 0.0
    @test s1.source_region == 7
    @test s1.injection_batch == 2
    @test s1.crossing_count == 3
    @test s1.first_crossing_time == 1.0
    @test s1.last_crossing_time == 3.0
    @test s1.reflection_flag == true

    s2 = provenance_summary(log, id2)
    @test s2.species == 3                       # species recovered from the id high bits
    @test s2.crossing_count == 0
    @test isnan(s2.first_crossing_time)
    @test s2.reflection_flag == false
    @test s2.max_kinetic_energy == 0.0
end

@testset "§9.4 maximum kinetic energy tracked over a run; log survives reordering" begin
    T = Float64
    ps = ParticleSet{1,T}(3)
    load_lattice_1d!(ps, 0.0, 1.0)
    set_density_weight!(ps, 1.0, FourierGrid((3,), (1.0,)))
    assign_global_particle_ids!(ps, 1; species = 0)
    ps.v[1] .= [1.0, 2.0, 0.5]                  # m = 1 ⇒ ke = ½ v²
    log = ParticleProvenanceLog()
    record_max_kinetic_energy!(log, ps, 0.0)    # ke = 0.5, 2.0, 0.125
    ps.v[1] .= [0.2, 1.0, 3.0]
    record_max_kinetic_energy!(log, ps, 1.0)    # ke = 0.02, 0.5, 4.5
    # per-particle running maximum over the two snapshots
    @test provenance_summary(log, ps.id[1]).max_kinetic_energy ≈ 0.5
    @test provenance_summary(log, ps.id[2]).max_kinetic_energy ≈ 2.0
    @test provenance_summary(log, ps.id[3]).max_kinetic_energy ≈ 4.5

    # the id-keyed log is invariant if the particles are physically reordered:
    # reverse the SoA arrays and confirm each id still resolves to the same record.
    reverse!(ps.id)
    reverse!(ps.v[1])
    @test provenance_summary(log, ps.id[1]).max_kinetic_energy ≈ 4.5   # was particle 3
    @test provenance_summary(log, ps.id[3]).max_kinetic_energy ≈ 0.5   # was particle 1

    # whole-log batch reduction covers exactly the logged ids
    summ = provenance_summary(log)
    @test length(summ) == 3
    @test summ[global_particle_id(2; species = 0)].max_kinetic_energy ≈ 2.0
end
