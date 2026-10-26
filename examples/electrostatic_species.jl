# A cold pair plasma: both charge signs move in the same periodic Poisson field.
using HybridPlasmaPIC

grid = FourierGrid((64,), (2π,))
species = [ParticleSet{1,Float64}(512; q = charge, m = 1) for charge in (-1, 1)]
for ps in species
    load_lattice!(ps, (0.0,), grid.L, (512,))
    set_density_weight!(ps, 1.0, grid)
    ps.x[1] .+= ps.q .* 0.001 .* sin.(ps.x[1])
end
state = ElectrostaticPIC(grid, species) # :require rejects nonneutral physical charge
init_espic!(state, species)
for _ = 1:100
    step_espic!(state, species, 0.01)
end
println(electrostatic_diagnostics(state, species))
# Optional checkpoint: save_checkpoint("pair-plasma.checkpoint", state, species)
