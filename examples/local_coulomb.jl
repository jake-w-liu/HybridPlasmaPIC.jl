# Run with: julia --project=. examples/local_coulomb.jl
# Stand-alone periodic collision step in the package's normalized hybrid units.
using HybridPlasmaPIC
using Random

units = PlasmaUnits(n0 = 1e6, B0 = 5e-9, mi = 1.6726e-27)
grid = FourierGrid((8,), (8.0,))
species = (ParticleSet{1,Float64}(256; q = 1, m = 1), ParticleSet{1,Float64}(128; q = 2, m = 4))
rng = MersenneTwister(24)
for ps in species
    rand!(rng, ps.x[1])
    ps.x[1] .*= grid.L[1]
    for component in ps.v
        randn!(rng, component)
        component ./= sqrt(ps.m)
    end
    set_density_weight!(ps, 1.0, grid)
end
workspace = CoulombCollisionWorkspace(species, grid)
collide_coulomb!(
    species,
    grid,
    0.01;
    coulomb_log = 10,
    coulomb_normalization(units)...,
    workspace,
    rng,
)
println("Advanced two species through one cell-local Coulomb substep.")
# Keep the evolved RNG alongside particles for exact continuation. Existing
# simulation checkpoints do not persist this caller-owned RNG, and this example
# does not select a collision/field splitting for a staggered time integrator.
