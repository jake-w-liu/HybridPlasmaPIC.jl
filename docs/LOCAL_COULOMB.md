# Local Coulomb collisions

`collide_coulomb!(species, grid, dt; ...)` advances one nonrelativistic collision
substep on a periodic `FourierGrid`. It accepts ordered tuples or vectors of
independent CPU `ParticleSet`s, all with the grid's dimension and Float32/64
precision. Self and cross-species pairs are selected within each cell. The
existing `collide_coulomb!(ps, gcoeff, dt)` method remains the prescribed-rate,
homogeneous single-population operator.

## Set the units explicitly

The required keywords are `coulomb_log`, `epsilon0` and `weight_to_number`.
`coulomb_log` is a nonnegative scalar or symmetric species-by-species matrix;
zero entries disable channels. Neutral species also have no active channels.
`weight_to_number` converts each stored macro-weight into a real-particle count.
`epsilon0` uses the same units as the supplied masses, charges, lengths and time.
The positive `u_floor` regularizes relative speed and defaults to 1e-3 in the
particle velocity units. Rescale it explicitly when changing unit systems.

For the package's hybrid normalization, use `coulomb_normalization(units)` with
the same `PlasmaUnits` as the fields and particles. Its return value supplies
`epsilon0` and `weight_to_number`; it does not choose a Coulomb logarithm or
velocity floor. Use `set_density_weight!` with density normalized by `units.n0`.
In one or two dimensions, the omitted directions have unit normalized extent;
physical cell volume is `prod(grid.dx)*d_i^3`. The helper therefore fixes the
reference transverse length/area used to interpret each macro-weight.

The executable [two-species example](../examples/local_coulomb.jl) shows the
complete setup. Its call is:

```julia
workspace = CoulombCollisionWorkspace(species, grid)
collide_coulomb!(species, grid, 0.01;
    coulomb_log = 10,
    coulomb_normalization(units)...,
    workspace, rng)
```

Allocate the workspace once and reuse it. Particle-count changes resize its
buffers; species count and grid cell count must still match. Use a separate
workspace and RNG for each concurrent collision call.

## Physical and numerical contract

The pair variance uses the reduced physical mass, cell volume, weighted partner
density and each reused partner's pairing count. Odd self populations use
half-strength triplets; zero-weight particles never pair. Equal macro weights
conserve pair momentum and energy to roundoff. Unequal-weight acceptance gives
statistical conservation and must be assessed across seeds or ensembles.

The shared Takizuka–Abe angle law switches to isotropic scattering when its
variance exceeds unity. The physical prefactor uses mathematical π and retains
original masses, cell widths and integer partner/reuse counts. Binary scaling
avoids intermediate rate overflow/underflow; a certified comparison resolves
the sharp regime boundary. Resolve the collision timestep and check relaxation
against timestep/cell/population refinement: finite-step accuracy does not
follow from conservation alone. The operator does not select the Coulomb
logarithm, validate a weak-coupling model for a physical application, or account
for relativistic collision physics.

Widened, compensated and exact evaluation preserve the stored sampled expression
where ordinary products or cancellation lose range. Exact completion can
allocate. Finite particle inputs are required; an accepted proposal that cannot
be represented as a finite particle velocity raises `ArgumentError`.

## State, failures and coupled evolution

Validation precedes sampling and workspace staging. Velocities are staged until
all active channels finish, then published for every population. Validation and
sampling failures leave particle state unchanged; RNG and scratch may have
advanced. Positions, weights, identifiers, tags, charges and masses are unchanged. The result is the
same species collection supplied by the caller.

This is a stand-alone operator. A coupled field integrator must choose an
explicit collision splitting and apply it at consistent particle time levels.
The package does not yet provide that composed integrator or store this external
RNG in simulation checkpoints. For exact continuation, retain the evolved RNG
alongside the ordered particle state and preserve the execution environment;
recreate the derived workspace. GPU and distributed collision execution remain
separate implementation work.

The [implementation tracker](LOCAL_MULTISPECIES_COULOMB_PLAN.md) records the
verification evidence and remaining composed-evolution/backend/restart gates.
