# Electromagnetic PIC with local Coulomb collisions

Status: source implementation; runtime acceptance is deferred. No tests
were added or run for this slice. Numerical, stochastic, compiler, formatter and
resource behavior have not been measured. This is a bounded CAP-010 capability;
it does not complete general collision/field composition.

`CollisionalEMPICStepper` composes the generic synchronous `EMPIC` species
collection with the existing cell-local `collide_coulomb!`. It supports CPU
Float32/64 in one to three spatial dimensions, all three velocity components,
ordinary Vector particle leaves, and arbitrary finite signed charges with
positive masses and nonnegative weights. Every population is mobile. The inner
state must have `n_sub=1`, `relativistic=false`, zero background charge,
`neutrality=:require` and no species schedule owner. Opted-in scheduled states
are rejected even with all-fast rows and `n_sub=1`; the wrapper's common-fast
collision phase and schema 1 remain unchanged. NGP, CIC and TSC retain their
existing source/gather rules.
Nearest rounding and gradual underflow are required; no arithmetic mode is
changed. The wrapper borrows the inner state exclusively.

The public operations are:

```julia
composed = CollisionalEMPICStepper(inner, species;
    dt, rng, coulomb_log,
    epsilon0=N0, weight_to_number=N0, u_floor)
init!(composed, species)
step!(composed, species, dt)
diagnostic = electromagnetic_diagnostics(composed, species)
save_checkpoint(path, composed, species)
archive_run(archive_path, composed, species; rng_seed)
```

Here `inner`, the ordered populations, fixed positive resolved `dt`, explicit
Xoshiro or MersenneTwister RNG, collision logarithms and normalization are supplied
by the caller. A scalar or symmetric species-by-species matrix supplies
`coulomb_log`; zero entries disable channels. The wrapper deeply copies the RNG
and logarithms. Particle objects, order, counts, charges, masses, field/collision
configuration and dt are fixed. Initialization requires physical full-time
velocities and seeded physical E/B, resets the physical clock, and preserves the
current RNG stream. Reinitializing a previously advanced inner state requires
the caller to supply physical initial velocities and fields again.

## Normalization and physical scope

Both `epsilon0` and `weight_to_number` are required and must resolve to the same
positive grid-precision value. This is a consistency condition, not a collision
rate fitted to the field step. With reference scales L₀, T₀, M₀, Q₀, the local
collision owner uses εc=εSI M₀L₀³/(Q₀²T₀²) and real count
`w_real=N0*w`, where `N0=weight_to_number`. The force-compatible electric-field
unit is M₀L₀/(Q₀T₀²). Maxwell's source coefficient is then N0/εc, and the field
energy in the stored particle-energy units carries εc/N0. The existing EMPIC
equations and diagnostic use unit coefficients, so this composition requires
εc=N0. Physical-count weights use the pair 1,1; density-normalized weights may
use another matching pair. A `coulomb_normalization(PlasmaUnits(...))` result
cannot be substituted unless its reference units also satisfy this EMPIC
normalization. D<3 retains the collision owner's unit reference transverse extent.

The collision kernel retains its original masses, signed charges, timestep,
widths, count conversion and pair-rate factors. Equal macro-weights inherit its
pairwise momentum/energy conservation to roundoff. Unequal macro-weights inherit
its conservation in expectation. Finite collision intervals still have sampling
and discretization error. No relativistic collision substitution, bath, fixed
background energy reservoir, source injection, particle growth, independent
species schedule or device/distributed evolution is supplied.

## Ordering, publication and diagnostics

For each fixed h, the canonical field candidate advances positions to n+1,
velocities and B to n+1/2, and E to n+1. On the first step it performs the existing
physical-state priming. Local Coulomb scattering then acts over h at the new
positions while velocities retain that common half-time label. This declared
field-then-collision Lie split is first order, with no symmetric-split or general
temporal-order theorem. The next field step uses the scattered velocities.

Collisions leave positions, charges, masses, weights and provenance unchanged.
The stored current remains the actual preceding drift source, and charge
snapshots retain the same discrete continuity relation. No postcollision current
refresh or second Maxwell update is performed. The canonical Fourier/Nyquist and
range/conditioning limits remain. Diagnostics delegate the canonical stored-time
particle moments and shared field energy, with nominal half-time offsets and the
completed source clock/count. The nonrelativistic momentum reduction retains its
existing range limitations; this composition adds no separate moment formula.
Diagnostics do not synchronize velocities or project energy. No
exact fully discrete energy, momentum or Gauss-storage law is claimed.

The field map, collisions and RNG advance on private candidates before one
prevalidated publication. A failed input, field, arithmetic or sampling stage
leaves public fields, particles, clock, initialization/phase, source ledger and
live RNG unchanged. Private FFT, field/particle/RNG candidates and collision
scratch may change and are reset before use. All-zero logarithms skip collision
and RNG reset/copy; the field candidate/publication owners are the ordinary EMPIC
owners. Zero dt validates and preserves the completed boundary without priming
or RNG draws. Other nonzero converted dt values must equal the configured h.

Standard Array/Vector storage identities remain. Resource exhaustion,
OS-protected memory, custom setter side effects and concurrent external mutation
are outside the ordinary publication guarantee. Calls require independent
states, grids, populations and RNGs. Directly advancing or changing the borrowed
inner state outside the wrapper invalidates its source/phase contract. External
particle or field edits retain the existing EMPIC requirement for consistent
sources and declared velocity phases; the wrapper adds no injection-current model.

## Restart and remaining work

Distinct `:collisional_empic` schema 1 requires integer `split_policy=1`,
`collision_position_policy=:new_integer` and
`collision_velocity_policy=:carried_fast_half`. It contains the canonical generic
EMPIC collection payload verbatim, the complete advanced RNG/algorithm/runtime,
fixed collision settings and completed source ledger. Only initialized ready
boundaries save. The original field/current/charge arrays restore without a
projection or recomputation. Full RNG state, not seed provenance, determines
stochastic continuation in the same Julia runtime/platform.

Restore validates and detaches the complete incoming graph before publication.
Destinations require matching fixed population counts and ordinary Vector
storage; saved weights, ids/tags, phase space and physical fields replace their
values in place. Field/collision configuration, normalization and RNG algorithm
must match. Fresh matching wrappers can load without `init!` afterward. Old
deterministic EMPIC and CollisionalHybrid payloads cannot substitute for this
schema; the shared checksum/archive container remains version 4.

Source/caller/inverse review is separate from runtime acceptance. Existing tests,
compiler assertions and tolerances remain unchanged. Collision-rate/relaxation,
coupled convergence, rollback/retry, advanced-RNG continuation and complete-call
cost remain unverified. Broader CAP-010 physics and backends remain open.

The separate [actual-pencil collisional EM source](DISTRIBUTED_COLLISIONAL_EMPIC.md)
now composes the distributed model with the same local scheduler on admitted
owned cells. It preserves actual drift-current/field history and full per-rank
RNG/physical state in one matching-layout generation. Local counts may change
through migration. Runtime, statistics and restart acceptance remain deferred;
the serial schema above and the distributed schema are distinct.
