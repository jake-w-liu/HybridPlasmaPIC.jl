# Electrostatic PIC with local Coulomb collisions

Status: source implementation; runtime acceptance is deferred. No tests
were added or run. Numerical, stochastic, compiler, formatter and complete-call
resource behavior remain unverified. This is a bounded CAP-010 composition;
general collision/field composition remains open.

`CollisionalElectrostaticPIC` composes the generic synchronous `ElectrostaticPIC`
collection with the existing cell-local `collide_coulomb!`. It supports CPU
Float32/64 in D=1–3, all three velocity components, ordinary Vector particle
leaves, arbitrary finite signed charges, positive masses and nonnegative weights.
All populations are mobile. The inner collection uses its canonical CIC and
periodic spectral Poisson owners, zero stored background and `neutrality=:require`.
Nearest rounding and gradual underflow are required and never changed internally.

```julia
composed = CollisionalElectrostaticPIC(inner, species;
    dt, rng, coulomb_log,
    epsilon0=N0, weight_to_number=N0, u_floor)
init!(composed, species)
step!(composed, species, dt)
diagnostic = electrostatic_diagnostics(composed, species)
save_checkpoint(path, composed, species)
archive_run(archive_path, composed, species; rng_seed)
```

The caller supplies an unprimed generic collection, ordered populations, fixed
positive resolved h, explicit Xoshiro or MersenneTwister, scalar or symmetric
species-by-species logarithms, and collision normalization. The wrapper deeply
copies RNG/logs. Population objects/order/counts/q/m, field/collision settings and
h are fixed. Evolution leaves weights, ids and tags unchanged. Typed empty
collections and empty populations retain the underlying owners' behavior.

## Normalization and physical scope

Both `epsilon0` and `weight_to_number` are required and must resolve to the same
positive grid-precision value. With reference scales L0,T0,M0,Q0, the collision
owner uses epsilon_c=epsilonSI*M0*L0³/(Q0²*T0²) and real count `w_real=N0*w`,
where N0=weight_to_number. The electric unit required by the particle force is
M0*L0/(Q0*T0²). Poisson then has source coefficient N0/epsilon_c, and field energy
in the stored particle-energy units carries epsilon_c/N0. The canonical ES
equations and energy use coefficient 1, so this composition requires epsilon_c=N0.
Physical-count weights use 1,1 in a consistent normalized unit system; another
matching pair can represent density-normalized weights. An arbitrary
`coulomb_normalization(PlasmaUnits(...))` pair must not be substituted unless its
reference units also satisfy this solver's normalization. D<3 retains the local
collision owner's unit reference transverse extent.

The original charge/mass, density deposition, force, coordinate, Poisson and
collision-rate owners remain canonical. Equal macro-weights inherit pairwise
momentum/energy conservation to roundoff; unequal weights inherit conservation
in expectation. Finite collision intervals still have sampling/discretization
error. No fixed bath, compensating reservoir, relativistic collision replacement,
source injection, particle growth, subcycling, or device/distributed evolution
is supplied. The periodic zero-mode, neutrality-tolerance, CIC/spectral and
existing range/conditioning limits remain.

## Time levels, failure and ownership

Initialize with physical full-time velocities. `init!` recomputes the self-field,
resets priming and the wrapper clock, and preserves the current RNG stream.
Reinitialization after evolution requires the caller to supply physical initial
velocities again; it does not recenter carried velocities.

Each positive h first forms the canonical ES kick/drift/Poisson candidate:
x/E/ne are at the new integer boundary and v is at the common preceding half
level. Coulomb scattering over h uses those new positions and retains the
carried velocity label. This field-then-collision Lie split is first order; no
symmetric-split or general accuracy/stability theorem is claimed. The next field
step uses the scattered velocities. Collision changes neither x nor charge,
so there is no second Poisson solve. ES stores no drift-current snapshot.

`source_time` and `source_count` are the wrapper's sole physical clock and
completed collision ledger. The inner ES state has no independent clock. Each
positive operation preflights finite strictly advancing stored time and count
range. Time accumulates with the canonical stored addition; it is not rebuilt as
count*h. Diagnostics reuse canonical particle moments and shared field energy
once, report the stored half-time offset and wrapper clock/ledger, and do not
synchronize particles or project energy. The existing nonrelativistic momentum
reduction's range limitations remain. No exact fully discrete energy law is
claimed.

Fields, particles, priming/history, clock and advanced RNG stage before one
prevalidated publication. An expected input, arithmetic or sampling error leaves
public state and live RNG unchanged. Private particle/Poisson/gather/FFT,
collision and candidate RNG scratch may change and resets before reuse. Zero dt
validates a completed boundary without priming or draws; another positive
converted dt must equal h. All-zero logarithms skip collisions and both RNG
reset/copy operations. Ordinary fixed-size public array identities remain.

The inner state, grid and populations are borrowed exclusively. A direct extra
same-h inner advance after priming leaves its phase flags unchanged and **cannot
be detected** from current inner metadata. Do not advance or refresh it outside
the wrapper. External phase-space/weight edits retain the underlying ES need for
a consistent old field and correctly labelled velocities; there is no injection
owner. Independent states/grids/populations are required for concurrency.
Readonly Fourier metadata, physical/candidate fields, gathers, collision scratch,
particles and writable RNG buffers must keep their declared independent roles.
Resource exhaustion, OS-protected storage, custom setters and concurrent external
mutation are outside the ordinary publication guarantee.

## Restart

Distinct `:collisional_espic` schema 1 requires integer `split_policy=1`,
`collision_position_policy=:new_integer` and
`collision_velocity_policy=:carried_half`. It stores the canonical generic ES
collection payload verbatim, sole source clock/count, fixed collision settings,
and complete advanced RNG/algorithm/runtime. Only initialized ready boundaries
save. A seed provenance record does not replace advanced stochastic state.

Restore validates and detaches the entire incoming graph before publication.
Destinations require matching fixed counts/q/m, field/collision configuration,
normalization, ordinary Vector storage and RNG algorithm/runtime. Saved weights,
ids/tags, phase space and physical E/ne replace values in place. Old destination
physical/clock/inner-phase values may be replaced; source metadata must pass all
clock/phase/history checks. Fresh matching wrappers can load without another
`init!`. Fields restore verbatim, without a Poisson projection or density refresh.
Derived collision/FFT/candidate RNG scratch resets before use. Deterministic ES,
CollisionalHybrid and CollisionalEMPIC payloads cannot substitute; the checksum/
archive container remains version 4.

Source and caller review is separate from acceptance. Existing tests, compiler
assertions and tolerances remain unchanged; no tests were added or run. Actual
coupled convergence, collision relaxation, rollback/retry, advanced-RNG restart,
formatting/JET and complete-call cost remain unverified.
