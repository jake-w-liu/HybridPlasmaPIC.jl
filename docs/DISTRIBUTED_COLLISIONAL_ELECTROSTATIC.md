# Actual-pencil electrostatic PIC with local Coulomb collisions

This source candidate composes the existing rank-owned periodic electrostatic
PIC with canonical local nonrelativistic Coulomb scattering. Runtime, multi-rank,
statistical, restart, compiler, formatting and resource acceptance are unverified.
No tests were added or run for this implementation. This is a bounded CAP-010
family; distributed collisional hybrid/EM evolution remains separate work.

Load HybridPlasmaPIC, MPI and PencilFFTs. Every rank calls
`distributed_collisional_electrostatic_pic(inner, species; dt, rng, coulomb_log,
epsilon0, weight_to_number, u_floor=1e-3)` with its actual
`distributed_electrostatic_pic` inner model and ordered local populations. The
returned extension-owned `DistributedCollisionalElectrostaticPIC` supports
`init!(state, species)`, `step!(state, species, dt)`,
`electrostatic_diagnostics(state, species)`, `save_mpi_checkpoint(dir, state)` and
`load_mpi_checkpoint!(state, dir)` collectively.

Use a 3D fully periodic CPU Float32/64 grid, all-mobile ordinary Vector particle
leaves, zero background, `neutrality=:require` and explicit `grid_dc=:project`.
The same ordered species, q/m, grid/layout and collision configuration are fixed.
Empty local populations are valid; migration replaces complete component leaves
while retaining ParticleSet owners. Particle counts are not a frozen binding.
Finite signed charges, positive masses and finite nonnegative weights are required.

Both normalization keywords are required and must resolve to the same positive
working-precision value. The stored Poisson coefficient is
`weight_to_number/epsilon0`, and the stored field-energy coefficient is its inverse.
Equality retains the canonical unit coefficients. Collision rates use the original
resolved factors and the original global cell widths. A generic PlasmaUnits
collision normalization need not satisfy this solver's reduced normalization.
The scalar or symmetric species-by-species Coulomb logarithms are copied; their
resolved values, `u_floor` and the positive resolved fixed interval cannot change.
Nearest rounding and gradual underflow are required, never set internally.

## Phase, collision cells and RNG ownership

Initialize with full-time velocities. Canonical private first-step priming shifts
them backwards by half an interval, then the shared ES mover advances positions
and solves charge/E once. Collisions act over the fixed interval at the **new
integer-time positions and carried half-time velocities**. They change no charge,
position, weight or provenance, so no repeated Poisson solve is performed. The
next ES kick consumes the scattered carried half directly. This is a field-then-
collision first-order Lie split, without a symmetric-split or exact discrete total-
energy claim. Unequal weights retain the canonical conservation-in-expectation
law; equal-weight pair conservation is limited by stored arithmetic.

Bins cover only the rank's actual logical owned ranges, independent of physical
array permutation. The canonical periodic collision floor determines global
cells; global dx determines collision volume/rates. No global particle or bin
array is allocated. Before a step publishes, every positive-weight particle must
have its canonical collision cell on its current migrated CIC-anchor rank. If
extreme stored geometry violates this condition, the operation rejects; a universal
CIC/floor ownership equivalence is not asserted. Accepted complete cell ownership
allows the unchanged serial self/cross/sampling scheduler to run locally. Tiny
cells/grids and periodic seams retain that scheduler's existing cell definition.

Each rank supplies an explicit stock Xoshiro or MersenneTwister. Construction
deeply detaches live and candidate RNG state. Ranks agree on algorithm/runtime,
not stream state; callers choose streams, with no automatic rank-seeding or
independence guarantee. Active collisions stage against the candidate stream.
All-zero logarithms skip collisions and both RNG copies, with no draws. On the
common admitted domain, this branch keeps the canonical ES map. Zero interval
validates without priming, copying or drawing; a nonzero value that converts to
zero rejects. Every positive interval must equal the configured resolved dt and
advance the stored clock; unresolved half intervals reject.

The inner time/step Refs are the sole physical clock. Wrapper `source_time` and
`source_count` must equal their actual stored values; time is not count*dt.
Reinitialization resets physical phase and this ledger while preserving the live
RNG stream. The wrapper exclusively borrows the inner/grid/communicator and
populations. Bare-inner advance is detected by ledger mismatch; same-clock edits
or refreshes cannot all be detected and are unsupported. This is not an injection
or source-refresh API. All ranks must use the same collective call order.

## Transaction and complete restart

Particles, fields, inner clock/phase, source ledger and advanced live RNG publish
after every rank passes private field/collision/final preparation. Ordinary stock
copies and fully prepared particle-leaf rebinding follow that gate; no public
resize, validation or stochastic draw remains. Recoverable local errors preserve
public physical state and live RNG. Private FFT/particle/collision/candidate-RNG
scratch can change and is reset/reused. Busy acquisition releases only flags
acquired by the failed call. MPI failure, fatal resource exhaustion, protected
memory, concurrent mutation/use and rank-divergent calls are excluded.

`pencil_collisional_electrostatic_manifest.ser` has a distinct schema1 model and
required split, collision-position, carried-velocity and canonical-owned-cell
policies. One checksummed generation stores complete E/rho and ordered local
particles, actual phase/clock/source ledger and the **full per-rank stock RNG**.
It is not an independent ES save followed by an RNG sidecar. Private collision,
FFT, migration and candidate-RNG scratch, plans and communicators are omitted.

Restore requires the same rank count, actual ownership/permutation, precision,
q/m, fixed parameters, RNG algorithm and Julia runtime/platform. Counts may change,
including empty populations. Incoming physical fields restore verbatim without
migration, Poisson or phase reconstruction. Incoming full RNG values, graph
independence, source ledger, checksums and physical neutrality validate before
publication. Old destination physical/clock/phase/RNG values may be stale; valid
ordinary writable storage and fixed configuration are still required. Both live
and reusable candidate stream receive the saved complete stream after the gate.

Use a shared filesystem and trusted Serialization data with exclusive directory
access. Earlier generation/write failures preserve the prior manifest; a reported
superseded-cleanup failure can accompany a usable newly published checkpoint.
There is no crash-durability, changed-layout/rank stream repartition or cross-runtime
stream-continuation promise. See [the physical ES implementation](DISTRIBUTED_ELECTROSTATIC_IMPLEMENTATION.md)
and [its checkpoint contract](DISTRIBUTED_ELECTROSTATIC_CHECKPOINT.md).
