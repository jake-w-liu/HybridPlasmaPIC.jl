# Actual-pencil electromagnetic PIC with local Coulomb collisions

This source candidate composes rank-owned periodic electromagnetic PIC with
canonical local nonrelativistic Coulomb scattering. No tests were added or run;
runtime, statistics, multi-rank rollback/restart, compiler, format and cost/scaling
acceptance remain unverified. It is a bounded CAP-010 family, not completion of
general source/collision, schedule or backend composition.

Load HybridPlasmaPIC, MPI and PencilFFTs. Every rank calls
`distributed_collisional_empic(inner, species; dt, rng, coulomb_log, epsilon0,
weight_to_number, u_floor=1e-3)` on its actual `distributed_empic` model. The returned
extension-owned `DistributedCollisionalEMPIC` supports collective
`init!(state,species)`, `step!(state,species,dt)`,
`electromagnetic_diagnostics(state,species)`, `save_mpi_checkpoint(dir,state)` and
`load_mpi_checkpoint!(state,dir)`.

The scope is fully periodic3D CPU Float32/64, all-mobile nonrelativistic CIC and
one common fixed positive resolved particle/field interval, strictly below the
inner spectral Maxwell CFL. Use zero background, `neutrality=:require`, explicit
`grid_dc=:project`, ordinary independent Vector particle leaves, finite signed q,
positive m and finite x/v/nonnegative weights. There is no |v|<c admission in this
nonrelativistic model. Ordered species owners and q/m are fixed; local counts and
component leaves can change under migration. Locally empty species are valid.
No fixed population, source injection, species schedule, relativity or device
collision substitution is supplied.

Both normalization keywords are required, finite positive and exactly equal after
canonical grid-precision conversion. With stored macro weights and force-compatible
E units, Maxwell's source coefficient is `weight_to_number/epsilon0`, and field
energy in stored particle units carries its inverse. Equality retains the existing
unit coefficients. Canonical collisions use the original resolved factors and
global cell widths; arbitrary PlasmaUnits collision normalization need not satisfy
this electromagnetic solver's admitted units. The scalar/symmetric matrix Coulomb
logs, u_floor, c/grid/source configuration and h are fixed. Nearest rounding and
gradual underflow are required without changing the caller's mode.

## Carried histories and collision ownership

Initialize with physical full-time velocities and seeded E/B. Canonical init
privately wraps/migrates, forms charge, corrects resolved longitudinal Gauss modes,
retains seeded transverse/null E and B, zeros J and copies the density endpoints.
It resets clock/phase and count bindings while preserving the live RNG stream.
Reinit permits new local counts but requires caller-supplied full-time fields and
velocities; it does not undo an old half phase or supply an injection current.

Every positive step first completes the canonical private EM prime/Boris/drift,
raw midpoint current, represented continuity correction and Maxwell candidate.
E and positions reach the new integer time; B and velocities retain the nominal
negative half-interval offset. Local Coulomb then scatters these carried-half
velocities at the **new positions** over h. This declared field-then-collision
first-order Lie split is not synchronized full-time scattering or a symmetric
split, and supplies no general temporal-accuracy or exact discrete energy theorem.

Scattering changes no x, charge, mass, weight, provenance or field. The eleven
E/B/J/rho_n/rho_np1 arrays retain the completed EM candidate. J is the actual
preceding drift source, rather than a moment redeposited from scattered velocities.
The density endpoints and their represented continuity residual retain that source
history. No postcollision current correction, Gauss repair or second Maxwell solve
runs. The next field step consumes the scattered carried half directly.
Diagnostics report postcollision kinetic/momentum, actual carried E/B energy,
phase offsets and original/grid charge/volume/resolved residual ledgers, with split
and source-clock metadata. They neither synchronize velocities nor imply exact
energy, momentum, continuity or Gauss storage. Unequal macro-weights inherit
conservation in expectation; equal-weight pair conservation is limited by rounding.

Bins cover each rank's actual owned logical cell ranges, independent of parent
array permutation. The canonical periodic collision floor and original global dx
control bins and physical rates. Positive-weight particles must have that cell
on the current migrated CIC-anchor rank before publication; unsupported extreme
stored-geometry ownership disagreement rejects. No universal floor/CIC equivalence
or silent omitted-pair behavior is asserted. Tiny grids/seams retain the existing
collision definition, and no global field/population/bin array is replicated.

Each rank supplies an explicit stock Xoshiro or MersenneTwister. Both live and
candidate streams are deeply owned; algorithms/runtime agree, states need not.
The caller chooses rank streams, with no implicit rank seed or independence claim.
Active scattering uses the private stream. All-zero logs skip collision and both
RNG copies/draws; the common-admitted field map is the ordinary EM map. True zero
dt validates without priming, draws or copies. Nonzero dt that converts to zero
rejects; every positive converted interval must equal the configured h.

The inner time/step Refs are the sole physical clock. Wrapper source_time/count
must equal those actual values, not count*h. Ready operations enforce current
published count/qm/owner bindings. The wrapper exclusively borrows inner/grid/
particles/communicator; bare-inner advance is detected by its source ledger,
while same-count/same-clock external edits cannot all be detected and are unsupported.
All ranks must call in the same order. Concurrent/rank-divergent use is excluded.

## Publication and complete restart

All ranks prepare collision/EM fields, complete particle replacements/count rows,
RNG and metadata before the final gate. Only prevalidated stock copies, ordinary
physical/cache/phase stores and typed source-ledger assignments then publish.
Recoverable local failures preserve physical particles/fields/clock/history/counts
and live RNG. Private FFT/candidate/midpoint/collision/RNG work may change and resets
before reuse. Busy flags release only when acquired by that call. Fatal resource
exhaustion, MPI loss, protected memory, custom setters and concurrent mutation are
outside the ordinary-store contract.

The distinct schema1 `pencil_collisional_empic_manifest.ser` requires fixed split,
new-position/carried-half and canonical-owned-cell policies. One checksummed
generation stores complete nested EM physical state (eleven fields, full ordered
particle leaves/provenance/qm and actual initialized/primed/last_dt/time/step), fixed
configuration/normalization/logs/source ledger and the **full per-rank live RNG**.
It is not an independent EM checkpoint plus an RNG sidecar. Candidate RNG, derived
count caches, collision/midpoint/FFT work, plans and communicator are omitted.

Restore requires matching ranks, actual logical ranges/permutation/parent shapes,
precision, fixed q/m/c/configuration, RNG algorithm and Julia runtime/platform.
Counts may change, including empty local species. Old destination physical values,
particle q/m, count-cache entries, clock/phase and RNG values/metadata can be stale;
valid ordinary storage and immutable cached configuration remain required.
Incoming full fields/phase/current constraints, particles, original source/RNG/
metadata links and physical neutrality validate before detachment/publication.
Both stream objects receive saved full state; public owner/field/cache/Ref identities
persist, and particle leaves rebind as in migration. No current redeposit, field
projection, B recentering, priming, migration or Maxwell solve occurs on load.

Use shared filesystem, trusted Serialization and exclusive directory access.
Closed rank files precede atomic checksummed manifest publication. Earlier failures
preserve the prior manifest; superseded-cleanup failure can accompany a usable new
checkpoint. No crash durability, changed-layout RNG repartition, serial-layout bit
identity or cross-runtime stream continuation is promised. Shared collision/RNG/
generation owners also serve the [collisional ES family](DISTRIBUTED_COLLISIONAL_ELECTROSTATIC.md);
see the [EM model](DISTRIBUTED_EMPIC_IMPLEMENTATION.md) and
[EM checkpoint](DISTRIBUTED_EMPIC_CHECKPOINT.md) for its underlying map and limits.
