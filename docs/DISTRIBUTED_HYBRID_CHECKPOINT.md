# Complete actual-pencil isothermal hybrid checkpoint

This is the complete matching-layout persistence source candidate for
`DistributedHybrid`, paired with the actual-pencil isothermal hybrid implementation.
No tests were added or run. Runtime, numerical continuation, multi-rank rollback,
compiler, format, liveness and resource acceptance remain unverified.

Every rank calls `save_mpi_checkpoint(dir, model)` on the same shared directory,
or `load_mpi_checkpoint!(model, dir)` into a matching configured model. Calls must
follow the same collective order and exclude concurrent grid/model/directory use.
The caller owns the communicator, FFT plan and filesystem lifetime. Serialization
inputs must be trusted; this format is not an untrusted-data decoder.

`pencil_hybrid_manifest.ser` is a distinct checksummed schema1 manifest with
`:distributed_isothermal_hybrid` model identity and explicit integrator, moment,
represented derivative, source-risk, layout and field-phase policies. It cannot
silently load an electrostatic, electromagnetic or serial hybrid checkpoint.
Physics and configuration retain their exact stored scalar types/values, including
isothermal Te, original and converted nfloor, original absent η/ηH/de2, q/m,
global n/L/dx, rank count and logical ownership ranges.

Each rank stores all eighteen public physical scalar arrays:

- B and E, three components each;
- density n and reciprocal density ninv;
- ion velocity ui, current J and pressure gradient gradp, three components each;
- electron pressure pe.

Complete ordered local particle x/v/weight/id/tag/q/m leaves and
`initialized`, `primed`, `last_dt`, `time`, `step` accompany those fields. Positions
and fields are at the nominal integer boundary, while advanced particle velocities
retain the saved `-last_dt/2` offset. Exposed endpoint-predictor moments and carried
E are saved directly. Restore does not derive them from stored half-time velocities.
At initialized step0, supplied velocities and phase remain full time. Time is not
reconstructed as step*dt, and magnetic RK subcycle NB is not persistent history.

Save requires initialized ready state, actual owner/count/qm bindings, valid
phase, finite fields and density at least nfloor. No density clipping, neutrality
condition, Poisson solve, field projection or exact-divB condition is introduced.
The physical hybrid charge/volume/divB diagnostics remain separate from checkpoint
integrity; the manifest does not imply exact Gauss, continuity or energy conservation.

The format omits derived count caches, floor-count scratch, RK stages, candidate
fields, particle/midpoint/predictor/gather buffers, layout lookup, FFT plans and
communicator. It owns no RNG or stochastic continuation state. The saved particle
leaves determine local counts when restoring; a species can become empty locally
and counts can differ from destination counts. The nonempty ordered species-slot
list and exact cached q/m still match the destination.

Restore requires matching precision, immutable physics/configuration, actual rank
count, owned logical ranges, parent dimensions and pencil permutation. Destination
physical values, particle q/m, clock/phase, local counts and count-cache entries may
be stale. Valid ordinary storage, cache lengths, field/particle/metadata/FFT
independence and cached configuration remain required. Current ready-state/count
admission is not called to validate incoming replacements.

Before publication every rank checks manifest/rank checksums, strict keys and
types, fixed policies, full incoming phase, independent finite eighteen-field
storage, density bounds, complete positive-ion particles and actual CIC ownership.
Incoming links to destination physical/private/metadata/particle/FFT storage reject
before detachment; incoming fields cannot alias particles or each other. A detached
copy validates again. Same-type complete replacement leaves, q/m/count rows and
all eighteen source/destination array references prepare before the final all-rank
gate. Ordinary leaf/scalar/field/cache/Ref stores then publish the boundary. Model
owners, public field arrays, cache Vectors and clock/phase Ref identities persist;
particle component arrays can rebind as in normal migration.

No deposit, curl, Ohm, predictor, recentering, priming, migration, Poisson or field
repair runs on restore. Old state is preserved on recoverable preparation failure;
private scratch may change during later normal operations and is not a snapshot.
The supplied file values define the carried boundary rather than a reconstructed
minimum field state.

Each save uses a unique owned generation directory. Closed rank files are hashed
before an atomic checksummed manifest publishes that generation. Existing matching
configuration and old rank checksums validate before replacement. Earlier failures
preserve the prior manifest and clean only canonical owned generation entries;
superseded-cleanup failure reports that the new checkpoint is already usable.
Symlink/path rules and error composition reuse the seven shared pencil checkpoint
mechanics. MPI loss, fatal resource exhaustion, protected memory, concurrent
mutation/use and crash durability are outside this ordinary-store contract.
Changed decomposition/rank count or backend conversion is not supplied.

See the [model implementation guide](DISTRIBUTED_HYBRID_IMPLEMENTATION.md) for the
physical map, original-product Faraday admission boundary and excluded physics.
