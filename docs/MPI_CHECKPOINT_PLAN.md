# Distributed checkpoint publication and collection restart

Status: implementation in progress after `d2c3da7`; no completion claim.

## Verified baseline

`/tmp/hybrid-mpi-save-baseline-probe.log` shows a serializer exception truncating
an existing rank file from 834 to 486 bytes while its manifest remains unchanged.
The prior checkpoint no longer deserializes. The same probe confirms that ordered
population collections have no applicable MPI save or load method.

## Contract

- All ranks participate with the same grid, ordered species parameters and
  integrator phase. Each rank retains its own population counts and provenance.
- Use the existing serial collection representation and restore preflight.
  Validate every rank before changing any restore destination or publishing files.
- A checkpoint directory belongs to one grid/layout and has exclusive access
  during save or load. Concurrent readers/writers are outside this interface.
- Write new generations under a unique owned subdirectory. Publish their manifest
  only after every rank file has closed successfully, using same-filesystem
  `Base.Filesystem.rename`. Julia 1.12's implementation directly calls the OS
  rename operation; `mv(force=true)` has a delete/copy fallback and is unsuitable.
- New MPI schema 2 records the generation. Read existing schema 1 flat-file
  checkpoints. Serial checkpoint schemas remain unchanged.
- Recoverable failures before publication retain the prior manifest and files.
  Delete only files owned by the attempted generation. Cleanup failures preserve
  both the primary error and the cleanup error.
- After publication, remove superseded checkpoint files. If that cleanup fails,
  report the failure with the new checkpoint already published and usable.
  No crash/power-loss durability or rollback of subsequent allocation failures
  is claimed; interrupted processes can leave an unpublished generation.

## Acceptance

- Scalar schema 1 reading and schema 2 replacement/restart.
- Ordered collection exact continuation with changing timestep/subcycles,
  migration, unequal counts, empty populations and independent writable views.
- One-rank validation, serializer, filesystem, malformed later-population and
  non-resizable destination failures; all peers report errors and remain usable.
- Prior checkpoint bytes and restart survive failed publication; successful
  repeated saves remove superseded generations.
- Both Julia 1.12.7 and 1.13.1; real MPI, affected serial callers, concrete JET,
  formatter and final diff checks.

The tests distinguish plausible incomplete repairs: atomically replacing only
an already-corrupted manifest cannot preserve the prior rank files; validating
only the local population strands peers or partially restores them; concatenating
populations loses ordered charge/mass/provenance; ignoring cleanup errors hides
unpublished files or a partially removed old generation.

## Verification in progress

The immutable v1 candidate repairs the original corruption probe: all 834 bytes
of the previous rank file remain intact after the serializer error, the manifest
is unchanged, and the old checkpoint deserializes. Log:
`/tmp/hybrid-mpi-save-fixed-probe.log`.

The first broad Julia 1.12 run passed 715 assertions and had one fixture error:
a concrete Vector of vector-backed ParticleSets rejected insertion of a
view-backed ParticleSet. The fixture now uses the intended common ParticleSet
supertype. The first two-rank Julia 1.13 run passed 360 continuation assertions
per rank, then reached that fixture error; it was stopped through its launch
session (exit 130). Neither run is a final full-pass claim. Final tests also add
failed-cleanup error preservation and cross-rank manifest/species mismatch cases.

A confined v1 probe also verified that replacing a generation directory with a
symlink caused cleanup to remove the linked directory's rank file. Log:
`/tmp/hybrid-mpi-generation-symlink-probe.log`. Replacement preflight and cleanup
now reject generation symlinks; the regression requires the prior manifest and
external rank-file bytes to survive. This guard is pending final verification.
