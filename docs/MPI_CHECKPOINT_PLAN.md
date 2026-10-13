# Distributed checkpoint publication and collection restart

Status: implemented and verified on Julia 1.12.7 and 1.13.1. Included in this
checkpoint commit; the wider project audit remains active.

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

## Verification

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

A confined candidate probe also verified that replacing a generation directory with a
symlink caused cleanup to remove the linked directory's rank file. Log:
`/tmp/hybrid-mpi-generation-symlink-probe.log`. Replacement preflight and cleanup
now reject generation symlinks; the regression requires the prior manifest and
external rank-file bytes to survive. The repaired probe passed
(`/tmp/hybrid-mpi-generation-symlink-fixed112.log`), and both-version real-MPI
regressions now cover this ownership boundary.

The final production source is frozen in the v3 and v4 snapshots. V4 additionally
enables the valid-but-different manifest test that final review found missing
from the case list; production source is unchanged.

| Check | Verified result and log |
| --- | --- |
| Full real two-rank MPI suite, Julia 1.12.7 | Pass, including 360 collection restart checks per rank and 48,049 collection evolution checks per rank; `/tmp/hybrid-mpi-generation-v3-fullmpi112.log` |
| Collection checkpoint suite, real two-rank Julia 1.13.1 | 360 continuation, 28 restore-failure and 48 publication checks per rank; `/tmp/hybrid-mpi-generation-v3-mpi113.log` |
| Final restore/publication cases, real two-rank both versions | 30 restore-failure and 48 publication checks per rank; `/tmp/hybrid-mpi-generation-v4-faults-mpi112.log`, `/tmp/hybrid-mpi-generation-v4-faults-mpi113.log` |
| Final single-rank failure cases, Julia 1.13.1 | 58 checks; `/tmp/hybrid-mpi-generation-v4-faults-self113.log` |
| Earlier adjacent serial checkpoint/MPI suite, Julia 1.12.7 | 736 checks, before the additional symlink and manifest-disagreement regressions; `/tmp/hybrid-mpi-generation-v2-serial112.log` |
| Final formatting | All six changed Julia files pass JuliaFormatter 1.0.62 without edits; `/tmp/hybrid-mpi-generation-final-format.log` |

Package JET retains the ceiling of 30: the final source reports 26 on Julia 1.12
and 25 on Julia 1.13, compared with baseline 24/22. All 24 concrete save/load
calls per version execute successfully. Their two additional reports concern
`buf` in MPI.jl's `bcast`: `collective.jl:89–104` assigns it in both branches of
the same immutable root flag before use. Both branches execute in the real-MPI
tests. Two save Serialization reports also occur on the scalar baseline.
The extra latest-Julia package-only report infers an empty grid-size tuple and
`Vector{Union{}}`; the public FourierGrid constructor requires a nonempty size
tuple. A native 1.13 probe confirms no empty two-argument constructor and checks
Int/Float32/Float64 descriptor arrays for every supported dimension (12 checks).
No report suppression or gate increase was added. Logs:
`/tmp/hybrid-mpi-generation-v3-jet112.log`,
`/tmp/hybrid-mpi-generation-v3-jet113.log`,
`/tmp/hybrid-mpi-generation-grid-types113.log`.

The warmed resource probe covers Float32/64, 1D–3D, a scalar 64-particle population
and collections of 64/17/0 particles on 4-cell-per-axis grids. Across five
replacements per case, 120 checks confirm one active generation with only its
rank file and no accumulating old generations. Save allocations are
24,640–43,072 B; restore allocations are 35,648–82,304 B. Fixture disk usage is
3,067–9,867 B. These local measurements are not cluster scaling evidence. Log:
`/tmp/hybrid-mpi-generation-v3-resources-warm112.log`.
