# Periodic seam consistency

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: pushed separately as `e8925e8` after local acceptance on Julia 1.12.7
and 1.13.2. [Exact-commit CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39893248890)
passes all six jobs. Checked again 2026-10-29; the historical local results
below keep their actual runtime versions.

## Verified behavior

The public witness in
`../HybridPlasmaPIC-audit-evidence/2026-10-23/migration-seam/witness.jl`
retains 144 passes and 72 failures on both Julia 1.12.7 and 1.13.2. It covers
Float32/64, every axis in 1D–3D, negative half-ulp and subnormal inputs,
ordinary interior controls, particle counts and lost-particle statistics.
The frozen migration V3 candidate reproduces the same failures on both runtimes.

For these inputs, `mod(x, 1)` rounds to exactly one. `apply_periodic!` folds that
upper endpoint to zero, as required by its `[lo, hi)` contract and existing
`test_boundary_edges.jl` regression. Direct `rank_of_position` instead assigns
the raw coordinate to the last slab. Migration retains a coordinate equal to
the excluded upper endpoint. The public outputs and the separate implementations
in `boundaries/periodic.jl`, `parallel/domain_decomposition.jl` and
`particles/migration.jl` verify this inconsistency.

## Repair contract

Use one private periodic-coordinate owner for wrapping and folding a rounded
upper endpoint to the lower endpoint. Apply it to the established public
boundary operation, rank classification and private migration preparation.
Retain existing interval/finite-input validation, nonperiodic loss behavior,
particle payloads and object identities. GROWTH-002's prepare/publish ordering
and collective error handling must remain intact.

Inspect every caller before editing. The canonical helper accepts the validated
interval and its already-computed width so per-particle loops retain their
current loop-invariant arithmetic. Ordinary coordinates must keep their bits;
the intended change is the existing public half-open seam rule. No broader
rank-arithmetic or integrator change is part of this repair.

## Acceptance

- [x] Reproduce the public inconsistency on both supported runtimes.
- [x] Implement the shared coordinate owner in a separate candidate.
- [x] Check Float32/64 and 1D–3D seams, nonzero lower bounds, mixed periodic
  layouts, ordinary bit compatibility, unchanged payloads and storage rollback.
- [x] Run real 2/4/8-rank seam routing and subsequent migration controls.
- [x] Check neighboring boundaries, rank lookup, migration, deposition and
  MPI/checkpoint behavior; compare complete compiler reports and allocations.
- [x] Commit and push separately as `e8925e8`.
- [x] `e8925e8`: all six exact-SHA CI jobs pass, as linked above.

The candidate at `../HybridPlasmaPIC-periodic-seam-20261023` layers this change
over frozen GROWTH-002 V3. Its original witness passes 216/216 assertions on both
runtimes. The seven changed files pass 14 formatter AST/idempotence checks;
27 fully concrete compiler calls have no reports on each runtime. All six
focused real-MPI commands exit successfully and check payloads, half-open
coordinates, nonperiodic loss and subsequent valid migration. Source and log
hashes are retained under `migration-seam` in the persistent evidence directory.

The expanded neighboring suite passes 54,480 assertions on each runtime.
The complete single-thread package compiler comparison passes, retaining all
30/27 parent reports and group order. The ordinary-output/allocation matrix
passes 667,584 assertions in each of its four commands. All 108 complete output
and allocation samples match the parent per runtime. The matrix
covers Float32/64, 1D–3D, empty/16/1,024 particles per rank, periodic/mixed
layouts, random interior positions and whole-box offsets. Rank operations must
retain zero warmed allocation; public boundary and migration allocations and
complete output hashes must match frozen V3. No elapsed speedup is claimed.

The first resource harness incorrectly required zero allocation from the public
boundary wrapper. Both unchanged-parent runs retain 667,044 passes and 108
failures, all at that assertion. Independent empty-particle probes pass 12
checks per runtime and measure identical allocation in `apply_periodic!` and
its unchanged `_validated_open_interval` owner: 640 bytes per axis on minimum
Julia and 384 on latest. The corrected, separately retained harness requires
the same measured validation allocation and exact parent/candidate allocation
equality. Production source and numerical assertions are unchanged.

Final review verifies that the MPI storage/collective protocol and migration
prepare/publish logic remain intact. The seven owned files are frozen;
348 other committed Julia files match parent `7c6e1e4`. All accepted native
observations, log hashes and comparisons are checked in
`migration-seam/final-verification.json`. The full review also names the tests
that reject rank-only, storage-only and nonperiodic-folding incomplete fixes.
