# Moment storage validation

Status: implemented and verified on Julia 1.12.7 and 1.13.1 in immutable
snapshots; the repair is included with this evidence.

The Julia 1.13.1 reproduction in
`/tmp/hybrid-serial-moment-alias-probe.log` shows `momentum!`,
`compute_moments!`, and `compute_moments_multi!` overwriting particle weights
when those weights are passed as scratch. The resulting bulk velocity is
`(2, 6, 24)` for a particle whose velocity is `(2, 3, 4)`.

At the reproduced baseline, `momentum!` validated only scratch length. The
single-population wrapper deposited density before reaching that check, and the
collection wrapper cleared outputs before depositing populations. Validation
now precedes the first write in each public operation.

## Implementation contract

- Grid outputs and grid scratch must have the expected dimensions and independent
  storage. Writable buffers must not overlap particle input storage.
- Particle scratch must have the required vector type, element type, axes and
  length, and must not overlap grid buffers or any population's input storage.
- Check every population and buffer before modifying outputs, particles or
  scratch. Preserve independent writable views and sequential reuse of a scratch
  vector when that reuse cannot overwrite inputs or accumulated outputs.
- Reuse the existing particle-storage alias predicates. Consolidate common
  validation in the moment layer and route serial and MPI moment entry points
  through it. MPI input failures must be reported collectively before deposition.
- Keep number-density, charge-weighted density, normalization and valid arithmetic
  unchanged. Document buffer requirements in the public moment docstrings.

## Verification

1. Reproduce weight, position, velocity, output and cross-population scratch
   overlaps. Assert rejection before every destination or input changes.
2. Exercise late failures: the last momentum component and a later population.
   This rejects a patch that validates only the currently processed component.
3. Cover all shapes, 1D–3D, Float32/64, empty populations, independent views and
   valid shared scratch reused sequentially. Compare with independent sums and
   the existing valid-path output.
4. Test real MPI with a failure on one rank and successful reuse afterward.
   A local exception must not strand peers in the reduction.
5. Run affected pressure/current, hybrid/CAM-CL, collection and MPI regressions,
   concrete compiler checks and both supported Julia version gates. Retain the
   existing `test_physics2.jl` allocation limit and measure validation overhead.

A length-only guard, a local-only MPI check, or copying aliased scratch after an
output has already been cleared would each leave part of the defect unresolved.

## Current evidence

- Final focused tests: 1,956 assertions on each Julia version. Real two-rank
  tests: 264 assertions per rank on each version, including failed-input
  recovery and valid scalar/collection moments.
- The earlier broad 1.12 run passed 2,145 neighboring assertions and 1,446
  assertions in the new file, but exited with one fixture error: position views
  used different index types, violating the ParticleSet constructor contract.
  The corrected final file passes all 1,956 assertions above; production files
  are byte-identical between those snapshots. The existing allocation gate in
  `test_physics2.jl` passed without modification.
- Package JET: 24 reports on Julia 1.12 and 22 on Julia 1.13 against the unchanged
  ceiling of 30. All 36 concrete moment calls per version have zero reports;
  72 analysis/runtime checks pass on each. The final 1.12 check uses a separately
  resolved Julia 1.12 tooling environment with JET 0.11.6.
- A matched old/new comparison passes 288 checks. All sampled grid outputs and
  scratch results agree exactly; all 18 warmed scalar/collection/pressure cases
  allocate 0 B. Alternating shared-machine samples give paired latency ratios
  of 0.80–1.24; these measurements do not establish a speedup.

Logs: `/tmp/hybrid-moment-storage-v2-focused112.log`,
`/tmp/hybrid-moment-storage-v2-focused113.log`,
`/tmp/hybrid-moment-storage-v1-mpi112.log`,
`/tmp/hybrid-moment-storage-v2-mpi113.log`,
`/tmp/hybrid-moment-storage-v1-regressions112.log`,
`/tmp/hybrid-moment-storage-v1-jet113.log`,
`/tmp/hybrid-moment-storage-v2-jet112.log`, and
`/tmp/hybrid-moment-storage-v2-comparison113.log`.
