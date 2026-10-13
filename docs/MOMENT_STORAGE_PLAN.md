# Moment storage validation

Status: planned; implementation has not started.

The Julia 1.13.1 reproduction in
`/tmp/hybrid-serial-moment-alias-probe.log` shows `momentum!`,
`compute_moments!`, and `compute_moments_multi!` overwriting particle weights
when those weights are passed as scratch. The resulting bulk velocity is
`(2, 6, 24)` for a particle whose velocity is `(2, 3, 4)`.

`momentum!` currently validates only scratch length. The single-population
wrapper deposits density before reaching that check, and the collection wrapper
clears outputs before depositing populations. Validation must precede the first
write in each public operation.

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
