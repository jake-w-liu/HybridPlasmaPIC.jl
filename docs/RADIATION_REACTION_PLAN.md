# Radiation reaction validation and arithmetic

Status: input validation repair verified on Julia 1.12.7 and 1.13.1 and pushed
as `fdca823`. Extreme-range arithmetic remains a separate confirmed repair;
the full radiation audit is open.

## Verified gaps

Both Julia 1.12.7 and 1.13.1 accept a finite Float64 electric or magnetic field
of magnitude `1e100` for Float32 particles. Conversion makes it infinite;
the resulting NaN invariant silently skips the update. The repaired public boundary
rejects nonfinite converted fields before mutation, using the same field
validator as the particle pusher. Coefficient and timestep validation also
require finite values in particle precision, including no-op calls. Negative
source values that narrow to signed zero remain invalid.

Independent 2,048-bit expressions also verify incorrect finite cooling at
extreme light speeds, field magnitudes and coefficient/timestep products.
Twelve Float32/64 examples return either unchanged velocity or zero where the
complete documented exponential substep has a finite intermediate result.
Each native-version probe passes 48 range-witness and four narrowing-witness
assertions. Logs: `/tmp/hybrid-radiation-baseline112.log` and
`/tmp/hybrid-radiation-baseline113.log`.

## Scope and acceptance

The first slice validates inputs once before the public loop. It retains the
signature, existing velocity preflight and valid zero-step behavior. Test every
field component, all supported dimensions, Float16/32/64, wider real inputs,
invalid values during no-ops, unchanged particles and ordinary numerical
equivalence. Run the complete forcing suite, concrete JET, allocation checks
and formatting on both supported Julia versions before committing.

The first candidate passes 1,659 forcing assertions on each Julia version, but
adds 2,880 B to ordinary calls. The shared `_validated_uniform_field` constructs
component-name strings before testing finiteness. The final repair moves that
construction into its failure branch, preserving conversion order and error
text. This extends the canonical-owner slice to `boris.jl`; its existing pusher
callers and radiation allocation/state comparisons pass as recorded below.

The final v2 implementation passes 1,810 assertions across 22 forcing, particle
and long-orbit testsets on each Julia version. All 27 concrete compiler calls
per version have zero reports, with 81 analysis/runtime checks passing. Package
JET retains 26/25 reports against the unchanged ceiling of 30; the existing MPI
report classification is in [MPI_CHECKPOINT_PLAN.md](MPI_CHECKPOINT_PLAN.md).
Logs: `/tmp/hybrid-radiation-inputs-v2-regressions112.log`,
`/tmp/hybrid-radiation-inputs-v2-regressions113.log`,
`/tmp/hybrid-radiation-inputs-v2-jet112.log`,
`/tmp/hybrid-radiation-inputs-v2-jet113.log`.

All 27 warmed Float16/32/64, 1D–3D, 0/1/256-particle cases retain the baseline's
exact velocity arrays and 0 B allocation (82 comparison checks). An additional
54 checks pass on both baseline and final shared validators, preserving error
text/order and converted tuple precision. Logs:
`/tmp/hybrid-radiation-inputs-v2-comparison112.log`,
`/tmp/hybrid-radiation-field-errors-baseline112.log`,
`/tmp/hybrid-radiation-field-errors-fixed112.log`.
The final post-test source edit changes only the docstring to describe the
frozen-rate approximation and remaining range limitation; the function body
matches the tested v2 snapshot exactly.

The arithmetic slice must retain the complete radiation invariant and optical
depth before materialization. Verify cooling direction, no heating, subluminal
outputs, zero-force states, extreme cancellation and finite neighboring values
against independent high-precision expressions. Measure ordinary and exceptional
costs; do not silently skip NaN invariants or equate an overflowing intermediate
with physically complete damping.

This API implements the dominant drag term, evaluated in uniform fields. Its
exponential substep freezes the initial drag coefficient; it is not an exact
solution of the nonlinear continuous equation. The input-validation slice now
documents that distinction and the remaining numeric-range limitation; temporal
convergence remains an acceptance test for the arithmetic repair. The wider classical/quantum
model distinction is described in the primary
[Smilei radiation documentation](https://smileipic.github.io/Smilei/Understand/radiation_loss.html).
