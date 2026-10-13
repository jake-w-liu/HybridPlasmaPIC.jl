# Radiation reaction validation and arithmetic

Status: input validation repair verified on Julia 1.12.7 and 1.13.1 and pushed
as `fdca823`. The arithmetic candidate passes serial numerical gates on both
versions, but a reproduced concurrent BigFloat failure remains unresolved.
The candidate is uncommitted and the radiation audit remains open.

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

## Arithmetic candidate under verification

The candidate evaluates a velocity factor in [0,1], avoiding the momentum
round trip that can heat weakly damped particles. Range-sensitive force
invariants use exact stored inputs; the complete damping factor is evaluated
with higher precision and logarithms. Finite-range, weak-drag, retained-precision
and first-order synchrotron convergence tests pass **16,299 assertions on each
Julia version**. Combined forcing, particles and relativistic-boundary checks
pass **30,906 per version**. Source SHA256:
`7f50d7b08bd7ea71bfbeabce2c3791e18fdfa965fb2759023966dcb4bd2dd7d6`.
Logs: `/tmp/hybrid-radiation-arithmetic-v2-focused112.log`,
`/tmp/hybrid-radiation-arithmetic-v2-focused113.log`,
`/tmp/hybrid-radiation-arithmetic-v2-neighbors112.log`,
`/tmp/hybrid-radiation-arithmetic-v2-neighbors113.log`.

Both versions pass 36 concrete JET calls with zero reports (108 assertions);
package counts remain 26/25. Measured ordinary calls retain 0 B across 27
configurations. Exceptional single-particle cases allocate 17,984–37,048 B.
Seven alternating local timing pairs give new/old ratios of 1.38–3.06, including
the identical velocity reset in both sides. Further cost reduction is needed;
the machine was shared. Logs: `/tmp/hybrid-radiation-arithmetic-v2-jet112.log`,
`/tmp/hybrid-radiation-arithmetic-v2-jet113.log`,
`/tmp/hybrid-radiation-arithmetic-v2-resources112.log`.

VERIFIED: the previous public arithmetic increases the exact stored velocity
norm under weak positive drag for Float16/32/64 on both versions (six witness
assertions each). Logs: `/tmp/hybrid-radiation-weak-baseline112.log` and
`/tmp/hybrid-radiation-weak-baseline113.log`. A first candidate also lost a tiny
nonzero force when its cancellation bound underflowed; the corrected selector
requires that bound to retain normal range. Six original-candidate witnesses
and 54 final regression checks retain that distinction.

**Open verification failure:** concurrent calls on Julia 1.12 intermittently
return almost unchanged BigFloat velocities for a tiny-field case. A captured
53-bit input returns about `(0.5, 0.25, 0)` instead of the independent
`(0.1769635, 0.0884818, 0)` result; serial replay of the identical saved inputs
matches the reference. Eighty repeated warm task batches pass, but fresh
processes reproduce the failure. Passing repeats do not close it. Logs:
`/tmp/hybrid-radiation-context-fresh112-4.log`,
`/tmp/hybrid-radiation-context-capture112-16.log`,
`/tmp/hybrid-radiation-context-replay112.log`. Native reduction below establishes
the cause; the ordinary installed-runtime failure remains open.

A finer captured trace verifies that the wrong intermediate is `log(rwide)`:
for stored `rwide=0.6875`, Julia 1.12.7 concurrently returned approximately
-45.87345. Exact invariants and optical depth were correct in that trace.
Log: `/tmp/hybrid-radiation-context-trace-v3-112-23.log`; captured operands:
`/tmp/hybrid-radiation-context-trace-v3-53-RoundingMode{:Up}().bin`.
This trace identifies the failing operation. The native reduction below
separately establishes its cause; passing repeated warm calls alone does not
establish a repair.

The reduction now reproduces the wrong logarithm using only Julia Base,
without importing this package (`/tmp/hybrid-radiation-base-concurrency112-24.log`).
Direct concurrent MPFR constant calls also return incorrect pi values on both
Julia 1.12.7 and 1.13.1, each using unpatched MPFR 4.2.2 with shared caches:
`/tmp/hybrid-mpfr-cache-concurrency112-1.log` and
`/tmp/hybrid-mpfr-cache-concurrency113-1.log`. The native cache-lock transition
is now verified by an instrumented C reproducer and a matched one-line repair.
Installed runtimes remain unchanged.
[Evidence and remaining integration](MPFR_CONCURRENCY_PLAN.md) distinguish the
temporary repair from ordinary Julia execution; no package-local workaround
has been adopted.

A branchless product-predicate experiment preserves 7,803 represented-operand
checks and ordinary measured 0 B allocation, but its nine alternating local
pairs have mixed new/candidate median ratios of 0.93–1.10. It has not been
promoted: these shared-machine measurements establish no consistent benefit.
Log: `/tmp/hybrid-radiation-guard-benchmark-v2-112.log`. A preceding benchmark
harness failed its allocation assertion and is excluded from timing evidence.
A separate forced-inlining prototype also has mixed local ratios (0.87–1.05),
including a slower large Float32 batch; it remains unpromoted. Log:
`/tmp/hybrid-radiation-inline-benchmark112.log`. Neither experiment changes the
frozen candidate or closes its recorded performance and runtime gates.
