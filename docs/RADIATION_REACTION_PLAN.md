# Radiation reaction validation and arithmetic

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

Recorded status: input validation (RR-001) is pushed as `fdca823`. The normwise
arithmetic repair of RR-002 and RR-003 [below](#normwise-arithmetic-repair-rr-002-rr-003)
replaces the held BigFloat candidate and passes the forcing, particle,
relativistic-boundary, Aqua and JET gates on Julia 1.12.7 and 1.13.2.
Exact-SHA CI is pending.

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

## Normwise arithmetic repair (RR-002, RR-003)

Float16/32/64 particles apply the frozen-rate substep as `v ← F·v`,
`F = 1/√(1 + r(e^{2α} − 1)) ∈ [0, 1]`, `r = 1 − |v|²/c²`, `α = K dt γ χ²`
(algebraically the momentum update), so no stored speed component increases.
The contract is the project's normwise standard: results agree with the exact
substep for the stored inputs to a few ulps of `|v|`
(`src/particles/radiation_reaction.jl`).

1. **Native evaluation** (Float16 particles in Float32): `r` from `c² − |v|²`,
   adding the exact rounding residuals of the squares, their sum and `c²` when
   `β² > 1/2`; `P = (v·E/c)²` squared at the force scale; `χ² = |w|² − P`;
   `X = r(e^{2α} − 1)` (`exp − 1` costs at most 2.7u of `F` in absolute terms,
   so `expm1` is not needed). The guard is an a-posteriori bound
   `errα ≤ G·u·(16|w|² + 8P + 12.1|E|² + 18.2|v|²|B|² + 8 floatmin)`,
   `G = K dt/√r`, accepted when `r ≥ 64u`, `errα ≤ 1/64` and
   `1.09·|dF/dα|·errα ≤ 8 eps(T)`. Non-finite intermediates fail the
   comparison; loop invariants require a normal `K·dt` and
   `4 floatmin/eps ≤ c² ≤ floatmax/4`.
2. **Fallbacks.** Float16 and Float32 particles re-run the same evaluation in
   Float64 (inputs convert exactly and no Float32 input can over- or underflow
   there). Float64 particles retry with a tight bound built from the computed
   products. Remaining particles use exponent-scaled double-double arithmetic:
   power-of-two scaling of `c`, `v` and the fields, `r` from an exact expansion
   near the light cone, `χ² = |w⊥|² + r|E∥|²` as nonnegative terms, and `α`
   as a double-double times `2^e`. `F = 0` for `α ≥ 4096`, `F = 1` for
   `α < 2^-60`, and a logarithmic form once `1 + X > 2^50` underflows to a
   subnormal or zero. Unresolved Lorentz gaps or invariant errors amplified
   beyond the certified budget reconstruct the expression from original stored
   operands with task-local exact integers. Native/DD paths retain zero
   allocation; this cold completion allocates. Hardware particles use no BigFloat.
3. **Other element types** (BigFloat) keep the previous momentum update at the
   ambient precision and rounding mode, bit for bit.

A native result is within `(8 + 4)·eps(T)·|v|` a priori (the bounds
over-estimate the observed error several-fold). The scaled certificate must
include the amplification by `K dt γ`; a small normalized force alone is not
negligible. The [force-span repair and error derivation](RADIATION_FORCE_SPAN_PLAN.md)
replace the earlier Lorentz-gap clamp and complete amplified cancellation from
exact original operands. The derivation states its normalized-DD, IEEE and libm
assumptions; it is not a machine-checked proof.

**Measured design inputs.** The brief's guard (native `r`, fallback below
`r < √eps(T)`) was rejected by measurement: against a 1024-bit reference the
native formula reached 28–38 ulps of `|v|` for `r ∈ [0.01, 1]` and 1,175
(Float32) / 6.5×10⁷ (Float64) ulps inside the accepted range. Budgets of 4, 16
and 64 eps were measured before choosing 8 (64 reached 8.8 ulps). A first
version squared `v·E` before dividing by `c²`, so `(v·E)²` underflowed for
`c = 2^-300`; all 78 resulting matrix failures were repaired by squaring at the
force scale.

**Tests.** `test/test_radiation_arithmetic.jl`, included from
`test_forcing.jl` (39,948 assertions per version): an extreme-scale matrix
with `c`, `|E|`, `|B|` at 2^±5 / 2^±40 / 2^±300 (Float16/32/64), `K dt` scaled
to the force and at independent extremes, both field orientations, all axes
and signs, D = 1–3; light-cone, three-component near-null and field-along-`v`
cases at optical depths around `ln γ`; random ordinary, relativistic and
near-E×B-drift draws; weak drag never heating; exact zero force; the held
cancellation witnesses; BigFloat identity with the previous update; frozen-rate
convergence; zero allocation and task independence. Every case is compared
with an independent 8192-bit momentum expression, normwise at 8 ulps of `|v|`.
The same file asserts per type that the previous arithmetic fails part of the
matrix (RR-002) and heats part of the weak-drag draws (RR-003).

**Results** (Julia 1.12.7 and 1.13.2, identical counts). `test_forcing.jl`
passes 41,607 assertions, `test_particles.jl` 147,
`test_relativistic_validation.jl` 12,801 and Aqua 9. Package JET reports
are 29 / 27 with report sets identical to the parent `d46029e`. Concrete JET over
31 entry-point and tier signatures reports nothing (`report_call`, and
`report_opt` for hardware floats). The largest observed error over the test
cases is 0.50 / 1.28 / 1.80 ulps of `|v|` (Float16/32/64). The previous
arithmetic exceeds 8 ulps or is nonfinite in 204 / 264 / 312 of 2,160 matrix
entries and heats 137 / 115 / 126 of 1,000 weak-drag draws. Every line of
`radiation_reaction.jl` executes under the new test file.

**Performance** (Apple M5, shared machine at load 23–55; warmed calls on
10⁵ particles; minimum ns per particle, alternating new and parent processes;
0 B allocated in every case). Ratios are medians of per-round ratios over three
rounds (Julia 1.12) and five rounds (1.13.2). A faster core state seen by only
one process of a round skewed a few rounds, with per-round extremes of 1.10 and
1.68; no round is excluded.

| Case | Float32 1.12 / 1.13.2 | Float64 1.12 / 1.13.2 |
| --- | --- | --- |
| E, B generic, weak drag (acceptance) | 1.05 / 1.02 | 1.03 / 1.00 |
| pure B, weak drag (acceptance) | 1.04 / 1.02 | 1.04 / 1.01 |
| relativistic, speeds 0.8–0.999 c (all compensated) | 1.19 / 1.17 | 1.19 / 1.15 |
| strong drag, `α ≈ 1` per step (5.9% tier 1) | 1.35 / 1.32 | 1.38 / 1.34 |

In-process, the native kernel costs the same for strong and weak drag, so the
strong-drag excess comes from the 5.9% of particles routed to tier 1. A
64-particle blocked loop that defers those particles was measured and not
adopted (weak drag slower, strong drag about 1 ns faster under load).

## Superseded BigFloat candidate

The held candidate below was never committed; the normwise repair above
replaces it, uses no BigFloat for hardware floats and is therefore unaffected
by NUM-001. Its record is retained.

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
