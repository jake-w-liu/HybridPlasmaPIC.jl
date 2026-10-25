# Coulomb sampled-angle range repair

Status: the numerical repair is pushed as `2a5cc12` after final local validation;
exact-SHA CI and fixed latency acceptance remain open. Baseline is `a18f09d`, including the exact-null
repair. This slice does not accept the full local multispecies feature.

## Pre-edit brief

- Mode: repair the prescribed-rate Coulomb sampler and velocity reconstruction.
- Requirement: preserve representable nonzero physical kicks when an
  intermediate variance, tangent, squared tangent or unit-direction product
  underflows. Preserve exact null rotations, RNG ordering, public validation
  and transactional velocity publication. The complete sampled expression,
  rather than a midpoint reconstruction, owns exceptional range arithmetic.
- Baseline evidence: both Julia 1.12.7 and 1.13.2 fail four of eight checks in
  `coulomb-null/legacy-variance-underflow-v1.jl`; source/commands/logs and
  independent BigFloat references are recorded in
  `coulomb-null/variance-underflow-review-v1.json`. These failures occur in both
  `839c109` and the exact-null repair; that repair did not cause them.
- Root cause: `collisions.jl` forms `var = numerator / denominator`, then
  `sqrt(var)`. The variance may round to zero while its square root and final
  velocity kick remain representable. Reconstructing midpoint plus/minus a
  rounded rotated vector also loses small increments. The latter class must
  receive independent failing witnesses before being marked repaired.
- Canonical owners: `collisions.jl` owns the legacy public API, validation,
  pairing/rejection and scratch publication. The frozen local-collision V10
  candidate has a shared `_coulomb_scattered_relative` sampler and a complete
  physical-update expression. Extract those mathematical owners once; do not
  introduce a competing scattering law or alter the frozen V10 worktree.
- Affected surfaces: source/include order, private sampler/pair callers,
  permanent collision tests, local-candidate composition, README, changelog
  and the shared implementation tracker. No public API, restart schema,
  backend or collision-model expansion is part of this repair.
- Compatibility: public signatures, RNG events, error behavior and null
  identities remain. Retain every ordinary comparison profile and report
  changed output bits. Stable increments intentionally change historical
  floating-point trajectories; reproducibility within the repaired version
  remains required. Validate ordinary results against an independent physical
  update and conservation/rate invariants, not only against the old midpoint
  implementation. Never retain an erased kick for bitwise compatibility. Native ordinary
  Float32/64 batches must retain zero allocations; generic and exceptional
  costs are measured and reported separately.
- Strategy: retain the original finite rate factors through exceptional
  sampling. Use mantissa/exponent arithmetic before the square root, keeping a
  binary-scaled tangent through the final physical expression. Ordinary
  sampling preserves its operation order. Complete-angle arithmetic may use
  the existing widened/compensated/exact owners when their documented bounds
  apply; this is a numerical evaluation path, not a silent degraded result.
  Use one common expression for all such paths. A threshold may select more
  accurate arithmetic; it must never turn a nonzero angle into an identity.
- Verification: retain before/after commands, hashes and logs for both
  runtimes. Cover Float16/32/64 and ordinary-precision BigFloat, all coordinate
  orientations, signs, unequal weights, pair/triplet populations, tiny and zero
  draws, coefficient/timestep product underflow/overflow, variance boundaries,
  near-backscatter and cancellation. Check componentwise independent references,
  conservation where applicable, scratch/error/RNG contracts, ordinary full
  outputs and allocations, concrete/full JET, Aqua, neighbor tests, actual
  package test sandboxes and whole-tree formatting. Recheck the composed local
  candidate before its eventual acceptance. Fixed latency and exceptional
  resource acceptance remain explicit gates, not inferred from untimed tests.

## Wrong-patch controls

1. Replacing a zero variance with a minimum positive constant gives the wrong
   physical rate; vary coefficient, speed and Gaussian draw independently.
2. Computing a better square root but still subtracting rotated vectors or
   rebuilding midpoints loses oblique and cancellation-sensitive increments;
   compare each final component with an independent complete expression.
3. Returning for every small angle hides finite kicks and changes draw order;
   retain exact-null, nonzero-small-angle, invalid-azimuth and rejection tests.
4. Copying a second sampler makes local and legacy laws diverge; inspect all
   callers and test a separately composed local candidate against the common
   mathematical owner.

## Tracker

- [x] Reproduce the erased-kick variance witness on both runtimes.
- [x] Trace public inputs, pair/triplet calls, normalized weights, scratch and
  the existing local complete-angle owner.
- [x] Extend independent failing witnesses to coefficient products, oblique
  increments and neighboring angle boundaries before production edits. Both
  runtimes retain 14,836 passes and 2,444 failed component checks, with no test
  errors (`coulomb-angle/before-v1-*`).
- [x] Implement one shared scaled-angle and complete-update path.
- [x] Verify numerical, ordinary, RNG/error, allocation and compiler contracts.
- [x] Review exact final source for the bounded legacy numerical repair.
- [ ] Reassemble and reverify the separate local-collision candidate.
- [x] Commit/push the bounded repair individually as `2a5cc12`.
- [ ] Complete all six exact-SHA CI jobs and full compiler-tree review.
- [ ] Complete fixed latency acceptance and remaining CAP-001 integration.

The frozen local sampler reconstructs exceptional tangents from logarithms;
its existing complete-angle tests use a 4,096-epsilon relative bound. Whether
that path satisfies this slice's stricter componentwise checks is untested.
Run a direct physical-rate witness before reusing it unchanged. No global
accuracy or optimization claim follows from this plan.

## Prototype review and compatibility decision

The first prototype retains 17,088 passes and 192 failures on each runtime.
All remaining failures are BigFloat rate-product profiles: their variance is
representable at BigFloat's wide exponent range, but midpoint reconstruction
still loses component accuracy. `diagnose-v1-112.toml` records every case.
The failure requires a stable increment update, not a hardware-type exception
or a larger error bound.

The initial plan overconstrained compatibility by requiring historical output
bits except for individually proved wrong results. The checked public contract
describes the physical rotation; `test/test_collisions.jl` verifies repeated
execution with the same RNG on the same implementation. It does not promise
historical-version bit patterns. The existing local-collision plan explicitly
reserves stable-increment migration for a separately reviewed numerical fix
with correctness, compatibility and cost acceptance. This is that repair.
Accordingly, use the common physical increment for all nonidentity updates,
retain the entire old/new differential inventory as diagnostic evidence and
add independent ordinary physical oracles. No previously failed gate is being
reported as passed; the baseline and first prototype remain retained.

## Verified correction to the large-midpoint regression

The first package-sandbox run fails the old `test_collisions.jl` assertion
that the pair `(9e307, 1e308)` remains unchanged in all velocity components.
The independent `large-midpoint-oracle-v1.jl` replays the same random events
and evaluates `theta = 2atan(draw/sqrt(abs(g)^3))` and the axis-aligned rotation
with 4,096-bit trigonometry. Both runtimes pass 72 checks over three seeds and
three signed/neighboring velocity pairs. The parent erases all transverse
kicks; the candidate agrees componentwise with the physical oracle within
the already declared 32-epsilon bound. The original sandbox failures and
source hashes remain in `coulomb-angle/sandbox-v1-*`.

The regression now checks the unchanged longitudinal velocities and each
nonzero transverse kick against this independent expression. This corrects
an invalid no-motion expectation; it does not loosen the numerical tolerance
or skip the near-range case. Original-factor regime-boundary replay also
confirms 42 passes and six failures on `a18f09d`, on both runtimes;
`regime-parent-v1-*` retains those records. The candidate's focused matrix
passes all 17,328 checks on both runtimes. Full acceptance is still open.

## Ordinary cancellation repair brief

The independent 1,024-bit Rodrigues probe records 88 final-component relative
differences on each runtime (`ordinary-physical-oracle-v1-*`). A second probe
separates the mathematical rotation from the native sampled angle and geometry:
an independently assembled matrix evaluates the latter expression at high
precision. It finds 76 components beyond the unchanged 32-epsilon relative
arithmetic bound, all admitted to the ordinary update. The complete path has
no failed component in this matrix. All 12,288 normwise physical checks pass
per runtime; these do not replace either retained componentwise measurement.

The current cancellation screen only selects complete evaluation after losing
half the floating-point significand. That permits substantial amplification
of ordinary rounding error before the screen is reached. The next bounded
candidate screens subtraction once it loses a leading bit, in both the
increment and final velocity owners, then reuses the existing complete
expression. This is a conservative dispatch screen, not a universal proof of
the final error bound. No angle, rate, rejection, identity, tolerance or oracle
changes accompany it. Rerun the complete sampled-expression matrix, expanded
held-out inputs, original range/identity/rollback tests and every ordinary
allocation profile; latency and exceptional costs remain acceptance gates.

Native sampled quantities and exact mathematical geometry are distinct
reference levels. The stage arithmetic contract uses the former. The broader
physical probe retains its component-relative failures separately, alongside
the normwise physical result. Do not claim a component-relative bound against
an infinitely precise rotation from these stage-arithmetic checks.

The first decomposition accidentally recomputed BigFloat's perpendicular
radius inside its higher-precision block. Hardware references were unaffected;
the 22 BigFloat candidates required a corrected replay. The corrected matrix
captures every sampled geometric operand before increasing precision. It passes
all 24,576 arithmetic and physical-norm checks on each runtime with the stronger
screen. An executed wrong-patch control restoring only the old screens fails
73 of those assertions on each runtime, including 54 hardware-float and 19
BigFloat component failures. Source, command, failure and correction records
are retained as `decomposition-fixture-correction-v1.json`,
`ordinary-decomposition-v3-*` and `mutation-old-cancellation-v1-*`.
The stricter mathematical-rotation component comparison still records 45
differences per runtime; these belong to angle/geometry rounding, rather than
failure of the sampled-expression arithmetic gate. No global relative-error
bound is claimed for components arbitrarily close to zero.

## Subnormal pair-weight factor

**VERIFIED:** `subnormal-weight-ratio-before-v1-*` records 47 passes and 25
failed public component checks on each runtime. The six-particle fixtures
pair two equal subnormal weights and leave four unit-weight particles at rest.
The native `wmax/wbar` quotient is subnormal; multiplying it by a large finite
rate returns to the normal range after relative information has already been
lost. An independent 4,096-bit original-factor reference detects the wrong
kick in Float16, Float32 and Float64.

Before this correction, the sampler's ordinary-path predicate checked only the
final numerator, denominator and coefficient product. It now also requires a
normal weight quotient; otherwise the existing scaled-factor owner evaluates
the original inputs. This preserves the law and existing resource
targets. Add the public witness permanently and rerun the unchanged range,
ordinary, error, allocation and compiler gates. The separate population-mean
defect remains tracked in [its own plan](COULOMB_WEIGHT_MEAN_PLAN.md).

The isolated predicate correction passes all 72 witness assertions on both
runtimes (`subnormal-weight-ratio-prototype-v1-*`). It now appears in the
canonical sampler and permanent range suite. The original failure records
remain unchanged. `final-v2-*` freezes and verifies this revised source.

## Final validation record

Before the weight-quotient correction, the final V1 source passes 54,883
package-sandbox assertions, 12,096 ordinary checks across 576 profiles and
40 concrete compiler signatures on each runtime. All ordinary profiles retain
zero warmed allocations and matching RNG/permutation state. Full historical
velocities and scratch are retained; 216 profiles change floating-point bits
under the explicit stable-increment migration.

The complete package report comparison preserves 29/27 reports. One private
pair-function entry signature changed with the removal of logarithmic inputs;
its entire callee chain, failure signature and message remain identical.
All other complete trees match after checked unchanged-source line mapping.
`compiler-v1/explicit-signature-review-v1.json` records that distinction;
the strict identical-tree comparison's failure is retained rather than hidden.

The permanent sampled-expression matrix passes 26,624 assertions per runtime.
A separate held-out seed range passes another 106,496 checks on each runtime.
The final V2 sandbox reruns the permanent matrix with the new 72-assertion
weight-quotient regression. Exceptional allocation measurements and fixed
latency acceptance are separate from ordinary zero-allocation checks.

Final V2 passes all 54,955 sandbox assertions, 12,096 ordinary checks and
40 concrete signatures on each runtime. The complete compiler review retains
the same one private-entry distinction and otherwise identical full trees.
`final-v2-review.json` verifies all source, log and full-state artifact hashes.
The separate exceptional resource inventory passes 720 state-integrity checks
over 72 profiles per runtime (six range families, four precisions, and 2/32/128
particles). All three warmed allocation samples agree within each profile;
the largest measured cost is 38,728 bytes per pair. Every row and population
scaling measurement remains in `exceptional-resources-v1-review.json`.

Fixed latency acceptance remains open under the shared 1.5 screening ceiling.
The numerical repair can be reviewed and pushed independently of that larger
optimization gate, as with the preceding exact-identity repair; no latency or
complete-Coulomb-correctness claim is made. Population-mean repair and local
multispecies composition remain separately tracked work.
