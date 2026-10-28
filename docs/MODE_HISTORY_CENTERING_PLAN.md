# Mode-history regression centering repair

Revision-scoped implementation record. The pre-edit findings and candidate
acceptance states below are historical; see the [audit tracker](DEEP_DEBUG_TRACKER.md)
for current integration and remaining work.

The correction landed as `4f26e1d`; [all six exact-SHA jobs pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/40206440280)
(checked 2026-10-29). The archived complete compiler review records unchanged
24 / 21-report forests after checked source-line mappings and 790,608 /
790,609 assertions on actual Julia 1.12.7 / 1.13.2. The original `839c109`
run remains failed.

## Pre-edit brief

- Mode: repair the canonical validation-only least-squares owner.
- Requirement: a constant complex history has zero fitted rate, residual and
  conditional standard error. Adding a representable constant observation
  offset preserves slopes and residual statistics. Keep every existing fit
  tolerance, input guard and phase/window convention.
- Evidence: `839c109` CI run `40194825415` has five green jobs and twelve
  Julia 1.12 failures at `test/test_mode_history.jl:23`; log-amplitude RMS is
  about 1.34e-13 against the unchanged 512-epsilon bound. The full failed log
  archive is retained. Independent constant-history witnesses fail 110 of 210
  checks on both local runtimes; 42 of 70 histories have a false residual,
  reaching 1.5916157281026244e-12. Source, commands, logs and per-input values
  are retained in `mode-fit-centering/before-v1-*` and
  `constant-public-before-v1-*` in the audit evidence directory.
- Cause **VERIFIED**: `_mode_line_fit` forms `sum(values)/n` at the absolute
  observation scale, then reconstructs `center+slope*x`. For each constant
  witness, the recorded mean error accounts for its residual. The exact
  constant-history OLS slope and residual are zero. The centered calculation
  must not reintroduce the absolute offset when evaluating residuals.
- Owner and callers: `validation/mode_history.jl` owns both full-window and
  half-window phase/log-magnitude fits. Case 34 and the provisional case 35
  consume these diagnostics; production solvers and raw histories are unchanged.
- Compatibility: retain signatures, return fields, guards and mathematical OLS
  semantics. Improved diagnostic floating-point bits are intentional. Replay
  fit-derived artifacts; do not retain incorrect residuals for bit equality.
- Strategy: subtract the first observation before computing its mean,
  covariance and residual. Restore the offset only for the returned center.
  Use the existing array reduction; no new dependency or special constant case.
- Verification: both original synthetic/dense-OLS/alias/long-history suites,
  independent constant histories over sizes around reduction boundaries,
  exactly representable translated noisy observations, views and reversed
  input; unchanged allocation scaling and concrete compiler checks; replay
  retained production histories; actual package sandbox and formatting; commit,
  push and all six exact-SHA CI jobs with full compiler-log review.

## Wrong patches

1. Loosening the RMS tolerance leaves constant histories with invented noise;
   exact-zero public witnesses reject it.
2. Correcting only the mean still rounds reconstructed fitted observations at
   the absolute scale; translated noisy histories check all residual statistics.
3. Special-casing constants leaves general offset sensitivity; independent
   dense least squares and translated nonconstant observations reject it.

## Progress

- [x] Preserve full CI failure logs and reproduce through an independent route.
- [x] Trace canonical full/half-window owner and affected artifacts.
- [x] Implement centered observation arithmetic and permanent regressions.
- [x] Complete numerical, resource/compiler and production-history replay gates.
- [x] Commit and push `4f26e1d`; six exact-SHA jobs and complete compiler review pass.

The separate [Coulomb population-mean repair](COULOMB_WEIGHT_MEAN_PLAN.md)
also landed, with its own exact-SHA acceptance; this centering repair does not
change that physical contract or its retained witnesses.

## Validation results

Both Julia 1.12.7 and 1.13.2 pass all 3,050 fitter assertions, including the
unchanged original synthetic tolerances, 280 exact-zero constant-history checks
and 128 translated-noise checks. The independent public witness passes all
210 checks after failing 110 before repair. Four concrete compiler signatures
and 64 existing workspace/input checks pass; the two-vector-plus-512-byte
allocation ceiling is unchanged.

Three executed wrong-patch controls retain the parent, correct only the mean,
or special-case constants. Each fails 31 translated-statistics comparisons on
each runtime. These failures verify that the permanent regression distinguishes
the full centered expression from a constant-only or mean-only patch.

Replaying all 40 frozen production histories passes 681 fixed-gate checks on
each runtime. A separate Python OLS/unwrapping implementation agrees on all
360 fitted statistics per runtime. Every raw-history hash remains unchanged;
this is a diagnostic replay, not a claim that new solver trajectories were run.
Formatter 1.0.62 passes AST/idempotence checks and leaves both tested Julia files
byte-identical. `focused-v1-*`, `mutation-v1-*`, `replay-v1-*` and
`independent-replay-v1.json` retain commands, inputs and results.

Both actual package test sandboxes pass 12,997 assertions across fitter,
production initialization/history, frequency extraction and configured Aqua
checks. All 94 `src`/`ext` Julia files match the parent byte for byte; the
validation helper is not part of the package's compilation surface. Its own
concrete compiler checks pass. Complete exact-SHA compiler review is accepted
as recorded above; the failed historical run remains failed evidence.
