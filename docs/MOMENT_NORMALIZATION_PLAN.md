# Complete per-species first moments

Status: MOMENT-001 is pushed as `240f77c` and its complete GitHub run is green
on actual Julia 1.12.7 and 1.13.2. Native compiler checks also pass with
test-defined shape methods loaded. Remaining arithmetic boundaries follow below.

## Scope and pre-edit brief

- Mode: repair range-sensitive deposition and normalization.
- Acceptance: retain representable complete per-species density, momentum and
  current when an intermediate product, accumulation or volume loses range.
  Keep ordinary arithmetic, generic raw deposition, scratch contents and input
  validation behavior.
- Baseline: both native versions reproduce 12 volume and 24 numerator/stencil/
  current witness assertions. An eight-assertion control rejects predividing
  weights: doing so overflows before a finite CIC split. See the
  [baseline evidence](CHARGE_DEPOSITION_CONSERVATION_PLAN.md#confirmed-prerequisite-cell-volume-range).
- Cause: **VERIFIED** at `deposit_scalar!`, `momentum!` and `current!`: raw
  products/sums and cell volume are rounded before subsequent normalization or
  species charge can restore range. Both EMPIC current paths separately form
  the same products; the legacy longitudinal mean has an independent reduction.
- Invariants: read finite original operands for complete-expression replay;
  never recover from already-rounded scratch. Validate before modifying public
  moment destinations. Particle arrays stay unchanged. Generic raw deposition
  and ordinary finite cancellation retain their existing rounded contract.
- Owners: one `_deposit_scalar_kernel!` shared by serial/threaded deposition,
  `_finish_moment!` in `src/coupling/moments.jl`, and the existing EMPIC current
  owners in `src/models/full_pic.jl`. Existing exact conversion/rounding helpers
  remain authoritative. No new public API or package dependency.
- Compatibility: public moments divide then multiply by charge; EMPIC multiplies
  by charge/volume. Retain both ordinary orders. Work buffers still contain the
  last weighted velocity, including overflow. Unsupported symbolic coefficients
  and nonfinite operands retain visible legacy arithmetic.
- Strategy: specialize range tracking out of generic raw deposition; checked
  moment calls retain eligibility and product/accumulation/narrowing flags.
  Worker flags have private storage and are reduced after workers join. Eligible
  range loss triggers one rational cell replay per component, followed by the
  existing rounding owner. Track each individual cell width before its product.
- Verification: exact rational/rounding-cell oracles, supported custom shapes,
  retained BigFloat coefficient precision, symbolic/nonfinite controls, views,
  preflight/no-mutation tests, all dimensions and native precisions, scratch and
  EM callers, compiler/Aqua/format gates, ordinary baseline/resource comparisons,
  both native versions and GitHub CI.

The checked/unchecked private policy keeps one stencil owner. The replay is
bounded to once per component and has O(grid-cell) accumulator storage, with
operand-dependent BigInt payloads. Exceptions propagate. It is not a retry,
clamp, skipped input or second public implementation.

Tests reject five incomplete repairs: normalizing weights before stencil
splitting; repairing only volume; forming momentum before a rescuing charge;
checking only the final stored value; and ignoring loss inside worker partials.

## Local verification

The isolated source is `/tmp/hybrid-moment-normalization-v1`; baseline/final
hashes and the six-file patch accompany its handoff. The final arithmetic with
symbolic/BigFloat eligibility passes **1,953 permanent assertions and 64
diagnostic assertions on each native version**. The 1.12 run precedes a two-line
worker-state presence guard; the final 1.13 run includes that guard and uses
four threads. Logs: `/tmp/hybrid-moment-normalization-quality112.log` and
`/tmp/hybrid-moment-normalization-final-quality113.log`.

The final guard passes on 1.12 with **1,956 storage/preflight, 754 threaded,
49 worker-status, 42 concrete JET, 18 ordinary baseline-comparison and seven
eligibility-specialization assertions**. Package JET and Aqua pass on each
version; all 22/21 terminal diagnostic messages match their baseline with
multiplicity. Logs: `...final-quality{112,113}.log` and
`...diagnostics{112,113}.log`. The additional fixture-loaded runs pass four
fixture checks and all three JET checks on each version, exit 0. All 22/21
terminal messages still match baseline; logs are `...fixture-quality{112,113}.log`
and `...fixture-diagnostics{112,113}.log`.

Earlier source passes **54,445 affected assertions on 1.12** and **54,908 on
1.13**, including hybrid, electrostatic and both EM families. The first 1.12
harness completed 52,345 assertions before a relative-include harness error;
the corrected include-context run passes the remaining 2,100. Those logs remain
`...final112.log`, `...remaining112.log` and `...final113.log`; they are not
misrepresented as post-eligibility/guard runs.

All 18 final ordinary configurations match archived `a2154cf` bit for bit.
Serial generic deposition, density, momentum and current allocate 0 B.
Threaded generic deposition uses 2,816–19,456 B versus 2,944–19,776 B at baseline;
threaded density uses 2,944–19,584 B. Final shared-host minimum-time ratios are
0.774–1.215 for generic serial, 0.773–1.857 for generic threaded, and 0.839–2.086
for normalized moments. These measurements establish neither a speedup nor a
portable overhead ceiling. Earlier cold current cases with 16/256 particles
and 8/64/512 cells allocate 0.198–15.198 MB; they are allocation totals, not
peak RSS. Complete samples are retained in the resource and final-quality logs.

The eligibility audit found two candidate regressions before promotion:
symbolic π cannot be converted to `Rational{BigInt}`, and `issubnormal` does not
accept BigFloat. The final predicates preserve symbolic legacy behavior and
support stored BigFloat coefficients with independent narrowing/rounding-cell
tests. A subsequent package gate found optional worker-status indexing through
`Nothing`; the presence guard removes all four added reports without disabling
analysis. Earlier failed diagnostics remain available in the handoff.

## Remaining boundaries

- Complete EMPIC charge and current combinations are separate repairs, now
  pushed and awaiting their own GitHub acceptance: [MOMENT-002](EMPIC_CHARGE_COMBINATION_PLAN.md)
  and [MOMENT-003](EMPIC_CURRENT_COMBINATION_PLAN.md).
- Hybrid bulk-velocity ratios, other charge-weighted species combinations and
  MPI cross-rank reductions retain their separate complete-expression boundaries.
  Those consumers are not certified by this per-species repair.
- Pressure's second moments/covariance and generic raw deposition retain their
  separate range boundaries. Ordinary finite cancellation stays rounded.
- The public FourierGrid constructor supports Float32/64. BigFloat coverage here
  concerns retained stencil coefficients; third-party float coefficients and
  stateful custom shapes are not given unverified guarantees.

- [x] Native focused, neighboring, preflight and resource gates.
- [x] Native package/concrete compiler and Aqua checks before fixture extensions.
- [x] Package analysis after the new test-defined shape methods are loaded.
- [x] Promote the exact verified patch and review the final diff; source/test
  hashes match the candidate, and both new test includes remain registered.
- [x] Commit/push `240f77c` and observe complete GitHub acceptance.

[GitHub run 38741178196](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/38741178196) passes all six jobs. Package suites pass
357,258 assertions on 1.12.7 and 357,259 on 1.13.2; package JET counts remain
22/21. Formatting and all three MPI ranks also pass. [Job-log evidence](CI_RELIABILITY_PLAN.md#acceptance-run-at-240f77c).
