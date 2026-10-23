# Legacy longitudinal current modal range

Status: 2026-10-23; fresh defect reproduced on Julia 1.12.7 and 1.13.2.
This bounded repair is separate from the open forward/inverse transform and
arbitrary-grid conditioning work in `CONTINUITY_DIAGNOSTIC_RANGE_PLAN.md`.

## Pre-edit brief

- Mode: canonical arithmetic repair in the legacy EMPIC1D current solve.
- Requirement: finite represented density Fourier coefficients, finite nonzero
  derivative wavenumber and finite nonzero timestep must retain a representable
  completed longitudinal current coefficient. Keep ordinary operation order.
- Baseline: both native runtimes pass 24 defect-characterization checks at
  exact `2993f5a`. Four-point cosine density increments produce `[NaN,-Inf,NaN,Inf]`
  or all-zero currents although the independent complete-expression references
  are finite and nonzero. Ordinary controls and source-state checks pass.
- **VERIFIED cause:** `src/models/full_pic.jl:481` materializes `cb[m]/dt`
  before dividing by `ikm`. The recorded first quotient overflows or underflows.
  The already accepted `_spectral_exact_projection` produces the finite correct
  coefficient for those same stored operands.
- Contract: preserve the forward/inverse FFTs, null-mode masking, physical DC
  current, deposit chronology, source arrays, nonfinite behavior and public API.
  Do not claim this slice repairs unrepresentable Fourier coefficients or
  inverse-transform intermediates.
- Canonical owners: retain the legacy two-quotient fast path in a small modal
  helper beside existing spectral projection helpers in `electrostatic.jl`;
  delegate range-sensitive completion to `_spectral_exact_projection`.
  `_esirkepov_Jx!` supplies the nonzero derivative mode and original timestep.
- Compatibility: no new projection algebra, normalization, solver state,
  checkpoint schema, dependency, or general EMPIC arithmetic change.
- Strategy: inspect both component quotients at both stages using existing
  range predicates. Eligible exceptional inputs use the shared exact owner;
  ordinary and ineligible inputs retain the prior expression and result.
- Verification: independent rational modal rounding, both signs/components,
  overflow/underflow/opposite factor order, true final overflow and zero;
  independent four-point current oracles, DC/Nyquist and source-state controls;
  both complete EM suites, concrete/full-package compiler checks and formatter.

## Cost and acceptance gates

Compare complete legacy-current outputs and allocations against the exact parent
over Float32/64, empty and populated species, mobile/immobile ions, ordinary
signed charge changes and grid lengths 4, 5, 16, 128 and 1024. Warm each path,
alternate paired timing order and retain every sample. Require identical ordinary
outputs and allocation bytes; candidate/baseline median ratio must be <=1.5
for each ordinary profile. Component fast-path calls must remain allocation-free.
Exceptional calls retain numerical and source-state checks and report measured
allocations without presenting those observations as general upper bounds.
Any failed gate remains a failed candidate; no retry-until-green acceptance.

## Progress

- [x] Reproduce both range failures and ordinary controls on both runtimes.
- [x] Implement the shared-owner modal completion and permanent regression in isolation.
- [x] Pass numerical, neighboring, compiler, format and resource gates.
- [x] Review the exact source, preserve held changes, commit and push separately.
- [ ] Require all six exact-SHA GitHub jobs.

Fixing only infinities misses the zero-result witness. Multiplying `dt*k` first
can itself overflow or underflow before the completed quotient. Routing every
ordinary case through a differently ordered projection can change established
bits. The regression and differential gates must distinguish all three.

Persistent evidence:
`../HybridPlasmaPIC-audit-evidence/2026-10-23/legacy-current-mode`, including
`witness.jl`, both baseline logs and `baseline-commands.json`.

## Retained first candidate

V1 passes 1,438 focused assertions on Julia 1.12.7. Its first complete-current
cost matrix retains 193 passes and seven latency failures across 40 profiles;
all ordinary output, state and allocation comparisons pass. The seven observed
ratios are 1.5163–1.7606 against the unchanged 1.5 limit. Latest-Julia cost and
the subsequent integration gates were not run for this rejected candidate.
Its exact sources, manifest, driver and logs remain in `rejected-v1`.

The generated `v1-current112.llvm:90` contains a call to the new modal helper
inside the Fourier-mode loop. V2 explicitly inlines that helper without changing
its arithmetic. Whether this meets the cost target is an untested hypothesis
until the same paired matrix is run on the revised source. The generated call
alone does not establish a timing cause.

V2 also fails the same minimum-Julia cost gate: 197 passes and three latency
failures (1.5590–1.6932). All state, output and byte comparisons still pass;
`rejected-v2` retains its complete source and measured matrix. It is not promoted.
V3 combines the four independent, pure range predicates with Boolean `|` rather
than short-circuit `||`, preserving their truth value and every floating-point
operation. The eligibility check still guards the exact completion. Its cost
effect remains to be measured against the unchanged gate.

V3 retains 197 passes and three minimum-Julia latency failures, so it is also
unaccepted and archived in `rejected-v3`. V4 groups the two successive quotient
magnitudes belonging to each source component: their minimum must be normal
unless that original component is zero. Both complex intermediates must still
be finite. This removes redundant component-zero tests while preserving the
two-stage range condition. The four original predicates provide a differential
guard check over hardware edge values and random stored bit patterns before
the same public numerical and cost gates run. No arithmetic or target changes.

V4 passes 433,206 guard-decision comparisons and 1,438 numerical assertions on
each runtime. Its minimum-Julia cost matrix nevertheless retains 199 passes
and one latency failure: Float32, 1,024 cells, 16 particles and mobile ions
measure 1.6226 times the parent median. All outputs, states and allocation bytes
match. The complete matrix and exact sources remain in `rejected-v4`; the
latest cost lane and subsequent full neighbors/compiler gates have not run.
The unchanged 1.5 gate remains unsatisfied.

## V5 direct component guards

V5 replaces `min(abs(first), abs(second)) >= floatmin(T)` with the conjunction
of the two direct comparisons. For finite operands these predicates are
equivalent; the existing finite guards reject nonfinite intermediates before
accepting the ordinary path. The zero-source exemptions, eligibility guard,
quotient arithmetic and exact completion remain unchanged. This is a proposed
cost reduction, not an established timing cause or measured improvement.
Repeat the independent decision matrix and numerical regressions on both
runtimes, then the unchanged 40-profile cost gate once on this source. Retain
all failures; only a passing candidate may proceed to full neighboring gates.

V5 passes the decision/numerical matrices again on both runtimes and all 200
ordinary output, state, allocation and latency assertions on each runtime.
All 40 profiles per runtime retain identical outputs and zero measured bytes
in all 120 warmed samples per owner. Maximum observed median ratios are
1.4821 on Julia 1.12.7 and 1.4626 on 1.13.2, within the unchanged 1.5 gate.
`cost-commands-v5.json`, both complete TOML matrices and `source-v5.json` retain
the frozen source and every sample. Full neighboring/compiler/exceptional gates
remain pending; this candidate is not yet promoted.

## Final local acceptance

Both frozen V5 native lanes pass all four final commands: 37,916 complete
electromagnetic-neighbor assertions, six concrete compiler checks, the complete
package compiler comparison, and 96 exceptional output/state assertions per
runtime. Package trees preserve every baseline group and its order (30 reports
on Julia 1.12.7, 27 on 1.13.2) after checked source-line shifts and generated
type-variable names. Existing reports remain open; none were suppressed.

The 16 exceptional profiles per runtime retain all 48 warmed allocation samples:
18,448–24,128 bytes on minimum Julia and 17,248–22,480 on latest Julia. These
are observations for the recorded four-point inputs, not general resource bounds.
The ordinary 40-profile gates and zero-allocation results remain as recorded
above. All 273 frozen source/test/extension files still match `source-v5.json`.
`neighbor-commands-v5-*.json`, `compiler-comparison.json`, and
`exceptional-review-v5.json` retain the final command, source and log evidence.
Pushed separately as `7b50b51`; all 358 committed Julia blobs match the validated
candidate and held hashes are unchanged. [All six exact-SHA CI jobs pass](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39905400737).
Complete logs verify 616,222 / 616,223 passing assertions on Julia 1.12.7 /
1.13.2. Full compiler trees match the reviewed metadata-repair CI baseline
(23 / 22 reports). All six log hashes and the downloaded archive hash are
retained in `ci/39905400737/complete-log-review.json` within the evidence root.
