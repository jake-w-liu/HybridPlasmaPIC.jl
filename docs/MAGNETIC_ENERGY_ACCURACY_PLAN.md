# Magnetic energy accuracy — ENERGY-010

Status (2026-10-23): reconstructed V5 and its composition onto `7be1916`
pass their final local gates on Julia 1.12.7 and 1.13.2. The repair is pushed
separately as `5c2796b`; all 357 committed Julia files match the candidate and
held hashes remain unchanged. Its [exact-commit CI](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/actions/runs/39897058507)
is pending; formatting is green. Complete-budget composition is a separate open repair.

The earlier `/tmp` artifacts below are unavailable; their historical counts
do not validate the reconstructed source. New evidence is retained at
`../HybridPlasmaPIC-audit-evidence/2026-10-23/field-accuracy`.

## Current reconstruction and integration contract

The shared Float32 owner certifies a wide complete-integral reduction and uses
bounded exact square bins when both outward bounds do not round to the same
value. It serves magnetic, legacy/general electromagnetic and inertia energy;
electromagnetic normalization retains `E + c²B`, while inertia retains linear
`de2`. Float64 ordinary arithmetic and the public total's association remain
unchanged. The independent exact rational oracle checks final rounding cells,
including midpoint tails, zero and overflow.

The final isolated V5 source passes 18,862 focused assertions, 756 cost
assertions across 164 profiles, 12 mutant checks, 50 concrete compiler checks
and 210 concurrency assertions per runtime. Maximum ordinary median ratios are
1.322575 and 1.244453, within the unchanged 1.5 gate. Ordinary component calls
allocate zero; public-budget allocations match the parent. The largest measured
large magnetic midpoint completion allocates 51,600 bytes, below 256 KiB.
All earlier numerical-harness and latency failures remain retained; no timing
cause is attributed to external load.

Before composition, root checked all nine owned hashes, eight final command
logs, the complete patch hash and 343 unchanged Julia files against `b251cb7`.
The current-parent integration will preserve the accepted interval, migration
and seam changes, insert the permanent accuracy test in the existing runner,
then run the same numerical/resource/concrete gates and complete package
compiler comparisons on the exact composed source. Source ownership and
held radiation hashes remain explicit promotion gates.

## Current composed acceptance

The frozen current-parent candidate passes 21,313 focused assertions on each
runtime, including the accepted interval-validation tests. Each runtime also
passes 756 resource assertions across 164 profiles, 12 independent mutant
checks, 50 concrete compiler checks and 210 concurrency assertions. All 18
formatting/parsed-semantics checks pass. Complete single-thread package compiler
trees and their ordering match `7be1916`: 30 reports on Julia 1.12.7 and 27 on
1.13.2; no diagnostics were suppressed.

All 160 ordinary profiles satisfy the unchanged 1.5 median-ratio gate. Maximum
Float32/Float64 ratios are 1.396129/1.191699 on minimum Julia and
0.931504/1.357780 on latest Julia. Component allocations remain zero, and
public-budget allocations match the parent. All four large midpoint profiles
meet the 256 KiB allocation gate; the measured maximum is 51,600 bytes. These
are paired measurements, not a claim that every call is faster.

Root verified all 11 native command/log hashes, 272 frozen source/test/extension
files and 348 unchanged committed Julia files. The independent review traces
the UInt128 bin bound, nonnegative wide-reduction error bound, exact coefficient
and spacing retention, rounding-cell oracle, midpoint and nonfinite paths,
unchanged Float64 owners and the complete final diff. The new test remains
registered alongside all accepted migration, seam and interval tests.

Persistent integration evidence is at
`../HybridPlasmaPIC-audit-evidence/2026-10-23/field-accuracy-combined`:
`commands.json`, `compiler-comparison.json`, `resource-and-source-review.json`,
the full native logs and per-profile TOML records. Earlier failed candidates
remain in `field-accuracy`; their failures were not removed or reclassified.

## Verified evidence

**VERIFIED:** on exact source `ffdb522aa858622a645c9ace104726ab54782493`,
`magnetic_energy` loses percent-level accuracy for ordinary large Float32
fields. Both Julia 1.12.7 and 1.13.2 produce the same results:

| Grid | Stored Bz | Public energy | Exact rational integral | Relative error |
| --- | --- | ---: | ---: | ---: |
| 128³ | Uniform 1.25 | 1,593,835.5 | 1,638,400 | 2.72000% |
| 128³ | Alternating 1 and 1.25 | 1,323,507.5 | 1,343,488 | 1.48721% |
| 160³ | Uniform 1.25 | 3,092,971.5 | 3,200,000 | 3.34464% |
| 160³ | Alternating 1 and 1.25 | 2,572,787.5 | 2,624,000 | 1.95170% |

Bx and By are zero; each stored cell spacing is one. References use exact
integer counts of the stored field levels and rational squares, independently
of the production reduction. All 32³/80³ Float32 controls and all eight
Float64 controls are exact. Each native run passes 88 assertions describing
these witnesses and controls; these are defect characterization checks, not
acceptance of the inaccurate result.

**VERIFIED reduction cause:** the public results equal a standalone sequential
Float32 square sum times one half in every measured profile.
`src/integrators/extrapolated_leapfrog.jl:1288–1298` accumulates this sum in
the field type; its documented range certificate retains that ordinary sum.
The exact reference shows accumulated rounding error even though all values
remain finite and normal. Range certification alone does not bound diagnostic
accuracy.

Evidence directory: `/tmp/hybrid-magnetic-accuracy-review-20261020`.
`probe.jl` contains the independent oracle and reduction control;
`baseline-source-hashes.json` pins the source; `results.json` records both
successful commands; `accuracy-summary.json` records all 32 profiles and log
hashes. No source change accompanied this assessment.

## Shared owners and compatibility

Source inspection at `ffdb522` identifies three public field reductions:
hybrid `magnetic_energy` at `extrapolated_leapfrog.jl:1288`, legacy
`em_field_energy(::EMPIC1D)` at `full_pic.jl:1573`, and general
`em_field_energy(::EMPIC)` at `full_pic.jl:2218`. All delegate normalization
and range completion to `_em_field_energy_finish` at `full_pic.jl:1515`;
its rational fallback shares `_field_square_sum_exact` at line 1481.
`_electron_inertia_energy` at `conservation.jl:10` has a separate ordinary
reduction with a linear coefficient, sharing the exact square sum.
These are inspected dependencies; the present numerical witnesses cover
hybrid magnetic energy only.

The ordinary-output checks in `test/test_hybrid_field_energy_range.jl:143`,
`test/test_em_energy.jl:173` and `test/inertia_energy/range.jl:87` intentionally
require the prior arithmetic's bits. For any owner changed by the accuracy
repair, its Float32 check must instead compare with the independent complete
integral. Preserve each allocation, input-ownership, nonfinite and Float64
compatibility assertion. Budget tests and benchmarks must use the newly
accepted component owner without disguising the component's own cost.

## Implementation contract and remaining gates

- [x] Reproduce against exact main on minimum and latest stable Julia.
- [x] Map shared hybrid, legacy/general electromagnetic and inertia quadratic
  reductions, coefficient/volume normalization and existing range fallbacks.
- [ ] Derive an accumulation error bound and choose a shared Float32 reduction
  that meets it. Float64 products are a candidate because each Float32 square
  is representable there; summation and normalization still need a complete
  error analysis. No grid-size cutoff or fixture-specific branch.
- [ ] Implement in the canonical owner, retaining validation, nonfinite
  behavior, zero handling and complete-input range repair. Improving ordinary
  Float32 bits is an intentional contract change; Float64 behavior must be
  separately checked. Update the range-only compatibility documentation.
- [ ] Add permanent independent rational oracles for uniform, alternating and
  irregular fields, coefficient/volume extremes, cancellation in budget totals,
  output midpoints, overflow and subnormal boundaries. Retain all previous
  field/inertia/PIC and full-budget checks.
- [ ] Measure latency and allocations on both native runtimes over small and
  large 1D–3D grids, including the four failing profiles. Set explicit numerical
  and cost gates before implementation; preserve failed measurements.
- [ ] Check concrete and complete compiler reports, native integration and
  concurrency, formatting and migration. Commit/push the accepted change
  separately to `main` and verify all six jobs for that exact SHA.

No accuracy implementation or performance improvement is claimed yet.

## Exact Float32 reduction prototype — pre-edit brief

The first isolated experiment targets the shared exact square reduction.
It does not modify the frozen complete-budget candidate.

Finite Float32 values have a significand of at most 24 bits. Squaring that
integer needs at most 48 bits. The stored exponent gives 254 possible square
scales; subnormals share the smallest normal's scale. With at most three
Int64-sized field arrays, each exponent group's integer sum is below 2^113,
so UInt128 bins can retain it exactly. Combining those bins at scale 2^-298
then constructs the rational square sum. These format/count derivations must
be checked against independently converted rationals and synthetic capacity
limits before integration. Reject unsupported inputs explicitly in the
prototype; inspect existing nonfinite behavior before changing a shared owner.

If this primitive passes, the intended public path first accumulates Float32
squares in Float64 and propagates an outward error bound through the coefficient
and cell volume. It accepts a result only when both endpoints round to the same
Float32 value. Uncertain cases use the exact binned reduction and existing exact
conversion. Float64 and other stored precisions retain their existing owner.
No public change is authorized by prototype results alone.

Acceptance targets are correctly rounded finite Float32 integrals, exact
Float64 baseline compatibility, zero ordinary warmed allocations, and paired
ordinary latency no greater than 1.5 times the prior owner over the declared
small/large matrices. Exceptional resource measurements must include large
exact-midpoint fields: bounded bin storage must avoid per-cell Rational work.
Keep nonfinite/signed-zero behavior, all existing range tests, full compiler
comparisons and native/CI integration. A failed target remains a failed
candidate until an explained implementation change is tested.

The isolated reduction at `/tmp/hybrid-field-accuracy-design-v2` passes 6,238
assertions on each Julia 1.12.7/1.13.2: every stored exponent group, signed
significand boundaries, random arrays in 1D–3D, rational equality, unchanged
inputs, empty arrays and the UInt128 capacity bound. V1's original test
mistakenly doubled the maximum significand in one capacity assertion; each
run retained 6,237 passes and that failure. V2 corrects only that assertion;
the reduction source is byte-identical and all assertions remain.
`resources.jl` passes 76 checks on each runtime for uniform, alternating and
exact-midpoint 32³–160³ fields. All 24 measured profiles allocate
48,088–48,304 bytes per exact reduction, below the predeclared 256 KiB gate.
Largest measured median times are 0.2145/0.2822 seconds on minimum/latest
Julia. These are absolute measurements; a paired public benchmark remains
required. `resource-results.json` and `resource-summary.json` retain commands,
profiles and log hashes. Public accuracy integration remains open.

The wide-accumulation prototype at
`/tmp/hybrid-field-accuracy-wide-v1-20261020` passes its numerical checks on
both Julia 1.12.7 and 1.13.2; all large-resource gates also pass.
Its brief derives a complete-integral bound for Float32 squares, exact stored
coefficients and one to three spacings; the verification driver compares the
bound and final rounding with independent per-cell rational integrals.
The driver also covers exponent groups, mixed scales, midpoint tails,
nonfinite delegation, ownership and ordinary allocations. Execution began after
the timing controls and thread-clock diagnostic. Each runtime passes 10,382
assertions: 10,206 independent enclosure/rounding checks, 170 midpoint/nonfinite/
zero checks and six ordinary allocation checks. The ordinary profiles allocate
zero bytes. This is prototype evidence; public integration and its full
acceptance remain open.

Its resource driver also covers uniform, alternating and exact-midpoint
32³–160³ fields. The predeclared complete-integral targets are correct rounding
for every result, zero warmed allocations for certified ordinary results,
and at most 256 KiB for each exact-midpoint evaluation. Each runtime passes
148 resource assertions across 12 profiles. Certified
ordinary evaluations allocate zero bytes; exact-midpoint evaluations allocate
51,464–51,632 bytes, below 256 KiB. Maximum measured median elapsed times are
0.2467/0.2564 seconds on minimum/latest Julia. `summary.json` retains all
commands, profiles and log hashes. These absolute observations do not satisfy
the still-pending paired public latency gate.
